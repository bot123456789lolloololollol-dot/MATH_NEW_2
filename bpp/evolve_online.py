"""Evolution of ONLINE placement policies.

Fitness = mean over training streams of (bins_policy - bins_BF)/L1(stream).
Zero == ties Best Fit; negative == beats it on average (paired per stream).
Island model + immigrants, same stack-machine representation with causal
features.
"""
import argparse
import json
import os
import time
import numpy as np

from .online import pack_online
from .gp import ProgSpec, rand_prog, rand_consts, mutate, crossover, \
    program_to_str
from .streams import STREAM_GENS, make_stream_set


def tournament(pop, fits, k, rng):
    idx = rng.integers(0, len(pop), size=k)
    return pop[idx[np.argmin([fits[i] for i in idx])]]


def signature(p):
    return p.ops.tobytes()


def breed(pop, fits, rng, p_cross=0.5, mut_rate=0.25):
    r = rng.random()
    if r < p_cross:
        ch = crossover(tournament(pop, fits, 4, rng),
                       tournament(pop, fits, 4, rng), rng)
        if rng.random() < 0.5:
            ch = mutate(ch, rng, rate=mut_rate * 0.5)
        return ch
    return mutate(tournament(pop, fits, 4, rng), rng, rate=mut_rate)


# ---------------------------------------------------------------- worker side
_W = {}


def _init_worker(streams):
    _W["streams"] = streams


def _eval_one(payload):
    ops, consts = payload[:2]
    streams = _W["streams"]
    tot = 0.0
    for s, cap, bf, l1 in streams:
        nb = pack_online(s, cap, ops, consts)
        tot += (nb - bf) / l1
    return tot / len(streams)


def evaluate_pop(pool, progs):
    payloads = [(p.ops, p.consts) for p in progs]
    if pool is None:
        return [_eval_one(x) for x in payloads]
    return list(pool.map(_eval_one, payloads, chunksize=1))


def run(out_prefix="experiments/onrun1", islands=4, pop_per_island=80,
        gens=120, n_streams=200, seed=7, elitism=4, parsimony=3e-4,
        migrate_every=15, immigrant_frac=0.15, log_every=10,
        max_stream_len=800):
    from concurrent.futures import ProcessPoolExecutor
    rng = np.random.default_rng(seed)
    t0 = time.time()

    sgen = np.random.default_rng(seed + 11)
    streams = []
    for nm, s, cap in make_stream_set(sgen, n_streams):
        if len(s) > max_stream_len:
            keep = rng.choice(len(s), size=max_stream_len, replace=False)
            s = s[keep]
        sd = np.asarray(s, dtype=np.int64)
        bf = int(np.max(np.zeros(1)))  # placeholder replaced below
        streams.append((sd, cap))      # baselines computed in workers? no ->
        streams[-1] = (sd, cap)

    # precompute BF & L1 here (fast enough once)
    from .online import online_bf
    # JIT warm-up in the PARENT so workers load the on-disk cache instead of
    # racing to compile/write it themselves (Windows file-lock crashes).
    _w = np.array([3, 7], dtype=np.int64)
    pack_online(_w, 10, np.array([[0, 1]], dtype=np.int32), np.zeros(8))
    online_bf(_w, 10)
    prepared = []
    for sd, cap in streams:
        bf = online_bf(sd, cap)
        l1 = max(1, -(-int(sd.sum()) // cap))
        prepared.append((sd, int(cap), bf, l1))

    if (os.cpu_count() or 4) > 1:
        pool = ProcessPoolExecutor(max_workers=os.cpu_count() - 1,
                                   initializer=_init_worker,
                                   initargs=(prepared,))
    else:
        pool = None

    def penal(f, p):
        return f + parsimony * p.ops.shape[0]

    pops, fitss = [], []
    for _ in range(islands):
        pop = [ProgSpec(rand_prog(rng), rand_consts(rng))
               for _ in range(pop_per_island)]
        raw = evaluate_pop(pool, pop)
        pops.append(pop)
        fitss.append([penal(f, p) for f, p in zip(raw, pop)])

    history = []
    for g in range(gens):
        for isl in range(islands):
            pop, fits = pops[isl], fitss[isl]
            n_off = int(0.5 * len(pop))
            children = []
            while len(children) < n_off:
                ch = breed(pop, fits, rng)
                children.append(ch)
            children += [ProgSpec(rand_prog(rng), rand_consts(rng))
                         for _ in range(max(1, int(immigrant_frac * n_off)))]
            craw = evaluate_pop(pool, children)
            cfits = [penal(f, p) for f, p in zip(craw, children)]
            merged = list(zip(pop, fits)) + list(zip(children, cfits))
            merged.sort(key=lambda pf: pf[1])
            uniq, used = [], set()
            for p, f in merged:
                sig = signature(p)
                if sig in used:
                    continue
                used.add(sig)
                uniq.append((p, f))
                if len(uniq) >= pop_per_island:
                    break
            pops[isl] = [p for p, _ in uniq]
            fitss[isl] = [f for _, f in uniq]

        if (g + 1) % migrate_every == 0 and islands > 1:
            champs = [pops[i][0] for i in range(islands)]
            for i in range(islands):
                donor = champs[(i + 1) % islands]
                worst = int(np.argmax(fitss[i][elitism:])) + elitism
                pops[i][worst] = donor
                fitss[i][worst] = penal(evaluate_pop(pool, [donor])[0], donor)

        allf = [f for fs in fitss for f in fs]
        flat = int(np.argmin(allf))
        bi, bj = flat // pop_per_island, flat % pop_per_island
        champion = pops[bi][bj]
        history.append({"gen": g, "best": float(allf[flat]),
                        "median": float(np.median(allf)),
                        "len_champ": int(champion.ops.shape[0])})
        if g % log_every == 0 or g == gens - 1:
            print(f"[gen {g:3d}] best={allf[flat]:+.6f} "
                  f"median={np.median(allf):+.6f} "
                  f"elapsed={time.time()-t0:.0f}s "
                  f"champ={program_to_str(champion)[:70]}", flush=True)

    if pool is not None:
        try:
            pool.shutdown(wait=True, cancel_futures=True)
        except OSError:
            pass  # benign Windows queue-feeder race at interpreter exit
    bestf, bi, bj = min((f, i, j) for i, fs in enumerate(fitss)
                        for j, f in enumerate(fs))
    best = pops[bi][bj]

    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    np.savez_compressed(out_prefix + "_champion.npz", ops=best.ops,
                        consts=best.consts)
    with open(out_prefix + "_log.json", "w") as fh:
        json.dump({"seed": seed, "islands": islands,
                   "pop_per_island": pop_per_island, "gens": gens,
                   "n_streams": n_streams, "parsimony": parsimony,
                   "best_fitness": float(bestf), "history": history,
                   "champion_str": program_to_str(best)}, fh, indent=2)
    print("champion:", program_to_str(best))
    print("paired-vs-BF fitness:", bestf)
    return best


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="experiments/onrun1")
    ap.add_argument("--islands", type=int, default=4)
    ap.add_argument("--pop", type=int, default=80)
    ap.add_argument("--gens", type=int, default=120)
    ap.add_argument("--streams", type=int, default=200)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    run(args.out, args.islands, args.pop, args.gens, args.streams, args.seed)
