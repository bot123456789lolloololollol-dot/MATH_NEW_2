"""Round-2 evolution: island model + PAIRED fitness vs min(FFD,BFD).

Fitness = mean_i (H_i - base_i)/lb_i - parsimony*len.
A program equivalent to the best classical heuristic scores exactly 0;
only strictly-better packings push it negative.  Islands exchange elites
periodically; fresh immigrants per generation counter diversity collapse.
"""
import argparse
import json
import os
import time
import numpy as np

from .gp import (ProgSpec, rand_prog, rand_consts, mutate, crossover,
                 program_to_str)
from .instances import GENERATORS, make_mixed_dataset
from .evaluate import ParallelEvaluator


def tournament(pop, fits, k, rng):
    idx = rng.integers(0, len(pop), size=k)
    best = idx[np.argmin([fits[i] for i in idx])]
    return pop[best]


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


def run_evolution2(out_prefix="experiments/run2", islands=4, pop_per_island=80,
                   gens=150, offspring_frac=0.5, n_train=300, seed=2026,
                   elitism=4, parsimony=3e-4, migrate_every=15,
                   train_names=None, log_every=10, immigrant_frac=0.15):
    rng = np.random.default_rng(seed)
    t0 = time.time()

    dsgen = np.random.default_rng(seed + 7)
    names = train_names or list(GENERATORS.keys())
    dataset = make_mixed_dataset(dsgen, n_train, names)
    ev = ParallelEvaluator(dataset, mode="paired")

    def penal(f, p):
        return f + parsimony * p.ops.shape[0]

    pops, fitss = [], []
    for _ in range(islands):
        pop = [ProgSpec(rand_prog(rng), rand_consts(rng))
               for _ in range(pop_per_island)]
        raw = ev.evaluate(pop)
        pops.append(pop)
        fitss.append([penal(f, p) for f, p in zip(raw, pop)])

    history = []
    gbest_fit = min(min(fs) for fs in fitss)

    for g in range(gens):
        for isl in range(islands):
            pop, fits = pops[isl], fitss[isl]
            order = np.argsort(fits)
            n_off = int(offspring_frac * len(pop))
            children = []
            while len(children) < n_off:
                ch = breed(pop, fits, rng)
                children.append(ch)
            n_imm = max(1, int(immigrant_frac * n_off))
            children += [ProgSpec(rand_prog(rng), rand_consts(rng))
                         for _ in range(n_imm)]
            craw = ev.evaluate(children)
            cfits = [penal(f, p) for f, p in zip(craw, children)]
            pool = list(zip(pop, fits)) + list(zip(children, cfits))
            pool.sort(key=lambda pf: pf[1])
            uniq, used = [], set()
            for p, f in pool:
                sig = signature(p)
                if sig in used:
                    continue
                used.add(sig)
                uniq.append((p, f))
                if len(uniq) >= pop_per_island:
                    break
            pops[isl] = [p for p, _ in uniq]
            fitss[isl] = [f for _, f in uniq]

        # migration: exchange island champions
        if (g + 1) % migrate_every == 0 and islands > 1:
            champs = [pops[i][0] for i in range(islands)]
            for i in range(islands):
                donor = champs[(i + 1) % islands]
                tgt = pops[i]
                # replace worst non-elite
                worst = int(np.argmax(fitss[i][elitism:])) + elitism
                pops[i][worst] = donor
                fitss[i][worst] = penal(ev.evaluate([donor])[0], donor)

        cur_best = min(min(fs) for fs in fitss)
        if cur_best < gbest_fit - 1e-12:
            gbest_fit = cur_best
        allf = [f for fs in fitss for f in fs]
        flat = int(np.argmin(allf))
        bi, bj = flat // pop_per_island, flat % pop_per_island
        champion = pops[bi][bj]
        history.append({"gen": g, "best": float(cur_best),
                        "median": float(np.median(allf)),
                        "len_champ": int(champion.ops.shape[0])})
        if g % log_every == 0 or g == gens - 1:
            print(f"[gen {g:3d}] best={cur_best:+.6f} "
                  f"median={np.median(allf):+.6f} "
                  f"elapsed={time.time()-t0:.0f}s "
                  f"champ={program_to_str(champion)[:70]}", flush=True)

    ev.close()
    allf = [(f, i, j) for i, fs in enumerate(fitss) for j, f in enumerate(fs)]
    bestf, bi, bj = min(allf)
    best = pops[bi][bj]

    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    np.savez_compressed(out_prefix + "_champion.npz", ops=best.ops,
                        consts=best.consts)
    with open(out_prefix + "_log.json", "w") as fh:
        json.dump({
            "seed": seed, "islands": islands, "pop_per_island": pop_per_island,
            "gens": gens, "n_train": n_train, "parsimony": parsimony,
            "train_names": names, "mode": "paired",
            "best_fitness": float(bestf),
            "history": history,
            "champion_str": program_to_str(best),
        }, fh, indent=2)
    print("champion:", program_to_str(best))
    print("paired fitness:", bestf)
    return best


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="experiments/run2")
    ap.add_argument("--islands", type=int, default=4)
    ap.add_argument("--pop", type=int, default=80)
    ap.add_argument("--gens", type=int, default=150)
    ap.add_argument("--train", type=int, default=300)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()
    run_evolution2(args.out, args.islands, args.pop, args.gens,
                   n_train=args.train, seed=args.seed)
