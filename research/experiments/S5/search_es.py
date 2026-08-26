"""Parallel evolutionary minimization of comparator networks.

Modes:
  sort    : minimize lexicographic (size, depth) of a SORTING network on n
            wires, verified by exhaustive bitset simulation.
  select  : minimize size of a k-SELECTION network (bottom-k inputs land on
            wires 0..k-1, unordered; tail unconstrained), verified by numpy
            exhaustive simulation.

Workers run independent ES chains (different RNG seeds); improvements are
appended to per-worker JSONL logs and merged by the coordinator.
"""
import argparse
import glob
import json
import os
import random
import sys
import time

import netlib


# ------------------------------------------------------------------ fitness

class SortJudge:
    def __init__(self, n):
        self.n = n
        _, em = netlib.masks_for(n)
        self.em = em

    def valid(self, net):
        w = netlib.evaluate(net, self.n)
        return all(w[k] == self.em[k] for k in range(self.n))

    def error(self, net):
        return netlib.error_bits(net, self.n)

    def metrics(self, net):
        return (len(net), netlib.depth(net))


class SelectJudge:
    """bottom-k of the input must appear (as a set) on wires 0..k-1."""

    def __init__(self, n, k):
        import numpy as np
        self.np = np
        self.n, self.k = n, k
        total = 1 << n
        p = np.arange(total, dtype=np.uint32)
        self.inp = ((p[:, None] >> np.arange(n, dtype=np.uint32)[None, :])
                    & 1).astype(np.uint8)          # (P, n)
        srt = np.sort(self.inp, axis=1)
        self.want_lo = np.sort(srt[:, :k], axis=1)  # sorted bottom-k rows

    def valid(self, net):
        np = self.np
        w = self.inp.copy()
        for (i, j) in net:
            wi, wj = w[:, i].copy(), w[:, j].copy()
            np.minimum(wi, wj, out=w[:, i])
            np.maximum(wi, wj, out=w[:, j])
        got = np.sort(w[:, :self.k], axis=1)
        return bool(np.array_equal(got, self.want_lo))

    def error(self, net):
        np = self.np
        w = self.inp.copy()
        for (i, j) in net:
            wi, wj = w[:, i].copy(), w[:, j].copy()
            np.minimum(wi, wj, out=w[:, i])
            np.maximum(wi, wj, out=w[:, j])
        got = np.sort(w[:, :self.k], axis=1)
        return int((got != self.want_lo).sum())

    def metrics(self, net):
        d = netlib.depth(net)
        return (len(net), d)


def make_judge(mode, n, k=None):
    return SortJudge(n) if mode == "sort" else SelectJudge(n, k)


# ------------------------------------------------------------------ seeds

def load_seed_nets(mode, n, k, seed_dir="seeds"):
    nets = []
    # inject previously-verified champions when present
    champ = os.path.join("verified", f"SEL_{n}_{k}_best.json")
    if os.path.exists(champ):
        with open(champ) as fh:
            d = json.load(fh)
        if d.get("n") == n and d.get("k") == k:
            nets.append([tuple(e) for e in d["network"]])
    if mode == "sort":
        for f in glob.glob(os.path.join(seed_dir, "Sorters_Sort_%d_*.json" % n)):
            with open(f) as fh:
                d = json.load(fh)
            nets.append([tuple(e) for e in d["nw"]])
    else:
        for f in glob.glob(os.path.join(seed_dir, "Median_Median_%d_*.json" % n)):
            with open(f) as fh:
                d = json.load(fh)
            nw = [tuple(e) for e in d["nw"]]
            kk = (n - 1) // 2
            if k == kk:
                nets.append(nw)   # median nets directly match this k
            # truncate a full-sorter view: median net + leftover is not a
            # sorter; instead derive selectors from sorting seeds below.
    if mode == "sort":
        import baselines
        nets.append(baselines.batcher_oem_net(n))
        nets.append(baselines.insertion_net(n))
    else:
        # baseline: take best known sorter at n (if any) and use its prefix
        # property? A full sorter IS a k-selector (wires 0..k-1 hold sorted
        # bottom-k). Use smallest available sorter seed or Batcher.
        import baselines
        best_sorter = None
        sizes = []
        for f in glob.glob(os.path.join(seed_dir, "Sorters_Sort_%d_*.json" % n)):
            sizes.append(f)
        if sizes:
            fbest = min(sizes, key=lambda f: json.load(open(f))["L"])
            with open(fbest) as fh:
                d = json.load(fh)
            best_sorter = [tuple(e) for e in d["nw"]]
        else:
            best_sorter = baselines.batcher_oem_net(n)
        nets.append(best_sorter)
    return nets


# ------------------------------------------------------------------ mutations

def mutate(net, n, rng):
    m = rng.randrange(7)
    L = len(net)
    g = list(net)
    if m == 0 and L > 1:                       # delete random CE
        del g[rng.randrange(L)]
    elif m <= 2:                               # insert random CE
        i = rng.randrange(n - 1)
        j = rng.randint(i + 1, n - 1)
        g.insert(rng.randrange(L + 1), (i, j))
    elif m == 3 and L > 1:                     # move CE
        r = g.pop(rng.randrange(L))
        g.insert(rng.randrange(L + 1), r)
    elif m == 4:                               # retarget one endpoint
        r = rng.randrange(L)
        i, j = g[r]
        if rng.random() < 0.5:
            i = max(0, min(j - 1, i + rng.choice((-2, -1, 1))))
        else:
            j = min(n - 1, max(i + 1, j + rng.choice((-2, -1, 1))))
        if i < j:
            g[r] = (i, j)
    elif m == 5 and L > 2:                     # swap two positions
        a, b = rng.sample(range(L), 2)
        g[a], g[b] = g[b], g[a]
    else:                                      # reverse a segment
        if L > 3:
            a = rng.randrange(L - 1)
            b = rng.randint(a + 1, min(L, a + 8))
            g[a:b] = reversed(g[a:b])
    return g


def prune_fast(net, judge, max_rounds=2):
    cur = list(net)
    for _ in range(max_rounds):
        changed = False
        r = 0
        while r < len(cur):
            cand = cur[:r] + cur[r + 1:]
            if judge.valid(cand):
                cur = cand
                changed = True
            else:
                r += 1
        if not changed:
            break
    return cur


# ------------------------------------------------------------------ repair


def remove_repair(net, judge, n, rng, tries=26):
    """Delete a random CE, then greedily insert CEs that reduce the error.
    Returns a net that may still be invalid (caller judges by fitness)."""
    g = list(net)
    r = rng.randrange(len(g))
    del g[r]
    e = judge.error(g)
    if e == 0:
        return g
    pos = min(r, len(g))
    for _ in range(tries):
        i = rng.randrange(n - 1)
        j = rng.randint(i + 1, n - 1)
        cand = g[:pos] + [(i, j)] + g[pos:]
        ec = judge.error(cand)
        if ec < e:
            g, e = cand, ec
            pos = min(pos + rng.choice((0, 1)), len(g))
            if e == 0:
                break
    return g


def crossover(a, b, rng):
    """Prefix/suffix splice of two nets (positions independent)."""
    ca = rng.randrange(len(a) + 1)
    cb = rng.randrange(len(b) + 1)
    return list(a[:ca]) + list(b[cb:])


def transplant(a, b, rng, seg=4):
    """Copy a contiguous segment of b into a at a random position."""
    if len(b) < 1:
        return list(a)
    s = rng.randrange(len(b))
    t = min(len(b), s + rng.randint(1, seg))
    pos = rng.randrange(len(a) + 1)
    return list(a[:pos]) + list(b[s:t]) + list(a[pos:])


def fitness(judge, net):
    """Lexicographic (error, size): invalid nets are stepping stones."""
    return (judge.error(net), len(net))


# ------------------------------------------------------------------ ES chain

def es_chain(worker_id, mode, n, k, minutes, seed, out_dir,
             init_nets, target=None):
    rng = random.Random(seed)
    judge = make_judge(mode, n, k)
    t_end = time.time() + minutes * 60
    log_path = os.path.join(out_dir, f"worker_{worker_id}.jsonl")
    logf = open(log_path, "a")

    # champion pool: distinct valid nets as (metrics, net), Pareto-sorted
    pool = []
    seen = set()
    for net in init_nets:
        if judge.valid(net):
            key = tuple(net)
            if key not in seen:
                seen.add(key)
                pool.append((judge.metrics(net), list(net)))
    if not pool:
        return
    pool.sort(key=lambda x: x[0])

    # current individual for the (error, size) walk: start from best valid
    cur = list(pool[0][1])
    cur_f = fitness(judge, cur)

    best_m = pool[0][0]
    eval_count = 0

    while time.time() < t_end:
        eval_count += 1
        u = rng.random()
        if u < 0.25:
            child = remove_repair(cur, judge, n, rng)
        elif u < 0.55:
            child = mutate(cur, n, rng)
        elif u < 0.85 and len(pool) > 1:
            # recombination: splice or transplant from another pool member
            other = pool[rng.randrange(len(pool))][1]
            child = crossover(cur, other, rng) if rng.random() < 0.5 \
                else transplant(other, cur, rng)
            if rng.random() < 0.6:
                child = remove_repair(child, judge, n, rng, tries=14)
        elif u < 0.92:
            # double delete + repair
            child = list(cur)
            for _ in range(2):
                if len(child) > 2:
                    del child[rng.randrange(len(child))]
            child = remove_repair(child, judge, n, rng, tries=30)
        else:
            # restart walk from a random pool member (keeps diversity)
            cur = list(pool[rng.randrange(len(pool))][1])
            cur_f = fitness(judge, cur)
            continue
        f = fitness(judge, child)
        if f <= cur_f:
            if f < cur_f or rng.random() < 0.5:
                cur, cur_f = child, f
            if f[0] == 0:
                m = judge.metrics(child)
                key = tuple(child)
                if m < best_m or (m not in [p[0] for p in pool]):
                    if key not in seen:
                        seen.add(key)
                        pool.append((m, list(child)))
                        pool.sort(key=lambda x: x[0])
                        del pool[40:]
                    if m < best_m:
                        best_m = m
                        logf.write(json.dumps(
                            {"t": time.time(), "size": m[0],
                             "depth": m[1], "net": list(child)}) + "\n")
                        logf.flush()
                        if target and m[0] <= target:
                            break
        if eval_count % 4000 == 0 and pool:
            pb = prune_fast(pool[0][1], judge)
            if len(pb) < len(pool[0][1]):
                key = tuple(pb)
                if key not in seen:
                    seen.add(key)
                    pool[0] = (judge.metrics(pb), pb)
                    pool.sort(key=lambda x: x[0])
                    best_m = min(best_m, pool[0][0])
                    logf.write(json.dumps(
                        {"t": time.time(), "size": pool[0][1].__len__(),
                         "depth": pool[0][0][1], "net": pool[0][1],
                         "pruned": True}) + "\n")
                    logf.flush()
    # final flush of top-5 distinct
    flushed = set()
    for m, net in sorted(pool, key=lambda x: x[0])[:8]:
        key = tuple(net)
        if key not in flushed:
            flushed.add(key)
            logf.write(json.dumps({"t": time.time(), "size": len(net),
                                   "depth": netlib.depth(net),
                                   "final": True, "net": net}) + "\n")
    logf.close()


# ------------------------------------------------------------------ driver

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("sort", "select"), required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--k", type=int, default=None)
    ap.add_argument("--minutes", type=float, default=20)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--seed", type=int, default=20260826)
    ap.add_argument("--target", type=int, default=None,
                    help="stop a worker early when size <= target")
    ap.add_argument("--out", default="runs")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    init_nets = load_seed_nets(args.mode, args.n, args.k)
    print(f"mode={args.mode} n={args.n} k={args.k} "
          f"seeds={len(init_nets)} workers={args.workers}", flush=True)

    import multiprocessing as mp
    procs = []
    for w in range(args.workers):
        pargs = (w, args.mode, args.n, args.k, args.minutes,
                 args.seed + 1000 * w, args.out, init_nets, args.target)
        p = mp.Process(target=es_chain, args=pargs)
        p.start()
        procs.append(p)
    for p in procs:
        p.join()
    print("done", flush=True)


if __name__ == "__main__":
    main()
