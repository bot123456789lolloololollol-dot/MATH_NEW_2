"""Steady-state evolutionary search for placement rules.

Every individual is evaluated on the full training set each generation
(affordable at ~2-4 s/generation on a 20-core machine), avoiding the
subset-vs-full fitness comparability problem.
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


def run_evolution(out_prefix="experiments/run1", pop_size=300, gens=80,
                  offspring_per_gen=120, n_train=240, seed=12345,
                  elitism=6, parsimony=4e-4,
                  train_names=None, log_every=5):
    rng = np.random.default_rng(seed)
    t0 = time.time()

    dsgen = np.random.default_rng(seed + 7)
    names = train_names or list(GENERATORS.keys())
    dataset = make_mixed_dataset(dsgen, n_train, names)
    ev = ParallelEvaluator(dataset)

    def penalized(f, p):
        return f + parsimony * p.ops.shape[0]

    pop = [ProgSpec(rand_prog(rng), rand_consts(rng)) for _ in range(pop_size)]
    raw_fits = ev.evaluate(pop)
    fits = [penalized(f, p) for f, p in zip(raw_fits, pop)]

    history = []
    best_i = int(np.argmin(fits))
    best, best_fit = pop[best_i], fits[best_i]
    print(f"[init] best={best_fit:.5f}  {program_to_str(best)[:110]}", flush=True)

    for g in range(gens):
        children = []
        while len(children) < offspring_per_gen:
            r = rng.random()
            if r < 0.55:
                ch = crossover(tournament(pop, fits, 5, rng),
                               tournament(pop, fits, 5, rng), rng)
            elif r < 0.9:
                ch = mutate(tournament(pop, fits, 5, rng), rng)
            else:
                ch = ProgSpec(rand_prog(rng), rand_consts(rng))
            if parsimony * ch.ops.shape[0] > 0.05:
                continue  # reject monsters early
            children.append(ch)
        craw = ev.evaluate(children)
        cfits = [penalized(f, p) for f, p in zip(craw, children)]

        pool = list(zip(pop, fits)) + list(zip(children, cfits))
        pool.sort(key=lambda pf: pf[1])
        uniq, used = [], set()
        for p, f in pool:
            sig = signature(p)
            if sig in used:
                continue
            used.add(sig)
            uniq.append((p, f))
            if len(uniq) >= pop_size:
                break
        pop = [p for p, _ in uniq]
        fits = [f for _, f in uniq]

        gbest = int(np.argmin(fits))
        improved = fits[gbest] < best_fit - 1e-12
        if improved:
            best_fit = fits[gbest]
            best = pop[gbest]
        history.append({"gen": g, "best": float(min(fits)),
                        "median": float(np.median(fits)),
                        "len_best": int(best.ops.shape[0])})
        if g % log_every == 0 or g == gens - 1:
            print(f"[gen {g:3d}] best={min(fits):.5f} "
                  f"median={np.median(fits):.5f} "
                  f"elapsed={time.time()-t0:.0f}s", flush=True)

    ev.close()

    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    np.savez_compressed(out_prefix + "_champion.npz", ops=best.ops,
                        consts=best.consts)
    with open(out_prefix + "_log.json", "w") as fh:
        json.dump({
            "seed": seed, "pop_size": pop_size, "gens": gens,
            "n_train": n_train, "parsimony": parsimony,
            "train_names": names,
            "best_fitness": float(best_fit),
            "history": history,
            "champion_str": program_to_str(best),
        }, fh, indent=2)
    print("champion:", program_to_str(best))
    print("fitness:", best_fit)
    return best, best_fit


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="experiments/run1")
    ap.add_argument("--pop", type=int, default=300)
    ap.add_argument("--gens", type=int, default=80)
    ap.add_argument("--offspring", type=int, default=120)
    ap.add_argument("--train", type=int, default=240)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()
    run_evolution(args.out, args.pop, args.gens, args.offspring, args.train,
                  args.seed)
