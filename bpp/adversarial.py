"""Adversarial instance search: find instances where candidate H wastes
maximally more bins than FFD/BFD ensemble.

Method: stochastic hill-climbing over multisets of integer sizes (fixed n),
perturbing one item size at a time, objective = bins(H) - min(bins_FFD,bins_BFD)
(normalized by LB to compare across sizes of instance).  Used both as an
adversarial test of champions and as a diagnostic tool.
"""
import numpy as np

from .core import ffd, bfd, l2_lower_bound
from .gp import pack_gp


def gap_objective(sizes, cap, ops, consts):
    sd = np.sort(sizes)[::-1].copy()
    h = pack_gp(sd, cap, ops, consts)[0]
    base = min(ffd(sd, cap), bfd(sd, cap))
    lb = max(1, l2_lower_bound(sizes, cap))
    return (h - base) / lb


def hill_climb_adversarial(ops, consts, cap=1000, n=60, iters=400,
                           seed=0, init="uniform", temp_schedule=None):
    rng = np.random.default_rng(seed)
    if init == "uniform":
        s = rng.integers(1, cap + 1, size=n).astype(np.int64)
    elif init == "ffd_trap":
        reps = n // 3
        s = np.concatenate([np.array([cap // 2 + 3, cap // 4 + 6, cap // 4 + 3])
                            * (cap // 1000 or 1) for _ in range(reps)]).astype(np.int64)[:n]
        s = np.clip(s, 1, cap)
    else:
        s = np.asarray(init, dtype=np.int64).copy()

    cur = gap_objective(s.astype(np.int32), cap, ops, consts)
    T0 = cap // 8
    best_s, best_v = s.copy(), cur
    for t in range(iters):
        T = max(1, int(T0 * (1 - t / iters)))
        s2 = s.copy()
        k = int(rng.integers(1, 4))
        idx = rng.choice(n, size=k, replace=False)
        s2[idx] = np.clip(s2[idx] + rng.integers(-T, T + 1, size=k), 1, cap)
        v = gap_objective(s2.astype(np.int32), cap, ops, consts)
        if v >= cur:
            s, cur = s2, v
            if v > best_v:
                best_v, best_s = v, s2.copy()
    return best_s.astype(np.int32), best_v


def random_restart_adversarial(ops, consts, restarts=8, cap=1000, n=60,
                               iters=300, seed0=0):
    best = None
    bv = -1e9
    for r in range(restarts):
        s, v = hill_climb_adversarial(ops, consts, cap=cap, n=n, iters=iters,
                                      seed=seed0 + r,
                                      init="uniform" if r % 3 else "ffd_trap")
        if v > bv:
            bv, best = v, s
    return best, bv
