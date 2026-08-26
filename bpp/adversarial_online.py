"""Adversarial stream search for ONLINE policies: find streams where the
champion uses many more bins than Best Fit (paired gap normalized by L1).

Hill-climbing over integer streams: perturb k item sizes per proposal,
accept improvements of objective = (bins_H - bins_BF)/L1.
"""
import numpy as np

from .online import pack_online, online_bf


def gap(sizes, cap, ops, consts):
    h = pack_online(sizes, cap, ops, consts)
    b = online_bf(sizes, cap)
    l1 = max(1.0, sizes.sum() / cap)
    return (h - b) / l1


def hill_climb_stream(ops, consts, cap=1000, n=200, iters=250, seed=0,
                      init=None):
    rng = np.random.default_rng(seed)
    if init is None:
        s = rng.integers(1, cap + 1, size=n).astype(np.int64)
    else:
        s = np.asarray(init, dtype=np.int64).copy()
        n = len(s)
    cur = gap(s, cap, ops, consts)
    best_s, best_v = s.copy(), cur
    for t in range(iters):
        s2 = s.copy()
        T = max(2, int(cap * 0.15 * (1 - t / iters)))
        k = int(rng.integers(1, 5))
        idx = rng.choice(n, size=k, replace=False)
        s2[idx] = np.clip(s2[idx] + rng.integers(-T, T + 1, size=k), 1, cap)
        v = gap(s2, cap, ops, consts)
        if v >= cur:
            s, cur = s2, v
            if v > best_v:
                best_v, best_s = v, s2.copy()
    return best_s, float(best_v)


def multi_restart(ops, consts, restarts=6, cap=1000, n=200, iters=200,
                  seed0=0):
    bs, bv = None, -1e18
    for r in range(restarts):
        s, v = hill_climb_stream(ops, consts, cap, n, iters, seed0 + r)
        if v > bv:
            bv, bs = v, s
    return bs, bv
