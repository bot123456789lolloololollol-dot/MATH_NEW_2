"""Tuned simple-policy baselines for online packing.

TBF(tau)  -- 'tight Best Fit' in the spirit of FunSearch's reported discovery:
place like Best Fit only when the tightest available fit is <= tau*cap;
otherwise place into the feasible bin whose POST-placement residual is
LARGEST (worst-fit behaviour) so as not to create small residuals.
PFB      -- perfect-fit first (exact fill priority), else Best Fit.
Both are causal online policies; tau is tuned on TRAINING streams only and
frozen for evaluation.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def tbf(sizes, cap, tau):
    n = sizes.shape[0]
    resids = np.empty(n, dtype=np.int64)
    nb = 0
    thresh = tau * cap
    for i in range(n):
        s = int(sizes[i])
        bb = -1
        br = cap + 1          # tightest fit
        wb = -1
        wr = -1               # largest post-residual
        for b in range(nb):
            r = resids[b]
            if r >= s:
                f = r - s
                if f < br:
                    br = f
                    bb = b
                if f > wr:
                    wr = f
                    wb = b
        if bb < 0:
            nb += 1
            resids[nb - 1] = cap - s
        elif br <= thresh:
            resids[bb] -= s
        else:
            resids[wb] -= s
    return nb


@njit(cache=True)
def pfb_bf(sizes, cap):
    n = sizes.shape[0]
    resids = np.empty(n, dtype=np.int64)
    nb = 0
    for i in range(n):
        s = int(sizes[i])
        pb = -1
        bb = -1
        br = cap + 1
        for b in range(nb):
            r = resids[b]
            if r >= s:
                f = r - s
                if f == 0:
                    pb = b
                    break
                if f < br:
                    br = f
                    bb = b
        if pb >= 0:
            resids[pb] = 0
        elif bb >= 0:
            resids[bb] -= s
        else:
            nb += 1
            resids[nb - 1] = cap - s
    return nb


def tune_tbf(streams, taus=(0.01, 0.02, 0.05, 0.08, 0.12, 0.18, 0.25)):
    """Return tau minimizing total bins across training streams."""
    best_tau, best_total = None, None
    for tau in taus:
        tot = 0
        for s, cap in streams:
            tot += tbf(np.asarray(s, dtype=np.int64), cap, tau)
        if best_total is None or tot < best_total:
            best_total, best_tau = tot, tau
    return best_tau, best_total


def tune_harmonic(streams, ms=range(2, 17)):
    from .online import harmonic_count
    best_m, best_total = None, None
    for m in ms:
        tot = sum(harmonic_count(s, cap, m) for s, cap in streams)
        if best_total is None or tot < best_total:
            best_total, best_m = tot, m
    return best_m, best_total
