"""Regime-Switching Best Fit (RSBF).

Online policy: behave exactly like Best Fit while collecting causal
statistics; when the past-item multiset gives strong evidence of a
structural regime in which an alternative any-fit variant empirically wins,
switch placement rule.  Modes are deterministic; switching is hysteresis-free
and monotone (once switched, stays).

Mode detection (all computed from PAST items only):
  triplet regime : fraction of past items in (cap/4, cap/2] >= theta_mid AND
                   fraction < cap/4 <= eps_small
  tail regime    : fraction of past items > 9cap/10 >= theta_tail  -> plain BF
                   (TBF loses here; keep BF)
Everything else -> BF.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def rsbf(sizes, cap, warmup, theta_mid, eps_small, tau_num, tau_den):
    """tau = tau_num/tau_den * cap defines 'tight' fit for TBF mode."""
    n = sizes.shape[0]
    resids = np.empty(n, dtype=np.int64)
    nb = 0
    hist_mid = 0        # past items in (cap/4, cap/2]
    hist_small = 0      # past items <= cap/4
    seen = 0
    thresh = tau_num * cap // tau_den
    mode = 0            # 0 = BF, 1 = TBF

    for i in range(n):
        s = int(sizes[i])
        if i >= warmup and mode == 0:
            if seen > 0 and hist_mid >= theta_mid * seen \
                    and hist_small <= eps_small * seen:
                mode = 1
        bb = -1
        br = cap + 1
        wb = -1
        wr = -1
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
        elif mode == 0:
            resids[bb] -= s
        else:
            if br <= thresh:
                resids[bb] -= s
            else:
                resids[wb] -= s
        # update causal stats
        if s > cap // 4:
            if s <= cap // 2:
                hist_mid += 1
        else:
            hist_small += 1
        seen += 1
    return nb


def tune_rsbf(train_streams, warmups=(20, 50, 100),
              thetas=(0.80, 0.90, 0.97), epss=(0.05, 0.15),
              taus=(0.02, 0.05)):
    best, best_cfg = None, None
    for w in warmups:
        for th in thetas:
            for e in epss:
                for t in taus:
                    tot = sum(rsbf(np.asarray(s, dtype=np.int64), c, w, th, e,
                                   int(round(t * 1000)), 1000)
                              for s, c in train_streams)
                    if best is None or tot < best:
                        best, best_cfg = tot, (w, th, e, t)
    return best_cfg, best
