"""Property checks for online champions.

anyfit_violations: does the policy ever open a new bin while a feasible
open bin exists?  Any-Fit policies inherit classical worst-case bounds.

restricted_class_stats: behavior on streams with all items > cap/3
(at most two items per bin) -- the class where Best Fit admits a clean
exact analysis.
"""
import numpy as np
from numba import njit

from .online import NF_ON


@njit(cache=True, boundscheck=False)
def _simulate_track_anyfit(sizes, cap, ops, consts):
    """Returns (nbins, violations). violation = opened new bin although some
    open bin had residual >= s."""
    n = sizes.shape[0]
    resids = np.empty(n + 1, dtype=np.int64)
    loads = np.empty(n + 1, dtype=np.int64)
    opened = np.empty(n + 1, dtype=np.int64)
    BLK = 32
    nblk = (cap // BLK) + 2
    cnt = np.zeros(cap + 2, dtype=np.int64)
    blk = np.zeros(nblk, dtype=np.int64)
    nb = 0
    placed = 0
    acc_past = 0
    min_past = cap + 1
    max_past = 0
    X = np.empty((NF_ON, n + 1), dtype=np.float64)
    stack = np.empty((12, n + 1), dtype=np.float64)
    violations = 0

    for i in range(n):
        s = int(sizes[i])
        m = nb + 1
        prog = i / n
        waste = (nb * cap - placed) / cap
        feasible_exists = False
        for b in range(nb):
            if resids[b] >= s:
                feasible_exists = True
                break
        for b in range(m):
            if b < nb:
                r = resids[b]
                is_new = 0.0
                bidx = b / m
                age = i - opened[b]
            else:
                r = cap
                is_new = 1.0
                bidx = 1.0
                age = 0
            fit = r - s
            perfect = 1.0 if fit == 0 else 0.0
            peq = 0.0
            ple = 0.0
            if 0 <= fit <= cap:
                fi = fit
                peq = float(cnt[fi])
                kb = fi // BLK
                for k in range(kb):
                    ple += blk[k]
                for v in range(kb * BLK, fi + 1):
                    ple += cnt[v]
            dead = 1.0 if (fit > 0 and ple == 0.0 and i > 0) else 0.0
            lc = -np.log2(s / cap)
            if lc < 0.0:
                lc = 0.0
            if lc > 8.0:
                lc = 8.0
            X[0, b] = float(r)
            X[1, b] = float(fit)
            X[2, b] = perfect
            X[3, b] = loads[b] / cap
            X[4, b] = r / cap
            X[5, b] = s / cap
            X[6, b] = prog
            X[7, b] = peq
            X[8, b] = ple
            X[9, b] = (acc_past / i / cap) if i > 0 else 0.0
            X[10, b] = max_past / cap
            X[11, b] = (min_past / cap) if min_past <= cap else 0.0
            X[12, b] = nb / max(1, n)
            X[13, b] = waste
            X[14, b] = is_new
            X[15, b] = bidx
            X[16, b] = age / max(1, n)
            X[17, b] = dead
            X[18, b] = (ple / i) if i > 0 else 0.0
            X[19, b] = lc
        sp = 0
        valid = True
        p = 0
        nops = ops.shape[0]
        while p < nops:
            op = ops[p, 0]
            arg = ops[p, 1]
            p += 1
            if sp >= 12:
                valid = False
                break
            if op == 0:
                for b in range(m):
                    stack[sp, b] = X[arg, b]
                sp += 1
            elif op == 1:
                cv = consts[arg % 8]
                for b in range(m):
                    stack[sp, b] = cv
                sp += 1
            elif op <= 7:
                if sp < 2:
                    valid = False
                    break
                for b in range(m):
                    y = stack[sp - 1, b]
                    x0 = stack[sp - 2, b]
                    if op == 2:
                        z = x0 + y
                    elif op == 3:
                        z = x0 - y
                    elif op == 4:
                        z = x0 * y
                    elif op == 5:
                        d = y
                        if -1e-12 < d < 1e-12:
                            d = 1e-12
                        z = x0 / d
                    elif op == 6:
                        z = x0 if x0 < y else y
                    else:
                        z = x0 if x0 > y else y
                    stack[sp - 2, b] = z
                sp -= 1
            else:
                if sp < 1:
                    valid = False
                    break
                for b in range(m):
                    y = stack[sp - 1, b]
                    if op == 8:
                        z = abs(y)
                    elif op == 9:
                        z = -y
                    else:
                        d = y
                        if -1e-12 < d < 1e-12:
                            d = 1e-12
                        z = 1.0 / d
                    stack[sp - 1, b] = z
        if not valid or sp < 1:
            best_b = nb
            bf = cap + 1
            for b in range(nb):
                f = resids[b] - s
                if f >= 0 and f < bf:
                    bf = f
                    best_b = b
        else:
            for b in range(nb):
                if resids[b] < s:
                    stack[sp - 1, b] = np.nan
            best_b = -1
            best_score = 0.0
            have = False
            for b in range(m):
                sc = stack[sp - 1, b]
                if sc != sc:
                    continue
                if not have or sc > best_score:
                    have = True
                    best_score = sc
                    best_b = b
            if best_b < 0:
                best_b = nb
        if best_b == nb and feasible_exists:
            violations += 1
        if best_b == nb:
            resids[nb] = cap - s
            loads[nb] = s
            opened[nb] = i
            nb += 1
        else:
            resids[best_b] -= s
            loads[best_b] += s
        placed += s
        cnt[s] += 1
        blk[s // BLK] += 1
        acc_past += s
        if s > max_past:
            max_past = s
        if s < min_past:
            min_past = s
    return nb, violations


def anyfit_violations(ops, consts, streams):
    tot = 0
    for s, cap in streams:
        _, v = _simulate_track_anyfit(np.asarray(s, dtype=np.int64),
                                      int(cap), ops, consts)
        tot += v
    return tot
