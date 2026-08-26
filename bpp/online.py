"""ONLINE bin packing: causal GP policies + classical baselines.

Items arrive one at a time in stream order; placement is irrevocable.
Features are CAUSAL: current item, open bins, and statistics of PAST items
(maintained via a sqrt-decomposition histogram -> O(sqrt(cap)) update/query).
"""
import numpy as np
from numba import njit

STACK_DEPTH = 12
N_CONSTS = 8

# opcodes (mirrored from gp.py to avoid cross-imports at njit time)
PUSH_FEAT, PUSH_CONST = 0, 1
ADD, SUB, MUL, DIV, MIN, MAX, ABS, NEG, INV = 2, 3, 4, 5, 6, 7, 8, 9, 10

FEAT_NAMES_ON = [
    "resid",            # 0
    "fit",              # 1  resid - s
    "perfect",          # 2
    "load_frac",        # 3
    "resid_frac",       # 4
    "size_frac",        # 5  current item / cap
    "progress",         # 6  i / n
    "past_eq_fit",      # 7  #past items with size == fit
    "past_le_fit",      # 8  #past items with size <= fit
    "mean_past_frac",   # 9
    "max_past_frac",    # 10
    "min_past_frac",    # 11
    "nbins_frac",       # 12 nbins / n
    "waste_frac",       # 13 (nbins*cap - placed) / cap
    "is_new_bin",       # 14
    "bin_idx_frac",     # 15
    "age_frac",         # 16 (i - opened) / n
    "dead_resid",       # 17 fit>0 and no past item size <= fit
    "le_fit_ratio",     # 18 past_le_fit / i
    "logclass",         # 19 -log2(size/cap), clipped [0,8]
]
NF_ON = len(FEAT_NAMES_ON)


@njit(cache=True, boundscheck=False)
def pack_online(sizes, cap, ops, consts):
    """Simulate GP policy on stream; returns nbins used."""
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
    stack = np.empty((STACK_DEPTH, n + 1), dtype=np.float64)

    for i in range(n):
        s = int(sizes[i])
        m = nb + 1
        prog = i / n
        waste = (nb * cap - placed) / cap

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

        # ---- run program once over candidates (vectorized over bins)
        sp = 0
        valid = True
        p = 0
        nops = ops.shape[0]
        while p < nops:
            op = ops[p, 0]
            arg = ops[p, 1]
            p += 1
            if sp >= STACK_DEPTH:
                valid = False
                break
            if op == PUSH_FEAT:
                for b in range(m):
                    stack[sp, b] = X[arg, b]
                sp += 1
            elif op == PUSH_CONST:
                cv = consts[arg % N_CONSTS]
                for b in range(m):
                    stack[sp, b] = cv
                sp += 1
            elif op <= MAX:   # binary
                if sp < 2:
                    valid = False
                    break
                for b in range(m):
                    y = stack[sp - 1, b]
                    x0 = stack[sp - 2, b]
                    if op == ADD:
                        z = x0 + y
                    elif op == SUB:
                        z = x0 - y
                    elif op == MUL:
                        z = x0 * y
                    elif op == DIV:
                        d = y
                        if -1e-12 < d < 1e-12:
                            d = 1e-12
                        z = x0 / d
                    elif op == MIN:
                        z = x0 if x0 < y else y
                    else:
                        z = x0 if x0 > y else y
                    stack[sp - 2, b] = z
                sp -= 1
            else:             # unary
                if sp < 1:
                    valid = False
                    break
                for b in range(m):
                    y = stack[sp - 1, b]
                    if op == ABS:
                        z = abs(y)
                    elif op == NEG:
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

    return nb


# ------------------------------------------------------------------ baselines
@njit(cache=True)
def online_nf(sizes, cap):
    nb = 0
    resid = 0
    for i in range(sizes.shape[0]):
        s = int(sizes[i])
        if s > resid:
            nb += 1
            resid = cap - s
        else:
            resid -= s
    return nb


@njit(cache=True)
def online_ff(sizes, cap):
    n = sizes.shape[0]
    resids = np.empty(n, dtype=np.int64)
    nb = 0
    for i in range(n):
        s = int(sizes[i])
        bb = -1
        for b in range(nb):
            if resids[b] >= s:
                bb = b
                break
        if bb < 0:
            bb = nb
            nb += 1
            resids[bb] = cap
        resids[bb] -= s
    return nb


@njit(cache=True)
def online_bf(sizes, cap):
    n = sizes.shape[0]
    resids = np.empty(n, dtype=np.int64)
    nb = 0
    for i in range(n):
        s = int(sizes[i])
        bb = -1
        br = cap + 1
        for b in range(nb):
            r = resids[b]
            if r >= s and r < br:
                br = r
                bb = b
        if bb < 0:
            nb += 1
            resids[nb - 1] = cap - s
        else:
            resids[bb] -= s
    return nb


@njit(cache=True)
def online_wf(sizes, cap):
    n = sizes.shape[0]
    resids = np.empty(n, dtype=np.int64)
    nb = 0
    for i in range(n):
        s = int(sizes[i])
        bb = -1
        br = -1
        for b in range(nb):
            r = resids[b]
            if r >= s and r > br:
                br = r
                bb = b
        if bb < 0:
            nb += 1
            resids[nb - 1] = cap - s
        else:
            resids[bb] -= s
    return nb


def harmonic_count(sizes, cap, m):
    """Harmonic-Greedy(m): class k holds sizes in (cap/(k+2), cap/(k+1)];
    each class packed independently & greedily into bins of capacity
    cap/(k+1).  Online-causal, deterministic."""
    s = np.asarray(sizes, dtype=np.int64)
    total = 0
    for k in range(m):
        hi = cap // (k + 1)
        lo = cap // (k + 2) + 1
        mask = (s >= lo) & (s <= hi)
        cnt_k = int(mask.sum())
        if cnt_k == 0:
            continue
        sizes_k = np.sort(s[mask])[::-1]
        open_r = []
        binsk = 0
        for sz in sizes_k:
            placed = False
            for j in range(len(open_r)):
                if open_r[j] >= sz:
                    open_r[j] -= sz
                    placed = True
                    break
            if not placed:
                open_r.append(hi - sz)
                binsk += 1
        total += binsk
    return total
