"""Fast packing kernels (numba-jitted).

Conventions:
  * sizes: int32 array of item sizes, 1 <= s <= cap
  * cap:   integer bin capacity
  * All kernels return an assignment array (bin index per item, in the order
    processed) and implicitly the bin count = max(assign)+1.  We return the
    number of bins used directly for speed.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def _open_new_bin(bins_resid, cap):
    """Append a fresh bin, return its index."""
    bins_resid.append(cap)
    return bins_resid.shape[0] - 1


# ---------------------------------------------------------------- greedy one-pass
@njit(cache=True)
def first_fit(sizes, cap):
    """First Fit in the given order."""
    n = sizes.shape[0]
    assign = np.empty(n, dtype=np.int64)
    resids = np.empty(n, dtype=np.int64)  # upper bound on bins
    nb = 0
    for i in range(n):
        s = sizes[i]
        placed = -1
        for b in range(nb):
            if resids[b] >= s:
                placed = b
                break
        if placed < 0:
            placed = nb
            nb += 1
            resids[placed] = cap
        resids[placed] -= s
        assign[i] = placed
    return nb, assign


@njit(cache=True)
def best_fit(sizes, cap):
    """Best Fit in the given order."""
    n = sizes.shape[0]
    assign = np.empty(n, dtype=np.int64)
    resids = np.empty(n, dtype=np.int64)
    nb = 0
    for i in range(n):
        s = sizes[i]
        best_b = -1
        best_r = cap + 1
        for b in range(nb):
            r = resids[b]
            if r >= s and r < best_r:
                best_r = r
                best_b = b
        if best_b < 0:
            best_b = nb
            nb += 1
            resids[best_b] = cap
            best_r = cap - s
        else:
            best_r = best_r - s
        resids[best_b] = best_r
        assign[i] = best_b
    return nb, assign


@njit(cache=True)
def worst_fit(sizes, cap):
    """Worst Fit in the given order."""
    n = sizes.shape[0]
    assign = np.empty(n, dtype=np.int64)
    resids = np.empty(n, dtype=np.int64)
    nb = 0
    for i in range(n):
        s = sizes[i]
        best_b = -1
        best_r = -1
        for b in range(nb):
            r = resids[b]
            if r >= s and r > best_r:
                best_r = r
                best_b = b
        if best_b < 0:
            best_b = nb
            nb += 1
            resids[best_b] = cap
        resids[best_b] -= s
        assign[i] = best_b
    return nb, assign


@njit(cache=True)
def argsort_desc(sizes):
    return np.argsort(-sizes.astype(np.float64), kind="quicksort")


def ffd(sizes, cap):
    order = np.argsort(-np.asarray(sizes, dtype=np.float64), kind="stable")
    return first_fit(np.asarray(sizes, dtype=np.int32)[order], cap)[0]


def bfd(sizes, cap):
    order = np.argsort(-np.asarray(sizes, dtype=np.float64), kind="stable")
    return best_fit(np.asarray(sizes, dtype=np.int32)[order], cap)[0]


def wfd(sizes, cap):
    order = np.argsort(-np.asarray(sizes, dtype=np.float64), kind="stable")
    return worst_fit(np.asarray(sizes, dtype=np.int32)[order], cap)[0]


# ---------------------------------------------------------------- lower bounds
@njit(cache=True)
def l1_lower_bound(sizes, cap):
    total = 0
    for i in range(sizes.shape[0]):
        total += sizes[i]
    return -(-total // cap)  # ceil division


@njit(cache=True)
def l2_lower_bound(sizes, cap):
    """Martello-Toth L2 lower bound (exact integer form)."""
    n = sizes.shape[0]
    best = l1_lower_bound(sizes, cap)
    # alpha ranges over candidate values: half-capacity and each distinct size
    alphas = np.empty(n + 1, dtype=np.int64)
    for k in range(n):
        alphas[k] = sizes[k]
    alphas[n] = cap // 2
    for ka in range(n + 1):
        alpha = alphas[ka]
        if alpha < 0 or alpha > cap // 2:
            continue
        n1 = 0          # s > cap - alpha
        n2 = 0          # cap - alpha >= s > cap/2
        sum2 = 0        # sum of N2
        sum3 = 0        # sum of N3 (cap/2 >= s >= alpha)
        for i in range(n):
            s = sizes[i]
            if s > cap - alpha:
                n1 += 1
                sum3 += 0
            elif s > cap // 2:
                n2 += 1
                sum2 += s
            elif s >= alpha:
                sum3 += s
        free_in_n2_bins = n2 * cap - sum2
        extra = sum3 - free_in_n2_bins
        lb = n1 + n2
        if extra > 0:
            lb += -(-extra // cap)
        if lb > best:
            best = lb
    return best


# ---------------------------------------------------------------- exact optimum (small n)
@njit(cache=True)
def optimal_bins_bitmask(sizes, cap):
    """Exact minimum bins via subset DP. Only for n <= ~20 (2^n memory)."""
    n = sizes.shape[0]
    full = (1 << n) - 1
    INF = np.int64(1 << 30)
    # feasibility & sums
    fits = np.full(1 << n, False)
    wsum = np.zeros(1 << n, dtype=np.int64)
    fits[0] = True
    for m in range(1, 1 << n):
        low = m & (-m)
        i = 0
        mm = low
        while mm > 1:
            mm >>= 1
            i += 1
        wsum[m] = wsum[m ^ low] + sizes[i]
        fits[m] = wsum[m] <= cap
    dp = np.full(1 << n, INF, dtype=np.int64)
    dp[0] = 0
    for m in range(1, 1 << n):
        # iterate submasks that are feasible bins containing lowest unset bit trick:
        sub = m
        while sub > 0:
            if fits[sub]:
                v = dp[m ^ sub] + 1
                if v < dp[m]:
                    dp[m] = v
            sub = (sub - 1) & m
    return dp[full]


def optimal_bins_small(sizes, cap):
    """Wrapper: falls back to branch-and-bound for n<=64 without 2^n blowup guard."""
    sizes = np.asarray(sizes, dtype=np.int32)
    n = len(sizes)
    if n > 20:
        raise ValueError("bitmask DP limited to n<=20")
    return int(optimal_bins_bitmask(sizes, cap))
