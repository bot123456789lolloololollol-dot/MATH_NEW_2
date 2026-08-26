#!/usr/bin/env python3
"""Algorithms for P||Cmax. All take an integer processing-time array p (unsorted)
and m; all return (assign, cmax) where assign is a list of m lists of original
job indices, and cmax is the makespan (int).

Implementations follow PREREGISTERED.md:
  B1 LPT       - Graham (1969), SIAM J. Appl. Math. 17(2):416-429.
  B2 MULTIFIT  - Coffman, Garey, Johnson (1978), SIAM J. Comput. 7(1):1-17.
                 (deviation: initial bracket [LB, sum p], up to 40 bisections)
  B3 LPT+LS    - LPT start; steepest-descent local search to local optimum over
                 {relocation (jump), pairwise interchange (swap)} involving >=1
                 critical machine (loses no strictly-improving move: any strict
                 improvement must reduce EVERY critical machine's load).
  C1 ROL       - bounded-lookahead rollout LPT (our mechanism; see SPEC.md).
Deterministic: no randomness anywhere; documented tie-breaks only.
"""
from typing import List, Tuple

import numpy as np


# ---------------------------------------------------------------- helpers ----
def cmax_of(assign: List[List[int]], p: np.ndarray) -> int:
    return int(max((sum(int(p[j]) for j in mach) for mach in assign), default=0))


def lb_bound(p: np.ndarray, m: int) -> int:
    """max(max_j p_j, ceil(sum_j p_j / m))"""
    return int(max(int(np.max(p)), -(-int(np.sum(p)) // m)))


def _argmin_load(loads: List[int]) -> int:
    best, bi = None, 0
    for i, l in enumerate(loads):
        if best is None or l < best:
            best, bi = l, i
    return bi


def _order_desc(p: np.ndarray) -> np.ndarray:
    """Stable descending order of job indices (ties by original index asc)."""
    return np.argsort(-p, kind="stable")


# ------------------------------------------------------------------ B1 LPT ---
def lpt(p: np.ndarray, m: int) -> Tuple[List[List[int]], int]:
    order = _order_desc(p)
    loads = [0] * m
    assign: List[List[int]] = [[] for _ in range(m)]
    for j in order.tolist():
        i = _argmin_load(loads)
        loads[i] += int(p[j])
        assign[i].append(int(j))
    return assign, max(loads)


# ------------------------------------------------------------ B2 MULTIFIT ---
def ffd_feasible(p_sorted_desc: List[int], m: int, cap: int) -> bool:
    """First-fit-decreasing packing into m bins of capacity `cap`."""
    loads = [0] * m
    for d in p_sorted_desc:
        for i in range(m):
            if loads[i] + d <= cap:
                loads[i] += d
                break
        else:
            return False
    return True


def multifit(p: np.ndarray, m: int) -> Tuple[List[List[int]], int]:
    order = _order_desc(p)
    pdesc = [int(p[j]) for j in order.tolist()]
    lo = max(pdesc[0], -(-sum(pdesc) // m))  # LB: infeasible below LB
    hi = sum(pdesc)                          # trivially feasible
    cap = hi
    while lo < hi:
        mid = (lo + hi) // 2
        if ffd_feasible(pdesc, m, mid):
            hi = mid
        else:
            lo = mid + 1
        cap = hi
    # rebuild explicit assignment at final capacity (first-fit-decreasing)
    loads = [0] * m
    assign: List[List[int]] = [[] for _ in range(m)]
    for idx, d in enumerate(pdesc):
        for i in range(m):
            if loads[i] + d <= cap:
                loads[i] += d
                assign[i].append(int(order[idx]))
                break
        else:
            raise AssertionError("multifit reconstruction failed")
    return assign, max(loads)


# ------------------------------------------------------------- B3 LPT + LS ---
def ls_no_improving_move(assign: List[List[int]], p: np.ndarray) -> bool:
    """Independent naive checker (no critical-machine restriction): verifies no
    single relocation or pairwise interchange strictly reduces the makespan."""
    cmax = cmax_of(assign, p)
    m = len(assign)
    locs = {}
    for i, mach in enumerate(assign):
        for j in mach:
            locs[int(j)] = i
    loads = [sum(int(p[j]) for j in mach) for mach in assign]
    n = len(p)
    # relocations
    for j in range(n):
        src = locs[j]
        for dst in range(m):
            if dst == src:
                continue
            if max(l + (int(p[j]) if i == dst else (-int(p[j]) if i == src else 0))
                   for i, l in enumerate(loads)) < cmax:
                return False
    # interchanges
    for j1 in range(n):
        a = locs[j1]
        for j2 in range(j1 + 1, n):
            b = locs[j2]
            if a == b:
                continue
            na = loads[a] - int(p[j1]) + int(p[j2])
            nb = loads[b] - int(p[j2]) + int(p[j1])
            ok = True
            for i, l in enumerate(loads):
                v = na if i == a else (nb if i == b else l)
                if v >= cmax:
                    ok = False
                    break
            if ok:
                return False
    return True


def lpt_ls(p: np.ndarray, m: int, max_passes: int = 100000) -> Tuple[List[List[int]], int]:
    assign, _ = lpt(p, m)
    loads = [sum(int(p[j]) for j in mach) for mach in assign]
    pnp = np.asarray(p, dtype=np.int64)
    pmembers = [np.asarray(mach, dtype=np.int64) for mach in assign]

    def refresh_members():
        for i in range(m):
            pmembers[i] = np.asarray(assign[i], dtype=np.int64)

    for _pass in range(max_passes):
        cmax = max(loads)
        crit = [i for i in range(m) if loads[i] == cmax]
        best = None  # (val, kind_order, c, i, j_pos, k_pos_or_None)

        for c in crit:
            jc = pmembers[c]
            if jc.size == 0:
                continue
            pj = pnp[jc]
            for i in range(m):
                if i == c:
                    continue
                ji = pmembers[i]
                base = max(
                    (loads[t] for t in range(m) if t != c and t != i),
                    default=-1,
                )
                # --- relocations j: c -> i
                vr = np.maximum(np.maximum(base, loads[c] - pj),
                                loads[i] + pj)
                idxs = np.nonzero(vr < cmax)[0]
                if idxs.size:
                    pos = int(idxs[np.argmin(vr[idxs])])
                    val = int(vr[pos])
                    cand = (val, 0, c, i, pos, None)
                    if best is None or cand < best:
                        best = cand
                # --- interchanges j in c <-> k in i
                if ji.size == 0:
                    continue
                pk = pnp[ji]
                nc = loads[c] - pj[:, None] + pk[None, :]
                ni = loads[i] - pk[None, :] + pj[:, None]
                val = np.maximum(np.maximum(nc, ni), base)
                mask = val < cmax
                if not mask.any():
                    continue
                flat = val.ravel()
                fpos = int(np.argmin(flat))
                rpos, kpos = divmod(fpos, ji.size)
                cand = (int(flat[fpos]), 1, c, i, int(rpos), int(kpos))
                if best is None or cand < best:
                    best = cand

        if best is None:
            break  # local optimum w.r.t. full neighborhood (see docstring)

        val, kind, c, i, rpos, kpos = best
        if val >= cmax:  # safety: never accept non-improving
            break
        if kind == 0:
            j = int(assign[c][rpos])
            assign[c].pop(rpos)
            assign[i].append(j)
            loads[c] -= int(p[j])
            loads[i] += int(p[j])
        else:
            j = int(assign[c][rpos])
            k = int(assign[i][kpos])
            assign[c][rpos] = k
            assign[i][kpos] = j
            loads[c] += int(p[k]) - int(p[j])
            loads[i] += int(p[j]) - int(p[k])
        refresh_members()

    return assign, max(loads)


# ------------------------------------------------------- C1 ROL (candidate) ---
def _simulate_lpt_finish(loads: List[int], rest: List[int], m: int) -> int:
    """Final makespan if remaining jobs (already sorted desc) run under plain LPT."""
    cur = loads[:]
    for d in rest:
        best, bi = None, 0
        for i in range(m):
            l = cur[i]
            if best is None or l < best:
                best, bi = l, i
        cur[bi] += d
    return max(cur)


def rol(p: np.ndarray, m: int, K_mult: int) -> Tuple[List[List[int]], int]:
    """Bounded-lookahead rollout LPT. The first K = K_mult*m jobs (largest ones)
    are placed by trying every distinct-load machine, completing the schedule
    with plain LPT, and keeping the placement minimizing the rollout makespan.
    Remaining jobs: plain LPT. Tie-breaks: smaller resulting load, then lower
    machine index. Deterministic."""
    n = len(p)
    order = _order_desc(p)
    ps = [int(x) for x in p[order]]
    orig = order.tolist()
    loads = [0] * m
    assign: List[List[int]] = [[] for _ in range(m)]
    K = min(K_mult * m, n)
    for idx in range(K):
        d = ps[idx]
        rest = ps[idx + 1:]
        # dedupe symmetric machine choices by current load
        seen: dict = {}
        for mi in range(m):
            if loads[mi] not in seen:
                seen[loads[mi]] = mi
        choices = sorted(seen.values())
        best = None  # (rollout_cmax, resulting_load, machine)
        for mi in choices:
            nl = loads[mi] + d
            loads[mi] = nl
            val = _simulate_lpt_finish(loads, rest, m)
            loads[mi] -= d
            cand = (val, nl, mi)
            if best is None or cand < best:
                best = cand
        _, _, mi = best
        loads[mi] += d
        assign[mi].append(int(orig[idx]))
    for idx in range(K, n):
        mi = _argmin_load(loads)
        loads[mi] += ps[idx]
        assign[mi].append(int(orig[idx]))
    return assign, max(loads)


ALGORITHMS = {
    "lpt": lambda p, m: lpt(p, m),
    "multifit": lambda p, m: multifit(p, m),
    "lpt_ls": lambda p, m: lpt_ls(p, m),
}
