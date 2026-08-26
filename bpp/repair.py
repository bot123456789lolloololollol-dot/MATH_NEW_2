"""Monotone repair heuristics on top of any constructive packing.

Move ("unload-and-repack"): choose r bins with smallest load, empty them,
redistribute their items into kept bins' residual space (first-fit, desc
order), pack leftovers among themselves via BFD in fresh bins; accept iff the
total bin count strictly decreases.

Every accepted move strictly decreases the count => termination + output is
never worse than input (dominance is structural).  Started from FFD it
inherits Dosa's tight bound FFD <= 11/9 OPT + 6/9.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def _bfd_assign_sizes(sizes, cap):
    """BFD over a size array -> (nbins, bin_id_per_input_position)."""
    m = sizes.shape[0]
    order = np.argsort(-sizes.astype(np.float64))
    resids = np.empty(m, dtype=np.int64)
    out = np.empty(m, dtype=np.int64)
    nb = 0
    for t in range(m):
        i = order[t]
        s = sizes[i]
        best_b = -1
        best_r = cap + 1
        for b in range(nb):
            r = resids[b]
            if r >= s and r < best_r:
                best_r = r
                best_b = b
        if best_b < 0:
            resids[nb] = cap - s
            nb += 1
            out[i] = nb - 1
        else:
            resids[best_b] -= s
            out[i] = best_b
    return nb, out


@njit(cache=True)
def ffd_assign(sizes, cap):
    """FFD -> (nbins, assign) for items in given order."""
    n = sizes.shape[0]
    sd = np.sort(sizes)[::-1].copy()
    resids = np.empty(n, dtype=np.int64)
    assign = np.empty(n, dtype=np.int64)
    nb = 0
    for i in range(n):
        s = sd[i]
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
        assign[i] = bb
    return nb, assign


@njit(cache=True)
def improve_ffd(sizes, cap, max_rounds=300, r_max=6, patience=40):
    """FFD start + iterated unload-r-smallest-bins repack.

    Returns (nbins_final, nbins_start); guarantees nbins_final <= nbins_start.
    """
    n = sizes.shape[0]
    sd = np.sort(sizes)[::-1].copy()
    cur, assign = ffd_assign(sizes, cap)
    nb_start = cur

    loads = np.zeros(2 * n + 2, dtype=np.int64)
    fails = 0
    rounds = 0

    while rounds < max_rounds and fails < patience and cur > 1:
        rounds += 1
        for b in range(cur):
            loads[b] = 0
        for i in range(n):
            loads[assign[i]] += sd[i]

        improved = False
        for r in range(1, min(r_max, cur - 1) + 1):
            # r smallest-load bins (ties -> lowest id)
            sel = np.full(r, -1, dtype=np.int64)
            for k in range(r):
                best_b = -1
                best_l = np.int64(1 << 60)
                for b in range(cur):
                    taken = False
                    for t in range(k):
                        if sel[t] == b:
                            taken = True
                            break
                    if not taken and loads[b] < best_l:
                        best_l = loads[b]
                        best_b = b
                sel[k] = best_b

            in_sel = np.zeros(cur, dtype=np.bool_)
            for k in range(r):
                in_sel[sel[k]] = True

            pool_idx = np.empty(n, dtype=np.int64)
            npool = 0
            for i in range(n):
                if in_sel[assign[i]]:
                    pool_idx[npool] = i
                    npool += 1
            pool_idx = pool_idx[:npool]

            nk = cur - r
            keep_ids = np.empty(nk, dtype=np.int64)
            keep_res = np.empty(nk, dtype=np.int64)
            j = 0
            for b in range(cur):
                if not in_sel[b]:
                    keep_ids[j] = b
                    keep_res[j] = cap - loads[b]
                    j += 1

            porder = np.argsort(-sd[pool_idx].astype(np.float64))
            trial_res = keep_res.copy()
            target_slot = np.empty(npool, dtype=np.int64)
            leftover_pos = np.empty(npool, dtype=np.int64)  # positions in pool
            nl = 0
            for t in range(npool):
                ii = pool_idx[porder[t]]
                s = sd[ii]
                placed = -1
                for b in range(nk):
                    if trial_res[b] >= s:
                        trial_res[b] -= s
                        placed = b
                        break
                target_slot[t] = placed
                if placed < 0:
                    leftover_pos[nl] = t
                    nl += 1

            extra = 0
            lbin = None
            if nl > 0:
                lsizes = np.empty(nl, dtype=np.int64)
                for z in range(nl):
                    lsizes[z] = sd[pool_idx[porder[leftover_pos[z]]]]
                extra, lbin = _bfd_assign_sizes(lsizes, cap)

            cand = nk + extra
            if cand < cur:
                for t in range(npool):
                    ii = pool_idx[porder[t]]
                    if target_slot[t] >= 0:
                        assign[ii] = keep_ids[target_slot[t]]
                for z in range(nl):
                    ii = pool_idx[porder[leftover_pos[z]]]
                    assign[ii] = cur - r + lbin[z]
                cur = cand
                improved = True
                break
        if improved:
            fails = 0
        else:
            fails += 1

    return cur, nb_start


@njit(cache=True)
def improve_from_bfd(sizes, cap, max_rounds=300, r_max=6, patience=40):
    """Same repair started from BFD. Returns (final, start)."""
    n = sizes.shape[0]
    sd = np.sort(sizes)[::-1].copy()
    cur, assign = _bfd_assign_sizes(sd, cap)
    nb_start = cur
    # reuse improve loop by running same logic inline (duplicate for speed)
    loads = np.zeros(2 * n + 2, dtype=np.int64)
    fails = 0
    rounds = 0
    while rounds < max_rounds and fails < patience and cur > 1:
        rounds += 1
        for b in range(cur):
            loads[b] = 0
        for i in range(n):
            loads[assign[i]] += sd[i]
        improved = False
        for r in range(1, min(r_max, cur - 1) + 1):
            sel = np.full(r, -1, dtype=np.int64)
            for k in range(r):
                best_b = -1
                best_l = np.int64(1 << 60)
                for b in range(cur):
                    taken = False
                    for t in range(k):
                        if sel[t] == b:
                            taken = True
                            break
                    if not taken and loads[b] < best_l:
                        best_l = loads[b]
                        best_b = b
                sel[k] = best_b
            in_sel = np.zeros(cur, dtype=np.bool_)
            for k in range(r):
                in_sel[sel[k]] = True
            pool_idx = np.empty(n, dtype=np.int64)
            npool = 0
            for i in range(n):
                if in_sel[assign[i]]:
                    pool_idx[npool] = i
                    npool += 1
            pool_idx = pool_idx[:npool]
            nk = cur - r
            keep_ids = np.empty(nk, dtype=np.int64)
            keep_res = np.empty(nk, dtype=np.int64)
            j = 0
            for b in range(cur):
                if not in_sel[b]:
                    keep_ids[j] = b
                    keep_res[j] = cap - loads[b]
                    j += 1
            porder = np.argsort(-sd[pool_idx].astype(np.float64))
            trial_res = keep_res.copy()
            target_slot = np.empty(npool, dtype=np.int64)
            leftover_pos = np.empty(npool, dtype=np.int64)
            nl = 0
            for t in range(npool):
                ii = pool_idx[porder[t]]
                s = sd[ii]
                placed = -1
                for b in range(nk):
                    if trial_res[b] >= s:
                        trial_res[b] -= s
                        placed = b
                        break
                target_slot[t] = placed
                if placed < 0:
                    leftover_pos[nl] = t
                    nl += 1
            extra = 0
            lbin = None
            if nl > 0:
                lsizes = np.empty(nl, dtype=np.int64)
                for z in range(nl):
                    lsizes[z] = sd[pool_idx[porder[leftover_pos[z]]]]
                extra, lbin = _bfd_assign_sizes(lsizes, cap)
            cand = nk + extra
            if cand < cur:
                for t in range(npool):
                    ii = pool_idx[porder[t]]
                    if target_slot[t] >= 0:
                        assign[ii] = keep_ids[target_slot[t]]
                for z in range(nl):
                    ii = pool_idx[porder[leftover_pos[z]]]
                    assign[ii] = cur - r + lbin[z]
                cur = cand
                improved = True
                break
        if improved:
            fails = 0
        else:
            fails += 1
    return cur, nb_start


def repair_best(sizes, cap):
    """min(improve_ffd, improve_from_bfd) counts (still >= each start)."""
    s = np.asarray(sizes, dtype=np.int32)
    a = improve_ffd(s, int(cap))[0]
    b = improve_from_bfd(s, int(cap))[0]
    return min(a, b)
