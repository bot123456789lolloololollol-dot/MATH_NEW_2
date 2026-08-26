"""Additional classical baselines: DJD / DJD-threshold family (jitted).

DJD (Djang & Finch, 1998): process items largest-first; open a bin, place the
largest remaining item, then keep adding the largest item that fits until the
bin load reaches a threshold (here: a fraction of capacity or "all it can"),
then close the bin.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def djd(sizes, cap, mode):
    """mode 0: fill until nothing fits (greedy full bin)
       mode k>0: close bin once load >= k*cap/4 (k=1..4)"""
    n = sizes.shape[0]
    taken = np.zeros(n, dtype=np.bool_)
    order = np.argsort(-sizes.astype(np.float64))
    nbins = 0
    remaining = n
    while remaining > 0:
        nbins += 1
        load = 0
        # seed with largest untaken
        for t in range(n):
            i = order[t]
            if not taken[i]:
                taken[i] = True
                load += sizes[i]
                remaining -= 1
                break
        # greedily add largest fitting until threshold / nothing fits
        progress = True
        while progress:
            progress = False
            for t in range(n):
                i = order[t]
                s = sizes[i]
                if not taken[i] and load + s <= cap:
                    if mode > 0 and load >= mode * cap / 4.0:
                        break
                    taken[i] = True
                    load += s
                    remaining -= 1
                    progress = True
                    break
            if mode == 0:
                continue
    return nbins


def djd_full(sizes, cap):
    return int(djd(np.asarray(sizes, dtype=np.int32), int(cap), 0))


def djd_half(sizes, cap):
    return int(djd(np.asarray(sizes, dtype=np.int32), int(cap), 2))
