"""Evaluation harness: fitness of programs on datasets, parallel runners."""
import os
import numpy as np
from concurrent.futures import ProcessPoolExecutor

from .core import ffd, bfd, wfd, l1_lower_bound, l2_lower_bound
from .gp import pack_gp

# ------------------------------------------------------------------ datasets


def prep_dataset(dataset):
    """Pre-sort instances desc once; return list of (sizes_desc, cap)."""
    out = []
    for sizes, cap in dataset:
        s = np.sort(np.asarray(sizes, dtype=np.int32))[::-1].copy()
        out.append((s, int(cap)))
    return out


def eval_program_on_dataset(ops, consts, dataset):
    """Return (mean_excess_ratio, total_bins). excess = (bins - L2)/L2."""
    tot = 0.0
    bins_sum = 0
    for s, cap in dataset:
        nb, _ = pack_gp(s, cap, ops, consts)
        lb = l2_lower_bound(s, cap)
        tot += (nb - lb) / lb
        bins_sum += nb
    return tot / len(dataset), bins_sum


def eval_baselines_on_instance(sizes, cap):
    s = np.asarray(sizes, dtype=np.int32)
    sd = np.sort(s)[::-1].copy()
    res = {
        "ffd": ffd(sd, cap),
        "bfd": bfd(sd, cap),
        "wfd": wfd(sd, cap),
        "l1": int(l1_lower_bound(s, cap)),
        "l2": int(l2_lower_bound(s, cap)),
    }
    return res


def prep_dataset_paired(dataset):
    """Pre-sort instances desc and precompute paired baseline info.

    Returns list of (sizes_desc, cap, base_bins, lb) where
    base_bins = min(FFD, BFD) on that instance and lb = L2 lower bound.
    """
    out = []
    for sizes, cap in dataset:
        s = np.asarray(sizes, dtype=np.int32)
        sd = np.sort(s)[::-1].copy()
        f = ffd(sd, cap)
        b = bfd(sd, cap)
        lb = max(1, int(l2_lower_bound(s, cap)))
        out.append((s, int(cap), min(f, b), lb))
    return out


def eval_paired_improvement(ops, consts, paired):
    """Mean over instances of (H - base)/lb.  0 == ties best classical;
    negative == strictly better; positive == worse."""
    tot = 0.0
    for sd, cap, base, lb in paired:
        nb, _ = pack_gp(sd, cap, ops, consts)
        tot += (nb - base) / lb
    return tot / len(paired)


# ------------------------------------------------------------------ parallel eval
_WORKER_DATA = {}


def _init_worker(train_pkl):
    _WORKER_DATA["ds"] = train_pkl
    _WORKER_DATA["paired"] = None


def _eval_one(payload):
    ops, consts, idx = payload
    ds = _WORKER_DATA["ds"]
    sub = ds if idx is None else ds[idx]
    v, _ = eval_program_on_dataset(ops, consts, sub)
    return v


def _eval_paired_one(payload):
    ops, consts, idx = payload
    paired = _WORKER_DATA["paired"]
    if idx is not None:
        paired = [paired[i] for i in idx]
    return eval_paired_improvement(ops, consts, paired)


def _init_worker_paired(paired_pkl):
    _WORKER_DATA["paired"] = paired_pkl


class ParallelEvaluator:
    """Evaluates populations over a fixed dataset using persistent processes.

    mode="excess":  fitness = mean((H-L2)/L2)
    mode="paired":  fitness = mean((H-min(FFD,BFD))/L2)  -- 0 for ties
    """

    def __init__(self, dataset, n_workers=None, mode="paired"):
        self.dataset = prep_dataset(dataset)
        self.paired = prep_dataset_paired(dataset) if mode == "paired" else None
        self.mode = mode
        self.n_workers = n_workers or max(1, (os.cpu_count() or 4) - 1)
        self._pool = None

    def _get_pool(self):
        if self._pool is None:
            if self.mode == "paired":
                self._pool = ProcessPoolExecutor(
                    max_workers=self.n_workers,
                    initializer=_init_worker_paired,
                    initargs=(self.paired,),
                )
            else:
                self._pool = ProcessPoolExecutor(
                    max_workers=self.n_workers,
                    initializer=_init_worker,
                    initargs=(self.dataset,),
                )
        return self._pool

    def evaluate(self, progs, subset_idx=None):
        payloads = [(p.ops, p.consts, subset_idx) for p in progs]
        if self.n_workers <= 1:
            fn = _eval_one if self.mode != "paired" else _eval_paired_one
            return [fn(x) for x in payloads]
        fn = _eval_one if self.mode != "paired" else _eval_paired_one
        return list(self._get_pool().map(fn, payloads, chunksize=1))

    def close(self):
        if self._pool is not None:
            self._pool.shutdown()
            self._pool = None
