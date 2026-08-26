"""Independent exhaustive verifier using numpy simulation over all 2^n
input patterns -- deliberately implemented differently from netlib's bigint
bitset evaluator so that the two cannot share an implementation bug.

For sorting: applies compare-exchanges to every pattern simultaneously and
checks (a) output is non-decreasing along wires, (b) output is a permutation
of the input (multiset preservation, not assumed even though comparators
preserve multisets by construction).

For selection of rank r: checks some wire equals the r-th smallest value on
every pattern, plus optional full-multiset checks.
"""
import numpy as np


def simulate_all(net, n, chunk_bits=20):
    """Yield (inputs, outputs) chunks; inputs/outputs shape (chunk, n)
    uint8, covering ALL 2^n patterns exactly once."""
    total = 1 << n
    wires = np.arange(total, dtype=np.uint64)
    # inputs[p, i] = bit i of p
    start = 0
    csz = 1 << chunk_bits
    while start < total:
        end = min(start + csz, total)
        p = np.arange(start, end, dtype=np.uint64)
        inp = ((p[:, None] >> np.arange(n, dtype=np.uint64)[None, :]) & 1
               ).astype(np.uint8)
        w = inp.copy()
        for (i, j) in net:
            wi = w[:, i]
            wj = w[:, j]
            mn = np.minimum(wi, wj)
            mx = np.maximum(wi, wj)
            w[:, i] = mn
            w[:, j] = mx
        yield inp, w
        start = end


def check_sorting(net, n):
    """True iff network sorts ALL 2^n zero-one inputs (independently)."""
    for inp, out in simulate_all(net, n):
        # multiset preservation
        if not np.array_equal(np.sort(inp, axis=1), np.sort(out, axis=1)):
            return False
        # sortedness
        if not np.all(out[:, :-1] <= out[:, 1:]):
            return False
    return True


def rank_values(inp, r):
    """r-th smallest per row."""
    return np.sort(inp, axis=1)[:, r]


def check_median(net, n, wire=None):
    """Check whether some wire always carries the lower-median
    (rank floor((n-1)/2)) of the inputs. If wire given, check that one."""
    r = (n - 1) // 2
    for inp, out in simulate_all(net, n):
        rv = rank_values(inp, r)
        if wire is not None:
            if not np.array_equal(out[:, wire], rv):
                return False
        else:
            if not np.any(np.all(out == rv[:, None], axis=0)):
                return False
    return True


def check_half_selection(net, n, h=None):
    """Check whether wires 0..hh-1 hold the hh smallest values (unordered),
    i.e. a half-selection/median-group network."""
    if h is None:
        h = n // 2
    for inp, out in simulate_all(net, n):
        want = np.sort(inp, axis=1)[:, :h]
        got = np.sort(out[:, :h], axis=1)
        if not np.array_equal(want, got):
            return False
        rest_ok = np.sort(inp, axis=1)[:, h:]
        got2 = np.sort(out[:, h:], axis=1)
        if not np.array_equal(rest_ok, got2):
            return False
    return True


def load_seed(path):
    """Load a SorterHunter JSON: 'nw' is a flat list of [i, j] pairs."""
    import json
    with open(path) as f:
        d = json.load(f)
    nw = d["nw"]
    if not all(isinstance(e, list) and len(e) == 2 and
               all(isinstance(x, int) for x in e) for e in nw):
        raise ValueError(f"unrecognized nw structure in {path}")
    return d["N"], [tuple(e) for e in nw], d


if __name__ == "__main__":
    import sys
    import netlib
    path = sys.argv[1]
    n, nw, d = load_seed(path)
    print(f"N={n} L={len(nw)} D={netlib.depth(nw)} "
          f"(file says L={d['L']} D={d['D']})")
    print("bitset verifier:", netlib.is_sorting(nw, n))
    print("numpy verifier:", check_sorting(nw, n))
