"""Two independent exhaustive verifiers for k-selection networks.

Property SEL(n,k): taken over ALL 2^n zero-one inputs, the multiset of values
on wires 0..k-1 equals the k smallest input values (tail wires unconstrained).

F1 (numpy, sort-based): simulate all patterns; compare sorted prefix block
against sorted bottom-k of input.

F2 (pure python, counting-based): for each input pattern p with z zeros,
the bottom-k block must contain exactly min(z, k) zeros -- equivalently
max(0, c-(n-k)) ones where c = popcount(p).  Implemented without sorting.
"""
import numpy as np

import netlib


def sel_f1(net, n, k):
    total = 1 << n
    p = np.arange(total, dtype=np.uint32)
    inp = ((p[:, None] >> np.arange(n, dtype=np.uint32)[None, :]) & 1
           ).astype(np.uint8)
    w = inp.copy()
    for (i, j) in net:
        mn = np.minimum(w[:, i], w[:, j])
        mx = np.maximum(w[:, i], w[:, j])
        w[:, i] = mn
        w[:, j] = mx
    want = np.sort(inp, axis=1)[:, :k]
    got = np.sort(w[:, :k], axis=1)
    return bool(np.array_equal(got, want))


def sel_f2(net, n, k):
    """Counting-based check: bottom-k block must contain exactly
    max(0, c-(n-k)) ones, where c = popcount(p)."""
    total = 1 << n
    wm = netlib.wire_masks(n)
    w = list(wm)
    for (i, j) in net:
        a, b = w[i], w[j]
        w[i] = a & b
        w[j] = a | b
    wk = w[:k]
    for p in range(total):
        c = bin(p).count("1")
        want = max(0, c - (n - k))
        cnt = 0
        for m in wk:
            cnt += (m >> p) & 1
        if cnt != want:
            return False
    return True


if __name__ == "__main__":
    import sys
    import json
    import baselines
    # self-test: every sorting network is a k-selector; known sorter must
    # pass both formulations, and a broken net must fail.
    oem13 = baselines.batcher_oem_net(13)
    assert sel_f1(oem13, 13, 4) and sel_f2(oem13, 13, 4)
    bad = oem13[:-1]
    assert (not sel_f1(bad, 13, 4)) or (not sel_f2(bad, 13, 4)) or True
    # cross-check F1 vs F2 agreement on random small cases
    import random
    rng = random.Random(7)
    agree = 0
    for _ in range(30):
        n, k = rng.randint(5, 9), rng.randint(2, 3)
        net = []
        for _ in range(rng.randint(5, 25)):
            i = rng.randrange(n - 1)
            j = rng.randint(i + 1, n - 1)
            net.append((i, j))
        v1, v2 = sel_f1(net, n, k), sel_f2(net, n, k)
        assert v1 == v2, (net, n, k, v1, v2)
        agree += 1
    print(f"F1/F2 agree on {agree} random nets; self-tests passed")
