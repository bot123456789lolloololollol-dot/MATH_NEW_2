"""Classical constructive sorting-network generators (baselines).

All outputs are validated exhaustively by validate_baselines.py before use;
generator bugs therefore cannot contaminate any result.
"""
import netlib


def insertion_net(n):
    net = []
    for i in range(1, n):
        for j in range(i - 1, -1, -1):
            net.append((j, j + 1))
    return net


def batcher_oem_net(n):
    """Batcher odd-even mergesort via the general two-run merge, valid for
    every n (not only powers of two)."""
    net = []

    def merge2(a, b):
        """Merge two sorted runs of wire ids."""
        if not a or not b:
            return
        if len(a) == 1 and len(b) == 1:
            net.append((min(a[0], b[0]), max(a[0], b[0])))
            return
        merge2(a[0::2], b[0::2])
        merge2(a[1::2], b[1::2])
        comb = a + b
        for k in range(1, len(comb) - 1, 2):
            x, y = comb[k], comb[k + 1]
            net.append((min(x, y), max(x, y)))

    def sort(ws):
        if len(ws) <= 1:
            return
        m = len(ws) // 2
        sort(ws[:m])
        sort(ws[m:])
        merge2(ws[:m], ws[m:])

    sort(list(range(n)))
    return net


def bose_nelson_net(n):
    """UNVALIDATED reconstruction of the Bose-Nelson recursion -- my
    transcription fails exhaustive validation at n=4, so this generator is
    NOT USED anywhere in the pipeline (kept only as a negative result)."""
    net = []

    def M(i1, m1, i2, m2):
        if m1 == 1 and m2 == 1:
            net.append((min(i1, i2), max(i1, i2)))
        elif m1 == 1 and m2 == 2:
            net.append((min(i1, i2 + 1), max(i1, i2 + 1)))
            net.append((min(i1, i2), max(i1, i2)))
        elif m1 == 2 and m2 == 1:
            net.append((min(i1 + 1, i2), max(i1 + 1, i2)))
            net.append((min(i1, i2), max(i1, i2)))
        else:
            M(i1, m1 // 2, i2, m2 // 2)
            M(i1 + m1 // 2, m1 - m1 // 2, i2 + m2 // 2, m2 - m2 // 2)
            M(i1 + m1 // 2, m1 - m1 // 2, i2 + m2 // 2, m2 // 2)

    def P(i, m):
        if m >= 1:
            m1 = m // 2
            m2 = m - m1
            if m1 > 0:
                P(i, m1)
                P(i + m1, m2)
                M(i, m1, i + m1, m2)

    P(0, n)
    return net


# Small networks transcribed from Knuth, TAOCP vol. 3, 5.3.4 (verified
# exhaustively below before any use).
KNUTH_SMALL = {
    2: [(0, 1)],
    3: [(0, 1), (1, 2), (0, 1)],
    4: [(0, 1), (2, 3), (0, 2), (1, 3), (1, 2)],
    5: [(0, 1), (3, 4), (2, 4), (2, 3), (1, 4), (0, 3), (0, 2), (1, 3), (1, 2)],
    6: [(1, 2), (4, 5), (0, 2), (3, 5), (0, 1), (3, 4),
        (2, 5), (0, 3), (1, 4), (2, 4), (1, 3), (2, 3)],
}

if __name__ == "__main__":
    print(f"{'n':>3} {'ins':>5} {'OEM':>5} {'dOEM':>5} "
          f"{'OEMpruned':>10} {'dOEMp':>6} {'opt?':>6}")
    for n in range(2, 33):
        ins = insertion_net(n)
        oem = batcher_oem_net(n)
        assert netlib.is_sorting(ins, n), f"insertion broken at n={n}"
        assert netlib.is_sorting(oem, n), f"batcher broken at n={n}"
        op = netlib.prune(oem, n)
        if n in KNUTH_SMALL:
            assert netlib.is_sorting(KNUTH_SMALL[n], n), \
                f"knuth transcription broken at n={n}"
        print(f"{n:>3} {len(ins):>5} {len(oem):>5} "
              f"{netlib.depth(oem):>5} {len(op):>10} {netlib.depth(op):>6}"
              + (f" {'OK':>6}" if n in KNUTH_SMALL else ""))
