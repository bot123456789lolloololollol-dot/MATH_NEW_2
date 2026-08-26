"""Baseline sizes for k-selection networks.

For each preregistered (N, k):
  B1: smallest known SORTER on N wires (any sorter is a k-selector)
  B2: B1 pruned under the selector property (drop unneeded CEs)
  B3: SorterHunter median net for N when its k matches
All values exhaustively verified with both selector formulations.
"""
import glob
import json

import numpy as np  # noqa: F401  (verifier dependency)

import netlib
import verify_selector as vs


def best_sorter(n):
    cands = []
    for f in glob.glob(f"seeds/Sorters_Sort_{n}_*.json"):
        with open(f) as fh:
            d = json.load(fh)
        cands.append((d["L"], [tuple(e) for e in d["nw"]]))
    import baselines
    cands.append((len(basels := baselines.batcher_oem_net(n)), basels))
    return min(cands, key=lambda x: x[0])


def main():
    targets = [(10, 4), (10, 5), (11, 5), (12, 5), (12, 6), (13, 6),
               (14, 6), (14, 7), (15, 7), (16, 7), (16, 8)]
    out = {}
    print(f"{'N':>3} {'k':>2} {'B1(sorter)':>11} {'B2(pruned)':>11} "
          f"{'B3(median)':>10} {'B1depth':>8} {'B2depth':>8}")
    for n, k in targets:
        L1, s1 = best_sorter(n)
        assert vs.sel_f1(s1, n, k) and vs.sel_f2(s1, n, k)
        # greedy prune under selector judge
        s2 = list(s1)
        changed = True
        while changed:
            changed = False
            r = 0
            while r < len(s2):
                cand = s2[:r] + s2[r + 1:]
                if vs.sel_f1(cand, n, k):
                    s2 = cand
                    changed = True
                else:
                    r += 1
        assert vs.sel_f1(s2, n, k) and vs.sel_f2(s2, n, k)
        b3 = "-"
        if k == (n - 1) // 2:
            fs = glob.glob(f"seeds/Median_Median_{n}_*.json")
            if fs:
                fb = min(fs, key=lambda f: json.load(open(f))["L"])
                with open(fb) as fh:
                    d = json.load(fh)
                m = [tuple(e) for e in d["nw"]]
                assert vs.sel_f1(m, n, k) and vs.sel_f2(m, n, k)
                b3 = len(m)
        out[f"{n}_{k}"] = {"B1": L1, "B2": len(s2),
                           "B2_depth": netlib.depth(s2), "B3": b3}
        print(f"{n:>3} {k:>2} {L1:>11} {len(s2):>11} {str(b3):>10} "
              f"{netlib.depth(s1):>8} {netlib.depth(s2):>8}")
    with open("select_baselines.json", "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
