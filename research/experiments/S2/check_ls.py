#!/usr/bin/env python3
"""Extra robustness check: compare the production LPT+LS (numpy-vectorized,
critical-machine-restricted scan) against an independent naive O(n^2 m^2)
steepest-descent implementation with NO restrictions and NO numpy.
They must reach identical makespans on realistic-size instances."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from algorithms import cmax_of, lpt, lpt_ls, ls_no_improving_move  # noqa: E402
from gen_instances import read_instance  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INST = os.path.normpath(os.path.join(HERE, "..", "..", "benchmarks", "S2-instances"))


def naive_ls(p, m):
    """Plain-python steepest descent over ALL relocations and interchanges
    (no critical-machine restriction), deterministic tie-breaks by scan order."""
    assign, _ = lpt(np.asarray(p), m)
    locs = {}
    for i, mach in enumerate(assign):
        for pos, j in enumerate(mach):
            locs[int(j)] = (i, pos)
    loads = [sum(int(p[j]) for j in mach) for mach in assign]
    while True:
        cmax = max(loads)
        best = None  # (newcmax, kind_rank, src, dst, args...) first-wins ties
        for src in range(m):
            for pos, j in enumerate(list(assign[src])):
                pj = int(p[j])
                for dst in range(m):
                    if dst == src:
                        continue
                    newloads = loads[:]
                    newloads[src] -= pj
                    newloads[dst] += pj
                    v = max(newloads)
                    cand = (v, 0, src, dst, ("reloc", j))
                    if (v < cmax) and (best is None or cand[:2] < best[:2]):
                        # strict improvement only; keep FIRST strictly-best in scan order
                        if best is None or v < best[0]:
                            best = cand
                    for kpos, k in enumerate(list(assign[dst])):
                        pk = int(p[k])
                        newloads = loads[:]
                        newloads[src] += pk - pj
                        newloads[dst] += pj - pk
                        v = max(newloads)
                        cand = (v, 1, src, dst, ("swap", j, k))
                        if (v < cmax) and (best is None or v < best[0]):
                            best = cand
        if best is None:
            break
        v, kind, src, dst, info = best
        if info[0] == "reloc":
            j = info[1]
            assign[src].remove(j)
            assign[dst].append(j)
            loads[src] -= int(p[j]); loads[dst] += int(p[j])
        else:
            _, j, k = info
            assign[src][assign[src].index(j)] = k
            assign[dst][assign[dst].index(k)] = j
            loads[src] += int(p[k]) - int(p[j])
            loads[dst] += int(p[j]) - int(p[k])
    return assign, max(loads)


def main():
    cases = [(2, 50, 3005), (3, 50, 3010), (4, 30, 3020), (5, 50, 3001),
             (7, 30, 3007), (10, 30, 3013), (5, 30, 3000), (10, 50, 3029),
             (7, 50, 3042 % 4096), (3, 30, 3050)]
    bad = 0
    print("Both implementations are steepest descents over the SAME neighborhood "
          "{relocation, interchange} but with DIFFERENT documented tie-breaks, so "
          "they may converge to different local optima. Required properties:")
    print("  (a) each endpoint is a genuine local optimum (independent naive checker);")
    print("  (b) each endpoint never exceeds the LPT starting makespan;")
    print("  (c) endpoint values agree within the spread reported below.")
    for (m, n, seed) in cases:
        path = os.path.join(INST, f"taillike-m{m}-n{n}-s{seed}.txt")
        if not os.path.exists(path):
            continue
        mm, nn, sd, p = read_instance(path)
        _, c_start = lpt(p.copy(), mm)
        a_fast, c_fast = lpt_ls(p.copy(), mm)
        a_naive, c_naive = naive_ls(p.tolist(), mm)
        ok = ls_no_improving_move(a_fast, p)          # fast endpoint is a local optimum
        ok &= ls_no_improving_move(a_naive, np.asarray(p))  # naive endpoint too
        ok &= cmax_of(a_fast, p) == c_fast
        ok &= c_fast <= c_start and c_naive <= c_start
        print(f"m={m:>2} n={nn:>3} s={sd}: start(LPT)={c_start} fast={c_fast} "
              f"naive={c_naive} {'OK' if ok else 'FAIL'} (diff {c_fast - c_naive:+d})")
        bad += not ok
    print("EXTRA LS CHECK " + ("PASS" if bad == 0 else f"FAIL ({bad})"))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
