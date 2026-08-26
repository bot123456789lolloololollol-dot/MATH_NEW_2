"""Sanity checks: kernels vs brute-force optimum; L1/L2 validity."""
import sys, time
import numpy as np

sys.path.insert(0, r"C:\Users\Deves\OneDrive\Desktop\MATH_2_NEW")
from bpp.core import (first_fit, best_fit, worst_fit, ffd, bfd, wfd,
                      l1_lower_bound, l2_lower_bound, optimal_bins_small)

rng = np.random.default_rng(0)

# warm up jit
ffd(np.array([3, 5, 5, 4], dtype=np.int32), 10)
bfd(np.array([3, 5, 5, 4], dtype=np.int32), 10)

fails = 0
t0 = time.time()
for trial in range(3000):
    n = int(rng.integers(1, 13))
    cap = int(rng.integers(5, 40))
    s = rng.integers(1, cap + 1, size=n).astype(np.int32)
    opt = optimal_bins_small(s, cap)
    lb1 = l1_lower_bound(s, cap)
    lb2 = l2_lower_bound(s, cap)
    ff = first_fit(s, cap)[0]
    bf = best_fit(s, cap)[0]
    wf = worst_fit(s, cap)[0]
    f = ffd(s, cap)
    b = bfd(s, cap)
    w = wfd(s, cap)
    ok = (lb1 <= opt and lb2 <= opt and lb1 <= lb2
          and opt <= min(ff, bf, f, b) and opt <= max(wf, w))
    if not ok:
        fails += 1
        print("FAIL", s, cap, "opt", opt, "l1", lb1, "l2", lb2,
              "ff", ff, "bf", bf, "wf", wf, "ffd", f, "bfd", b, "wfd", w)
        if fails > 5:
            break
print(f"3000 trials vs exact optimum in {time.time()-t0:.1f}s, failures={fails}")

# L2 tightness spot-check: known case where L1 < OPT == L2
s = np.array([6, 6, 6, 6, 6, 6, 10, 10, 10, 10, 10, 10], dtype=np.int32)
print("L1:", l1_lower_bound(s, 18), "L2:", l2_lower_bound(s, 18),
      "OPT:", optimal_bins_small(s, 18))

# speed probe: how many FFD evals/sec on n=200 instances
s = rng.integers(1, 1001, size=200).astype(np.int32)
t0 = time.time()
R = 20000
for _ in range(R):
    ffd(s, 1000)
dt = time.time() - t0
print(f"FFD(n=200): {R/dt:,.0f} evals/sec single-core")
