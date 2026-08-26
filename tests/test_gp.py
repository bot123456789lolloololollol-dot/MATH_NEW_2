"""Smoke-test the GP interpreter: hand-coded BF/WF programs + random progs."""
import sys, time
import numpy as np

sys.path.insert(0, r"C:\Users\Deves\OneDrive\Desktop\MATH_2_NEW")
from bpp.core import ffd, bfd
from bpp.gp import (pack_gp, rand_prog, rand_consts, ProgSpec,
                    program_to_str, PUSH_FEAT, PUSH_CONST, ADD, SUB, MUL,
                    DIV, MIN, MAX, ABS, NEG, INV)

rng = np.random.default_rng(1)


def make_prog(instrs, consts=None):
    ops = np.array([x for pair in instrs for x in pair], dtype=np.int32).reshape(-1, 2)
    if consts is None:
        consts = rand_consts(rng)
    return ProgSpec(ops, consts)


# Best Fit program: minimize fit -> score = -fit  (virtual bin has fit=cap-s, large)
bf = make_prog([[PUSH_FEAT, 1], [NEG, 0]])
# Worst Fit: maximize fit
wf = make_prog([[PUSH_FEAT, 1]])
# FFD-like: prefer fullest bin that fits => score = load_frac (ties->lowest idx)
ff_like = make_prog([[PUSH_FEAT, 3]])


def sort_desc(s):
    return np.sort(np.asarray(s, dtype=np.int32))[::-1].copy()


fails = 0
for trial in range(500):
    n = int(rng.integers(2, 60))
    cap = int(rng.integers(20, 200))
    s = rng.integers(1, cap + 1, size=n).astype(np.int32)
    sd = sort_desc(s)
    # reference best fit on desc order:
    from bpp.core import best_fit
    bf_ref = best_fit(sd.astype(np.int32), cap)[0]
    bf_gp = pack_gp(sd, cap, bf.ops, bf.consts)[0]
    wf_ref_gp = pack_gp(sd, cap, wf.ops, wf.consts)[0]
    if bf_gp != bf_ref:
        fails += 1
        print("BF mismatch", trial, "gp", bf_gp, "ref", bf_ref)
        if fails > 4:
            break
print("BF-program matches kernel best_fit on 500 instances:", fails == 0)

# random programs run without crash and give sane counts
t0 = time.time()
nbins_rand = []
for _ in range(200):
    p = ProgSpec(rand_prog(rng), rand_consts(rng))
    n = int(rng.integers(10, 300))
    cap = int(rng.integers(100, 1000))
    s = rng.integers(1, cap + 1, size=n).astype(np.int32)
    sd = sort_desc(s)
    nb, _ = pack_gp(sd, cap, p.ops, p.consts)
    assert nb > 0 and nb <= len(s)
    nbins_rand.append(nb)
print(f"200 random programs OK ({time.time()-t0:.1f}s incl compile)")

# speed: one program over many instances
s = rng.integers(1, 1001, size=300).astype(np.int32)
sd = sort_desc(s)
p = make_prog([[PUSH_FEAT, 1], [PUSH_FEAT, 7], [MUL, 0], [PUSH_FEAT, 0], [DIV, 0],
               [PUSH_FEAT, 17], [SUB, 0], [PUSH_FEAT, 8], [MAX, 0]],
              rand_consts(rng))
pack_gp(sd[:50], 1000, p.ops, p.consts)  # warm
t0 = time.time()
R = 2000
for _ in range(R):
    pack_gp(sd, 1000, p.ops, p.consts)
dt = time.time() - t0
print(f"GP evals n=300: {R/dt:,.0f} evals/sec single-core")

print("example program:", program_to_str(p))
