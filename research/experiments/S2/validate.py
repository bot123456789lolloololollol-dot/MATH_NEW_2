#!/usr/bin/env python3
"""Validation gates for S2 (must pass before any evaluation run).

1. LPT on Graham's hand instance p=[3,3,2,2,2], m=2 -> Cmax = 7 (OPT=6, LB=6).
   MULTIFIT on the same instance -> minimal feasible capacity 6.
2. All algorithms: assignments valid (each job exactly once), Cmax recomputable.
3. Brute-force optimal for tiny instances: every algorithm's Cmax >= OPT; LS ends
   at a local optimum verified by an independent naive checker.
4. Determinism: repeated runs give identical results.
Writes outputs/validation_report.txt
"""
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from algorithms import (lb_bound, lpt, lpt_ls, ls_no_improving_move, cmax_of,
                        multifit, rol)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)

lines = []
def log(s=""):
    print(s)
    lines.append(str(s))


def brute_opt(p, m):
    n = len(p)
    best = None
    for choice in itertools.product(range(m), repeat=n):
        loads = [0] * m
        for j, mi in enumerate(choice):
            loads[mi] += int(p[j])
        v = max(loads)
        if best is None or v < best:
            best = v
    return best


def check_valid(assign, p, m, tag):
    allj = sorted(j for mach in assign for j in mach)
    assert allj == list(range(len(p))), f"{tag}: assignment not a permutation: {allj[:10]}..."
    assert len(assign) == m, f"{tag}: wrong machine count"
    return True


rng = np.random.default_rng(7)
fails = 0

# ---- Gate 1: hand-checkable instance ---------------------------------------
p_hand = np.array([3, 3, 2, 2, 2], dtype=np.int64)
m_hand = 2
a, c_lpt = lpt(p_hand, m_hand)
log(f"[gate1] LPT p=[3,3,2,2,2] m=2 -> Cmax={c_lpt} expect 7 "
    f"{'PASS' if c_lpt == 7 else 'FAIL'}")
fails += c_lpt != 7
a_mf, c_mf = multifit(p_hand, m_hand)
log(f"[gate1] MULTIFIT same -> Cmax={c_mf} expect 6 "
    f"{'PASS' if c_mf == 6 else 'FAIL'}")
fails += c_mf != 6
a_ls, c_ls = lpt_ls(p_hand, m_hand)
log(f"[gate1] LPT+LS same -> Cmax={c_ls} (here LS reaches OPT={brute_opt(p_hand, m_hand)} "
    f"via one swap; LS remains heuristic in general - see gate23)")
check_valid(a_ls, p_hand, m_hand, "LS-hand")
ok_naive = ls_no_improving_move(a_ls, p_hand)
log(f"[gate1] LS local-optimum verified by independent naive checker: {ok_naive} "
    f"{'PASS' if ok_naive else 'FAIL'}")
fails += not ok_naive
lbh = lb_bound(p_hand, m_hand)
log(f"[gate1] LB(hand)={lbh} expect 6 {'PASS' if lbh == 6 else 'FAIL'}")
fails += lbh != 6

# ---- Gates 2+3: random tiny instances vs brute force ------------------------
n_bad = 0
for t in range(200):
    n = int(rng.integers(4, 8))
    mm = int(rng.integers(2, 4))
    pt = rng.integers(1, 20, size=n).astype(np.int64)
    opt = brute_opt(pt, mm)
    for name, fn in (("lpt", lambda x, y: lpt(x, y)),
                     ("multifit", lambda x, y: multifit(x, y)),
                     ("lpt_ls", lambda x, y: lpt_ls(x, y)),
                     ("rol", lambda x, y: rol(x, y, K_mult=2))):
        assign, c = fn(pt.copy(), mm)
        try:
            check_valid(assign, pt, mm, name)
            rc = cmax_of(assign, pt)
            assert rc == c, f"{name}: reported {c} != recomputed {rc}"
            assert c >= opt, f"{name}: Cmax {c} < OPT {opt}"
            if name == "lpt_ls":
                assert ls_no_improving_move(assign, pt), "LS not locally optimal"
        except AssertionError as e:
            n_bad += 1
            log(f"[gate23] FAIL t={t} {name}: {e}")
log(f"[gate23] 200 random tiny instances x 4 algorithms vs brute force: "
    f"{'ALL PASS' if n_bad == 0 else str(n_bad) + ' FAILURES'}")
fails += n_bad

# ---- Gate 4: determinism ----------------------------------------------------
pdet = np.random.default_rng(123).integers(1, 100, size=60).astype(np.int64)
det_ok = True
for name, fn in (("lpt", lambda x, y: lpt(x, y)),
                 ("multifit", lambda x, y: multifit(x, y)),
                 ("lpt_ls", lambda x, y: lpt_ls(x, y)),
                 ("rol", lambda x, y: rol(x, y, K_mult=2))):
    r1 = fn(pdet.copy(), 5)
    r2 = fn(pdet.copy(), 5)
    same = (r1[1] == r2[1]) and ([sorted(x) for x in r1[0]] == [sorted(x) for x in r2[0]])
    log(f"[gate4] determinism {name}: {'PASS' if same else 'FAIL'} (Cmax={r1[1]})")
    det_ok &= same
fails += not det_ok

log()
log("VALIDATION " + ("PASS" if fails == 0 else f"FAIL ({fails})"))
with open(os.path.join(OUT, "validation_report.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")
sys.exit(0 if fails == 0 else 1)
