"""
SAT-based exact solver for monochromatic-solution-free coloring problems.
Independent implementation cross-checking solver.py's backtracking search.

Encoding:
  variable x(i,c) = "value/residue i has color c"   (i 0-based internally)
  - at-least-one color per i; at-most-one via pairwise clauses
  - for every forbidden solution triple and color c:
        NOT(x(i,c) AND x(j,c) AND x(l,c))
  - symmetry breaking: residue/value 0 fixed to color 0 (safe by relabeling)

Solver: CaDiCaL (via python-sat). Deterministic given fixed instance.
"""

import sys
from pysat.solvers import Cadical153
from pysat.formula import IDPool

from solver import (interval_constraints, cyclic_constraints,
                    filter_constraint)


def sat_coloring(nvars, triples, t):
    """Return list color[i] or None if UNSAT."""
    pool = IDPool()
    cnf = []
    for i in range(nvars):
        cnf.append([pool.id(("x", i, c)) for c in range(t)])
        for c in range(t):
            for d in range(c + 1, t):
                cnf.append([-pool.id(("x", i, c)), -pool.id(("x", i, d))])
    # symmetry breaking: element 0 has color 0
    if nvars > 0:
        cnf.append([-pool.id(("x", 0, c)) for c in range(1, t)])
    for (a, b, c3) in triples:
        for col in range(t):
            cnf.append([-pool.id(("x", a, col)),
                        -pool.id(("x", b, col)),
                        -pool.id(("x", c3, col))])
    s = Cadical153(bootstrap_with=cnf)
    if s.solve():
        model = set(l for l in s.get_model() if l > 0)
        inv = {v: k for k, v in pool.obj2id.items()}
        color = [None] * nvars
        for l in model:
            key = inv[l]
            if key[0] == "x":
                color[key[1]] = key[2]
        return color
    return None


def interval_coloring_sat(N, k, t, conv="all"):
    tr = [(x, y, z) for (x, y, z) in interval_constraints(N, k)
          if filter_constraint(x + 1, y + 1, z + 1, conv)]
    return sat_coloring(N, tr, t)


def cyclic_coloring_sat(n, k, t, conv="all"):
    tr = [(x, y, z) for (x, y, z) in cyclic_constraints(n, k)
          if filter_constraint(x, y, z, conv)]
    return sat_coloring(n, tr, t)


def rado_number_sat(k, t, conv="all", cap=100000):
    N = 1
    while True:
        col = interval_coloring_sat(N, k, t, conv)
        if col is None:
            return N
        N += 1
        assert N <= cap


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "rado":
        kmax, t = int(sys.argv[2]), int(sys.argv[3])
        conv = sys.argv[4] if len(sys.argv) > 4 else "all"
        out = []
        for k in range(2, kmax + 1):
            r = rado_number_sat(k, t, conv)
            print(f"{k},{r}", flush=True)
    elif mode == "cyclic":
        k, t = int(sys.argv[2]), int(sys.argv[4] if len(sys.argv) > 4 else 2)
        conv = sys.argv[5] if len(sys.argv) > 5 else "weak"
        lo, hi = int(sys.argv[3]), int(sys.argv[3])
        row = []
        import json
        data = {}
        kmax = k
        hi = int(sys.argv[3])
        for kk in range(2, kmax + 1):
            S = [n for n in range(2, hi + 1)
                 if cyclic_coloring_sat(n, kk, t, conv) is not None]
            data[kk] = S
            print(kk, S, flush=True)
        json.dump(data, open(f"cyclic_{conv}_t{t}.json", "w"))
