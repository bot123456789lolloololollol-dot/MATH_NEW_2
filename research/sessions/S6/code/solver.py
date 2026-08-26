"""
Exact solvers for monochromatic-solution-free colorings.

Problems:
  A) Interval generalized Schur: does a t-coloring of {1..N} avoid monochromatic
     solutions to x + y = k*z (x,y >= 1)?
  B) Cyclic generalized Schur: does a t-coloring of Z_n avoid monochromatic
     solutions to x + y = k*z (mod n)?

Method: deterministic backtracking over variables in a fixed order, with
symmetry breaking over color labels (first occurrences must appear in label
order -- any valid coloring can be relabeled to satisfy this, so this prunes
nothing but relabelings) and incremental constraint checking.

Soundness/completeness of "unsat" answers: the search enumerates all
assignments modulo color relabeling; if it fails, no valid coloring exists.
Deterministic; no randomness; no heuristic pruning beyond symmetry and
conflict detection.

CLI prints CSV rows.
"""

import sys
from sys import setrecursionlimit

setrecursionlimit(100000)


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def egcd(a, b):
    if b == 0:
        return a, 1, 0
    g, p, q = egcd(b, a % b)
    return g, q, p - (a // b) * q


def modinv(a, m):
    g, x, _ = egcd(a % m, m)
    assert g == 1, "not invertible"
    return x % m


# ---------------------------------------------------------------- constraints

def interval_constraints(N, k):
    """Triples (x, y, z) as 0-based variable indices of values in {1..N},
    original x <= y, with x + y = k*z."""
    out = []
    for z in range(1, N + 1):
        s = k * z
        for x in range(1, s // 2 + 1):
            y = s - x
            if y <= N:
                out.append((x - 1, y - 1, z - 1))
    return out


def cyclic_constraint_map(n, k):
    """For each unordered pair-rep (x, y) with x <= y: list of z with k*z == x+y.
    Returns dict[(x,y)] -> list[z]."""
    cons = {}
    for x in range(n):
        for y in range(x, n):
            s = (x + y) % n
            g = gcd(k, n)
            zs = []
            if s % g == 0:
                ng = n // g
                step = ng
                z0 = (modinv(k // g, ng) * ((s // g) % ng)) % ng
                zs = [(z0 + i * step) % n for i in range(g)]
            cons[(x, y)] = zs
    return cons


def cyclic_constraints(n, k):
    out = []
    cm = cyclic_constraint_map(n, k)
    for (x, y), zs in cm.items():
        for z in zs:
            out.append((x, y, z))
    return out


# ---------------------------------------------------------------- solver core

class CSP:
    """Backtracking solver.  Variables 0..n-1, colors 0..t-1.
    constraints: list of tuples (variable indices); forbidden: all equal."""

    def __init__(self, nvars, t, constraints, var_order=None):
        self.n = nvars
        self.t = t
        self.order = var_order or list(range(nvars))
        # per-variable constraint incidence
        self.inc = [[] for _ in range(nvars)]
        for ci, c in enumerate(constraints):
            for v in set(c):
                self.inc[v].append(ci)
        self.constraints = [tuple(sorted(set(c))) for c in constraints]
        self.color = [None] * nvars

    def solve(self):
        """Return a valid coloring list or None."""
        self.color = [None] * self.n
        return self._dfs(0, -1)

    def _dfs(self, oi, maxcolor):
        if oi == len(self.order):
            return self.color[:]
        v = self.order[oi]
        hi = min(maxcolor + 1, self.t - 1)
        for c in range(hi + 1):
            self.color[v] = c
            if self._consistent(v):
                r = self._dfs(oi + 1, max(maxcolor, c))
                if r is not None:
                    return r
            self.color[v] = None
        return None

    def _consistent(self, v):
        col = self.color
        cv = col[v]
        for ci in self.inc[v]:
            c = self.constraints[ci]
            ok = True
            same = True
            for u in c:
                cu = col[u]
                if cu is None:
                    same = False
                    break
                if cu != cv:
                    same = False
                    break
            if same:
                return False
        return True


def interval_coloring(N, k, t):
    """Return coloring of {1..N} (list, index i = value i+1) avoiding mono
    x+y=kz, or None."""
    cons = interval_constraints(N, k)
    return CSP(N, t, cons).solve()


def cyclic_coloring(n, k, t):
    """Return coloring of Z_n (list index = residue 0..n-1) avoiding mono
    x+y=kz mod n, or None."""
    cons = cyclic_constraints(n, k)
    return CSP(n, t, cons).solve()


# ---------------------------------------------------------------- conventions

def filter_constraint(x, y, z, convention):
    """Return True if this solution triple should be FORBIDDEN as a
    monochromatic solution under the given convention."""
    if convention == "all":
        return True
    if convention == "nontrivial":
        # exclude the all-equal diagonal solution x=y=z
        return not (x == y == z)
    if convention == "weak":
        # only pairwise-distinct solutions count
        return len({x, y, z}) == 3
    raise ValueError(convention)


def cyclic_coloring_conv(n, k, t, convention="all"):
    cons = [(x, y, z) for (x, y, z) in cyclic_constraints(n, k)
            if filter_constraint(x, y, z, convention)]
    return CSP(n, t, cons).solve()


# ---------------------------------------------------------------- Rado numbers

def rado_number(k, t, cap=20000):
    """Smallest N such that every t-coloring of {1..N} has mono x+y=kz.
    Assumes such N <= cap.  Returns N."""
    N = 1
    while interval_coloring(N, k, t) is not None:
        N += 1
        if N > cap:
            raise RuntimeError("cap exceeded")
    return N


# ---------------------------------------------------------------- CLI / data

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "validate"

    if mode == "validate":
        print("# validation against known values")
        print("RR_2(x+y=z), expect 5:", rado_number(1, 2))
        print("RR_3(x+y=z), expect 14:", rado_number(1, 3))

    elif mode == "cyclic":
        # args: kmax nmax t [convention]
        kmax = int(sys.argv[2]); nmax = int(sys.argv[3]); t = int(sys.argv[4])
        conv = sys.argv[5] if len(sys.argv) > 5 else "all"
        print(f"# convention={conv}, t={t}; row per k: n=2..{nmax}")
        print("k,colorable(n=2..%d)" % nmax)
        for k in range(2, kmax + 1):
            row = "".join(str(int(cyclic_coloring_conv(n, k, t, conv) is not None))
                          for n in range(2, nmax + 1))
            print(f"{k},{row}")

    elif mode == "interval":
        kmax = int(sys.argv[2]); t = int(sys.argv[3])
        for k in range(2, kmax + 1):
            try:
                print(f"RR_{t}(x+y={k}z) =", rado_number(k, t))
            except RuntimeError:
                print(f"RR_{t}(x+y={k}z) = >cap")
