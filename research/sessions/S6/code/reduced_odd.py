"""
Reduced exact solver for ODD k, weak convention, 2 colors.

Construction family ("class-halfplane + multiple-sequence"):
    c(v) = 1 if v mod k in [1, m]        (m = k//2, non-multiples)
    c(k*q) = u_q                          (multiples, free choice)

THEOREM-CANDIDATE (reduction): for odd k >= 5, a coloring of [N] avoiding
mono distinct solutions of x+y=kz exists IFF such a structured coloring
exists.  (Direction <= is the content; direction => is trivial.)

This module decides, for given odd k and N, whether some u_1..u_{N//k}
satisfies all reduced constraints:
    for 1 <= i < j <= N//k, skipping degenerate pairs (i+j == k*i or == k*j):
       not (u_i == u_j == M(i+j)),   M(s) = 1 iff s mod k in [1,m],
       and if s = i+j is itself a multiple of k, M(s) := u_{s/k}.
"""

from pysat.solvers import Cadical153
from pysat.formula import IDPool


def reduced_sat(k, N):
    """Odd k. True iff some u works for prefix [N]. Returns (sat, u or None)."""
    m = k // 2
    qmax = N // k
    if qmax == 0:
        return True, []
    pool = IDPool(); cnf = []
    U = {q: pool.id(("u", q)) for q in range(1, qmax + 1)}
    def M(s):
        r = s % k
        if r == 0:
            return U[s // k] if s // k <= qmax else None
        return True if 1 <= r <= m else False
    for i in range(1, qmax + 1):
        for j in range(i + 1, qmax + 1):
            if j == (k - 1) * i or i == (k - 1) * j:
                continue
            ms = M(i + j)
            if ms is True:
                cnf.append([-U[i], -U[j]])
            elif ms is False:
                cnf.append([U[i], U[j]])
            else:
                t = U[(i + j) // k]
                cnf.append([-U[i], -U[j], -t])
                cnf.append([U[i], U[j], t])
    s = Cadical153(bootstrap_with=cnf)
    if not s.solve():
        return False, None
    mdl = set(l for l in s.get_model() if l > 0)
    u = [1 if U[q] in mdl else 0 for q in range(1, qmax + 1)]
    return True, u
