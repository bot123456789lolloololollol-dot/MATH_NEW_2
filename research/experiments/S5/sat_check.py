"""SAT-based symbolic verifier for sorting networks (Campaign C1).

Encodes "there EXISTS an input on which the network fails to sort" as CNF:
input bits are variables, each compare-exchange contributes Tseitin-defined
min/max gates, and a disjunction over adjacent-output violations
(out_i=1 AND out_{i+1}=0) is asserted.  The encoding is equisatisfiable to
the negation of sortness by construction, hence:
    solver SAT  <=>  network does NOT sort (counterexample extracted)
    solver UNSAT <=> network sorts ALL inputs (exhaustive certificate)

This is a third verification path, independent of both the bigint bitset
simulator and the numpy simulator, and its cost is O(L*n) clauses rather
than O(2^n) memory/time.
"""
from pysat.solvers import Cadical153
from pysat.formula import IDPool


def sat_sorts(net, n, want_model=False):
    """Return True iff net sorts all inputs.  With want_model, also return
    a violating input assignment when the network fails."""
    pool = IDPool()
    cur = [pool.id(("v", i)) for i in range(n)]
    clauses = []

    def AND_gate(y, u, v):
        clauses.extend([[-u, -v, y], [u, -y], [v, -y]])

    def OR_gate(y, u, v):
        clauses.extend([[u, v, -y], [-u, y], [-v, y]])

    for (i, j) in net:
        u, v = cur[i], cur[j]
        mn = pool.id(("mn", i, j, len(clauses)))
        mx = pool.id(("mx", i, j, len(clauses)))
        AND_gate(mn, u, v)
        OR_gate(mx, u, v)
        cur[i] = mn
        cur[j] = mx

    viol = []
    for i in range(n - 1):
        a = pool.id(("bad", i))
        AND_gate(a, cur[i], -cur[i + 1])   # a <-> out_i AND NOT out_{i+1}
        viol.append(a)
    clauses.append(viol)                   # some adjacent inversion exists

    s = Cadical153(bootstrap_with=clauses)
    r = s.solve()
    model = s.get_model() if (r and want_model) else None
    s.delete()
    if r and want_model:
        val = {abs(l): (l > 0) for l in model}
        inp = tuple(1 if val[pool.id(("v", i))] else 0 for i in range(n))
        return False, inp
    return r is False, None


if __name__ == "__main__":
    import sys
    import json
    import netlib
    import baselines

    # self-consistency battery
    print("n  net                 bitset   SAT     expect")
    cases = []
    for n in range(2, 9):
        oem = baselines.batcher_oem_net(n)
        cases.append((n, oem, True))
        cases.append((n, oem[:-1], None))          # truncated: usually broken
        rng = __import__("random").Random(n)
        rnd = []
        for _ in range(n * (n - 1) // 4):
            i = rng.randrange(n - 1)
            j = rng.randint(i + 1, n - 1)
            rnd.append((i, j))
        cases.append((n, rnd, None))
        cases.append((n, baselines.KNUTH_SMALL[n] if n in baselines.KNUTH_SMALL
                      else oem, True))
    mismatches = 0
    for n, net, expect in cases:
        b = netlib.is_sorting(net, n)
        s, _ = sat_sorts(net, n)
        flag = "" if expect is None or b == expect else " UNEXPECTED"
        if expect is not None and b != expect:
            mismatches += 1
        if b != s:
            mismatches += 1000
        print(f"{n:<2} {str(net[:3])+'...':<19} {str(b):<7} {str(s):<7}"
              f"{'':6}{flag}")
    print("MISMATCHES:", mismatches)
