"""Exact decision procedure: does a k-selection network on n wires with at
most L compare-exchange elements exist?  SAT-encoded.

Model: L slots; each slot selects exactly one "action" via one-hot variable
y[t,c] over candidate comparators c=(i,j), i<j, plus an IDLE action.  Any
network with L' <= L real comparators is representable (pad the sequence
with IDLE slots), so
    SAT  at L  =>  exists selector with <= L CEs   (witness decoded & verified)
    UNSAT at L =>  no selector with <= L CEs       (certified lower bound)

Per input pattern p (ALL 2^n enumerated -> complete), wire values are
constants at slot entry; slot semantics are encoded with Tseitin gates and
gated AND/OR terms; the selection property is enforced per pattern as a
cardinality constraint on the bottom-k block:
      #{i<k: out_i(p)=1} = max(0, c_p-(n-k)),  c_p = popcount(p)  [constant]

Sorting mode (--property sort) replaces that constraint by pairwise output
order + is validated against literature-proven optimal sorter sizes.

Symmetry breaking used (soundness argued in report): pairs restricted to
i<j; nothing else -- no completeness-unsafe constraints.
"""
from pysat.solvers import Cadical153
from pysat.formula import IDPool, CNF
from pysat.card import CardEnc, EncType


def build_and_encode(n, k, L, property="select", verbose=False):
    pool = IDPool()
    cnf = CNF()
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    actions = pairs + [None]                       # None = idle
    A = len(actions)
    y = {}
    for t in range(L):
        lits = []
        for ai in range(A):
            v = pool.id(("y", t, ai))
            y[(t, ai)] = v
            lits.append(v)
        cnf.extend(CardEnc.atleast(lits=lits, bound=1, vpool=pool).clauses)
        cnf.extend(CardEnc.atmost(lits=lits, bound=1, vpool=pool).clauses)
    # symmetry breaking (both preserve completeness):
    # (1) no two adjacent identical comparators -- deleting the first of a
    #     duplicate pair leaves behavior identical, so minimal nets never
    #     contain them;
    # (2) IDLE slots form a suffix -- deleting an interior idle slot yields
    #     an equivalent network.
    for ai in range(A - 1):
        for t in range(1, L):
            cnf.append([-y[(t, ai)], -y[(t - 1, ai)]])
    idle_ai = A - 1
    for t in range(L - 1):
        cnf.append([-y[(t, idle_ai)], y[(t + 1, idle_ai)]])

    total = 1 << n

    def AND2(a, b):
        v = pool.id(("and", a, b))
        cnf.append([-a, -b, v]); cnf.append([a, -v]); cnf.append([b, -v])
        return v

    def AND3(a, b, c):
        v = pool.id(("and3", a, b, c))
        cnf.append([-a, -b, -c, v]); cnf.append([a, -v])
        cnf.append([b, -v]); cnf.append([c, -v])
        return v

    def ORlist(lits_):
        if len(lits_) == 1:
            return lits_[0]
        acc = lits_[0]
        for x in lits_[1:]:
            nv = pool.id(("or", acc, x))
            cnf.append([acc, x, -nv]); cnf.append([-acc, nv]); cnf.append([-x, nv])
            acc = nv
        return acc

    for p in range(total):
        cur = []
        for i in range(n):
            b = (p >> i) & 1
            v = pool.id(("cst", p, i))
            cnf.append([v] if b else [-v])
            cur.append(v)
        for t in range(L):
            terms = [[] for _ in range(n)]     # disjuncts adding 1 to wire w
            touched = [[] for _ in range(n)]   # action indices touching wire w
            for ai in range(A - 1):
                i, j = actions[ai]
                yv = y[(t, ai)]
                terms[i].append(AND3(yv, cur[i], cur[j]))
                terms[j].append(ORlist([AND2(yv, cur[i]), AND2(yv, cur[j])]))
                touched[i].append(yv)
                touched[j].append(yv)
            nxt = []
            for w in range(n):
                if touched[w]:
                    selw = ORlist(touched[w])
                    # ns <-> NOT selw
                    ns = pool.id(("ns", t, p, w))
                    cnf.append([ns, selw]); cnf.append([-ns, -selw])
                    keep = AND2(ns, cur[w])
                else:
                    keep = cur[w]
                nxt.append(ORlist([keep] + terms[w]))
            cur = nxt
        if property == "sort":
            for i in range(n - 1):
                cnf.append([-cur[i], cur[i + 1]])
        else:
            ones_needed = max(0, bin(p).count("1") - (n - k))
            block = cur[:k]
            if ones_needed > 0:
                cnf.extend(CardEnc.atleast(
                    lits=list(block), bound=ones_needed, vpool=pool).clauses)
            cnf.extend(CardEnc.atmost(
                lits=list(block), bound=max(0, ones_needed),
                vpool=pool).clauses)
    return cnf, pool, actions


def decode(model, pool, actions, L):
    val = {abs(x): (x > 0) for x in model}
    net = []
    for t in range(L):
        for ai, a in enumerate(actions[:-1]):
            v = pool.id(("y", t, ai))
            if val.get(v, False):
                net.append(a)
                break
    return net


def decide(n, k, L, property="select", time_limit=None):
    """Return ('SAT', net) or ('UNSAT', None) or ('UNKNOWN', None)."""
    cnf, pool, actions = build_and_encode(n, k, L, property)
    s = Cadical153(bootstrap_with=cnf.clauses)
    if time_limit:
        s.conf_budget(time_limit)
    r = s.solve_limited()
    if r is True:
        return "SAT", decode(s.get_model(), pool, actions, L)
    if r is False:
        return "UNSAT", None
    return "UNKNOWN", None


if __name__ == "__main__":
    import sys
    import time
    import verify_selector as vs
    import netlib

    prop = sys.argv[1] if len(sys.argv) > 1 else "validate"
    if prop == "validate":
        # validate sorting-mode encoder against literature-proven optima
        proven = {3: 3, 4: 5, 5: 9}          # Knuth TAOCP 5.3.4
        ok = True
        for n, s_opt in proven.items():
            t0 = time.time()
            r_hi, net = decide(n, None, s_opt, "sort")
            r_lo, _ = decide(n, None, s_opt - 1, "sort")
            good = (r_hi == "SAT" and r_lo == "UNSAT")
            ok &= good
            print(f"n={n}: SAT@{s_opt}={r_hi}, UNSAT@{s_opt-1}={r_lo} "
                  f"[{time.time()-t0:.1f}s] {'OK' if good else 'FAIL'}")
            if r_hi == "SAT":
                assert netlib.is_sorting(net, n), "decoded witness invalid!"
        print("encoder validation:", "PASSED" if ok else "FAILED")
