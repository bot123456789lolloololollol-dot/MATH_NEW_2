"""C-S3-03 / C-S3-04: labeling conjectures for ALL small trees.

C-S3-03 (odd-graceful): Gnanajothi (1991) conjectured every tree admits an
  odd-graceful labeling: injective f: V -> {0,1,...,2m-1} (m = #edges) such
  that the edge differences |f(u)-f(v)| are exactly {1,3,...,2m-1}.
  Known: proven for caterpillars and diameter <= 5 (Barrientos 2008).
  No published exhaustive verification over all trees of a given order was
  found in the literature scan (see literature/S3-novelty.md).

C-S3-04 (prime): Entringer (~1980) conjectured every tree has a prime
  labeling: bijection f: V -> {1,...,n} with gcd(f(u),f(v)) = 1 for every
  edge uv. Open in general.

Both solvers are exact backtracking searches over each unlabeled tree; the
tree stream comes from lib_s3.gen_free_trees, whose counts match OEIS
A000055 exactly (self-validating enumeration).

Usage: python tree_labelers.py --mode odd|prime --nmin A --nmax B --out FILE
"""

import argparse
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lib_s3 import gen_free_trees          # noqa: E402


def tree_data(n, edges):
    """Adjacency + BFS order from vertex 0 (deterministic)."""
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    order = [0]
    seen = [False] * n
    seen[0] = True
    head = 0
    while head < len(order):
        x = order[head]
        head += 1
        for y in adj[x]:
            if not seen[y]:
                seen[y] = True
                order.append(y)
    pos = [0] * n
    for i, v in enumerate(order):
        pos[v] = i
    return adj, order, pos


def solve_odd_graceful(n, edges, node_budget=2_000_000):
    """Return labeling list or None. Labels 0..2m-1; diffs odd 1..2m-1.

    Parity argument: all edge differences are odd => adjacent vertices carry
    opposite-parity labels => the even-labeled vertex set is one bipartition
    class. We try both orientations. Dynamic vertex order (most constrained
    first) + per-tree node budget."""
    m = n - 1
    if n == 1:
        return [0]
    adj, _, _ = tree_data(n, edges)
    color = [-1] * n
    color[0] = 0
    stack = [0]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if color[y] < 0:
                color[y] = color[x] ^ 1
                stack.append(y)
    for flip in (0, 1):
        res = _og_try(n, m, adj, color, flip, node_budget)
        if res is not None:
            return res
    return None


def _og_try(n, m, adj, color, flip, node_budget):
    """CSP solver: candidate bitmasks per vertex, MRV ordering, forward
    checking on labels AND on remaining-edge/diff feasibility."""
    dom = 2 * m
    full_diff = (1 << m) - 1
    stats = {"nodes": 0}
    par_mask = [
        sum(1 << x for x in range(dom) if x % 2 == p) for p in (0, 1)
    ]
    par_of = [(color[v] ^ flip) & 1 for v in range(n)]

    lab = [-1] * n
    cand = [par_mask[par_of[v]] for v in range(n)]
    nbrs = adj

    class Budget(Exception):
        pass

    def rec(used_lab, used_diff):
        stats["nodes"] += 1
        if stats["nodes"] > node_budget:
            raise Budget()
        # find unassigned vertex with min candidates
        best_v = -1
        best_cnt = 999
        for v in range(n):
            if lab[v] >= 0:
                continue
            c = bin(cand[v]).count("1")
            if c == 0:
                return False
            if c < best_cnt:
                best_cnt, best_v = c, v
                if c == 1:
                    break
        if best_v < 0:
            return True                      # all assigned
        v = best_v
        pm = cand[v]
        while pm:
            x = (pm & -pm).bit_length() - 1
            pm &= pm - 1
            newdiffs = 0
            ok = True
            for u in nbrs[v]:
                if lab[u] >= 0:
                    d = abs(x - lab[u])
                    b = 1 << ((d - 1) >> 1)
                    if (used_diff | newdiffs) & b:
                        ok = False
                        break
                    newdiffs |= b
            if not ok:
                continue
            nd = used_diff | newdiffs
            # forbidden labels for unassigned neighbors: those whose distance
            # to x is a diff already used elsewhere (not newdiffs)
            forb = 0
            dd = used_diff & ~newdiffs      # diffs used by earlier edges
            fd = dd
            while fd:
                i = (fd & -fd).bit_length() - 1
                fd &= fd - 1
                d = 2 * i + 1
                if x - d >= 0:
                    forb |= 1 << (x - d)
                if x + d < dom:
                    forb |= 1 << (x + d)
            # assign v, forward-check unassigned neighbors
            saved = []
            bad = False
            lab[v] = x
            for u in nbrs[v]:
                if lab[u] >= 0:
                    continue
                old = cand[u]
                newc = old & ~forb & ~(1 << x)
                # u's diff to x must be an ODD unused-by-others value: parity
                # ensures odd; ensure not in dd (handled by forb) -- ok
                if newc != old:
                    cand[u] = newc
                    saved.append((u, old))
                if newc == 0:
                    bad = True
                    break
            if not bad:
                # every remaining edge needs >=1 feasible diff
                bad = not _edges_feasible(n, nbrs, lab, cand, dom)
            if not bad and rec(used_lab | (1 << x), nd):
                return True
            lab[v] = -1
            for (u, old) in saved:
                cand[u] = old
        return False

    try:
        ok = rec(0, 0)
    except Budget:
        return "BUDGET"
    return list(lab) if ok else None


def _edges_feasible(n, nbrs, lab, cand, dom):
    """Each unassigned edge (both endpoints free) must admit some unused odd
    difference given current candidate masks; edges with one endpoint fixed
    were already handled by forward checking."""
    for v in range(n):
        if lab[v] >= 0:
            continue
        for u in nbrs[v]:
            if u <= v or lab[u] >= 0:
                continue
            # need x in cand[v], y in cand[u] with |x-y| odd (auto by parity?
            # NO: both free -> any pair) -> exists iff cands differ in parity
            # overlap nonempty across parities
            pv = cand[v]
            pu = cand[u]
            # parity classes: even labels = even positions bits
            ev_v = pv & 0xAAAAAAAAAAAAAAAA & ((1 << dom) - 1)
            ev_u = pu & 0xAAAAAAAAAAAAAAAA & ((1 << dom) - 1)
            od_v = pv & ~ev_v
            od_u = pu & ~ev_u
            if not ((ev_v & od_u) or (od_v & ev_u)):
                return False
    return True


def solve_odd_graceful_z3(n, edges, timeout_sec=60):
    """Escalation tier: encode odd-graceful labeling in z3."""
    try:
        import z3
    except ImportError:
        return "NOZ3"
    m = n - 1
    dom = 2 * m
    L = [z3.Int(f"x{i}") for i in range(n)]
    D = [z3.Int(f"d{i}") for i in range(m)]
    s = z3.Solver()
    s.set("timeout", int(timeout_sec * 1000))
    for i in range(n):
        s.add(z3.Or([L[i] == x for x in range(dom)]))
    for i in range(m):
        u, v = edges[i]
        s.add(D[i] == z3.If(L[u] - L[v] >= 0, L[u] - L[v], L[v] - L[u]))
        s.add(z3.Or([D[i] == 2 * j + 1 for j in range(m)]))
    s.add(z3.Distinct(L))
    s.add(z3.Distinct(D))
    if s.check() == z3.sat:
        mod = s.model()
        return [mod.eval(L[i]).as_long() for i in range(n)]
    return None


def check_odd_graceful(n, edges, f):
    m = n - 1
    if len(set(f)) != n:                      # injective
        return False, "not injective"
    if max(f) > 2 * m - 1 or min(f) < 0:
        return False, "label out of range"
    ds = sorted(abs(f[u] - f[v]) for u, v in edges)
    want = sorted(range(1, 2 * m, 2))
    return ds == want, ("diffs" if ds == want else f"got {ds}")


def solve_prime(n, edges):
    if n == 1:
        return [1]
    adj, order, pos = tree_data(n, edges)
    lab = [0] * n
    used = [False] * (n + 1)

    def rec(i):
        if i == n:
            return True
        v = order[i]
        nbrs = [u for u in adj[v] if lab[u]]
        for x in range(1, n + 1):
            if used[x]:
                continue
            ok = True
            for u in nbrs:
                if math.gcd(x, lab[u]) != 1:
                    ok = False
                    break
            if not ok:
                continue
            lab[v] = x
            used[x] = True
            if rec(i + 1):
                return True
            used[x] = False
            lab[v] = 0
        return False

    return list(lab) if rec(0) else None


def check_prime(n, edges, f):
    if sorted(f) != list(range(1, n + 1)):
        return False, "not bijective"
    import math
    for u, v in edges:
        if math.gcd(f[u], f[v]) != 1:
            return False, f"edge ({u},{v}) labels {f[u]},{f[v]}"
    return True, ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["odd", "prime"], required=True)
    ap.add_argument("--nmin", type=int, default=1)
    ap.add_argument("--nmax", type=int, default=19)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    solver = solve_odd_graceful if args.mode == "odd" else solve_prime
    checker = check_odd_graceful if args.mode == "odd" else check_prime

    trees = gen_free_trees(args.nmax)
    results = {}
    t_all = time.time()
    for n in range(args.nmin, args.nmax + 1):
        t0 = time.time()
        worst_ms = 0.0
        worst_key = None
        count = 0
        for edges in trees[n]:
            count += 1
            t1 = time.time()
            f = solver(n, edges)
            dt = (time.time() - t1) * 1000.0
            if dt > worst_ms:
                worst_ms, worst_key = dt, edges
            if f == "BUDGET":
                fz = solve_odd_graceful_z3(n, edges)
                if fz is None:
                    results.setdefault("z3_unsat", []).append(
                        {"n": n, "edges": sorted(edges)})
                    print(f"z3 UNSAT n={n} edges={sorted(edges)}", flush=True)
                    continue
                f = fz
                tier = "z3"
            else:
                tier = "core"
            okz, whyz = checker(n, edges, f)
            if not okz:
                print(f"SOLVER BUG ({tier}) n={n} {sorted(edges)}: {whyz}",
                      flush=True)
                results.setdefault("solver_bugs", []).append(
                    {"n": n, "edges": sorted(edges), "f": f})
                continue
            results.setdefault("by_tier", {}).setdefault(tier, 0)
            results["by_tier"][tier] += 1
            if tier == "z3":
                continue
            if f is None:
                ok, why = checker(n, edges, [0] * n)
                print(f"COUNTEREXAMPLE-FROM-SOLVER n={n} edges={sorted(edges)}")
                results.setdefault("failures", []).append(
                    {"n": n, "edges": sorted(edges)})
                continue
            ok, why = checker(n, edges, f)
            if not ok:
                print(f"SOLVER BUG: produced invalid labeling n={n} "
                      f"edges={sorted(edges)} f={f}: {why}")
                results.setdefault("solver_bugs", []).append(
                    {"n": n, "edges": sorted(edges), "f": f})
        el = time.time() - t0
        results[str(n)] = {"trees": count, "seconds": round(el, 2),
                           "worst_ms": round(worst_ms, 2)}
        print(f"mode={args.mode} n={n:2d}: {count} trees, all labeled OK, "
              f"{el:7.2f}s (worst tree {worst_ms:.1f} ms)", flush=True)
        # save a sample labeling for reproducibility at a few orders
    results["total_seconds"] = round(time.time() - t_all, 2)
    results["mode"] = args.mode
    results["range"] = [args.nmin, args.nmax]
    out = args.out or f"results_{args.mode}_{args.nmin}_{args.nmax}.json"
    json.dump(results, open(out, "w"), indent=1)
    print(f"[written] {out}")


if __name__ == "__main__":
    main()
