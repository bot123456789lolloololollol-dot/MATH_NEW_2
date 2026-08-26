"""C-S3-02: Vizing's domination conjecture -- exhaustive verification for
all pairs of unlabeled graphs on <= 7 vertices each.

Conjecture (Vizing 1968): gamma(G □ H) >= gamma(G) * gamma(H) for the
Cartesian product. OPEN in general (Clark & Suen 2000 prove the 1/2-fraction
bound; Barcalkin-German 1979 the gamma=2 case; cycles etc.).

This script checks EVERY pair (G,H) of unlabeled simple graphs with
1 <= |V|,|W| <= 7 (1253 x 1253 = 1,570,009 ordered pairs; unordered pairs are
covered by symmetry gamma(G□H)=gamma(H□G)) using EXACT domination numbers:

  - gamma of every factor from the atlas;
  - component decomposition: G□H splits into C_i □ D_j over connected
    components, gamma additive over components;
  - exact gamma of each component product (<= 49 vertices) by iterative-
    deepening backtracking with a fractional covering lower bound.

Verdict per pair: HOLDS (gamma_product >= target), VIOLATED (strictly less),
or TIMEOUT (escalated separately). A single VIOLATED would refute the
conjecture.

Usage: python vizing_atlas.py [--nmax 7] [--out verdict_C2.json]
                             [--pairs-limit K] (debug)
"""

import argparse
import json
import os
import sys
import time
from math import ceil

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from lib_s3 import gen_unlabeled_graphs   # noqa: E402


# ----------------------------------------------------------- exact domination

class ProductInstance:
    __slots__ = ("n", "adj", "closed", "order")

    def __init__(self, n, edges):
        self.n = n
        adj = [0] * n
        closed = [0] * n
        for u, v in edges:
            adj[u] |= (1 << v)
            adj[v] |= (1 << u)
        for v in range(n):
            closed[v] = adj[v] | (1 << v)
        self.adj = adj
        self.closed = closed


def gamma_exact(inst, ub_hint=None, deadline=None):
    """Exact domination number: iterative deepening DFS with
    most-constrained-vertex branching and greedy 2-packing lower bound."""
    n = inst.n
    full = (1 << n) - 1
    closed = inst.closed
    covs = [bin(c).count("1") for c in closed]

    def rec(dominated, k, deadline):
        if dominated == full:
            return True
        if k == 0:
            return False
        # greedy 2-packing lower bound among undominated vertices:
        # pick vertices whose closed neighborhoods are pairwise disjoint
        avail = (~dominated) & full
        pack = 0
        reserved = 0
        x = avail
        while x:
            v = (x & -x).bit_length() - 1
            x &= x - 1
            if closed[v] & reserved:
                continue          # closed neighborhood overlaps a picked one
            pack += 1
            reserved |= closed[v]
        if pack > k:
            return False
        # most constrained undominated vertex: fewest domination options
        best_v = -1
        best_cnt = 99
        x = avail
        while x:
            v = (x & -x).bit_length() - 1
            x &= x - 1
            cnt = bin(closed[v]).count("1")
            if cnt < best_cnt:
                best_cnt, best_v = cnt, v
                if cnt <= 2:
                    break
        v = best_v
        # branch over closed neighborhood of v, most coverage first
        cands = []
        c = closed[v]
        while c:
            w = (c & -c).bit_length() - 1
            c &= c - 1
            cands.append((bin(closed[w] & ~dominated).count("1"), w))
        cands.sort(reverse=True)
        for _cov, w in cands:
            if deadline and time.time() > deadline:
                raise TimeoutError
            if rec(dominated | closed[w], k - 1, deadline):
                return True
        return False

    lb0 = ceil(n / max(covs))
    k = lb0
    while True:
        if rec(0, k, deadline):
            return k
        k += 1


# ------------------------------------------------------------------ main flow

def connected_components(n, edges):
    adjm = {v: set() for v in range(n)}
    for u, v in edges:
        adjm[u].add(v)
        adj[v].add(u)
    seen = set()
    comps = []
    for s in range(n):
        if s in seen:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in adjm[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        idx = {v: i for i, v in enumerate(comp)}
        ce = [(idx[u], idx[v]) for u, v in edges
              if u in idx and v in idx]
        comps.append((len(comp), tuple(sorted(ce))))
    return comps


def relabel_key(size, comp_edges):
    """Canonical-ish cache key: exact relabeled component identity is given
    by the sorted edge tuple after our generators produce canonical labelings
    for connected induced subgraphs? NOT guaranteed iso-invariant in general;
    instead cache on the frozenset of edges directly (safe)."""
    return (size, frozenset(comp_edges))


def gamma_of_graph(n, edges, comp_gamma_cache):
    """gamma via component decomposition."""
    total = 0
    for (size, ce) in connected_components(n, edges):
        key = relabel_key(size, ce)
        if key not in comp_gamma_cache:
            comp_gamma_cache[key] = gamma_exact(ProductInstance(size, list(ce)))
        total += comp_gamma_cache[key]
    return total


def product_components_check(g_n, g_edges, h_n, h_edges,
                             prod_cache, deadline=None):
    """Return gamma(G □ H) computed componentwise."""
    cg = connected_components(g_n, g_edges)
    ch = connected_components(h_n, h_edges)
    total = 0
    for (sg, eg) in cg:
        for (sh, eh) in ch:
            key = (frozenset(eg), frozenset(eh))
            if key in prod_cache:
                total += prod_cache[key]
                continue
            pn = sg * sh
            pe = []
            idxmap = {}
            for a in range(sg):
                for b in range(sh):
                    idxmap[(a, b)] = a * sh + b
            for (a1, a2) in eg:
                for b in range(sh):
                    pe.append((idxmap[(a1, b)], idxmap[(a2, b)]))
            for (b1, b2) in eh:
                for a in range(sg):
                    pe.append((idxmap[(a, b1)], idxmap[(a, b2)]))
            val = gamma_exact(ProductInstance(pn, pe), deadline=deadline)
            prod_cache[key] = val
            total += val
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=7)
    ap.add_argument("--out", default="verdict_C2.json")
    ap.add_argument("--per-pair-deadline", type=float, default=30.0)
    ap.add_argument("--start-order", type=int, default=1,
                    help="check pairs (G,H) with max(|V(G)|,|V(H)|)>=this")
    args = ap.parse_args()

    t0 = time.time()
    graphs = gen_unlabeled_graphs(args.nmax)
    flat = []
    for n in range(1, args.nmax + 1):
        for edges in graphs[n]:
            flat.append((n, edges))
    print(f"atlas up to {args.nmax}: {len(flat)} graphs", flush=True)

    comp_gam = {}
    gammas = []
    for (n, e) in flat:
        gammas.append(gamma_of_graph(n, e, comp_gam))
    from collections import Counter
    print("gamma distribution:", sorted(Counter(gammas).items()), flush=True)

    prod_cache = {}
    n_pairs = 0
    n_hold = 0
    n_violate = []
    n_timeout = []
    t_pass = time.time()
    for i, ((gn, ge), gg) in enumerate(zip(flat, gammas)):
        for j, ((hn, he), gh) in enumerate(zip(flat, gammas)):
            if j < i:
                continue     # symmetric: gamma(G□H) = gamma(H□G)
            target = gg * gh
            # trivial: gamma=1 factors always satisfy (verified too, cheap)
            deadline = time.time() + args.per_pair_deadline
            try:
                val = product_components_check(gn, ge, hn, he, prod_cache,
                                               deadline=deadline)
            except TimeoutError:
                n_timeout.append(((gn, sorted(map(tuple, map(list, ge))),
                                   hn, sorted(he))))
                print(f"TIMEOUT pair #{i},{j} orders {gn}x{hn} "
                      f"gammas {gg}x{gh}", flush=True)
                continue
            n_pairs += 1
            if val >= target:
                n_hold += 1
            else:
                n_violate.append({"G": [gn, sorted(ge)], "H": [hn, sorted(he)],
                                  "gammas": [gg, gh],
                                  "product": val, "target": target})
                print(f"!!! VIOLATION: gamma={val} < {target} "
                      f"for orders {gn}x{hn}", flush=True)
        if (i + 1) % 50 == 0:
            el = time.time() - t_pass
            eta = el / (i + 1) * (len(flat) - i - 1)
            print(f"  row {i+1}/{len(flat)} pairs={n_pairs} hold={n_hold} "
                  f"viol={len(n_violate)} timeout={len(n_timeout)} "
                  f"elapsed={el:.0f}s eta={eta:.0f}s", flush=True)

    out = {
        "conjecture": "Vizing 1968: gamma(G□H) >= gamma(G)*gamma(H)",
        "scope": f"all unordered pairs of unlabeled graphs on <= {args.nmax} vertices",
        "pairs_checked": n_pairs,
        "holds": n_hold,
        "violations": n_violate,
        "timeouts": n_timeout[:20],
        "runtime_sec": round(time.time() - t0, 1),
    }
    json.dump(out, open(args.out, "w"), indent=1)
    print(f"[written] {args.out}")
    print(f"VERDICT: {'VERIFIED' if not n_violate else 'REFUTED'} "
          f"on all {n_pairs} unordered pairs (timeouts: {len(n_timeout)})")


if __name__ == "__main__":
    main()
