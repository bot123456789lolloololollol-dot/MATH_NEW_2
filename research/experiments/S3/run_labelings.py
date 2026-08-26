"""S3 unified exhaustive labeling runner (retry session).

Four published OPEN conjectures, finite consequence per bound N:
"every unlabeled free tree on n <= N vertices admits labeling L".

  odd   : Gnanajothi 1991 odd-graceful  f: V->{0..2m-1} inj., edge |diffs| = odd 1..2m-1
  prime : Entringer ~1980               f: V->{1..n} bij., gcd(f(u),f(v))=1 on edges
  anti  : Hartsfield-Ringel 1990        g: E->{1..m} bij., vertex sums pairwise distinct
  nprime: Ryan 2014 neighborhood-prime  f: V->{1..n} bij., gcd of N(v) labels = 1 for deg>=2

Determinism: no randomness anywhere; tree stream from lib_s3.gen_free_trees
(OEIS A000055-validated). Checkers are written directly from the definitions
and are independent of the solvers. Every level result is appended to a JSONL
log immediately after the level finishes (crash-safe). A level counts toward a
verdict only if recorded "complete": true.

Usage:
  python run_labelings.py --mode odd|prime|anti|nprime --nmin A --nmax B \
      --deadline 240 --outdir DIR [--selftest]
"""
import argparse
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_s3 import gen_free_trees  # noqa: E402


def tree_data(n, edges):
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    order, seen = [0], [False] * n
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


# ======================================================================
# ODD-GRACEFUL  (parity reduction + MRV CSP + forward checking; Z3 tier)
# ======================================================================

def solve_odd(n, edges, node_budget=2_000_000):
    """Parity argument: all edge differences odd => adjacent vertices carry
    opposite parities => the even-labeled set is exactly one bipartition
    class; both orientations tried."""
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
    dom = 2 * m
    stats = {"nodes": 0}
    full = ((1 << dom) - 1)
    par_mask = [
        sum(1 << x for x in range(dom) if x % 2 == p) for p in (0, 1)
    ]
    lab = [-1] * n
    cand = [par_mask[(color[v] ^ flip) & 1] for v in range(n)]

    class Budget(Exception):
        pass

    def rec(used_diff):
        stats["nodes"] += 1
        if stats["nodes"] > node_budget:
            raise Budget()
        best_v, best_cnt = -1, 999
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
            return True
        v = best_v
        pm = cand[v]
        while pm:
            x = (pm & -pm).bit_length() - 1
            pm &= pm - 1
            newdiffs, ok = 0, True
            for u in adj[v]:
                if lab[u] >= 0:
                    d = abs(x - lab[u])
                    b = 1 << ((d - 1) >> 1)
                    if (used_diff | newdiffs) & b:
                        ok = False
                        break
                    newdiffs |= b
            if not ok:
                continue
            dd = used_diff & ~newdiffs          # earlier-used diffs
            forb = 0
            fd = dd
            while fd:
                i = (fd & -fd).bit_length() - 1
                fd &= fd - 1
                d = 2 * i + 1
                if x - d >= 0:
                    forb |= 1 << (x - d)
                if x + d < dom:
                    forb |= 1 << (x + d)
            saved = []
            bad = False
            lab[v] = x
            for u in adj[v]:
                if lab[u] >= 0:
                    continue
                old = cand[u]
                newc = old & ~forb & ~(1 << x)
                if newc != old:
                    cand[u] = newc
                    saved.append((u, old))
                if newc == 0:
                    bad = True
                    break
            if not bad:
                bad = not _edges_feasible(n, adj, lab, cand, full)
            if not bad and rec(used_diff | newdiffs):
                return True
            lab[v] = -1
            for (u, old) in saved:
                cand[u] = old
        return False

    try:
        ok = rec(0)
    except Budget:
        return "BUDGET"
    return list(lab) if ok else None


def _edges_feasible(n, nbrs, lab, cand, full):
    """Every free-free edge needs a pair of candidates of opposite parity."""
    ev = 0xAAAAAAAAAAAAAAAA & full      # bits at odd indices = odd labels
    for v in range(n):
        if lab[v] >= 0:
            continue
        for u in nbrs[v]:
            if u <= v or lab[u] >= 0:
                continue
            pv, pu = cand[v], cand[u]
            if not (((pv & ev) & (pu & ~ev)) or ((pv & ~ev) & (pu & ev))):
                return False
    return True


def solve_odd_z3(n, edges, timeout_sec=30):
    import z3
    m = n - 1
    L = [z3.Int(f"x{i}") for i in range(n)]
    s = z3.Solver()
    s.set("timeout", int(timeout_sec * 1000))
    for i in range(n):
        s.add(z3.Or([L[i] == x for x in range(2 * m)]))
    for (u, v) in edges:
        d = z3.Int(f"d{u}_{v}")
        s.add(d == z3.If(L[u] - L[v] >= 0, L[u] - L[v], L[v] - L[u]))
        s.add(z3.Or([d == 2 * j + 1 for j in range(m)]))
    s.add(z3.Distinct(L))
    if s.check() == z3.sat:
        mod = s.model()
        return [mod.eval(L[i]).as_long() for i in range(n)]
    return None


def check_odd(n, edges, f):
    m = n - 1
    if len(set(f)) != n:
        return False, "not injective"
    if max(f) > 2 * m - 1 or min(f) < 0:
        return False, "label out of range"
    ds = sorted(abs(f[u] - f[v]) for u, v in edges)
    want = list(range(1, 2 * m, 2))
    return (ds == want), ("diffs" if ds != want else "")


# ======================================================================
# PRIME
# ======================================================================

def solve_prime(n, edges, node_budget=2_000_000):
    """MRV CSP over labels 1..n with bitmask forward checking:
    assigning label x removes from every unassigned neighbor's domain all
    labels NOT coprime to x (precomputed coprime masks)."""
    if n == 1:
        return [1]
    adj, _, _ = tree_data(n, edges)
    full = (1 << (n + 1)) - 2                    # bits 1..n
    cop = [0] * (n + 1)
    for x in range(1, n + 1):
        m = 0
        for y in range(1, n + 1):
            if math.gcd(x, y) == 1:
                m |= 1 << y
        cop[x] = m
    lab = [-1] * n
    cand = [full for _ in range(n)]
    stats = {"nodes": 0}

    class Budget(Exception):
        pass

    def rec(unassigned):
        stats["nodes"] += 1
        if stats["nodes"] > node_budget:
            raise Budget()
        best_v, best_cnt = -1, 10 ** 9
        for v in range(n):
            if lab[v] >= 0:
                continue
            c = bin(cand[v]).count("1")
            if c < best_cnt:
                best_cnt, best_v = c, v
                if c == 1:
                    break
        if best_v < 0:
            return True
        v = best_v
        pm = cand[v]
        while pm:
            x = (pm & -pm).bit_length() - 1
            pm &= pm - 1
            saved = []
            bad = False
            lab[v] = x
            keep = cop[x]
            for u in adj[v]:
                if lab[u] >= 0:
                    continue
                old = cand[u]
                newc = old & keep
                if newc != old:
                    cand[u] = newc
                    saved.append((u, old))
                if newc == 0:
                    bad = True
                    break
            if not bad and rec(unassigned - 1):
                return True
            lab[v] = -1
            for (u, old) in saved:
                cand[u] = old
        return False

    try:
        found = rec(n)
    except Budget:
        return "BUDGET"
    return list(lab) if found else None


def check_prime(n, edges, f):
    if sorted(f) != list(range(1, n + 1)):
        return False, "not bijective"
    for u, v in edges:
        if math.gcd(f[u], f[v]) != 1:
            return False, f"gcd>1 on ({u},{v})"
    return True, ""


# ======================================================================
# ANTIMAGIC
# ======================================================================

def _edge_bfs_order(n, edges):
    adj, order, pos = tree_data(n, edges)
    keyed = []
    for i, (u, v) in enumerate(edges):
        keyed.append((max(pos[u], pos[v]), min(pos[u], pos[v]), i))
    keyed.sort()
    return [i for (_a, _b, i) in keyed], adj


class AntiBudget(Exception):
    pass


def solve_anti(n, edges, node_budget=4_000_000):
    """Backtrack over edges (BFS order) assigning unused labels 1..m.

    Soundness of pruning: the ONLY pruning is (a) a vertex whose incident
    edges are all labeled gets its FINAL sum, which must differ from every
    other final sum. No interval-based rejection is used: an unfinished
    vertex merely *can* reach a range of sums, and reachability alone never
    implies a collision. Hence no valid labeling is ever discarded.
    Leaves finalize early (most vertices of a tree are leaves), which gives
    strong incremental filtering."""
    m = n - 1
    if n == 1:
        return []                     # single vertex: vacuously distinct
    if n == 2:
        return None                   # K2: both endpoints get the same sum
    eorder, adj = _edge_bfs_order(n, edges)
    inc_count = [len(adj[v]) for v in range(n)]
    elab = [0] * m                    # indexed by ORIGINAL edge id
    vsum = [0] * n
    done_ct = [0] * n
    finals = {}
    stats = {"nodes": 0}

    def rec(idx, unused):
        stats["nodes"] += 1
        if stats["nodes"] > node_budget:
            raise AntiBudget()
        if idx == m:
            return len(finals) == n and len(set(finals.values())) == n
        ei = eorder[idx]
        u, v = edges[ei]
        for x in unused:
            done_ct[u] += 1
            done_ct[v] += 1
            vsum[u] += x
            vsum[v] += x
            newly_final = []
            for w in (u, v):
                if done_ct[w] == inc_count[w]:
                    newly_final.append(w)
            ok = True
            seen_new = set()
            for w in newly_final:
                s = vsum[w]
                if s in finals.values() or s in seen_new:
                    ok = False
                    break
                seen_new.add(s)
            if ok:
                elab[ei] = x
                for w in newly_final:
                    finals[w] = vsum[w]
                if rec(idx + 1, [z for z in unused if z != x]):
                    return True
                for w in newly_final:
                    del finals[w]
                elab[ei] = 0
            done_ct[u] -= 1
            done_ct[v] -= 1
            vsum[u] -= x
            vsum[v] -= x
        return False

    try:
        ok = rec(0, list(range(1, m + 1)))
    except AntiBudget:
        return "BUDGET"
    if not ok:
        return None
    return list(elab)


def check_anti(n, edges, g):
    m = n - 1
    if sorted(g) != list(range(1, m + 1)):
        return False, "not a permutation of 1..m"
    sums = [0] * n
    for i, (u, v) in enumerate(edges):
        sums[u] += g[i]
        sums[v] += g[i]
    return len(set(sums)) == n, ("vertex sums not distinct" 
                                 if len(set(sums)) != n else "")


# ======================================================================
# NEIGHBORHOOD-PRIME
# ======================================================================

def _primes_to(n):
    sieve = [True] * (n + 1)
    ps = []
    for p in range(2, n + 1):
        if sieve[p]:
            ps.append(p)
            for q in range(p * p, n + 1, p):
                sieve[q] = False
    return ps


def solve_nprime(n, edges, node_budget=2_000_000):
    if n <= 2:
        return None                    # conjecture starts at n >= 3
    adj, order, _ = tree_data(n, edges)
    lab = [0] * n
    used = [False] * (n + 1)
    primes = _primes_to(n)
    g = [1] * n                        # running gcd of assigned nbr labels
    rem = [len(adj[v]) for v in range(n)]
    stats = {"nodes": 0}

    class Budget(Exception):
        pass

    def rec(i):
        stats["nodes"] += 1
        if stats["nodes"] > node_budget:
            raise Budget()
        if i == n:
            return True
        v = order[i]
        for x in range(1, n + 1):
            if used[x]:
                continue
            lab[v] = x
            used[x] = True
            touched = [(w, g[w]) for w in adj[v] if len(adj[w]) >= 2]
            ok = True
            for w, old in touched:
                ng = math.gcd(old, x)
                g[w] = ng
                rem[w] -= 1
                if ng > 1 and rem[w] == 0:
                    ok = False
                elif ng > 1:
                    # if every still-unused label shares a prime p|ng with all
                    # options, gcd can never drop to 1 along remaining nbrs
                    for p in primes:
                        if ng % p:
                            continue
                        if all(used[y] or y % p == 0 for y in range(1, n + 1)):
                            ok = False
                            break
            if ok and rec(i + 1):
                return True
            for w, old in touched:
                g[w] = old
                rem[w] += 1
            used[x] = False
            lab[v] = 0
        return False

    try:
        found = rec(0)
    except Budget:
        return "BUDGET"
    return list(lab) if found else None


def check_nprime(n, edges, f):
    if sorted(f) != list(range(1, n + 1)):
        return False, "not bijective"
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    for v in range(n):
        if len(adj[v]) >= 2:
            gg = 0
            for u in adj[v]:
                gg = math.gcd(gg, f[u])
            if gg != 1:
                return False, f"vertex {v} neighborhood gcd {gg}"
    return True, ""


# ======================================================================
# modes registry
# ======================================================================

MODES = {
    "odd":    (solve_odd,  check_odd),
    "prime":  (solve_prime, check_prime),
    "anti":   (solve_anti, check_anti),
    "nprime": (solve_nprime, check_nprime),
}

# levels outside a conjecture's domain (K2 excluded by Hartsfield-Ringel;
# neighborhood-prime defined for n >= 3; n=1 odd-graceful has m=0 and an
# empty label range, treated as vacuous)
SKIP_LEVELS = {"odd": {1}, "prime": set(), "anti": {2}, "nprime": {1, 2}}


def corrupt_variants(mode, n, edges, witness):
    """Provably-invalid mutations of a valid witness (for checker tests)."""
    out = []
    if mode == "anti":
        m = len(witness)
        if m >= 2:
            bad = list(witness)
            bad[1] = bad[0]                     # duplicate edge label
            out.append(("dup_edge_label", bad))
        bad = list(witness)
        bad[0] = witness[0] + 10 * (n + 1)       # out of range
        out.append(("out_of_range", bad))
    else:
        bad = list(witness)
        if mode == "odd":
            bad[0] = bad[1] if n > 1 else bad[0] + 1   # duplicate label
        else:
            mx = max(bad)
            bad[0] = mx + 7                          # breaks bijection/range
        out.append(("dup_or_range", bad))
    return out


def selftest():
    """Checkers must ACCEPT solver witnesses and REJECT corrupted variants;
    includes tiny known cases."""
    cases = {
        "odd":    (5, [(0, 1), (0, 2), (0, 3), (0, 4)]),     # star S4
        "prime":  (5, [(0, 1), (0, 2), (0, 3), (0, 4)]),
        "anti":   (6, [(0, 1), (0, 2), (0, 3), (1, 4), (2, 5)]),
        "nprime": (6, [(0, 1), (0, 2), (0, 3), (1, 4), (2, 5)]),
    }
    for mode, (n, edges) in cases.items():
        solve, check = MODES[mode]
        w = solve(n, edges)
        assert w not in (None, "BUDGET"), f"{mode}: no witness on tiny case"
        ok, why = check(n, edges, w)
        assert ok, f"{mode}: checker rejected own witness: {why}"
        assert check(n, edges, list(w)[::-1])[0] in (True, False)  # runs
        for tag, bad in corrupt_variants(mode, n, edges, w):
            okb, whyb = check(n, edges, bad)
            assert not okb, \
                f"{mode}: checker ACCEPTED corrupted variant {tag}: {bad}"
        print(f"selftest {mode}: witness ok, "
              f"{len(corrupt_variants(mode, n, edges, w))} corruptions rejected")
    assert solve_anti(2, [(0, 1)]) is None, "K2 must be non-antimagic"
    assert solve_odd(1, []) == [0]
    assert solve_prime(1, []) == [1]
    assert check_nprime(3, [(0, 1), (0, 2)],
                        solve_nprime(3, [(0, 1), (0, 2)]))[0]
    print("selftest passed")


# ======================================================================
# runner
# ======================================================================

def run_mode(mode, nmin, nmax, deadline_s, outdir):
    solve, check = MODES[mode]
    os.makedirs(outdir, exist_ok=True)
    log_path = os.path.join(outdir, f"log_{mode}.jsonl")
    t_start = time.time()
    print(f"[{mode}] generating free trees to n={nmax} ...", flush=True)
    trees = gen_free_trees(nmax)
    print(f"[{mode}] generation done in {time.time()-t_start:.1f}s",
          flush=True)
    for n in range(nmin, nmax + 1):
        if time.time() - t_start > deadline_s:
            print(f"[{mode}] deadline before level n={n}; stopping",
                  flush=True)
            return
        if n in SKIP_LEVELS.get(mode, ()):
            with open(log_path, "a") as fh:
                fh.write(json.dumps(
                    {"mode": mode, "n": n, "skipped": True,
                     "note": "outside conjecture domain"}) + "\n")
            continue
        t0 = time.time()
        cnt = 0
        failures, bugs, budget_hits, z3_used = [], [], 0, 0
        worst_ms, worst_key = 0.0, None
        sample = None
        aborted = False
        for edges in trees[n]:
            cnt += 1
            ta = time.time()
            w = solve(n, edges)
            dtms = (time.time() - ta) * 1000.0
            if dtms > worst_ms:
                worst_ms, worst_key = dtms, sorted(edges)
            if w == "BUDGET":
                if mode == "odd":
                    wz = solve_odd_z3(n, edges)
                    if wz is None:
                        failures.append({"edges": sorted(edges),
                                         "note": "z3 unsat/timeout"})
                        print(f"CANDIDATE-COUNTEREXAMPLE n={n} "
                              f"{sorted(edges)}", flush=True)
                        continue
                    z3_used += 1
                    w = wz
                else:
                    budget_hits += 1
                    continue
            if w is None:
                failures.append({"edges": sorted(edges), "note": "unsat"})
                print(f"COUNTEREXAMPLE-FROM-SOLVER n={n} {sorted(edges)}",
                      flush=True)
                continue
            ok, why = check(n, edges, w)
            if not ok:
                bugs.append({"edges": sorted(edges), "w": w, "why": why})
                print(f"SOLVER-BUG n={n} {sorted(edges)}: {why}", flush=True)
                continue
            if sample is None:
                sample = {"edges": sorted(edges), "labeling": w}
            remaining = deadline_s - (time.time() - t_start)
            if remaining <= 0 and cnt < len(trees[n]):
                aborted = True
                break
        clean = (not aborted) and (not failures) and (not bugs) \
            and budget_hits == 0
        rec = {
            "mode": mode, "n": n,
            "complete": (not aborted),
            "clean": clean,
            "trees": cnt, "seconds": round(time.time() - t0, 2),
            "worst_ms": round(worst_ms, 2),
            "failures": failures[:20], "solver_bugs": bugs[:20],
            "budget_hits": budget_hits, "z3_escalations": z3_used,
            "sample_witness": sample,
        }
        with open(log_path, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        status = "OK" if clean else "ATTENTION"
        print(f"[{mode}] n={n:2d} {'COMPLETE' if not aborted else 'ABORTED'}"
              f" {cnt} trees {rec['seconds']}s worst {worst_ms:.0f}ms "
              f"[{status}]", flush=True)
        if failures or bugs:
            return


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=list(MODES) + ["all"])
    ap.add_argument("--nmin", type=int, default=1)
    ap.add_argument("--nmax", type=int, default=19)
    ap.add_argument("--deadline", type=float, default=240,
                    help="wall-clock budget seconds for this invocation")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        if args.mode is None:
            sys.exit(0)
    modes = list(MODES) if args.mode in (None, "all") else [args.mode]
    for md in modes:
        run_mode(md, args.nmin, args.nmax, args.deadline, args.outdir)
