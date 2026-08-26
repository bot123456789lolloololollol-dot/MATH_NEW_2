"""
solver.py — Reference (cross-check) exact solver for maximum Sidon sets in Z_n.

Same problem and same search protocol as solver.js, written independently in
Python. Used for (a) cross-implementation agreement tests with the JS primary
solver, (b) small-range full runs.

Strict convention: all sums a+b (unordered, doubles included) distinct in Z_n
(OEIS A260999 definition).  Search-side formulation: all C(k,2) circular
distances min((b-a) mod n,(a-b) mod n) distinct AND none equal to n/2 when n
is even.
"""
import argparse, json, sys, time

sys.setrecursionlimit(100000)


def upper_bound_k(n: int) -> int:
    slots = (n - 1) // 2
    k = 1
    while ((k + 1) * k) // 2 <= slots:
        k += 1
    return k


def solve(n: int, k: int, deadline_ms: int):
    half = n // 2
    used = [0] * (half + 1)
    elems = [0] * k          # elems[0] = 0 fixed by translation symmetry
    state = {"nodes": 0, "timeout": False}
    deadline = time.monotonic() * 1000 + deadline_ms
    witness = []

    def dfs(t: int, last: int) -> bool:
        state["nodes"] += 1
        if state["nodes"] % 4096 == 0 and time.monotonic() * 1000 > deadline:
            state["timeout"] = True
            raise TimeoutError
        if t == k:
            witness[:] = elems[:]
            return True
        if (n - 1) - last < k - t:      # counting prune (completeness-safe)
            return False
        added = []
        ok_all = False
        for x in range(last + 1, n):
            ok = True
            del added[:]
            for i in range(t):
                d = x - elems[i]
                if d > half:
                    dd = n - d
                else:
                    dd = d
                if n % 2 == 0 and dd == half:
                    ok = False
                    break
                if used[dd]:
                    ok = False
                    break
                used[dd] = 1
                added.append(dd)
            if ok:
                elems[t] = x
                if dfs(t + 1, x):
                    ok_all = True
                    break
            for s in added:
                used[s] = 0
            if state["timeout"]:
                raise TimeoutError
        if not ok_all:
            # ensure cleanup of last failed candidate already done above
            pass
        else:
            pass
        return ok_all

    try:
        found = dfs(1, 0)
        return {"status": "SAT" if found else "UNSAT",
                "witness": witness[:] if found else None,
                "nodes": state["nodes"]}
    except TimeoutError:
        return {"status": "TIMEOUT", "witness": None, "nodes": state["nodes"]}


def solve_n(n: int, deadline_ms: int):
    rec = {"n": n, "ub": upper_bound_k(n), "levels": [],
           "k_max": None, "witness": None}
    for k in range(rec["ub"], 0, -1):
        r = solve(n, k, deadline_ms)
        lvl = {"k": k, "status": r["status"], "nodes": r["nodes"]}
        if r["status"] == "SAT":
            lvl["witness"] = r["witness"]
            rec["k_max"] = k
            rec["witness"] = r["witness"]
        rec["levels"].append(lvl)
        if r["status"] in ("SAT", "TIMEOUT"):
            break
    rec["certified"] = (rec["k_max"] is not None
                        and all(l["status"] != "TIMEOUT" for l in rec["levels"]))
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="lo", type=int, default=2)
    ap.add_argument("--to", dest="hi", type=int, default=40)
    ap.add_argument("--deadline-per-level-ms", type=int, default=120000)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    results = []
    for n in range(max(2, args.lo), args.hi + 1):
        rec = solve_n(n, args.deadline_per_level_ms)
        results.append(rec)
        lv = " ".join(f"{l['k']}:{l['status']}[{l['nodes']}n]" for l in rec["levels"])
        print(f"n={n} ub={rec['ub']} k_max={rec['k_max']} "
              f"certified={rec['certified']} | {lv}", flush=True)
        if args.out:
            with open(args.out, "w") as f:
                json.dump(results, f, indent=1)
    print("DONE", len(results), "cases")


if __name__ == "__main__":
    main()
