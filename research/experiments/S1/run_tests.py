"""
run_tests.py — Test battery for the S1 Sidon/MGR exact-value pipeline.

Test 1  checker corruption selftest (check_sidon.py --selftest logic).
Test 2  THIRD independent implementation: naive exhaustive search over all
        C(n,k) subsets (itertools) for small n; establishes ground-truth
        k_max without any symmetry assumptions at all. Compare with both
        solvers' outputs for n in [2..16].
Test 3  cross-implementation agreement: solver.py vs solver.js on [2..N1]
        (values AND witnesses up to translation/negation-free equality of
        canonical form: both solvers use identical deterministic order, so
        raw witness lists must match exactly).
Outputs a PASS/FAIL summary; exit nonzero on any failure.

Usage:
  python run_tests.py [--agree-to 36] [--naive-to 16] [--js-deadline-ms 60000]
Raw output also written to outputs/tests_raw.json.
"""
import argparse
import itertools
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_sidon            # noqa: E402
import solver as pysolver     # noqa: E402


def naive_kmax(n):
    """Ground truth by full enumeration: max k such that some k-subset of Z_n
    is Sidon (strict). No symmetry reductions whatsoever. Monotonicity used
    to stop: subsets of Sidon sets are Sidon, so the first failing k is final."""
    best = 0
    universe = range(n)
    for k in range(1, n + 1):
        found = False
        for combo in itertools.combinations(universe, k):
            residues = list(combo)
            ok1, _ = check_sidon.check_f1_sums(n, residues)
            if ok1:
                found = True
                break
        if found:
            best = k
        else:
            # larger k only harder once one size fails? NOT guaranteed monotone
            # a priori, but Sidon-ness is hereditary (subsets of Sidon are
            # Sidon), so first failing k means no larger works.
            break
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--agree-to", type=int, default=36)
    ap.add_argument("--naive-to", type=int, default=16)
    ap.add_argument("--js-deadline-ms", type=int, default=60000)
    args = ap.parse_args()

    os.makedirs(os.path.join(HERE, "outputs"), exist_ok=True)
    report = {"tests": []}
    failures = []

    # ---- Test 1: checker selftest ---------------------------------------
    bad = check_sidon.selftest()
    t1 = {"name": "checker_corruption_selftest", "ok": not bad, "detail": bad}
    report["tests"].append(t1)
    print(("PASS" if t1["ok"] else "FAIL"), "- checker corruption selftest")
    failures += bad

    # ---- Test 2: naive ground truth vs python solver --------------------
    mism = []
    for n in range(2, args.naive_to + 1):
        gt = naive_kmax(n)
        rec = pysolver.solve_n(n, deadline_ms=120000)
        if not rec["certified"] or rec["k_max"] != gt:
            mism.append({"n": n, "naive": gt, "solver_py": rec["k_max"],
                         "certified": rec["certified"]})
        else:
            st, det = check_sidon.verify_case({"n": n, "set": rec["witness"],
                                               "expected_k": gt})
            if st != "PASS":
                mism.append({"n": n, "checker_rejected_witness": det})
    t2 = {"name": "naive_vs_pysolver", "range": f"2..{args.naive_to}",
          "ok": not mism, "mismatches": mism}
    report["tests"].append(t2)
    print(("PASS" if t2["ok"] else "FAIL"),
          f"- naive enumeration vs python solver on 2..{args.naive_to}")
    failures += [str(m) for m in mism]

    # ---- Test 3: python vs node agreement --------------------------------
    out_json = os.path.join(HERE, "outputs", "agreement_js.json")
    proc = subprocess.run(
        ["node", os.path.join(HERE, "solver.js"),
         "--from", "2", "--to", str(args.agree_to),
         "--deadline-per-level-ms", str(args.js_deadline_ms),
         "--out", out_json],
        capture_output=True, text=True, timeout=540)
    js_ok = proc.returncode == 0
    js = {r["n"]: r for r in json.load(open(out_json))}
    mism3 = []
    for n in range(2, args.agree_to + 1):
        pyrec = pysolver.solve_n(n, deadline_ms=args.js_deadline_ms)
        jrec = js.get(n)
        if not jrec or not pyrec["certified"] or not jrec["certified"]:
            mism3.append({"n": n, "why": "uncertified/missing"})
            continue
        if pyrec["k_max"] != jrec["k_max"]:
            mism3.append({"n": n, "py": pyrec["k_max"], "js": jrec["k_max"]})
        elif sorted(pyrec["witness"]) != sorted(jrec["witness"]):
            mism3.append({"n": n, "why": "different witnesses (unexpected: "
                                        "identical deterministic search)"})
    t3 = {"name": "py_vs_node_agreement", "range": f"2..{args.agree_to}",
          "node_ran": js_ok, "ok": js_ok and not mism3, "mismatches": mism3}
    report["tests"].append(t3)
    print(("PASS" if t3["ok"] else "FAIL"),
          f"- python vs node agreement on 2..{args.agree_to}")
    failures += [str(m) for m in mism3]

    with open(os.path.join(HERE, "outputs", "tests_raw.json"), "w") as f:
        json.dump(report, f, indent=1)

    if failures:
        print("OVERALL: FAIL")
        sys.exit(1)
    print("OVERALL: PASS")


if __name__ == "__main__":
    main()
