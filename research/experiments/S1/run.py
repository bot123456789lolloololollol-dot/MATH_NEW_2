"""
run.py — Orchestrator for the S1 certified exact-table experiment.

Phases:
  solve    : node solver.js over [LO..HI] (chunked by caller), raw JSON out
  verify   : (a) independent checker validates every witness
             (b) compare k_max vs OEIS A260999 b-file (third_party/b260999.txt)
             (c) audit Cariboni witness sets (third_party/w260999.txt)
                 under the strict checker for all n <= HI
             (d) recompute inverse terms a(k) = min{n : f(n)>=k} and compare
                 with OEIS A004136 prefix
Usage:
  python run.py solve --lo 2 --hi 90 --tag chunkA [--deadline-ms 300000]
  python run.py verify --hi 120 --results outputs/main_chunkA.json,...
Exit codes: 0 ok; 1 mismatches found; 2 usage.
"""
import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_sidon  # noqa: E402


def phase_solve(lo, hi, tag, deadline_ms):
    out = os.path.join(HERE, "outputs", f"main_{tag}.json")
    cmd = ["node", os.path.join(HERE, "solver.js"),
           "--from", str(lo), "--to", str(hi),
           "--deadline-per-level-ms", str(deadline_ms),
           "--out", out]
    print("RUN:", " ".join(cmd), flush=True)
    r = subprocess.run(cmd, timeout=560)
    return r.returncode


def load_b_file(path):
    vals = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            vals[int(parts[0])] = int(parts[1])
    return vals


def load_cariboni_witnesses(path):
    """Parse 'a(n)=k from {0,1,3,...}' style lines."""
    import re
    out = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = re.match(r"a\((\d+)\)\s*=\s*\d+\s*from\s*\{(.*)\}", line)
            if not m:
                continue
            n = int(m.group(1))
            S = sorted({int(t.strip()) % n for t in m.group(2).split(",")})
            out[n] = S
    return out


def phase_verify(hi, result_files):
    merged = {}
    for rf in result_files:
        path = os.path.join(HERE, "outputs", rf)
        for rec in json.load(open(path)):
            merged[rec["n"]] = rec
    ns = sorted(n for n in merged if 2 <= n <= hi)
    report = {"verified_ns": len(ns), "uncertified": [], "witness_fails": [],
              "bfile_mismatches": [], "bfile_missing": [],
              "cariboni_audit": {"checked": 0, "fails": [], "size_mismatch": []},
              "inverse_a004136": {"recomputed": {}, "oeis": [1, 3, 7, 13, 21, 31,
                                                             48, 57, 73, 91, 120,
                                                             133, 168, 183, 255,
                                                             255, 273, 307],
                                  "mismatches": []}}

    # (a)+(b): checker on every witness + compare with b-file
    cases = []
    for n in ns:
        rec = merged[n]
        if not rec["certified"]:
            report["uncertified"].append(n)
            continue
        cases.append({"n": n, "set": rec["witness"], "expected_k": rec["k_max"]})
    npass = 0
    for case in cases:
        st, det = check_sidon.verify_case(case)
        if st != "PASS":
            report["witness_fails"].append({"n": case["n"], "status": st,
                                            "detail": det})
        else:
            npass += 1
    report["checker_passed"] = npass

    bvals = load_b_file(os.path.join(HERE, "third_party", "b260999.txt"))
    for n in ns:
        rec = merged[n]
        if not rec["certified"]:
            continue
        if n in bvals and bvals[n] != rec["k_max"]:
            report["bfile_mismatches"].append({"n": n, "ours": rec["k_max"],
                                               "bfile": bvals[n]})
        if n not in bvals:
            report["bfile_missing"].append(n)

    # (c): strict-checker audit of Cariboni's published witnesses
    try:
        wit = load_cariboni_witnesses(
            os.path.join(HERE, "third_party", "w260999.txt"))
        for n, S in sorted(wit.items()):
            if n < 2 or n > hi:
                continue
            report["cariboni_audit"]["checked"] += 1
            st, det = check_sidon.verify_case({"n": n, "set": S,
                                               "expected_k": bvals.get(n)})
            if st != "PASS":
                report["cariboni_audit"]["fails"].append(
                    {"n": n, "status": st, "detail": det})
            elif n in bvals and len(S) != bvals[n]:
                report["cariboni_audit"]["size_mismatch"].append(
                    {"n": n, "set_size": len(S), "bfile": bvals[n]})
    except FileNotFoundError:
        report["cariboni_audit"]["note"] = "witness file missing"

    # (d): inverse sequence a(k)=min{n: f(n)>=k}; valid because [2..hi]
    # fully certified contiguously
    f_of_n = {n: merged[n]["k_max"] for n in ns if merged[n]["certified"]}
    ak = {}
    for n in sorted(f_of_n):
        for k in range(1, f_of_n[n] + 1):
            if k not in ak:
                ak[k] = n
    report["inverse_a004136"]["recomputed"] = ak
    for k, v in ak.items():
        idx = k - 1
        if idx < len(report["inverse_a004136"]["oeis"]) and \
                report["inverse_a004136"]["oeis"][idx] != v:
            report["inverse_a004136"]["mismatches"].append(
                {"k": k, "ours": v, "oeis": report["inverse_a004136"]["oeis"][idx]})

    out = os.path.join(HERE, "outputs", f"verification_upto_{hi}.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=1)

    clean = (not report["uncertified"] and not report["witness_fails"]
             and not report["bfile_mismatches"]
             and not report["cariboni_audit"]["fails"]
             and not report["cariboni_audit"]["size_mismatch"]
             and not report["inverse_a004136"]["mismatches"])
    print(json.dumps({k: (v if k != "cariboni_audit" else
                          {kk: vv for kk, vv in v.items()}) for k, v in
                      report.items()}, indent=1)[:3000])
    print(f"\nVERIFY {'CLEAN' if clean else 'ISSUES FOUND'} -> {out}")
    return 0 if clean else 1


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="phase", required=True)
    s = sub.add_parser("solve")
    s.add_argument("--lo", type=int, required=True)
    s.add_argument("--hi", type=int, required=True)
    s.add_argument("--tag", required=True)
    s.add_argument("--deadline-ms", type=int, default=300000)
    v = sub.add_parser("verify")
    v.add_argument("--hi", type=int, required=True)
    v.add_argument("--results", nargs="+", required=True)
    args = ap.parse_args()
    if args.phase == "solve":
        sys.exit(phase_solve(args.lo, args.hi, args.tag, args.deadline_ms))
    sys.exit(phase_verify(args.hi, args.results))


if __name__ == "__main__":
    main()
