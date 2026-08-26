"""
check_sidon.py — Independent checker for Sidon sets in Z_n (STRICT convention).

Deliberately written against the DEFINITION, sharing no code or formulation
tricks with the solvers:

  F1 (primary, definitional): A is Sidon in Z_n iff all sums
      a_i + a_j (i <= j, doubles included) are pairwise distinct modulo n.
      This is exactly the OEIS A260999 definition ("distinct sums of any two
      elements").

  F2 (secondary, difference form): all unordered-pair circular distances
      min((b-a) mod n, (a-b) mod n) are pairwise distinct AND, when n is
      even, none equals n/2.  (F1 <=> F2: see SPEC.md Lemma 1.)

A case passes only if BOTH formulations accept. Any inconsistency between F1
and F2 is reported as ERROR (would indicate a misconception).

Usage:
  python check_sidon.py CASES.json        # {"cases":[{"n":int,"set":[ints], "expected_k":opt}]}
  python check_sidon.py --selftest        # corruption battery; exits nonzero on failure
Exit code 0 iff all cases PASS.
"""
import json
import sys


def check_f1_sums(n, residues):
    """Definitional check: sums a+b (i<=j) all distinct mod n."""
    k = len(residues)
    seen = set()
    for i in range(k):
        for j in range(i, k):
            s = (residues[i] + residues[j]) % n
            if s in seen:
                return False, f"F1 sum collision at value {s}"
            seen.add(s)
    return True, "F1 ok"


def check_f2_distances(n, residues):
    """Difference-form check: circular distances distinct, no n/2 pair."""
    seen = set()
    k = len(residues)
    for i in range(k):
        for j in range(i + 1, k):
            d = (residues[j] - residues[i]) % n
            dd = min(d, n - d)
            if n % 2 == 0 and dd == n // 2:
                return False, f"F2 half-modulus pair ({residues[i]},{residues[j]})"
            if dd in seen:
                return False, f"F2 distance collision at {dd}"
            seen.add(dd)
    return True, "F2 ok"


def verify_case(case):
    """Returns (status, detail). status in {'PASS','FAIL','ERROR'}."""
    n = case["n"]
    S = case["set"]
    detail = []
    # basic hygiene
    if n < 2:
        return "FAIL", "n<2"
    if len(S) == 0:
        return "FAIL", "empty set"
    residues = [x % n for x in S]
    if len(set(residues)) != len(residues):
        return "FAIL", "duplicate residues"
    if "expected_k" in case and case["expected_k"] != len(residues):
        return "FAIL", (f"claimed k={case['expected_k']} but set has "
                        f"{len(residues)} distinct residues")
    ok1, d1 = check_f1_sums(n, residues)
    ok2, d2 = check_f2_distances(n, residues)
    if ok1 != ok2:
        return "ERROR", f"formulation mismatch: {d1} / {d2} -- investigate!"
    if not ok1:
        return "FAIL", f"{d1}; {d2}"
    return "PASS", f"{d1}; {d2}; k={len(residues)}"


def selftest():
    """Corruption battery: every deliberately broken object must FAIL,
    every genuine object must PASS. Returns list of failures."""
    failures = []

    def expect(status, case, why):
        got, detail = verify_case(case)
        if got != status:
            failures.append(f"{why}: expected {status}, got {got} ({detail})")

    # --- genuine objects must PASS --------------------------------------
    expect("PASS", {"n": 7, "set": [0, 1, 3]}, "classic (7,3,1) singer set")
    expect("PASS", {"n": 13, "set": [0, 1, 3, 9], "expected_k": 4},
           "(13,4,1) perfect difference set")
    expect("PASS", {"n": 6, "set": [0, 1]}, "trivial 2-set in Z_6")
    expect("PASS", {"n": 16, "set": [0, 1, 4, 6], "expected_k": 4},
           "size-4 Sidon set in Z_16")

    # --- corrupted objects must FAIL ------------------------------------
    expect("FAIL", {"n": 7, "set": [0, 1, 3, 3]}, "duplicate residue")
    expect("FAIL", {"n": 5, "set": [0, 1, 2]},
           "sum collision 0+2 = 1+1 in Z_5")
    expect("FAIL", {"n": 8, "set": [0, 4]}, "half-modulus pair in Z_8")
    expect("FAIL", {"n": 12, "set": [0, 5, 10]},
           "distance collision in Z_12 (circular distance 5 twice)")
    expect("FAIL", {"n": 13, "set": [0, 1, 3], "expected_k": 4},
           "wrong claimed k")
    # weak-Sidon object that violates STRICT convention (A260998-style):
    # {0,1,2,4} in Z_6 has distinct sums of *distinct* pairs but 2+2 = 0+4.
    expect("FAIL", {"n": 6, "set": [0, 1, 2, 4]},
           "weak-only set must fail strict checker")

    # sanity: re-verify the two nontrivial PASS fixtures explicitly
    for fix in ({"n": 16, "set": [0, 1, 4, 6]}, {"n": 13, "set": [0, 1, 3, 9]}):
        st, det = verify_case(fix)
        if st != "PASS":
            failures.append(f"spot check {fix} unexpectedly {st}: {det}")

    return failures


def main():
    if "--selftest" in sys.argv:
        bad = selftest()
        if bad:
            print("SELFTEST FAILED:")
            for b in bad:
                print("  -", b)
            sys.exit(1)
        print("SELFTEST PASSED: all genuine objects PASS, all corrupted objects FAIL.")
        sys.exit(0)

    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)

    with open(sys.argv[1]) as f:
        data = json.load(f)
    npass = nfail = nerr = 0
    worst = None
    for case in data["cases"]:
        st, det = verify_case(case)
        tag = f"n={case['n']} k={len(case['set'])}"
        if st == "PASS":
            npass += 1
        elif st == "FAIL":
            nfail += 1
            print(f"FAIL {tag}: {det}")
        else:
            nerr += 1
            print(f"ERROR {tag}: {det}")
        if st != "PASS":
            worst = worst or tag
    print(f"SUMMARY: {npass} PASS, {nfail} FAIL, {nerr} ERROR "
          f"of {npass+nfail+nerr} cases")
    sys.exit(0 if (nfail + nerr) == 0 else 1)


if __name__ == "__main__":
    main()
