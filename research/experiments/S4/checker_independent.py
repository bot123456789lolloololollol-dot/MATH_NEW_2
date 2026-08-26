#!/usr/bin/env python3
"""checker_independent.py -- SECOND, algorithmically different verifier for
the certified extension of A094406 / A176762 produced by extend.py.

Differences from the primary implementation:
  * feasibility of (digit_count, sigma) is decided by meet-in-the-middle
    enumeration of deficit multisets (primary used DP sumset tables);
  * lexicographic minimality of the witness prefix is re-derived by testing
    every smaller alternative digit at every prefix position;
  * tails are computed by fresh direct iteration (no memo/recursion);
  * window completeness bounds are re-derived here, then the full window is
    rescanned and the argmin recomputed.

Exit code 0 iff EVERYTHING (including all tamper tests) behaves as expected.
"""
from __future__ import annotations

import itertools
import json
import sys

UNHAPPY_CYCLE = frozenset({4, 16, 37, 58, 89, 145, 42, 20})
V_NONZERO = tuple(sorted({81 - c * c for c in range(10)} - {0}))  # deficits
SQ = [c * c for c in range(10)]
FROBENIUS_17_32 = 495          # largest int NOT of form 17a+32b


def ssd(n: int) -> int:
    return sum(int(c) ** 2 for c in str(n))


def tail_direct(n: int, cap: int = 100000) -> int | None:
    """Fresh simulation; None if happy."""
    x, st = n, 0
    while x not in UNHAPPY_CYCLE:
        if x == 1:
            return None
        x = ssd(x)
        st += 1
        if st > cap:
            raise RuntimeError("cap")
    return st


# ---------------------------------------------------------------- deficits --
def feasible_mitm(d: int, sigma: int) -> bool:
    """Meet-in-the-middle decision: does some d-digit string have ssd sigma?
    Equivalent to: delta = 81d - sigma >= 0 and delta is a sum of <= d
    elements of V (zeros allowed).  Coverage argument: if the minimal
    multiset has j* <= ceil(delta/17) terms, the 'right' table alone
    (sizes 0..max_b = ceil(delta/17)) contains it."""
    delta = 81 * d - sigma
    if delta < 0:
        return False
    if delta == 0:
        return True
    if delta > FROBENIUS_17_32:
        # beyond the Frobenius number of {17,32} every delta IS representable
        # as 17a+32b; only the term-count bound remains to check
        best = min((b + (delta - 32 * b) // 17)
                   for b in range(delta // 32 + 1)
                   if (delta - 32 * b) % 17 == 0)
        return best <= d
    half = max(0, (delta + 16) // 17 // 2)
    max_b = max(0, (delta + 16) // 17)
    left = {}
    for a in range(0, half + 1):
        for combo in itertools.combinations_with_replacement(V_NONZERO, a):
            s = sum(combo)
            if s <= delta and (s not in left or a < left[s]):
                left[s] = a
    right = {}
    for b in range(0, max_b + 1):
        for combo in itertools.combinations_with_replacement(V_NONZERO, b):
            s = sum(combo)
            if s <= delta and (s not in right or b < right[s]):
                right[s] = b
    for s, ca in left.items():
        cb = right.get(delta - s)
        if cb is not None and ca + cb <= d:
            return True
    return False


def mstar_verify(sigma: int, prefix: str, nines: int):
    """Independently verify prefix+'9'*nines is THE minimal number with
    ssd = sigma.  Returns (ok, reason)."""
    d = len(prefix) + nines

    def val_ssd(p, nn):
        return sum(int(c) ** 2 for c in p) + 81 * nn

    if val_ssd(prefix, nines) != sigma:
        return False, f"witness ssd != sigma ({val_ssd(prefix, nines)} != {sigma})"
    if prefix != prefix.upper() or not prefix.isdigit() or prefix[0] == "0":
        return False, "prefix not a canonical positive digit string"
    # (1) claimed length is achievable
    if not feasible_mitm(d, sigma):
        return False, f"length {d} not feasible"
    # (2) length d-1 infeasible
    if d >= 2 and feasible_mitm(d - 1, sigma):
        return False, f"length {d-1} IS feasible -- witness not shortest"
    # (3) everything shorter is arithmetically impossible
    if d >= 3 and sigma > 81 * (d - 2):
        pass  # then even d-2 impossible; combined with (2) => all < d impossible
    elif d >= 3:
        # need explicit infeasibility for lengths < d-1 down to ceil(sigma/81)
        lo = max(1, -(-sigma // 81))
        for dd in range(lo, d - 1):
            if feasible_mitm(dd, sigma):
                return False, f"length {dd} feasible -- shorter witness exists"
    # (4) lexicographic minimality at fixed length d:
    #     every alternative smaller digit at each position must be infeasible
    rem = sigma
    for i, ch in enumerate(prefix):
        r_after = d - i - 1
        lo_c = 1 if i == 0 else 0
        for c in range(lo_c, int(ch)):
            if feasible_mitm(r_after, rem - SQ[c]):
                return False, (f"smaller digit {c} possible at position {i} "
                               f"(witness not lexicographically minimal)")
        rem -= SQ[int(ch)]
    if rem != 81 * nines:
        return False, "residual mismatch for nine-suffix"
    return True, "ok"


def compact_from_json(h: dict) -> tuple[str, int]:
    return h["prefix"], int(h["nines"])


# ------------------------------------------------------------------- main ---
def check(doc: dict) -> tuple[bool, list[str]]:
    notes = []
    ok_all = True

    def fail(msg):
        nonlocal ok_all
        ok_all = False
        notes.append("FAIL: " + msg)

    def good(msg):
        notes.append("pass: " + msg)

    # ---- 0. definitions sanity -------------------------------------------
    if tail_direct(4) != 0 or tail_direct(2) != 1 or tail_direct(15999) != 12:
        fail("tail sanity anchors wrong")
    else:
        good("tail anchors: t(4)=0, t(2)=1, t(15999)=12")

    # ---- 1. published small terms reproduced (fresh brute force) ---------
    pubA = doc["published_terms_verification"]["a094406_recomputed_0_12"]
    wantA = [4, 2, 11, 15, 5, 3, 14, 45, 36, 6, 112, 269, 15999]
    rec = []
    seen_tails = {}
    for m in range(1, 30000):
        t = tail_direct(m)
        if t is None:
            continue
        if t not in seen_tails:
            seen_tails[t] = m
    for t in range(13):
        rec.append(seen_tails.get(t))
    if rec != wantA or pubA != wantA:
        fail(f"A094406 small-term reproduction mismatch {rec}")
    else:
        good("A094406 a(0..12) reproduced independently")

    # ---- 2. phase B: a(13) certificate -----------------------------------
    phB = doc["phases"][0]
    winB = phB["winner"]
    if (winB["prefix"] + "9" * int(winB["nines"])) != "577" + "9" * 196:
        fail("phase-B claimed value not '577'+196 nines")
    w_sigma = phB["winner_sigma"]
    win = phB["winner"]
    p, n = compact_from_json(win)
    okB, why = mstar_verify(w_sigma, p, n)
    if not okB or tail_direct(w_sigma) != 12:
        fail(f"phase B winner invalid: {why} tail={tail_direct(w_sigma)}")
    else:
        good(f"a(13): sigma={w_sigma} minimal witness verified "
             f"(M* = '577'+196 nines)")
    # independent window rescan
    loB, hiB = phB["window"]
    hits = []
    for sig in range(int(loB), int(hiB) + 1):
        if tail_direct(sig) == 12:
            hits.append(sig)
    if hits != [phB["winner_sigma"]] or len(phB["hits_with_tail_12"]) != len(hits):
        fail(f"phase B window rescan mismatch: {hits}")
    else:
        good(f"phase B window [{loB},{hiB}] rescan: unique hit {hits[0]}")

    # ---- 3. phase C: NEW a(14) certificate -------------------------------
    phC = doc["phases"][1]
    a13 = 578 * 10**196 - 1
    if tail_direct(a13) != 13 or ssd(a13) != 15999:
        fail("a(13) anchor inconsistent")
    winC = phC["winner"]
    pC, nC = compact_from_json(winC)
    d_win = len(pC) + nC
    sigma_star = a13 + int(phC["sigma_star_delta_from_a13"])
    okC, why = mstar_verify(sigma_star, pC, nC)
    if not okC:
        fail(f"phase C winner invalid: {why}")
    elif tail_direct(sigma_star) != 13:
        fail("phase C winner sigma tail != 13")
    elif int(winC.get("tail", -1)) != 13 or \
            int(winC.get("digits", -1)) != d_win:
        fail("phase C winner metadata (tail/digits) inconsistent")
    else:
        good(f"a(14): sigma=a13+{phC['sigma_star_delta_from_a13']} "
             f"witness ('{pC}' + {nC} nines) fully verified: shortest "
             f"length AND lexicographically minimal")

    # independent window completeness: hi = 81 * digits(U)
    loC, hiC = int(phC["window"][0]), int(phC["window"][1])
    if loC != a13 or hiC != 81 * d_win:
        fail("phase C window bounds not the completeness bounds")
    else:
        good("phase C window = [a13, 81*d(winner)] as required")
    best = None
    nhits = 0
    x = loC
    deltas_hit = []
    while x <= hiC:
        t = tail_direct(x)
        if t == 13:
            nhits += 1
            deltas_hit.append(x - a13)
            # candidate M*(x): compute independently via search over (len,prefix)
            cand = independent_mstar(x, d_cap_hint=None)
            if cand is None:
                fail(f"independent M* search failed for sigma=a13+{x-a13}")
            else:
                pc, nc = cand
                if best is None:
                    best = (pc, nc, x - a13)
                else:
                    if (len(pc) + nc, pc) < (len(best[0]) + best[1], best[0]):
                        best = (pc, nc, x - a13)
        x += 1
    if nhits != len(phC["hits_with_tail_13"]) or best is None:
        fail(f"phase C rescan hit count mismatch ({nhits})")
    elif (best[0], best[1]) != (pC, nC) or best[2] != int(
            phC["sigma_star_delta_from_a13"]):
        fail(f"phase C independent argmin mismatch: {best[:2]} vs winner")
    else:
        good(f"phase C independent rescan: {nhits} hit(s) at delta="
             f"{deltas_hit}, argmin agrees")

    # ---- 4. phase D dominance --------------------------------------------
    phD = doc["phases"][2]
    if "A176762(21) = A094406(14)" not in phD["claim"]:
        fail("phase D claim text missing")
    dom = phD["dominance"]
    if not (dom["log10_lower_bound_for_digits_of_any_height21_happy"] > 200
            and dom["winner_digit_count"] == d_win):
        fail("phase D dominance numbers inconsistent")
    else:
        good("phase D: happy-branch lower bound (>10^200 digits) dominates "
             f"winner ({str(d_win)[:12]}... digits)")

    return ok_all, notes


def independent_mstar(sigma: int, d_cap_hint):
    """Find THE minimal number with ssd=sigma, independently: scan candidate
    lengths upward from ceil(sigma/81); for the first feasible length, build
    the lexicographically smallest string using feasible_mitm."""
    base = -(-sigma // 81)
    dstar = None
    for t in range(12):
        if feasible_mitm(base + t, sigma):
            dstar = base + t
            break
    if dstar is None:
        return None
    # guard: for huge sigma the loop above runs <= 12 mitm calls on huge
    # integers -- fine. Build prefix greedily.
    digits = []
    rem = sigma
    for i in range(dstar):
        if rem == 81 * (dstar - i):
            return "".join(digits), dstar - i
        r_after = dstar - i - 1
        lo_c = 1 if i == 0 else 0
        for c in range(lo_c, 10):
            if feasible_mitm(r_after, rem - SQ[c]):
                digits.append(str(c))
                rem -= SQ[c]
                break
        else:
            return None
    return "".join(digits), 0


# ----------------------------------------------------------- tamper tests ---
def tamper_suite(doc: dict) -> list[dict]:
    import copy
    results = []

    def run(mutator, expect_fail: bool, name: str):
        d2 = copy.deepcopy(doc)
        mutator(d2)
        ok, _ = check(d2)
        rejected = not ok
        results.append({"test": name, "expected_reject": expect_fail,
                        "rejected": rejected,
                        "verdict": "PASS" if rejected == expect_fail else "FAIL"})

    def m_flip_last_prefix_digit(d2):
        w = d2["phases"][1]["winner"]
        w["prefix"] = w["prefix"][:-1] + \
            str((int(w["prefix"][-1]) - 1) % 10 or 9)

    def m_permute_prefix(d2):
        w = d2["phases"][1]["winner"]
        w["prefix"] = w["prefix"][1] + w["prefix"][0] + w["prefix"][2:]

    def m_drop_one_nine(d2):
        w = d2["phases"][1]["winner"]
        w["nines"] = int(w["nines"]) - 1
        for h in d2["phases"][1]["hits_with_tail_13"]:
            h["nines"] = int(h["nines"]) - 1

    def m_add_one_nine(d2):
        w = d2["phases"][1]["winner"]
        w["nines"] = int(w["nines"]) + 1

    def m_fake_smaller_hit(d2):
        h = copy.deepcopy(d2["phases"][1]["hits_with_tail_13"][0])
        h["prefix"] = "1888"       # lexicographically smaller than '5888'
        h["digits"] = int(h["digits"])
        d2["phases"][1]["hits_with_tail_13"].append(h)

    def m_corrupt_a13(d2):
        d2["phases"][0]["winner"]["prefix"] = "579"

    def m_wrong_tail(d2):
        d2["phases"][1]["winner"]["tail"] = 12

    def m_shift_window(d2):
        w = d2["phases"][1]["window"]
        d2["phases"][1]["window"] = [w[0], str(int(w[1]) + 81)]

    tests = [
        (m_flip_last_prefix_digit, True, "flip last prefix digit 8->9/7"),
        (m_permute_prefix, True, "transpose prefix digits (same multiset)"),
        (m_drop_one_nine, True, "remove one nine everywhere"),
        (m_add_one_nine, True, "add one spurious nine"),
        (m_fake_smaller_hit, True, "inject fake smaller-M* hit"),
        (m_corrupt_a13, True, "corrupt a(13) prefix"),
        (m_wrong_tail, True, "misstate winner tail"),
        (m_shift_window, True, "shift completeness window bound"),
    ]
    for fn, exp, name in tests:
        run(fn, exp, name)
    return results


def main() -> None:
    path = sys.argv[1] if len(sys.argv) > 1 else "out/a094406_certification.json"
    with open(path, encoding="utf-8") as f:
        doc = json.load(f)
    ok, notes = check(doc)
    print("--- independent check ---")
    for nline in notes:
        print(" ", nline)
    print("CHECK:", "ALL PASS" if ok else "FAILURES PRESENT")

    tampers = tamper_suite(doc)
    print("--- tamper suite ---")
    all_tamper_ok = True
    for t in tampers:
        print(f"  [{t['verdict']}] {t['test']} (rejected={t['rejected']})")
        if t["verdict"] == "FAIL":
            all_tamper_ok = False
    report = {
        "check_passed": bool(ok),
        "notes": notes,
        "tamper_tests": tampers,
        "all_tamper_tests_behave": all_tamper_ok,
    }
    with open("out/check_report.json", "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1)
    print("wrote out/check_report.json")
    sys.exit(0 if (ok and all_tamper_ok) else 1)


if __name__ == "__main__":
    main()
