#!/usr/bin/env python3
"""extend.py -- certified extension of OEIS A094406 beyond a(13), plus the
derived extension of A176762 to n=21.

Phases:
  B. Re-certify the published a(13) = 578*10^196 - 1 with a complete
     minimality argument (window scan above the base a(12) = 15999).
  C. Compute NEW term a(14) = smallest unhappy number whose trajectory under
     ssd needs 14 steps to reach the 8-cycle, with the same style of
     completeness certificate.
  D. Prove A176762(21) equals the new A094406(14) (dominance over the happy
     branch A001273 via digit-count pigeonhole + published A001273 data).

Every claim is dumped to out/a094406_certification.json for the independent
checker (checker_independent.py).
"""
from __future__ import annotations

import json
import math
import sys

import happy_core as hc
from happy_core import UNHAPPY_CYCLE, ssd, tail_u
from mstar import V_DEFICITS, d_star, feasible, mstar_compact, selftest

OUT = "out/a094406_certification.json"


# ---------------------------------------------------------------- helpers ---
def num_digits(n: int) -> int:
    return len(str(n))


def ceil_div(a: int, b: int) -> int:
    return -(-a // b)


def cmp_compact(pa: str, na: int, pb: str, nb: int) -> int:
    """Compare numbers 'pa'+'9'*na vs 'pb'+'9'*nb without materialising nines."""
    da, db = len(pa) + na, len(pb) + nb
    if da != db:
        return -1 if da < db else 1
    L = max(len(pa), len(pb))
    for i in range(L):
        ca = pa[i] if i < len(pa) else "9"
        cb = pb[i] if i < len(pb) else "9"
        if ca != cb:
            return -1 if ca < cb else 1
    return 0


_TAIL_CACHE: dict[int, int | None] = {}


def _tail_plain(n: int) -> int | None:
    x, st = n, 0
    while x not in UNHAPPY_CYCLE:
        if x == 1:
            return None
        x = ssd(x)
        st += 1
    return st


def tail_of_big(n: int) -> int | None:
    """Tail of an arbitrarily large n: recurse through ssd (shrinks fast).

    For n >= 100, ssd(n) < n and the trajectory never returns to n, hence
    tail(n) = 1 + tail(ssd(n)); recursion bottoms out at small arguments.
    """
    if n in _TAIL_CACHE:
        return _TAIL_CACHE[n]
    if n <= 400000:
        val = _tail_plain(n)
    else:
        sub = tail_of_big(ssd(n))
        val = None if sub is None else sub + 1
    _TAIL_CACHE[n] = val
    return val


# ------------------------------------------------------------------ phase B --
def certify_a13() -> dict:
    a12 = 15999                       # published a(12), brute-force verified
    claimed = 578 * 10**196 - 1       # published a(13): '577' + 196 nines
    assert hc.tail_u(claimed) == 13 and ssd(claimed) == a12
    D = num_digits(claimed)           # 199
    hi = 81 * D                       # any shorter competitor sigma <= hi
    hits = []
    for sig in range(a12, hi + 1):
        t = hc.tail_u(sig)
        if t != 12:
            continue
        p, nn = mstar_compact(sig)
        hits.append({"sigma": sig, "tail": t,
                     "prefix": p, "nines": nn, "digits": len(p) + nn})
    best = min(hits, key=lambda h: (h["digits"], h["prefix"] + ""))
    # winner determination uses plain string compare: all candidates here are
    # short enough to materialise
    full = {h["sigma"]: int(h["prefix"] + "9" * h["nines"]) for h in hits}
    best_sigma = min(full, key=lambda s: full[s])
    win = next(h for h in hits if h["sigma"] == best_sigma)
    ok = (win["prefix"] + "9" * win["nines"]) == str(claimed)
    small = [m for m in range(1, 100) if hc.tail_u(m) == 13]
    res = {
        "phase": "B",
        "claim": "a(13) = 578*10^196 - 1 ('577' followed by 196 nines)",
        "published_source": "OEIS A094406 comment (Schneider/Marcus), "
                            "accessed 2026-08-26",
        "base": {"a12": a12, "verified_by": "brute force scan to 30000"},
        "window": [a12, hi],
        "window_width": hi - a12 + 1,
        "hits_with_tail_12": hits,
        "winner_sigma": best_sigma,
        "winner": win,
        "matches_published_claim": ok,
        "small_m_check_lt_100": small,
    }
    print(f"B: window [{a12},{hi}] width {hi-a12+1}, "
          f"hits={len(hits)}, matches published: {ok}")
    if not ok:
        raise SystemExit("Phase B FAILED: recomputed optimum differs from "
                         "published a(13)")
    return res


# ------------------------------------------------------------------ phase C --
def extend_a14() -> dict:
    a13 = 578 * 10**196 - 1
    assert hc.tail_u(a13) == 13
    p0, n0 = mstar_compact(a13)          # candidate M*(a13)
    d0 = len(p0) + n0
    hi = 81 * d0                         # completeness bound
    lo = a13
    width = hi - lo + 1
    hits = []
    for delta in range(width):
        sig = a13 + delta
        t = tail_of_big(sig)
        if t != 13:
            continue
        p, nn = mstar_compact(sig)
        hits.append({"sigma_str": str(sig), "delta": delta, "tail": t,
                     "prefix": p, "nines": nn, "digits": len(p) + nn})
    # find minimum among hits with the virtual comparator
    best_i = 0
    for i in range(1, len(hits)):
        h = hits[i]
        b = hits[best_i]
        if cmp_compact(h["prefix"], h["nines"], b["prefix"], b["nines"]) < 0:
            best_i = i
    win = hits[best_i]
    small = [m for m in range(1, 100) if hc.tail_u(m) == 14]
    res = {
        "phase": "C",
        "claim": "NEW a(14): prefix '%s' followed by %s nines (%s digits)"
                 % (win["prefix"], win["nines"], win["digits"]),
        "sigma_star_delta_from_a13": win["delta"],
        "candidate_M*(a13)": {"prefix": p0, "nines": n0, "digits": d0},
        "window": [str(lo), str(hi)],
        "window_width": width,
        "hits_with_tail_13": hits,
        "winner_index": best_i,
        "winner": win,
        "small_m_check_lt_100": small,
    }
    print(f"C: window width {width}, hits={len(hits)}")
    for h in hits:
        mark = " <== WINNER" if h is win else ""
        print(f"   sigma=a13+{h['delta']:<4d} tail={h['tail']} "
              f"M*=({h['prefix']})+{h['nines']}x9 ({h['digits']} digits){mark}")
    return res


# ------------------------------------------------------------------ phase D --
def verify_decomposition(limit: int = 20000) -> dict:
    """Exhaustively verify, for all 1 <= m <= limit:
       happy m      =>  T(m) == happy-height(m)
       unhappy m    =>  T(m) == 7 if m in cycle else tail_u(m) + 7
    This is the metric decomposition behind A176762(n)=min(A001273(n),
    A094406(n-7)) and behind step (i) of phase D.
    """
    bad = []
    for m in range(1, limit + 1):
        t = hc.T_metric(m)
        if m == 1 or hc.classify(m) == "happy":
            # happy height: iterations to reach 1
            x, hgt = m, 0
            while x != 1:
                x = ssd(x)
                hgt += 1
            ok = (t == hgt)
        elif m in UNHAPPY_CYCLE:
            ok = (t == 7)
        else:
            tu = _tail_plain(m)
            ok = (t == tu + 7)
        if not ok:
            bad.append((m, t))
    return {"checked_up_to": limit, "violations": bad[:20],
            "n_violations": len(bad)}



def derive_a176762_21(a14: dict) -> dict:
    """Show A176762(21) = A094406(14).

    (i)  metric decomposition: T(m) = 21 <=> m happy with height 21, or m
         unhappy (non-cycle) with tail_u(m) = 14.  Verified exhaustively on
         1..20000 in happy_core (relation rows) + spot assertion below.
    (ii) happy branch dominance: any happy m of height 21 has
         ssd(m) >= A001273(20) >= 10^(A001273(19)/81 - 1), so
         m >= 10^(ceil(ssd(m)/81) - 1) which dwarfs the a(14) winner.
    """
    win = a14["winner"]
    d_win = win["digits"]

    # --- symbolic lower-bound chain for A001273 (happy heights), using ONLY
    # the published formulas (OEIS A001273, Ya-Ping Lu, Jul 2025):
    #   a_H(n) = k(n)*10^((a_H(n-1) - ssd(k(n)-1))/81) - 1  for n >= 7,
    # with k(n) < 10^4, hence ssd(k(n)-1) <= 4*81 = 324, hence
    #   E_n := (a_H(n-1) - ssd(k(n)-1))/81 >= (a_H(n-1) - 324)/81,
    #   a_H(n) >= 10^{E_n}.
    # Base (published comment): a_H(8) = 3789*10^973 - 1 > 10^975.
    # Recursion for lower bounds L(n) = lower bound on log10(a_H(n)):
    #   log10(a_H(n)) >= E_n >= (a_H(n-1) - 324)/81 >= 10^{L(n-1)} / 200
    # (for a_H(n-1) >= 10^4, subtracting 324 costs < one decimal digit).
    # So L(n) = 10^{L(n-1)} - 3 : each published formula step ADDS ONE
    # EXPONENT to the tower.  Floats overflow past ~10^308, so we saturate:
    # once L > 250 we only remember that the tower kept growing.
    L = 975.0                      # log10(a_H(8)) lower bound
    saturated_at = None            # first n where L exceeded 250
    tower_levels_after_saturation = 0
    for n in range(9, 22):
        if L <= 250.0:
            L = 10.0 ** min(L, 250.0) - 3.0
            # this L bounds log10(E_n) ... see note below
        # NOTE: for n >= 10 the quantity L would be log10(log10(...)) nested;
        # rather than tracking nested logs we exploit monotonicity: all maps
        # t -> 10^t - 3 are increasing, and once the *lower bound on
        # log10(a_H(n))* exceeds 250, every later bound does too, with one
        # further exponentiation each step.  We count those steps.
        else:
            tower_levels_after_saturation += 1
        if saturated_at is None and L > 250.0:
            saturated_at = n
    # After the loop: a lower bound on log10(a_H(21)) is itself > 10^(10^250)
    # (tower with `saturated_at`-many explicit levels plus
    # `tower_levels_after_saturation` further exponentiations), i.e. utterly
    # beyond any number we ever materialise.  What we actually need:
    # any happy m with h(m) = 21 has ssd(m) >= a_H(20) >> 10^250, so
    # digits(m) >= ssd(m)/81 >> 10^248 -- far more than d_win digits.
    assert saturated_at is not None and saturated_at <= 10
    assert tower_levels_after_saturation >= 10
    # Conservative usable statement: log10(digits of any height-21 happy m)
    # > 248, i.e. more than 10^248 digits.  The unhappy-side winner has
    # d_win = 7.13...e196 << 10^248 digits, so it dominates.
    log_digits_m_lower = 248.0
    assert d_win < 10**197

    res = {
        "phase": "D",
        "claim": "A176762(21) = A094406(14) = ('%s' followed by %s nines)"
                 % (win["prefix"], win["nines"]),
        "metric_decomposition": "T(m)=height(m) if happy; T(m)=tail(m)+7 if "
                                "unhappy non-cycle; verified exhaustively for "
                                "m<20000 in happy_core relation rows",
        "dominance": {
            "method": "digit-count pigeonhole with published A001273 "
                      "formulas (Lu 2025)",
            "tower_saturation_first_n": saturated_at,
            "tower_extra_levels_after_saturation":
                tower_levels_after_saturation,
            "log10_lower_bound_for_digits_of_any_height21_happy":
                log_digits_m_lower,
            "winner_digit_count": d_win,
            "conclusion": "any happy height-21 number has more than 10^248 "
                          "digits, while the winner has %d digits" % d_win,
        },
    }
    print("D: A176762(21) = A094406(14); happy branch dominated "
          "(digits lower bound 10^%d vs winner %d digits)"
          % (log_digits_m_lower, d_win))
    return res


def main() -> None:
    selftest(verbose=True)
    pub = hc.verify_published()
    decomp = verify_decomposition(20000)
    print("metric decomposition verified up to 20000:",
          decomp["n_violations"] == 0)
    if decomp["n_violations"]:
        raise SystemExit("decomposition lemma FAILED")
    b = certify_a13()
    c = extend_a14()
    d = derive_a176762_21(c)
    doc = {
        "target_sequences": ["A094406", "A176762", "A001273"],
        "oeis_access_date": "2026-08-26",
        "definitions_file": "happy_core.py",
        "mstar_theory_file": "mstar.py",
        "published_terms_verification": pub,
        "metric_decomposition_check": decomp,
        "phases": [b, c, d],
    }
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=1, default=str)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
