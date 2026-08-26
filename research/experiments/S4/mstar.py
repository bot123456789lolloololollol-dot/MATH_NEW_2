#!/usr/bin/env python3
"""mstar.py -- exact computation of M*(sigma), the smallest positive integer
whose decimal digit-square-sum equals sigma.

Theory (full write-up in REPORT.md / SPEC.md):
  * A d-digit number n = sum c_i 10^i has ssd(n) = sum c_i^2 <= 81 d.
    Define the DEFICIT delta(d, sigma) = 81 d - sigma = sum (81 - c_i^2).
    Per-digit deficits v_c = 81 - c^2 form V = {81,80,77,72,65,56,45,32,17,0}.
    So sigma is achievable with exactly d digits  <=>  delta >= 0 and delta is
    a sum of AT MOST d elements of V (zeros allowed for padding digits 9...9
    wait -- zeros in V correspond to digit 9; padding with digit 9 adds
    deficit 0 and one digit, hence "at most d nonzero-deficit terms").
  * d*(sigma) = min{ d : feasible }, and any shorter number beats any longer
    one numerically, so M*(sigma) has exactly d*(sigma) digits.
  * Among d-digit strings with ssd = sigma, the lexicographically smallest is
    the numerically smallest; it is built greedily left-to-right choosing the
    smallest next digit whose residual is still feasible.
  * Greedy termination: as soon as residual == 81 * remaining_digits, ALL
    remaining digits are forced to be 9, so the suffix is nines and the whole
    number is compactly (prefix string, number of trailing nines).
"""
from __future__ import annotations

V_DEFICITS = sorted({81 - c * c for c in range(10)})          # [0,17,...,81]
SQ_DIGIT = [c * c for c in range(10)]

# ---------------------------------------------------------------- sumsets ---
_DELTA_CAP = 4000        # table covers deltas 0..CAP
_J_MAX = 48              # table covers sums of <= J_MAX elements of V

_D_TABLE: list[set[int]] = [set() for _ in range(_J_MAX + 1)]
_D_TABLE[0].add(0)
for _j in range(1, _J_MAX + 1):
    prev = _D_TABLE[_j - 1]
    cur = _D_TABLE[_j]
    for x in prev:
        for v in V_DEFICITS:
            y = x + v
            if y <= _DELTA_CAP:
                cur.add(y)


def _min_terms(delta: int) -> int | None:
    """Smallest j <= _J_MAX with delta in D_j, else None."""
    if delta > _DELTA_CAP:
        raise ValueError("delta out of table range")
    for j in range(_J_MAX + 1):
        if delta in _D_TABLE[j]:
            return j
    return None


def feasible(d: int, sigma: int) -> bool:
    """Can exactly d digits (first digit >= 1 irrelevant here) sum to sigma?"""
    if d <= 0:
        return sigma == 0
    delta = 81 * d - sigma
    if delta < 0:
        return False
    if delta > _DELTA_CAP:
        raise ValueError("delta out of range")
    j = _min_terms(delta)
    return j is not None and j <= d


def d_star(sigma: int, extra_scan: int = 10) -> int:
    """Minimal digit count of a number with ssd == sigma."""
    base = -(-sigma // 81)  # ceil(sigma/81)
    for t in range(extra_scan + 1):
        if feasible(base + t, sigma):
            return base + t
    raise AssertionError(f"no d* found within +{extra_scan} for sigma={sigma}")


def mstar_compact(sigma: int) -> tuple[str, int]:
    """Return (prefix_digits, nines_count): number = prefix + '9'*nines."""
    if sigma <= 0:
        raise ValueError("sigma must be >= 1")
    d = d_star(sigma)
    digits: list[str] = []
    rem = sigma
    for i in range(d):
        if rem == 81 * (d - i):               # all remaining digits forced to 9
            return "".join(digits), d - i
        r_after = d - i - 1
        if rem == 81 * r_after:               # rest is all nines, forced
            return "".join(digits), r_after
        lo = 1 if i == 0 else 0
        for c in range(lo, 10):
            if feasible(r_after, rem - SQ_DIGIT[c]):
                digits.append(str(c))
                rem -= SQ_DIGIT[c]
                break
        else:
            raise AssertionError("greedy stuck (should never happen)")
    assert rem == 0
    return "".join(digits), 0


def mstar_int(sigma: int) -> int:
    """Full integer value (only for small results)."""
    p, nines = mstar_compact(sigma)
    return int(p + "9" * nines) if nines else int(p)


# ------------------------------------------------------- brute validation ---
def brute_mstar(sigma: int, limit: int = 400000) -> int | None:
    m = 1
    while m <= limit:
        if ssd_fast(m) == sigma:
            return m
        m += 1
    return None


def ssd_fast(n: int) -> int:
    s = 0
    while n:
        q, r = divmod(n, 10)
        s += r * r
        n = q
    return s


def selftest(sigma_limit: int = 350, verbose: bool = True) -> dict:
    """Cross-validate mstar against literal brute force for small sigma."""
    bad = []
    checked = 0
    for sigma in range(1, sigma_limit + 1):
        got = mstar_int(sigma)
        want = brute_mstar(sigma)
        checked += 1
        if got != want:
            bad.append((sigma, got, want))
    if verbose:
        print(f"mstar selftest: {checked} sigmas checked, "
              f"{len(bad)} mismatches")
    if bad:
        for row in bad[:10]:
            print("  MISMATCH sigma,got,want =", row)
        raise SystemExit("mstar selftest FAILED")
    return {"checked": checked, "mismatches": 0}


if __name__ == "__main__":
    selftest()
    print("examples:")
    for s in (15999, 16031, 42, 120):
        p, n = mstar_compact(s)
        print(f"  sigma={s}: prefix={p!r} nines={n} d={len(p)+n}")
