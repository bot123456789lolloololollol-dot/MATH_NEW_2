#!/usr/bin/env python3
"""happy_core.py -- core metrics for the sum-of-squares-of-digits (ssd) map.

Definitions (match OEIS A003132, A094406, A176762):
  ssd(n)            = sum of squares of the decimal digits of n.
  UNHAPPY_CYCLE C   = (4, 16, 37, 58, 89, 145, 42, 20).
  tail_u(n)         = for unhappy n: number of ssd-iterations until the
                      trajectory first hits a member of C.  Cycle members
                      have tail 0.  Happy numbers have no tail (None).
  T(n)              = (#distinct values of the trajectory n -> ssd(n) -> ...
                      up to (excluding) the first repeated value) - 1.
                      This is the A176762 metric ("steps to reach a cycle"):
                      fixed point 1 has T=0; a member of the 8-cycle has T=7.

Key lemma used everywhere (proved in REPORT.md): for every n >= 100,
ssd(n) < n, hence for n >= 100:
      tail_u(n) = 1 + tail_u(ssd(n))   (n unhappy)
      T(n)      = 1 + T(ssd(n))
"""
from __future__ import annotations

import sys

UNHAPPY_CYCLE = frozenset({4, 16, 37, 58, 89, 145, 42, 20})


def ssd(n: int) -> int:
    """Sum of squares of decimal digits of n (A003132)."""
    return sum(int(c) ** 2 for c in str(n))


def classify(n: int) -> str:
    """'happy' / 'unhappy' by iterating ssd until 1 or a cycle member."""
    x = n
    seen = set()
    while x != 1 and x not in UNHAPPY_CYCLE:
        if x in seen:
            raise RuntimeError("cycle detection failure")
        seen.add(x)
        x = ssd(x)
    return "happy" if x == 1 else "unhappy"


def tail_u(n: int) -> int | None:
    """Steps from unhappy n until first member of the 8-cycle; None if happy."""
    x, steps = n, 0
    seen = set()
    while x not in UNHAPPY_CYCLE:
        if x == 1:
            return None
        if x in seen:
            raise RuntimeError("runaway trajectory")
        seen.add(x)
        x = ssd(x)
        steps += 1
        if steps > 10_000_000:
            raise RuntimeError("step cap exceeded")
    return steps


def T_metric(n: int) -> int:
    """A176762 metric: (#distinct trajectory values until first repetition)-1.

    Walk n -> ssd(n) -> ... adding values to `seen`; stop when a value is
    revisited.  For a fixed point (1) or an 8-cycle member this correctly
    traverses the entire closed loop before stopping.
    """
    x, seen = n, set()
    while x not in seen:
        seen.add(x)
        x = ssd(x)
    return len(seen) - 1


def brute_force_records(metric, limit: int, max_index: int) -> list[int]:
    """Smallest m>=1 with metric(m)==n for n=0..max_index (scan to limit)."""
    best = [-1] * (max_index + 1)
    for m in range(1, limit + 1):
        v = metric(m)
        if v is None or v > max_index:
            continue
        if best[v] == -1:
            best[v] = m
    return best


# Published terms (OEIS, accessed 2026-08-26).
PUBLISHED_A176762 = [1, 10, 13, 23, 19, 7, 356, 4, 2, 11, 15, 5, 3, 14, 45,
                     36, 6, 112, 269, 15999]
PUBLISHED_A094406 = [4, 2, 11, 15, 5, 3, 14, 45, 36, 6, 112, 269, 15999]


def verify_published(verbose: bool = True) -> dict:
    """Reproduce all published small terms of A176762 and A094406."""
    out = {}
    rec_T = brute_force_records(T_metric, limit=30000, max_index=19)
    ok_T = rec_T == PUBLISHED_A176762
    out["a176762_recomputed"] = rec_T
    out["a176762_match"] = ok_T

    def tail_or_none(m: int):
        return tail_u(m)

    rec_U = brute_force_records(tail_or_none, limit=30000, max_index=12)
    ok_U = rec_U[:13] == PUBLISHED_A094406
    out["a094406_recomputed_0_12"] = rec_U[:13]
    out["a094406_match"] = ok_U

    # relation A176762(n) = min(A001273(n), A094406(n-7)) checked on data:
    rel = []
    for n in range(7, 20):
        u = PUBLISHED_A094406[n - 7] if n - 7 < len(PUBLISHED_A094406) else None
        rel.append((n, PUBLISHED_A176762[n], u,
                    (u is not None and PUBLISHED_A176762[n] == u)))
    out["relation_rows_n,a176762,a094406(n-7),equal"] = rel
    out["relation_holds_on_published_range"] = all(r[3] for r in rel)

    if verbose:
        print("A176762 published terms reproduced:", ok_T)
        print("A094406 published terms reproduced:", ok_U)
        print("min-formula identity holds on published range:",
              out["relation_holds_on_published_range"])
    if not (ok_T and ok_U):
        raise SystemExit("FAILED to reproduce published terms")
    return out


if __name__ == "__main__":
    import json
    res = verify_published()
    with open("out/published_terms.json", "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1, default=str)
    print("wrote out/published_terms.json")
