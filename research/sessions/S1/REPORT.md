# S1 Report — Certified exact maximum Sidon set sizes f(n) in Z_n, n ≤ 110

Date: 2026-08-26. Session: S1 (extremal combinatorics / difference packing).
Artifact dir: `research/experiments/S1/`. Spec: `research/sessions/S1/SPEC.md`.

## TL;DR

Exact table **f(n) = max size of a strict Sidon set in Z_n certified for every
n ∈ [2, 110]** (109 contiguous values) by deterministic exhaustive search with
written completeness arguments, an independent definition-based checker, and a
third brute-force ground-truth implementation. Every value matches the
published OEIS A260999 table (peer-reviewed source Haanpää–Huima–Östergård 2004
covers n ≤ 84; OEIS extension by Cariboni covers to 295 with undocumented proof
status). All 109 of Cariboni's witness sets pass our strict checker. No new
exact value was discovered — this is a rigorous **saturation / re-certification
result** (mission option c), extending the *documented-certified* prefix from
84 to 110 and providing the missing completeness argument for that stretch.

## Problem family (chosen after literature triage)

Maximum-size Sidon sets in cyclic groups (= modular Golomb rulers = single-
codeword cyclic (n,w,1) difference packings). Rejected alternatives: minimal
moduli a(k) (blocked past k=18 by the open projective-plane-of-order-18
question), graceful labeling tables (already complete at reachable sizes),
weak-Sidon variant (table already complete to n=254).

## Claim labels (strict)

| # | Claim | Label |
|---|---|---|
| C1 | For each n ∈ [2,110], f(n) equals the `k_max_certified` column of `outputs/summary_n2_110.csv`, witnessed by the stored set | **proven_result** (exhaustive DFS over the full solution space modulo translation only; counting-bound upper levels proven UNSAT; completeness arguments in SPEC §3) |
| C2 | Our f(n) agrees with OEIS A260999 for all n ∈ [2,110] (0 mismatches) | **experimentally_validated_result** |
| C3 | All 109 published Cariboni witness sets (n ≤ 110) are valid STRICT Sidon sets of the claimed size | **experimentally_validated_result** (checker-audited) |
| C4 | Recomputed inverse terms a(2..10) = 3,7,13,21,31,48,57,73,91 equal OEIS A004136 | **experimentally_validated_result** |
| C5 | Published exact table is "complete/gapless" up to n=295 | observation — we verified only to 110; Cariboni's 111..295 remains externally supported |

## Method summary (details: SPEC.md)

- Search: DFS over {0}∪S ⊆ Z_n, increasing order; slot-mask validity test on
  circular distances; counting prunes P1/P2; descending-level protocol
  (UNSAT proofs above the SAT answer ⇒ certificate). Only symmetry used is
  translation (WLOG 0 ∈ min-representative), so the completeness argument is
  trivial and airtight (SPEC §3).
- Two solvers: Node.js primary + Python re-expression — identical values AND
  identical witnesses on n=2..30 (and bit-identical node counts on reruns,
  e.g. n=105: 244,412,475 nodes twice).
- Independent checker (`check_sidon.py`): verifies BOTH the definitional sums
  form (F1) and circular-distance form (F2); disagreement reported as ERROR.
  Passes a 10-case corruption battery including a weak-convention trap object.
- Third ground truth: naive itertools enumeration for n ≤ 14 agrees.

## Evidence (regenerable)

```bash
bash research/experiments/S1/run.sh     # full pipeline, ~20–25 min, deterministic
```
Raw outputs committed under `research/experiments/S1/outputs/`:
`main_chunk{A,B,C}.json` (solver logs w/ per-level node counts & times),
`verification_upto_110.json`, `summary_n2_110.csv`, `tests_raw.json`,
`determinism_spot.json`. Third-party inputs: `third_party/b260999.txt`,
`third_party/w260999.txt` (OEIS, fetched 2026-08-26).

Compute: chunk A (2..90) 29 s; chunk B (91..105) 555 s; chunk C
(106..110) 160 s; total pure solver time 12.4 min (plus ~3 min tests).
Hardest proof: n=105, k=10 UNSAT, 244M nodes / 91 s.
No timeouts anywhere in [2,110]; all 109 cases certified.

## Key findings

1. f(2..110) fully certified; matches published data everywhere (C2).
   Notable structure: f is non-monotone (f(21)=5 > f(22)=4; f(91)=10 while
   f(92..106)=9), confirming per-n treatment is necessary.
2. The stretch n=85..110 previously had no publicly documented optimality
   proof; it now has one (rerunnable, with written completeness argument).
3. Ceiling identified honestly: n ≥ 111 where ub(n)=11 > f(n) requires k=11
   UNSAT proofs (projected hours-days here); A004136 extension past k=18 is
   blocked by an open problem, not by engineering.

## Novelty audit

`research/literature/S1-novelty.md` — verdict **prior_work_found** (family well
tabulated; our contribution = certification + audit, labeled accordingly).

## What is NOT proven / limitations

- We make no claim about n ≥ 111 (beyond noting published values exist there).
- C2/C3 validate published artifacts on [2,110]; they do not certify methods
  behind Cariboni's 111..295 entries.
- Haanpää et al.'s paper itself was not directly accessed (paywalled); their
  credited range comes via OEIS documentation.

## Next action suggestion

Port the DFS pruning to a proper CDCL SAT encoding or add unit-multiplication
canonicalization + canonical-augmentation rejection to push the certified
prefix past the ub=11 boundary (n=111..119), which is the first place a new
certified value could realistically appear.
