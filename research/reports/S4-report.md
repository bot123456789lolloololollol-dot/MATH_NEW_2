# S4 Session Report — Computational Number Theory (2026-08-26)

**Deliverable:** two certified sequence extensions with complete minimality proofs,
machine-checkable certificates, and an independent second-algorithm checker.

## Headline results
1. **A094406(14) — new term, proven.**
   a(14) = "5888" followed by N nines where N is the explicit 197-digit integer
   7135802469…358022 (total digit count N+4; equivalently a(14) = 5889·10^N − 1).
   Published range ended at a(13) = 578·10^196 − 1 (OEIS, rev. 2026-08-22).
2. **A094406(13) re-certified** — first written completeness proof for the published value
   (exhaustive window scan [15999, 16119] + exact deficit-multiset optimality).
3. **A176762(21) = A094406(14)** — new derived value for the mixed cycle-reach record
   sequence, via exhaustive metric-decomposition verification and a digit-count dominance
   proof over the happy branch (published A001273 data consumed as cited input).
   Also derived: A176762(20) = A094406(13).

## Why the method matters (naive search fails)
The answer has ~7.14·10^196 digits — no enumeration of candidate numbers is possible.
Instead: tail(m)=n forces ssd(m) ≥ a(n−1); m ≥ M*(ssd(m)); any better competitor satisfies
ssd(m) ≤ 81·digits(U). That collapses minimality to an exhaustive scan of 108 big integers
plus an exact O(polylog) characterization of minimal preimages via per-digit deficits
V = {0,17,…,81}.

## Verification standard applied
* Primary implementation (`extend.py`, `mstar.py`): DP sumset tables, greedy lex-min
  construction, brute-force cross-validation on σ ≤ 350 (350/350), published-term
  reproduction for both sequences (all 33 terms match OEIS).
* Independent checker (`checker_independent.py`): meet-in-the-middle deficit feasibility,
  fresh simulations, re-derived bounds, lexicographic-minimality recheck — ALL PASS.
* Adversarial tamper suite: 8/8 corrupted certificate variants correctly rejected
  (including a same-multiset digit transposition that passes ssd checks but fails lex-minimality).
* Claim labels: `proven_result` (see sessions/S4/REPORT.md §4 for the one inherited input).

## Artifacts
* `research/experiments/S4/` — code + raw outputs + `run.sh` (one command, <2 min).
* `research/sessions/S4/REPORT.md`, `research/sessions/S4/SPEC.md`.
* `research/candidates/C-S4-01.json`, `research/candidates/C-S4-02.json`.
* `research/literature/S4-novelty.md` — novelty audit, sources dated 2026-08-26.

## Honest limitations
a(15) is beyond explicit certification here (its nine-count would itself have ~7·10^196
digits); general web search engines blocked automation, so novelty rests on OEIS state,
entry-cited literature, and arXiv. Rejected targets (distinct-cubes/squares exact-count
records) documented in REPORT.md §1.
