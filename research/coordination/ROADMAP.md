# Global Roadmap — v1 (2026-08-26)

## Phase 0 — Infrastructure (done)
Repo structure, protocol, scoring rubric, five mission briefs.

## Phase 1 — Parallel exploration (running)
Five sessions, distinct domains, hard quality gates:
- S1 Extremal combinatorics: certified exact values for small-case packing/difference problems.
- S2 Combinatorial optimization: scheduling heuristics vs honest baselines (preregistered seeds).
- S3 Adversarial math: counterexample hunting against open conjectures at reachable sizes.
- S4 Computational number theory: extend hard integer sequences beyond published range w/ certificates.
- S5 Circuit/logic optimization: SAT-based synthesis of sorting networks; validate tooling on known optima, attack open cases.

## Phase 2 — Triage & promotion
Coordinator scores all candidates; top 1–2 promoted to adversarial verification round.

## Phase 3 — Adversarial round
Independent attacker agents: counterexample search, certificate fuzzing, hidden-assumption audit,
benchmark-leakage check, deeper literature search for equivalent formulations.

## Phase 4 — Independent reproduction
Fresh agent reimplements the verifier/spec from the candidate's SPEC section only.

## Phase 5 — Final report + archival
`reports/FINAL_REPORT.md`, discoveries vs failed archived honestly, leaderboard finalized, commits tagged.

## Kill criteria
Any candidate whose novelty audit finds direct prior work, whose evidence is not regenerable, or
whose gain is within noise goes to `/failed` with post-mortem. No zombie candidates.
