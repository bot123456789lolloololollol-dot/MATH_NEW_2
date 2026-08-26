# S2 Novelty Audit — bounded-lookahead rollout LPT (ROL) for P||Cmax

Auditor: Session S2. All queries run on 2026-08-26 (EDT). Search indexes used:
DBLP publication search API, arXiv API, Semantic Scholar API (rate-limited; partial),
Mojeek web search (DuckDuckGo bot-blocked, Bing timed out), Wikipedia.

## What we searched and found

| Query / check | Index | Date | Result |
|---|---|---|---|
| "lookahead makespan" | DBLP | 2026-08-26 | 0 hits |
| "rollout scheduling parallel" | DBLP | 2026-08-26 | 1 hit: "An efficient rollout algorithm for unrelated parallel machine scheduling with random rework", Int. J. Prod. Res., 2025 — DIFFERENT problem (unrelated machines, stochastic rework); shows the rollout concept exists in scheduling, but not for P||Cmax dispatch |
| "lookahead machine scheduling" | DBLP | 2026-08-26 | 12 hits, ALL online/semi-online scheduling where "lookahead" means extra information about FUTURE jobs (e.g., "Semi-Online Scheduling on Two Identical Parallel Machines with Initial-Lookahead Information", 2024) — different meaning from our offline bounded-lookahead simulation |
| all:"makespan" AND (all:"rollout" OR all:"one-step lookahead") | arXiv API | 2026-08-26 | 0 entries |
| "lookahead scheduling identical parallel machines makespan" | Semantic Scholar | 2026-08-26 | rate-limited (HTTP 429); no results retrieved |
| "COMBINERULE" scheduling | DBLP + Mojeek | 2026-08-26 | 0 relevant hits; rule name not locatable in any index we can reach |
| "LRM" parallel-machines dispatch rule | Mojeek | 2026-08-26 | no relevant hits |
| "beam search" scheduling heuristics | DBLP | 2026-08-26 | hits exist but for single-machine early/tardy problems (2008-2012), NOT P||Cmax |
| "Bounds multiprocessing timing anomalies" (baseline citation check) | DBLP | 2026-08-26 | CONFIRMED: R.L. Graham, SIAM J. Appl. Math., 1969 |
| "Coffman bin-packing multiprocessor scheduling" (baseline citation check) | DBLP | 2026-08-26 | CONFIRMED: E.G. Coffman Jr., M.R. Garey, D.S. Johnson, SIAM J. Comput., 1978; Wikipedia (accessed 2026-08-26) confirms MULTIFIT name and 13/11 approximation bound |
| "Lee multiprocessor MULTIFIT" | DBLP | 2026-08-26 | CONFIRMED: "Multiprocessor scheduling: combining LPT and MULTIFIT", Discrete Applied Mathematics, 1988 — closest-in-spirit prior art (combining LPT with a bin-packing search) |

## Verdict on our mechanism

Mechanism: for each of the first K largest jobs, try every distinct-load machine,
complete the schedule with plain LPT (policy simulation), commit the placement that
minimizes the simulated final makespan; finish with plain LPT.

- CLOSEST PRIOR ART:
  1. Rollout algorithms as a generic class (policy-simulation improvement of greedy
     heuristics) — established concept in the literature; the 2025 IJPR paper above uses
     rollout thinking for an unrelated-machines variant.
  2. Lee & Massey (1988): hybridizing LPT with bin-packing search — same goal
     (improving LPT) via a different mechanism.
  3. Beam-search scheduling heuristics — partial-decision exploration, but applied to
     single-machine weighted problems, and using bounding functions rather than greedy
     completion rollouts.
- NO prior art found for a K-bounded LPT-rollout dispatch heuristic specifically for
  P||Cmax in the indexes reachable to us. However, the general mechanism class is well
  known, so this is an INCREMENT (novel instantiation), not a new principle.
- Search limitations honestly stated: Semantic Scholar was rate-limited; Google Scholar
  was not directly queryable from this environment; paywalled full texts were not
  scanned. A reviewer with better access may find e.g. textbook exercises or minor
  workshop papers using "lookahead LPT".

Novelty status: audited_clean_with_caveats — no blocking prior art found for the exact
mechanism; mechanism class known; claim labels kept at "experimentally_validated_result"
for empirical statements only.
