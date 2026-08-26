# Session Mission Briefs — Phase 1

Five independent sessions, distinct domains, no overlap. Each session: implement, benchmark,
attack own work honestly, novelty-audit against literature, and write up with strict claim labels.

---

## S1 — Extremal combinatorics: certified exact values
**Domain:** combinatorial constructions / difference packing.
**Mission:** Find a parameterized extremal-combinatorics family whose published tables have gaps at
sizes reachable by exhaustive search on this machine (candidates to evaluate, pick one):
cyclic/modular Golomb rulers (difference triangle sets) for small n and m marks; maximum Sidon
sets in Z_n for n where tables are incomplete; graceful/harmonious labeling existence for small
families; perfect difference families. Build a deterministic backtracking/SAT search with justified
symmetry breaking. Deliver either (a) a new exact value with completeness argument + checker, or
(b) a new best-known construction + checker, or (c) an honest saturation report (table already
complete → document and fall back within family).
**Gate:** any "exact/optimal" claim requires the completeness argument in writing; every object ships
with an independent checker script.

## S2 — Combinatorial optimization: scheduling heuristics, honest baselines
**Domain:** P||Cmax makespan minimization.
**Mission:** Implement Taillard-instance generation (published generator, seeded), a strong baseline
(LPT + pairwise-swap improvement; plus a known strong dispatch rule from literature), then attempt a
genuinely different mechanism (e.g., bounded lookahead, or a small parameterized rule family tuned by
preregistered search). Preregister instance seeds AND evaluation protocol in `PREREGISTERED.md`
BEFORE running comparisons; ≥30 runs per config; report mean±std gaps vs LPT lower bound... (use
known lower bounds max(p_max, sum/p)); paired significance test. Novelty audit against known rules
(LPT, COMBINERULE, LRM, CDJS, etc.). A null result is acceptable if documented rigorously.
**Gate:** no seed may be chosen after seeing results; all raw outputs committed.

## S3 — Adversarial math: counterexample hunt
**Domain:** SAT/constraint solving applied to open conjectures.
**Mission:** Assemble ≥5 precise published open conjectures (or recently-settled ones worth auditing)
that have finite verifiable consequences at reachable sizes (graph theory conjectures with small
counterexample candidates; digit/number properties; tiling/coloring claims). For each: formalize,
brute force/SAT-solve to the largest reachable bound, record verdict:
REFUTED (counterexample found + verified independently) / VERIFIED_UP_TO_N / INFEASIBLE.
Any refutation must be reproduced with a second independent implementation before it is claimed.
Check literature first — do not waste effort on already-settled conjectures.
**Gate:** verdict files per conjecture; counterexamples ship with checker scripts.

## S4 — Computational number theory: extend hard sequences
**Domain:** number-theoretic constructions/records.
**Mission:** Pick 1–3 integer sequences whose published/OEIS terms end at a value reachable with a
smarter-than-naive search on this machine (e.g., minimal numbers with exactly k representations of
some form; least number needing k steps under a deterministic process; maximal-determinant records
for small orders; Egyptian-fraction extremes). Verify exactly where the published range ends
(OEIS b-file, literature), extend it with certified computation (checker script per term), and
document why naive search would have failed. Cross-check each new term against OEIS/Literature for
concurrent discoveries.
**Gate:** every new term has a machine-checkable certificate; OEIS cross-check dated in writeup.

## S5 — Circuit/logic optimization: sorting networks via SAT
**Domain:** circuit/logic synthesis.
**Mission:** Implement a SAT encoding of size-minimal sorting networks (0-1 encoding of comparators +
permutation constraints; validate against literature-known optimal sizes for small n, e.g., n≤8
exactly as a tooling sanity check). Then attempt the smallest genuinely open case reachable
(literature check required for current best-known upper/lower bounds — cite them). Deliverable ladder:
(a) new smaller network (upper-bound improvement, certificate = explicit network + checker),
(b) UNSAT-based lower-bound improvement (requires completeness argument for the encoding),
(c) validated open-source encoder + narrowed bound on a variant (selection/merging networks),
(d) honest negative result documenting solver scaling limits with data.
**Gate:** networks checked by an independent brute-force verifier over all 0-1 inputs; all bounds cited
with sources and dates.

---

## Shared output contract (all sessions)

Write into the repo (and nowhere else outside your listed touchpoints):
1. `research/sessions/S{n}/REPORT.md` — findings with claim labels + regenerable commands.
2. `research/sessions/S{n}/SPEC.md` — precise spec of the key artifact so a stranger can reimplement.
3. `research/candidates/C-S{n}-01.json` (from TEMPLATE.json) — one record per serious candidate.
4. `research/experiments/S{n}/` — code + raw outputs + run scripts.
5. `research/literature/S{n}-novelty.md` — novelty audit with sources and access dates.
6. `research/leaderboard/session_S{n}.json` — self-assessment scores with one-line justifications.

Rules: no git commands; timebox experiments (≤20 min wall-clock each unless essential); deterministic
seeds; Python (`python` = 3.11) is available; Node 24 available; NO gcc.
