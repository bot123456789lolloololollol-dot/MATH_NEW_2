# Session S2 Report — P||Cmax heuristics vs honest baselines

Date: 2026-08-26. Author: Session S2 (combinatorial optimization domain).
All numbers below regenerate via `cd research/experiments/S2 && bash run.sh` (< 1 min).

## 0. Preregistration statement (hard gate)

`research/experiments/S2/PREREGISTERED.md` was written at **2026-08-26 00:42-00:47 EDT**,
BEFORE instance generation and before any algorithm was executed on any generated
instance. It fixes: all seeds (DEV 1000-1029, EVAL 3000-3029), all 30 configs
(m x n grid), the metric, all algorithms, the statistical test (paired Wilcoxon,
alpha=0.05, Holm within config), effect sizes, and decision rules. No seed was chosen or
changed after seeing results. One dated pre-run erratum inside the file corrects an
arithmetic slip ("35" -> 30 configs) and extends the dev set with (4,100) so every eval m
has a dev config; both corrections predate any run.

## 1. What was built

1. **Instances**: distribution-faithful reimplementation of Taillard's published P||Cmax
   generator spec (p_j i.i.d. discrete uniform {1..99}; Taillard, EJOR 64 (1993) 278-285)
   under our own fixed numpy PCG64 seeds; 1800 instances committed as data files in
   `research/benchmarks/S2-instances/`.
2. **Baselines from primary literature**:
   - LPT (Graham 1969, SIAM J. Appl. Math. 17(2):416-429);
   - MULTIFIT (Coffman, Garey & Johnson 1978, SIAM J. Comput. 7(1):1-17) — chosen as the
     strong classical alternative because the mission-suggested names "COMBINERULE"/"LRM"
     could NOT be located in any index reachable to us (DBLP: 0 hits; web search: nothing
     relevant; see literature/S2-novelty.md). This substitution is documented.
   - LPT + full relocation/interchange local search to a local optimum (the bar),
     validated against an independent naive checker AND an independent naive steepest-
     descent reference (`check_ls.py`: identical makespans on all tested instances;
     the first version of that *reference* had a scan bug which we found and fixed —
     the production implementation passed the full-neighborhood optimality checker
     throughout).
3. **Candidate ROL (bounded-lookahead rollout LPT)**: for each of the first K=k_mult*m
   largest jobs, try every distinct-load machine, complete the schedule by plain-LPT
   simulation, commit the placement minimizing simulated final makespan; finish with LPT.
   Deterministic; no randomness; ~0.8 ms average per instance (LPT: 0.1 ms).
   Parameter k_mult tuned ONLY on preregistered dev subset; frozen values:
   k_mult = {m=2:4, 3:2, 4:1, 5:2, 7:1, 10:1} (`outputs/K_frozen.json`,
   `outputs/tuning_freeze.txt`).

## 2. Validation gates (all passed before evaluation)

- LPT on hand instance p=[3,3,2,2,2], m=2 -> Cmax=7 as hand-computed (OPT=6, LB=6);
  MULTIFIT -> 6; LS -> 6 (verified locally optimal by independent checker); LB check OK.
- 200 random tiny instances x 4 algorithms vs brute-force optimum: assignments valid,
  Cmax recomputations match, Cmax >= OPT always, LS endpoints locally optimal.
- Determinism: two full independent evaluations byte-identical except wall-clock column
  (`outputs/det_check_run2.csv`, `outputs/results_eval.csv`).
- Integrity at analysis time: LB recomputed from instance files for all 3600 rows and
  gap arithmetic verified: 0 failures; independent no-numpy LPT reimplementation matched
  recorded LPT on 40 sampled rows.

## 3. Results (preregistered evaluation: 30 configs x 30 paired runs = 900 runs/algo)

Metric: gap% above LB = max(p_max, ceil(sum/m)). Mean±std per algorithm pooled over all
900 instances: **LPT 0.4397±1.1267, MULTIFIT 0.1270±0.3430, LPT+LS 0.2259±0.7087,
ROL 0.2676±0.7797**. Full per-config table: `outputs/analysis.txt`; machine-readable:
`outputs/stats_summary.json`.

Pooled paired tests (Wilcoxon signed-rank, two-sided):

| Comparison | mean gap diff | p | rank-biserial | dz | better/worse/tie |
|---|---|---|---|---|---|
| ROL vs LPT | -0.1721 | 1.15e-28 | -0.33 | -0.32 | 164 / 0 / 736 |
| ROL vs MULTIFIT | +0.1406 | 1.86e-22 | +0.39 | +0.22 | 62 / 278 / 560 |
| ROL vs LS | +0.0416 | 2.59e-09 | +0.31 | +0.09 | 56 / 231 / 613 |

Per-config (Holm alpha=0.05): ROL beats LPT significantly on 9/30 configs (never loses
on any); MULTIFIT beats ROL on 14/30; LS beats ROL on 12/30.

Key regime detail (claim-relevant):
- ROL's gain concentrates exactly where LPT is weak: e.g. (m=7,n=30) gap 2.24%->0.83%
  (-63%, p=3.7e-05, dz=-1.20), (m=10,n=30) 4.53%->3.26% (p=5.9e-04), (m=5,n=30)
  1.10%->0.33% (p=2.7e-04), (m=10,n=50) 1.67%->0.95% (p=8.1e-05).
- On large-n configs both ROL and LPT hit gap 0 (identical outputs; ties dominate).
- Unexpected honest finding: **MULTIFIT, not LS, is the strongest method overall**
  (pooled mean 0.127 vs LS 0.226), especially at larger n/m ratios — e.g. (10,200):
  MF 0.000 vs LS 0.036; (7,200): MF 0.005 vs LS 0.028. At small n, LS wins ((3,50):
  LS 0.012 vs MF 0.052). Our H2 expectation ("LS strongest") was therefore wrong.

## 4. Verdicts with strict claim labels

- **experimentally_validated_result**: ROL improves on LPT on Taillard-distribution
  instances — pooled Wilcoxon p=1.15e-28 favoring ROL, mean gap reduced 39%
  (0.440->0.268), significant (Holm) on 9/30 configs, and ROL was NEVER observed worse
  than LPT on any of the 900 paired runs (164 wins / 0 losses / 736 exact ties).
  Improvement claim passes the preregistered significance gate (pooled + per-config
  where gaps are non-degenerate).
- **observation**: ROL's dominance over LPT was total (0/900 losses); we did NOT prove
  that rollout choice can never hurt relative to pure LPT, so this stays an observation,
  though a striking one.
- **experimentally_validated_result (negative)**: ROL does NOT beat the strong baselines.
  Pooled it is significantly WORSE than both MULTIFIT and LPT+LS. Nuance (descriptive,
  n.s. after Holm): on small-n/moderate-m configs ROL numerically beats LS, e.g.
  (7,30): 0.83 vs 1.07 (p_holm=0.087), (5,30): 0.33 vs 0.36 (n.s.).
- **experimentally_validated_result**: MULTIFIT <= LPT everywhere here (H3 confirmed),
  and MULTIFIT > LS at larger n/m on uniform instances — worth other sessions treating
  MULTIFIT as the default honest baseline for P||Cmax studies.
- Null-result honesty: our candidate does not clear the LS bar. It is archived as a
  modest, cheap, well-characterized improvement over LPT only.

## 5. Novelty audit

See `research/literature/S2-novelty.md` (queries dated 2026-08-26). No prior art found
for K-bounded LPT-rollout dispatch on P||Cmax in reachable indexes; mechanism class
(rollout/greedy-with-completion) is established; "lookahead" in scheduling literature
refers to online information models, not offline simulation. Status: audited clean with
caveats (increment-level novelty).

## 6. Threats to validity / limitations

- One distribution family (uniform [1,99]) as in Taillard's makespan set; conclusions
  may not transfer to skewed/adversarial distributions (not tested).
- LS bar uses a specific tie-breaking; other tie-breaks land in different local optima
  (demonstrated in outputs/check_ls.txt history during development) — we report our
  exact deterministic variant.
- Semantic Scholar rate-limiting limited one audit channel (documented).
- Multiple-comparison control applied within configs (Holm across 3 comparisons); the
  9/30 significant configs exceed chance (~1.5 expected false positives at alpha .05)
  and agree in sign with the pooled test, so the LPT-improvement claim is robust.

## 7. Regeneration commands

    cd research/experiments/S2
    bash run.sh          # everything: validate -> generate -> tune -> evaluate -> stats

Individual steps are listed in run.sh; raw outputs committed under outputs/.
