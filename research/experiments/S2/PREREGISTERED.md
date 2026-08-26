# S2 PREREGISTRATION — P||Cmax heuristic comparison

Written: 2026-08-26 00:42 EDT (before any comparison experiment was run; before instance
generation and before any algorithm was executed on any generated instance).
This file fixes seeds, configurations, metrics, algorithms, and the statistical protocol.
No seed may be added or changed after results are seen. Any deviation must be flagged as
post-hoc in REPORT.md.

## Problem

P||Cmax: schedule n independent jobs with integer processing times p_j on m identical
parallel machines to minimize makespan Cmax.

## Instance generation (fixed before runs)

Distribution-faithful reimplementation of Taillard's published generator description:
processing times are i.i.d. discrete uniform integers in {1,...,99}
(source: E. Taillard, "Benchmarks for basic scheduling problems", EJOR 64 (1993) 278-285;
his paper specifies the uniform [1,99] distribution for the P||Cmax benchmarks. We do not
reproduce his exact PRNG streams — his files use his own LCG state which is not fully
specified in the paper — we reproduce the distribution he specifies, with our own fixed
seeds via numpy PCG64: `numpy.random.default_rng(seed).integers(1, 100, size=n)`).

- Configurations (ALL of these, no others): m in {2,3,4,5,7,10} x n in {30,50,100,200,500}.
- DEV_SEEDS = [1000..1029] (30 seeds) — used ONLY for tuning candidate parameters.
- EVAL_SEEDS = [3000..3029] (30 seeds) — used ONLY for the final evaluation.

> ERRATUM (2026-08-26 00:47 EDT, before instance generation and before ANY algorithm was
> executed on any generated instance): the original text said "35 configs"; the correct
> count for the SAME preregistered sets is 6 x 5 = 30 configurations. No configuration,
> seed, metric, or algorithm was changed. Additionally, the dev subset originally omitted
> an m=4 configuration, leaving K*(m=4) undefined; we extend D by (4,100) so that every
> evaluation m value {2,3,4,5,7,10} has a dev configuration. Both corrections are made
> before any run; nothing below is interpreted from results.
- Instance file name: `taillike-m{m}-n{n}-s{seed}.txt` in `research/benchmarks/S2-instances/`.
- Instances are committed data files; all algorithms read exactly these files. The same
  instance file is fed to every algorithm (paired design).

## Metric

gap% = 100 * (Cmax_alg - LB) / LB, with
LB = max( max_j p_j , ceil(sum_j p_j / m) ).
Lower is better. All comparisons are PARED per instance (same instance across algorithms).

## Algorithms

Baselines (implemented from primary literature):
- B1 LPT — Graham, R.L. (1969), "Bounds on multiprocessing timing anomalies",
  SIAM J. Appl. Math. 17(2):416-429. Sort jobs non-increasing p; assign each to least-loaded
  machine (tie: lowest machine index).
- B2 MULTIFIT — Coffman, E.G., Garey, M.R., Johnson, D.S. (1978), "An application of
  bin-packing to multiprocessor scheduling", SIAM J. Comput. 7(1):1-17. Binary search on bin
  capacity c over bracket [LB, sum p]; feasibility test = first-fit-decreasing packing into m
  bins of capacity c; iterate until convergence or 40 halvings. (Original paper's tighter
  initial upper bound replaced by the safe bracket above; documented deviation.)
- B3 LPT+LS — LPT followed by local search to a local optimum over the union neighborhood
  {single-job relocation (jump), pairwise inter-machine interchange (swap)}, best-improvement
  (steepest descent), deterministic tie-breaks (first found in fixed scan order). A move is
  accepted iff it strictly reduces the current makespan. Scan restricted to moves involving
  >=1 critical machine (a machine whose load equals Cmax); this loses NO strictly improving
  move because any strict improvement must strictly reduce every critical machine's load.
  This is "the bar to beat".

Candidate (our mechanism):
- C1 ROL ("bounded-lookahead rollout LPT"): jobs sorted non-increasing. For each of the first
  K jobs (K >= 1), tentatively place it on each machine choice (machine choices deduplicated
  by identical current load; ties broken by lowest index), complete the schedule by running
  plain LPT on ALL remaining jobs from that tentative partial state, and evaluate the final
  rollout Cmax; commit the placement minimizing the rollout Cmax (tie-break: smaller resulting
  load of chosen machine, then lower machine index). The remaining n-K jobs are finished with
  plain LPT (no rollout). Mechanism class: simulation-based decision making (rollout policy
  evaluation), distinct from static priority rules (LPT), bin-packing search (MULTIFIT), and
  neighborhood descent (LS).
- Parameter K tuned ONLY on DEV_SEEDS x dev configs D = {(2,30),(3,50),(4,100),(5,100),
  (7,200),(10,500)} [erratum: (4,100) added pre-run, see above],
  grid K in {m, 2m, 4m}. Selection rule: for each m with an eval configuration, minimize
  mean gap over the dev configuration(s) of that m; tie -> smaller K.
  K* frozen before ANY eval-seed run. Eval uses K*(m) per configuration.

## Statistical protocol (fixed)

- Primary endpoint: ROL vs LPT paired gap difference per configuration.
- Secondary endpoints: ROL vs LS, ROL vs MULTIFIT, and MULTIFIT vs LPT (descriptive).
- Test: two-sided paired Wilcoxon signed-rank test (scipy.stats.wilcoxon, zero_method='wilcox'),
  alpha = 0.05, Holm-Bonferroni correction within each configuration across the 3 comparisons
  {ROL-LPT, ROL-MULTIFIT, ROL-LS}. Sensitivity check: paired t-test reported alongside.
- Effect sizes: matched-pairs rank-biserial correlation r_rb and Cohen's d_z =
  mean(diff)/std(diff, ddof=1). Also pooled-across-all-configs test per comparison.
- Decision rule: claim "improvement" only if Holm-adjusted p < 0.05 AND mean gap difference
  favors ROL; otherwise label "within noise" / null result.
- n = 30 paired runs per configuration (EVAL_SEEDS). No exclusions; all runs reported,
  including failures/exceptions if any occur.

## Validation gates (must pass before eval)

- LPT on hand instance p=[3,3,2,2,2], m=2 must give Cmax=7 (Graham's classic example;
  OPT=6, LB=6, gap=16.67%). LS reaches its own local optimum there (documented value).
- Independent brute-force optimal (n<=8, m<=3) cross-check: algorithms produce valid
  schedules (all jobs assigned once) and correct Cmax recomputation from assignments.
- LS output must satisfy an independent naive O(n^2 m^2) no-improving-move checker.
- Determinism: full pipeline run twice on 3 spot configs must give byte-identical outputs.

## Falsifiable hypotheses (stated now, before results)

- H1 (primary): ROL(K*) has significantly smaller gap than LPT on the majority of the 30
  configs (Holm-adjusted p<0.05, favorable sign) and pooled across configs.
- H2 (openly expected possible null): ROL does NOT beat LPT+LS; LS is expected to be the
  strongest method on these uniform instances.
- H3: MULTIFIT's mean gap <= LPT's mean gap on most configs.

A null outcome on H1 will be reported honestly as such (archive-worthy negative).
