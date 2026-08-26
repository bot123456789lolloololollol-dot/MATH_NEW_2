# S2 SPEC — Bounded-lookahead rollout LPT (ROL) for P||Cmax, with honest baselines

This spec is sufficient for a stranger to reimplement everything. All code lives in
`research/experiments/S2/`; one command (`bash run.sh`) regenerates every number.

## 1. Problem

P||Cmax: n independent jobs, integer processing times p_j, m identical parallel machines;
minimize the maximum machine load (makespan).

## 2. Instances

Distribution-faithful reimplementation of Taillard's published P||Cmax generator
specification: p_j i.i.d. discrete uniform on {1,...,99}
(source: E. Taillard, "Benchmarks for basic scheduling problems", EJOR 64 (1993) 278-285).
His original PRNG stream is not specified in the paper; we use numpy PCG64:

    rng = numpy.random.default_rng(seed)
    p = rng.integers(1, 100, size=n)   # dtype int64, values in {1..99}

Configs: m in {2,3,4,5,7,10} x n in {30,50,100,200,500} (30 configs).
Seeds (preregistered before any run, see PREREGISTERED.md):
DEV = 1000..1029 (tuning only), EVAL = 3000..3029 (final evaluation only).
Files: `research/benchmarks/S2-instances/taillike-m{m}-n{n}-s{seed}.txt`;
line 1 "m n seed", line 2 the n integers. Instances are committed data.

## 3. Metric

gap% = 100 * (Cmax_alg - LB) / LB, LB = max(max_j p_j, ceil(sum_j p_j / m)).
Comparisons are paired per instance.

## 4. Algorithms

### B1 LPT (Graham 1969, SIAM J. Appl. Math. 17(2):416-429)
Sort jobs by non-increasing p (ties: ascending original index). Place each job on the
least-loaded machine (ties: lowest machine index).

### B2 MULTIFIT (Coffman, Garey & Johnson 1978, SIAM J. Comput. 7(1):1-17)
Binary search the bin capacity c in [max(p_max, ceil(S/m)), S], S = sum p:
feasibility(c) = first-fit-decreasing packing of all jobs into m bins of capacity c
(bins scanned in index order). Iterate until lo == hi (<= ~40 bisections); final Cmax =
smallest feasible c found. Deviation from the paper: we use this safe bracket instead of
the paper's tighter initial upper bound (behaviorally equivalent after convergence).

### B3 LPT+LS ("the bar")
Start from LPT. Repeat until no improving move exists:
  - Let Cmax be current makespan; a move is improving iff it strictly reduces Cmax.
  - Neighborhood: single-job relocation (jump) and pairwise inter-machine interchange.
  - Scan restricted to moves involving at least one critical machine (load == Cmax);
    this loses NO strictly improving move because any strict improvement must strictly
    reduce EVERY critical machine's load, so moves touching no critical machine cannot help.
  - Accept the best move found (steepest descent). Tie-breaks: smaller resulting makespan,
    then relocation before interchange, then lower source machine, lower target machine,
    then first position in stored job order (numpy argmin first-occurrence semantics).
Termination is guaranteed (strictly decreasing integer makespan). Implementation is
numpy-vectorized per (critical machine, target machine) pair; validated against an
independent naive O(n^2 m^2) unrestricted steepest descent (check_ls.py: identical
makespans on all tested instances) and against an independent no-improving-move checker
(algorithms.ls_no_improving_move, used in validate.py on 200 random tiny instances).

### C1 ROL — bounded-lookahead rollout LPT (our mechanism)
Sort jobs non-increasing (ties: ascending original index), ps[0..n-1]. Let K = k_mult*m.
For idx = 0 .. K-1 (job duration d = ps[idx], remaining suffix rest = ps[idx+1..]):
  - For each machine mi whose CURRENT load value is distinct (machines with equal loads
    are symmetric; among equals keep the lowest index), tentatively set loads[mi] += d,
    simulate plain LPT on `rest` from these loads (each remaining job goes to the current
    least-loaded machine, ties to lowest index), obtaining final makespan val; undo.
  - Commit placement minimizing the tuple (val, loads[mi]+d, mi).
Remaining jobs ps[K..]: plain LPT. No randomness anywhere.
Mechanism class: simulation-based decision making (greedy rollout policy evaluation),
distinct from static priority rules (LPT), packing search (MULTIFIT), and neighborhood
descent (LS). Cost O(K * m * n * m) worst case (~milliseconds at our sizes).
Parameter protocol: k_mult in {1,2,4} tuned ONLY on DEV seeds x D =
{(2,30),(3,50),(4,100),(5,100),(7,200),(10,500)}; rule: per m, minimize mean dev gap,
tie -> smaller k_mult. Frozen result (outputs/K_frozen.json): {2:4, 3:2, 4:1, 5:2, 7:1, 10:1}.
EVAL seeds were never touched during tuning.

## 5. Evaluation protocol (preregistered)

See PREREGISTERED.md (written before any run; contains one dated pre-run erratum about
config count arithmetic and adding (4,100) to the dev set; no seed/metric/test changes).
Per config: 30 paired runs (EVAL seeds). Two-sided paired Wilcoxon signed-rank test
(scipy zero_method='wilcox'; degenerate all-zero-difference case reported as p=1),
Holm-Bonferroni within each config across {ROL-LPT, ROL-MULTIFIT, ROL-LS}; paired t-test
as sensitivity; effect sizes: matched-pairs rank-biserial r_rb (positive = first-listed
algorithm worse) and Cohen's dz.

## 6. Reproduction

    cd research/experiments/S2 && bash run.sh

Requires Python 3.11+, numpy, scipy. Everything is deterministic; two independent full
evaluations produced identical outputs modulo wall-clock timing column (verified,
outputs/det_check_run2.csv vs results_eval.csv).
