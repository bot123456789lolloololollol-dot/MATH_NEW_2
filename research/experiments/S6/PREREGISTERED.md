# S6 Preregistration — Equation Discovery Study

Written **before** any experiment was run (2026-08-26). Seeds, metrics, sample sizes,
decision thresholds, and analysis procedures below are fixed in advance; deviations are
allowed only if explicitly flagged as post-hoc in the report.

## Global seed policy

- Master seed: `20260826`.
- Experiment *k* uses seeds `1000*k + [0..29]` wherever multiple stochastic runs are needed.
- All RNG via `numpy.random.default_rng(seed)`. No seed is chosen after looking at results.

## Systems under study (fixed)

| id | system | params | purpose |
|---|---|---|---|
| SYS1 | damped linear oscillator | k=2.0, c=0.3 | sanity recovery |
| SYS2 | Lotka–Volterra | a=1.0, b=0.1, c=1.1, d=0.4 | nonlinear recovery |
| SYS3 | Lorenz-63 | σ=10, β=8/3, ρ=28 | chaotic recovery + horizon analysis |
| SYS4 | nonlinear pendulum | g/L=1.0 (dimensionless) | period-law generalization |
| SYS5 | two-body gravity | μ varied | Kepler law + eccentricity-invariance |
| SYS6 | series RLC | R,L,C grids | hidden parameter relations |
| SYS7 | conservative Duffing x''=x−x³ | — | invariant discovery |
| SYS8 | damped Duffing (δ=0.3) / forced variants | — | invariant negative control |

## Data-collection protocol (fixed)

- Integrator: `scipy.integrate.solve_ivp`, RK45, rtol=1e-10, atol=1e-12.
- Sampling: uniform grid; dt = 1e-3 (SYS3), 2e-3 (SYS1/2), 5e-3 (SYS4), 1e-2 (SYS5/6).
- Noise: additive Gaussian on states, σ expressed as fraction of per-channel state std;
  levels {0, 0.005, 0.01, 0.02, 0.05}; derivatives from central differences
  (Savitzky–Golay window 11, polyorder 3 when σ>0).
- Train/held-out split: distinct initial conditions and/or distinct parameter values,
  never a random row split of the same trajectory.

## Metrics (fixed definitions)

- Coefficient error: `max_i |ĉ_i − c_i| / max_j |c_j|` over true-active terms.
- Support recovery: Jaccard similarity of active-term sets.
- Rollout error: relative Frobenius error of discovered-vs-true model integrated from an
  unseen IC over a fixed horizon (SYS1/2: T=10; SYS3: T=2).
- Predictability horizon t*: smallest t with normalized trajectory deviation > 0.25.
- Invariant quality: after best affine fit to ground-truth H, residual std along held-out
  trajectories divided by range(H); detection criterion for "invariant exists":
  smallest-to-largest singular-value ratio < 1e-8 AND spectral gap σ[-2]/σ[-1] > 100.
- Function-fit NMSE: MSE divided by Var(y) on held-out data.

## Decision thresholds (fixed before runs)

- "Equation recovered" ⇔ coefficient error < 1e-6 (clean) or < 10×(noise floor bound),
  support Jaccard = 1, AND held-out rollout error < 5% over the fixed horizon.
- "Law generalizes beyond discovery regime" ⇔ held-out-regime NMSE < 5% using the exact
  discovered symbolic form (no refitting on test data).
- "No law found" (honest-failure verdict) ⇔ held-out NMSE > 1e-2 after exhaustive search;
  this verdict must be reported rather than a spurious fit whenever triggered.
- Statistical comparisons (exp09): ≥30 paired runs per method, report mean ± std and
  two-sided Wilcoxon signed-rank p-value at α=0.05.

## Baselines for comparison (fixed)

Plain least squares on full library (OLS-full), Ridge (α=1e-3 on standardized columns),
Lasso (α picked by BIC over 20-point log grid), and genetic-programming symbolic regression
(population 300, generations 80, tournament 3, p_cx=0.7, subtree-mutation 0.2).

## Planned experiments

exp01 known-law recovery · exp02 noise sweeps · exp03 Occam/equivalence tests ·
exp04 pendulum period extrapolation · exp05 Kepler III + e-invariance · exp06 RLC relations ·
exp07 invariant discovery + negative controls · exp08 adversarial suite (noise-only,
out-of-class truth y=tanh(x), confounded causal test, chaos horizon) · exp09 baseline study.

## Post-hoc deviations

- 2026-08-26 (during library smoke-testing, before any recorded experiment): derivative
  operator upgraded from 2nd-order (`np.gradient`) to 4th-order central stencils with
  trimmed edges. Reason: O(dt^2) truncation bias projected onto spurious library terms
  at the 1e-6 level; the 4th-order operator removes the bias rather than loosening the
  recovery criterion. No result files existed when this change was made.
- 2026-08-26 (bring-up): stencil upgraded once more to 7-point 6th order after the
  O(dt^4) residue was observed to leave ~1e-7 "dust" terms in Lorenz; model selection
  changed from BIC to chronological validation-split RSS for the same reason (BIC's iid
  residual assumption is violated by correlated FD error). Both changes affect all
  experiments symmetrically; recovery criteria unchanged.
- 2026-08-26 (bring-up): support-recovery tolerance clarified as scale-relative
  |c| < 1e-6 * max|c_true| per system (an earlier absolute tolerance misclassified
  unit-coefficient terms in Lorenz). Subsequently replaced entirely by a z=5
  significance test (|c| > 5 standard errors from the OLS covariance on the selected
  support), which is the correct noise-floor rule under P2 of THEORY.md; a fixed
  tolerance either counts fitted noise as active or hides real small terms.
- 2026-08-26 (bring-up, adversarial A2 only): verdict requires (i) holdout NMSE < 1e-2,
  (ii) compactness <= 25 nodes (positive-case truths have <= 17), and (iii) domain
  transfer: NMSE < 1e-2 on [3.5,6] u [-6,-3.5], outside the training interval.
  Rationale discovered during bring-up: a 16-node expression fits tanh(2x) in-domain
  without representing its mechanism; without (iii) curve fitting would pass as law.


