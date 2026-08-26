# S6 SPEC — Equation-Discovery Pipeline (reimplementation guide)

Audience: a stranger reimplementing the key artifacts from this description alone.
Environment used: Python 3.11.15, numpy 2.4.3, scipy 1.17.1, sympy 1.14.0,
scikit-learn 1.9.0, matplotlib 3.11.1 (Windows x64). No compiled extensions.

## A. Data generation (ground truth never shown to learners)

Integrate ODEs with `scipy.integrate.solve_ivp`, RK45, rtol=1e-10, atol=1e-12,
uniform sampling via `t_eval`. Systems and parameters are listed in
`experiments/S6/PREREGISTERED.md` (table SYS1..SYS8).

Period measurements:
- Pendulum: event detection on upward zero-crossings of theta (`direction=+1`);
  period = mean gap between consecutive crossing times. Cross-check against
  `4*ellipk(sin^2(theta0/2))`.
- Two-body: periapsis passages = events y=0 with direction=+1 (periapsis placed on
  +x axis, CCW start); period = mean gap. NOTE the start point registers an event
  at t~0; use mean of consecutive gaps, never first-to-last difference.
- RLC ring: peak spacing pi/omega_d; decay rate from linear fit of ln|peak amps|
  vs time, restricted to samples above 1e-10 * max|i|. Simulate at >=60 samples
  per fastest period.

## B. Derivatives

Seven-point 6th-order central stencil:
    dX = (X[6:] - 9X[5:-1] + 45X[4:-2] - 45X[2:-4] + 9X[1:-5] - X[:-6])/(60 dt)
aligned so dX[j] corresponds to state X[j+6] after trimming three more rows on each
side (`differentiate()` returns the aligned pair). With measurement noise: Savitzky-
Golay filter (window 11, polyorder 3, axis 0) before differentiating.
Rationale: lower-order stencils leave truncation residue that projects onto spurious
library terms at the 1e-6..1e-7 level (observed during bring-up).

## C. Sparse dynamics discovery (SINDy-style)

1. Library Theta(x): [1, all monomials up to degree p, optional sin/cos of chosen
   variables], column labels human-readable ("x*z", "sin(x)").
2. STLSQ: standardize columns by their norms; iterate {least squares -> zero out
   entries below threshold tau} until support stable; final refit on RAW columns of
   the surviving support (unbiased coefficients).
3. Threshold selection: chronological split of regression rows, first 75% / last 25%.
   For each tau on a 40-point geometric grid (std(dX)*1e-7 .. std(dX)*10):
   fit on train part only; score = RSS on validation part / number of validation
   entries; pick min score, ties -> fewer terms; refit chosen support on ALL rows.
   (BIC fails here because FD residual is correlated noise.)
4. Support significance: a term is reported active iff |c_ij| > 5 * SE(c_ij) in some
   target column, SE from the OLS covariance on the selected support.
5. Rollout validation: integrate the discovered RHS from unseen initial conditions;
   relative Frobenius error vs ground-truth rollout over a fixed horizon.

## D. Genetic-programming symbolic regression

Programs: binary trees over {add, sub, mul, protected div, sin, cos, clipped exp,
log|x|+eps guarded, sqrt|x|} x {variables, ephemeral constants}. Complexity = node
count. Fitness cost = n*ln(MSE_train) + kappa*nodes*ln(n), kappa=1 (MDL/BIC-style).
Search: population 300-800, generations 80-200, tournament(3) selection, elitism 5%,
crossover 0.7 (subtree swap), subtree mutation 0.2, depth cap ~8. Deterministic via
numpy default_rng(seed). Output BOTH the best-cost model and the full Pareto front
(nodes vs NMSE); final model choice may use held-out data for MDL scoring (exp03).
Known limitation (documented): vanilla GP reliably finds exact minimal forms for
trigonometric/log identities but can miss short root-form expressions like c/sqrt(g)
within these budgets.

## E. Conserved-quantity discovery

Features Phi(s): monomials up to degree r (constant term dropped).
Constraint matrix rows: Phi(S_k+1) - Phi(S_k) for consecutive samples, stacked over
MULTIPLE trajectories covering distinct energy levels (a single closed trajectory
makes extra spurious null directions -- verified failure mode).
Column-scale A by std of each column (NO centering: the constraint must stay
homogeneous); take SVD; the smallest right singular vector v gives c = v/sd_col;
F(s) = Phi(s) . c.
Detection test that an invariant exists at all:
    sigma_min/sigma_max < 1e-8  AND  spectral gap sigma_{r-1}/sigma_{r+...} > 100.
Conservation quality on held-out trajectories: max within-trajectory span of F
divided by cross-trajectory spread of mean(F) (never divide by within-trajectory
range of a conserved quantity -- it is ~0 by definition and ill-conditioned).
Compare to ground truth via affine fit a*H+b pooled across energies.

## F. Decision rules (preregistered)

- "equation recovered": coefficient error < 1e-6 (clean), significant support equals
  truth, held-out rollout error < 5%.
- "law generalizes": discovered symbolic form evaluated without refitting on
  held-out regime data, NMSE < 5%.
- "no law found" (honest failure): holdout NMSE > 1e-2 OR fails compactness (<=25
  nodes) OR fails domain transfer (NMSE < 1e-2 required outside training interval).
- Statistical comparisons: >=30 paired seeds, Wilcoxon signed-rank, alpha=0.05.

## G. File map & reproduction

    experiments/S6/
      src/{systems,sindy,symreg,invariant,metrics,plotting,common,pde}.py
      exp0{1..9}_*.py        # one script per experiment, writes results/*.json + figures/*.png
      exp09b_gplearn.py      # external GP baseline (requires pip install gplearn)
      exp10_pde.py           # PDE discovery: Burgers + KdV (round 2)
      exp11_control.py       # controller synthesized from discovered model (round 2)
      exp12_modes.py         # hidden modal structure in coupled oscillators (round 3)
      exp13_integral.py      # integral-form identification at high noise (round 3)
      run_all.py             # regenerates EVERYTHING: python run_all.py
      check_claims.py        # asserts REPORT.md headline numbers vs results/*.json
      PREREGISTERED.md       # seeds/protocol fixed before runs + deviation log
      PREREGISTERED_R2.md    # round-2 protocol + deviations
      PREREGISTERED_R3.md    # round-3 protocol
      tests/test_core.py     # pytest tests/ -q   (15 tests total incl. test_round2.py)
    sessions/S6/REPORT.md    # findings with claim labels (incl. round 2/3 addenda)
    sessions/S6/SPEC.md      # this file
    sessions/S6/proofs/THEORY.md, certificate_duffing.py

All RNG through numpy.random.default_rng with seeds from the preregistration files.

## H. Round-2/3 method notes (beyond sections A-F)

PDE discovery (exp10): spectral derivatives on a periodic grid; time derivative by
4th-order central stencil over snapshot spacing; library {1,u,u2,u_x,u_xx,u_xxx,
u u_x,u^2 u_x,u u_xx}; STLSQ as section C. CRITICAL requirements learned the hard way:
(1) data must contain more than one traveling wave (single waves make the library
collinear -- P9); (2) verify stencil signs against sin(3t) before trusting results;
(3) rollouts of discovered-vs-truth PDEs must use IDENTICAL relaxed tolerances.
Control loop (exp11): identify f_hat on bounded-angle excited data (input channel
subtracted using the KNOWN gain); feedback linearization u = -f_hat - kp th - kd om,
kp=4, kd=4; compare against identical controller built from true f.
Integral form (exp13): target X(t+W)-X(t) vs cumulative-sum-integrated library,
W=25 samples; never differentiate noisy states.
