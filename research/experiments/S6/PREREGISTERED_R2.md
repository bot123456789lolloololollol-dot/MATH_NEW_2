# S6 Preregistration — ROUND 2 addendum

Written before any round-2 experiment was run (2026-08-26 ~03:52 local). Extends
`PREREGISTERED.md`; all global rules (seeds, claim vocabulary, regenerability) carry over.

## R2-exp10 — PDE discovery from simulation data (numerical modeling)

- Systems (fixed): viscous Burgers `u_t = -u u_x + nu u_xx` (nu=0.05, periodic domain
  x in [0,2pi], IC u0(x)=sin(x)·(1+cos? fixed below)); KdV `u_t = -6 u u_x - u_xxx`
  (periodic, soliton IC). Integrator: method='RK45' rtol=1e-10/atol=1e-12 on a fine
  grid, or 'DOP853'; spatial derivatives by spectral (FFT) differentiation on the
  periodic grid -- learner may also use them (spectral derivs are standard numerics,
  not ground-truth leakage).
- Library (fixed): {1, u, u^2, u_x, u_xx, u_xxx, u*u_x, u^2*u_x, u*u_xx}.
- Recovery criterion (fixed): coefficient error < 1e-8 clean on active terms, support
  exact via z=5 significance, AND held-out prediction: discovered PDE integrated from a
  NEW initial condition matches truth over fixed horizon with rel L2 < 5%.
- Noise: additive Gaussian on u at sigma in {0, 1e-4, 1e-3} of max|u| (PDE regression
  is far more noise-sensitive than ODE case; scope limited accordingly).
- Seeds: 1000*k+i policy; Burgers/KdV sims deterministic (no RNG except noise).

## R2-exp11 — control loop built on the discovered model (control systems)

- System (fixed): driven pendulum theta'' = -(g/L) sin(theta) - c theta' + u,
  g/L=1, c=0.15. Learner sees only (theta, omega, u) data with u = filtered random
  excitation (preregistered seed), identifies f(theta, omega) in library
  {1, th, om, th*om, sin(th)? NO -- library {1, th, om, th^2, th*om, om^2} WITHOUT sin}
  plus KNOWN input gain b=1 (actuation physics assumed known; disclosed).
  Expectation: polynomial library CANNOT capture sin(th) exactly -> this tests control
  under model mismatch on a bounded domain near theta=0 where sin(th)~th is decent.
- Controller (fixed): feedback linearization u = -f_hat(theta,omega) + PD(theta_ref -
  theta) with PD gains kp=4, kd=2*sqrt(kp) (critically damped), computed ONLY from the
  discovered f_hat; applied to the TRUE plant in closed-loop sim.
- Metrics: settling time to |error|<0.05 rad from initial condition theta0=0.5, omega0=0;
  IAE (integral absolute error) over 10 s; same controller built from true f as oracle
  comparison. Success: IAE(discovered)/IAE(oracle) < 1.25 AND stable.
- Adversarial: repeat with initial condition theta0=2.5 (far from linearization-valid
  region) -- expect degraded but reported performance; no success claimed there.

## R2-exp09b — external-library GP baseline

- Install gplearn (pure Python); same three function classes as exp03, 30 paired seeds;
  metrics: holdout NMSE, node complexity of final expression; Wilcoxon vs our GP-MDL.
  Hypothesis: gplearn comparable accuracy, larger expressions (no MDL front).

## R2-checks — evidence-chain hardening

- `check_claims.py`: re-reads committed results/*.json and asserts every headline number
  quoted in sessions/S6/REPORT.md within stated tolerance; prints PASS/FAIL table.

## Post-hoc deviations (round 2)

- exp10 KdV (2026-08-26, after first run): single-soliton initial condition replaced by
  a TWO-soliton condition and snapshots refined from 101 to 401 with a 4th-order
  temporal stencil. Reason: one traveling wave makes all advective/derivative library
  terms nearly collinear (the strong-rank condition of THEORY P1 fails -- a genuine,
  instructive failure mode, archived here), and the fast soliton made O(dt^2) time-
  differentiation bias dominate. Data-excitation fix, criteria unchanged.
- exp10 rollout horizon for KdV shortened to 0.8 (two-soliton collision completes by
  then on the periodic domain; avoids wrap-around ambiguity in the comparison).
- exp11 (2026-08-26, after first run): excitation clipping tightened +1.5/-1.5 ->
  +0.4/-0.4 and library extended to degree 3. Reason archived: the original torque
  spun the pendulum through 32 radians, where the polynomial premise is vacuous
  (with sin(th) added to the library the SAME run recovers -0.9998*sin(th),
  -0.1499*om -- pipeline correct, design wrong); bounded-angle operation is the
  regime the experiment intended to probe.

