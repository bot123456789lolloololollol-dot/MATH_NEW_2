# S6 Preregistration — ROUND 4 addendum (micro-study)

Written before implementation/run (2026-08-26 ~06:17 local).

## R4-exp14 — does constant optimization fix the root-form miss?

- Motivation: exp03/exp04 documented that vanilla GP misses short root-form
  expressions (e.g., c/sqrt(g)) within budget while finding other exact forms.
  Standard remedy in the SR literature (used by PySR et al.): periodic numerical
  optimization of ephemeral constants. We test exactly this one mechanism.
- Tasks (fixed): T1 y = c/sqrt(x), x~U[0.5,2], c=6.2832 (the exp04 miss);
  T2 y = c/(1+x^2), x~U[-2,2], c=2.0 (Lorentzian -- rational root form);
  n=500 train + 300 held-out per seed.
- Methods (fixed):
    base : SymbolicRegressor(population=600, generations=150, kappa=1)
    opt  : identical + constant refinement every 10 generations on the top-10
           cost individuals, 25 iterations of cyclic coordinate descent with
           backtracking step (numpy-only; no external optimizers).
- Endpoints (fixed): success rate over 30 seeds (holdout NMSE < 1e-6 AND numeric
  equivalence to ground truth within 1e-6 on fresh grid), median holdout NMSE,
  median nodes; paired two-sided Wilcoxon on holdout NMSE, alpha=0.05.
- Either outcome is publishable within the session: positive -> tool upgraded and
  exp04 note amended; negative -> vanilla search vindicated at this scale, note kept.
