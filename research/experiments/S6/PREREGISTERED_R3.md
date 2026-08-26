# S6 Preregistration — ROUND 3 addendum

Written before any round-3 experiment was run (2026-08-26 ~05:20 local).

## R3-exp12 — hidden modal structure in coupled oscillators

- System (fixed): two coupled linear oscillators
    q1'' = -k1 q1 - kc (q1 - q2)/m1 ;  q2'' = -k2 q2 - kc (q2 - q1)/m2,
  parameters k1=4, k2=9, kc=1.5, m1=m2=1 (analytic normal-mode frequencies available
  in closed form via the 4x4 system matrix eigenvalues).
- Data: TWO training trajectories from generic ICs (seeded); learner sees states only.
- Discovery: identify x' = A x (linear SINDy degenerate case) by least squares +
  significance support; then compute eigenvalues/eigenvectors of discovered A.
- Claims tested (fixed criteria):
  C-i: A recovered entrywise within 1e-8 clean;
  C-ii: discovered eigenfrequency pairs match analytic sqrt(eig(A)) within 1e-6 rel;
  C-iii: eigenvectors (mode shapes) match analytic up to scale within 1e-8 after sign
         alignment -- this is "hidden structure" not present as explicit terms in any
         textbook description of the data;
  C-iv: held-out trajectory from unseen IC predicted with rel err < 1e-6.
- Noise extension: sigma in {0, 0.005, 0.01}; report degradation honestly.

## R3-exp13 — integral-form identification against the sigma=5% failure

- Motivation: exp02 documented that derivative-based STLSQ fails at LV sigma=0.05
  (rollout divergence). Known alternative family (integral/weak form) avoids
  differentiating noisy states by integrating the equation against test functions.
- Implementation (fixed): for LV, integrate dX over sliding windows of length W=25
  samples via cumulative sums: target becomes X(t+W)-X(t) ≈ Theta-integrated c.
  Same library; same validation-split threshold selection; z=5 significance.
- Protocol: sigma in {0.02, 0.05}, 30 paired seeds per level; primary endpoint =
  held-out rollout error < 5% (preregistered recovery criterion); comparison metric =
  rollout error ratio integral/derivative. Report whichever way it goes; a negative
  result (no rescue) is an acceptable outcome and will be labeled as such.
- Theory note: P12 in THEORY.md sketches why integration acts as a low-pass filter
  with noise variance shrinking as O(W^-1) relative to signal, versus O(dt^-2)
  AMPLIFICATION for differentiation.

## R3-misc

- SPEC.md extended to cover round-2/round-3 artifacts.
- After all experiments: full `python run_all.py` regression pass + `check_claims.py`
  must be green before the final push (guards against late edits breaking stored
  results).

## Post-hoc deviations (round 3)

- exp03 (2026-08-26, found by the regression pass): per-seed RNG derivation replaced
  `hash((case, seed))` -- whose string hash is randomized per Python process -- with
  `zlib.crc32(case)*1000+seed`. exp03 results were regenerated; qualitative rates
  unchanged in kind (S1 4/5, S2 3/5, S3 5/5 equivalent forms); output verified
  byte-identical across consecutive runs.
