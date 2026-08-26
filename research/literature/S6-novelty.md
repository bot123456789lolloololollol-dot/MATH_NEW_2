# S6 Novelty Audit — Equation Discovery Study

Audited 2026-08-26 (all sources accessed this date). Verdicts use the lab claim
vocabulary; "blocking prior work" means the idea as stated is already published.

## Method components

| Component | Closest prior work | Verdict |
|---|---|---|
| Sparse regression discovery of ODEs from data (STLSQ + library) | Brunton, Proctor & Kutz, PNAS 113(15):3932 (2016), arXiv:1509.03580; PySINDy (Kaptanoglu et al., JOSS 2022) | **Known** — we implement and validate it; no novelty claimed |
| PDE discovery from snapshot data | Rudy, Brunton, Proctor & Kutz, Science 357(6355):940 (2017), "Data-driven discovery of partial differential equations" (PDE-FIND); Schaeffer, SIAM J. Sci. Comput. 39(3):A1234 (2017) | **Known** — exp10 is a replication; our archived traveling-wave-degeneracy failure (P9) is an instructive illustration of their stated sampling requirements, not new theory |
| STLSQ convergence/threshold theory | Zhang & Schaeffer, SIMADS 17(3):948 (2019), arXiv:1805.06445 | **Known** — our P1/P2 restates their noiseless analysis in entrywise form |
| GP symbolic regression | Koza 1992 (MIT Press); Schmidt & Lipson, Science 324:81 (2009); Cranmer, arXiv:2305.01582 (PySR) | **Known** — our GP is a small teaching-grade implementation of a classic idea |
| Conserved quantities via nullspace/SVD on feature increments F(s_{k+1})−F(s_k)=0 | **AI Poincaré**: Liu & Tegmark, PRL 126, 180604 (2021), arXiv:2011.04698 (same increment formulation, solved by LP); Liu et al., PRE 109, L023301 (2024), arXiv:2305.19525 (SVD nullspace on feature matrix); Ha & Jeong, PRR 3, L042035 (2021) (neural); Zhu et al., PRE 108, L022301 (2023) (neural deflation) | **Blocking prior work for the method itself.** Our exp07 is therefore labeled a *replication + calibration* of this known approach class on two new test systems (Duffing double-well, planar central force with non-polynomial energy), NOT a novel method |
| Symbolic zero-derivative certificate (∇F·f̂ ≡ 0 machine-checked) | Not found in this exact certificate form in items above (Schmidt–Lipson report recovered Hamiltonians but no symbolic identity checker artifact) | Small methodological nicety; labeled an engineering contribution, not a discovery |
| Validation-split threshold selection against correlated FD error | Standard cross-validation practice; not found specifically applied to SINDy threshold dust in cited works | Engineering refinement, disclosed in PREREGISTERED deviations |
| Significance-test support detection (\|c\| > z·SE) | Standard OLS inference; its use as the SINDy support criterion is unremarkable | Engineering |
| Predictability-horizon-aware validation for discovered chaotic models | Pathak et al., Chaos 27, 121102 (2017) established valid-prediction-time/Lyapunov-time normalization; Gilpin, NeurIPS 2021 (arXiv:2110.05266) benchmarks horizons in Lyapunov units | **Known framing** — our A4 twin-experiment decomposition (attributing horizon differences to log(ε-ratio)/λ) is a textbook consequence, used here as evaluation hygiene |

## Physics relationships recovered

Kepler's third law (+ eccentricity invariance), the amplitude expansion of the
nonlinear pendulum period, RLC damping/dispersion formulas, Duffing Hamiltonian,
angular momentum conservation: all classical, proven results — see THEORY.md P5–P8.
The scientific value is the demonstration that the *discovery pipeline* recovers them
under preregistered held-out/noise/adversarial protocols, together with honest failure
modes (GP's failure to find exact root forms; detectability boundary at tiny damping;
tanh domain-transfer failure). These demonstrations are reproducible artifacts, not new
physics claims.

## Overall novelty position

No component of this session is claimed as new science. The deliverable is:
(1) independent replication of the SINDy/invariant-discovery program with full
regenerability; (2) calibrated negative controls and adversarial gates (generalization
gate, compactness + domain-transfer bars) documented so others can adopt them;
(3) proofs (P1–P8) collected with certificates tying discovered equations to known law.
Per PROTOCOL scoring anchors this is Novelty 1–2 (known result / incremental),
which we self-report accordingly.
