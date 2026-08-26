# Automated algorithm discovery for bin packing

Research project: discover new or materially improved algorithms/heuristics
by **executable candidate generation + objective benchmarking** (no
prose-only invention).  Problem domain: 1-dimensional bin packing, both
offline and online, chosen because candidates can be evaluated millions of
times per hour and because strong baselines with known worst-case bounds
exist.

## Headline results so far

### Offline line (negative result + micro-findings)

1. **No cheap rule beats the classical ensemble.** Across ~1.6M evaluations:
   (a) a full genetic-programming search over placement scoring functions
   (20 features incl. exact-fit counters and remaining-multiset statistics),
   (b) hand-designed lookahead rules, and (c) monotone repair moves
   (unload-smallest-bins + redistribute + BFD; plus a global low-load-rebuild
   variant) — NOTHING beat `min(FFD, BFD)` on any instance family at greedy
   runtimes.  FFD/BFD solutions are locally optimal under these move classes
   on all 900 standard-suite instances (Falkenauer u/t, Scholl I & hard,
   Wäscher, hard28).  This is consistent with (and sharpens) the published
   record, where improvements over FFD always add machinery (MFFD reductions,
   DJD exact-fill search, metaheuristics).
2. **Monotone-repair wrapper with proven safety** (`bpp/repair.py`): output is
   provably never worse than its start (Theorem 1 in docs/THEORY.md), hence
   inherits Dósa's tight FFD bound 11/9·OPT+6/9 — but empirically it never
   fires on standard instances.  Committed as an honest negative result.
3. **EFF vs FFD incomparability**: exact-fit-priority First Fit wins vs FFD
   on ~1.7% of random instances and loses on ~0.018% (90:1 ratio, 2M-trial
   study); neither dominates.  Both directions counterexamples committed.

### Online line — RESULT: Regime-Switching Best Fit (RSBF)

Campaign 2 first established that the space of fixed any-fit scoring rules
around Best Fit is near-saturated: GP over 20 causal features plateaued at
BF-equivalence; target-residual and tight-or-first families tie or lose;
fixed-threshold TBF(tau) trades wins on mid-band regimes against losses on
heavy-tailed streams.  That trade-off motivated the discovered construction:

**RSBF** runs exactly Best Fit while causally tracking the past-item size
distribution; when it evidences a mid-band ("triplet") regime, it switches
permanently to TBF(tau).  Across 4,800 paired evaluation streams over 16
distributions (incl. OR-Library Falkenauer u/t and held-out capacities):

- significant wins vs BF on all three structured regimes
  (best p = 3.5e-08, incl. a capacity never seen during tuning),
- bit-for-bit BF behaviour everywhere else (zero significant regressions),
- provably identical to BF pre-trigger; Any-Fit worst-case bounds retained,
- no adversarial stream found (12-restart hill-climb) that costs even one
  extra bin.

Details: docs/RESULTS.md, docs/THEORY.md §3, experiments/RSBF_table.md.

## Repository map

    bpp/core.py            jitted offline kernels (FF/BF/WF/FFD/BFD/L1/L2/exact DP)
    bpp/gp.py              offline GP interpreter (verified == best_fit)
    bpp/online.py          ONLINE causal GP interpreter + NF/FF/BF/WF/Harmonic
    bpp/repair.py          monotone unload-and-repack (dominance theorem)
    bpp/baselines.py       DJD variants
    bpp/suites.py          loaders: Falkenauer/OR-Library, Scholl, Wäscher, hard28
    bpp/streams.py         online stream generators (+held-out families)
    bpp/tbf.py             tuned baselines TBF(τ), PFB, Harmonic tuning
    bpp/evolve.py          offline GP evolution (round 1)
    bpp/evolve2.py         offline islands + paired fitness (round 2)
    bpp/evolve_online.py   online policy evolution (main line)
    bpp/adversarial.py     adversarial instance search (offline)
    bpp/rsbf.py            Regime-Switching Best Fit + tuning
    bpp/final_eval.py      full paired evaluation protocol (tuning/eval seed split)
    bpp/adversarial_online.py  adversarial stream search (online)
    bpp/analysis.py        any-fit tracking / policy property checks
    bpp/stats_tools.py     Wilcoxon, Cohen's dz, bootstrap CI
    docs/THEORY.md         statements + proofs + references
    benchmarks/            standard suites (downloaded from OR-Library / ESICUP)
    experiments/           champions (*.npz), logs, evaluation JSON

## Reproducing

    python tests/sanity.py                 # kernels vs exact optimum (3000 cases)
    python tests/test_gp.py                # interpreter == Best Fit equivalence
    python -m bpp.suites                   # load all 900 standard instances
    python -m bpp.evolve_online --out experiments/onrunX --seed S
    python -m bpp.ood_online experiments/onrunX_champion.npz

## Environment notes

Windows, Python 3.11, numba JIT; multiprocessing capped by BPP_WORKERS env
var.  All experiments use fixed seeds; every champion is stored as
(ops, consts) arrays and is fully deterministic.
