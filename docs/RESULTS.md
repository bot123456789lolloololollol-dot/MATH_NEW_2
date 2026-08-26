# RESULTS — Regime-Switching Best Fit (RSBF) and the discovery campaigns

## Headline

**RSBF (Regime-Switching Best Fit)** is a new online bin-packing policy that
(i) is *provably identical* to Best Fit until a causally-computed trigger
fires, and (ii) after the trigger fires on streams whose past-item
distribution matches the "mid-band/triplet" regime, switches to the
tight-fit-threshold variant TBF(tau).  Across 4,800 paired evaluation
streams drawn from 16 distributions (12 synthetic families incl. two tuned
only for baselines + OR-Library Falkenauer u/t streams + Weibull tails +
4 held-out families), RSBF:

* **significantly beats Best Fit** exactly on the structured mid-band regimes:
  - `triplet_like`      102W/156T/42L, Wilcoxon p = 3.3e-06
  - `orlib_t`           43W/243T/14L,  p = 1.2e-04   (real benchmarks)
  - held-out `triplet_like_cap333` (capacity never seen in tuning):
    117W/132T/51L, p = 3.5e-08
* is **bit-for-bit BF elsewhere** (0–2 differing streams out of 300 per
  family; zero significant regressions anywhere), and
* resists a dedicated hill-climbing adversarial search: 12 restarts could not
  find ANY stream where RSBF uses even one extra bin vs BF.

Effect sizes are small but real: ~0.1–0.2% fewer bins on affected families,
consistent with the theory that Best Fit's average case is near-optimal on
i.i.d.-like streams and its losses concentrate in structured regimes.

**Scale validation**: on 40 fresh streams of n=10,000 items (triplet_like),
RSBF saves 169 bins vs BF (+0.099%, Wilcoxon p = 3.8e-06); on uniform
n=20,000 and discrete_3 n=10,000 it is exactly BF (0 differing streams) —
no false-trigger cost grows with stream length.

**Tuning stability**: re-running the hyper-parameter search under three
different tuning seeds yields warmup=20 and theta_mid=0.8 every time, with
(eps, tau) jittering inside {(0.05,0.02), (0.05,0.05), (0.15,0.05)}; all
variants in this neighbourhood reproduce the significant-regime-win /
zero-regression pattern.

Full table: experiments/RSBF_table.md; raw paired stats:
experiments/RSBF_results.json.  Tuning (seed 555) and evaluation (seed
20260826) use disjoint stream draws; hyper-parameters (warmup=20,
theta_mid=0.8, eps_small=0.15, tau=0.05) were chosen on tuning draws only.

## What RSBF is

```
state: residual list; causal counters of past items by size band
loop over arriving items s:
    if not yet switched and #past >= warmup and
         frac_past( cap/4 < size <= cap/2 ) >= theta_mid  and
         frac_past( size <= cap/4 )        <= eps_small:
         switch mode -> TBF(tau)          # monotone, no hysteresis
    place s:
        mode BF : tightest feasible bin, else new bin          (= Best Fit)
        mode TBF: if tightest fit <= tau*cap -> that bin;
                  else -> feasible bin with LARGEST post-residual
```

Properties proved / verified (docs/THEORY.md):
1. **Pre-trigger equivalence to BF** — by construction (identical code path).
2. **Post-trigger any-fit** — TBF always places into some feasible bin when
   one exists, so the whole run remains an Any-Fit algorithm and inherits
   the classical AF worst-case family bounds.
3. **Trigger safety** — firing requires >=80% of all past items inside a
   width-(cap/4) band and <=15% small items; on every tested non-matching
   family the trigger either never fires or fires without changing any bin
   count (0–2 streams of 300 differ, all wins/ties).
4. **Adversarial robustness** — hill-climbing over integer streams finds no
   regression instance.

## Why this is a fair contribution

* FunSearch (Nature 2023) and EoH (ICML 2024) evolved single fixed priority
  functions and report empirical-only gains vs BF/FF.  RSBF instead changes
  the *policy architecture* (causal regime detection + switching), comes with
  significance testing across 16 distributions incl. real benchmark streams
  and held-out capacities, an explicit safety analysis, and an adversarial
  probe.  The mechanism is simple enough to adopt in practice (O(1) extra
  state, few comparisons per item).
* We also contribute the negative results below, which justify WHY wins here
  are small and where future search should look.

## Negative results (both campaigns)

### Offline constructive/repair space is saturated (campaign 1)
* GP over placement-scoring functions (20 features, exact-fit & remaining-
  multiset statistics), island model, paired fitness vs min(FFD,BFD):
  converged to baseline-equivalent programs; nothing beat the classical
  ensemble on any of 900 standard instances or generated families.
* Hand-designed lookahead rules (exact-fit priority, dead-residual avoidance,
  cnt-fit weighting): all tie or lose.
* Monotone repair (unload-smallest-bins+redistribute+BFD, provably never
  worse than start; Theorem 1 in docs/THEORY.md): ZERO improving moves found
  on ALL standard instances — min(FFD,BFD) packings are locally optimal under
  this move class; a global low-load-rebuild variant also fails.
* Consistent with literature (MFFD/DJD/metaheuristics needed for progress).

### Online scoring-rule space around BF is near-saturated (campaign 2)
* GP with causal features (incl. past-histogram statistics), seeded with
  strong templates (BF/WF/exact-priority/dead-resid/TOF(K)): plateaued at
  BF-equivalence on aggregate fitness.
* Target-residual policies argmax -|fit - tau*cap|: strictly worse for
  tau > 0 everywhere tested.
* Tight-or-first TOF(K): ties BF on mid-band families.
* Fixed-threshold TBF(tau): wins up to 0.3% on mid-band regimes BUT loses
  similar amounts on heavy-tail/uniform families — the trade-off that
  MOTIVATES RSBF's switching design.
* EFF vs FFD offline micro-study: incomparable algorithms (90:1 win ratio).

## Reproduction

    python -m bpp.final_eval            # regenerates table + JSON from seeds
    python -m bpp.final_eval_adversarial  # (see experiments/) adversarial probe
