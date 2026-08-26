# Lab Protocol

## Discovery lifecycle

```
Idea → implementation → benchmark → adversarial testing → improvement
     → novelty search → proof → independent reproduction → final evaluation
```

A candidate may enter `/discoveries` only after: (a) novelty audit finds no blocking prior work,
(b) proof artifact exists at the claimed level, (c) an agent that did not author it reproduces the
key evidence from the written spec alone.

## Evidence standards

- **Exact/optimal claims**: require a deterministic, rerunnable exhaustive search or a certificate
  (explicit object + independent checker script), plus a written completeness argument
  (why the search space reduction/symmetry breaking loses no solutions).
- **Best-known/record claims**: require explicit object + checker + statement of what is NOT proven.
- **Algorithmic improvement claims**: fixed seeds preregistered before results are seen; ≥30 runs;
  report mean ± std and a paired significance test; baseline implemented from its primary source.
- Every experiment dir must contain `run.sh` (or `run.py`) with the exact command to regenerate outputs.

## Scoring rubric (each factor 0–5, product = total, max 3125)

| Factor | 5 | 3 | 1 |
|---|---|---|---|
| Novelty | No prior art found after serious audit; fills a stated gap | Increment on known results | Known result |
| Usefulness | Practitioners/table-keepers would adopt | Niche but real value | Toy |
| Provability | Proof complete or clearly within reach now | Partial (bounds/constructions only) | Far out of reach |
| Reproducibility | One-command rerun from spec alone; independently reproduced | Rerunnable by authors | Not rerunnable |
| Improvement | New exact value / >10% robust gain / new record | Modest gain (~2–10%) | ≤ noise |

## Claim vocabulary

`observation` (noticed, unverified) / `hypothesis` (stated falsifiably, weakly checked) /
`conjecture` (checked in many cases, unproven) / `experimentally_validated_result` (all claims
regenerable, adversarially tested) / `proven_result` (rigorous proof incl. exhaustive-search
completeness arguments) / `formally_verified_result` (machine-checked proof, e.g. Lean 4).

## Coordinator decision log

Material decisions (promotions, kills, mission changes) go to `coordination/DECISIONS.md`.
