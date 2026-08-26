# Session S6-MATH Report — Mathematical Discovery Researcher

Date: 2026-08-26. Mission: discover a genuinely new mathematical theorem via the
loop known-result → modification → conjecture → massive counterexample search →
proof → verification. Note: another concurrent session also used label "S6"
(physics/SINDy); this session's artifacts use the S6M prefix where global.

## TL;DR

Two candidates pursued through the full pipeline:

1. **C-S6M-001 (KILLED).** Independently discovered exact two-color Rado
   numbers for x+y=kz with pairwise-distinct variables — piecewise quadratic
   formula for k ≥ 8 with sporadic exceptions. Killed by novelty audit: A.
   Towell published the identical theorem (Zenodo 19372727 / OEIS A394445,
   April 2026). Post-mortem + independent confirmation archived. The kill is a
   process success: no rediscovery was labeled as discovery.

2. **C-S6M-002 (ACTIVE, headline).** New problem family and first-ever
   per-modulus classification: which cyclic groups Z_n admit 2-colorings
   avoiding monochromatic pairwise-distinct solutions of x+y ≡ kz (mod n)?
   - **Conjecture (validated):** for n ≥ 12 exactly three residue classes are
     good: k ≡ 0; k = n/2 (n even); k ≡ ±n/3 (9 | n). Finite exceptional
     tables for n ≤ 11; k=2 sub-case = Grier 2012 (cited, not claimed).
   - **Proven:** all three infinite construction families (Theorems A, B, C in
     sessions/S6/proofs/THEORY_C-S6M-002.md) — explicit colorings with full
     elementary proofs.
   - **Verified:** exhaustive Good(n,k) determination for every n ∈ [12,60]
     and every residue k (~1700 instances, zero mismatches); n ≤ 200 sweeps
     per fixed k ∈ [2,12]; three independent solver implementations agree;
     constructions verified to n ≈ 600–1350.
   - **Open:** the badness direction in full generality (roadmap §7 of theory
     doc). Status label honestly: experimentally_validated_result, not yet
     proven_result.
   - Novelty audit #2 (81 tool calls): no prior work for k ≥ 3.

## Process narrative

- Recon: fresh lab repo, protocol/roadmap present; Python 3.11 available,
  sympy absent, Lean absent, gh absent (no remote push possible — commits local).
- Tooling: built DFS backtracking solver + CaDiCaL SAT encoding + naive brute
  force (three independent implementations); validated on classical anchors
  S(2)=5, S(3)=14 before trusting anything.
- C-S6M-001 arc: computed strict Rado table k≤30 → found piecewise-quadratic
  law → verified k=8..30 → analyzed extremal colorings (residue-class structure
  mod k, odd-k reduction to a 1-D constraint problem) → novelty audit killed it
  (identical theorem published 5 months earlier). Salvaged as tool validation +
  independent confirmation.
- C-S6M-002 arc: pivoted to the cyclic variant (audit-clean gap). Data hunting:
  transposed the problem (bad residues per n) → conjectured gcd-shaped
  classification → exhaustive verification n=12..60 → designed and proved three
  construction families (half-split; two-block with parity case-split; layered
  column-constant with tuned special columns after two falsified design
  attempts) → adversarial audits on novelty and correctness.
- Failure honesty: naive divisor-lifting lemma falsified by data (Z_17 vs Z_34);
  EH-sumset badness lemma refuted by certified instances (Z_11,k=10); two
  candidate ±n/3 constructions falsified by verification before the third
  closed. All recorded.

## Regenerability (key commands)

All under research/sessions/S6/code:
- `python solver.py validate`                       # S(2)=5, S(3)=14 anchors
- `python sat_solver.py rado 30 2 weak`             # C-S6M-001 table (k<=30)
- cyclic sweeps + classification sweep + construction verification:
  see sessions/S6/data/*.json and the transcripts embedded in commits; scripts:
  solver.py, sat_solver.py, reduced_odd.py (all deterministic).

## Score self-assessment for C-S6M-002

novelty 4 (audited clean; crowded neighborhood), usefulness 3 (table-keepers /
Ramsey-on-groups niche), provability 4 (constructions proven; converse open but
tractable), reproducibility 5 (deterministic, triple-implemented), improvement 4
(first determination of an entire family). Product 960.

## Handoff to coordinator

- Promote C-S6M-002 to adversarial round: attack the badness direction at
  n ∈ [61, 200] and the small-n tables; attempt independent reproduction from
  THEORY_C-S6M-002.md spec alone.
- Push remotes absent; all evidence committed locally (see git log S6M).
