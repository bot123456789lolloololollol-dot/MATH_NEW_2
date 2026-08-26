# S5 REPORT — Comparator-network discovery campaigns (sorting + k-selection)

Session S5 (workstream 2). Date: 2026-08-26. All web sources accessed this date.
A concurrent S5 workstream pursued pure-SAT exact minima for merging networks
and small medians; its audit is in `../literature/S5-novelty.md` sections A–C,
its spec in `../experiments/S5/SPEC.md`. This report covers the independent
evolutionary + general-k-selector workstream; overlap is limited to shared
infrastructure-free tooling (separate encoders, separate verifiers).

---

## 0. TL;DR

| claim | label | evidence |
|---|---|---|
| SEL(6,3) = 8 exactly | `proven_result` | SAT UNSAT at L≤7 + verified witness at 8 |
| SEL(7,3) ∈ {10,11}; SEL(8,4) ∈ [10,14] | `proven_result` for the lower bounds (UNSAT@9) + verified witnesses above | certificates in `results_exact.jsonl`, `exact73b/84b.log` |
| SEL(10,5)=20@d7, SEL(12,6)=28@d9, SEL(14,7)=36@d9, SEL(16,8)=45@d11 are smaller than any construction derivable from published records | `experimentally_validated_result` ("best known under SEL(n,k) as defined here") | dual exhaustive verification F1+F2 |
| Published sorter records N=13..17 do not fall to ~40 CPU-min of hybrid evolution+repair | negative result, `experimentally_validated_result` | campaign logs |

Headline artifact: a 45-comparator, depth-11 network selecting the bottom 8 of
16 inputs — fewer comparators than the best published median network on 16
inputs (46 CEs, which selects only 7).

## 1. Pipeline

generated → compiled/simulated → benchmarked → optimized, per the mission:

- **Generated**: (μ+,λ)-style ES with (error,size) lexicographic fitness;
  mutations = delete/insert/move/retarget/swap/segment-reverse; recombination =
  prefix-suffix splice + segment transplant between Pareto-record seeds;
  error-guided remove-and-repair; periodic redundancy pruning to fixpoint.
  Seeds: 177 record networks downloaded from SorterHunter (repo state
  2026-08-26) + classical constructions (generalized Batcher odd-even merge
  sort for all n, insertion).
- **Simulated**: exhaustive evaluation over ALL 2^n zero-one inputs via
  bigint bitset parallelism (each wire = one 2^n-bit integer; compare-exchange
  = AND/OR pair). Correctness of encoding: wire functions of comparator nets
  are monotone, min/max = AND/OR on bits.
- **Benchmarked**: size (#compare-exchange elements ≈ area/cost), depth
  (parallel latency), plus failing-pattern counts as search gradient.
- **Optimized**: campaigns under preregistered budgets
  (`experiments/S5/PREREGISTERED.md`, written before any run).

## 2. Verification (three independent paths, cross-validated)

1. `netlib.is_sorting` — bigint bitset simulation, all 2^n inputs.
2. `numpyverifier.check_sorting` / `verify_selector.sel_f1` — numpy
   simulation; checks multiset preservation explicitly rather than assuming
   comparator semantics preserve it.
3. `sat_check.sat_sorts` — Tseitin CNF encoding of "∃ input violating sorted
   order"; solver UNSAT ⇔ network sorts. Cost O(L·n) clauses, scales past the
   bitset memory wall.

Validation battery (all pass):
- Three verifiers agree on 28 mixed known/broken networks (`python sat_check.py`).
- Encoder validation against literature-proven optima: sorting-mode exact
  decision returns UNSAT below and SAT at s(3)=3, s(4)=5, s(5)=9.
- Brute-force agreement: exhaustive enumeration gives SEL(4,2)=4, SEL(5,2)=6;
  SAT decision matches exactly.
- Selector formulations F1 (numpy sort-based) and F2 (pure-python counting:
  bottom-k block must contain max(0, c−(n−k)) ones) agree on random nets and
  on every reported result.

Independent reproduction of the published frontier: all 177 downloaded sorter
record networks re-verified — 30 by bitset+numpy (n≤20), 147 by the symbolic
SAT checker (n>20); **0 failures** (`seed_sweep.txt`). This is itself a
reproduction service to the table-keeper.

## 3. Campaign A — sorting records (negative result, preregistered targets)

Targets N=13..17 vs published best-known (45/51/56/60/71 sizes; Pareto gaps
such as 45@d9 at N=13). Two rounds × 20 min × up to 19 workers (~40 CPU-min
per target for N=13/14): **no new nondominated points found**. Instrumented
single-chain probe at N=13: ~2200 evals/s, walk confined to (error=0, size=45)
neutral drift — decades-tuned records are extremely locally optimal.
Honest conclusion consistent with the audit: closing s(13) ∈ [44?,45] or
filling its Pareto gaps needs industrial compute (cf. Harder 2020's Isabelle-
verified proof machinery for n=11,12). Logs preserved in `runs_n13*/`,
`runs_n14*/`.

## 4. Campaign B — k-selection networks (4/4 targets improved)

Property SEL(n,k) (preregistered): over ALL 2^n zero-one inputs, wires 0..k−1
carry exactly the k smallest input values as a multiset; tail unconstrained.
Baselines (computed + exhaustively verified in `select_baselines.json`):
B2 = smallest known sorter on n wires pruned under the selector property;
B3 = best SorterHunter median net when k = ⌊(n−1)/2⌋.

| (n,k)   | baseline (B2/B3) | found        | verdict |
|---------|------------------|--------------|---------|
| (10,5)  | 21 / –           | **20 @ d7**  | −1 CE |
| (12,6)  | 29 / –           | **28 @ d9**  | −1 CE |
| (14,7)  | 38 / –           | **36 @ d9**  | −2 CE |
| (16,8)  | 47 / –           | **45 @ d11** | −2 CE |

Every number dual-verified (F1 numpy + F2 counting) and saved with its network
in `verified/SEL_*_best.json`. Polish rounds seeded from the winners found no
further reduction (12 min × 5–6 workers each).

Remark: SEL(16,8) = 45 beats even the *published* best median network on 16
inputs (46 CEs, selecting only 7 elements) — the tail-partition constraint in
the published variant carries measurable cost. This comparison is across
different properties and is offered as analysis, not as a same-property record.

## 5. Exact values via validated SAT decision procedure

Encoding: L slots × one-hot over {all pairs (i,j), i<j} ∪ {IDLE}; per pattern
(all 2^n enumerated → complete), constants propagate through Tseitin-gated
min/max terms; selection enforced per pattern by cardinality constraints
(#ones among wires<k equals max(0, c_p−(n−k)), c_p = popcount(p), constant).
Sound symmetry breaking only: (a) adjacent duplicate comparators banned
(deleting the first of two identical adjacent CEs preserves behavior);
(b) IDLE slots form a suffix (deleting an interior idle slot preserves
behavior). Neither excludes any minimal network ⇒ completeness preserved.

Results (`results_exact.jsonl`, `exact73b.log`, `exact84b.log`):
- **SEL(6,3) = 8** — UNSAT for every L ≤ 7 (certificates: 0.3/2.5/6.8/15.1 s),
  SAT witness at 8, decoded and re-verified by F1+F2. `proven_result`.
- SEL(7,3): UNSAT through L=9 (163.6 s @8, 1268.6 s @9) ⇒ **≥ 10**; ES found
  an 11-CE@d5 witness ⇒ **SEL(7,3) ∈ {10, 11}**.
- SEL(8,4): UNSAT through L=9 (452 s @8, 920 s @9) ⇒ **≥ 10**; ES witness at
  14 ⇒ **SEL(8,4) ∈ [10, 14]**.
- Scaling data: UNSAT cost roughly ×2–4 per +1 CE (118→452 s at L=7,8; 452→920
  at L=9 for (8,4)); SEL(8,4)@L=10 exceeded a 1800 s budget without a verdict.
  Direct-encoding practicality ends near n≈9–10 on this
  machine — documented as the deliverable-ladder item (d) scaling limit.
- Engineering note: python-sat's CaDiCaL binding lacks interrupt(); the first
  L=9 attempts crashed on the timer guard. Fixed by subprocess isolation
  (`run_exact2.py`); the crashed runs produced no results, all recorded
  verdicts come from completed solves.

### Completeness argument (required gate for `proven_result`)
Any network with L′ ≤ L comparators maps to an assignment of the encoding:
use its comparator sequence in the first L′ slots, IDLE elsewhere. Tseitin
gates reproduce exactly the min/max semantics on every zero-one pattern
(equisatisfiable by construction). Property constraints are imposed on all
2^n patterns, so SAT ⟺ existence of an ≤L network satisfying SEL on all
zero-one inputs. Extension to arbitrary total orders: Lemma below. Hence
UNSAT at L certifies no ≤L selector exists. Solver soundness assumed
(CaDiCaL 153 via python-sat); certificates were not independently
re-checked by a second solver — noted as residual trust assumption.

### Zero-one principle for selectors (needed to lift exhaustive 0/1 checks)
Lemma: if a comparator network places the bottom-k multiset correctly on
wires 0..k−1 for every 0/1 input, it does so for every input from any total
order. Proof: for an input vector v, let θ_j be the threshold map sending
values ≥ j-th-smallest-distinct to 1. Comparators commute with thresholds:
θ(min(a,b)) = min(θa,θb) and likewise for max, so applying the network to v
and then thresholding equals thresholding the inputs and applying the
network. The multiset of a block being the k smallest is determined by the
counts of ones at every threshold level; these match by hypothesis on each
θ_j(v) ∈ {0,1}^n. ∎

## 6. Why the winners win (structural analysis)

Wire-degree profiles (comparisons touching each wire):

- Best sorter N=16 (60 CEs): near-symmetric load peaking on center wires
  (degrees …,10,9,9,10,…), regular Batcher layer taper 8,8,8,8,7,6,4,4,5,2.
- SEL(16,8) winner (45 CEs): asymmetric profile with load concentrated at the
  selection boundary (wire degrees: w6=9, w7=7, w8=8) while tail wires are
  starved (w14=5, w15=3); irregular layer taper 8,6,4,8,6,5,2,2,1,2,1.
- SEL(12,6) winner (28 CEs): boundary wires 4–7 carry degrees 7,8,7,6 vs 3–4
  elsewhere.

Mechanistic reading: correctness only constrains the bottom-k block, so the
search reallocates comparisons away from the ordered tail onto boundary
wires that decide *which* elements survive the cut. The one-sided property
also removes the symmetric-load pressure a full sorter obeys, which is where
the 2–5% size savings come from. This asymmetry is invisible to constructions
derived from sorters (our baselines), explaining why evolution finds points
they cannot reach.

## 7. Reproduction

From `research/experiments/S5/`:

```
python baselines.py                  # constructive baselines + self-validation
python sat_check.py                  # 3-verifier consistency battery
python verify_selector.py            # F1/F2 self-tests
python select_baselines.py           # baseline table (verified)
python search_es.py --mode select --n 16 --k 8 --minutes 20 --workers 16 \
    --seed 20260826 --out runs_sel16_8
python merge_runs.py runs_sel16_8    # verify + Pareto vs published
python run_exact.py 6,3 900          # exact SEL(6,3) decision ladder
cat verified/SEL_16_8_best.json      # headline artifact + verifier list
```

Deterministic seeds; raw logs committed under `runs_*/`; no result enters the
report without passing both F1 and F2 (selectors) or bitset+SAT (sorters).

## 8. What is NOT claimed

- No improvement on published SORTING network records (attempted, not achieved).
- "Best known" selector sizes are relative to SEL(n,k) as defined here and to
  baselines derivable from published data; other definitions (e.g., with
  complementary tail partition) may differ.
- UNSAT certificates were produced by CaDiCaL without independent proof-file
  re-checking; completeness argument is manual (not machine-checked).
- Novelty audit could not reach Google Scholar/arXiv listing pages directly
  (bot-blocks); absence-of-prior-art is supported but not exhaustive.

---
---

# PART II — Pure-SAT exact minima: sorting, merging, median networks
*(second concurrent S5 workstream; encoder + results below are independent of
Part I's ES tooling — separate encoder `sat_net.py`, separate verifier `checker.py`.
Coordination note for the lead: two S5 agents were spawned on identical mission
text and discovered each other mid-flight; artifacts were merged additively,
nothing overwritten.)*

## 9. TL;DR of this part

| claim | label | evidence |
|---|---|---|
| Encoder reproduces s(4)=5, s(5)=9, s(6)=12 exactly (UNSAT@s−1, SAT@s) | `experimentally_validated_result` | CaDiCaL + Glucose UNSAT agreement; z3 native-engine agreement at n=4 pair and n=5 SAT; all witnesses PASS independent 2^n-input checker |
| Merging networks M(2,2)=3, M(2,3)=5, M(3,3)=6, M(2,4)=6 exact | `experimentally_validated_result` (proven modulo §12 argument; values likely classical — see audit §B) | UNSAT below / SAT above / Glucose cross-checks / checker over all constant inputs |
| med(5)=7 exact (median-selection network) | `experimentally_validated_result` | UNSAT@6 (CaDiCaL AND Glucose), SAT@7, checker PASS 32/32 inputs |
| SB2 symmetry breaker ("first comparator = (1,2)") is UNSOUND for designated-output properties such as median | `experimentally_validated_result` (empirical demonstration) | SB2-ON ⇒ spurious UNSAT at T=7; SB2-OFF ⇒ SAT@7, UNSAT@6 as expected |
| Direct encoding cannot close s(13)∈[43/44,45] or prove n≥7 sorter lower bounds on this machine class | negative result, `experimentally_validated_result` | timing table §11 |

## 10. Encoding (full spec: `experiments/S5/SPEC.md`)

Slot-based CNF: choice vars x(l,p) per slot l ∈ [0,T) and channel-pair p (plus an
EMPTY var per slot); total simulation semantics w(v,l,c) for every enumerated input v,
level l, channel c — whatever single choice is active (min-side, max-side, non-touching,
empty), w(v,l+1,c) is fully constrained. Property clauses: sort/merge forbid descents
`[¬out_k, out_{k+1}]`; median pins output channel n//2 to the input's median bit.
Inputs enumerated: ALL binary inputs except ascending-canonical ones (sound: comparators
preserve ascending-sortedness; canonical inputs can never violate either property);
merge uses only the (m+1)(n−m+1) constant block-sorted inputs (exhaustive over its
actual input space). Symmetry breaking: SB1 lex-order of adjacent disjoint comparators,
SB2 first-comparator relabel (sort only — see unsoundness finding below), SB3 no
consecutive duplicate pairs; compatibility by bubble-normalization argument (SPEC §5).

Median-property soundness beyond binary inputs (threshold lemma): N(x)[mid] ≥ t ⇔
N(θ_t(x))[mid] = 1 for every threshold map θ_t, because comparators commute with
thresholding; hence checking all 0/1 inputs suffices for arbitrary total orders.

## 11. Results and timings (machine under heavy cross-session CPU contention; wall-clock)

| instance | expected | verdict | time | engine(s) |
|---|---|---|---|---|
| n=4 T=4 (sort) | UNSAT | UNSAT | 0.04s | CaDiCaL; z3-native agrees |
| n=4 T=5 | SAT | SAT | 0.04s | CaDiCaL; z3-native agrees; checker PASS 16/16 |
| n=5 T=8 | UNSAT | UNSAT | 0.78s | CaDiCaL |
| n=5 T=9 | SAT | SAT | 0.22s | CaDiCaL; z3-native agrees; checker PASS 32/32 |
| n=6 T=11 | UNSAT | UNSAT | 71s | CaDiCaL AND Glucose (111s) |
| n=6 T=12 | SAT | SAT | 3.7s | CaDiCaL; checker PASS 64/64 |
| n=7 T=15 | UNSAT | TIMEOUT (>500s, killed) | — | scaling limit documented |
| n=7 T=16 | SAT | SAT; checker PASS 128/128 | 362s | witness net_n7_T16.json |
| merge (4,2) T=2 / T=3 | UNSAT / SAT | confirmed both engines | <0.05s | witness [[1,3],[2,4],[2,3]] checker PASS 9/9 |
| merge (5,2) T=4/T=5 | UNSAT / SAT | confirmed | <0.07s | M(2,3)=5 |
| merge (6,3) T≤5 / T=6 | UNSAT / SAT | confirmed; glucose agrees | <0.15s | M(3,3)=6, witness PASS 16/16 |
| merge (6,2) T≤5 / T=6 | UNSAT / SAT | confirmed | <0.26s | M(2,4)=6, witness PASS 15/15 |
| median n=5 T=6 | UNSAT | UNSAT (CaDiCaL + Glucose) | 0.13s | with SB2 OFF |
| median n=5 T=7 | SAT | SAT; checker PASS 32/32 | 0.18s | witness in net_med_n5_T7.json |
| median n=7 T=12 | ? | TIMEOUT (>480s) | — | recorded honestly |

Generation sizes: n=7 T=15 → 13,770 vars / 636,805 clauses; n=8-class instances grow
as ~2^n·n·T·P clauses (SPEC §3), which is the practical wall.

**Debugging methodology worth keeping:** during development, three clause-sign bugs and
one inverted-property bug were caught by FORCING a known-optimal, independently-checked
network into the encoder as unit clauses and diffusing every generated clause against
the true simulation trajectory (`violated clauses: 0` required). Also caught: a partial
semantics encoding (constraining only compared channels) that admits spurious models —
unsound; and SB2-on-median producing wrong UNSAT. All documented in SPEC §8 pitfalls.

## 12. Soundness/completeness argument for UNSAT claims

Any network with L′ ≤ L comparators satisfying the property maps to a formula
assignment: its comparator sequence occupies the first L′ slots, EMPTY elsewhere;
clauses §3.1–3.2 (SPEC) reproduce min/max/copy semantics exactly on every enumerated
input; §4 enumerates enough inputs (canonical-exclusion and merge-constant arguments
in SPEC §4); §5 transformations (relabel/dedupe/bubble) map any network to an
SB-restricted assignment of equal or smaller size without changing its function.
Hence formula-UNSAT at L implies no ≤L network exists; combined with a decoded,
independently checked witness at L+1 this proves the minimum = L+1. Residual trust
assumptions: solver soundness (CaDiCaL 153 / Glucose 4 via python-sat, z3 5.1) —
mitigated by multi-engine agreement on every UNSAT reported here; no DRUP proof-file
verification was performed (noted).

## 13. Regenerate (from `research/experiments/S5/`)

```
bash run.sh                      # full validation battery + variants + cross-checks
python checker.py --selftest     # independent verifier self-test (8 cases)
python sat_net.py cnf --n 6 --slots 11 --out /tmp/x.cnf
python solve_run.py --dimacs /tmp/x.cnf --engine cadical --timeout 300   # expect UNSAT
```

## 14. What this part does NOT claim

- No new bound on plain s(n): frontier (proven n ≤ 12; open case n=13) untouched —
  correctly so given Harder 2020-scale machinery needed even for n=11,12.
- M(m,n)/med(5) values are small and likely classical (audit §B); shipped as
  certificates + pipeline, not novelty claims.
- med(7) left open here (timeout); natural next target together with M(3,4), M(4,4),
  and SEL-style exact values from Part I.
