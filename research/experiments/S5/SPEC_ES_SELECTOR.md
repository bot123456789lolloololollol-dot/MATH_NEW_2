# S5 SPEC (workstream 2) — ES engine + exact selector SAT encoding

Companion to REPORT.md; specifies this workstream's artifacts precisely enough
for independent reimplementation. The pure-SAT merging/median encoder by the
other S5 workstream is specified separately in `SPEC.md`.

## A. Definitions

- Network on n wires: sequence of compare-exchange elements (i,j), i<j;
  element maps wire values (a,b) → (min,max) onto wires (i,j) respectively.
- Size = number of elements. Depth = greedy layering: an element's layer is
  1 + max(layer of its two wires before it).
- Sorting: outputs ascending for every input.
- SEL(n,k): for EVERY zero-one input pattern p, the multiset of values on
  wires 0..k−1 equals the k smallest inputs of p (tail unconstrained).

## B. Exhaustive verification semantics

All claims are checked over ALL 2^n patterns (complete), then lifted to all
total orders by the threshold-commutation lemma (REPORT.md §5).

Three independent implementations must agree:
1. bigint bitset (`netlib.py`): wire i ↔ 2^n-bit integer; CE = AND/OR pair;
   sortedness ⇔ wire k equals precomputed threshold mask T_k.
2. numpy (`numpyverifier.py`, `verify_selector.sel_f1`): explicit per-pattern
   arrays; multiset preservation CHECKED, not assumed.
3. SAT symbolic (`sat_check.py`): variables = input bits; Tseitin min/max
   gates per element; assert ∃ adjacent inversion; UNSAT ⇔ sorts.

Selector check F2 (`verify_selector.sel_f2`): independent counting formulation
— bottom-k block must contain exactly max(0, c−(n−k)) ones where c=popcount(p).

## C. Exact decision encoding (`exact_selector.py`)

Decision problem: ∃ network of ≤ L elements satisfying SEL(n,k)?

- Variables y[t,a]: slot t ∈ [0,L) selects action a ∈ {(i,j): i<j} ∪ {IDLE},
  exactly one per slot (CardEnc atleast+atmost).
- Per pattern p ∈ [0,2^n): wire entry values are unit-clause constants.
- Slot semantics per touched wire w: next(w) = keep(w) ∨ selected-contributions;
  keep = ¬selany ∧ cur[w]; low contribution = y∧cur[i]∧cur[j]; high =
  y∧(cur[i]∨cur[j]); all Tseitin (AND2/AND3/OR chains).
- Constraint per pattern: cardinality over output block literals:
  #ones(wires 0..k−1) = max(0, c_p−(n−k)) (atleast+atmost with equal bounds).
- Symmetry breaking (soundness argued in REPORT §5; both preserve every
  minimal network): no two adjacent slots select the same comparator;
  IDLE selections form a suffix.
- Sort-mode variant (for validation): replace block constraint by pairwise
  output order clauses; validated: UNSAT below / SAT at s(3)=3, s(4)=5, s(5)=9,
  matching literature-proven optima; brute-force match at SEL(4,2)=4, SEL(5,2)=6.

Driver `run_exact.py`: linear ladder from L=2; SAT → decode witness (one-hot
readout), re-verify via F1+F2 before recording; UNSAT recorded as bound;
per-L wall-clock interrupt guard (threading.Timer → solver.interrupt).

## D. ES engine (`search_es.py`)

- Fitness: lexicographic (error, size); error = # mismatched (pattern,block)
  entries for selectors, mismatched (pattern,wire) bits for sorters. Invalid
  intermediates are accepted when they reduce error ⇒ mobility across
  Pareto points.
- Operators: delete / insert random CE / move / retarget ±2 / swap /
  segment-reverse; crossover = prefix-suffix splice between pool members;
  transplant = copy 1..4-CE segment into current walk; remove-and-repair =
  delete random CE then greedily insert up to 30 candidate CEs keeping the
  best error reduction (early exit at 0).
- Pool: distinct valid nets, Pareto-sorted by (size, depth); neutral drift
  accepted with prob 0.5; periodic fixpoint pruning of pool leader.
- Parallelism: independent worker chains (distinct RNG seeds) appending JSONL
  improvement logs; merge step re-verifies everything (`merge_runs.py`).

## E. Artifacts

- `verified/SEL_{n}_{k}_best.json` — champions incl. full networks.
- `verified/campaign_B_winners.json` — the four improvements.
- `results_exact.jsonl` — exact-search ladder log incl. witnesses.
- `seed_sweep.txt` — reproduction status of all 177 downloaded record nets.
- `PREREGISTERED.md` — budgets/seeds/success criteria, written pre-run.
