# PREREGISTERED — S5 comparator-network campaigns

Written before any search was run (2026-08-26). Fixed now; not to be edited
after results are observed. Any deviation will be flagged in REPORT.md.

## Verifiers (fixed in advance)

- Primary: `netlib.is_sorting` — bigint bitset simulation over ALL 2^n
  zero-one inputs (zero-one principle links this to all total orders).
- Independent: `numpyverifier.check_sorting` — structurally different numpy
  simulation; additionally checks multiset preservation instead of assuming
  it. A result is claimed only if BOTH pass in fresh processes.

## Campaign A — sorting-network upper bounds

Targets (fixed): N in {13, 14, 15, 16, 17}.
Objective: lexicographic (size, depth). Reference best-known (Dobbelaere,
fetched 2026-08-26):

| N | best size | best depth | Pareto gap candidates |
|---|---|---|---|
| 13 | 45 (proven optimal) | 9 @46 | 45@9 unknown |
| 14 | 51 | 9 @52 | 51@9 unknown |
| 15 | 56 | 9 @57 | 56@9 unknown |
| 16 | 60 (Green 1969) | 9 proven @61 | 60@9 unknown |
| 17 | 71 | 10 @74, 11 @72 | 71@11, 71@10 unknown |

Success criteria (declared in advance):
- HIT-1: any network strictly dominating a published (size,depth) point.
- HIT-2: any unpublished nondominated point (e.g. 45@9) = new design.
- Otherwise: honest report of achieved (size,depth) vs baselines.

Budget: <= 20 min wall-clock per (N, seed-family) run; up to 4 runs per N;
seeds = 20260826 + run index. Workers = 16 processes.

## Campaign B — k-selection networks

Property: bottom-k inputs on wires 0..k-1 (unordered), tail unconstrained;
verified exhaustively by numpy checker (and cross-checked by a second
formulation: prefix-threshold consistency, see verify_selector.py).
Targets (fixed): (N,k) in {(10,4),(10,5),(11,5),(12,5),(12,6),(13,6),(14,6),
(14,7),(15,7),(16,7),(16,8)}.
Baselines: truncated full sorter of best-known size; SorterHunter median nets
(bottom-k with k=floor((N-1)/2)); greedy half-cleaner construction.
Success: size < best baseline for that (N,k).

## Campaign C — SAT

C1: symbolic counterexample checker (CNF over input vars + Tseitin gates);
validate against known-optimal s(n) for n<=8 by UNSAT on random non-sorting
nets (must be SAT) and on known sorting nets (must be UNSAT).
C2: CEGIS synthesis for selectors at small N; report any size < Campaign B
best. Lower-bound claims ONLY from solver UNSAT with completeness argument;
otherwise labeled conjecture/heuristic evidence.
