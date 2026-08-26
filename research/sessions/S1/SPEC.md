# SPEC — Certified exact maximum Sidon set sizes f(n) in Z_n, n = 2..110

Artifact: `research/experiments/S1/` (solvers, independent checker, raw outputs).
Purpose: a stranger must be able to reimplement the verifier from this document alone.

## 1. Problem

For integer n ≥ 2, define

    f(n) = max{ |A| : A ⊆ Z_n, all sums a+b (a,b ∈ A) pairwise distinct in Z_n }
           (a,b unordered with doubles included; i.e., i ≤ j)

This is the STRICT Sidon / modular-Golomb-ruler convention = OEIS **A260999**
("maximal size of a subset of Z_n with distinct sums of any two elements").
Its inverse a(k) = min{n : f(n) ≥ k} is OEIS **A004136**.

Deliverable: exact certified values of f(n) for every n in [2, 110], together
with a witness set per n, produced by deterministic exhaustive search.

## 2. Equivalent formulations (both machine-checked, must agree)

- **F1 (definitional):** the multiset {a+a} ∪ {a+b : a<b} has k(k+1)/2 values,
  all distinct mod n.
- **F2 (search form):** the k(k−1)/2 *circular distances*
  dd(a,b) = min((b−a) mod n, (a−b) mod n) are pairwise distinct AND, when n is
  even, no pair has dd = n/2.

**Lemma 1 (F1 ⟺ F2).** If two sums collide nontrivially, a+b ≡ c+d, then
a−c ≡ d−b gives equal ordered differences from distinct ordered pairs unless
{a,b}={c,d}; conversely a−b ≡ c−d implies a+d ≡ b+c. For doubles: a+a ≡ b+b
forces a ≡ ±b, and a ≡ −b (a≠b) happens iff n even and dd(a,b)=n/2. Hence F1
fails exactly when some circular distance repeats or some pair sits at n/2.
∎ (Both directions used by the checker: F1 and F2 are verified independently;
any disagreement is reported as ERROR, not FAIL.)

**Counting upper bound.** Circular distances take values in {1,…,⌊(n−1)/2⌋}
(the value n/2 at even n is unusable by Lemma 1). There are ⌊(n−1)/2⌋ usable
slots and C(k,2) pairs, so

    ub(n) = max{k : k(k−1)/2 ≤ ⌊(n−1)/2⌋},   f(n) ≤ ub(n).

This bound needs no search and is a proven upper bound.

## 3. Search (solver.js primary; solver.py independent re-expression)

Decision procedure `solve(n,k)`: does a k-element strict-Sidon set exist?

Space enumerated: sets {0} ∪ S, S ⊆ {1,…,n−1}, elements chosen in INCREASING
order; depth-first; candidates x tested in ascending order.

**Symmetry reductions and completeness arguments (why no solution class is lost):**

1. **Translation.** If A ⊆ Z_n is Sidon and t ∈ A, then A−t is Sidon and
   contains 0 (translation permutes Z_n, preserving sums-distinctness). Every
   solution class therefore contains a representative containing 0; we search
   exactly those. *Complete.*
2. **Ordered enumeration.** With 0 ∈ A fixed as minimum, every subset of
   {1,…,n−1} corresponds to exactly one branch of the DFS; increasing-order
   selection enumerates each subset once. *Complete.*
   (No multiplication-by-unit canonicalization and no reflection reduction are
   applied — these would be sound but are deliberately omitted so the
   completeness argument stays trivial. Cost: constant-factor slowdown only.)

**Insertion validity test (necessary and sufficient incrementally):** adding x
to current set T is allowed iff for every a ∈ T, dd(x,a) ∉ UsedSlots and,
for even n, dd(x,a) ≠ n/2. By Lemma 1 this enforces exactly the defining
condition on the growing set.

**Pruning rules (all preserve completeness — each only removes provably dead
branches):**

- P1 counting prune: at depth t with maximum element m, if (n−1) − m < k − t,
  too few residues remain to place k−t more marks. *Sound.*
- P2 static bound: levels k > ub(n) are never searched (Section 2). *Sound.*

No dominance, memoization, or learned clauses are used. The search explores
the FULL tree up to P1/P2 cuts, so "UNSAT" outcomes are exhaustive proofs.

**Protocol per n:** for k = ub(n), ub(n)−1, …: run solve(n,k); stop at first
SAT. f(n) = that k. Certification requires: every level ABOVE the answer
returned UNSAT (no timeout), and the answer level returned SAT. Any timeout ⇒
n marked uncertified and excluded from claims (none occurred for n ≤ 110).

**Determinism:** no randomness, no maps/dicts in the hot path (arrays only),
fixed iteration order; reruns reproduce identical node counts and witnesses
(verified: n=92 → 38,251,495 nodes twice; n=105 → 244,412,475 nodes twice).

## 4. Independent checker (check_sidon.py)

Given (n, set S, expected_k), verifies from scratch:
residues distinct mod n; |S| = expected_k; F1 via O(k²) sum-multiset check;
F2 via O(k²) distance check; both must accept, disagreement ⇒ ERROR.
The checker shares NO code with either solver and uses the definitional
formulation, guarding against solver-side misconceptions.

**Corruption battery (`--selftest`, must always pass):** genuine objects
{(7,{0,1,3}), (13,{0,1,3,9}), (16,{0,1,4,6})} PASS; deliberately broken
objects each FAIL: duplicate residue; sum collision (5,{0,1,2}); half-modulus
pair (8,{0,4}); repeated circular distance (12,{0,5,10}); wrong claimed k;
and a WEAK-Sidon object valid under the A260998 convention but not ours
(6,{0,1,2,4}, where 2+2 = 0+4) — proving the checker distinguishes the strict
convention from the weak one.

## 5. Third implementation (ground truth at small n)

`run_tests.py::naive_kmax` brute-forces ALL subsets via itertools (no symmetry
assumptions beyond heredity of Sidon-ness) for n ≤ 14 and agrees with both
solvers; python vs node agree on values and witnesses for n = 2..30.

## 6. Verification against published data (outputs/verification_upto_110.json)

For all 109 certified n ∈ [2,110]:
- witness passes independent checker (109/109);
- f(n) equals OEIS A260999 b-file value (0 mismatches);
- Cariboni's published witness set (third_party/w260999.txt) passes the strict
  checker with size equal to the b-file value (109/109 audited, 0 failures);
- inverse terms recomputed a(2..10) = 3,7,13,21,31,48,57,73,91 match OEIS
  A004136 exactly (k=1 differs only because OEIS counts the trivial group Z_1,
  outside our n ≥ 2 domain).

## 7. Reproduction

```
bash research/experiments/S1/run.sh          # everything, ~20–25 min total
node research/experiments/S1/solver.js --only 61            # single case
python research/experiments/S1/check_sidon.py --selftest     # checker battery
```
Raw outputs: `outputs/main_chunk{A,B,C}.json`, `outputs/tests_raw.json`,
`outputs/verification_upto_110.json`, `outputs/summary_n2_110.csv`,
`outputs/agreement_js.json`, `outputs/determinism_spot.json`.

## 8. Known limits

- n ≥ 111 with f(n) < ub(n) (e.g., 111..119 where ub=11 but f=10) need UNSAT
  proofs at k=11, projected hours-days in this implementation; NOT attempted.
- f is NOT monotone in n (e.g., f(21)=5 > f(22)=4): each n is solved
  independently; monotonicity is never assumed.
