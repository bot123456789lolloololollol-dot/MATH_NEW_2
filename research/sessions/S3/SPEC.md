# S3 SPEC — Adversarial sweep of four tree-labeling conjectures at reachable sizes

Session: S3 (retry). Date: 2026-08-26.
Domain: adversarial mathematics — hunting counterexamples to published open
conjectures by exhaustive computation at reachable sizes.

## Triage (status-checked BEFORE compute, see literature/S3-novelty.md)

Web tooling available this session: WebFetch only (no WebSearch). Sources used:
Wikipedia "Graceful labeling", arXiv API queries (2026-08-26).

| # | Conjecture | Published status found | Prior exhaustive verification found? | Decision |
|---|---|---|---|---|
| C1 | Odd-graceful: every tree is odd-graceful (Gnanajothi 1991) | OPEN; proven for caterpillar forests and diameter <= 5 (Barrientos 2008, arXiv:0807.3738) | none found on arXiv | RUN |
| C2 | Prime: every tree has a prime labeling (Entringer ~1980) | OPEN; still listed as conjecture in Rao survey arXiv:2006.03801 | none found on arXiv | RUN |
| C3 | Antimagic: every connected graph != K2 is antimagic, in particular every tree (Hartsfield-Ringel 1990) | OPEN "even for trees" (arXiv:1905.06595, 2506.15221); Meng 2026 machine-checked only odd-order trees with exactly two degree-2 vertices up to 25 | none covering ALL trees up to a bound | RUN |
| C4 | Neighborhood-prime: every tree of order >= 3 admits a neighborhood-prime labeling (Ryan 2014; Cloys-Fox arXiv:1801.01802 prove caterpillars/spiders/firecrackers) | OPEN in general per those abstracts | none found | RUN |
| — | Graceful tree conjecture | OPEN but verified by computer to 27 (Aldred-McKay 1998), 29 (Horton), claimed 35 (Fang project) | YES, far beyond our budget | SKIP |
| — | Harmonious trees (Graham-Sloane 1980) | OPEN but Fang 2011 (arXiv) verified ALL trees <= 31 nodes | YES, beyond our budget | SKIP |

## Finite consequences tested (exact formalization)

For each labeling L and bound N: **every unlabeled free tree T with |V(T)| = n,
1 <= n <= N, admits an L-labeling**; verdict VERIFIED_UP_TO_N or REFUTED.

- odd-graceful(T): injective f: V -> {0,...,2m-1} (m = |E| = n-1) with edge
  differences {|f(u)-f(v)|} = {1,3,...,2m-1}.
- prime(T): bijective f: V -> {1,...,n} with gcd(f(u),f(v)) = 1 on every edge.
- antimagic(T): bijection g: E -> {1,...,m} such that vertex sums
  s(v) = sum of g(e) over edges e incident to v are pairwise distinct.
- neighborhood-prime(T): bijective f: V -> {1,...,n} with
  gcd{f(u) : u in N(v)} = 1 for every v with deg(v) >= 2. (n >= 3.)

## Method

1. Complete enumeration of unlabeled free trees via `lib_s3.gen_free_trees`
   (centroid construction); counts validated against OEIS A000055 in
   `lib_s3.self_test` (passed previously for n <= 16; re-run each session).
   Completeness argument: every free tree has 1 or 2 centroids; rooting at a
   centroid gives either (a) branches all of size <= (n-1)/2 around one root,
   or (b) two half-size rooted halves joined at their roots; the generators
   enumerate exactly these multisets/pairs, so no tree is missed and none
   duplicated (seen-set dedup).
2. Exact backtracking solvers per labeling (bitmask candidate domains,
   MRV ordering, forward checking; details in run_labelings.py docstrings).
   Any solver output that fails the INDEPENDENT checker function aborts the
   run with a logged solver-bug record (checker is written from the definition,
   not from the solver).
3. Escalation tier for odd-graceful: Z3 encoding (z3-solver installed OK).
4. Runs are chunked: one process per mode with internal wall-clock deadline
   (< 300 s) and incremental JSONL logging after each completed order n, so a
   killed process loses nothing. A level n counts toward a verdict only if
   recorded `"complete": true`.
5. Checker corruption self-test: each checker must REJECT mutated valid
   labelings (label swaps, value perturbations, duplicate labels) and ACCEPT
    the solver's output; test runs before any compute.

## Claim discipline

Expected verdicts are `experimentally_validated_result` of the form
"VERIFIED_UP_TO_N" (complete enumeration + exact solvers + independent
checkers + completeness argument above), NOT proofs of the conjectures.
A single genuine UNSAT from the exact solver would be a counterexample and
would require independent reimplementation before any REFUTED claim.
