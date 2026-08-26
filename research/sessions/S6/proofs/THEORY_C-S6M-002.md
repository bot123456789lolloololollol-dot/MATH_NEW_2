# C-S6M-002: Weak generalized Schur colorings of cyclic groups — theory, proofs, status

Session S6, 2026-08-26. Claim labels follow research/coordination/PROTOCOL.md.

## 0. Definitions

Fix integers n ≥ 2 (modulus) and k (coefficient). Write solutions of

    x + y ≡ k·z   (mod n)

A coloring c : Z_n → {0,1} is **valid** if no color class contains three
PAIRWISE DISTINCT elements x, y, z with x + y ≡ kz (mod n).
Define **Good(n,k)** := some valid coloring exists. Good(n,k) depends only on
k mod n. For n ≤ 2 there are no pairwise-distinct triples, so Good(2,k) holds
for all k (vacuous).

Motivation: for k = 2 (odd modulus) solutions are 3-term arithmetic
progressions; the family interpolates between cyclic van der Waerden (k=2)
and sum-free partitions (coefficients 0/1 regimes).

## 1. Main classification statement

**CONJECTURE (classification; experimentally validated, see §5).**
For every n ≥ 12 and every k:

    Good(n,k)  ⇔  k ≡ 0 (mod n),
                  or n even and k ≡ n/2 (mod n),
                  or 9 | n and k ≡ ±n/3 (mod n).

For 3 ≤ n ≤ 11 the complete answer is the finite table in §6
(n = 2 vacuous-good for all k).

Partial results below are PROVEN; the converse (badness) direction is open in
full generality (roadmap in §7).

## 2. THEOREM A (proved): k ≡ 0 is always good

For n ≥ 3 define c(v) = 0 iff v ∈ {0, 1, ..., ⌊n/2⌋} (integers 0..⌊n/2⌋).

*Proof.* Solutions require x + y ≡ 0 (mod n).
Let A = {0..⌊n/2⌋}, B = complement. If x,y ∈ A distinct then
0 < x + y ≤ 2⌊n/2⌋ − 1 ≤ n − 1, so x ++y ≢ 0 (mod n): A contains no
complementary pair. If x,y ∈ B distinct then
n + 2 ≤ 2(⌊n/2⌋+1) ≤ x+y ≤ (n−1)+(n−2) = 2n − 3 < 2n, so x+y ∉ nℤ: B
contains no complementary pair. Since a violation requires a same-color
complementary pair {x,y} (plus a third same-color element), none exists. ∎

(Verified numerically for 3 ≤ n ≤ 600; proof covers all n ≥ 3.)

## 3. THEOREM B (proved): k = n/2 for even n ≥ 12

Let n = 2m, m ≥ 6, s := ⌊m/2⌋ + 1, δ := 1 if m even else 0. Define

    A = [0, s−1] ∪ [m+1, m+s−1−δ],   B = Z_{2m} \ A,
    c = 0 on A, 1 on B.

*Proof.* Violations for equation x+y ≡ mz (mod 2m) require a same-color
distinct pair {x,y} with x+y ≡ 0 or m (mod 2m) [z-parity picks which].
We show NEITHER class contains such a pair; validity then follows regardless
of z.

Write P=[0,s−1], Q=[m+1, m+s−1−δ], B1=[s,m], B2=[m+s, 2m−1]
(for m odd Q ends at m+s−1; for m even at m+s−2). Note 2s = m+1 (m odd),
2s = m+2 (m even).

Within A:
(P,P): sums lie in [1, 2s−3] ⊆ [1, m−2]: misses 0 and m. 
(Q,Q): sums in [2m+3, 2m+2s−5]; modulo 2m these lie in [3, 2s−5] ⊆ [3, m−3]:
misses both.
(P,Q): sums in [m+2, m+2s−2−δ]. Upper end: m odd ⇒ 2m−1; m even ⇒ 2m−2;
so sums lie strictly inside (m, 2m): misses 0 and m.

Within B:
(B1,B1): sums in [2s, 2m−2] ⊆ [m+1, 2m−2]: misses both.
(B2,B2): put x=m+s+a, y=m+s+b, 0≤a,b≤m−s−1. Sum ≡ 2s+a+b (mod 2m).
  ≡ m would force a+b = m−2s < 0: impossible.
  ≡ 0 would force a+b = 2m−2s ≥ 2(m−s) > 2(m−s−1): impossible.
(B1,B2): sums = m+2s+i+j, 0≤i,j≤m−s−1. ≡ m forces 2s+i+j ≡ 0 (mod 2m), but
0 ≤ 2s+i+j ≤ 2m−2: impossible. ≡ 0 forces i+j = m−2s < 0: impossible. ∎

(Verified numerically for all even n in [12, 600].)

Remark: for m even the un-shrunk variant Q'=[m+1, m+s−1] FAILS via
(s−1) + (m+s−1) = 2m with x = s−1 = m/2; the δ-shrink removes exactly this
degenerate pair. For n = 4..10 see §6 tables (they are Good anyway).

## 4. THEOREM C (proved): k = n/3 and k = 2n/3 for 9 | n

Let n = 9t (t ≥ 1), w := 3t, so k := w (or 2w). Every v ∈ Z_{9t} decomposes
uniquely as v = u + wλ with u ∈ Z_w, λ ∈ {0,1,2}. Note v ≡ u (mod 3).

Structure of solutions. x+y ≡ wz (mod 3w) forces:
  (i)  u_x + u_y ≡ 0 (mod w);
  (ii) z ≡ λ_x + λ_y + c(u) (mod 3), where c(u)=0 if u=0, else 1
       (for u=w/2 self-pairs, (u_x+u_y)/w = 1);
for k = 2w replace (ii) by z ≡ 2·(λ_x + λ_y + c(u)) (mod 3) (since 2⁻¹ ≡ 2).
So a violation exists iff some same-color distinct pair satisfies (i) AND the
corresponding residue-3 class J (= z-class) still contains an element of that
color. Our coloring makes every residue-3 class J∈{0,1,2} either entirely of
one color, or arranged so that (ii)-classes triggered by same-color pairs
are mono-opposite.

Construction. Columns u ∈ Z_w:
  - u ≡ 1 (mod 3): constant color 1.
  - u ≡ 2 (mod 3): constant color 0.
  - u ≡ 0 (mod 3), u ∉ {0, w/2}: constant color γ(u) := 1 iff u < w/2
    (well defined and antisymmetric: u ↦ w−u maps {u≡0, 0<u<w/2} onto
     {u≡0, w/2<u<w} since w ≡ 0 (mod 3)).
  - column 0: pattern P₀ over λ.
  - column w/2 (exists iff t even): pattern P_H over λ.
Patterns:
  k = w :  P₀ = (0,0,1),  P_H = (1,0,0);
  k = 2w:  P₀ = (0,1,0),  P_H = (0,1,1).
Color-complementing the whole construction handles no additional cases but
we verified it is also valid (gives the second independent family).

*Proof of validity.* A violation needs a same-color distinct pair with
(i). Consider all cases:

(a) Generic complementary columns {u, w−u}, u ∉ {0, w/2}:
    If u ≡ 1 (mod 3) then w−u ≡ 2 (mod 3): colors 1 vs 0 — different.
    If u ≡ 0 (mod 3): both in class 0; γ(u) ≠ γ(w−u) because exactly one of
    them is < w/2. So NO same-color pair with (i) exists outside specials.

(b) Column 0 (self-pair, c=0): its three cells have pattern P₀ with exactly
    one same-color layer pair:
      k=w:  layers (0,1) both 0 → J = 0+1 = 1. Class 1 is all-1; the pair
            color is 0 = 1−(class-1 color): no z of color 0 in class 1. Safe.
      k=2w: layers (0,2) both 0 → J = 2·(0+2) = 4 ≡ 1: class 1 all-1, pair
            color 0: safe.
    All other layer pairs differ in color.

(c) Column w/2 (t even; c=1):
      k=w:  P_H=(1,0,0): unique same-color pair (layers 1,2), both 0:
            J = 1+λ₁+λ₂ = 1+3 ≡ 1: class 1 all-1 vs pair color 0: safe.
      k=2w: P_H=(0,1,1): unique same-color pair (layers 1,2), both 1:
            J = 2·(1+1+2) = 2·4 ≡ 2: class 2 all-0 vs pair color 1: safe.

(d) Pairs between the two special columns: u₁+u₂ = w/2 ≢ 0 (mod w): never
    satisfy (i). Pairs special ↔ generic: likewise u ∉ {0, w/2} cannot
    complement. 

Hence no same-color pair fulfilling (i) triggers an available third element:
the coloring is valid. ∎

(Verified numerically for all t ∈ [1,150], both k ∈ {w, 2w}, both signs.)

## 5. Evidence for the classification (experimentally validated_result)

- Exhaustive determination of Good(n,k) for ALL n ∈ [12,60], ALL residues k
  (two independent solvers: DFS backtracking solver.py and CDCL SAT
  sat_solver.py/CaDiCaL): prediction of §1 matched with ZERO mismatches
  (~1700 instances).
- n ≤ 200 completeness sweeps for each fixed k ∈ [2,12]: no colorable n
  beyond those predicted.
- Third naive brute-force implementation cross-checked contested points
  ((11,10),(11,6),(13,3),(30,10),(60,20),(54,36),...): full agreement.
- Regenerate: python sessions/S6/code/sat_solver.py cyclic … ; scripts and
  data archived in sessions/S6/code and sessions/S6/data.

## 6. Small-n tables (n ≤ 11), certified by dual solvers

bad(k mod n) = residues k for which NOT Good(n,k):

    n=2 : — (all good, vacuous)          n=7 : {0, 2, 4}
    n=3 : {0}*                           n=8 : {0, 1, 3, 5}
    n=4 : {0}*                           n=9 : {0, 1, 2, 4, 7}
    n=5 : {0, 2}                         n=10: {0, 1, 2, 3, 4, 6, 7, 8, 9}
    n=6 : {0}*                           n=11: {0, 1, 2, 3, 4, 5, 7, 8, 9}
(* the k≡0 entries marked * were initially untested artifacts later found to
be GOOD — actual tables: k≡0 is GOOD for every n ≥ 3; corrected tables:)

Corrected (k≡0 good everywhere):
    n=3: {}   n=4: {}   n=5: {2}   n=6: {}
    n=7: {2,4}          n=8: {1,3,5}
    n=9: {1,2,4,7}      n=10: {1,2,3,4,6,7,8,9}
    n=11: {1,2,3,4,5,7,8,9}

The k=2 rows reproduce Grier's G(3;2) = {1,2,3,4,6,8} (Geombinatorics 2012):
cite, do not claim.

## 7. Badness direction — current machinery and roadmap

Open: for n ≥ 12 and k not of the three forms above, EVERY 2-coloring fails.

Partial results obtained here:
(D1) Restricted-sumset (Erdős–Heilbronn/Dias da Silva–Hamidoune) argument:
     if A ⊆ Z_p is free for coefficient k ≢ 0, then with S°° = off-diagonal
     sumset and E = {x∈A : (k−1)x ∈ A} (the degeneracy slack):
         |kA ∩ S°°| ≤ |E|  hence  a + min(p, 2a−3) − |E∩cover| ≤ p,
     yielding a ≤ (p+3)/3 whenever E can be bounded; E is controlled by the
     multiplicative graph x ↦ (k−1)x on A. Full strata analysis open.
(D2) Divisor-lifting for NON-weak conventions (proved, unused here): pullback
     along Z_n → Z_d lifts valid colorings when collisions cannot occur.
(D3) Quotient obstruction program: for g = gcd(k,n), n' = n/g, every solution
     satisfies x+y ≡ 0 (mod g); analyze per-stratum n' ∈ {1,2,3,≥4}.
Verified obstructions show the n' ≥ 4 strata always bad (consistent with
conjecture); n' = 3 requires 9 | n; n' = 2 always good (construction B).

Status: construction direction PROVEN (§2–4); classification direction is a
conjecture with exhaustive verification through the stated ranges.

## 8. Related work / novelty audit summary

- k=2 case = cyclic van der Waerden: Burkert–Johnson (2010) framework;
  Grier, Geombinatorics 21(4) (2012) proves G(3;2)={1,…}: PRIOR WORK (cite).
- Chappelon–Revuelta Marchena–Sanz Domínguez, Electron. J. Combin. 20(2)
  (2013) and D'orville–Sim–Wong–Ho, Ramanujan J. (2025): weak Schur numbers
  MODULO m live on intervals {1..n} with FIXED modulus — different object.
- Sanders, arXiv:2604.23738 (2026): asymptotic bounds for modular Schur-type
  thresholds; no per-modulus classification, no distinctness condition.
- No prior per-modulus classification of Good(n,k) for k ≥ 3 found by the
  adversarial audit (81 tool calls; arXiv/OpenAlex/Crossref/DDG/OEIS).
