# SPEC — S4 artifact: certified extension of OEIS A094406 / A176762

Purpose: a stranger must be able to REIMPLEMENT verification from this spec alone and
confirm the two new certified values. Everything here is self-contained; no trust in our
scripts is required.

## 0. Definitions
* `ssd(n)` = sum of squares of the decimal digits of n (OEIS A003132).
* Unhappy cycle C = {4, 16, 37, 58, 89, 145, 42, 20}; 1 is a fixed point.
* `tail(n)`: iterate x ← ssd(x); if n's orbit reaches 1, tail is undefined (happy);
  otherwise tail(n) = number of steps until an element of C is first reached.
  Members of C have tail 0.  [A094406 metric]
* `T(n)`: number of distinct values of the trajectory n → ssd(n) → … up to (excluding)
  the first repeated value, minus 1.  [A176762 metric]

## 1. Facts you must verify yourself (all cheap)
F1. ssd(n) < n for all n ≥ 100 (finite case check d=3; inequality 81d < 10^(d−1) for d ≥ 4).
    Consequently tail(n) = 1 + tail(ssd(n)) and T(n) = 1 + T(ssd(n)) for n ≥ 100.
F2. Brute force over m = 1..30000 reproduces the published terms:
    A094406 a(0..12) = 4, 2, 11, 15, 5, 3, 14, 45, 36, 6, 112, 269, 15999
    A176762 a(0..19) = 1, 10, 13, 23, 19, 7, 356, 4, 2, 11, 15, 5, 3, 14, 45, 36, 6,
                       112, 269, 15999
    ("smallest m whose metric equals k", first occurrence per metric value).
F3. Metric decomposition: T(m) = height(m) for happy m; T(m) = 7 on C; T(m) = tail(m)+7
    otherwise. Verify exhaustively for all m < 20000.

## 2. The minimal-preimage function M*(σ)
M*(σ) := smallest positive integer with ssd = σ.

Per-digit deficits v(c) = 81 − c² for digits c = 0..9 give V = {0,17,32,45,56,65,72,77,80,81}.
KEY CHARACTERIZATION: σ is the ssd of some string of exactly d digits ⟺
Δ := 81·d − σ ≥ 0 AND Δ is a sum of at most d elements of V.
Decide this exactly by enumerating multisets of nonzero deficits (Δ ≤ ~110 in all uses
here — trivially enumerable; for Δ > 495 use that every Δ = 17a+32b beyond the Frobenius
number F(17,32) = 17·32−17−32).

Then:
* d*(σ) = min feasible d.  (Any shorter length beats any longer length numerically.)
* M*(σ)'s digit string is built greedily left→right: place the smallest digit (≥1 in the
  first position) such that the residual σ′ is feasible with the remaining positions.
* EARLY EXIT: once residual = 81·(remaining positions), every remaining digit is forced
  to be 9.  Hence M*(σ) is always of the form prefix + "9"*N with a SHORT prefix.
Sanity anchor: M*(15999) = "577" followed by 196 nines = 578·10^196 − 1 (199 digits),
because Δ(199) = 81·199 − 15999 = 120 decomposes over V\{0} ONLY as {56,32,32}
(digits {5,7,7}); Δ(198) = 39 has NO decomposition; arrangements of the multiset are
lexicographically minimized by ascending order "577".
Cross-validate your M* against literal brute force for every σ ∈ [1,350] before trusting it.

## 3. Window completeness theorem (the minimality engine)
For n ≥ 14:  a(n) = min{ M*(σ) : tail(σ) = n−1 }   (∞-safe: only unhappy σ count).
Proof obligations:
(i)  any m with tail(m)=n has σ=ssd(m) with tail(σ)=n−1 (F1) and m ≥ M*(σ);
(ii) any m with ssd(m)=σ needs ≥ ⌈σ/81⌉ digits ⇒ m ≥ 10^(⌈σ/81⌉−1);
(iii) therefore any competitor better than the best candidate U (#digits D) has
     σ ≤ 81·D.  Combined with σ ≥ a(n−1) (definition of a(n−1)), scan the FINITE window
     [a(n−1), 81·D] exhaustively, compute tail for each integer (one ssd step drops it to
     ≤ 81·199 ≈ 16100, then simulate), and minimize exact M* over hits.
Also check no m < 100 has the target tail (direct simulation of 1..99).

## 4. The two certified values
### 4.1 A094406 a(13) (published value, re-certified)
a(12) = 15999 (verified in F2). Candidate U = "577"+196 nines, D = 199, window [15999, 16119].
Exhaustive scan: exactly one σ in the window has tail 12 — σ = 15999 itself.
So a(13) = M*(15999) = 578·10^196 − 1.  ✔ matches the OEIS comment (no prior proof existed).

### 4.2 A094406 a(14) — NEW
σ₀ := a(13) = 578·10^196 − 1; verify tail(σ₀) = 13 directly (ssd(σ₀) = 15999 → 12 more steps
to C). U = M*(σ₀): Δ₀ = 81⌈σ₀/81⌉ − σ₀ = 26, which has NO decomposition over V\{0};
Δ₀ + 81 = 107 decomposes ONLY as {56,17,17,17} → digits {5,8,8,8} or {45,45,17} → {6,6,8};
"5888" < "668" lexicographically ⇒ U = "5888" + "9"*N with N := (σ₀ − 217)/81 + ... precisely
N = d* − 4 where d* = ⌈σ₀/81⌉ + 1, i.e.

N  =  71358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358022

(check: 217 + 81·N = σ₀ exactly; total digit count D = N+4 ends in …026.)
Window [σ₀, 81·D], width 108: exhaustive rescan finds EXACTLY ONE σ with tail 13, namely
δ = 0 (σ = σ₀).  Hence **a(14) = U = 5889·10^N − 1**, and the certificate is complete.

### 4.3 A176762(21) — NEW derived value
By F3, A176762(21) = min(A001273(21), A094406(14)).  Using only published A001273 data
(a_H(8) = 3789·10^973 − 1 > 10^975; recursion a_H(n) = k(n)·10^((a_H(n−1)−ssd(k(n)−1))/81) − 1
with k(n) < 10⁴ so ssd(k(n)−1) ≤ 324), twelve recursion steps make log₁₀ a_H(20) tower past
all representation; pigeonhole (ii) then forces EVERY happy height-21 number to have more
than 10^248 digits, while a(14) has ≈ 7.14·10^196 digits.  So A176762(21) = a(14).
(Same argument one level down gives A176762(20) = a(13).)

## 5. What to check programmatically (recommended test list)
1. F2 brute-force term reproduction; F3 decomposition sweep.
2. M* vs brute force for σ ≤ 350.
3. tail anchors: tail(4)=0, tail(2)=1, tail(15999)=12, tail(578·10^196 − 1)=13,
   tail(M*(15999)) = 13, ssd chains as above.
4. Feasibility characterisation: Δ(198,15999)=39 infeasible, Δ(199,15999)=120 ↔ {56,32,32};
   Δ(d*−1, σ₀)=26 infeasible, Δ(d*,σ₀)=107 ↔ {56,17,17,17}|{45,45,17}.
5. Lex-minimality: for each prefix position i and each digit c < prefix[i] (c ≥ 1 at i=0),
   assert residual infeasible with remaining positions.
6. Full window rescans of §4.1/§4.2 with independently written tail code.
7. Tamper tests: corrupting the prefix (digit flip OR permutation), ±1 nine-count, fake
   smaller hit, wrong metadata tail, shifted window bound must ALL be rejected.

Our implementations: primary = `mstar.py` + `extend.py` (DP sumset tables);
independent = `checker_independent.py` (meet-in-the-middle deficit enumeration, fresh
simulation, re-derived bounds, tamper suite). Raw certificates: `out/a094406_certification.json`.
