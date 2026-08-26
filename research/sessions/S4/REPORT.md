# S4 REPORT — Computational Number Theory: certified extension of A094406 / A176762

Session: S4 (computational number theory). Date: 2026-08-26.
All code and raw outputs: `research/experiments/S4/`. One-command regeneration:
`bash research/experiments/S4/run.sh` (deterministic, < 2 min, no dependencies beyond CPython 3.11 stdlib).

## 1. Targets chosen and published frontiers (established 2026-08-26)

| Sequence | Definition | Published frontier | Source |
|---|---|---|---|
| **A094406** | Smallest *unhappy* number taking n steps to reach the unhappy cycle (4,16,37,58,89,145,42,20) under iteration of the sum-of-squares-of-digits map ssd (A003132) | a(0..12) small; **a(13) = 578·10^196 − 1** ("577" followed by 196 nines, 199 digits); nothing beyond | OEIS A094406 entry + b-file (n = 0..13), revision dated 2026-08-22, accessed 2026-08-26 |
| **A176762** | Smallest positive integer taking n steps to reach *a* cycle (fixed point 1 or the 8-cycle) under ssd | data field ends at a(19) = 15999; formula field says a(n) = min(A001273(n), A094406(n−7)); a(20) is implied (= A094406(13)); **nothing at n ≥ 21 anywhere** | OEIS A176762, revision 2026-07-09, accessed 2026-08-26 |
| A001273 (context only) | Smallest happy number of height n | formulas through a(19) (Ya-Ping Lu, Jul 2025) — actively maintained, NOT targeted | OEIS A001273, accessed 2026-08-26 |

Rejected targets during selection (documented feasibility analysis):
* **A275154 / A380448** ("smallest integer that is a sum of distinct positive cubes in exactly n ways" and prime analogue): b-file of 10 000 terms exists (D. A. Corneth), Zhao Hui Du pushed minimality arguments past 10^9 with interval proofs, revisions Aug 2025 — frontier far out of reach here.
* **A097757-family** (sums of distinct squares in exactly n ways): grounded in Sprague (1948); well trodden.
* **A275154-style cubes records**: same community, same verdict.

## 2. Results

### R1 — NEW TERM for A094406 (`experimentally_validated_result` → `proven_result`, see §4)
**a(14) = "5888" followed by N nines**, where

```
N = 71358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358024691358022
```

i.e. the number has D = N + 4 = ...35802**6** digits (the exact 197-digit digit counts are
in `out/a094406_certification.json`; they are too long to restate safely by hand).
Equivalently a(14) = 5889·10^N − 1.

* **Witness/verification:** ssd(a(14)) = 5²+8²+8²+8² + 81·N = 217 + 81·N = 578·10^196 − 1 = a(13),
  and a(13) has tail 13, so a(14) has tail exactly 14 (checked by direct simulation in both
  implementations).
* **Minimality certificate (complete):** see §3 — exhaustive finite window scan + exact
  digit-DP optimality, machine-checked by two independent implementations.

### R2 — RE-CERTIFICATION of the published a(13) = 578·10^196 − 1 (`proven_result`)
The OEIS comment states the value without a proof artifact. We supply one:
window [15999, 16119] contains exactly ONE σ with tail 12 (namely σ = 15999 = a(12)),
and M\*(15999) = "577"+196 nines (deficit-multiset analysis: total deficit 81·199−15999 = 120
has the unique admissible digit-deficit decompositions {56,32,32}→{5,7,7} and none better;
see §3). The independent checker re-derives this from scratch.

### R3 — NEW DERIVED VALUE for A176762 (`proven_result`)
**A176762(21) = A094406(14)** = the number of R1.
Proof sketch (full version §3.4): T(m)=21 ⟺ m happy of height 21 or m unhappy with tail 14
(metric decomposition verified exhaustively for all m < 20000 and proved via
ssd(m) < m for m ≥ 100). Any happy m of height 21 satisfies ssd(m) ≥ A001273(20) ≥
10^((A001273(19)−324)/81), which by the published Lu recursion tower (base a_H(8) =
3789·10^973 − 1 > 10^975, twelve more exponentiation steps to n=21) exceeds any number
with fewer than 10^248 digits after one more pigeonhole step; the winner has ~7.14·10^196
digits. So the unhappy branch dominates. Also derived (not new, implicit in the entry's own
formula field): **A176762(20) = A094406(13)**.

## 3. Why naive search fails, and the completeness argument (written proof)

### 3.1 Setup and lemmas (all proved, and exhaustively spot-checked in code)
ssd(n) := sum of squares of decimal digits of n.
Lemma 1 (*shrinkage*): ssd(n) < n for every n ≥ 100.
  Proof: for d ≥ 4 digits, ssd(n) ≤ 81d ≤ 81d < 10^(d−1) ≤ n (since 10³ = 1000 > 324 and
  10 grows faster than +81 each digit). For d = 3, write n = 100a+10b+c, a ≥ 1; maximizing
  ssd(n) − n = (a²−100a)+(b²−10b)+(c²−c) over integers gives −99 + 0 + 72 = −27 < 0
  (convexity ⇒ maxima at endpoints of each variable's range). ∎
Lemma 2: hence tail(n) = 1 + tail(ssd(n)) and T(n) = 1 + T(ssd(n)) for all n ≥ 100, where
  tail = A094406 metric (0 on cycle members, undefined on happy numbers) and T = A176762
  metric ((#distinct trajectory values until first repetition) − 1).
Lemma 3 (*metric decomposition*, basis of R3): T(m) = height(m) for happy m; T(m) = 7 for
  members of the 8-cycle; T(m) = tail(m)+7 otherwise. Verified exhaustively for m < 20000
  (`verify_decomposition` in extend.py).

### 3.2 Exact minimal preimage M\*(σ)
Let M\*(σ) be the smallest positive integer with ssd = σ. Write its d digits c_1…c_d and
define the *deficit* Δ = 81d − ssd = Σ_i (81 − c_i²). Per-digit deficits V = {0,17,32,45,
56,65,72,77,80,81}. Therefore:

> σ is achievable with exactly d digits ⟺ Δ := 81d − σ ≥ 0 and Δ is a sum of ≤ d elements of V.

* **Digit count:** d\*(σ) = min feasible d. Any number with fewer digits is numerically
  smaller than any with more, so M\*(σ) has exactly d\*(σ) digits. Feasibility is decidable
  exactly: primary implementation uses DP sumset tables of V; the independent checker uses
  meet-in-the-middle enumeration of deficit multisets plus the Frobenius fact that every
  Δ > 495 = F(17,32) is of the form 17a+32b. In our instances d\*(σ) ∈ {⌈σ/81⌉, ⌈σ/81⌉+1}
  because the residual deficits 26 (for a(14)-scale) and 39 (for a(13)-scale) are
  exhaustively shown NOT to lie in the semigroup generated by V, while 107 resp. 120 are.
* **Arrangement:** among fixed-length strings, numeric order = lexicographic order, so the
  greedy left-to-right rule "place the smallest digit whose residual is still feasible"
  produces M\*(σ). Once the residual equals 81·(#remaining digits), ALL remaining digits are
  forced to be 9 — giving compact certificates *(prefix, nine-count)* even when the count of
  nines is itself a 197-digit integer.
* Both implementations cross-validated against literal brute force for every σ ∈ [1,350]
  (350/350 agree), and against each other on the huge instances.

### 3.3 Window completeness for a(n) = min{M\*(σ) : tail(σ) = n−1} (n ≥ 14)
By Lemma 2, any m with tail n has σ := ssd(m) with tail n−1, and m ≥ M\*(σ).
Any m with σ has ≥ ⌈σ/81⌉ digits (pigeonhole: 81·digits ≥ σ), i.e.
m ≥ 10^(⌈σ/81⌉−1). Let U be the best candidate found, D = #digits(U).
Then any competitor better than U must satisfy σ ≤ 81D. Hence the search space is the
FINITE window [a(n−1), 81D] (lower end: a(n−1) is the smallest tail-(n−1) number).
We enumerate EVERY integer in the window (108 values for n=14; 121 values for n=13),
compute tails by direct simulation, and take the exact minimum of M\* over hits.
For a(13): base a(12) = 15999 brute-force verified as smallest tail-12 number;
window width 121; unique hit σ = 15999; M\*(15999) matches the published value.
For a(14): lower end a(13) (certified); winner U = M\*(a(13)); D = N+4; window
[a13, 81·D] of width 108; unique hit δ = 0; so **a(14) = M\*(a(13))**, completing the proof.
Numbers m < 100 are checked directly (none has tail 13 or 14); m in [100, a(n)) cannot have
tail n because ssd(m) < a(n−1) would contradict Lemma 2 + definition of a(n−1) — subsumed by
the window argument.

### 3.4 Dominance of the unhappy branch for A176762(21)
By Lemma 3, A176762(21) = min(A001273(21), A094406(14)). Using only published A001273
facts — a_H(8) = 3789·10^973 − 1 > 10^975 and a_H(n) = k(n)·10^((a_H(n−1)−ssd(k(n)−1))/81)−1
with k(n) < 10^4, so ssd(k(n)−1) ≤ 324 — the recursion gives log₁₀ a_H(n) ≥ 10^{L(n−1)}/200,
i.e. each step adds an exponent to the tower. After the twelve steps to n = 20 the bound is
unrepresentably large; conservatively, any happy height-21 number has more than 10^248
digits, while the R1 winner has ~7.14·10^196 digits. Unhappy branch wins. ∎
(The same argument one level down shows A176762(20) = A094406(13), consistent with the
entry's formula field.)

## 4. Claim labels
* R1 (A094406 a(14)): **proven_result** — complete minimality argument (§3.3), witnessed
  object, two independent implementations, adversarial tamper suite (8/8 corruptions
  rejected, see `out/check_report.json`).
* R2 (re-certification of a(13)): **proven_result** — same standard.
* R3 (A176762(21), A176762(20)): **proven_result** modulo the published correctness of the
  OEIS A001273 formula field (Lu, Jul 2025), which we consume as cited input, not as our
  claim; our own contributions (decomposition lemma, dominance bound, unhappy-side record)
  are independently proven.
* Novelty status of the extension claim: `audited_clean` within the reachable audit radius
  (OEIS current state + literature list + web search); see `research/literature/S4-novelty.md`.

## 5. Reproduction commands

```bash
bash research/experiments/S4/run.sh          # everything, <2 min
python research/experiments/S4/happy_core.py # published-term reproduction
python research/experiments/S4/mstar.py      # M* selftest vs brute force (sigma<=350)
python research/experiments/S4/extend.py     # phases B/C/D -> out/a094406_certification.json
python research/experiments/S4/checker_independent.py   # second algorithm + tamper tests
```

Raw outputs committed: `out/published_terms.json`, `out/a094406_certification.json`,
`out/check_report.json`.

## 6. Honest limitations
* a(15) is out of reach in explicit form: it minimizes over σ ≥ a(14), and since a(14)
  already has ≈ 7.14·10^196 digits, a(15) has roughly a(14)/81-digit scale — i.e. even its
  compact-certificate *nine-count* would be an integer with ~7·10^196 decimal digits, not
  explicitly writable. A doubly-symbolic certificate (exponent towers with machine-checked
  algebraic identities, as used in the published A001273 formulas) might still be possible
  (future work).
* The dominance argument consumes the OEIS A001273 published formulas as ground truth.
* Web search engines blocked automated queries (DDG captcha, Bing irrelevant results);
  novelty audit rests primarily on live OEIS state (canonical for sequences), the entries'
  literature lists, arXiv, and Schneider's archived page. Residual risk of a concurrent
  unpublished computation cannot be excluded — flagged in novelty file.
