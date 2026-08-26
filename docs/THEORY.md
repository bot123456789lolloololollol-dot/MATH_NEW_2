# Theory notes and provable properties

## 1. Offline line: structural dominance of the monotone-repair wrapper

**Setting.** One-dimensional bin packing, capacity c, item sizes integers in
[1, c]. Let FFD denote First-Fit-Decreasing and BFD Best-Fit-Decreasing.
Dósa (COCOA 2007) proved the tight bound

    FFD(I) <= (11/9) OPT(I) + 6/9        for every instance I,

and the bound is attained for every k by instances with OPT = 6k built from
blocks {1/2+e, 1/4+e, 1/4-2e} and {1/4+2e, 1/4+2e, 1/4-2e, 1/4-2e}.

**Algorithm A(r_max, patience).** Start from any constructive solution P
(e.g. min over {FFD, BFD}). Repeat: choose the r bins of smallest load,
empty them, re-insert their items into the residual capacity of the kept
bins in descending order (first-fit), pack the unplaced leftovers among
themselves by BFD into fresh bins; accept the move iff the total number of
open bins strictly decreases. Stop after `patience` consecutive rejections or
`max_rounds` moves.

**Theorem 1 (monotonicity / dominance).** Algorithm A terminates, and its
output uses at most as many bins as its starting solution. Consequently,
started from FFD it satisfies A(I) <= FFD(I) <= (11/9)OPT(I) + 6/9, and
started from min{FFD,BFD} it satisfies A(I) <= min{FFD(I),BFD(I)} pointwise
on every instance.

*Proof.* Each accepted move replaces r open bins plus the fresh-bin packing
of the leftovers by kept bins plus fewer-than-r bins in total; the move is
accepted only when the count strictly decreases, so the bin count is a
strictly decreasing integer-valued Lyapunov function and termination is
immediate. The final count is therefore <= the initial count. The worst-case
claims follow by composing with Dósa's bound (and BFD's asymptotic 11/9).
QED.

**Remark (honesty about scope).** Theorem 1 gives *safety*, not *progress*:
empirically (see experiments/repair_probe) the accepted-move count is ZERO on
all 900 standard-suite instances and all generated families tested, i.e.
min(FFD, BFD) solutions are already locally optimal under this move class.
This matches the literature: published improvements over FFD-class packings
either use reductions before packing (MFFD, Johnson & Garey 1985,
71/60·OPT+1), exact-fill combinatorial search (DJD), or heavy metaheuristic
search (Falkenauer's grouping GA; Alvim et al. 2004 hybrid improvers). Our
negative result is stronger than prior statements in one respect: it shows
even a *global* rebuild variant (emptying ALL bins below a load threshold)
fails to escape, so the local-optimum basin is robust, not an artifact of
small neighborhoods.

## 2. Exact-fit-priority variants: incomparability with FF (empirical)

Define EFF (exact-fit-first first fit): place each item into the lowest-index
bin whose residual equals the item size exactly, if any; otherwise lowest-index
feasible bin. Brute-force search over 2M random instances (n<=40, cap<=120,
both arrival orders) found:

    EFF < FF on ~1.67% of instances
    EFF > FF on ~0.018% of instances   (ratio ≈ 90:1)
    ties otherwise.

So neither algorithm dominates the other (counterexamples both ways are
committed under experiments/eff_ff_divergence), but EFF is a nearly-free
drop-in variant that wins 90x more often than it loses on random instances.
No claim of worst-case improvement is made; both are Any-Fit algorithms and
inherit the Johnson et al. Any-Fit asymptotic bounds.

## 3. Online line: analysis targets for discovered policies

The online GP policies are functions of causal features only. Two rigorous
properties will be established for the final champion P*:

**(a) Any-Fit lemma.** If P* never opens a new bin while a feasible open bin
exists (checked programmatically by instrumented simulation), then classical
arguments give CR(P*) <= 2, matching the Any-Fit family bound; if additionally
P*'s choice among feasible bins coincides with Best Fit whenever the tightest
fit is unique-minimum, the Dósa-Sgall style weight function for BF (tight
ratio floor(1.7)·OPT) can be attempted as a template.

**(b) Restricted-class dominance.** On streams where every item exceeds c/3
(bins hold at most two items), online packing reduces to greedy matching;
for this class we prove exact characterizations of the champion's behavior
vs Best Fit (see docs/proof_restricted_class.md).

**(c) Finite-state verification (planned).** For small capacities c, the
joint state (multiset of open-bin residuals, past-size histogram summary the
program actually reads) is finite up to relabeling; the champion's decision
table can then be exhaustively explored to compute its exact competitive
curve against adversarial continuations, giving machine-checked worst-case
bounds for small c. This is implemented in bpp/verify_smallcap.py.

## References

- G. Dósa. The tight bound of First Fit Decreasing bin-packing algorithm is
  FFD(I) <= (11/9)OPT(I)+6/9. COCOA 2007, LNCS.
- D.S. Johnson, M.R. Garey. A 71/60 theorem for bin packing. J. Complexity
  1(1):65-106, 1985.
- G. Dósa, J. Sgall. Optimal analysis of Best Fit bin packing. ICALP 2014.
  (BF tightly floor(1.7)·OPT.)
- M.R. Garey, R.L. Graham, J.D. Ullman / D.S. Johnson 1973-74: Any-Fit family
  asymptotic ratio <= 17/10 for FF/BF; AF class results in Johnson's MIT
  thesis (1973).
- S.F. Wu? -- (Any-Fit <= 2 statement: see Johnson thesis ch. for AF bounds.)
- B. Romera-Paredes et al. Mathematical discoveries from program search with
  large language models (FunSearch). Nature 625, 2023. [Online BPP priority
  evolution; empirical only, no proofs.]
- A.F. Alvim, C.C. Ribeiro, F. Glover, D.J. Aloise. A hybrid improvement
  heuristic for the one-dimensional bin packing problem. J. Heuristics
  10(2):205-229, 2004. [Load-redistribution operators, no guarantees.]
