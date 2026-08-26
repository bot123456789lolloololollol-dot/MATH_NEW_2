# S5 Novelty Audit — SAT-based size-optimal comparator networks
Auditor: Session S5. Access dates: 2026-08-26 (all web sources).

## A. Plain sorting networks s(n) — frontier as verified today

| source | date accessed | content |
|---|---|---|
| https://en.wikipedia.org/wiki/Sorting_network | 2026-08-26 | Exact minimal sizes proven for n = 1..12; best-known sizes n=13..16 are 45/51/56/60; depths proven optimal n <= 16 (Bundala & Zavodny 2014); 0-1 principle statement. |
| https://bertdobbelaere.github.io/sorting_networks.html | 2026-08-26 (page last updated Nov 7, 2025) | Living table of best-known sizes/depths to n=32; marks proven optimality per value; n=13 size listed as bounded in 44..45 (lower bound via Van Voorhis-style induction), n=14..16 bounded 47/48+, 51+, 55+. |
| arXiv API query (submittedDate desc) | 2026-08-26 | Jannis Harder, "An Answer to the Bose-Nelson Sorting Problem for 11 and 12 Channels", arXiv:2012.04400 v3 (submitted Dec 8, 2020): proves s(11)=35, s(12)=39 by generalized Van Voorhis DP + Isabelle/HOL-verified checker. |
| arXiv:1405.5754 v3 | 2026-08-26 | Codish, Cruz-Filipe, Frank, Schneider-Kamp, "Twenty-Five Comparators is Optimal when Sorting Nine Inputs (and Twenty-Nine for Ten)" (2014): SAT + symmetry proof of s(9)=25, s(10)=29. |
| arXiv:1411.6408, arXiv:1707.08725 | 2026-08-26 | "Sorting Networks: the End Game" (2014) front/back-end pruning; Frasinaru & Raschip 2017 subsumption speedups — confirm the machinery class used for proven optima. |
| Bundala & Zavodny, "Optimal-Size Sorting Networks" (LPAR 2014) / Wikipedia corroboration | 2026-08-26 | Depth optimality d(n) for n<=16; SAT-based methodology ancestor of ours. |

**Conclusion:** smallest genuinely open case is **s(13) in [44?,45]** (Wikipedia's fetch reported a lower bound of 43 from Van Voorhis-style induction anchored at S(11)=35; Dobbelaere's page suggests 44; both sources agree the upper bound is 45 and that the case is open). Closing it required industrial-scale compute (cf. Harder 2020 for n=11,12) — out of scope for this session's budget; we document scaling limits instead (see REPORT.md).

## B. Variant prior art (merging networks, median/selection networks)

Search attempts on 2026-08-26:
- arXiv API `abs:"selection networks"` etc.: only ML architecture hits; no classical selection-network optimality papers surfaced.
- DuckDuckGo HTML: bot-blocked (CAPTCHA), unusable.
- Bing (two phrasings): no relevant hits for merging-network minimality; for median networks the results were generic statistics pages, but the retrieval model's own summary volunteered "the well-known result that a median-of-5 network needs 7 comparators", consistent with med(5)=7 being FOLKLORE-KNOWN. We therefore treat med(5)=7 as PRIOR ART (folklore), not a novel result of ours.
- No authoritative modern table of exact minimal merging-network sizes M(m,n) was located within the search budget. Classical constructions (Batcher odd-even mergers, Knuth TAOCP vol. 3 §5.3.4) provide upper bounds; we could not verify whether M(2,2)=3, M(2,3)=5, M(3,3)=6, M(2,4)=6 have been published as PROVEN minima before. They are small values likely known or derivable; we label them accordingly (experimentally_validated_result here; novelty unknown/partially-known).

Honest assessment: this session did NOT find (within its search budget) published SAT-based exact-minima proofs for these specific variant parameter points; but given their size, absence of evidence is weak. Claims are framed as certificates + reproducible pipeline, not as record claims.

## C. What is genuinely attributable to this session

1. A validated-from-scratch open SAT encoder (DIMACS, solver-agnostic) reproducing known s(n) optima for n<=6 (n=7,8 partially; see REPORT) — an engineering artifact with documented debugging methodology (trajectory clause-diffing).
2. An empirical soundness finding: the standard "relabel first comparator to (1,2)" symmetry breaker is UNSOUND for median-selection variants (designated-output property) — demonstrated by wrong UNSAT at med(5)=7 with SB2 ON vs correct SAT with SB2 OFF. Instructive for anyone reusing sorting-network symmetry breakers on selection properties.
3. Certified exact values (SAT UNSAT below + explicit network above + independent brute-force verification): M(2,2)=3, M(2,3)=5, M(3,3)=6, M(2,4)=6, med(5)=7 (+ any further values recorded in REPORT.md).
4. Scaling data for the plain direct encoding up to n=7/8 (documented limits).

---

## D. Second-pass audit addendum (independent S5 workstream, 2026-08-26 ~02:00)

A second S5 run pursued evolutionary campaigns + a general-k selector exact
search (complementary to sections A-C above). Additional sources checked:

| source | date | relevance |
|---|---|---|
| bertdobbelaere.github.io/median_networks.html | 2026-08-26 | best-known median sizes N=3..20; proven ONLY N=3, N=9 (=19, Smith96); convention verified empirically = bottom floor((N-1)/2) onto wires 0..k-1 |
| github.com/bertdobbelaere/SorterHunter (544 network JSONs) | 2026-08-26 | record networks downloaded, used as seeds AND independently re-verified (see seed_sweep.txt) |
| Google/Bing/DuckDuckGo queries for selection-network exact tables | 2026-08-26 | no published proven-optimal or best-known table for the general k-selector property SEL(n,k) found |

Novelty positions claimed by THIS workstream:
1. Proven-optimal sizes for SEL(n,k) at small n (SAT UNSAT certificates +
   verified witnesses): no prior art found. Values recorded in
   experiments/S5/results_exact.jsonl.
2. Improved upper bounds for SEL(n,k) at (10,5),(12,6),(14,7),(16,8) vs
   derived baselines (pruned record sorters / adapted median nets).
3. Negative result: published sorter records at N=13..17 did not fall to
   ~40 CPU-minutes of hybrid evolution+repair search (consistent with the
   maturity documented in section A).
