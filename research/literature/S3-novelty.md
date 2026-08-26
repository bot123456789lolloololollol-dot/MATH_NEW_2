# S3 novelty / status audit — 2026-08-26

Tooling note: this session had **WebFetch only** (no WebSearch tool). Sources
consulted: English Wikipedia ("Graceful labeling") and the arXiv API
(export.arxiv.org), queries run 2026-08-26. Where a claim below rests only on
abstract text fetched the same day, it is marked [arXiv abstract].

## Items RUN by S3

### C1. Odd-graceful tree conjecture (Gnanajothi 1991)
- Definition and conjecture: every tree with m edges admits an injective
  labeling f: V -> {0,...,2m-1} whose edge differences are exactly
  {1,3,...,2m-1}. Attributed to R. Gnanajothi, "Topics in graph theory",
  PhD thesis, Madurai Kamaraj University 1991 (as cited by Barrientos).
- Barrientos, *Odd-graceful labelings of trees of diameter 5*,
  arXiv:0807.3738 (2008): proves every tree of diameter <= 5 is odd-graceful,
  and caterpillar forests. [arXiv abstract]
- Wang-Xu-Yao arXiv:1512.08634 (2015): improper variants; no exhaustive check.
- **No published exhaustive verification "all trees up to N vertices" found**
  in the arXiv query (2026-08-26). => an exact verified-up-to-N statement at
  any N >= 6 (diameter <= 5 already covers all trees on <= ? vertices... no:
  diameter bound does NOT cover all small trees; e.g. some n=9+ trees have
  diameter >= 6) fills a gap.
- Caveat: we could not rule out verification in journals not indexed by
  arXiv (e.g., electronic journals of labeling theory); label stays
  experimentally_validated_result regardless.

### C2. Prime tree conjecture (Entringer, ~1980)
- Every tree on n vertices admits a bijection f: V -> {1..n} with gcd = 1
  across every edge. Still open per Rao's 2020 survey arXiv:2006.03801
  (lists "Entringer's prime tree conjecture" among open conjectures).
  [arXiv abstract]
- Known proven families per that survey + Fox-Spaeth arXiv:2601.14535 (2026,
  total prime labelings: snakes, books, prisms, prime trees...) and
  Cloys-Fox arXiv:1801.01802 (neighborhood-prime variants). [arXiv abstracts]
- Seoud/Youssef-type papers reportedly verify small orders in journal-only
  venues; none found via arXiv query. **No all-trees-up-to-N computational
  verification found** in our scan.

### C3. Antimagic tree conjecture (Hartsfield-Ringel 1990)
- Hartsfield-Ringel, *Pearls in Graph Theory* (1990): every connected graph
  other than K2 is antimagic (edge labels 1..m, distinct vertex sums).
- Status: "the conjecture remains open, even for trees" per arXiv:1905.06595
  (2019) and arXiv:2506.15221 (2025). [arXiv abstracts]
- Proven classes: trees with at most one degree-2 vertex (Kaplan-Lev-Roditty
  2009; Liang-Wong-Zhu 2014); caterpillars (Lozano-Mora-Seara-Tey,
  arXiv:1812.06715, Mediterr. J. Math. 2021); spiders/double spiders;
  odd-order trees with <= two degree-2 vertices (Meng, 2026, arXiv:2307.16836
  follow-ups). Meng 2026 machine-checked constructions on all 2,245,070
  odd-order trees up to 25 **with exactly two degree-2 vertices** only.
  [arXiv abstracts]
- **No exhaustive check of ALL trees up to any order bound found** =>
  verified-up-to-N over complete tree levels is new as far as our scan goes.

### C4. Neighborhood-prime tree conjecture (Ryan ~2014)
- Bijection f: V -> {1..n} such that gcd of the labels in N(v) is 1 for every
  v with deg(v) >= 2; conjectured to exist for every tree of order >= 3
  (S. Ryan et al., ~2014; see Cloys & Fox, arXiv:1801.01802, "Neighborhood
  prime labelings of trees", which prove caterpillars, spiders, firecrackers,
  polygonal snakes, books). [arXiv abstract]
- **No exhaustive all-trees-up-to-N verification found.**

## Items SKIPPED after status check (verified beyond our reach)

- Graceful tree conjecture (Ringel 1963 / Rosa 1967): OPEN but computationally
  verified for ALL trees with <= 27 vertices (Aldred & McKay 1998),
  <= 29 (Horton Honours thesis), claimed <= 35 (Fang's Graceful Tree
  Verification Project, 2010) — source: Wikipedia "Graceful labeling",
  retrieved 2026-08-26. Out of budget.
- Harmonious tree conjecture (Graham-Sloane 1980): Fang, *New Computational
  Result on Harmonious Trees* (arXiv, 2011/2012): every tree with at most
  **31 nodes** is harmonious. Also asymptotically settled by Montgomery-
  Pokrovskiy-Sudakov (2018, rainbow embeddings) and bounded-degree case fully
  resolved by Muyesser-Pokrovskiy (2025). Out of budget.

## Novelty verdict for S3 outputs

For each of C1-C4 the specific artifact "exact backtracking verification that
every unlabeled free tree on n <= N vertices admits the labeling, with
per-order witness samples and OEIS-validated complete enumeration" was NOT
found in the scanned literature; the nearest prior art per item is listed
above (class proofs, not exhaustive sweeps). Novelty status: audited_clean
for the verified-up-to-N statements themselves, with the caveat that
journal-only venues were outside scan reach (WebFetch-only session).

Dates: all queries executed 2026-08-26.
