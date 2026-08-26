# S4 Novelty Audit — extension of A094406 / A176762 (digit-map cycle-reach records)

Audit date: **2026-08-26**. Auditor: session S4.

## Claim under audit
1. A094406(14) = "5888" followed by N nines (N = 197-digit integer in
   `research/experiments/S4/out/a094406_certification.json`), with completeness proof.
2. Re-certification (first written minimality proof) of published A094406(13).
3. A176762(21) = the number in (1), and A176762(20) = A094406(13).

## Sources checked (all accessed 2026-08-26 unless noted)

| Source | What was checked | Finding |
|---|---|---|
| OEIS A094406, `https://oeis.org/search?q=id:A094406&fmt=json` | full entry + revision time | Entry rev. **2026-08-22**; data a(0..12); comment gives a(13) = 577999…999 (199 digits) "too large to include"; b-file (`/A094406/b094406.txt`) covers 0..13, marked "synthesized from sequence entry". **No a(14), no proof artifact.** |
| OEIS A176762, same API | full entry | Rev. **2026-07-09** (Hasler, definitional comment). Data ends a(19)=15999. Formula: a(n)=min(A001273(n), A094406(n−7)). Nothing at n ≥ 20 is displayed; n=20 follows mechanically; nothing anywhere about n=21. |
| OEIS A001273, same API | full entry + literature links | Formulas through a(19) credited to Ya-Ping Lu (**Jul 26–27 2025**) — actively maintained; NOT targeted by us. Literature: Cai & Zhou, Rocky Mountain J. Math. 38 (2008) 1921–1926; Grundman & Teeple, Fib. Quart. 41 (2003) 301–306; Havermann blog "Big and Happy" (2010); Lapointe arXiv:1904.12032; Mei & Read-McFarland arXiv:1511.01441. |
| arXiv:1904.12032 abstract (Lapointe) | scope check via WebFetch | Smallest happy numbers per height **in bases 2 and 3**, recursive relations; decimal-base records and unhappy numbers out of scope. No conflict. |
| OEIS full-text search `15999 happy unhappy steps cycle` | sibling-sequence sweep | Only A094406 and A176762 contain the frontier term; no third sequence tracking these records further. |
| Earlier S4 probe: A275154/A380448/A097757 families | target-feasibility screen | Distinct-cubes/squares exact-count records are actively maintained through 2025 with b-files of 10⁴ terms — frontiers unreachable here; documented as rejected targets in sessions/S4/REPORT.md. |
| Walter Schneider's "Unhappy Numbers" page (web.archive.org capture of wschnei.de, 2004) | historical origin of a(13)-scale values | Cited by OEIS as source; page predates any a(14). |
| General web search (DuckDuckGo html/lite endpoints; Bing via WebFetch) | hobbyist/record pages for happy-number step records | DDG returned bot-block pages (anomaly detection); Bing results were irrelevant (careers pages). **Automated general web search unavailable this session** — residual audit gap flagged below. |

## Verdict
`audited_clean` within the reachable radius: the canonical registry (OEIS) — including both
target entries revised as recently as Aug 2026 — stops exactly where our new terms begin,
and the cited literature does not go further. The most likely collision risk is the active
OEIS contributor group around A001273/A176762 (M. F. Hasler, Ya-Ping Lu, D. Corneth);
nothing they have posted extends to our claims as of 2026-08-26.

## Caveats (honesty)
* Automated general-purpose web search failed (bot blocking); novelty rests on OEIS state +
  entry-cited literature + arXiv rather than an open-web sweep.
* Concurrent discovery before OEIS submission cannot be excluded; if contributing, cite
  access date 2026-08-26.
