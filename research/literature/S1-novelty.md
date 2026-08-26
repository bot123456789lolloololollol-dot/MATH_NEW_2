# S1 Novelty Audit — maximum Sidon sets in Z_n / modular Golomb rulers

Session: S1 (extremal combinatorics / difference packing)
Audit date: 2026-08-26. Tools: OEIS (curl), arXiv API, publisher pages via OEIS links.
WebSearch tool unavailable in this session; OEIS/arXiv fetched directly over HTTPS.
Status after audit: **prior_work_found** for the family and its published tables;
our contribution is positioned as independent re-certification + saturation, not a
new exact value.

## What the family is

- Strict Sidon sets in cyclic groups Z_n ("distinct sums of any two elements"),
  equivalently modular Golomb rulers / (n,w,1) single-codeword cyclic difference
  packings. Two OEIS views: per-modulus max f(n) = A260999; inverse minimal
  modulus a(k) = A004136. Weak variant (sums of distinct elements only):
  A260998 / A004135; integer-interval variant: A004133.

## Published state found (with access dates)

| Source | Content | Access date |
|---|---|---|
| OEIS A260999 + b-file (n=1..295) + witness file | values n≤84 credited to Haanpää–Huima–Östergård (Nov 2000); n=85..295 by F. Cariboni (Dec 2017–Jan 2018); proof method for 85..295 undocumented on the page | 2026-08-26 (oeis.org/A260999, b260999.txt, a260999_2.txt) |
| H. Haanpää, A. Huima, P. R. J. Östergård, "Sets in Z_n with Distinct Sums of Pairs", Discrete Appl. Math. 138 (2004) 99–106, doi:10.1016/S0166-218X(03)00273-7 | exact computation of maximal sizes (published range to n=84 per OEIS credit) | 2026-08-26 (via OEIS refs; paywalled full text not accessed) |
| OEIS A004136 | minimal-modulus sequence, terms n=1..18 only, keyword "more"; a(19)=343 would require projective plane of order 18 (open) — confirms extension beyond k=18 is blocked by an open problem | 2026-08-26 |
| M. Buratti, D. R. Stinson, arXiv:2007.01908 (2020) "New Results on Modular Golomb Rulers, OOC and Related Structures" | settles existence of modular Golomb rulers of order k for all k ≤ 11 (minimal-modulus direction) | 2026-08-26 (arXiv abstract) |
| D. M. Gordon, arXiv:2408.16721 (2024) "Modular Golomb rulers and almost difference sets" | existence questions via almost difference sets; λ=0 case = MGR | 2026-08-26 (arXiv abstract) |
| B. Bajnok, arXiv:1705.07444 "Additive Combinatorics: A Menu of Research Problems" p.162 | survey context for Sidon-in-Z_n problems | 2026-08-26 (via A004136 link) |
| R. L. Graham, N. J. A. Sloane, "On Additive Bases and Harmonious Graphs", SIAM J. Alg. Disc. Meth. 1 (1980) 382–404 | origin of v_alpha/v_delta additive-basis tables | 2026-08-26 (OEIS refs) |

## Assessment

1. The exact values f(n) for n ≤ 84 are **established literature** (Haanpää et al.).
   Re-proving them is reproduction, not novelty.
2. The range n = 85..295 rests on Cariboni's OEIS contribution: witness sets are
   published but the entry does not document how optimality was proven (search
   method, completeness argument, or solver). Our run independently certifies
   85..110 with a written completeness argument; that is a *verification service*
   to the table-keeper, not a new value. No discrepancy found ⇒ honest saturation
   outcome within reach of our machine (mission option c).
3. No gap suitable for a NEW exact value exists at sizes reachable here:
   - extending f(n) past ~111 requires UNSAT proofs one level above the counting
     bound where published data says f stays below ub (hours-days scale here);
   - extending A004136 past k=18 is blocked by the open existence of a projective
     plane of order 18.
4. Adjacent families checked and deliberately NOT entered (would duplicate known
   work or exceed budget): graceful-labeling max-edge tables (A004137, well
   tabulated), weak-Sidon tables (A260998 complete to n=254), optimal Golomb
   rulers A003022 (k≤27 known in literature).

Conclusion: novelty_status = **prior_work_found** (family well covered to n=295
in OEIS; peer-reviewed to n=84). Any claim we make must be labeled as
re-certification/saturation, NOT new discovery.

## Negative checks performed

- OEIS search "modular golomb ruler" (4 hits, all above) — no untabulated variant.
- arXiv queries: all:"modular Golomb ruler" (15 newest), all:"Sidon sets in Z_n"
  (0 hits), all:"Sidon sets"+"cyclic group" (API returned empty/parse issues;
  noted as tool limitation, mitigated by OEIS reference harvesting which cites
  the primary literature).
