# S6-MATH Novelty Audit Log

(Session "S6-MATH" = mathematical discovery session; distinct from the other
session that also claimed the S6 label for a physics/SINDy study.)

## Audit 1 — Candidate C-S6M-001 (strict/distinct-variable Rado numbers, x+y=kz, 2 colors)
Date: 2026-08-26. Agent: general-purpose subagent, 50 tool calls
(arXiv API, Crossref, OpenAlex, Semantic Scholar, OEIS, Bing/DDG).

VERDICT: PRIOR ART FOUND — CANDIDATE KILLED.
- A. Towell, "Distinct-Variable Rado Numbers for x+y=kz: Computation,
  Closed-Form Formula, and Complete Proof", Zenodo DOI 10.5281/zenodo.19372727
  (2026-04): identical table (k ≤ 500), identical parity-split closed form for
  k ≥ 8 ((k²+2k+2)/2 even k, (k²+3k)/2 odd k), identical exceptions k ≤ 7,
  complete claimed proof.
- OEIS A394445 (Towell, approved 2026-03): same sequence + formula verbatim.
- Non-strict family classical: OEIS A100542; Burr–Loo; Harborth–Maasberg,
  Discrete Math. 197–198 (1999) 397–407; Jones–Schaal, Discrete Math. 289
  (2004) 63–69.

Our independent computation is retained as validation evidence for the tooling
and as an independent confirmation of Towell's table.
Post-mortem: research/sessions/S6/postmortems/C-001.md;
candidate record: research/candidates/C-S6M-001.json.

## Audit 2 — Candidate C-S6M-002 (per-modulus classification of weak 2-colorings
of Z_n avoiding mono distinct x+y ≡ kz)
Date: 2026-08-26. Agent: general-purpose subagent, 81 tool calls.

VERDICT: NOVELTY SURVIVES for k ≥ 3; the k = 2 case is prior work and must be
cited, not claimed.

KILLED as novelty (sub-case):
- Grier, "On the cyclic van der Waerden numbers", Geombinatorics 21(4) (2012):
  G(3;2) = {1,2,3,4,6,8} — exactly our k=2 colorable set (minus n=1), within
  Burkert–Johnson's W_c framework (Soifer ed., Progr. Math. 285, Birkhäuser
  2010, pp. 97–113).

NOT FOUND anywhere (arXiv API, OpenAlex, Crossref, DDG; OEIS partially blocked):
any per-modulus classification of 2-colorability of Z_n avoiding monochromatic
pairwise-distinct solutions of x+y ≡ kz for fixed k ≥ 3 or across the family.

Adjacent but DIFFERENT objects (examined in detail):
- Abbott–Wang, Proc. AMS 67 (1977) 11–16: interval {1..n} colored, modulus n+1,
  coefficient 1 only.
- Chappelon–Revuelta Marchena–Sanz Domínguez, "Modular Schur numbers",
  Electron. J. Combin. 20(2) (2013): weakly l-sum-free modulo m on initial
  segments {1,…,n} with FIXED modulus (exact values only m=1,2,3) — the dual
  "fixed small modulus" direction.
- D'orville–Sim–Wong–Ho, Ramanujan J. 67(4) (2025),
  doi:10.1007/s11139-025-01123-5: extends WS_m(k,l); same interval setting.
- Sanders, arXiv:2604.23738 (2026): asymptotic bounds for modular Schur-type
  thresholds over nonzero residues; no per-modulus classification, no
  distinctness condition; explicitly shows coefficient variants do not reduce
  to the a=1 theory.
- Zero-sum generalized Schur numbers (Robertson 2018; Robertson–Roy–Sarkar
  2018; Roy–Sarkar 2018): color-sum zero-sum condition on integer intervals —
  different condition entirely.
- Rainbow-free / anti-van der Waerden colorings of abelian groups: forbids
  RAINBOW solutions with ≥3 colors — complementary condition.

Caveats: OEIS search partially unreachable (HTTP 403 via fetchers); theses and
non-indexed preprints cannot be fully excluded by automated audit.
