#!/usr/bin/env bash
# S4 -- certified extension of OEIS A094406 / A176762 (sum-of-squares-of-digits map)
# Regenerates every artifact in out/ from scratch.  Deterministic; runtime < 2 min.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p out

python happy_core.py            # reproduce published terms A176762 a(0..19), A094406 a(0..12)
python extend.py                # certify a'(13), compute NEW a'(14), derive A176762(21)
python checker_independent.py   # independent second-algorithm verification + tamper suite

echo "ALL S4 EXPERIMENTS COMPLETE -- see out/published_terms.json,"
echo "out/a094406_certification.json, out/check_report.json"
