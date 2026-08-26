#!/usr/bin/env bash
# run.sh — reproduce every S1 experiment end-to-end (deterministic).
# Total wall time ~20-25 min on a desktop machine; each step < 10 min.
set -e
cd "$(dirname "$0")"

echo "== Step 0: checker corruption selftest =="
python check_sidon.py --selftest

echo "== Step 1: test battery (naive ground truth + cross-impl agreement) =="
python run_tests.py --agree-to 30 --naive-to 14 --js-deadline-ms 30000

echo "== Step 2: main exact-table computation, chunks A/B/C (n=2..110) =="
python run.py solve --lo 2   --hi 90  --tag chunkA --deadline-ms 300000
python run.py solve --lo 91  --hi 105 --tag chunkB --deadline-ms 300000
python run.py solve --lo 106 --hi 110 --tag chunkC --deadline-ms 420000

echo "== Step 3: verification phase (checker on witnesses, b-file compare,"
echo "          Cariboni witness audit, inverse A004136 compare) =="
python run.py verify --hi 110 \
  --results main_chunkA.json main_chunkB.json main_chunkC.json

echo "== Step 4: summary CSV =="
python - <<'EOF'
import json
merged = {}
for f in ("main_chunkA.json", "main_chunkB.json", "main_chunkC.json"):
    for r in json.load(open("outputs/" + f)):
        merged[r["n"]] = r
with open("outputs/summary_n2_110.csv", "w") as fo:
    fo.write("n,counting_ub,k_max_certified,witness_mod_n\n")
    for n in sorted(merged):
        r = merged[n]
        assert r["certified"], n
        fo.write('%d,%d,%d,"%s"\n' % (n, r["ub"], r["k_max"],
                                      ";".join(map(str, r["witness"]))))
print("wrote outputs/summary_n2_110.csv,", len(merged), "rows")
EOF

echo "== ALL DONE: expected verdict = verification CLEAN, no b-file mismatches =="
