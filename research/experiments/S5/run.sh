#!/usr/bin/env bash
# S5 experiment regeneration script.
# Reproduces: encoder validation (s(n) for n=4..6 + merge/median variant optima).
# Heavy runs (n=7 UNSAT@15, n=8 pair, med(7)) take minutes to hours; see REPORT.md
# for recorded timings. Run from this directory.
set -e

echo "== checker self-test (independent verifier) =="
python checker.py --selftest

echo "== sorting-network validation: UNSAT at s(n)-1, SAT at s(n) =="
for spec in "4 4" "4 5" "5 8" "5 9" "6 11" "6 12"; do
  set -- $spec; n=$1; T=$2
  python sat_net.py cnf --n $n --slots $T --out cnf_n${n}_T${T}.dimacs
  python solve_run.py --dimacs cnf_n${n}_T${T}.dimacs --engine cadical --timeout 300 \
      --model-out net_n${n}_T${T}.json --n $n --slots $T --variant sort || true
done
python checker.py --net net_n4_T5.json --property sort
python checker.py --net net_n5_T9.json --property sort
python checker.py --net net_n6_T12.json --property sort

echo "== n=7: SAT side of s(7)=16 (UNSAT@15 exceeds single-machine budget; see REPORT) =="
python sat_net.py cnf --n 7 --slots 16 --out cnf_n7_T16.dimacs
python solve_run.py --dimacs cnf_n7_T16.dimacs --engine cadical --timeout 3000 \
    --model-out net_n7_T16.json --n 7 --slots 16 --variant sort || true
python checker.py --net net_n7_T16.json --property sort

echo "== merging networks M(m,nm) exact optima (small cases) =="
# M(2,2)=3 : UNSAT@2 SAT@3 ; M(2,3)=5 ; M(3,3)=6 ; M(2,4)=6
for spec in "4 2 2" "4 2 3" "5 2 3" "5 2 4" "5 2 5" "6 3 3" "6 3 4" "6 3 5" "6 3 6" "6 2 4" "6 2 5" "6 2 6"; do
  set -- $spec; n=$1; m=$2; T=$3
  python sat_net.py cnf --n $n --slots $T --variant merge --m $m \
      --out cnf_merge_n${n}_m${m}_T${T}.dimacs
  python solve_run.py --dimacs cnf_merge_n${n}_m${m}_T${T}.dimacs --engine cadical \
      --timeout 300 || true
done
python solve_run.py --dimacs cnf_merge_n4_m2_T3.dimacs --engine cadical --timeout 60 \
    --model-out net_merge_n4_m2_T3.json --n 4 --slots 3 --variant merge || true
python checker.py --net net_merge_n4_m2_T3.json --property merge --m 2

echo "== median networks med(n): UNSAT@s-1, SAT@s (SB2 must be OFF) =="
for T in 6 7; do
  python sat_net.py cnf --n 5 --slots $T --variant median --no-sb2 \
      --out cnf_med_n5_T${T}.dimacs
  python solve_run.py --dimacs cnf_med_n5_T${T}.dimacs --engine cadical --timeout 300 || true
done
python solve_run.py --dimacs cnf_med_n5_T7.dimacs --engine cadical --timeout 60 \
    --model-out net_med_n5_T7.json --n 5 --slots 7 --variant median || true
python checker.py --net net_med_n5_T7.json --property median

echo "== second-engine cross-checks of key UNSATs =="
python solve_run.py --dimacs cnf_n5_T8.dimacs --engine glucose --timeout 300 || true
python solve_run.py --dimacs cnf_merge_n6_m3_T5.dimacs --engine glucose --timeout 300 || true
python solve_run.py --dimacs cnf_med_n5_T6.dimacs --engine glucose --timeout 300 || true

echo "== z3 independent-engine spot checks (smallest cases) =="
python solve_run.py --engine z3 --z3-native --n 4 --slots 4 --variant sort --timeout 120 || true
python solve_run.py --engine z3 --z3-native --n 4 --slots 5 --variant sort --timeout 120 || true
python solve_run.py --engine z3 --z3-native --n 5 --slots 8 --variant sort --timeout 600 || true

echo "Done. See REPORT.md for full results table and timings."
