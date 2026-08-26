#!/usr/bin/env bash
# S2 full pipeline. Run from this directory. Total wall-clock < 1 minute on a
# typical laptop (Python 3.11, numpy, scipy required).
set -euo pipefail
cd "$(dirname "$0")"

# 0. validation gates (hand instance, brute-force cross-check, determinism)
python validate.py

# extra robustness: vectorized LS vs independent naive steepest descent
python check_ls.py

# 1. generate committed benchmark instances (skips existing files; deterministic)
python gen_instances.py

# 2. dev-phase tuning on DEV_SEEDS x D (preregistered; see PREREGISTERED.md)
python run_eval.py --phase dev
python freeze_K.py | tee outputs/tuning_freeze.txt

# 3. final preregistered evaluation on EVAL_SEEDS x 30 configs
python run_eval.py --phase eval

# 4. statistics + integrity checks
python analyze.py

echo "Done. Raw outputs in outputs/:"
ls -la outputs/
