#!/usr/bin/env bash
# S3 retry session - exhaustive nprime verification, chunked (<=300s per invocation).
# Logs append to log_nprime.jsonl; a level counts only if "complete": true.
set -e
cd "$(dirname "$0")"
python ../run_labelings.py --selftest
python ../run_labelings.py --mode nprime --nmin 3 --nmax 19 --deadline 280 --outdir .
