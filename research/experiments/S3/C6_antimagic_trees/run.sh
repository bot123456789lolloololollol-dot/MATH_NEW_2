#!/usr/bin/env bash
# S3 retry session - exhaustive anti verification, chunked (<=300s per invocation).
# Logs append to log_anti.jsonl; a level counts only if "complete": true.
set -e
cd "$(dirname "$0")"
python ../run_labelings.py --selftest
python ../run_labelings.py --mode anti --nmin 3 --nmax 19 --deadline 280 --outdir .
