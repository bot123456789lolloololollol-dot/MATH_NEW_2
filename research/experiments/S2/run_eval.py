#!/usr/bin/env python3
"""S2 evaluation runner (deterministic, paired design).

--phase dev : runs on DEV_SEEDS x dev-configs D; ROL evaluated at every grid
              K_mult in {1,2,4}; writes outputs/results_dev.csv
--phase eval: reads outputs/K_frozen.json (written by freeze_K.py from dev data
              ONLY) and runs LPT / MULTIFIT / LPT+LS / ROL(K*) on
              EVAL_SEEDS x all 30 configs; writes outputs/results_eval.csv

Raw CSV columns:
config,m,n,seed,algo,K_mult,cmax,lb,gap_pct,time_ms,assign_sha16
"""
import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from algorithms import lb_bound, lpt, lpt_ls, multifit, rol  # noqa: E402
from gen_instances import DEV_SEEDS, EVAL_SEEDS, MS, NS, OUT, read_instance  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "outputs")
os.makedirs(OUTDIR, exist_ok=True)

DEV_CONFIGS = [(2, 30), (3, 50), (4, 100), (5, 100), (7, 200), (10, 500)]
K_GRID = [1, 2, 4]  # K = K_mult * m


def assign_sha(assign) -> str:
    h = hashlib.sha256()
    for i, mach in enumerate(assign):
        h.update((f"{i}:" + ",".join(map(str, sorted(mach))) + ";").encode())
    return h.hexdigest()[:16]


def run_one(algo: str, p: np.ndarray, m: int, k_mult: int):
    t0 = time.perf_counter()
    if algo == "lpt":
        assign, cmax = lpt(p, m)
    elif algo == "multifit":
        assign, cmax = multifit(p, m)
    elif algo == "lpt_ls":
        assign, cmax = lpt_ls(p, m)
    elif algo == "rol":
        assign, cmax = rol(p, m, k_mult)
    else:
        raise ValueError(algo)
    dt_ms = (time.perf_counter() - t0) * 1000.0
    return assign, int(cmax), dt_ms


def rows_for(phase: str):
    if phase == "dev":
        configs = DEV_CONFIGS
        seeds = DEV_SEEDS
    else:
        configs = [(m, n) for m in MS for n in NS]
        seeds = EVAL_SEEDS
    with open(os.path.join(OUTDIR, "K_frozen.json")) as f:
        kfrozen = json.load(f)
    plan = [("lpt", None), ("multifit", None), ("lpt_ls", None)]
    if phase == "dev":
        plan += [("rol", km) for km in K_GRID]
    else:
        plan += [("rol", None)]  # uses per-m frozen K*

    out_rows = []
    for (m, n) in configs:
        for seed in seeds:
            path = os.path.join(OUT, f"taillike-m{m}-n{n}-s{seed}.txt")
            mm, nn, sd, p = read_instance(path)
            assert (mm, nn, sd) == (m, n, seed), path
            lb = lb_bound(p, m)
            for algo, km in plan:
                if phase == "eval" and algo == "rol":
                    km = kfrozen[str(m)]
                assign, cmax, dt = run_one(algo, p, m, km if km else 1)
                gap = 100.0 * (cmax - lb) / lb
                out_rows.append(
                    f"{phase}-m{m}-n{n},{m},{n},{seed},{algo},"
                    f"{km if km else ''},{cmax},{lb},{gap:.6f},{dt:.3f},"
                    f"{assign_sha(assign)}"
                )
        print(f"[{phase}] done config m={m} n={n}", flush=True)
    return out_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=["dev", "eval"], required=True)
    args = ap.parse_args()
    header = ("config,m,n,seed,algo,K_mult,cmax,lb,gap_pct,time_ms,assign_sha16")
    rows = rows_for(args.phase)
    fn = os.path.join(OUTDIR, f"results_{args.phase}.csv")
    with open(fn, "w") as f:
        f.write(header + "\n")
        f.write("\n".join(rows) + "\n")
    print(f"wrote {len(rows)} rows -> {fn}")


if __name__ == "__main__":
    main()
