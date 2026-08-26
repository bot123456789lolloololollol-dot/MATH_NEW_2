#!/usr/bin/env python3
"""Freeze ROL parameter K*(m) from DEV results ONLY (preregistered selection
rule: per m, choose K_mult minimizing mean gap on the dev configuration of
that m; tie -> smaller K_mult). Writes outputs/K_frozen.json and prints the
tuning table. Must be run after --phase dev and before --phase eval."""
import csv
import json
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "outputs")

mean_gap = defaultdict(list)  # (m, k) -> [gaps]
with open(os.path.join(OUTDIR, "results_dev.csv")) as f:
    for row in csv.DictReader(f):
        if row["algo"] != "rol":
            continue
        mean_gap[(int(row["m"]), int(row["K_mult"]))].append(float(row["gap_pct"]))

ms = sorted({m for (m, _) in mean_gap})
frozen = {}
print(f"{'m':>3} {'K/m=1':>8} {'K/m=2':>8} {'K/m=4':>8} {'chosen K_mult':>14}")
for m in ms:
    stats = {k: sum(v) / len(v) for (mm, k), v in mean_gap.items() if mm == m}
    best = min(stats, key=lambda k: (stats[k], k))  # min mean gap; tie -> smaller K
    frozen[m] = best
    print(f"{m:>3} " + " ".join(f"{stats.get(k, float('nan')):>8.4f}" for k in (1, 2, 4))
          + f" -> K={best}*{m}={best*m}")
    for k in (1, 2, 4):
        v = mean_gap[(m, k)]
        assert len(v) == 30, f"expected 30 dev runs for m={m},K={k}, got {len(v)}"

with open(os.path.join(OUTDIR, "K_frozen.json"), "w") as f:
    json.dump({str(m): frozen[m] for m in ms}, f, indent=1)
print("frozen:", json.dumps(frozen))
