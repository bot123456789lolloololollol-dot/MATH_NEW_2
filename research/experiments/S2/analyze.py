#!/usr/bin/env python3
"""S2 statistical analysis (preregistered protocol).

Integrity checks:
  A. every row's LB matches max(max p, ceil(sum p/m)) recomputed from the
     committed instance file;
  B. gap_pct == 100*(cmax-lb)/lb;
  C. independent pure-python (no numpy) LPT reimplementation reproduces the
     recorded lpt cmax on a fixed sample of 40 rows.
Statistics per config and pooled:
  paired Wilcoxon signed-rank (two-sided, zero_method='wilcox'), Holm-Bonferroni
  across the 3 comparisons {rol vs lpt, rol vs multifit, rol vs lpt_ls} within
  each config; paired t-test sensitivity; effect sizes rank-biserial r_rb and
  Cohen's d_z. Writes outputs/stats_summary.json and outputs/analysis.txt.
"""
import csv
import json
import math
import os
from collections import defaultdict

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "outputs")
INST = os.path.normpath(os.path.join(HERE, "..", "..", "benchmarks", "S2-instances"))

ALGOS = ["lpt", "multifit", "lpt_ls", "rol"]
COMPARISONS = [("rol", "lpt"), ("rol", "multifit"), ("rol", "lpt_ls")]


def load_rows(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            r["m"], r["n"], r["seed"] = int(r["m"]), int(r["n"]), int(r["seed"])
            for k in ("cmax", "lb"):
                r[k] = int(r[k])
            r["gap"] = float(r["gap_pct"])
            rows.append(r)
    return rows


def check_integrity(rows):
    bad = 0
    for r in rows:
        with open(os.path.join(INST, f"taillike-m{r['m']}-n{r['n']}-s{r['seed']}.txt")) as f:
            parts = f.read().split()
        m, n = int(parts[0]), int(parts[1])
        p = list(map(int, parts[3:]))
        assert len(p) == n
        lb = max(max(p), -(-sum(p) // m))
        if lb != r["lb"]:
            bad += 1
            print(f"INTEGRITY FAIL lb row {r}")
        g = 100.0 * (r["cmax"] - r["lb"]) / r["lb"]
        if abs(g - r["gap"]) > 1e-6:
            bad += 1
            print(f"INTEGRITY FAIL gap row {r}")
    return bad


def indep_lpt(p, m):
    order = sorted(range(len(p)), key=lambda j: (-p[j], j))
    loads = [0] * m
    for j in order:
        i = min(range(m), key=lambda t: (loads[t], t))
        loads[i] += p[j]
    return max(loads)


def check_indep_lpt(rows):
    rng = np.random.default_rng(20260826)
    cand = [r for r in rows if r["algo"] == "lpt"]
    sample = [cand[i] for i in rng.choice(len(cand), size=40, replace=False)]
    bad = 0
    for r in sample:
        with open(os.path.join(INST, f"taillike-m{r['m']}-n{r['n']}-s{r['seed']}.txt")) as f:
            parts = f.read().split()
        p = list(map(int, parts[3:]))
        if indep_lpt(p, r["m"]) != r["cmax"]:
            bad += 1
            print(f"INDEP-LPT MISMATCH {r}")
    return bad


def paired(a_gaps, b_gaps):
    d = np.array(a_gaps) - np.array(b_gaps)
    out = {"mean_a": float(np.mean(a_gaps)), "mean_b": float(np.mean(b_gaps)),
           "std_a": float(np.std(a_gaps, ddof=1)), "std_b": float(np.std(b_gaps, ddof=1)),
           "mean_diff": float(np.mean(d)), "std_diff": float(np.std(d, ddof=1))}
    nz = d[d != 0]
    npos = int((d > 0).sum()); nneg = int((d < 0).sum())
    try:
        w = stats.wilcoxon(a_gaps, b_gaps, zero_method="wilcox", alternative="two-sided")
        out["wilcoxon_p"] = float(w.pvalue)
    except ValueError:  # all differences zero
        out["wilcoxon_p"] = 1.0
    if not math.isfinite(out["wilcoxon_p"]):
        # scipy returns nan for the degenerate all-zero-differences case;
        # convention: identical paired samples -> p = 1.0
        out["wilcoxon_p"] = 1.0
        out["degenerate_all_zero_diffs"] = True
    ranks = stats.rankdata(np.abs(d))
    rp, rn = ranks[d > 0].sum(), ranks[d < 0].sum()
    out["rank_biserial"] = float((rp - rn) / ranks.sum())  # >0 favors a smaller? see note
    out["dz"] = out["mean_diff"] / out["std_diff"] if out["std_diff"] > 0 else 0.0
    t = stats.ttest_rel(a_gaps, b_gaps)
    out["ttest_p"] = float(t.pvalue)
    if not math.isfinite(out["ttest_p"]):
        out["ttest_p"] = 1.0
    out["n_pos_diff"] = npos  # a>b count (a worse)
    out["n_neg_diff"] = nneg  # a<b count (a better)
    out["n_zero_diff"] = int(len(d) - npos - nneg)
    return out


def holm(pvals):
    """Holm-Bonferroni adjusted p-values (same order as input)."""
    m = len(pvals)
    idx = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    prev = 0.0
    for rank, i in enumerate(idx):
        v = min(1.0, (m - rank) * pvals[i])
        v = max(v, prev)
        adj[i] = v
        prev = v
    return adj


def main():
    rows = load_rows(os.path.join(OUTDIR, "results_eval.csv"))
    assert len(rows) == 3600, len(rows)

    print("== integrity ==")
    bad = check_integrity(rows)
    bad += check_indep_lpt(rows)
    print(f"integrity failures: {bad}")

    by = defaultdict(dict)
    for r in rows:
        by[(r["m"], r["n"])][r["algo"]] = by[(r["m"], r["n"])].get(r["algo"], [])
        by[(r["m"], r["n"])][r["algo"]].append(r["gap"])

    report = {"integrity_failures": bad, "configs": {}, "pooled": {}}
    txt = []
    txt.append(f"{'config':>12} | {'LPT mean±sd':>16} | {'MF mean±sd':>16} | "
               f"{'LS mean±sd':>16} | {'ROL mean±sd':>16} | ROL vs each (p_holm, r_rb, dz)")
    txt.append("-" * 150)
    wins = {c: 0 for c in COMPARISONS}
    sig_counts = {c: 0 for c in COMPARISONS}
    for cfg in sorted(by):
        gaps = {a: by[cfg][a] for a in ALGOS}
        for a in ALGOS:
            assert len(gaps[a]) == 30, (cfg, a, len(gaps[a]))
        res = {}
        padj = holm([paired(gaps[x], gaps[y])["wilcoxon_p"] for (x, y) in COMPARISONS])
        for (x, y), ph in zip(COMPARISONS, padj):
            pr = paired(gaps[x], gaps[y])
            pr["wilcoxon_p_holm"] = ph
            res[f"{x}_vs_{y}"] = pr
            if ph < 0.05:
                sig_counts[(x, y)] += 1
                if pr["mean_diff"] < 0:
                    wins[(x, y)] += 1  # significant improvement for x over y
        means = {a: f"{np.mean(gaps[a]):.4f}±{np.std(gaps[a], ddof=1):.4f}" for a in ALGOS}
        comp_str = "; ".join(
            f"{x}>{y}: p={res[f'{x}_vs_{y}']['wilcoxon_p_holm']:.2g}, "
            f"d={res[f'{x}_vs_{y}']['mean_diff']:+.4f}, r_rb={res[f'{x}_vs_{y}']['rank_biserial']:+.2f}, "
            f"dz={res[f'{x}_vs_{y}']['dz']:+.2f}"
            for x, y in COMPARISONS)
        txt.append(f"{cfg[0]:>3},{cfg[1]:>4}   | {means['lpt']:>16} | {means['multifit']:>16} | "
                   f"{means['lpt_ls']:>16} | {means['rol']:>16} | {comp_str}")
        report["configs"][f"{cfg[0]}-{cfg[1]}"] = res

    # pooled across all configs
    pooled = {}
    for x, y in COMPARISONS:
        xs, ys = [], []
        for cfg in sorted(by):
            xs.extend(by[cfg][x]); ys.extend(by[cfg][y])
        pr = paired(xs, ys)
        pr["wilcoxon_p_holm"] = pr["wilcoxon_p"]  # single family test
        pooled[f"{x}_vs_{y}"] = pr
    report["pooled"] = pooled
    report["sig_config_counts"] = {f"{x}_vs_{y}": {"significant": sig_counts[(x, y)],
                                                   "significant_and_favorable": wins[(x, y)],
                                                   "of": len(by)} for x, y in COMPARISONS}

    txt.append("")
    txt.append("== POOLED (900 paired runs) ==")
    for k, v in pooled.items():
        txt.append(f"{k}: mean {v['mean_b']:.4f}->{v['mean_a']:.4f}, diff {v['mean_diff']:+.5f}, "
                   f"Wilcoxon p={v['wilcoxon_p']:.3g}, r_rb={v['rank_biserial']:+.3f}, dz={v['dz']:+.3f}, "
                   f"W/L/T(better/worse/tie)={v['n_neg_diff']}/{v['n_pos_diff']}/{v['n_zero_diff']}")
    txt.append("")
    txt.append(f"Config counts (Holm p<0.05): {json.dumps(report['sig_config_counts'])}")

    out_txt = "\n".join(txt)
    print(out_txt)
    with open(os.path.join(OUTDIR, "analysis.txt"), "w") as f:
        f.write(out_txt + "\n")
    with open(os.path.join(OUTDIR, "stats_summary.json"), "w") as f:
        json.dump(report, f, indent=1)


if __name__ == "__main__":
    main()
