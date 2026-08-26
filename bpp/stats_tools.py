"""Paired statistical tests + effect sizes for policy comparisons."""
import numpy as np
from scipy import stats


def paired_report(a, b, label_a="A", label_b="B"):
    """a, b: per-instance metric (bins) arrays, same length. Lower is better."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    d = a - b
    nz = d[d != 0]
    out = {
        "n": int(len(d)),
        "mean_a": float(a.mean()), "mean_b": float(b.mean()),
        "total_a": float(a.sum()), "total_b": float(b.sum()),
        "wins_a": int((d < 0).sum()), "wins_b": int((d > 0).sum()),
        "ties": int((d == 0).sum()),
    }
    if len(nz) >= 5:
        w = stats.wilcoxon(a, b, zero_method="wilcox")
        out["wilcoxon_stat"] = float(w.statistic)
        out["wilcoxon_p"] = float(w.pvalue)
    else:
        out["wilcoxon_p"] = None
    if d.std() > 0:
        out["cohens_dz"] = float(d.mean() / d.std())
        t = stats.ttest_rel(a, b)
        out["ttest_p"] = float(t.pvalue)
    else:
        out["cohens_dz"] = 0.0
        out["ttest_p"] = None
    # bootstrap CI for mean difference
    rng = np.random.default_rng(0)
    if len(d) >= 10:
        boots = [float(rng.choice(d, size=len(d), replace=True).mean())
                 for _ in range(2000)]
        boots.sort()
        out["boot_ci95"] = [boots[50], boots[1949]]
    return out
