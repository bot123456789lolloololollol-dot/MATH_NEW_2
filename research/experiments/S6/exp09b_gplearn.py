"""Exp 09b — external-library GP baseline (round 2, preregistered R2).

gplearn 0.4.3 (a widely used, independent GP symbolic-regression implementation)
on the same three exp03 function classes, 30 paired seeds. Metrics: holdout NMSE
and expression size; paired Wilcoxon against our MDL-selected GP.
Hypothesis: comparable accuracy with substantially larger expressions (no MDL front).
"""
import sys
import numpy as np

sys.path.insert(0, ".")
from src import symreg as sr                                        # noqa: E402
from src.common import save_results, exp_seeds                      # noqa: E402
from scipy.stats import wilcoxon                                    # noqa: E402


def dataset(case, rng):
    if case == "S1_sincos":
        X = rng.uniform(-3, 3, (600, 1))
        f = lambda Z: np.sin(Z[:, 0]) * np.cos(Z[:, 0])
    elif case == "S2_rational":
        def sample(n):
            x = rng.uniform(0.0, 2.0, n)
            x = x[np.abs(x - 1.0) > 0.02][:n]
            while len(x) < n:
                more = rng.uniform(0.0, 2.0, n)
                more = more[np.abs(more - 1.0) > 0.02]
                x = np.concatenate([x, more])[:n]
            return x
        X = sample(600)[:, None]
        f = lambda Z: (Z[:, 0] ** 3 - 1) / (Z[:, 0] - 1)
    else:
        X = rng.uniform(1.05, 4.0, (600, 1))
        f = lambda Z: np.log(Z[:, 0] ** 2 - 1) - np.log(Z[:, 0] - 1)
    y = f(X)
    Xte = rng.uniform(X.min(), X.max(), (400, 1))
    return X, y, Xte, f(Xte)


def main():
    from gplearn.genetic import SymbolicTransformer  # noqa: E402
    out = {}
    for case in ("S1_sincos", "S2_rational", "S3_log"):
        ours_nmse, gp_nmse, gp_sizes = [], [], []
        for seed in exp_seeds(9)[:30]:
            rng = np.random.default_rng(seed + 555)
            Xtr, ytr, Xte, yte = dataset(case, rng)
            est = SymbolicTransformer(population_size=1000, generations=20,
                                      tournament_size=5, p_crossover=0.7,
                                      p_subtree_mutation=0.15, p_hoist_mutation=0.05,
                                      p_point_mutation=0.1, max_samples=0.9,
                                      function_set=("add", "sub", "mul", "div",
                                                    "sin", "cos", "log", "sqrt"),
                                      parsimony_coefficient=0.001,
                                      random_state=int(seed % 2**31), n_jobs=1)
            try:
                est.fit(Xtr, ytr)
                best = None
                for prog in est._programs[-1]:
                    if prog is None:
                        continue
                    p = prog.execute(Xte)
                    e = float(np.mean((p - yte) ** 2)) / max(np.var(yte), 1e-300)
                    if np.isfinite(e) and (best is None or e < best[0]):
                        best = (e, int(prog.length_))
                if best and np.isfinite(best[0]):
                    gp_nmse.append(min(best[0], 1e6))
                    gp_sizes.append(best[1])
                else:
                    gp_nmse.append(1e6); gp_sizes.append(-1)
            except Exception:
                gp_nmse.append(1e6); gp_sizes.append(-1)

            reg = sr.SymbolicRegressor(1, population=300, generations=80,
                                       seed=seed).fit(Xtr, ytr)
            tree, _ = min(reg.pareto_front_.items(),
                          key=lambda kv: len(ytr) * np.log(max(kv[1][1] * np.var(ytr),
                                                               1e-300)) + kv[0] * np.log(len(ytr)))[1]
            ours_nmse.append(sr.nmse(tree, Xte, yte))
        a = np.array(ours_nmse); b = np.array(gp_nmse)
        try:
            p = float(wilcoxon(a, b).pvalue)
        except ValueError:
            p = 1.0
        out[case] = {"ours_median_holdout_nmse": float(np.median(a)),
                     "gplearn_median_holdout_nmse": float(np.median(b)),
                     "gplearn_median_nodes": float(np.median(gp_sizes)),
                     "wilcoxon_p": p}
        print(f"{case}: ours med NMSE={np.median(a):.2e} | "
              f"gplearn med NMSE={np.median(b):.2e} "
              f"med nodes={np.median(gp_sizes):.0f} | Wilcoxon p={p:.3f}")

    save_results("exp09b_gplearn.json", out)


if __name__ == "__main__":
    main()
