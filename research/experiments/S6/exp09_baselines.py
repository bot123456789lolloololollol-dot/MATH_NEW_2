"""Exp 09 — Baseline comparison (preregistered): does sparsity-aware discovery beat
standard regressors, and does MDL-selected GP beat accuracy-maximizing GP?

Dynamics task: Lotka-Volterra at sigma=0.01 and 0.02, 30 paired seeds.
  Methods: STLSQ+BIC (ours), OLS-full-library, Ridge(alpha=1e-3), Lasso(BIC grid).
  Metrics: coefficient error, support Jaccard, rollout error; paired Wilcoxon.
Function task: y=sin(x)cos(x)+noise(1%), 30 seeds:
  GP-MDL selection vs GP-min-train-error vs degree-15 polynomial OLS.
Seeds and protocol fixed in PREREGISTERED.md before any run.
"""
import sys
import numpy as np
from scipy.stats import wilcoxon

sys.path.insert(0, ".")
from src import systems as sy, sindy as sd, symreg as sr, metrics as mt   # noqa: E402
from src.common import save_results, exp_seeds                            # noqa: E402


def dynamics_one(sigma, seed):
    rng = np.random.default_rng(seed)
    d = sy.lotka_volterra(t_end=25.0)
    Xn = d["X"] + rng.normal(size=d["X"].shape) * (sigma * d["X"].std(axis=0))
    Xs, dX = sy.differentiate(Xn, d["t"][1] - d["t"][0], smooth=True)
    Theta, nm, fn = sd.build_library_named(Xs, ["x", "y"], 2)
    methods = {
        "stlsq_bic": sd.fit_sindy(Theta, dX)["C"],
        "ols_full": sd.ols_full(Theta, dX),
        "ridge": sd.ridge(Theta, dX, alpha=1e-3),
        "lasso_bic": sd.lasso_bic(Theta, dX),
    }
    scale = 1.1
    C_true = np.zeros((len(nm), 2))
    for nmeq in ("x", "y", "x*y"):
        pass
    true_entries = {("x", 0): 1.0, ("x*y", 0): -0.1, ("y", 1): -1.1, ("x*y", 1): 0.4}
    for (nmeq, k), ct in true_entries.items():
        C_true[nm.index(nmeq), k] = ct
    at = np.zeros(len(nm), bool)
    for nmeq in ("x", "y", "x*y"):
        at[nm.index(nmeq)] = True
    true_rhs = lambda s: np.stack([s[..., 0] - .1 * s[..., 0] * s[..., 1],
                                   .4 * s[..., 0] * s[..., 1] - 1.1 * s[..., 1]], -1)
    out = {}
    t_eval = np.arange(0.0, 20.0 + 1e-9, 0.05)
    for mname, C in methods.items():
        cerr = mt.coeff_rel_err(C, C_true)
        act_hat = sd.significant_support(Theta, dX, C, z=5.0)
        jac = float(np.sum(act_hat.any(1) & at) / np.sum(act_hat.any(1) | at))
        rerr = mt.rollout_error(sd.DiscoveredModel(C, fn, 2).rhs, true_rhs,
                                [12.0, 4.0], t_eval) if np.all(np.isfinite(C)) else 1e9
        out[mname] = {"coeff_rel_err": min(cerr, 1e6), "support_jaccard": jac,
                      "rollout_rel_err": min(rerr, 1e9)}
    return out


def function_one(seed):
    rng = np.random.default_rng(seed)
    Xtr = rng.uniform(-3, 3, (600, 1)); Xte = rng.uniform(-3, 3, (400, 1))
    f = lambda Z: np.sin(Z[:, 0]) * np.cos(Z[:, 0])
    ytr = f(Xtr) + rng.normal(size=600) * 0.01 * np.std(f(Xtr))
    yte = f(Xte)

    reg = sr.SymbolicRegressor(1, population=300, generations=80, seed=seed % 7).fit(
        Xtr, ytr)
    best_mdl, best_acc = None, None
    for k, (tree, e_tr) in reg.pareto_front_.items():
        cost = len(ytr) * np.log(max(e_tr * np.var(ytr), 1e-300)) + k * np.log(len(ytr))
        if best_mdl is None or cost < best_mdl[0]:
            best_mdl = (cost, tree)
        if best_acc is None or e_tr < best_acc[0]:
            best_acc = (e_tr, tree)
    gp_mdl_nmse = sr.nmse(best_mdl[1], Xte, yte)
    gp_acc_nmse = sr.nmse(best_acc[1], Xte, yte)

    # degree-15 polynomial OLS baseline
    P_tr = np.vander(Xtr[:, 0], 16, increasing=True)
    P_te = np.vander(Xte[:, 0], 16, increasing=True)
    c, *_ = np.linalg.lstsq(P_tr, ytr, rcond=None)
    poly_nmse = float(np.mean((P_te @ c - yte) ** 2) / np.var(yte))
    return {"gp_mdl": gp_mdl_nmse, "gp_min_train": gp_acc_nmse,
            "poly15_ols": poly_nmse}


def main():
    dyn = {}
    for sigma in (0.01, 0.02):
        per_method = {}
        runs = [dynamics_one(sigma, s) for s in exp_seeds(9)]
        for m in ("stlsq_bic", "ols_full", "ridge", "lasso_bic"):
            per_method[m] = {k: (float(np.mean([r[m][k] for r in runs])),
                                 float(np.std([r[m][k] for r in runs])))
                             for k in ("coeff_rel_err", "support_jaccard",
                                       "rollout_rel_err")}
        # paired significance: stlsq vs others on coeff err and rollout
        sig_tests = {}
        for other in ("ols_full", "ridge", "lasso_bic"):
            for metric in ("coeff_rel_err", "rollout_rel_err"):
                a = [r["stlsq_bic"][metric] for r in runs]
                b = [r[other][metric] for r in runs]
                try:
                    p = float(wilcoxon(a, b).pvalue)
                except ValueError:
                    p = 1.0
                sig_tests[f"stlsq_vs_{other}_{metric}"] = p
        dyn[str(sigma)] = {"methods": per_method, "wilcoxon_p": sig_tests}
        print(f"sigma={sigma}:")
        for m, agg in per_method.items():
            print(f"   {m:10s}: coeff={agg['coeff_rel_err'][0]:.3e} "
                  f"jac={agg['support_jaccard'][0]:.3f} "
                  f"rollout={agg['rollout_rel_err'][0]:.3e}")

    fun_runs = [function_one(s) for s in exp_seeds(9)]
    fun = {k: (float(np.mean([r[k] for r in fun_runs])),
               float(np.std([r[k] for r in fun_runs]))) for k in fun_runs[0]}
    print("function task holdout NMSE:", fun)
    try:
        p_gp = float(wilcoxon([r["gp_mdl"] for r in fun_runs],
                              [r["gp_min_train"] for r in fun_runs]).pvalue)
    except ValueError:
        p_gp = 1.0
    save_results("exp09_baselines.json",
                 {"dynamics": dyn, "function": fun, "wilcoxon_gp_vs_acc": p_gp})


if __name__ == "__main__":
    main()
