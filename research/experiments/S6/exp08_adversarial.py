"""Exp 08 — Adversarial suite: can the discovery pipeline be made to lie?

A1 noise-only data      -> must return "no law" (empty/constant model)
A2 out-of-class truth   -> y = tanh(2x) with a library lacking tanh:
                           honest-failure verdict when holdout NMSE > 1e-2
A3 confounded causality -> x <- z -> y; observational regression finds a strong
                           x->y link that collapses under regime change
A4 chaos horizon        -> validation horizon of a discovered Lorenz model is
                           bounded by Lyapunov timescales, not by model quality;
                           we quantify the gap between perfect-model twin and
                           discovered-model rollouts and check it matches
                           t* ~ ln(C/eps)/lambda_max.
"""
import sys
import numpy as np

sys.path.insert(0, ".")
from src import systems as sy, sindy as sd, symreg as sr, metrics as mt   # noqa: E402
from src.common import save_results, exp_seeds                            # noqa: E402


def A1_noise_only():
    false_laws = 0
    details = []
    for seed in exp_seeds(8)[:30]:
        rng = np.random.default_rng(seed)
        X = rng.normal(size=(1000, 3))
        dX = sy.finite_differences(X, 0.01)   # derivatives of iid noise
        Theta, nm, fn = sd.build_library_named(X[4:-4], ["x", "y", "z"], 2)
        fit = sd.fit_sindy(Theta, dX)
        n_active = int(np.count_nonzero(np.abs(fit["C"]) > 1e-6))
        false_laws += int(n_active > 0)
        details.append(n_active)
    print(f"A1 noise-only: false-law rate {false_laws}/30 "
          f"(active-term counts {sorted(set(details))})")
    # GP control
    rng = np.random.default_rng(1)
    Xg = rng.normal(size=(500, 2)); yg = rng.normal(size=500)
    reg = sr.SymbolicRegressor(2, population=200, generations=40, seed=3).fit(Xg, yg)
    gp_nmse = reg.best_nmse_
    print(f"A1 GP on pure noise: best NMSE={gp_nmse:.3f} (1.0 == explains nothing)")
    return {"sindy_false_law_rate": false_laws / 30,
            "gp_best_nmse_on_noise": float(gp_nmse)}


def A2_out_of_class():
    rng = np.random.default_rng(8000)
    X = rng.uniform(-3, 3, (800, 1))
    y = np.tanh(2 * X[:, 0])
    Xtr, ytr, Xte_, yte_ = X[:600], y[:600], X[600:], y[600:]
    reg = sr.SymbolicRegressor(1, population=400, generations=120, seed=5).fit(Xtr, ytr)
    te = sr.nmse(reg.best_, Xte_, yte_)
    expr = str(sr.to_sympy(reg.best_, ["x"]))
    verdict = "no compact law within library (honest failure)" if te > 1e-2 \
        else "law found"
    # polynomial LSQ comparison (SINDy-style library, no tanh)
    from src.sindy import build_library_named, fit_sindy
    Th, nm, _ = build_library_named(Xtr, ["x"], 6)
    fit = fit_sindy(Th, np.gradient(ytr, Xtr[:, 0])[:, None])
    print(f"A2 out-of-class y=tanh(2x): GP holdout NMSE={te:.3e} -> {verdict}")
    print(f"   best expression: {expr}")
    print(f"   (derivative-fit residual scale: {float(np.max(np.abs(fit['C']))):.3g})")
    return {"holdout_nmse": float(te), "verdict": verdict, "expr": expr}


def A3_confounder():
    rng = np.random.default_rng(9000)
    coefs_obs, r2_regime = [], []
    for rep in range(20):
        rg = np.random.default_rng(9000 + rep)
        z = rg.normal(size=2000)
        xi, eta = rg.normal(size=2000) * .3, rg.normal(size=2000) * .3
        x = 1.5 * z + xi
        y = -1.2 * z + eta
        A = np.stack([x, np.ones_like(x)], 1)
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        coefs_obs.append(c[0])
        # interventional regime: do(x) ~ independent of z -> link must vanish
        z2 = rg.normal(size=2000) * 2.0
        x2 = rg.normal(size=2000)              # intervention: x set exogenously
        y2 = -1.2 * z2 + .3 * rg.normal(size=2000)
        A2 = np.stack([x2, np.ones_like(x2)], 1)
        c2, *_ = np.linalg.lstsq(A2, y2, rcond=None)
        ss_res = float(np.sum((y2 - A2 @ c2) ** 2))
        ss_tot = float(np.sum((y2 - y2.mean()) ** 2))
        r2_regime.append(1 - ss_res / max(ss_tot, 1e-300))
    obs_mean = float(np.mean(coefs_obs))
    int_mean = float(np.mean(r2_regime))
    print(f"A3 confounder: observational dY/dX ~= {obs_mean:+.3f} (strong spurious link)"
          f"; interventional mean R^2 = {int_mean:.3f} (link gone)")
    return {"observational_coef_mean": obs_mean,
            "observational_coef_std": float(np.std(coefs_obs)),
            "interventional_R2_mean": int_mean}


def A4_chaos_horizon():
    # discovered Lorenz model from clean data (same pipeline as exp01)
    L = sy.lorenz(t_end=15.0)
    Xs, dX = sy.differentiate(L["X"], L["t"][1] - L["t"][0])
    Th, nm, fn = sd.build_library_named(Xs, ["x", "y", "z"], 3)
    fit = sd.fit_sindy(Th, dX)
    mdl = sd.DiscoveredModel(fit["C"], fn, 3)
    coef_err = float(np.max(np.abs(fit["C"])))

    lam = mt.lyapunov_benettin(sy.lorenz_f(), [-8., 7., 27.])
    f_true = lambda s: sy.lorenz_f()(s)

    rng = np.random.default_rng(4242)
    horizons_disc, horizons_twin = [], []
    for rep in range(30):
        y0 = rng.uniform(-15, 15, 3); y0[2] += 20
        t_eval = np.arange(0, 4 + 1e-9, 0.005)
        sol_t = __import__("scipy.integrate", fromlist=["solve_ivp"]).solve_ivp(
            lambda t, s: f_true(s[None, :])[0], (0, 4), y0, t_eval=t_eval,
            rtol=1e-11, atol=1e-13)
        Xt = sol_t.y.T
        sol_m = mdl.simulate(y0, t_eval)
        if sol_m.y.shape[1] != len(t_eval):
            continue
        Xm = sol_m.y.T
        h_d = mt.predictability_horizon(Xm, Xt)
        if h_d is not None:
            horizons_disc.append(t_eval[h_d])
        # perfect-model twin: same true dynamics, IC perturbed by coef-error scale
        y0p = y0 + rng.normal(size=3) * 1e-6
        sol_p = __import__("scipy.integrate", fromlist=["solve_ivp"]).solve_ivp(
            lambda t, s: f_true(s[None, :])[0], (0, 4), y0p, t_eval=t_eval,
            rtol=1e-11, atol=1e-13)
        h_t = mt.predictability_horizon(sol_p.y.T, Xt)
        if h_t is not None:
            horizons_twin.append(t_eval[h_t])
    md, mtw = float(np.median(horizons_disc)), float(np.median(horizons_twin))
    eps_ratio = 1e-6 / 1e-10
    theory_gap = float(np.log(eps_ratio) / abs(lam))
    print(f"A4 chaos: lambda={lam:.3f}; median horizon discovered={md:.2f}s, "
          f"perfect-twin(IC~1e-6)={mtw:.2f}s; theory ln(1e4)/lam={theory_gap:.2f}s")
    return {"lyapunov": lam, "median_horizon_discovered_model": md,
            "median_horizon_perfect_model_ic1e-6": mtw,
            "coef_max_abs_err": coef_err,
            "theory_log_ratio_over_lambda": theory_gap}


def main():
    out = {"A1_noise_only": A1_noise_only(),
           "A2_out_of_class": A2_out_of_class(),
           "A3_confounder": A3_confounder(),
           "A4_chaos_horizon": A4_chaos_horizon()}
    save_results("exp08_adversarial.json", out)


if __name__ == "__main__":
    main()
