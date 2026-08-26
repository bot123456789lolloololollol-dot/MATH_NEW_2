"""Exp 01 — Can the system recover known governing equations exactly from clean data?

Systems: damped linear oscillator, Lotka-Volterra, Lorenz-63.
Learner sees only (states, finite-difference derivatives). Verdict per preregistered
criterion: coefficient error < 1e-6, support Jaccard = 1, held-out rollout < 5%.
"""
import sys
import numpy as np

sys.path.insert(0, ".")
from src import systems as sy, sindy as sd, metrics as mt          # noqa: E402
from src.common import save_results                                # noqa: E402


def identify(X, dt, var_names, poly_degree, noise_sigma=0.0, seed=0):
    if noise_sigma > 0:
        rng = np.random.default_rng(seed)
        scale = noise_sigma * X.std(axis=0)
        Xn = X + rng.normal(size=X.shape) * scale
        Xs, dX = sy.differentiate(Xn, dt, smooth=True)
    else:
        Xs, dX = sy.differentiate(X, dt)
    Theta, names, theta_fn = sd.build_library_named(Xs, var_names, poly_degree)
    fit = sd.fit_sindy(Theta, dX)
    return fit, names, sd.DiscoveredModel(fit["C"], theta_fn, X.shape[1])


def report(name, fit, names, eq_names, C_true, active_true, model, true_rhs, y0,
           t_horizon, dt_out=0.01):
    C = fit["C"]
    jaccard = mt.support_jaccard(C, active_true, 1e-8)
    cerr = mt.coeff_rel_err(C, C_true)
    t_eval = np.arange(0.0, t_horizon + 1e-9, dt_out)
    rerr = mt.rollout_error(model.rhs, true_rhs, y0, t_eval)
    terms = {eq: {names[i]: round(C[i, k], 10) for i in range(len(names))
                  if abs(C[i, k]) > 1e-8}
             for k, eq in enumerate(eq_names)}
    ok = (cerr < 1e-6) and (jaccard == 1.0) and (rerr < 0.05)
    rec = {"system": name, "terms": terms, "coeff_rel_err": cerr,
           "support_jaccard": jaccard, "rollout_rel_err": rerr,
           "bic": fit["bic"], "threshold": fit["threshold"],
           "recovered_per_preregistered_criterion": bool(ok)}
    print(f"--- {name}: coeff_err={cerr:.2e} jaccard={jaccard} "
          f"rollout={rerr:.2e} -> {'RECOVERED' if ok else 'NOT RECOVERED'}")
    for eq, tt in terms.items():
        print("   ", eq, "=", "  ".join(f"{c:+g}*{t}" for t, c in tt.items()))
    return rec


def main():
    out = []

    # SYS1 damped linear oscillator: x' = v ; v' = -2 x - 0.3 v
    d = sy.damped_oscillator(t_end=10.0)
    fit, names, mdl = identify(d["X"], d["t"][1] - d["t"][0], ["x", "v"], 2)
    Ct = np.zeros((len(names), 2)); at = np.zeros(len(names), bool)
    Ct[names.index("v"), 0] = 1.0
    Ct[names.index("x"), 1] = -2.0; Ct[names.index("v"), 1] = -0.3
    at[names.index("v")] = True; at[names.index("x")] = True
    out.append(report("damped_oscillator", fit, names, ["x_dot", "v_dot"], Ct, at, mdl,
                      lambda s: np.stack([s[..., 1], -2 * s[..., 0] - .3 * s[..., 1]], -1),
                      [0.3, -1.2], 10.0))

    # SYS2 Lotka-Volterra
    lv = sy.lotka_volterra(t_end=25.0)
    fit, names, mdl = identify(lv["X"], lv["t"][1] - lv["t"][0], ["x", "y"], 2)
    Ct = np.zeros((len(names), 2)); at = np.zeros(len(names), bool)
    Ct[names.index("x"), 0] = 1.0; Ct[names.index("x*y"), 0] = -0.1
    Ct[names.index("y"), 1] = -1.1; Ct[names.index("x*y"), 1] = 0.4
    for nm_ in ("x", "y", "x*y"):
        at[names.index(nm_)] = True
    out.append(report("lotka_volterra", fit, names, ["x_dot", "y_dot"], Ct, at, mdl,
                      lambda s: np.stack([s[..., 0] - .1 * s[..., 0] * s[..., 1],
                                          .4 * s[..., 0] * s[..., 1] - 1.1 * s[..., 1]], -1),
                      [12.0, 4.0], 20.0))

    # SYS3 Lorenz
    L = sy.lorenz(t_end=15.0)
    fit, names, mdl = identify(L["X"], L["t"][1] - L["t"][0], ["x", "y", "z"], 3)
    Ct = np.zeros((len(names), 3)); at = np.zeros(len(names), bool)
    Ct[names.index("y"), 0] = 10; Ct[names.index("x"), 0] = -10
    Ct[names.index("x"), 1] = 28; Ct[names.index("z"), 1] = -1; Ct[names.index("x*z"), 1] = -1
    Ct[names.index("x*y"), 2] = 1; Ct[names.index("z"), 2] = -8.0 / 3.0
    for nm in ("x", "y", "z", "x*z", "x*y"):
        at[names.index(nm)] = True
    f_true = lambda s: np.stack([10 * (s[..., 1] - s[..., 0]),
                                 28 * s[..., 0] - s[..., 1] - s[..., 0] * s[..., 2],
                                 s[..., 0] * s[..., 1] - (8 / 3) * s[..., 2]], -1)
    out.append(report("lorenz63", fit, names, ["x_dot", "y_dot", "z_dot"], Ct, at,
                      mdl, f_true, [6.0, -2.0, 18.0], 2.0))

    save_results("exp01_known_laws.json", {"runs": out})


if __name__ == "__main__":
    main()
