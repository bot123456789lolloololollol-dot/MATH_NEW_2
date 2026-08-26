"""Exp 02 — Noise robustness of equation discovery (preregistered sweep).

sigma in {0, 0.005, 0.01, 0.02, 0.05} x 30 seeds for Lotka-Volterra and Lorenz.
Metrics: coefficient error, support Jaccard (tolerance 5% of max |true coef|),
held-out rollout error. Savitzky-Golay smoothing + 4th-order derivatives when sigma>0.
"""
import sys
import numpy as np

sys.path.insert(0, ".")
from src import systems as sy, sindy as sd, metrics as mt          # noqa: E402
from src.common import save_results, exp_seeds                     # noqa: E402
from src.plotting import newfig, savefig                           # noqa: E402

SIGMAS = [0.0, 0.005, 0.01, 0.02, 0.05]


def one_run(system, sigma, seed):
    rng = np.random.default_rng(seed)
    if system == "lotka_volterra":
        d = sy.lotka_volterra(t_end=25.0)
        names, degree = ["x", "y"], 2
        true_rhs = lambda s: np.stack([s[..., 0] - .1 * s[..., 0] * s[..., 1],
                                       .4 * s[..., 0] * s[..., 1] - 1.1 * s[..., 1]], -1)
        active = ["x", "y", "x*y"]
        y0, horizon = [12.0, 4.0], 20.0
    else:
        d = sy.lorenz(t_end=15.0)
        names, degree = ["x", "y", "z"], 3
        true_rhs = lambda s: np.stack([10 * (s[..., 1] - s[..., 0]),
                                       28 * s[..., 0] - s[..., 1] - s[..., 0] * s[..., 2],
                                       s[..., 0] * s[..., 1] - (8 / 3) * s[..., 2]], -1)
        active = ["x", "y", "z", "x*z", "x*y"]
        y0, horizon = [6.0, -2.0, 18.0], 2.0

    X = d["X"]
    if sigma > 0:
        X = X + rng.normal(size=X.shape) * (sigma * X.std(axis=0))
    Xs, dX = sy.differentiate(X, d["t"][1] - d["t"][0], smooth=(sigma > 0))
    Theta, nm, theta_fn = sd.build_library_named(Xs, names, degree)
    fit = sd.fit_sindy(Theta, dX)
    C = fit["C"]

    scale = 1.1 if system == "lotka_volterra" else 28.0
    tol = 0.05 * scale
    if system == "lotka_volterra":
        tvals = {("x", 0): 1.0, ("x*y", 0): -0.1, ("y", 1): -1.1, ("x*y", 1): 0.4}
    else:
        tvals = {("y", 0): 10, ("x", 0): -10, ("x", 1): 28, ("z", 1): -1,
                 ("x*z", 1): -1, ("x*y", 2): 1, ("z", 2): -8 / 3}
    C_true = np.zeros_like(C)
    for (nmeq, k), ct in tvals.items():
        i = nm.index(nmeq)
        C_true[i, k] = ct
    cerr = mt.coeff_rel_err(C, C_true)
    at = np.zeros(len(nm), bool)
    for nmeq in active:
        at[nm.index(nmeq)] = True
    jac = mt.support_jaccard(C, at, tol)
    mdl = sd.DiscoveredModel(C, theta_fn, d["X"].shape[1])
    rerr = mt.rollout_error(mdl.rhs, true_rhs, y0,
                            np.arange(0.0, horizon + 1e-9, 0.01))
    return {"coeff_rel_err": cerr, "support_jaccard": jac, "rollout_rel_err": rerr}


def main():
    results = {}
    for system in ("lotka_volterra", "lorenz"):
        per_sigma = {}
        for sig in SIGMAS:
            runs = [one_run(system, sig, s) for s in exp_seeds(2)]
            agg = {k: (float(np.mean([r[k] for r in runs])),
                       float(np.std([r[k] for r in runs]))) for k in runs[0]}
            per_sigma[str(sig)] = {"mean_std": agg, "n": len(runs)}
            print(f"{system} sigma={sig}: coeff_err={agg['coeff_rel_err'][0]:.3e} "
                  f"jaccard={agg['support_jaccard'][0]:.3f} "
                  f"rollout={agg['rollout_rel_err'][0]:.3e}")
        results[system] = per_sigma

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for ax, system in zip(axes, ("lotka_volterra", "lorenz")):
        xs = SIGMAS
        for key, label in (("coeff_rel_err", "coefficient error"),
                           ("rollout_rel_err", "rollout error"),
                           ("support_jaccard", "support Jaccard")):
            ys = [results[system][str(s)]["mean_std"][key][0] for s in xs]
            es = [results[system][str(s)]["mean_std"][key][1] for s in xs]
            ax.errorbar(xs, ys, yerr=es, marker="o", capsize=3, label=label)
        ax.set_xscale("symlog", linthresh=1e-3)
        ax.set_yscale("log")
        ax.set_xlabel("noise sigma (fraction of state std)")
        ax.set_title(system)
        ax.legend(fontsize=8)
        ax.grid(alpha=.3)
    axes[0].set_ylabel("error (log)")
    savefig(fig, str(__import__("pathlib").Path("figures/exp02_noise.png")))

    # verdict: recovery preserved up to which noise level? (rollout<5%, jaccard=1)
    verdict = {}
    for system in results:
        ok_levels = []
        for s in SIGMAS:
            m = results[system][str(s)]["mean_std"]
            if m["support_jaccard"][0] >= 0.99 and m["rollout_rel_err"][0] < 0.05:
                ok_levels.append(s)
        verdict[system] = {"max_noise_with_recovery_mean": max(ok_levels) if ok_levels else None}
    save_results("exp02_noise.json", {"results": results, "verdict": verdict})
    print("verdict:", verdict)


if __name__ == "__main__":
    main()
