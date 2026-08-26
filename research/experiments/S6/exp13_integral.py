"""Exp 13 — integral-form identification vs the derivative-form sigma=5% failure
(round 3, preregistered R3).

Motivation: exp02 showed derivative-based STLSQ fails on Lotka-Volterra at
sigma=0.05 (rollouts diverge). Integrating the equation over sliding windows turns
the regression target into state DIFFERENCES and the features into time-integrated
library values (cumulative sums), so noisy states are never differentiated.

Protocol: LV, sigma in {0.02, 0.05} x 30 paired seeds; same library, same
validation-split selection and z=5 significance; primary endpoint = fraction of runs
with held-out rollout error < 5% (preregistered recovery criterion). A negative
result is acceptable and will be reported as such.
"""
import sys
import numpy as np

sys.path.insert(0, ".")
from src import systems as sy, sindy as sd                          # noqa: E402
from src.common import save_results, exp_seeds                      # noqa: E402

TRUE_RHS = lambda s: np.stack([s[..., 0] - .1 * s[..., 0] * s[..., 1],
                               .4 * s[..., 0] * s[..., 1] - 1.1 * s[..., 1]], -1)
ACTIVE = ("x", "y", "x*y")
TRUE_ENTRIES = {("x", 0): 1.0, ("x*y", 0): -0.1, ("y", 1): -1.1, ("x*y", 1): 0.4}
Y0, HORIZON = [12.0, 4.0], 20.0


def make_data(sigma, seed):
    rng = np.random.default_rng(seed)
    d = sy.lotka_volterra(t_end=25.0)
    Xn = d["X"] + rng.normal(size=d["X"].shape) * (sigma * d["X"].std(axis=0))
    return Xn, d["t"][1] - d["t"][0]


def derivative_fit(Xn, dt):
    Xs, dX = sy.differentiate(Xn, dt, smooth=True)
    Theta, nm, fn = sd.build_library_named(Xs, ["x", "y"], 2)
    fit = sd.fit_sindy(Theta, dX)
    return fit, Theta, dX, nm, fn


def integral_fit(Xn, dt, W=25):
    """X(t+W)-X(t) = (sum of Theta rows)*dt @ c ; sliding windows."""
    from src.sindy import build_library_named
    # smooth lightly then integrate
    from scipy.signal import savgol_filter
    Xs = savgol_filter(Xn, 11, 3, axis=0, mode="interp")
    Theta, nm, fn = build_library_named(Xs, ["x", "y"], 2)
    cum = np.vstack([np.zeros((1, Theta.shape[1])),
                     np.cumsum(Theta, axis=0)]) * dt          # len n+1
    n = len(Xs)
    lo = np.arange(0, n - W)
    hi = lo + W
    # integral of Theta over [t_lo, t_hi] ~ sum_{k=lo}^{hi-1} Theta_k * dt
    Thi = cum[hi] - cum[lo]
    target = Xs[hi] - Xs[lo]
    fit = sd.fit_sindy(Thi, target)
    return fit, Thi, target, nm, fn


def metrics(fit, Theta_eval, dX_eval, nm, fn):
    C = fit["C"]
    C_true = np.zeros_like(C)
    for (nmeq, k), ct in TRUE_ENTRIES.items():
        C_true[nm.index(nmeq), k] = ct
    cerr = float(np.max(np.abs(C[C_true != 0] - C_true[C_true != 0])) / 1.1)
    act_hat = sd.significant_support(Theta_eval, dX_eval, C, z=5.0).any(1)
    at = np.zeros(len(nm), bool)
    for nmeq in ACTIVE:
        at[nm.index(nmeq)] = True
    jac = float(np.sum(act_hat & at) / np.sum(act_hat | at))
    mdl = sd.DiscoveredModel(C, fn, 2)
    rerr = mt_rollout(mdl)
    return {"coeff_rel_err": cerr, "support_jaccard": jac, "rollout_rel_err": rerr}


def mt_rollout(mdl):
    from src import metrics as mt
    rerr = mt.rollout_error(mdl.rhs, TRUE_RHS, Y0,
                            np.arange(0.0, HORIZON + 1e-9, 0.05))
    if not np.isfinite(rerr):
        rerr = 1e6
    return min(rerr, 1e6)


def main():
    out = {}
    for sigma in (0.02, 0.05):
        rows = []
        for seed in exp_seeds(13)[:30]:
            Xn, dt = make_data(sigma, seed)

            fit_d, Th_tr, dX_tr, nm, fn = derivative_fit(Xn, dt)
            m_d = metrics(fit_d, Th_tr, dX_tr, nm, fn)

            fit_i, Thi, tgt, nm2, fn2 = integral_fit(Xn, dt)
            # evaluate significance on the integrated system itself
            m_i = metrics(fit_i, Thi, tgt, nm2, fn2)

            rows.append({"deriv": m_d, "integral": m_i})
        agg = {}
        for meth in ("deriv", "integral"):
            rr = [r[meth] for r in rows]
            agg[meth] = {
                "recovery_rate_rollout<5%": float(np.mean(
                    [r["rollout_rel_err"] < 0.05 for r in rr])),
                "median_coeff_err": float(np.median([r["coeff_rel_err"] for r in rr])),
                "median_jaccard": float(np.median([r["support_jaccard"] for r in rr])),
                "median_rollout": float(np.median([r["rollout_rel_err"] for r in rr]))}
        out[str(sigma)] = agg
        print(f"sigma={sigma}: derivative recovery {agg['deriv']['recovery_rate_rollout<5%']:.2f}, "
              f"integral recovery {agg['integral']['recovery_rate_rollout<5%']:.2f} | "
              f"median rollout {agg['deriv']['median_rollout']:.3f} vs "
              f"{agg['integral']['median_rollout']:.3f}")

    save_results("exp13_integral.json", out)


if __name__ == "__main__":
    main()
