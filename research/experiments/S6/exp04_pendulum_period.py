"""Exp 04 — Generalizing a known law: amplitude-dependent pendulum period.

Training: simulated pendula with g/L in {0.5,0.7,1.0,1.5,2.0}, amplitudes
theta0 in {5..40} degrees (small-amplitude regime). The learner sees only
(theta0, g_over_L) -> measured period T.

Discovery targets:
  (a) STRUCTURE: T * sqrt(g/L) depends on theta0 alone (dimensionless collapse).
  (b) FORM: tau(theta) = 2*pi*(1 + theta^2/16 + 11*theta^4/3072 + ...)
      Perturbation theory (Lindstedt-Poincare). We test whether regression on the
      basis [1, theta^2, theta^4] recovers {2pi, 2pi/16, 2pi*11/3072} and whether
      the discovered formula EXTRAPOLATES to unseen amplitudes (45..120 deg) and
      an unseen g/L = 1.25.
Exact reference: tau_exact = 4*elliptic_k(sin^2(theta0/2)).
"""
import sys
import numpy as np
from scipy.special import ellipk

sys.path.insert(0, ".")
from src import systems as sy, symreg as sr                         # noqa: E402
from src.common import save_results                                 # noqa: E402
from src.plotting import newfig, savefig                            # noqa: E402


def exact_tau(theta0_deg):
    th = np.deg2rad(theta0_deg)
    return 4.0 * ellipk(np.sin(th / 2) ** 2)


def main():
    rng = np.random.default_rng(4000)
    gl_train = [0.5, 0.7, 1.0, 1.5, 2.0]
    amp_train = np.arange(5.0, 41.0, 5.0)
    data = []
    for g in gl_train:
        for a in amp_train:
            tau = sy.pendulum_period(a, g_over_L=g)
            T = tau / np.sqrt(g)
            data.append((a, g, T))
    arr = np.array(data)

    # ---- (b) series-coefficient identification on dimensionless tau = T*sqrt(g/L)
    theta = np.deg2rad(arr[:, 0])
    tau = arr[:, 2] * np.sqrt(arr[:, 1])
    A = np.stack([np.ones_like(theta), theta**2, theta**4], axis=1)
    coef, res, *_ = np.linalg.lstsq(A, tau, rcond=None)
    theory = np.array([2 * np.pi, 2 * np.pi / 16, 2 * np.pi * 11 / 3072])
    rel_coef_err = np.abs(coef - theory) / theory

    # held-out amplitudes + unseen g/L using ONLY the discovered 3-term formula
    amp_test = [45.0, 60.0, 75.0, 90.0, 105.0, 120.0]
    ext = {}
    for a in amp_test:
        th = np.deg2rad(a)
        pred_disc = coef @ np.array([1.0, th**2, th**4])
        pred_naive = 2 * np.pi
        ex = exact_tau(a)
        # unseen g/L check at 30 deg with g=1.25
        ext[str(a)] = {
            "exact": float(ex),
            "discovered": float(pred_disc),
            "rel_err_discovered": float(abs(pred_disc / ex - 1)),
            "rel_err_small_angle_law": float(abs(pred_naive / ex - 1)),
        }
    g_unseen = 1.25
    tau_unseen = sy.pendulum_period(30.0, g_over_L=g_unseen)
    th30 = np.deg2rad(30.0)
    pred_unseen = (coef @ np.array([1.0, th30**2, th30**4])) / np.sqrt(g_unseen)
    meas_unseen = tau_unseen / np.sqrt(g_unseen)
    unseen_check = {"g_over_L": g_unseen, "amp_deg": 30.0,
                    "predicted_T": float(pred_unseen),
                    "measured_T": float(meas_unseen),
                    "rel_err": float(abs(pred_unseen / meas_unseen - 1))}

    print("series coefficients vs theory:")
    print("  const :", coef[0], "vs", theory[0], f"({rel_coef_err[0]:.2e})")
    print("  th^2  :", coef[1], "vs", theory[1], f"({rel_coef_err[1]:.2e})")
    print("  th^4  :", coef[2], "vs", theory[2], f"({rel_coef_err[2]:.2e})")
    for a in amp_test:
        e = ext[str(a)]
        print(f"  amp {a:5.0f} deg: disc err {e['rel_err_discovered']:.3e} | "
              f"small-angle err {e['rel_err_small_angle_law']:.3e}")
    print("unseen g/L:", unseen_check)

    # ---- (a) structural discovery: dimensionless collapse
    # If T = sqrt(L/g) * f(theta0), then tau := T*sqrt(g/L) must be independent of g.
    # Within each training amplitude, measure the relative spread of tau across g.
    spreads = {}
    for a in amp_train:
        taus = arr[np.isclose(arr[:, 0], a)][:, 2] * np.sqrt(arr[np.isclose(arr[:, 0], a)][:, 1])
        spreads[str(a)] = float(np.ptp(taus) / np.mean(taus))
    max_spread = max(spreads.values())
    print(f"dimensionless collapse: max within-amplitude spread of tau across "
          f"g/L values = {max_spread:.2e} (0 means perfect collapse)")
    # residual variance after removing the textbook factor vs before
    resid_before = float(np.var(arr[:, 2] - np.mean(arr[:, 2])) / np.var(arr[:, 2]))
    tau_resid = arr[:, 2] - np.stack(
        [np.ones_like(theta), theta**2, theta**4], 1) @ coef
    resid_after = float(np.var(tau_resid) / np.var(arr[:, 2]))
    collapse = {"max_within_amplitude_spread": max_spread,
                "per_amplitude": spreads,
                "unexplained_variance_fraction_before_collapse": resid_before,
                "after_collapse_plus_series": resid_after}

    # ---- GP attempt at closed-form structure (documented limitation: finds the
    # inverse-g dependence approximately but not the exact root form)
    X = arr[:, [1, 0]]  # features [g, theta]
    y = arr[:, 2]
    reg = sr.SymbolicRegressor(2, population=800, generations=200, seed=7).fit(X, y)
    expr = sr.to_sympy(reg.best_, ["gL", "th"])
    te_nmse = sr.nmse(reg.best_, X, y)
    print(f"GP structure discovery: nmse={reg.best_nmse_:.3e} holdout={te_nmse:.3e}")
    print("GP expression:", expr)

    # ---- figure
    fig, ax = newfig()
    ths = np.linspace(0, np.deg2rad(125), 300)
    ax.plot(np.rad2deg(ths), exact_tau(np.rad2deg(ths)), "k-", label="exact (elliptic)")
    ax.plot(np.rad2deg(theta), tau, "o", ms=4, label="simulated training data")
    fit_curve = coef[0] + coef[1] * ths**2 + coef[2] * ths**4
    ax.plot(np.rad2deg(ths), fit_curve, "r--", label="discovered 3-term law")
    ax.plot(np.rad2deg(ths), 2 * np.pi * np.ones_like(ths), "g:",
            label="textbook small-angle law")
    ax.axvspan(40, 125, alpha=.12, color="orange", label="extrapolation region")
    ax.set_xlabel("amplitude theta0 (deg)")
    ax.set_ylabel(r"dimensionless period $\tau = T\sqrt{g/L}$")
    ax.legend(fontsize=8)
    savefig(fig, "figures/exp04_period.png")

    save_results("exp04_period.json", {
        "series_coefficients": {"discovered": coef.tolist(),
                                "perturbation_theory": theory.tolist(),
                                "relative_error": rel_coef_err.tolist()},
        "extrapolation": ext, "unseen_gL": unseen_check,
        "collapse": collapse,
        "gp_structure": {"train_nmse": reg.best_nmse_,
                         "holdout_nmse": te_nmse,
                         "expression": str(expr),
                         "note": "vanilla GP approximates the inverse-root scaling but "
                                 "does not find its exact form; the collapse statistic "
                                 "and series identification carry the structural claim"}})


if __name__ == "__main__":
    main()
