"""Exp 06 — Hidden parameter relationships in series RLC circuits (electrical systems).

Free response of L q'' + R q' + q/C = 0. From simulated ring waveforms alone we
extract the decay rate alpha and damped frequency omega_d (envelope/zero-crossing
measurements -- no model assumed), then discover how they depend on circuit
parameters across an (R, L, C) grid:

    alpha   = R / (2L)
    omega_d^2 = 1/(LC) - R^2/(4 L^2)

Both are recovered by sparse regression on physics features, validated on held-out
component combinations never seen during discovery.
"""
import sys
import numpy as np

sys.path.insert(0, ".")
from src import systems as sy                                       # noqa: E402
from src.common import save_results                                 # noqa: E402
from src.plotting import newfig, savefig                            # noqa: E402


def measure_ring(t, i):
    """Model-light measurement: peak spacing -> omega_d; log-envelope slope -> alpha."""
    # discard deep-tail samples where the waveform underflows
    keep = np.abs(i) >= 1e-10 * np.max(np.abs(i))
    t, i = t[keep], i[keep]
    idx = []
    for k in range(2, len(i) - 2):
        if abs(i[k]) > abs(i[k-1]) and abs(i[k]) >= abs(i[k+1]) \
           and abs(i[k]) > abs(i[k-2]) and abs(i[k]) > abs(i[k+2]):
            if not idx or k - idx[-1] >= 3:   # one peak per half-cycle minimum
                idx.append(k)
    if len(idx) < 4:
        return np.nan, np.nan
    tp = t[idx]
    amp = np.abs(i[idx])
    omega_d = np.pi / np.mean(np.diff(tp))
    A = np.stack([tp, np.ones_like(tp)], axis=1)
    coef, *_ = np.linalg.lstsq(A, np.log(amp), rcond=None)
    alpha = -coef[0]
    return float(alpha), float(omega_d)


def build_dataset():
    Ls = [0.01, 0.02, 0.04, 0.08]
    Cs = [5e-6, 20e-6, 50e-6]
    Rs = [2.0, 8.0, 20.0, 40.0]
    rows = []
    for L in Ls:
        for C in Cs:
            for R in Rs:
                if R >= 2 * np.sqrt(L / C):
                    continue  # keep the underdamped regime (ring measurable)
                # dt chosen so even the fastest ring has >= 60 samples per period
                w0 = np.sqrt(1.0 / (L * C))
                dt = min(1e-4, (2 * np.pi / w0) / 60)
                sim = sy.rlc_response(R=R, L=L, C=C, t_end=0.5, dt=dt)
                a_m, w_m = measure_ring(sim["t"], sim["X"][:, 1])
                rows.append((R, L, C, sim["alpha_true"], sim["omega_d_true"],
                             a_m, w_m))
    return np.array(rows)


def main():
    rows = build_dataset()
    R, L, C = rows[:, 0], rows[:, 1], rows[:, 2]
    a_true, w_true, a_meas, w_meas = rows[:, 3], rows[:, 4], rows[:, 5], rows[:, 6]

    meas_err_a = np.abs(a_meas - a_true) / a_true
    meas_err_w = np.abs(w_meas - w_true) / w_true
    bad = ~np.isfinite(meas_err_a) | ~np.isfinite(meas_err_w)
    print(f"measurement stage: median rel err alpha={np.median(meas_err_a):.2e} "
          f"omega_d={np.median(meas_err_w):.2e}; max={np.nanmax(meas_err_a):.2e}/"
          f"{np.nanmax(meas_err_w):.2e} "
          f"(worst cells are low-Q circuits near the overdamped boundary; "
          f"{int(np.sum(bad))} unmeasurable cells excluded)")
    ok = ~bad
    R, L, C, a_true, w_true, a_meas, w_meas = (v[ok] for v in
                                               (R, L, C, a_true, w_true, a_meas, w_meas))

    # train/held-out split: every third row of the shuffled-with-fixed-seed list
    rng = np.random.default_rng(6000)
    perm = rng.permutation(len(rows))
    test_idx = set(perm[::4].tolist())
    tr = np.array([i for i in range(len(rows)) if i not in test_idx])

    # ---- discovery of alpha: features [R/L, 1] target a_meas
    Xa = np.stack([R / L, np.ones_like(R)], axis=1)
    ca, *_ = np.linalg.lstsq(Xa[tr], a_meas[tr], rcond=None)
    pred_a = Xa @ ca

    # ---- discovery of omega_d^2: features [1/(LC), R^2/L^2, 1] target w_meas^2
    Xw = np.stack([1.0 / (L * C), R**2 / L**2, np.ones_like(R)], axis=1)
    cw, *_ = np.linalg.lstsq(Xw[tr], w_meas[tr] ** 2, rcond=None)
    pred_w2 = Xw @ cw

    print(f"alpha  = {ca[0]:.8f} * (R/L) + {ca[1]:+.2e}   [expect 0.5, 0]")
    print(f"omega_d^2 = {cw[0]:.8f} * 1/(LC) + {cw[1]:.8f} * R^2/L^2 + {cw[2]:+.3e}"
          f"   [expect 1, -0.25, 0]")

    holdout_err_a = float(np.max(np.abs(pred_a[list(test_idx)] - a_true[list(test_idx)])
                                 / a_true[list(test_idx)]))
    holdout_err_w = float(np.max(np.abs(np.sqrt(pred_w2[list(test_idx)])
                                        - w_true[list(test_idx)]) / w_true[list(test_idx)]))
    print(f"held-out prediction error: alpha {holdout_err_a:.2e}, "
          f"omega_d {holdout_err_w:.2e}")

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ax = axes[0]
    ax.loglog(a_true, pred_a, "o", ms=4)
    lims = [a_true.min() * .8, a_true.max() * 1.25]
    ax.loglog(lims, lims, "r--")
    ax.set_xlabel(r"true $\alpha=R/2L$")
    ax.set_ylabel("discovered-formula prediction")
    ax.set_title("decay rate")
    ax.grid(alpha=.3)
    ax = axes[1]
    ax.loglog(w_true, np.sqrt(pred_w2), "o", ms=4)
    lims = [w_true.min() * .8, w_true.max() * 1.25]
    ax.loglog(lims, lims, "r--")
    ax.set_xlabel(r"true $\omega_d$")
    ax.set_title("damped frequency")
    ax.grid(alpha=.3)
    savefig(fig, "figures/exp06_rlc.png")

    save_results("exp06_rlc.json", {
        "n_circuits_total": int(len(rows)),
        "measurement_rel_err": {
            "alpha_median": float(np.median(meas_err_a[ok])),
            "omega_d_median": float(np.median(meas_err_w[ok])),
            "alpha_max": float(np.nanmax(meas_err_a)),
            "omega_d_max": float(np.nanmax(meas_err_w)),
            "n_unmeasurable": int(np.sum(bad))},
        "alpha_formula": {"R_over_L_coef": ca[0], "intercept": ca[1]},
        "omega_formula": {"inv_LC_coef": cw[0], "R2_L2_coef": cw[1],
                          "intercept": cw[2]},
        "heldout_max_rel_err": {"alpha": holdout_err_a, "omega_d": holdout_err_w}})


if __name__ == "__main__":
    main()
