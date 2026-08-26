"""Exp 12 — hidden modal structure in coupled oscillators (round 3, preregistered R3).

Two coupled linear oscillators; the learner sees only state trajectories. After
recovering the system matrix A by linear sparse regression, the EIGENSTRUCTURE of the
discovered A (normal-mode frequencies and mode shapes) is compared with the analytic
eigenstructure -- structure that appears nowhere as explicit terms in the data.
"""
import sys
import numpy as np

sys.path.insert(0, ".")
from src import systems as sy, sindy as sd                          # noqa: E402
from src.common import save_results                                 # noqa: E402

K1, K2, KC = 4.0, 9.0, 1.5


def coupled_rhs(s):
    q1, q2, v1, v2 = s
    return [v1, v2,
            -K1 * q1 - KC * (q1 - q2),
            -K2 * q2 - KC * (q2 - q1)]


def simulate(y0, t_end=30.0, dt=5e-3, sigma=0.0, seed=0):
    from scipy.integrate import solve_ivp
    t_eval = np.arange(0.0, t_end + 1e-9, dt)
    sol = solve_ivp(lambda t, s: coupled_rhs(s), (0, t_end), list(y0),
                    t_eval=t_eval, rtol=1e-11, atol=1e-13, method="DOP853")
    X = sol.y.T
    if sigma > 0:
        rng = np.random.default_rng(seed)
        X = X + rng.normal(size=X.shape) * (sigma * X.std(axis=0))
    return X, t_eval


def analytic_A():
    A = np.zeros((4, 4))
    A[0, 2] = A[1, 3] = 1.0
    A[2, 0] = -(K1 + KC); A[2, 1] = KC
    A[3, 1] = -(K2 + KC); A[3, 0] = KC
    return A


def spectrum(A):
    """Sorted (ascending) normal-mode frequencies omega = |Im(eig)|."""
    ev = np.linalg.eigvals(A)
    omegas = np.sort(np.abs(ev.imag))
    return ev, omegas


def mode_ratio(A, omega_target):
    """Displacement ratio |q2/q1| of the eigenvector whose frequency is closest."""
    evs, V = np.linalg.eig(A)
    idx = int(np.argmin(np.abs(np.abs(evs.imag) - omega_target)))
    w = V[:, idx]
    if abs(w[0]) < 1e-15:
        return np.inf
    return float(abs(w[1] / w[0]))


def main():
    out = {}
    A_true = analytic_A()
    ev_t, om_t = spectrum(A_true)

    for label, sigma, seed in (("clean", 0.0, 12000), ("sigma0.005", 0.005, 12001),
                               ("sigma0.01", 0.01, 12002)):
        T1, te = simulate([1.0, 0.3, 0.0, 0.0], sigma=sigma, seed=seed)
        T2, _ = simulate([0.0, 1.0, 0.7, -0.2], sigma=sigma, seed=seed + 100)
        dt = te[1] - te[0]
        Xs = np.vstack([sy.differentiate(T1, dt)[0], sy.differentiate(T2, dt)[0]])
        dX = np.vstack([sy.differentiate(T1, dt)[1], sy.differentiate(T2, dt)[1]])
        names = ["q1", "q2", "v1", "v2"]
        Theta, nm, fn = sd.build_library_named(Xs, names, poly_degree=1)
        fit = sd.fit_sindy(Theta, dX)                 # constant col present; zeroed by STLSQ
        order = [nm.index(n) for n in names]
        A_hat = fit["C"][order, :].T                  # x' = A_hat x

        entry_err = float(np.max(np.abs(A_hat - A_true)))
        _, om_h = spectrum(A_hat)
        freq_err = float(np.max(np.abs(om_h - om_t) / om_t))
        ratios_t = [mode_ratio(A_true, o) for o in om_t]
        ratios_h = [mode_ratio(A_hat, o) for o in om_t]
        shape_err = float(max(abs(a / b - 1) for a, b in zip(ratios_h, ratios_t)))

        from scipy.integrate import solve_ivp
        Xho, te_ho = simulate([0.4, -0.6, 0.5, 0.2], t_end=20.0)
        sol_m = solve_ivp(lambda t, s: A_hat @ s, (te_ho[0], te_ho[-1]), Xho[0],
                          t_eval=te_ho, rtol=1e-10, atol=1e-12)
        rel = float(np.linalg.norm(sol_m.y.T - Xho) / max(np.linalg.norm(Xho), 1e-300))

        rec = {"max_abs_entry_err": entry_err, "max_rel_eig_freq_err": freq_err,
               "mode_shape_ratio_err": shape_err, "holdout_traj_rel_err": rel}
        out[label] = rec
        print(f"{label}: |A-A*|_max={entry_err:.2e} eigfreq_err={freq_err:.2e} "
              f"mode-shape err={shape_err:.2e} holdout={rel:.2e}")

    save_results("exp12_modes.json", out)


if __name__ == "__main__":
    main()
