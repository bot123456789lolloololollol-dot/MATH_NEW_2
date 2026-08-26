"""Shared metrics for the S6 study. All definitions are in PREREGISTERED.md."""
import numpy as np
from scipy.integrate import solve_ivp


def coeff_rel_err(C_hat, C_true):
    """Max relative coefficient error over truly-active entries.

    C_hat, C_true: (n_terms, n_targets). Compares entrywise wherever
    C_true != 0; a coefficient placed on the wrong term therefore scores ~1.
    """
    mask = np.abs(C_true) > 0
    if not np.any(mask):
        return 0.0
    scale = float(np.max(np.abs(C_true[mask])))
    return float(np.max(np.abs(C_hat[mask] - C_true[mask])) / scale)


def support_jaccard(C_hat, active_true, tol):
    act = np.abs(C_hat).max(axis=1) > tol if C_hat.ndim > 1 else np.abs(C_hat) > tol
    inter = np.sum(act & active_true)
    union = np.sum(act | active_true)
    return float(inter / union) if union else 1.0


def rollout_error(f_model, f_true_vecf, y0, t_eval):
    """Relative Frobenius error between model rollout and true rollout from same IC."""
    sol_m = solve_ivp(lambda t, s: f_model(t, s), (t_eval[0], t_eval[-1]),
                      np.asarray(y0, float), t_eval=t_eval, rtol=1e-10, atol=1e-12)
    X = np.stack([f_true_vecf(np.asarray(y0)[None, :])[0]], 0)
    # integrate true system too
    sol_t = solve_ivp(lambda t, s: np.asarray(f_true_vecf(s[None, :])[0]),
                      (t_eval[0], t_eval[-1]), np.asarray(y0, float),
                      t_eval=t_eval, rtol=1e-10, atol=1e-12)
    if not (sol_m.success and sol_t.success):
        return np.inf
    Xm, Xt = sol_m.y.T, sol_t.y.T
    return float(np.linalg.norm(Xm - Xt) / max(np.linalg.norm(Xt), 1e-300))


def predictability_horizon(Xm, Xt, thresh=0.25):
    """First time index where normalized deviation exceeds thresh; None if never."""
    dev = np.linalg.norm(Xm - Xt, axis=1) / np.maximum(np.linalg.norm(Xt, axis=1), 1e-12)
    idx = np.where(dev > thresh)[0]
    return None if len(idx) == 0 else int(idx[0])


def lyapunov_benettin(f, y0, dt=0.01, n_steps=6000, delta0=1e-7, renorm_every=10,
                      seed=0):
    """Largest Lyapunov exponent via Benettin et al.; f maps state->dstate/dt."""
    rng = np.random.default_rng(seed)

    def rhs(t, s):
        return np.asarray(f(s[None, :])[0])
    y = np.asarray(y0, float).copy()
    d = np.asarray(y0, float) + delta0 * rng.normal(size=len(y0))
    lam, n_rounds = 0.0, n_steps // renorm_every
    from scipy.integrate import solve_ivp
    for _ in range(n_rounds):
        sY = solve_ivp(rhs, (0, dt * renorm_every), y, rtol=1e-11, atol=1e-13)
        sD = solve_ivp(rhs, (0, dt * renorm_every), d, rtol=1e-11, atol=1e-13)
        y = sY.y[:, -1]
        diff = sD.y[:, -1] - y
        dist = np.linalg.norm(diff)
        lam += np.log(dist / delta0)
        d = y + diff * (delta0 / dist)
    return float(lam / (n_rounds * dt * renorm_every))


def nmse(pred, y):
    return float(np.mean((np.asarray(pred) - np.asarray(y)) ** 2) / max(np.var(y), 1e-300))
