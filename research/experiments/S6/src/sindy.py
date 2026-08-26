"""SINDy-style equation discovery: candidate-library construction + STLSQ + model selection.

The learner sees only (X, dX) arrays. Nothing in this module imports ground truth.
"""
import itertools
import numpy as np


# ------------------------------------------------------------------ library
def build_library_named(X, var_names, poly_degree=2, trig_vars=()):
    """Candidate feature dictionary: 1, monomials up to poly_degree, sin/cos of chosen vars.

    Returns (Theta (n,p), names list, theta_fn callable mapping (m,d)->(m,p)).
    """
    n, d = X.shape
    """Same as build_library but with explicit variable names used in term labels."""
    n, d = X.shape

    def term_fn(exponent):
        def fn(S):
            out = np.ones(len(S))
            for j, e in enumerate(exponent):
                if e:
                    out = out * S[:, j] ** e
            return out
        return fn

    terms, names = [], []
    terms.append(term_fn((0,) * d)); names.append("1")
    for deg in range(1, poly_degree + 1):
        for expo in itertools.combinations_with_replacement(range(d), deg):
            exponent = [0] * d
            for j in expo:
                exponent[j] += 1
            label = "*".join(var_names[j] + (f"^{exponent[j]}" if exponent[j] > 1 else "")
                             for j in sorted(set(expo)))
            terms.append(term_fn(tuple(exponent))); names.append(label)
    for j in trig_vars:
        for fname, f in (("sin", np.sin), ("cos", np.cos)):
            def g(S, f=f, j=j):
                return f(S[:, j])
            terms.append(g); names.append(f"{fname}({var_names[j]})")

    Theta = np.column_stack([t(X) for t in terms])
    theta_fn = lambda S: np.column_stack([t(np.atleast_2d(S)) for t in terms])
    return Theta, names, theta_fn


# ------------------------------------------------------------------ STLSQ
def stlsq(Theta, dX, threshold, max_iter=25):
    """Sequentially thresholded least squares, per target column.

    Columns are standardized internally; the final active set is refit on raw columns
    so reported coefficients are unbiased. Returns coefficient matrix (p, n_targets)
    with zeros for inactive terms.
    """
    Theta = np.asarray(Theta, float)
    dX = np.asarray(dX, float)
    p = Theta.shape[1]
    norms = np.linalg.norm(Theta, axis=0)
    norms[norms == 0] = 1.0
    Thn = Theta / norms

    C = np.zeros((p, dX.shape[1]))
    for k in range(dX.shape[1]):
        y = dX[:, k]
        active = np.ones(p, bool)
        c = np.linalg.lstsq(Thn, y, rcond=None)[0]
        for _ in range(max_iter):
            small = np.abs(c) < threshold
            if np.all(small):
                active[:] = False
                c = np.zeros(p)
                break
            if np.any(small & active):
                active &= ~small
            c_new = np.zeros(p)
            c_new[active] = np.linalg.lstsq(Thn[:, active], y, rcond=None)[0]
            if np.array_equal(c_new != 0, c != 0) and np.allclose(c_new, c):
                c = c_new
                break
            c = c_new
        idx = np.where(active)[0]
        if len(idx):                      # unbiased refit on original scale
            C[idx, k] = np.linalg.lstsq(Theta[:, idx], y, rcond=None)[0]
    return C


def bic(C, Theta, dX):
    """Bayesian information criterion of a (possibly multi-target) linear model."""
    resid = dX - Theta @ C
    n = dX.size
    rss = float(np.sum(resid**2))
    k = int(np.count_nonzero(C))
    if rss <= 0:
        return -np.inf
    return n * np.log(rss / n) + k * np.log(n)


def fit_sindy(Theta, dX, thresholds=None):
    """STLSQ over a threshold grid; select by BIC. Returns dict with C, threshold, bic."""
    if thresholds is None:
        ymax = float(np.std(dX))
        thresholds = np.geomspace(ymax * 1e-7, ymax * 10, 40)
    best = None
    for thr in thresholds:
        C = stlsq(Theta, dX, thr)
        score = bic(C, Theta, dX)
        if best is None or score < best["bic"]:
            best = {"C": C, "threshold": thr, "bic": score}
    return best


# ------------------------------------------------------------------ baselines
def ols_full(Theta, dX):
    return np.linalg.lstsq(Theta, dX, rcond=None)[0]


def ridge(Theta, dX, alpha=1e-3):
    from sklearn.linear_model import Ridge
    mu, sd = Theta.mean(0), Theta.std(0)
    sd[sd == 0] = 1.0
    Ths = (Theta - mu) / sd
    out = np.zeros_like(ols_full(Theta, dX))
    for k in range(dX.shape[1]):
        m = Ridge(alpha=alpha, fit_intercept=False).fit(Ths, dX[:, k])
        out[:, k] = m.coef_ / sd
    return out


def lasso_bic(Theta, dX):
    from sklearn.linear_model import Lasso
    mu, sd = Theta.mean(0), Theta.std(0)
    sd[sd == 0] = 1.0
    Ths = (Theta - mu) / sd
    alphas = np.geomspace(1e-5, 1.0, 20)
    best = None
    for a in alphas:
        C = np.zeros_like(ols_full(Theta, dX))
        ok = True
        for k in range(dX.shape[1]):
            try:
                m = Lasso(alpha=a, fit_intercept=False, max_iter=20000).fit(Ths, dX[:, k])
            except Exception:
                ok = False
                break
            C[:, k] = m.coef_ / sd
        if not ok:
            continue
        score = bic(C, Theta, dX)
        if best is None or score < best[0]:
            best = (score, C)
    return best[1]


# ------------------------------------------------------------------ simulation
class DiscoveredModel:
    """Wraps fitted coefficients into a simulable RHS."""

    def __init__(self, C, theta_fn, d):
        self.C, self.theta_fn, self.d = C, theta_fn, d

    def rhs(self, t, s):
        th = self.theta_fn(np.asarray(s)[None, :])[0]
        return self.C.T @ th

    def simulate(self, y0, t_eval):
        from scipy.integrate import solve_ivp
        sol = solve_ivp(self.rhs, (t_eval[0], t_eval[-1]), np.asarray(y0, float),
                        t_eval=t_eval, rtol=1e-10, atol=1e-12, method="RK45")
        return sol
