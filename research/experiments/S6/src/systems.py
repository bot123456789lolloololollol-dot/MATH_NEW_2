"""Ground-truth dynamical systems for S6 equation-discovery study.

Every simulator returns a dict with:
  X : (n, d) float64 state samples on a uniform grid
  t : (n,)   time grid
  names : list of state-variable names
  f_true : human-readable RHS strings (for documentation only, never used by learners)
"""
import numpy as np
from scipy.integrate import solve_ivp

RTOL, ATOL = 1e-10, 1e-12


def _integrate(f, t0, t1, y0, dt, events=None, dense=False):
    t_eval = np.arange(t0, t1 + 0.5 * dt, dt)
    sol = solve_ivp(f, (t0, t1), np.atleast_1d(y0), method="RK45",
                    rtol=RTOL, atol=ATOL, t_eval=t_eval, events=events,
                    dense_output=dense)
    if not sol.success:
        raise RuntimeError("integration failed: " + sol.message)
    return sol


# ---------------------------------------------------------------- SYS1
def damped_oscillator(k=2.0, c=0.3, x0=1.5, v0=0.0, t_end=10.0, dt=2e-3):
    def f(t, s):
        return [s[1], -k * s[0] - c * s[1]]
    sol = _integrate(f, 0.0, t_end, [x0, v0], dt)
    return {"X": sol.y.T, "t": sol.t, "names": ["x", "v"],
            "f_true": ["v", "-2.0 x - 0.3 v"], "f": lambda s: np.stack(
                [s[..., 1], -k * s[..., 0] - c * s[..., 1]], axis=-1)}


# ---------------------------------------------------------------- SYS2
def lotka_volterra(a=1.0, b=0.1, c=1.1, d=0.4, x0=(10.0, 5.0), t_end=25.0, dt=2e-3):
    def f(t, s):
        x, y = s
        return [a * x - b * x * y, d * x * y - c * y]
    sol = _integrate(f, 0.0, t_end, list(x0), dt)
    return {"X": sol.y.T, "t": sol.t, "names": ["x", "y"],
            "f_true": ["1.0 x - 0.1 x y", "-1.1 y + 0.4 x y"], "f": lambda s: np.stack(
                [a * s[..., 0] - b * s[..., 0] * s[..., 1],
                 d * s[..., 0] * s[..., 1] - c * s[..., 1]], axis=-1)}


# ---------------------------------------------------------------- SYS3
def lorenz(sigma=10.0, beta=8.0 / 3.0, rho=28.0, x0=(-8.0, 7.0, 27.0),
           t_end=20.0, dt=1e-3):
    def f(t, s):
        x, y, z = s
        return [sigma * (y - x), x * (rho - z) - y, x * y - beta * z]
    sol = _integrate(f, 0.0, t_end, list(x0), dt)
    return {"X": sol.y.T, "t": sol.t, "names": ["x", "y", "z"],
            "f_true": ["10 (y-x)", "28 x - y - x z", "x y - (8/3) z"]}


def lorenz_f(sigma=10.0, beta=8.0 / 3.0, rho=28.0):
    def f(s):
        x, y, z = s[..., 0], s[..., 1], s[..., 2]
        return np.stack([sigma * (y - x), x * (rho - z) - y, x * y - beta * z], axis=-1)
    return f


# ---------------------------------------------------------------- SYS4
def pendulum_period(theta0_deg, g_over_L=1.0, n_periods_max=40.0):
    """Measure dimensionless period tau = T*sqrt(g/L) of the nonlinear pendulum.

    Period measured from successive upward zero-crossings of theta (event detection).
    Returns tau; starts at rest at theta0.
    """
    th0 = np.deg2rad(theta0_deg)

    def f(t, s):
        return [s[1], -g_over_L * np.sin(s[0])]

    def cross_up(t, s):
        return s[0]
    cross_up.direction = 1
    cross_up.terminal = False

    # generous time span: linear period is 2pi; slow large-amplitude motion ~ 4-6x
    t_end = min(2 * np.pi / np.sqrt(g_over_L) * n_periods_max, 400.0)
    sol = _integrate(f, 0.0, t_end, [th0, 0.0], 5e-3, events=cross_up)
    te = sol.t_events[0]
    if len(te) < 3:
        raise RuntimeError(f"period not found for theta0={theta0_deg}")
    T = te[2] - te[0]
    return T * np.sqrt(g_over_L)


# ---------------------------------------------------------------- SYS5
def two_body(mu=1.0, a=1.0, e=0.3, n_orbits=3, dt=1e-2, drag=None):
    """Planar two-body problem in relative coordinates, periapsis on +x axis.

    drag(t, s): optional velocity-dependent perturbation used only in falsification tests.
    Returns trajectory and measured period (via periapsis passages).
    """
    r0 = a * (1.0 - e)
    v0 = np.sqrt(mu * (1.0 + e) / (a * (1.0 - e)))

    def f(t, s):
        x, y, vx, vy = s
        r = np.hypot(x, y)
        acc = [-mu * x / r**3, -mu * y / r**3]
        if drag is not None:
            acc[0] += drag(t, s)[0]
            acc[1] += drag(t, s)[1]
        return [vx, vy, acc[0], acc[1]]

    def periapsis(t, s):
        return s[1]
    periapsis.direction = 1

    t_end = 2 * np.pi * np.sqrt(a**3 / mu) * (n_orbits + 0.5)
    sol = _integrate(f, 0.0, t_end, [r0, 0.0, 0.0, v0], dt, events=periapsis)
    te = sol.t_events[0]
    # the start point (periapsis, y=0, upward) registers as an event at t~0,
    # so average consecutive gaps instead of differencing first-to-last
    T_meas = float(np.mean(np.diff(te))) if len(te) >= 3 else np.nan
    return {"X": sol.y.T, "t": sol.t, "names": ["x", "y", "vx", "vy"],
            "T": T_meas, "T_theory": 2 * np.pi * np.sqrt(a**3 / mu)}


# ---------------------------------------------------------------- SYS6
def rlc_response(R=10.0, L=0.05, C=10e-6, q0=1e-3, t_end=0.05, dt=1e-4):
    """Series RLC capacitor discharge: L q'' + R q' + q/C = 0. State [q, i]."""
    def f(t, s):
        q, i = s
        return [i, (-q / C - R * i) / L]
    sol = _integrate(f, 0.0, t_end, [q0, 0.0], dt)
    return {"X": sol.y.T, "t": sol.t, "names": ["q", "i"],
            "alpha_true": R / (2 * L),
            "omega_d_true": np.sqrt(1.0 / (L * C) - (R / (2 * L)) ** 2)}


# ---------------------------------------------------------------- SYS7/8
def duffing(delta=0.0, x0=1.8, v0=0.0, t_end=50.0, dt=5e-3, forcing=0.0, omega=1.0):
    """Duffing oscillator x'' = x - x^3 - delta x' + F cos(omega t).

    delta=0: conservative double-well, H = v^2/2 + x^4/4 - x^2/2.
    """
    def f(t, s):
        x, v = s
        return [v, x - x**3 - delta * v + forcing * np.cos(omega * t)]
    sol = _integrate(f, 0.0, t_end, [x0, v0], dt)
    out = {"X": sol.y.T, "t": sol.t, "names": ["x", "v"]}
    if delta == 0.0 and forcing == 0.0:
        H = 0.5 * out["X"][:, 1]**2 + 0.25 * out["X"][:, 0]**4 - 0.5 * out["X"][:, 0]**2
        out["H_true"] = H
    return out


def central_force_2d(mu=1.0, r0=(1.0, 0.3), v0=(0.2, 1.1), t_end=30.0, dt=1e-2):
    """Generic central-force orbit (not necessarily closed); angular momentum conserved."""
    def f(t, s):
        x, y, vx, vy = s
        r = np.hypot(x, y)
        return [vx, vy, -mu * x / r**3, -mu * y / r**3]
    sol = _integrate(f, 0.0, t_end, list(r0) + list(v0), dt)
    return {"X": sol.y.T, "t": sol.t, "names": ["x", "y", "vx", "vy"]}


def differentiate(X, dt, smooth=False, order=6):
    """Aligned (states, derivatives) pair for system identification.

    Returns X_sub, dX where dX[i] approximates d/dt of X_sub[i].
    order=6 uses a seven-point stencil (O(dt^6)); its leading truncation term
    is ~3 orders of magnitude below the O(dt^4) one at our sampling rates,
    which matters because truncation residue otherwise projects onto spurious
    library terms during sparse regression.
    """
    if order == 6:
        dX = (X[6:] - 9 * X[5:-1] + 45 * X[4:-2] - 45 * X[2:-4]
              + 9 * X[1:-5] - X[:-6]) / (60 * dt)
        # dX covers original indices 3..n-4; drop 3 more on each side
        return X[6:-6], dX[3:-3]
    dX = finite_differences(X, dt, smooth=smooth,
                            order=order, trim=(order == 4))
    if order == 4:
        return X[4:-4], dX
    return X, dX


def finite_differences(X, dt, smooth=False, order=4, trim=True):
    """Central-difference derivative along axis 0.

    order=2: np.gradient (O(dt^2)); order=4: five-point stencil (O(dt^4)).
    With order=4 the two edge rows on each side are dropped unless trim=False
    (they would need contaminated one-sided formulas).
    """
    if smooth:
        from scipy.signal import savgol_filter
        X = savgol_filter(X, window_length=11, polyorder=3, axis=0, mode="interp")
    if order == 2:
        return np.gradient(X, dt, axis=0)
    if order != 4:
        raise ValueError("order must be 2 or 4")
    dX = (-X[4:] + 8 * X[3:-1] - 8 * X[1:-3] + X[:-4]) / (12 * dt)
    if trim:
        return dX[2:-2]
    out = np.empty_like(X)
    out[2:-2] = dX
    out[:2] = np.gradient(X, dt, axis=0)[:2]
    out[-2:] = np.gradient(X, dt, axis=0)[-2:]
    return out
