"""Unit tests for the S6 equation-discovery stack. Run: pytest tests/ -q"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import systems as sy, sindy as sd, symreg as sr, invariant as inv, metrics as mt  # noqa


def test_differentiate_sixth_order():
    t = np.arange(0, 8, 1e-3)
    X = np.sin(3 * t)[:, None]
    Xs, dX = sy.differentiate(X, 1e-3)
    assert np.max(np.abs(dX[:, 0] - 3 * np.cos(3 * t[6:-6]))) < 1e-10


def test_stlsq_recovers_oscillator_exactly():
    d = sy.damped_oscillator(t_end=8.0)
    Xs, dX = sy.differentiate(d["X"], d["t"][1] - d["t"][0])
    Th, names, fn = sd.build_library_named(Xs, ["x", "v"], 2)
    fit = sd.fit_sindy(Th, dX)
    C = fit["C"]
    assert abs(C[names.index("v"), 0] - 1) < 1e-8
    assert abs(C[names.index("x"), 1] + 2) < 1e-8
    assert abs(C[names.index("v"), 1] + 0.3) < 1e-8
    # support exactly {x, v}
    act = {names[i] for i in range(len(names)) if np.abs(C[i]).max() > 1e-9}
    assert act == {"x", "v"}


def test_discovered_model_rollout():
    d = sy.damped_oscillator(t_end=8.0)
    Xs, dX = sy.differentiate(d["X"], d["t"][1] - d["t"][0])
    Th, names, fn = sd.build_library_named(Xs, ["x", "v"], 2)
    fit = sd.fit_sindy(Th, dX)
    mdl = sd.DiscoveredModel(fit["C"], fn, 2)
    rerr = mt.rollout_error(mdl.rhs,
                            lambda s: np.stack([s[..., 1], -2 * s[..., 0] - .3 * s[..., 1]], -1),
                            [0.4, 1.0], np.linspace(0, 5, 100))
    assert rerr < 1e-6


def test_protected_operations_never_nan():
    rng = np.random.default_rng(0)
    X = rng.uniform(-5, 5, (500, 1))
    cases = {
        "div": ("div", ("var", 0), ("add", ("var", 0), ("const", 0.0))),
        "log": ("log", ("var", 0)),
        "sqrt": ("sqrt", ("var", 0)),
        "exp": ("exp", ("var", 0)),
    }
    for tree in cases.values():
        v = sr.evaluate(tree, X)
        assert np.all(np.isfinite(v)), tree[0]


def test_gp_finds_minimal_equivalent_form():
    rng = np.random.default_rng(0)
    X = rng.uniform(-3, 3, (400, 1))
    y = np.sin(X[:, 0]) * np.cos(X[:, 0])
    reg = sr.SymbolicRegressor(1, population=300, generations=100, seed=42).fit(X, y)
    assert reg.best_nmse_ < 1e-10
    assert sr.nodes(reg.best_) <= 9


def test_pendulum_period_matches_elliptic_integral():
    from scipy.special import ellipk
    for a_deg in (10.0, 30.0):
        tau_sys = sy.pendulum_period(a_deg)
        tau_ref = 4 * ellipk(np.sin(np.deg2rad(a_deg) / 2) ** 2)
        assert abs(tau_sys / tau_ref - 1) < 1e-7


def test_two_body_periods():
    for a, e, mu in ((1.0, 0.2, 1.0), (1.4, 0.45, 1.6)):
        o = sy.two_body(mu=mu, a=a, e=e, n_orbits=3)
        assert abs(o["T"] / o["T_theory"] - 1) < 1e-8


def test_invariant_recovery_and_negative_control():
    trajs = [sy.duffing(x0=1.8)["X"], sy.duffing(x0=-0.4, v0=0.9)["X"]]
    feat = lambda S: inv.poly_features(S, degree=4, var_names=["x", "v"])[0]
    res = inv.discover_invariant(trajs, feat)
    assert res["sigma_ratio"] < 1e-8 and res["spectral_gap"] > 100
    _, pnames = inv.poly_features(trajs[0], degree=4, var_names=["x", "v"])
    hdir = np.zeros(len(pnames))
    hdir[pnames.index("v^2")] = .5
    hdir[pnames.index("x^2")] = -.5
    hdir[pnames.index("x^4")] = .25
    cosang = abs(res["c"] @ hdir) / (np.linalg.norm(res["c"]) * np.linalg.norm(hdir))
    assert cosang > 1 - 1e-9
    held = [sy.duffing(x0=0.7, v0=-1.1)["X"], sy.duffing(x0=2.2)["X"]]
    assert inv.conservation_ratio(res["F"], held) < 1e-7
    # negative control: damping destroys the nullspace
    neg = inv.discover_invariant([sy.duffing(delta=0.3)["X"]], feat)
    assert not (neg["sigma_ratio"] < 1e-8 and neg["spectral_gap"] > 100)


def test_rlc_laws_recovered():
    from exp06_rlc import build_dataset
    rows = build_dataset()
    R, L, C, a_t, w_t, a_m, w_m = rows.T
    ok = np.isfinite(a_m) & np.isfinite(w_m)
    assert np.median(np.abs(a_m[ok] - a_t[ok]) / a_t[ok]) < 1e-2
    Xa = np.stack([R[ok] / L[ok], np.ones(ok.sum())], 1)
    ca, *_ = np.linalg.lstsq(Xa, a_m[ok], rcond=None)
    assert abs(ca[0] - 0.5) < 1e-2


def test_lyapunov_lorenz_near_literature():
    lam = mt.lyapunov_benettin(sy.lorenz_f(), [-8., 7., 27.], n_steps=4000)
    assert 0.6 < lam < 1.2


def test_generalization_gate_rejects_noise_laws():
    rng = np.random.default_rng(3)
    X = rng.normal(size=(1000, 3))
    Xs, dX = sy.differentiate(X, 0.01)
    Th, names, fn = sd.build_library_named(Xs, ["x", "y", "z"], 2)
    fit = sd.fit_sindy(Th[:700], dX[:700])
    tail = dX[700:] - Th[700:] @ fit["C"]
    nmse_tail = float(np.sum(tail**2) / np.sum((dX[700:] - dX[700:].mean(0))**2))
    assert nmse_tail > 0.9   # nothing explains unseen noise
