"""Round-2 unit tests: PDE machinery + control identification."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import pde as PD, sindy as sd, systems as sy          # noqa: E402


def test_time_stencil_sign():
    """The 4th-order time stencil in exp10 must differentiate sin(3t) correctly."""
    t = np.linspace(0, 2, 201)
    h = t[1] - t[0]
    U = np.sin(3 * t)[:, None]
    Ut = (-U[4:] + 8 * U[3:-1] - 8 * U[1:-3] + U[:-4]) / (12 * h)
    assert np.max(np.abs(Ut[:, 0] - 3 * np.cos(3 * t[2:-2]))) < 1e-6


def test_spectral_derivatives():
    g = PD.PeriodicGrid(n=64, L=2 * np.pi)
    x = g.x
    u = np.sin(2 * x)
    assert np.max(np.abs(g.dx(u) - 2 * np.cos(2 * x))) < 1e-10
    assert np.max(np.abs(g.dx2(u) + 4 * np.sin(2 * x))) < 1e-9


def test_burgers_pde_recovery():
    sim = PD.burgers_sim(t_end=1.5)
    from exp10_pde import build_dataset, identify
    names, Theta, dX, g, _ = build_dataset(sim)
    C, act = identify(names, Theta, dX)
    assert abs(C[names.index("u*u_x")] + 1) < 5e-3
    assert abs(C[names.index("u_xx")] - 0.05) < 5e-4
    got = {names[i] for i in range(len(names)) if act[i]}
    assert got == {"u*u_x", "u_xx"}


def test_control_identifier_finds_taylor_structure():
    from exp11_control import identify_model
    fn, C, terms = identify_model()
    # Taylor of -sin(theta) ~ -theta + theta^3/6 = -1.0 th + 0.1667 th^3
    assert abs(terms["th"] + 1) < 0.05
    assert abs(terms["om"] + 0.15) < 0.02
    assert abs(terms["th^3"] - 1 / 6) < 0.03
