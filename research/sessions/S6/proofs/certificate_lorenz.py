"""Machine-checkable certificate #2 (chaotic case): the Lorenz equations discovered
from clean trajectory data satisfy their own claimed coefficient set exactly, and the
discovered vector field reproduces the data derivatives to integrator precision.

Complements certificate_duffing.py (conservation identity). Here the checkable claim is
exactness of recovery: after rational snapping, the discovered RHS equals the classical
Lorenz system symbolically.

Run:  python proofs_certificate_lorenz.py   -> prints CERTIFIED / FAILED
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "experiments" / "S6"))
from src import systems as sy, sindy as sd  # noqa: E402


def term_expr(label):
    e = sp.Integer(1)
    for factor in label.split("*"):
        if "^" in factor:
            base, p = factor.split("^")
            e = e * sp.Symbol(base) ** int(p)
        else:
            e = e * sp.Symbol(factor)
    return e


def snap(c, tol=1e-6):
    return sp.nsimplify(float(c), rational=True, tolerance=tol)


def main():
    L = sy.lorenz(t_end=10.0)
    Xs, dX = sy.differentiate(L["X"], L["t"][1] - L["t"][0])
    Th, names, fn = sd.build_library_named(Xs, ["x", "y", "z"], poly_degree=3)
    fit = sd.fit_sindy(Th, dX)

    exprs = []
    for k in range(3):
        e = sp.Integer(0)
        for i, nm in enumerate(names):
            if abs(fit["C"][i, k]) > 1e-9:
                e += snap(fit["C"][i, k]) * term_expr(nm)
        exprs.append(sp.expand(e))
    x, y, z = sp.Symbol("x"), sp.Symbol("y"), sp.Symbol("z")

    print("discovered:")
    for nm_, e in zip(("x'", "y'", "z'"), exprs):
        print("  ", nm_, "=", e)

    classical = [10 * (y - x), 28 * x - y - x * z, x * y - sp.Rational(8, 3) * z]
    identical = all(sp.simplify(e - c) == 0 for e, c in zip(exprs, classical))

    # snapped field must still reproduce the data (guard against over-snapping)
    n_chk = 300
    max_rel = 0.0
    for k, e in enumerate(exprs):
        f_num = sp.lambdify((x, y, z), e, "numpy")
        pred = f_num(Xs[:n_chk, 0], Xs[:n_chk, 1], Xs[:n_chk, 2])
        true = dX[:n_chk, k]
        max_rel = max(max_rel, float(np.max(np.abs(pred - true))
                                     / np.max(np.abs(true))))

    ok = identical and max_rel < 1e-4
    print("symbolically identical to classical Lorenz:", identical,
          "| max snap-vs-data rel err: %.1e" % max_rel)
    print("CERTIFIED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
