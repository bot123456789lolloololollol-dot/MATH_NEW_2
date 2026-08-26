"""Machine-checkable certificate for P4b (THEORY.md):

The invariant F discovered from Duffing trajectories alone, together with the vector
field f_hat discovered from the same trajectories by sparse regression, satisfies

    dF/dt = grad F . f_hat  ==  0   as a symbolic polynomial identity.

Discovered floats are snapped to rationals (nsimplify); the snapped forms are verified
to still match the data before certification, so the certificate is about what was
actually learned, not about a hand-written model.

Run:  python proofs_certificate_duffing.py     -> prints CERTIFIED / FAILED
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "experiments" / "S6"))
from src import systems as sy, sindy as sd, invariant as inv  # noqa: E402


def term_expr(label):
    e = sp.Integer(1)
    for factor in label.split("*"):
        if "^" in factor:
            base, p = factor.split("^")
            e = e * sp.Symbol(base) ** int(p)
        else:
            e = e * sp.Symbol(factor)
    return e


def snap(coef, tol=1e-6):
    r = sp.nsimplify(float(coef), rational=True, tolerance=tol)
    return r


def main():
    # --- discover dynamics by sparse regression (data only)
    d = sy.duffing(x0=1.8)
    Xs, dX = sy.differentiate(d["X"], d["t"][1] - d["t"][0])
    Th, names, fn = sd.build_library_named(Xs, ["x", "v"], poly_degree=4)
    fit = sd.fit_sindy(Th, dX)
    C = fit["C"]

    fx = fv = sp.Integer(0)
    for i, nm in enumerate(names):
        if abs(C[i, 0]) > 1e-9:
            fx += snap(C[i, 0]) * term_expr(nm)
        if abs(C[i, 1]) > 1e-9:
            fv += snap(C[i, 1]) * term_expr(nm)

    # snapped field must still reproduce the data (guard against over-snapping)
    err_x = np.max(np.abs(np.array([
        float(fx.subs({sp.Symbol("x"): s[0], sp.Symbol("v"): s[1]})) for s in Xs[:500]])
        - dX[:500, 0])) / np.max(np.abs(dX[:500, 0]))
    err_v = np.max(np.abs(np.array([
        float(fv.subs({sp.Symbol("x"): s[0], sp.Symbol("v"): s[1]})) for s in Xs[:500]])
        - dX[:500, 1])) / np.max(np.abs(dX[:500, 1]))
    print("discovered x' =", fx, " (snap rel err %.1e)" % err_x)
    print("discovered v' =", fv, " (snap rel err %.1e)" % err_v)

    # --- discover invariant from trajectory increments only
    trajs = [sy.duffing(x0=1.8)["X"], sy.duffing(x0=-0.4, v0=0.9)["X"]]
    res = inv.discover_invariant(
        trajs, lambda S: inv.poly_features(S, degree=4, var_names=["x", "v"])[0])
    _, pnames = inv.poly_features(trajs[0], degree=4, var_names=["x", "v"])
    c = res["c"]
    c_norm = c / np.max(np.abs(c))          # scale-invariant ratios
    F = sp.Integer(0)
    for ci, nm in zip(c_norm, pnames):
        if abs(ci) < 1e-9:
            continue
        F += snap(ci) * term_expr(nm)
    print("discovered invariant (up to scale) F =", F)

    # --- certificate: dF/dt along the discovered field vanishes identically
    dFdt = sp.expand(sp.diff(F, sp.Symbol("x")) * fx + sp.diff(F, sp.Symbol("v")) * fv)
    print("dF/dt expands to:", dFdt)
    ok_dyn = err_x < 1e-4 and err_v < 1e-4
    print("CERTIFIED" if (dFdt == 0 and ok_dyn) else "FAILED")
    return 0 if (dFdt == 0 and ok_dyn) else 1


if __name__ == "__main__":
    sys.exit(main())
