"""Exp 10 — PDE discovery from simulation data (round 2, preregistered R2).

Systems: viscous Burgers (nu=0.05) and KdV on periodic domains.
Learner: spatiotemporal library {1, u, u^2, u_x, u_xx, u_xxx, u*u_x, u^2*u_x,
u*u_xx} regressed onto u_t by STLSQ with significance support; time derivative from
snapshot finite differences in t; spatial derivatives spectral.
Validation: coefficient exactness + held-out IC rollout of the discovered PDE vs truth.
Noise levels (Burgers only): sigma in {0, 1e-4, 1e-3} * max|u| with spectral lowpass.
Prior art to cite: Rudy et al., Science 2017 (data-driven PDE discovery / PDE-FIND).
"""
import sys
import numpy as np

sys.path.insert(0, ".")
from src import pde as PD, sindy as sd                               # noqa: E402
from src.common import save_results                                  # noqa: E402


def build_dataset(sim, noise_sigma=0.0, seed=0):
    g = sim["grid"]
    U = sim["U"].copy()
    scale = float(np.max(np.abs(U)))
    if noise_sigma > 0:
        rng = np.random.default_rng(seed)
        U = U + rng.normal(size=U.shape) * (noise_sigma * scale)
        U = g.lowpass(U, frac=0.75)
    t = sim["T"]
    dt = t[1] - t[0]
    # 4th-order central difference in time:
    # f'(t_m) ~ (-f_{m+2} + 8 f_{m+1} - 8 f_{m-1} + f_{m-2}) / (12 dt)
    Ut = (-U[4:] + 8 * U[3:-1] - 8 * U[1:-3] + U[:-4]) / (12 * dt)
    Um = U[2:-2]
    n_t, n_x = Um.shape
    feats = {
        "1": np.ones((n_t, n_x)),
        "u": Um,
        "u^2": Um**2,
        "u_x": g.dx(Um),
        "u_xx": g.dx2(Um),
        "u_xxx": g.dx3(Um),
        "u*u_x": Um * g.dx(Um),
        "u^2*u_x": Um**2 * g.dx(Um),
        "u*u_xx": Um * g.dx2(Um),
    }
    names = list(feats)
    Theta = np.column_stack([feats[k].ravel() for k in names])
    dX = Ut.ravel()[:, None]
    return names, Theta, dX, g, sim


def rollout_burgers(C_coef, names, g, nu_true, t_end, n=128):
    """Integrate discovered RHS from the true-model IC; compare against true rollout.

    Both integrations use identical (relaxed) tolerances -- the metric compares two
    models, not against integration precision.
    """
    idx_map = {k: names.index(k) for k in names}
    coef = {k: float(C_coef[idx_map[k]]) for k in names if abs(C_coef[idx_map[k]]) > 0}

    def rhs_from(coefs):
        def f(t, uh):
            u = np.real(np.fft.ifft(uh))
            ux = g.dx(u)
            vals = {"1": 1.0, "u": u, "u^2": u * u,
                    "u_x": ux, "u_xx": g.dx2(u), "u_xxx": g.dx3(u),
                    "u*u_x": u * ux, "u^2*u_x": u * u * ux,
                    "u*u_xx": u * g.dx2(u)}
            r = 0.0
            for k, c in coefs.items():
                r = r + c * vals[k]
            return np.fft.fft(r)
        return f

    if nu_true is not None:
        def truth(t, uh):
            u = np.real(np.fft.ifft(uh))
            return np.fft.fft(-u * g.dx(u) + nu_true * g.dx2(u))
        x = g.x
        u0 = np.sin(x) * (1 + np.cos(x))
    else:
        def truth(t, uh):
            u = np.real(np.fft.ifft(uh))
            return np.fft.fft(-6 * u * g.dx(u) - g.dx3(u))
        x = g.x
        u0 = 12.0 / np.cosh((x - 8.0) / 2.0) ** 2

    t_eval = np.linspace(0, t_end, 41)
    ivp = __import__("scipy.integrate", fromlist=["solve_ivp"])
    solT = ivp.solve_ivp(truth, (0, t_end), np.fft.fft(u0), t_eval=t_eval,
                         method="RK45", rtol=1e-8, atol=1e-9)
    solM = ivp.solve_ivp(rhs_from(coef), (0, t_end), np.fft.fft(u0), t_eval=t_eval,
                         method="RK45", rtol=1e-8, atol=1e-9)
    XT = np.real(np.fft.ifft(solT.y, axis=0))
    XM = np.real(np.fft.ifft(solM.y, axis=0))
    return float(np.linalg.norm(XM - XT) / np.linalg.norm(XT))


def identify(names, Theta, dX):
    fit = sd.fit_sindy(Theta, dX)
    C = fit["C"][:, 0]
    act = sd.significant_support(Theta, dX, fit["C"], z=5.0)[:, 0]
    return C, act


def main():
    out = {}

    # ---------------- Burgers
    sim = PD.burgers_sim()
    for sigma in (0.0, 1e-4, 1e-3):
        names, Theta, dX, g, s = build_dataset(sim, sigma, seed=1010 + int(sigma * 1e6))
        C, act = identify(names, Theta, dX)
        active = {names[i]: round(float(C[i]), 8) for i in range(len(names)) if act[i]}
        rerr = rollout_burgers(C, names, g, nu_true=0.05, t_end=1.5)
        cerr = max(abs(C[names.index("u*u_x")] + 1), abs(C[names.index("u_xx")] - 0.05)) \
            if ("u*u_x" in names and "u_xx" in names) else None
        out[f"burgers_sigma{sigma}"] = {"active_terms": active,
                                        "coeff_err_active": float(cerr),
                                        "holdout_rel_err_new_IC": rerr}
        print(f"Burgers sigma={sigma}: terms={active} coeff_err={cerr:.2e} "
              f"holdout(new IC) rel err={rerr:.3e}")

    # ---------------- KdV (clean)
    simK = PD.kdv_sim()
    names, Theta, dX, g, _ = build_dataset(simK)
    C, act = identify(names, Theta, dX)
    active = {names[i]: round(float(C[i]), 6) for i in range(len(names)) if act[i]}
    # held-out new IC (single soliton, different from the two-soliton training data);
    # horizon per R2 deviation note
    rerr = rollout_burgers(C, names, g, nu_true=None, t_end=0.8, n=g.n)
    cerr = abs(C[names.index("u*u_x")] + 6) if "u*u_x" in names else None
    cerr = max(cerr, abs(C[names.index("u_xxx")] + 1))
    out["kdv_clean"] = {"active_terms": active, "coeff_err_active": float(cerr),
                        "holdout_rel_err_new_IC": rerr}
    print(f"KdV: terms={active} coeff_err={cerr:.2e} holdout(new IC) rel err={rerr:.3e}")

    save_results("exp10_pde.json", out)


if __name__ == "__main__":
    main()
