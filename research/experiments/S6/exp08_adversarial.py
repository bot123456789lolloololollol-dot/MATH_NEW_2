"""Exp 08 — Adversarial suite: can the discovery pipeline be made to lie?

A1 noise-only data      -> any "law" must die at the preregistered generalization
                           gate (held-out one-step NMSE > 1e-2 on untouched rows).
                           Two datasets: iid noise, and a smooth OU-integrated curve.
A2 out-of-class truth   -> y = tanh(2x) with a library lacking tanh. Verdict rule:
                           "law found" iff holdout NMSE < 1e-2 AND compact (<=25 nodes,
                           post-hoc addition -- see PREREGISTERED deviations).
A3 confounded causality -> x <- z -> y; observational regression finds a strong x->y
                           link that collapses under intervention on x.
A4 chaos horizon        -> validation horizon of a discovered Lorenz model is bounded
                           by Lyapunov timescales: t*(eps) ~ ln(C*delta/eps)/lambda.
                           We measure horizons for the discovered model and for perfect
                           twins with IC perturbations 1e-6 / 1e-10 and check the
                           spacing matches theory.
"""
import sys
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, ".")
from src import systems as sy, sindy as sd, symreg as sr, metrics as mt   # noqa: E402
from src.common import save_results, exp_seeds                            # noqa: E402


def _fit_with_gate(Theta, dX, split):
    """Fit on rows[:split]; gate on untouched rows[split:]. Returns C, gate_nmse."""
    fit = sd.fit_sindy(Theta[:split], dX[:split])
    C = fit["C"]
    tail_d = dX[split:] - Theta[split:] @ C
    nmse_tail = float(np.sum(tail_d**2) / max(np.sum((dX[split:]
                      - dX[split:].mean(axis=0))**2), 1e-300))
    return C, nmse_tail


def A1_noise_only():
    raw_support_runs, gated_false = [], 0
    for seed in exp_seeds(8)[:30]:
        rng = np.random.default_rng(seed)
        # dataset 1: iid noise states
        X = rng.normal(size=(1000, 3))
        Xs, dX = sy.differentiate(X, 0.01)
        Theta, nm, fn = sd.build_library_named(Xs, ["x", "y", "z"], 2)
        C, gate = _fit_with_gate(Theta, dX, 700)
        n_active = int(np.count_nonzero(np.abs(C) > 1e-6))
        claimed = n_active > 0 and gate < 1e-2
        gated_false += int(claimed)
        raw_support_runs.append(n_active)
    print(f"A1 iid-noise: runs where STLSQ returns nonzero support: "
          f"{sum(r > 0 for r in raw_support_runs)}/30 "
          f"(counts {sorted(set(raw_support_runs))}); after generalization gate: "
          f"false laws = {gated_false}/30")

    # smooth-process control: integrated Ornstein-Uhlenbeck (physical-looking curve)
    rng = np.random.default_rng(777)
    dt = 0.01
    T = np.arange(0, 10 + 1e-9, dt)
    Xs_list = []
    for ch in range(3):
        ou = np.zeros(len(T))
        for k in range(1, len(T)):
            ou[k] = ou[k-1] + (-ou[k-1] * 0.5) * dt + rng.normal() * np.sqrt(dt) * .5
        Xs_list.append(ou)
    Xsm = np.stack(Xs_list, 1)
    Xa, dXa = sy.differentiate(Xsm, dt)
    Th2, _, _ = sd.build_library_named(Xa, ["x", "y", "z"], 2)
    C2, gate2 = _fit_with_gate(Th2, dXa, 700)
    n2 = int(np.count_nonzero(np.abs(C2) > 1e-6))
    print(f"A1 smooth-OU control: support={n2}, held-out gate NMSE={gate2:.3f} "
          f"-> claimed law: {bool(n2 > 0 and gate2 < 1e-2)}")

    reg = sr.SymbolicRegressor(2, population=200, generations=40, seed=3).fit(
        np.random.default_rng(1).normal(size=(500, 2)),
        np.random.default_rng(1).normal(size=500))
    print(f"A1 GP on pure noise: best NMSE={reg.best_nmse_:.3f} (1.0 == nothing)")
    return {"iid_raw_support_rate": float(np.mean(np.array(raw_support_runs) > 0)),
            "iid_gated_false_law_rate": gated_false / 30,
            "smooth_ou_support": n2, "smooth_ou_gate_nmse": float(gate2),
            "gp_best_nmse_on_noise": float(reg.best_nmse_)}


def A2_out_of_class():
    rng = np.random.default_rng(8000)
    X = rng.uniform(-3, 3, (800, 1))
    y = np.tanh(2 * X[:, 0])
    Xtr, ytr, Xte_, yte_ = X[:600], y[:600], X[600:], y[600:]
    reg = sr.SymbolicRegressor(1, population=400, generations=120, seed=5).fit(Xtr, ytr)
    te = sr.nmse(reg.best_, Xte_, yte_)
    nodes = sr.nodes(reg.best_)
    compact = nodes <= 25
    # domain-transfer check: does the discovered expression survive outside the
    # training interval? (tanh stays bounded; spurious fits typically blow up)
    Xex = np.concatenate([rng.uniform(3.5, 6.0, (300, 1)),
                          rng.uniform(-6.0, -3.5, (300, 1))])
    yex = np.tanh(2 * Xex[:, 0])
    with np.errstate(all="ignore"):
        ext_nmse = sr.nmse(reg.best_, Xex, yex)
    transfers = bool(np.isfinite(ext_nmse) and ext_nmse < 1e-2)
    law_found = bool(te < 1e-2 and compact and transfers)
    verdict = ("law found" if law_found else
               f"NOT a law: fits in-domain ({te:.1e} NMSE, {nodes} nodes) but "
               f"extrapolation NMSE={ext_nmse:.2e} -> domain-limited curve fit"
               if te < 1e-2 else
               f"no law (accuracy failure, NMSE={te:.1e})")
    print(f"A2 out-of-class y=tanh(2x): holdout NMSE={te:.3e}, nodes={nodes}, "
          f"extrapolation NMSE={ext_nmse:.3e} -> {verdict}")
    return {"holdout_nmse": float(te), "nodes": int(nodes),
            "compact_bar": 25, "extrapolation_nmse": float(ext_nmse)
            if np.isfinite(ext_nmse) else None,
            "transfers_out_of_domain": transfers, "verdict": verdict}


def A3_confounder():
    coefs_obs, r2_regime = [], []
    for rep in range(20):
        rg = np.random.default_rng(9000 + rep)
        z = rg.normal(size=2000)
        xi, eta = rg.normal(size=2000) * .3, rg.normal(size=2000) * .3
        x = 1.5 * z + xi
        y = -1.2 * z + eta
        A = np.stack([x, np.ones_like(x)], 1)
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        coefs_obs.append(c[0])
        z2 = rg.normal(size=2000) * 2.0
        x2 = rg.normal(size=2000)              # do(x): exogenous
        y2 = -1.2 * z2 + .3 * rg.normal(size=2000)
        A2 = np.stack([x2, np.ones_like(x2)], 1)
        c2, *_ = np.linalg.lstsq(A2, y2, rcond=None)
        ss_res = float(np.sum((y2 - A2 @ c2) ** 2))
        ss_tot = float(np.sum((y2 - y2.mean()) ** 2))
        r2_regime.append(1 - ss_res / max(ss_tot, 1e-300))
    obs_mean = float(np.mean(coefs_obs))
    int_mean = float(np.mean(r2_regime))
    print(f"A3 confounder: observational dY/dX ~= {obs_mean:+.3f} (strong spurious link)"
          f"; interventional mean R^2 = {int_mean:.4f} (link gone)")
    return {"observational_coef_mean": obs_mean,
            "observational_coef_std": float(np.std(coefs_obs)),
            "interventional_R2_mean": int_mean}


def A4_chaos_horizon():
    L = sy.lorenz(t_end=15.0)
    Xs, dX = sy.differentiate(L["X"], L["t"][1] - L["t"][0])
    Th, nm, fn = sd.build_library_named(Xs, ["x", "y", "z"], 3)
    fit = sd.fit_sindy(Th, dX)
    mdl = sd.DiscoveredModel(fit["C"], fn, 3)
    coef_err = float(mt.coeff_rel_err(fit["C"], _lorenz_true_C(nm)))

    lam = mt.lyapunov_benettin(sy.lorenz_f(), [-8., 7., 27.])
    fvec = sy.lorenz_f()

    def rollout(rhs, y0, t_eval):
        s = solve_ivp(rhs, (t_eval[0], t_eval[-1]), y0, t_eval=t_eval,
                      rtol=1e-11, atol=1e-13)
        return s.y.T if s.success and s.y.shape[1] == len(t_eval) else None

    rng = np.random.default_rng(4242)
    t_eval = np.arange(0, 25 + 1e-9, 0.01)
    h_disc, h_tiny, h_big = [], [], []
    for rep in range(30):
        y0 = rng.uniform(-15, 15, 3); y0[2] += 20
        Xt = rollout(lambda t, s: np.asarray(fvec(s[None, :])[0]), y0, t_eval)
        if Xt is None:
            continue
        Xm = mdl.simulate(y0, t_eval)
        if Xm.y.shape[1] == len(t_eval):
            hd = mt.predictability_horizon(Xm.y.T, Xt)
            if hd is not None:
                h_disc.append(t_eval[hd])
        for eps, store in ((1e-10, h_tiny), (1e-6, h_big)):
            yp = y0 + rng.normal(size=3) * eps
            Xp = rollout(lambda t, s: np.asarray(fvec(s[None, :])[0]), yp, t_eval)
            if Xp is not None:
                hb = mt.predictability_horizon(Xp, Xt)
                if hb is not None:
                    store.append(t_eval[hb])
    md = float(np.median(h_disc)) if h_disc else None
    m10 = float(np.median(h_tiny)) if h_tiny else None
    m6 = float(np.median(h_big)) if h_big else None
    # calibrate state-scale constant from the big-perturbation twin, predict others
    pred = {}
    if m6:
        C_const = np.exp(m6 * lam) * 1e-6
        for label, eps, meas in (("twin_1e-10", 1e-10, m10), ("discovered", None, md)):
            if meas is not None:
                e_eff = eps if eps else max(coef_err, 1e-12)
                pred[label] = float(np.log(C_const / e_eff) / lam)
    print(f"A4 chaos: lambda={lam:.3f} (literature 0.906); horizons medians: "
          f"discovered={md}, twin(1e-10)={m10}, twin(1e-6)={m6}")
    print(f"   theory-predicted horizons: {pred}")
    gap_meas = (m10 - m6) if (m10 and m6) else None
    gap_theory = float(np.log(1e-6 / 1e-10) / lam)
    print(f"   measured twin-horizon gap={gap_meas}s vs theory ln(1e4)/lam="
          f"{gap_theory:.2f}s")
    return {"lyapunov": lam, "median_horizon_discovered": md,
            "median_horizon_twin_1e-10": m10, "median_horizon_twin_1e-6": m6,
            "coef_max_rel_err": coef_err,
            "measured_gap_1e-10_vs_1e-6": gap_meas,
            "theory_gap_ln_ratio_over_lambda": gap_theory,
            "predicted_from_calibration": pred}


def _lorenz_true_C(names):
    Ct = np.zeros((len(names), 3))
    Ct[names.index("y"), 0] = 10; Ct[names.index("x"), 0] = -10
    Ct[names.index("x"), 1] = 28; Ct[names.index("y"), 1] = -1
    Ct[names.index("x*z"), 1] = -1
    Ct[names.index("x*y"), 2] = 1; Ct[names.index("z"), 2] = -8.0 / 3.0
    return Ct


def main():
    out = {"A1_noise_only": A1_noise_only(),
           "A2_out_of_class": A2_out_of_class(),
           "A3_confounder": A3_confounder(),
           "A4_chaos_horizon": A4_chaos_horizon()}
    save_results("exp08_adversarial.json", out)


if __name__ == "__main__":
    main()
