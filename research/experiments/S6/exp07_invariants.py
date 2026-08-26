"""Exp 07 — Conserved-quantity discovery from trajectories alone + calibrated rejection.

Positive cases:
  * conservative Duffing  x'' = x - x^3      -> H = v^2/2 + x^4/4 - x^2/2
  * planar central-force motion           -> angular momentum Lz = x*vy - y*vx
Method: polynomial feature dictionary; invariance constraint F(s_{k+1}) - F(s_k)=0
solved by SVD on stacked increments (multiple trajectories required -- see report).
Negative controls (must be REJECTED by the preregistered detection criterion
sigma_ratio < 1e-8 AND spectral gap > 100):
  * damped Duffing   delta = 0.3
  * forced Duffing   F = 0.5, omega = 1
Calibration sweep: dissipation delta in {0, 0.01, ..., 0.3} -> where does detection fail?
"""
import sys
import numpy as np

sys.path.insert(0, ".")
from src import systems as sy, invariant as inv                     # noqa: E402
from src.common import save_results                                 # noqa: E402
from src.plotting import newfig, savefig                            # noqa: E402

DETECTED = lambda r: bool(r["sigma_ratio"] < 1e-8 and r["spectral_gap"] > 100)


def duffing_case():
    train = [sy.duffing(x0=1.8)["X"], sy.duffing(x0=-0.4, v0=0.9)["X"]]
    feat = lambda S: inv.poly_features(S, degree=4)[0]
    res = inv.discover_invariant(train, feat)
    pool = [sy.duffing(x0=x0, v0=v0) for x0, v0 in
            [(1.8, 0.0), (-0.4, 0.9), (0.7, -1.1), (2.2, 0.0), (-1.9, 0.4)]]
    held = [sy.duffing(x0=x0, v0=v0)["X"] for x0, v0 in
            [(0.7, -1.1), (2.2, 0.0), (-1.9, 0.4)]]
    F = res["F"]
    cons_held = inv.conservation_ratio(F, held)
    Fp = np.concatenate([F(T["X"]) for T in pool])
    Hp = np.concatenate([T["H_true"] for T in pool])
    a, b, r2, resid = inv.affine_match(Fp, Hp)
    # analytic direction comparison in feature space
    _, pnames = inv.poly_features(pool[0]["X"], degree=4)
    hdir = np.zeros(len(pnames))
    hdir[pnames.index("v^2")] = 0.5
    hdir[pnames.index("x^2")] = -0.5
    hdir[pnames.index("x^4")] = 0.25
    cosang = float(abs(res["c"] @ hdir) / (np.linalg.norm(res["c"]) * np.linalg.norm(hdir)))
    return {"sigma_ratio": res["sigma_ratio"], "spectral_gap": res["spectral_gap"],
            "heldout_conservation_ratio": cons_held,
            "affine_R2_vs_true_H": r2,
            "resid_over_rangeH": resid,
            "cos_angle_to_H_direction": cosang,
            "detected": DETECTED(res)}


def central_force_case():
    train = [sy.central_force_2d(r0=(1.0, 0.3), v0=(0.2, 1.1))["X"],
             sy.central_force_2d(r0=(1.5, -0.2), v0=(-0.3, 0.8))["X"],
             sy.central_force_2d(r0=(0.8, 0.6), v0=(0.9, -0.4))["X"]]
    feat = lambda S: inv.poly_features(S, degree=2)[0]
    res = inv.discover_invariant(train, feat)
    _, pnames = inv.poly_features(train[0], degree=2)
    Ldir = np.zeros(len(pnames))
    Ldir[pnames.index("s0*s3")] = 1.0     # x * vy
    Ldir[pnames.index("s1*s2")] = -1.0    # y * vx
    cosang = float(abs(res["c"] @ Ldir) / (np.linalg.norm(res["c"]) * np.linalg.norm(Ldir)))
    held = [sy.central_force_2d(r0=(1.1, -0.5), v0=(0.4, 1.0))["X"],
            sy.central_force_2d(r0=(0.6, 0.9), v0=(-1.0, 0.3))["X"]]
    return {"sigma_ratio": res["sigma_ratio"], "spectral_gap": res["spectral_gap"],
            "heldout_conservation_ratio": inv.conservation_ratio(res["F"], held),
            "cos_angle_to_Lz_direction": cosang,
            "detected": DETECTED(res),
            "note": "energy is NOT polynomial (contains 1/r); only Lz is recoverable"}


def negative_controls():
    out = {}
    feat = lambda S: inv.poly_features(S, degree=4)[0]
    res_d = inv.discover_invariant([sy.duffing(delta=0.3)["X"]], feat)
    out["damped_delta0.3"] = {"sigma_ratio": res_d["sigma_ratio"],
                              "spectral_gap": res_d["spectral_gap"],
                              "correctly_rejected": not DETECTED(res_d)}
    def forced(t, s):  # integrate via generic call below
        pass
    simF = sy.duffing(forcing=0.5, omega=1.0, t_end=50.0)
    # forcing is time-dependent: conservation claim needs autonomous features;
    # we test whether the algorithm wrongly claims an invariant from state-only data
    res_f = inv.discover_invariant([simF["X"]], feat)
    out["forced_F0.5"] = {"sigma_ratio": res_f["sigma_ratio"],
                          "spectral_gap": res_f["spectral_gap"],
                          "correctly_rejected": not DETECTED(res_f)}
    return out


def calibration_sweep():
    feat = lambda S: inv.poly_features(S, degree=4)[0]
    rows = []
    for delta in (0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3):
        trajs = [sy.duffing(delta=delta, x0=1.8)["X"],
                 sy.duffing(delta=delta, x0=-0.4, v0=0.9)["X"]]
        res = inv.discover_invariant(trajs, feat)
        rows.append({"delta": delta, **{"sigma_ratio": res["sigma_ratio"],
                    "spectral_gap": res["spectral_gap"],
                    "detected": DETECTED(res)}})
        print(f"delta={delta:.3f}: sigma_ratio={rows[-1]['sigma_ratio']:.2e} "
              f"gap={rows[-1]['spectral_gap']:.2e} detected={rows[-1]['detected']}")
    return rows


def main():
    print("== Duffing conservative ==")
    du = duffing_case()
    for k, v in du.items():
        print(f"   {k}: {v}")
    print("== Central force ==")
    cf = central_force_case()
    for k, v in cf.items():
        print(f"   {k}: {v}")
    print("== Negative controls ==")
    nc = negative_controls()
    for k, v in nc.items():
        print(f"   {k}: {v}")
    print("== Dissipation calibration sweep ==")
    sweep = calibration_sweep()

    fig, ax = newfig()
    ds = [r["delta"] for r in sweep]
    ax.semilogy(ds, [max(r["sigma_ratio"], 1e-16) for r in sweep], "o-",
                label="sigma_min/sigma_max")
    ax.axhline(1e-8, color="red", ls="--", label="detection threshold")
    ax.set_xlabel("damping delta")
    ax.set_ylabel("null-space residual")
    ax.set_title("when is conservation still detectable?")
    ax.legend(fontsize=8)
    savefig(fig, "figures/exp07_calibration.png")

    save_results("exp07_invariants.json",
                 {"duffing": du, "central_force": cf,
                  "negative_controls": nc, "calibration": sweep})


if __name__ == "__main__":
    main()
