"""Evidence-chain checker: verifies every headline number quoted in
sessions/S6/REPORT.md against the committed results/*.json.

Run from experiments/S6:  python check_claims.py
Prints a PASS/FAIL table; exit code 1 on any FAIL.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
R = HERE / "results"

checks = []


def load(name):
    return json.loads((R / name).read_text())


def check(label, value, lo, hi):
    ok = lo <= value <= hi
    checks.append((label, value, f"[{lo}, {hi}]", ok))
    return ok


def main():
    e1 = load("exp01_known_laws.json")["runs"]
    by = {r["system"]: r for r in e1}
    check("exp01 osc coeff err", by["damped_oscillator"]["coeff_rel_err"], 0, 1e-9)
    check("exp01 LV coeff err", by["lotka_volterra"]["coeff_rel_err"], 0, 1e-9)
    check("exp01 lorenz coeff err", by["lorenz63"]["coeff_rel_err"], 0, 1e-9)
    for s in ("damped_oscillator", "lotka_volterra", "lorenz63"):
        check(f"exp01 {s} rollout", by[s]["rollout_rel_err"], 0, 5e-2)
        check(f"exp01 {s} jaccard", by[s]["support_jaccard"], 1.0, 1.0)

    e4 = load("exp04_period.json")
    rc = np.array(e4["series_coefficients"]["discovered"])
    th = np.array(e4["series_coefficients"]["perturbation_theory"])
    check("exp04 const rel err", abs(rc[0] / th[0] - 1), 0, 1e-5)
    check("exp04 th^2 rel err", abs(rc[1] / th[1] - 1), 0, 3e-3)
    check("exp04 collapse spread", e4["collapse"]["max_within_amplitude_spread"], 0, 1e-9)
    check("exp04 extrap@120 err", e4["extrapolation"]["120.0"]["rel_err_discovered"], 0, 0.03)
    check("exp04 unseen gL err", e4["unseen_gL"]["rel_err"], 0, 1e-5)

    e5 = load("exp05_kepler.json")
    c = e5["coefficients"]
    check("exp05 c1", c["log_a"], 1.4999999, 1.5000001)
    check("exp05 c2", c["log_mu"], -0.5000001, -0.4999999)
    check("exp05 |c_e|", max(abs(c["e"]), abs(c["e2"])), 0, 1e-6)
    check("exp05 heldout rms", e5["heldout_rms_resid"], 0, 1e-8)

    e6 = load("exp06_rlc.json")
    check("exp06 alpha coef", e6["alpha_formula"]["R_over_L_coef"], 0.498, 0.502)
    check("exp06 w2 invLC coef", e6["omega_formula"]["inv_LC_coef"], 0.999, 1.001)
    check("exp06 w2 R2L2 coef", e6["omega_formula"]["R2_L2_coef"], -0.251, -0.249)
    check("exp06 holdout omega_d", e6["heldout_max_rel_err"]["omega_d"], 0, 5e-3)

    e7 = load("exp07_invariants.json")
    check("exp07 duffing cos", e7["duffing"]["cos_angle_to_H_direction"], 1 - 1e-11, 1.0000000001)
    check("exp07 duffing resid/rangeH", e7["duffing"]["resid_over_rangeH"], 0, 1e-7)
    check("exp07 Lz cos", e7["central_force"]["cos_angle_to_Lz_direction"], 1 - 1e-11, 1.0000000001)
    check("exp07 neg damped rejected", int(e7["negative_controls"]["damped_delta0.3"]["correctly_rejected"]), 1, 1)
    check("exp07 neg forced rejected", int(e7["negative_controls"]["forced_F0.5"]["correctly_rejected"]), 1, 1)

    e8 = load("exp08_adversarial.json")
    check("exp08 gated false laws", e8["A1_noise_only"]["iid_gated_false_law_rate"], 0, 0)
    check("exp08 tanh in-domain nmse", e8["A2_out_of_class"]["holdout_nmse"], 0, 1e-3)
    ext = e8["A2_out_of_class"].get("extrapolation_nmse")
    if ext is not None:
        check("exp08 tanh extrapolation nmse", ext, 0.5, 10.0)
    check("exp08 confounder link", e8["A3_confounder"]["observational_coef_mean"], -1.0, -0.5)
    check("exp08 interventional R2", e8["A3_confounder"]["interventional_R2_mean"], 0, 0.01)
    lam = e8["A4_chaos_horizon"]["lyapunov"]
    check("exp08 lyapunov range", lam, 0.75, 1.15)

    e9 = load("exp09_baselines.json")
    stl = e9["dynamics"]["0.01"]["methods"]["stlsq_bic"]
    ols = e9["dynamics"]["0.01"]["methods"]["ols_full"]
    check("exp09 stlsq rollout sigma=.01", stl["rollout_rel_err"][0], 0, 5e-2)
    check("exp09 ours beats OLS rollout",
          1 if stl["rollout_rel_err"][0] < ols["rollout_rel_err"][0] else 0, 1, 1)

    width = max(len(c[0]) for c in checks) + 2
    fails = 0
    for label, val, rng_, ok in checks:
        print(f"{'PASS' if ok else 'FAIL'}  {label:<{width}} value={val!r:>24}  expected {rng_}")
        fails += (not ok)
    print(f"\n{len(checks)-fails}/{len(checks)} checks PASS")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
