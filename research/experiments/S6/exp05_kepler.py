"""Exp 05 — Kepler's third law from two-body simulations, plus an invariance result.

Training set: orbits over a grid of semi-major axis a, eccentricity e, and
gravitational parameter mu; the learner sees only (a, e, mu) -> measured period T.

Claims tested:
  C1 (law):  log T = alpha*log a + beta*log mu + gamma with alpha=1.5, beta=-0.5.
  C2 (invariance / "undocumented" absence-of-dependence): gamma_e = gamma_ee = 0 --
     the period does not depend on eccentricity. Discovering that a variable does
     NOT enter a law is as important as discovering the law itself.
  C3 (generalization): the discovered relation predicts held-out cells
     (unseen mu AND unseen e) to integration accuracy.
  C4 (falsifiability): under velocity drag the relation measurably breaks and the
     break is detected (no unfalsifiable "law").
"""
import sys
import numpy as np

sys.path.insert(0, ".")
from src import systems as sy                                       # noqa: E402
from src.common import save_results                                 # noqa: E402
from src.plotting import newfig, savefig                            # noqa: E402


def build_dataset():
    mus = [1.0, 1.3, 1.6]
    aa = np.linspace(0.6, 1.8, 7)
    es = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
    rows = []
    for mu in mus:
        for a in aa:
            for e in es:
                o = sy.two_body(mu=mu, a=a, e=e, n_orbits=4)
                rows.append((a, e, mu, o["T"], o["T_theory"]))
                assert abs(o["T"] / o["T_theory"] - 1) < 1e-8, "period measurement"
    return np.array(rows)


def fit_law(rows):
    # features: [1, log a, log mu, e, e^2]; target: log T
    X = np.stack([np.ones(len(rows)), np.log(rows[:, 0]), np.log(rows[:, 2]),
                  rows[:, 1], rows[:, 1] ** 2], axis=1)
    y = np.log(rows[:, 3])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    return coef, float(np.sqrt(np.mean(resid**2)))


def main():
    rows = build_dataset()
    # train: mu in {1.0, 1.6} and e <= 0.3 ; test: everything else
    train_mask = ((np.isin(rows[:, 2], [1.0, 1.6])) & (rows[:, 1] <= 0.30))
    coef, rms_train = fit_law(rows[train_mask])

    Xall = np.stack([np.ones(len(rows)), np.log(rows[:, 0]), np.log(rows[:, 2]),
                     rows[:, 1], rows[:, 1] ** 2], axis=1)
    resid_all = np.log(rows[:, 3]) - Xall @ coef
    rms_heldout = float(np.sqrt(np.mean(resid_all[~train_mask] ** 2)))

    print("discovered log T = c0 + c1*log(a) + c2*log(mu) + c3*e + c4*e^2")
    names = ["c0 (expect log 2pi)", "c1 (expect 1.5)", "c2 (expect -0.5)",
             "c3 (expect 0)", "c4 (expect 0)"]
    for nme, c in zip(names, coef):
        print(f"   {nme}: {c:+.8f}")
    print(f"train RMS(resid)={rms_train:.3e}  held-out RMS={rms_heldout:.3e}")

    # eccentricity-invariance significance: |c3| CI via bootstrap on train rows
    rng = np.random.default_rng(5000)
    tr = rows[train_mask]
    boot = []
    for _ in range(200):
        idx = rng.integers(0, len(tr), len(tr))
        boot.append(fit_law(tr[idx])[0])
    boot = np.array(boot)
    ci_c3 = np.percentile(boot[:, 3], [2.5, 97.5]).tolist()
    ci_c4 = np.percentile(boot[:, 4], [2.5, 97.5]).tolist()
    print(f"c3 95% CI: {ci_c3}  contains 0: {ci_c3[0] <= 0 <= ci_c3[1]}")
    print(f"c4 95% CI: {ci_c4}  contains 0: {ci_c4[0] <= 0 <= ci_c4[1]}")
    # practical equivalence (the CI at integration precision is far tighter than
    # any physical meaning of "zero"): |coef| below 1e-6 counts as no dependence
    e_invariant = bool(abs(coef[3]) < 1e-6 and abs(coef[4]) < 1e-6)
    print("eccentricity practically invariant (|c|<1e-6):", e_invariant)

    # ---- falsification: drag breaks Kepler and the break is detected
    def drag(t, s):
        return [-0.05 * s[2], -0.05 * s[3]]
    o = sy.two_body(mu=1.0, a=1.2, e=0.3, n_orbits=6, drag=drag)
    X = o["X"]; t = o["t"]
    r = np.hypot(X[:, 0], X[:, 1])
    mins = []
    for i in range(2, len(r) - 2):
        if r[i] < r[i-1] and r[i] < r[i+1] and r[i] <= r[i-2] and r[i] <= r[i+2]:
            if not mins or t[i] - mins[-1] > 0.3:   # one periapsis per orbit
                mins.append(t[i])
    drift, apparent = None, None
    if len(mins) >= 4:
        periods = np.diff(mins[:4])                 # first three orbital periods only
        apparent = float(np.mean(periods))
        drift = float((periods[-1] - periods[0]) / periods[0])
    kepler_pred = 2 * np.pi * np.sqrt(1.2**3 / 1.0)
    print(f"drag run: first-orbits mean period {apparent} vs clean-Kepler {kepler_pred}; "
          f"per-orbit relative drift={drift}")
    falsified = drift is not None and abs(drift) > 1e-3
    print("Kepler relation flagged violated under drag:", bool(falsified))

    # ---- figure: collapse plot T vs a^{3/2}/sqrt(mu); if C1+C2 hold, all data
    # falls on the identity line
    fig, ax = newfig()
    xax = rows[:, 0] ** 1.5 / np.sqrt(rows[:, 2])
    ax.loglog(xax[train_mask], rows[train_mask, 3], "o", ms=4, label="train")
    ax.loglog(xax[~train_mask], rows[~train_mask, 3], "x", ms=5, label="held out")
    lims = [xax.min() * 0.9, xax.max() * 1.1]
    ax.loglog(lims, lims, "r--", lw=1,
              label=r"identity: $T=a^{3/2}/\sqrt{\mu}$")
    ax.set_xlabel(r"$a^{3/2}/\sqrt{\mu}$")
    ax.set_ylabel("measured period T")
    ax.legend(fontsize=8)
    savefig(fig, "figures/exp05_kepler.png")

    save_results("exp05_kepler.json", {
        "coefficients": {"const": coef[0], "log_a": coef[1], "log_mu": coef[2],
                         "e": coef[3], "e2": coef[4]},
        "expected_const_log2pi": float(np.log(2 * np.pi)),
        "bootstrap_ci_e_coeff": ci_c3, "bootstrap_ci_e2_coeff": ci_c4,
        "train_rms_resid": rms_train, "heldout_rms_resid": rms_heldout,
        "eccentricity_practically_invariant": e_invariant,
        "eccentricity_ci_excludes_zero": bool(not (ci_c3[0] <= 0 <= ci_c3[1])),
        "note": "CIs are at integrator precision (~1e-10); practical-equivalence "
                "threshold 1e-6 is the meaningful invariance claim",
        "drag_falsification": {"apparent_period_first_orbits": apparent,
                               "kepler_prediction": kepler_pred,
                               "relative_per_orbit_drift": drift,
                               "flagged_violated": bool(falsified)}})


if __name__ == "__main__":
    main()
