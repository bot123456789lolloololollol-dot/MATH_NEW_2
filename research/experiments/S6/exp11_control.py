"""Exp 11 — Control loop synthesized from the DISCOVERED dynamics (round 2).

Plant (true, hidden from identifier): driven pendulum
    theta'' = -sin(theta) - c*omega + u,   g/L=1, c=0.15.
Identifier sees only (theta, omega, u) triples under preregistered filtered-random
excitation and fits f_hat(theta, omega) in a POLYNOMIAL library that does NOT contain
sin -- i.e., the identified model is knowingly imperfect away from theta ~ 0.
Input gain b=1 assumed known (actuation physics), disclosed in PREREGISTERED_R2.md.

Controller: exact-feedback-linearization computed ONLY from f_hat:
    u = -f_hat(theta, omega) - kp*theta - kd*omega     (regulate theta -> 0)
Oracle comparison: identical structure with the TRUE f.
Metrics (preregistered): IAE over 10 s from (0.5, 0); settling time |e|<0.05;
success iff stable and IAE(discovered)/IAE(oracle) < 1.25.
Adversarial: same from (2.5, 0) -- outside the region where poly approximates sin well;
result reported without a success claim.
"""
import sys
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, ".")
from src import systems as sy, sindy as sd                          # noqa: E402
from src.common import save_results                                 # noqa: E402

C_TRUE, G_OVER_L = 0.15, 1.0


def plant_rhs(t, s, u_fun):
    th, om = s
    return [om, -np.sin(th) - C_TRUE * om + u_fun(t)]


def excitation(seed, t_eval):
    """Band-limited random torque: OU process, clipped.

    Amplitude limited (+/-0.4) so identification keeps |theta| < ~1.2 rad: the
    experiment's premise is control under bounded-angle model mismatch, not under
    full rotations where a polynomial premise is vacuous (first run wound theta to
    32 rad -- see PREREGISTERED_R2 deviations).
    """
    rng = np.random.default_rng(seed)
    dt = t_eval[1] - t_eval[0]
    x = np.zeros(len(t_eval))
    for k in range(1, len(t_eval)):
        x[k] = x[k - 1] + (-x[k - 1] * 1.0) * dt + rng.normal() * np.sqrt(dt) * 1.2
    return np.clip(x, -0.4, 0.4)


def identify_model():
    t_end, dt = 25.0, 5e-3
    t_eval = np.arange(0.0, t_end + 1e-9, dt)
    u_seq = excitation(1111, t_eval)
    sol = solve_ivp(lambda t, s: plant_rhs(t, s, lambda tt: np.interp(tt, t_eval, u_seq)),
                    (0, t_end), [0.6, 0.0], t_eval=t_eval, rtol=1e-10, atol=1e-12,
                    method="DOP853")
    X = sol.y.T                                   # theta, omega
    assert np.max(np.abs(X[:, 0])) < 2.0, \
        "identification excursion left the bounded-angle regime"
    Xs, dX = sy.differentiate(X, dt)              # aligned; dX[:,1] ~ omega'
    us = u_seq[6:-6]
    target = (dX[:, 1] - us)[:, None]             # subtract KNOWN input channel
    # degree-3 library: sin(theta)'s leading nonlinearity is cubic
    Theta, names, fn = sd.build_library_named(Xs, ["th", "om"], poly_degree=3)
    fit = sd.fit_sindy(Theta, target)
    C = fit["C"][:, 0]
    act = sd.significant_support(Theta, target, fit["C"], z=5.0)[:, 0]
    terms = {names[i]: round(float(C[i]), 6) for i in range(len(names)) if act[i]}
    return fn, C, terms


def make_controller(fn, C):
    def fhat(th, om):
        S = np.array([[th, om]])
        return float((fn(S) @ C)[0])

    KP, KD = 4.0, 2 * np.sqrt(4.0)

    def u(t, th, om):
        return -fhat(th, om) - KP * th - KD * om
    return u


def closed_loop(u_ctrl, y0, t_end=10.0, dt_out=1e-3):
    t_eval = np.arange(0.0, t_end + 1e-9, dt_out)
    sol = solve_ivp(lambda t, s: plant_rhs(t, s, lambda tt: u_ctrl(tt, s[0], s[1])),
                    (0, t_end), list(y0), t_eval=t_eval, rtol=1e-10, atol=1e-12,
                    method="DOP853")
    th = sol.y[0]
    err = np.abs(th)
    iae = float(np.trapezoid(err, sol.t))
    below = err < 0.05
    settle = None
    for k in range(len(below)):
        if below[k:].all():
            settle = float(sol.t[k])
            break
    return {"IAE": iae, "settling_time": settle,
            "stable": bool(np.all(np.isfinite(err)) and np.max(err) < 10.0)}


def main():
    fn, C, terms = identify_model()
    print("identified f_hat terms:", terms)

    ctrl_disc = make_controller(fn, C)
    ctrl_oracle = make_controller(lambda S: np.stack(
        [-np.sin(S[..., 0]) - C_TRUE * S[..., 1]], axis=-1), np.array([1.0]))

    res = {}
    for label, y0 in (("near", (0.5, 0.0)), ("far_adversarial", (2.5, 0.0))):
        d = closed_loop(ctrl_disc, y0)
        o = closed_loop(ctrl_oracle, y0)
        ratio = d["IAE"] / max(o["IAE"], 1e-300)
        success = bool(d["stable"] and ratio < 1.25) if label == "near" else None
        res[label] = {"discovered": d, "oracle": o, "iae_ratio": ratio,
                      "preregistered_success": success}
        print(f"{label}: discovered IAE={d['IAE']:.4f} settle={d['settling_time']} "
              f"| oracle IAE={o['IAE']:.4f} settle={o['settling_time']} "
              f"| ratio={ratio:.3f} success={success}")

    save_results("exp11_control.json", {"f_hat_terms": terms, "runs": res})


if __name__ == "__main__":
    main()
