"""Conserved-quantity discovery from trajectories alone (no knowledge of the RHS).

Method: if a scalar function F(s) is conserved, then for consecutive samples of one
trajectory F(s_{k+1}) - F(s_k) = 0.  With F linear in a feature dictionary Phi,
this is a homogeneous linear system A c = 0 with
    A = [ Phi(S_1) - Phi(S_0) ; ... ]   (stacked over all training trajectories),
solved by SVD: the right singular vector with smallest singular value gives the
invariant; the ratio sigma_min/sigma_max and the spectral gap provide a built-in
detection test ("is there any conserved combination at all?").
"""
import numpy as np


def discover_invariant(traj_list, feat_fn):
    """traj_list: list of state arrays (n_i, d) sampled consecutively.

    Returns dict with coefficient vector c (normalized to unit norm), singular values,
    detection statistics, and the scalar field F(s) as a callable.
    """
    rows = []
    for S in traj_list:
        P = feat_fn(S)
        rows.append(P[1:] - P[:-1])
    A = np.vstack(rows)

    # column SCALING only (no centering: the constraint must stay homogeneous,
    # i.e. A c = 0 exactly -- centering would accept constant-increment functions)
    sd = A.std(axis=0)
    sd[sd < 1e-12] = 1.0
    As = A / sd

    U, sv, Vt = np.linalg.svd(As, full_matrices=False)
    c_std = Vt[-1]
    c = c_std / sd                       # back to raw feature space

    # residual statistics on the standardized system
    ratio_min = float(sv[-1] / sv[0])
    gap = float(sv[-2] / sv[-1]) if len(sv) > 1 else np.inf

    def F(S):
        return feat_fn(np.atleast_2d(S)) @ c

    return {"c": c, "singular_values": sv.tolist(), "sigma_ratio": ratio_min,
            "spectral_gap": gap, "F": F, "n_constraints": len(A)}


def invariant_drift(F, traj):
    """Relative drift of F along a held-out trajectory: std / range."""
    vals = F(traj)
    rng = float(np.max(vals) - np.min(vals))
    return float(np.std(vals) / max(rng, 1e-300))


def affine_match(F_vals, H_vals):
    """Best affine fit a*H+b ~ F; returns (a, b, R^2, resid_std_over_range)."""
    A = np.stack([H_vals, np.ones_like(H_vals)], axis=1)
    coef, *_ = np.linalg.lstsq(A, F_vals, rcond=None)
    pred = A @ coef
    ss_res = float(np.sum((F_vals - pred) ** 2))
    ss_tot = float(np.sum((F_vals - np.mean(F_vals)) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1e-300)
    resid_rel = float(np.std(F_vals - pred) / max(np.ptp(F_vals), 1e-300))
    return float(coef[0]), float(coef[1]), float(r2), resid_rel


def poly_features(X, degree=4, var_names=None):
    """Polynomial feature dictionary (with names), same construction as sindy."""
    import itertools
    n, d = X.shape
    names = []

    def term_fn(exponent):
        def fn(S):
            out = np.ones(len(S))
            for j, e in enumerate(exponent):
                if e:
                    out = out * S[:, j] ** e
            return out
        return fn

    terms, cols = [], []
    vn = var_names or [f"s{j}" for j in range(d)]
    for deg in range(0, degree + 1):
        for expo in itertools.combinations_with_replacement(range(d), deg):
            exponent = [0] * d
            for j in expo:
                exponent[j] += 1
            label = "*".join(vn[j] + (f"^{exponent[j]}" if exponent[j] > 1 else "")
                             for j in sorted(set(expo))) or "1"
            terms.append(term_fn(tuple(exponent)))
            names.append(label)
            cols.append(term_fn(tuple(exponent))(X))
    Phi = np.column_stack(cols)
    # drop constant column: it cannot contribute to differences but pollutes conditioning
    keep = ~(np.allclose(Phi, Phi[0]))
    return Phi[:, keep], [nm for nm, k in zip(names, keep) if k]
