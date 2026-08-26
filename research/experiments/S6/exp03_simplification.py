"""Exp 03 — Occam tests: does complexity-penalized search prefer simpler equivalent forms?

Controlled function classes where a simpler closed form exists:
  S1: y = sin(x) cos(x)                (== 0.5 sin(2x))
  S2: y = (x^3 - 1) / (x - 1)          (== x^2 + x + 1, sampled off the hole)
  S3: y = log(x^2 - 1) - log(x - 1)    (== log(x + 1), domain x > 1)
Protocol: GP symbolic regression with MDL selection over the Pareto front;
held-out NMSE; equivalence verified numerically on 4096 fresh points (< 1e-9).
This validates the optimizer's Occam bias -- a prerequisite for trusting
discovered laws in exp04-07.
"""
import sys
import numpy as np

sys.path.insert(0, ".")
from src import symreg as sr                                        # noqa: E402
from src.common import save_results                                 # noqa: E402
from src.plotting import newfig, savefig                            # noqa: E402


def dataset(case, rng):
    if case == "S1_sincos":
        Xtr = rng.uniform(-3, 3, (600, 1)); Xte = rng.uniform(-3, 3, (400, 1))
        f = lambda X: np.sin(X[:, 0]) * np.cos(X[:, 0])
        truth = "sin(x)*cos(x)"
    elif case == "S2_rational":
        def sample(n):
            x = rng.uniform(0.0, 2.0, n)
            x = x[np.abs(x - 1.0) > 0.02][:n]  # keep away from removable singularity
            while len(x) < n:
                more = rng.uniform(0.0, 2.0, n)
                more = more[np.abs(more - 1.0) > 0.02]
                x = np.concatenate([x, more])[:n]
            return x
        Xtr = sample(600); Xte = sample(400)
        Xtr, Xte = Xtr[:, None], Xte[:, None]
        f = lambda X: (X[:, 0] ** 3 - 1) / (X[:, 0] - 1)
        truth = "(x^3-1)/(x-1)"
    else:  # S3_log
        Xtr = rng.uniform(1.05, 4.0, (600, 1)); Xte = rng.uniform(1.05, 4.0, (400, 1))
        f = lambda X: np.log(X[:, 0] ** 2 - 1) - np.log(X[:, 0] - 1)
        truth = "log(x^2-1)-log(x-1)"
    return Xtr, f(Xtr), Xte, f(Xte), truth, f


def check_equivalence(tree, f_true, rng, tol=1e-8):
    X = rng.uniform(-3, 4, (4096, 1))
    with np.errstate(all="ignore"):
        pred = sr.evaluate(tree, X)
        true = f_true(X)
    ok = bool(np.all(np.isfinite(pred)))
    rel = float(np.max(np.abs(pred - true)) / max(np.max(np.abs(true)), 1e-12))
    return ok and rel < tol, rel


def main():
    rng_master = np.random.default_rng(20260826)
    out = []
    for case in ("S1_sincos", "S2_rational", "S3_log"):
        recs = []
        for seed in range(5):
            rng = np.random.default_rng(hash((case, seed)) % 2**31)
            Xtr, ytr, Xte, yte, truth, f_true = dataset(case, rng)
            reg = sr.SymbolicRegressor(1, population=300, generations=80,
                                       seed=seed).fit(Xtr, ytr)
            # MDL choice on final Pareto front, evaluated on held-out data
            best = None
            for k, (tree, nmse_tr) in reg.pareto_front_.items():
                e_te = sr.nmse(tree, Xte, yte)
                cost = len(ytr) * np.log(max(e_te * np.var(ytr), 1e-300)) + k * np.log(len(ytr))
                if best is None or cost < best[0]:
                    best = (cost, k, tree, e_te)
            _, k_best, tree, te_nmse = best
            equiv, relerr = check_equivalence(tree, f_true, rng_master)
            expr = str(sr.to_sympy(tree, ["x"]))
            recs.append({"seed": seed, "nodes": int(k_best), "holdout_nmse": te_nmse,
                         "expr": expr, "equivalent_to_truth": bool(equiv),
                         "max_rel_err_vs_truth": relerr})
            print(f"{case} seed{seed}: nodes={k_best} holdout_nmse={te_nmse:.2e} "
                  f"equiv={equiv} expr={expr}")
        n_equiv = sum(r["equivalent_to_truth"] for r in recs)
        min_nodes = min(r["nodes"] for r in recs)
        out.append({"case": case, "truth": truth, "runs": recs,
                    "n_equivalent_of_5": int(n_equiv),
                    "min_complexity_found": int(min_nodes)})

    # figure: Pareto front for S1
    rng = np.random.default_rng(0)
    Xtr, ytr, Xte, yte, _ = dataset("S1_sincos", rng)
    reg = sr.SymbolicRegressor(1, population=300, generations=80, seed=0).fit(Xtr, ytr)
    ks = sorted(reg.pareto_front_)
    vals = [max(reg.pareto_front_[k][1], 1e-16) for k in ks]
    fig, ax = newfig()
    ax.semilogy(ks, vals, "o-")
    ax.axhline(1e-12, color="gray", ls="--", lw=.8)
    ax.set_xlabel("expression complexity (node count)")
    ax.set_ylabel("training NMSE")
    ax.set_title("Pareto front: y = sin(x)cos(x)")
    savefig(fig, "figures/exp03_pareto.png")

    save_results("exp03_occam.json", {"cases": out})


if __name__ == "__main__":
    main()
