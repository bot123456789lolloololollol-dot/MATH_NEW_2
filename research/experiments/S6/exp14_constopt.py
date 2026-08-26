"""Exp 14 — does constant optimization fix the root-form miss? (round 4, R4).

A/B test: identical GP search with vs without periodic numerical refinement of
ephemeral constants (cyclic coordinate descent on the top-10 individuals every 10
generations). 30 paired seeds x 2 root-form tasks. See PREREGISTERED_R4.md.
"""
import sys
import numpy as np

sys.path.insert(0, ".")
from src import symreg as sr                                        # noqa: E402
from src.common import save_results, exp_seeds                      # noqa: E402
from scipy.stats import wilcoxon                                    # noqa: E402


def refine_constants(tree, X, y, iters=25, step0=0.3):
    """Cyclic coordinate descent over the tree's constants (in-place)."""
    positions = []

    def collect(t, pos=()):
        if t[0] == "const":
            positions.append(pos)
        else:
            for i, ch in enumerate(t[1:] if t[0] in sr.BINOPS + sr.UNOPS else ()):
                if isinstance(ch, tuple):
                    collect(ch, pos + (i + 1,))
    collect(tree)

    def cost():
        e = sr.nmse(tree, X, y)
        return min(e, 1e12)

    best = cost()
    step = step0
    for _ in range(iters):
        improved = False
        for pos in positions:
            node = list(sr._get(tree, pos))
            c0 = float(node[1])
            for direction in (+1.0, -1.0):
                cand = ("const", round(c0 * (1.0 + direction * step)
                                       + direction * step * 0.01, 10))
                tree2 = sr._replace(tree, pos, cand)
                e = min(sr.nmse(tree2, X, y), 1e12)
                if e < best:
                    best, tree = e, tree2
                    improved = True
                    break
        if not improved:
            step *= 0.5
            if step < 1e-6:
                break
    return tree


class SymbolicRegressorOpt(sr.SymbolicRegressor):
    """Identical search; refines constants of the top-10 elites every 10 gens."""

    def fit(self, X, y):
        n = len(y)
        logn = np.log(n)
        vars_avail = list(range(self.n_vars))
        rng = self.rng

        def cost_of(tr):
            e = sr.nmse(tr, X, y)
            mse = e * max(np.var(y), sr.EPS)
            return n * np.log(max(mse, 1e-300)) + self.kappa * sr.nodes(tr) * logn, e

        pop = [sr.random_tree(rng, vars_avail, self.max_depth)
               for _ in range(self.pop_size)]
        scored = [(t,) + cost_of(t) for t in pop]

        for gen in range(self.generations):
            scored.sort(key=lambda z: z[1])
            if scored[0][2] < self.early_nmse:
                break
            elite_n = max(2, self.pop_size // 20)
            newpop = [scored[i][0] for i in range(elite_n)]

            if gen % 10 == 9:
                refined = []
                for t in scored[:10]:
                    refined.append(refine_constants(t[0], X, y))
                newpop[:len(refined)] = refined

            def tour():
                idx = rng.integers(0, len(scored), self.k)
                return min((scored[i] for i in idx), key=lambda z: z[1])[0]

            while len(newpop) < self.pop_size:
                r = rng.random()
                a = tour()
                if r < self.p_cx:
                    child = sr.crossover(rng, a, tour())
                elif r < self.p_cx + self.p_mut:
                    child = sr.mutate(rng, a, vars_avail)
                else:
                    child = a
                if sr.depth(child) > self.max_depth + 3:
                    continue
                newpop.append(child)

            scored = [(t,) + cost_of(t) for t in newpop]

        scored.sort(key=lambda z: z[1])
        self.best_ = scored[0][0]
        self.best_cost_, self.best_nmse_ = scored[0][1], scored[0][2]
        front = {}
        for t, c, e in scored:
            k = sr.nodes(t)
            if k not in front or e < front[k][1]:
                front[k] = (t, e)
        self.pareto_front_ = dict(sorted(front.items()))
        return self


def task_data(task, seed):
    rng = np.random.default_rng(seed)
    if task == "T1_invsqrt":
        x = rng.uniform(0.5, 2.0, 500)
        f = lambda z: 6.2832 / np.sqrt(z)
    else:
        x = rng.uniform(-2.0, 2.0, 500)
        f = lambda z: 2.0 / (1.0 + z**2)
    xte = rng.uniform(x.min(), x.max(), 300)
    return x[:, None], f(x), xte[:, None], f(xte)


def equivalent(tree, f_true, lo, hi, tol=1e-6):
    Xg = np.linspace(lo + 0.02 * (hi - lo), hi - 0.02 * (hi - lo), 2048)[:, None]
    with np.errstate(all="ignore"):
        pred = sr.evaluate(tree, Xg)
        true = f_true(Xg)
    return bool(np.all(np.isfinite(pred))
                and np.max(np.abs(pred - true)) <= tol * max(np.max(np.abs(true)), 1e-9))


def main():
    results = {}
    for task, (lo, hi) in (("T1_invsqrt", (0.5, 2.0)), ("T2_lorentzian", (-2.0, 2.0))):
        rows = []
        for seed in exp_seeds(14)[:30]:
            Xtr, ytr, Xte, yte = task_data(task, seed)
            out = {}
            for mname, cls in (("base", sr.SymbolicRegressor),
                               ("opt", SymbolicRegressorOpt)):
                reg = cls(1, population=600, generations=150,
                          seed=seed % 2**31).fit(Xtr, ytr)
                nmse = sr.nmse(reg.best_, Xte, yte)
                succ = nmse < 1e-6 and equivalent(reg.best_,
                                                  (lambda z: 6.2832 / np.sqrt(z)) if task == "T1_invsqrt"
                                                  else (lambda z: 2.0 / (1 + z**2)),
                                                  lo, hi)
                out[mname] = {"holdout_nmse": min(nmse, 1e6), "success": bool(succ)}
            rows.append(out)
            print(f"{task} seed {seed%100}: base={out['base']['holdout_nmse']:.1e}"
                  f"({out['base']['success']}) opt={out['opt']['holdout_nmse']:.1e}"
                  f"({out['opt']['success']})")
        agg = {}
        for m in ("base", "opt"):
            rr = [r[m] for r in rows]
            agg[m] = {"success_rate": float(np.mean([r["success"] for r in rr])),
                      "median_holdout_nmse": float(np.median([r["holdout_nmse"] for r in rr]))}
        try:
            p = float(wilcoxon([r["base"]["holdout_nmse"] for r in rows],
                               [r["opt"]["holdout_nmse"] for r in rows]).pvalue)
        except ValueError:
            p = 1.0
        agg["wilcoxon_p"] = p
        results[task] = {"runs": rows, "aggregate": agg}
        print(f"{task}: success base={agg['base']['success_rate']:.2f} "
              f"opt={agg['opt']['success_rate']:.2f} | median NMSE "
              f"{agg['base']['median_holdout_nmse']:.2e} vs "
              f"{agg['opt']['median_holdout_nmse']:.2e} | p={p:.4f}")

    save_results("exp14_constopt.json", results)


if __name__ == "__main__":
    main()
