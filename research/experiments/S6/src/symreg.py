"""Compact deterministic genetic-programming symbolic regression.

Programs are nested tuples:
    ('op', child, ...)          ops: add, sub, mul, div*, sin, cos, exp*, log*, sqrt*
    ('var', j)                  j-th input variable
    ('const', value)            ephemeral mutable constant
(* = protected variants, safe under division-by-zero / log of non-positives.)

Model selection is MDL/BIC-style: cost = n*ln(MSE) + kappa * nodes * ln(n),
and a full Pareto front (nodes vs NMSE) is returned so the Occam choice is explicit.
"""
import numpy as np

EPS = 1e-9
BINOPS = ("add", "sub", "mul", "div")
UNOPS = ("sin", "cos", "exp", "log", "sqrt")


# ------------------------------------------------------------------ evaluation
def evaluate(tree, X):
    kind = tree[0]
    if kind == "var":
        return X[:, tree[1]].copy()
    if kind == "const":
        return np.full(len(X), float(tree[1]))
    if kind in BINOPS:
        a = evaluate(tree[1], X)
        b = evaluate(tree[2], X)
        if kind == "add":
            return a + b
        if kind == "sub":
            return a - b
        if kind == "mul":
            return a * b
        safe_b = np.where(np.abs(b) > EPS, b, 1.0)
        return np.where(np.abs(b) > EPS, a / safe_b, 1.0)
    if kind in UNOPS:
        a = evaluate(tree[1], X)
        if kind == "sin":
            return np.sin(a)
        if kind == "cos":
            return np.cos(a)
        if kind == "exp":
            return np.exp(np.clip(a, -20.0, 20.0))
        if kind == "log":
            return np.log(np.abs(a) + EPS)
        return np.sqrt(np.abs(a))
    raise ValueError(f"bad node {kind}")


def nodes(tree):
    if tree[0] in ("var", "const"):
        return 1
    return 1 + sum(nodes(c) for c in tree[1:])


def depth(tree):
    if tree[0] in ("var", "const"):
        return 1
    return 1 + max(depth(c) for c in tree[1:])


def _safe(v):
    return np.all(np.isfinite(v))


def nmse(tree, X, y):
    with np.errstate(all="ignore"):
        pred = evaluate(tree, X)
    if not _safe(pred):
        return 1e12
    r = float(np.mean((pred - y) ** 2)) / max(float(np.var(y)), EPS)
    return r if np.isfinite(r) else 1e12


# ------------------------------------------------------------------ genetics
def _rand_const(rng):
    choice = rng.random()
    if choice < 0.7:
        return round(float(rng.uniform(-3, 3)), 4)
    return round(float(10 ** rng.uniform(-2, 1)), 4)


def random_tree(rng, vars_avail, max_depth=4, p_const=0.25):
    if max_depth <= 0 or rng.random() < 0.15:
        if vars_avail and rng.random() < 1 - p_const:
            return ("var", int(rng.choice(vars_avail)))
        return ("const", _rand_const(rng))
    if rng.random() < 0.7:
        op = str(rng.choice(BINOPS))
        return (op,
                random_tree(rng, vars_avail, max_depth - 1, p_const),
                random_tree(rng, vars_avail, max_depth - 1, p_const))
    op = str(rng.choice(UNOPS))
    return (op, random_tree(rng, vars_avail, max_depth - 1, p_const))


def _subtree_positions(tree, prefix=()):
    out = [prefix]
    if tree[0] not in ("var", "const"):
        for i, ch in enumerate(tree[1:]):
            out.extend(_subtree_positions(ch, prefix + (i + 1,)))
    return out


def _get(tree, pos):
    for i in pos:
        tree = tree[i]
    return tree


def _replace(tree, pos, new):
    if not pos:
        return new
    kids = list(tree[1:])
    kids[pos[0] - 1] = _replace(kids[pos[0] - 1], pos[1:], new)
    return (tree[0],) + tuple(kids)


def _pick(rng, positions):
    return positions[int(rng.integers(0, len(positions)))]


def crossover(rng, a, b):
    pa = _pick(rng, _subtree_positions(a))
    pb = _pick(rng, _subtree_positions(b))
    return _replace(a, pa, _get(b, pb))


def mutate(rng, t, vars_avail, max_depth=4):
    r = rng.random()
    pos = _pick(rng, _subtree_positions(t))
    if r < 0.6:                                   # subtree replacement
        return _replace(t, pos, random_tree(rng, vars_avail, 2))
    node = list(_get(t, pos))
    if r < 0.8 and node[0] == "const":            # constant jitter
        node[1] = float(node[1]) * float(rng.uniform(0.5, 1.5)) + float(rng.normal(0, 0.1))
        return _replace(t, pos, ("const", round(node[1], 4)))
    if node[0] in BINOPS and len(node) == 3:      # turn binop into unop wrapper or drop
        if rng.random() < 0.5:
            return _replace(t, pos, (str(rng.choice(UNOPS)), (node[0],) + tuple(node[1:])))
        return _replace(t, pos, node[rng.integers(1, 3)])
    return t


def _depth_cap(t, cap=8):
    """Truncate overdeep offspring."""
    if depth(t) <= cap:
        return t
    return random_tree(np.random.default_rng(0), (), 2)  # tiny fallback


class SymbolicRegressor:
    def __init__(self, n_vars, population=300, generations=80, seed=0,
                 p_crossover=0.7, p_subtree_mut=0.2, tournament=3, kappa=1.0,
                 max_depth=5, early_nmse=1e-13):
        self.n_vars = n_vars
        self.pop_size = population
        self.generations = generations
        self.rng = np.random.default_rng(seed)
        self.p_cx, self.p_mut = p_crossover, p_subtree_mut
        self.k = tournament
        self.kappa = kappa
        self.max_depth = max_depth
        self.early_nmse = early_nmse

    def fit(self, X, y):
        n = len(y)
        logn = np.log(n)
        vars_avail = list(range(self.n_vars))

        def cost_of(tr):
            e = nmse(tr, X, y)
            mse = e * max(np.var(y), EPS)
            return n * np.log(max(mse, 1e-300)) + self.kappa * nodes(tr) * logn, e

        pop = [random_tree(self.rng, vars_avail, self.max_depth)
               for _ in range(self.pop_size)]
        scored = [(t,) + cost_of(t) for t in pop]
        best_hist = []
        self.pareto_ = {}

        for gen in range(self.generations):
            scored.sort(key=lambda z: z[1])
            if scored[0][2] < self.early_nmse:
                break
            elite_n = max(2, self.pop_size // 20)
            newpop = [scored[i][0] for i in range(elite_n)]

            def tour():
                idx = self.rng.integers(0, len(scored), self.k)
                return min((scored[i] for i in idx), key=lambda z: z[1])[0]

            while len(newpop) < self.pop_size:
                r = self.rng.random()
                a = tour()
                if r < self.p_cx:
                    b = tour()
                    child = crossover(self.rng, a, b)
                elif r < self.p_cx + self.p_mut:
                    child = mutate(self.rng, a, vars_avail)
                else:
                    child = a
                if depth(child) > self.max_depth + 3:
                    continue
                newpop.append(child)

            scored = []
            for t in newpop:
                c, e = cost_of(t)
                scored.append((t, c, e))

            front = {}
            for t, c, e in sorted(scored, key=lambda z: (nodes(z[0]), z[2])):
                k = nodes(t)
                if k not in front or e < front[k]:
                    front[k] = e
            self.pareto_[gen] = dict(front)
            best_hist.append(scored[0][2])

        scored.sort(key=lambda z: z[1])
        self.best_ = scored[0][0]
        self.best_cost_, self.best_nmse_ = scored[0][1], scored[0][2]

        # Pareto front at final generation: minimal NMSE per node count
        front = {}
        for t, c, e in scored:
            k = nodes(t)
            if k not in front or e < front[k][1]:
                front[k] = (t, e)
        self.pareto_front_ = {k: v for k, v in sorted(front.items())}
        return self


# ------------------------------------------------------------------ sympy bridge
def to_sympy(tree, var_names=None):
    import sympy as sp
    names = var_names or ["x0", "x1", "x2", "x3"]

    def conv(t):
        k = t[0]
        if k == "var":
            return sp.Symbol(names[t[1]])
        if k == "const":
            return sp.Float(t[1])
        if k == "add":
            return conv(t[1]) + conv(t[2])
        if k == "sub":
            return conv(t[1]) - conv(t[2])
        if k == "mul":
            return conv(t[1]) * conv(t[2])
        if k == "div":
            return conv(t[1]) / conv(t[2])
        if k == "sin":
            return sp.sin(conv(t[1]))
        if k == "cos":
            return sp.cos(conv(t[1]))
        if k == "exp":
            return sp.exp(conv(t[1]))
        if k == "log":
            return sp.log(sp.Abs(conv(t[1])) + EPS)
        if k == "sqrt":
            return sp.sqrt(sp.Abs(conv(t[1])))
        raise ValueError(k)
    return conv(tree)
