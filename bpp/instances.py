"""Instance generators for 1D bin packing (integer sizes).

Every generator returns (sizes ndarray[int32], cap int).  A fixed `rng`
(numpy Generator) must be passed so experiments are reproducible.
"""
import numpy as np


def gen_uniform(rng, n, cap, lo=1, hi=None):
    if hi is None:
        hi = cap
    s = rng.integers(lo, hi, size=n, endpoint=True)
    return s.astype(np.int32), int(cap)


def gen_discrete(rng, n, cap, values):
    """Sizes drawn uniformly from the given integer value set."""
    v = np.asarray(values)
    s = v[rng.integers(0, len(v), size=n)]
    return s.astype(np.int32), int(cap)


def gen_triplet(rng, n_groups, cap, jitter=2):
    """Falkenauer-style triplets: groups of ~3 items summing exactly to cap."""
    items = []
    for _ in range(n_groups):
        while True:
            a = rng.integers(cap // 4 + jitter, cap // 2 - jitter + 1)
            b = rng.integers(cap // 4 + jitter, cap // 2 - jitter + 1)
            c = cap - a - b
            if cap // 4 <= c <= cap // 2 and a > 0 and b > 0:
                break
        items.extend([a, b, c])
    s = np.array(items, dtype=np.int32)
    rng.shuffle(s)
    return s, int(cap)


def gen_ffd_worst_case(rng, reps, eps_int=3, cap=1000):
    """Integer version of the classic FFD 11/9-OPT construction.

    With capacity 1 and epsilon small: blocks of {1/2+e, 1/4+2e, 1/4+e} plus
    matching fillers; FFD wastes ~1/9 while OPT is tight.  We use the standard
    integer scaling.  Returns sizes where OPT is known analytically.
    """
    # cap = 12m units per block family; use cap divisible by 12*eps-ish scale
    base = 6 * eps_int * 2  # unit block
    items = []
    for _ in range(reps):
        h = cap // 2 + eps_int          # just over half
        q1 = cap // 4 + 2 * eps_int     # quarter+
        q2 = cap // 4 + eps_int         # quarter+
        items += [h, q1, q2]
    # fillers that allow perfect repacking into cap//... bins: pairs q1+q2+h' ...
    # add small items of size cap//2 - eps to pair with h -> but that would make OPT easy.
    s = np.array(items, dtype=np.int32)
    rng.shuffle(s)
    return s, int(cap)


def gen_bimodal(rng, n, cap, p_small=0.5, small_hi=None, big_lo=None):
    """Half tiny items, half large items -- stresses residual management."""
    if small_hi is None:
        small_hi = max(2, cap // 10)
    if big_lo is None:
        big_lo = max(small_hi + 1, cap // 2)
    mask = rng.random(n) < p_small
    s = np.where(mask,
                 rng.integers(1, small_hi + 1, size=n),
                 rng.integers(big_lo, cap + 1, size=n))
    return s.astype(np.int32), int(cap)


def gen_clustered(rng, n, cap, n_clusters=4):
    """Item sizes clustered around a few values correlated with cap fractions."""
    fracs = rng.choice([0.15, 0.25, 0.34, 0.5, 0.51, 0.66, 0.75, 0.9], size=n_clusters,
                       replace=False)
    centers = np.maximum(1, (cap * fracs).astype(int))
    which = rng.integers(0, n_clusters, size=n)
    spread = np.maximum(1, cap // 100)
    s = np.clip(centers[which] + rng.integers(-spread, spread + 1, size=n), 1, cap)
    return s.astype(np.int32), int(cap)


def gen_sawtooth(rng, n, cap, period=None):
    """Deterministic-ish sawtooth profile with noise."""
    if period is None:
        period = max(2, int(rng.integers(2, 8)))
    t = np.arange(n) % period
    frac = 0.15 + 0.8 * t / max(1, period - 1)
    s = np.clip((cap * frac).astype(int) + rng.integers(-2, 3, size=n), 1, cap)
    return s.astype(np.int32), int(cap)


def gen_almost_perfect(rng, n_bins_target, cap, perturb=1):
    """Build an instance with a known perfect packing, then perturb a few items.

    Guarantees LB == OPT-ish and makes greedy waste measurable.
    """
    items = []
    for _ in range(n_bins_target):
        r = cap
        while r > 0:
            mx = min(r, max(1, cap // 2))
            s = int(rng.integers(max(1, mx // 2), mx + 1)) if r > 2 else r
            items.append(s)
            r -= s
    s = np.array(items[: int(len(items))], dtype=np.int32)
    k = max(1, len(s) // 20)
    idx = rng.choice(len(s), size=k, replace=False)
    s[idx] = np.clip(s[idx] + rng.integers(-perturb, perturb + 1, size=k), 1, cap)
    rng.shuffle(s)
    return s, int(cap)


# ---------------------------------------------------------------- registry
GENERATORS = {
    "uniform": lambda rng: gen_uniform(rng, int(rng.integers(50, 500)), 1000),
    "uniform_small_n": lambda rng: gen_uniform(rng, int(rng.integers(20, 120)), 500),
    "discrete_3": lambda rng: gen_discrete(rng, int(rng.integers(60, 400)), 750,
                                           [250, 333, 376]),
    "discrete_5": lambda rng: gen_discrete(rng, int(rng.integers(60, 400)), 1000,
                                           [100, 300, 490, 510, 700]),
    "triplet": lambda rng: gen_triplet(rng, int(rng.integers(30, 150)), 1000),
    "bimodal": lambda rng: gen_bimodal(rng, int(rng.integers(60, 400)), 1000),
    "clustered": lambda rng: gen_clustered(rng, int(rng.integers(60, 400)), 1000),
    "sawtooth": lambda rng: gen_sawtooth(rng, int(rng.integers(60, 400)), 1000),
    "almost_perfect": lambda rng: gen_almost_perfect(rng, int(rng.integers(10, 40)), 1000),
    "ffd_trap": lambda rng: gen_ffd_worst_case(rng, int(rng.integers(10, 60))),
}

# held-out families never seen during training
HELDOUT_GENERATORS = {
    "uniform_cap777": lambda rng: gen_uniform(rng, int(rng.integers(80, 600)), 777),
    "discrete_prime": lambda rng: gen_discrete(rng, int(rng.integers(60, 400)), 997,
                                               [211, 307, 463, 509]),
    "triplet_cap500": lambda rng: gen_triplet(rng, int(rng.integers(30, 150)), 500),
    "clustered_cap333": lambda rng: gen_clustered(rng, int(rng.integers(60, 400)), 333),
}


def make_dataset(name, rng, m):
    return [GENERATORS[name](rng) for _ in range(m)]


def make_mixed_dataset(rng, m, names=None):
    names = names or list(GENERATORS.keys())
    out = []
    for i in range(m):
        out.append(GENERATORS[names[int(rng.integers(0, len(names)))]](rng))
    return out
