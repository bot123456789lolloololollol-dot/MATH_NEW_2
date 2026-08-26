"""Online stream generators (arrival order matters; NOT sorted)."""
import numpy as np


def _cap(rng):
    return int(rng.choice([100, 200, 500, 1000]))


def st_uniform(rng, n=None, cap=None):
    cap = cap or _cap(rng)
    n = n or int(rng.integers(100, 1200))
    return rng.integers(1, cap + 1, size=n).astype(np.int32), cap


def st_discrete(rng, n=None, cap=750, values=(250, 333, 376)):
    cap = cap or _cap(rng)
    n = n or int(rng.integers(100, 1200))
    v = np.asarray(values)
    return v[rng.integers(0, len(v), size=n)].astype(np.int32), int(cap)


def st_bimodal(rng, n=None, cap=1000):
    cap = cap or _cap(rng)
    n = n or int(rng.integers(100, 1200))
    mask = rng.random(n) < 0.5
    s = np.where(mask, rng.integers(1, cap // 10 + 1, size=n),
                 rng.integers(cap // 2 + 1, cap + 1, size=n))
    return s.astype(np.int32), int(cap)


def st_clustered(rng, n=None, cap=1000):
    cap = cap or _cap(rng)
    n = n or int(rng.integers(100, 1200))
    fracs = rng.choice([0.15, 0.25, 0.34, 0.5, 0.51, 0.66, 0.75, 0.9], size=4,
                       replace=False)
    centers = np.maximum(1, (cap * fracs).astype(int))
    which = rng.integers(0, 4, size=n)
    spread = max(1, cap // 100)
    return np.clip(centers[which] + rng.integers(-spread, spread + 1, size=n),
                   1, cap).astype(np.int32), int(cap)


def st_half_plus_eps(rng, n=None, cap=1000):
    """Alternating >half items and fillers: classic online killer."""
    cap = cap or _cap(rng)
    n = n or int(rng.integers(100, 1200))
    e = max(1, cap // 100)
    half_cap = int(n * 0.4)
    big = rng.integers(cap // 2 + e, cap // 2 + 8 * e, size=half_cap)
    small = rng.integers(1, cap // 2 - 8 * e, size=n - half_cap)
    s = np.concatenate([big, small])
    rng.shuffle(s)
    return np.clip(s, 1, cap).astype(np.int32), int(cap)


def st_triplet_like(rng, n=None, cap=1000):
    """Sizes in (cap/4, cap/2]: every bin fits <=3 items; arrivals random."""
    cap = cap or _cap(rng)
    n = n or int(rng.integers(100, 1200))
    e = max(1, cap // 200)
    s = rng.integers(cap // 4 + e, cap // 2 + 1, size=n)
    return s.astype(np.int32), int(cap)


def st_sawtooth(rng, n=None, cap=1000):
    cap = cap or _cap(rng)
    n = n or int(rng.integers(100, 1200))
    period = int(rng.integers(3, 10))
    t = np.arange(n) % period
    frac = 0.15 + 0.8 * t / max(1, period - 1)
    return np.clip((cap * frac).astype(int) + rng.integers(-2, 3, size=n),
                   1, cap).astype(np.int32), int(cap)


STREAM_GENS = {
    "uniform": st_uniform,
    "discrete_3": lambda rng: st_discrete(rng, values=(250, 333, 376)),
    "discrete_5": lambda rng: st_discrete(
        rng, cap=1000, values=(100, 300, 490, 510, 700)),
    "bimodal": st_bimodal,
    "clustered": st_clustered,
    "half_plus_eps": st_half_plus_eps,
    "triplet_like": st_triplet_like,
    "sawtooth": st_sawtooth,
}

HELDOUT_STREAMS = {
    "uniform_cap777": lambda rng: st_uniform(rng, cap=777),
    "discrete_prime": lambda rng: st_discrete(
        rng, cap=997, values=(211, 307, 463, 509)),
    "triplet_like_cap333": lambda rng: st_triplet_like(rng, cap=333),
    "clustered_cap555": lambda rng: st_clustered(rng, cap=555),
}


def make_stream_set(rng, m, names=None):
    names = names or list(STREAM_GENS.keys())
    out = []
    for _ in range(m):
        nm = names[int(rng.integers(0, len(names)))]
        out.append((nm,) + STREAM_GENS[nm](rng))
    return out


# ---------------------------------------------------------------- real-benchmark streams
_ORLIB_CACHE = {}


def _orlib_multisets():
    if not _ORLIB_CACHE:
        import os
        from .suites import load_falkenauer
        try:
            insts = load_falkenauer(classes=("u120", "u250", "t60", "t120"))
        except Exception:
            insts = []
        for i in insts:
            _ORLIB_CACHE.setdefault(i["name"].split("_")[0][:1], []).append(
                (i["sizes"], i["cap"]))
    return _ORLIB_CACHE


def st_orlib(rng, cls):
    ms = _orlib_multisets()[cls]
    sizes, cap = ms[int(rng.integers(0, len(ms)))]
    s = np.array(sizes, dtype=np.int32)
    rng.shuffle(s)
    return s, int(cap)


def st_weibull(rng, n=None, cap=1000, shape=1.5):
    cap = cap or 1000
    n = n or int(rng.integers(200, 2000))
    x = rng.weibull(shape, size=n)
    x = x / x.max()
    return np.clip((x * cap).astype(np.int32), 1, cap), int(cap)


if _orlib_multisets():
    STREAM_GENS["orlib_u"] = lambda rng: st_orlib(rng, "u")
    STREAM_GENS["orlib_t"] = lambda rng: st_orlib(rng, "t")
STREAM_GENS["weibull15"] = lambda rng: st_weibull(rng, shape=1.5)
STREAM_GENS["weibull30"] = lambda rng: st_weibull(rng, shape=3.0)
