"""Shared helpers: paths, seeds, result IO."""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
RESULTS = EXP / "results"
FIGURES = EXP / "figures"
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

MASTER_SEED = 20260826


def exp_seeds(exp_id, n=30):
    return [1000 * exp_id + i for i in range(n)]


def save_results(name, payload):
    def _default(o):
        if isinstance(o, (np.floating, np.integer)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(str(type(o)))
    p = RESULTS / name
    p.write_text(json.dumps(payload, indent=1, default=_default))
    print(f"[saved] {p}")
    return p


def add_src_to_path():
    for p in (str(HERE), str(EXP)):
        if p not in sys.path:
            sys.path.insert(0, p)
