#!/usr/bin/env python3
"""S2 instance generator: Taillard-distribution-faithful P||Cmax instances.

Source for distribution: E. Taillard, "Benchmarks for basic scheduling problems",
European Journal of Operational Research 64 (1993) 278-285. His P||Cmax benchmark
processing times are i.i.d. discrete uniform integers in [1, 99]. We reproduce that
distribution exactly (his original PRNG stream/state is not specified in the paper,
so we use our own fixed-seed numpy PCG64 streams; instances are committed as data).

Deterministic: p = default_rng(seed).integers(1, 100, size=n)  -> values in {1..99}.
File format: line 1 "m n seed"; line 2: n integers.
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "benchmarks", "S2-instances"))

MS = [2, 3, 4, 5, 7, 10]
NS = [30, 50, 100, 200, 500]
DEV_SEEDS = list(range(1000, 1030))   # tuning only (see PREREGISTERED.md)
EVAL_SEEDS = list(range(3000, 3030))  # final evaluation only


def gen_instance(m: int, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.integers(1, 100, size=n).astype(np.int64)


def write_instance(path: str, m: int, n: int, seed: int, p: np.ndarray) -> None:
    with open(path, "w") as f:
        f.write(f"{m} {n} {seed}\n")
        f.write(" ".join(map(str, p.tolist())) + "\n")


def read_instance(path: str):
    with open(path) as f:
        first = f.readline().split()
        m, n, seed = int(first[0]), int(first[1]), int(first[2])
        p = np.array(list(map(int, f.read().split())), dtype=np.int64)
    assert len(p) == n, f"{path}: expected {n} processing times, got {len(p)}"
    return m, n, seed, p


def main():
    os.makedirs(OUT, exist_ok=True)
    count = 0
    for m in MS:
        for n in NS:
            for phase, seeds in (("dev", DEV_SEEDS), ("eval", EVAL_SEEDS)):
                for s in seeds:
                    path = os.path.join(OUT, f"taillike-m{m}-n{n}-s{s}.txt")
                    if os.path.exists(path):
                        continue
                    write_instance(path, m, n, s, gen_instance(m, n, s))
                    count += 1
    print(f"wrote {count} instance files to {OUT}")


if __name__ == "__main__":
    sys.exit(main())
