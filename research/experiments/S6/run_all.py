"""Run the full S6 experiment suite. One command regenerates every result.

    python run_all.py
"""
import sys
import time
import traceback

sys.path.insert(0, ".")

EXPS = [
    "exp01_known_laws",
    "exp02_noise",
    "exp03_simplification",
    "exp04_pendulum_period",
    "exp05_kepler",
    "exp06_rlc",
    "exp07_invariants",
    "exp08_adversarial",
    "exp09_baselines",
    "exp10_pde",
    "exp11_control",
]


def main():
    only = sys.argv[1:] or EXPS
    for name in only:
        mod = __import__(name)
        t0 = time.time()
        print(f"\n===== {name} =====")
        try:
            mod.main()
            print(f"[ok] {name} in {time.time()-t0:.1f}s")
        except Exception:
            print(f"[FAIL] {name}")
            traceback.print_exc()
            sys.exit(1)


if __name__ == "__main__":
    main()
