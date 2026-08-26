"""Subprocess-isolated exact decision steps (CaDiCaL has no interrupt() in
python-sat; isolate instead of interrupt).

Usage: python run_exact2.py <n> <k> <startL> <maxL> <timeout_sec>
Writes one JSON line per step to results_exact.jsonl.
"""
import json
import subprocess
import sys
import time


CHILD = r'''
import json, sys
sys.path.insert(0, ".")
from exact_selector import decide
import verify_selector as vs
n, k, L = map(int, sys.argv[1:4])
r, net = decide(n, k, L)
out = {"res": r}
if r == "SAT":
    out["net"] = [list(e) for e in net]
    out["witness_valid"] = bool(vs.sel_f1(net, n, k) and vs.sel_f2(net, n, k))
json.dump(out, open(sys.argv[4], "w"))
'''


def main():
    n, k, lo, hi, tmo = map(int, sys.argv[1:6])
    logf = open("results_exact.jsonl", "a")
    with open("_exact_child.py", "w") as f:
        f.write(CHILD)
    for L in range(lo, hi + 1):
        t0 = time.time()
        try:
            subprocess.run([sys.executable, "_exact_child.py",
                            str(n), str(k), str(L), "_step_out.json"],
                           timeout=tmo, check=True)
            out = json.load(open("_step_out.json"))
            rec = {"n": n, "k": k, "L": L, "res": out["res"],
                   "secs": round(time.time() - t0, 1)}
            logf.write(json.dumps(rec) + "\n")
            logf.flush()
            print(f"SEL({n},{k}) L={L}: {out['res']} "
                  f"[{rec['secs']}s]", flush=True)
            if out["res"] == "SAT":
                rec2 = {"n": n, "k": k, "L": L,
                        "witness_valid": out.get("witness_valid"),
                        "net": out.get("net")}
                logf.write(json.dumps(rec2) + "\n")
                logf.flush()
                print(f"  witness valid={out.get('witness_valid')} "
                      f"net={out.get('net')}", flush=True)
                break
        except subprocess.TimeoutExpired:
            logf.write(json.dumps({"n": n, "k": k, "L": L,
                                   "res": "TIMEOUT",
                                   "secs": round(time.time() - t0, 1)}) + "\n")
            logf.flush()
            print(f"SEL({n},{k}) L={L}: TIMEOUT", flush=True)
            break


if __name__ == "__main__":
    main()
