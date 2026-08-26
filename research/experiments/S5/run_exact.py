"""Driver: search exact minimum sizes for k-selection networks.

For each (n,k): linear search upward from lo; records every SAT witness
(decoded network, re-verified by F1/F2 verifiers) and every UNSAT certificate.
Results append to results_exact.jsonl.
"""
import json
import sys
import time
from threading import Timer

from exact_selector import decide
import verify_selector as vs


def decide_guarded(n, k, L, prop, limit):
    out = {}

    def target():
        out["stop"] = True

    import exact_selector as ES
    cnf, pool, actions = ES.build_and_encode(n, k, L, prop)
    from pysat.solvers import Cadical153
    s = Cadical153(bootstrap_with=cnf.clauses)
    tm = Timer(limit, s.interrupt)
    t0 = time.time()
    tm.start()
    r = s.solve()
    tm.cancel()
    if r is True:
        net = ES.decode(s.get_model(), pool, actions, L)
        return "SAT", net, time.time() - t0
    if r is False:
        return "UNSAT", None, time.time() - t0
    return "UNKNOWN", None, time.time() - t0


def main():
    jobs = [(6, 3), (7, 3), (8, 4), (9, 4), (10, 5)]
    if len(sys.argv) > 1:
        jobs = [tuple(map(int, j.split(","))) for j in sys.argv[1:]]
    limit = float(sys.argv[-1]) if sys.argv[-1].replace('.', '', 1).isdigit() \
        else 900.0
    logf = open("results_exact.jsonl", "a")
    for (n, k) in jobs:
        lo = 2
        while True:
            r, net, dt = decide_guarded(n, k, lo, "select", limit)
            logf.write(json.dumps({"n": n, "k": k, "L": lo, "res": r,
                                   "secs": round(dt, 1)}) + "\n")
            logf.flush()
            print(f"SEL({n},{k}) L={lo}: {r} [{dt:.1f}s]", flush=True)
            if r == "SAT":
                ok = vs.sel_f1(net, n, k) and vs.sel_f2(net, n, k)
                logf.write(json.dumps({"n": n, "k": k, "L": lo,
                                       "witness_valid": bool(ok),
                                       "net": [list(e) for e in net]}) + "\n")
                logf.flush()
                print(f"  witness verified: {ok}, net={net}", flush=True)
                break
            if r == "UNKNOWN":
                print(f"  budget exhausted; stopping this job at L={lo}",
                      flush=True)
                break
            lo += 1


if __name__ == "__main__":
    main()
