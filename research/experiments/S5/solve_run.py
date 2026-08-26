#!/usr/bin/env python3
"""
S5 solve driver. Solves a DIMACS CNF with a chosen engine and wall-clock limit.

Engines: cadical | glucose | maple  (python-sat) ; z3 (native API build from
sat_net.Enc, used as an independent second engine for small instances).

Usage:
  python solve_run.py --dimacs f.dimacs --engine cadical --timeout 300 \
      [--model-out net.json --n N --slots T --variant sort]
Prints one line: RESULT engine dimacs timeout seconds verdict
Verdict in {SAT, UNSAT, TIMEOUT, ERROR}
"""
import argparse
import json
import sys
import time
import threading


def solve_pysat(path, engine, timeout):
    from pysat.formula import CNF
    from pysat.solvers import Cadical153, Glucose4, MapleCM
    cls = {"cadical": Cadical153, "glucose": Glucose4, "maple": MapleCM}[engine]
    t0 = time.time()
    cnf = CNF(from_file=path)
    tl = time.time() - t0
    s = cls(bootstrap_with=cnf)
    box = {"verdict": None, "model": None}

    def run():
        try:
            r = s.solve_limited(expect_interrupt=True)
            if r is True:
                box["verdict"] = "SAT"
                box["model"] = s.get_model()
            elif r is False:
                box["verdict"] = "UNSAT"
            else:
                box["verdict"] = "TIMEOUT"
        except Exception:
            box["verdict"] = "ERROR"

    th = threading.Thread(target=run, daemon=True)
    th.start()
    th.join(timeout)
    if th.is_alive():
        try:
            s.interrupt()
            th.join(30)
        except Exception:
            pass
        if box["verdict"] is None:
            box["verdict"] = "TIMEOUT"
    el = time.time() - t0
    return box["verdict"], box["model"], el + tl


def solve_z3_native(n, T, variant, m, sb1, sb2, sb3, timeout):
    import z3
    sys.path.insert(0, ".")
    from sat_net import Enc
    enc = Enc(n, T, variant, m, sb1=sb1, sb2=sb2, sb3=sb3)
    t0 = time.time()
    vs = [z3.Bool("v%d" % k) for k in range(enc.nvars + 1)]
    s = z3.Solver()
    s.set("timeout", int(timeout * 1000))
    for cl in enc.clauses:
        s.add(z3.Or([vs[abs(l)] if l > 0 else z3.Not(vs[abs(l)]) for l in cl]))
    r = s.check()
    el = time.time() - t0
    if r == z3.sat:
        mdl = s.model()
        model = []
        for k in range(1, enc.nvars + 1):
            val = mdl.eval(vs[k], model_completion=True)
            if z3.is_true(val):
                model.append(k)
        return "SAT", model, el
    if r == z3.unsat:
        return "UNSAT", None, el
    return "TIMEOUT", None, el


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dimacs")
    ap.add_argument("--engine", default="cadical",
                    choices=["cadical", "glucose", "maple", "z3"])
    ap.add_argument("--timeout", type=float, default=300.0)
    ap.add_argument("--model-out")
    ap.add_argument("--n", type=int)
    ap.add_argument("--slots", type=int)
    ap.add_argument("--variant", default="sort")
    ap.add_argument("--m", type=int)
    ap.add_argument("--z3-native", action="store_true",
                    help="ignore --dimacs; rebuild encoding natively in z3")
    ap.add_argument("--no-sb1", action="store_true")
    ap.add_argument("--no-sb2", action="store_true")
    ap.add_argument("--no-sb3", action="store_true")
    args = ap.parse_args()

    if args.engine == "z3" or args.z3_native:
        v, model, el = solve_z3_native(args.n, args.slots, args.variant, args.m,
                                       not args.no_sb1, not args.no_sb2,
                                       not args.no_sb3, args.timeout)
    else:
        v, model, el = solve_pysat(args.dimacs, args.engine, args.timeout)

    if v == "SAT" and args.model_out and model is not None:
        from sat_net import extract_model
        net = extract_model(model, args.n, args.slots, args.variant)
        with open(args.model_out, "w") as f:
            json.dump(net, f)
    print("RESULT %s %s %s %.2fs" %
          (v, args.engine, args.dimacs or ("z3native n%d T%d %s" %
           (args.n, args.slots, args.variant)), el))


if __name__ == "__main__":
    main()
