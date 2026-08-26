
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
