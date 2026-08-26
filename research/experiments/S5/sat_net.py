#!/usr/bin/env python3
"""
S5 encoder: existence of a comparator network with <= T comparators on n channels
satisfying a target property, encoded to CNF (DIMACS).

Variants:
  sort    : full sorting network (all 2^n binary inputs must emerge ascending)
  median  : designated output channel must carry the median bit of the input
            (ascending median position idx = n//2 for odd n)
  merge   : first m channels start ascending, rest start ascending; full output
            must be ascending. Only the (m+1)*(nm-m+1) constant input pairs are
            enumerated.

Soundness core (0-1 principle): a comparator network sorts all totally-ordered
inputs iff it sorts all binary inputs. We enumerate ALL binary inputs except the
(n+1) ascending-canonical ones 0^k 1^(n-k), which provably stay ascending through
any comparators (comparators preserve the ascending-sortedness invariant).

Semantics encoding is TOTAL: for every channel c and slot l, whatever single
choice (comparator covering c / comparator not covering c / empty) is active,
the next-level value w[l+1][c] is fully constrained (min / max / copy resp.).

Symmetry breaking (each preserves satisfiability-for-existence; see SPEC.md):
  --sb1 : adjacent slots holding disjoint comparators may be assumed lex-
          nondecreasing (disjoint comparators commute; bubble-sort argument).
  --sb2 : first slot's comparator assumed (channel1, channel2) (channel
          relabeling bijection).
  --sb3 : consecutive slots never hold the identical pair (second occurrence is
          a no-op; removal + re-bubbling terminates).
SB2/SB3/SB1 apply to sort/median. For merge only sb1/sb3 are offered (channel
roles block relabeling).

Usage:
  python sat_net.py cnf  --n 8 --slots 19 [--variant sort] [--m 2] \
      [--no-sb1 --no-sb2 --no-sb3] --out cnf.dimacs [--stats]
  python sat_net.py model --dimacs cnf.dimacs --n 8 --slots 19 --variant sort \
      [--m 2] > net.json        # extract comparator list from a SAT model line
"""
import argparse
import itertools
import json
import sys


def pairs_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


class Enc:
    def __init__(self, n, T, variant="sort", m=None, sb1=True, sb2=True, sb3=True):
        self.n, self.T, self.variant, self.m = n, T, variant, m
        self.sb1, self.sb2, self.sb3 = sb1, sb2, sb3
        self.pairs = pairs_of(n)
        self.P = len(self.pairs)
        self.clauses = []
        self._build_inputs()
        self._alloc_vars()
        self._encode()

    # ---- input enumeration ----------------------------------------------
    def _build_inputs(self):
        n = self.n
        if self.variant == "merge":
            m = self.m
            ins = []
            for a in range(m + 1):
                left = [1 if k >= m - a else 0 for k in range(m)]
                for b in range(n - m + 1):
                    right = [1 if k >= (n - m) - b else 0 for k in range(n - m)]
                    ins.append(tuple(left + right))
            self.inputs = ins
        elif self.variant == "sort":
            canon = set()
            for k in range(n + 1):
                canon.add(tuple([0] * (n - k) + [1] * k))
            self.inputs = [v for v in itertools.product([0, 1], repeat=n) if v not in canon]
        else:  # median
            canon = set()
            for k in range(n + 1):
                canon.add(tuple([0] * (n - k) + [1] * k))
            self.inputs = [v for v in itertools.product([0, 1], repeat=n) if v not in canon]

    # ---- variables --------------------------------------------------------
    def _alloc_vars(self):
        n, T, P = self.n, self.T, self.P
        self.nv_x0 = 1                       # x[l][p] : l*T? -> l*P+p+1
        self.x = lambda l, p: l * P + p + 1
        self.e = lambda l: T * P + l + 1     # empty-slot var
        off = T * P + T
        self.w = {}
        # w[(vi,l,c)]
        self.wid = lambda vi, l, c: off + (vi * (T + 1) + l) * n + c + 1
        self.nvars = off + len(self.inputs) * (T + 1) * n

    def C(self, lits):
        self.clauses.append(lits)

    # ---- encoding ---------------------------------------------------------
    def _encode(self):
        n, T, P = self.n, self.T, self.P
        pairs = self.pairs
        # slot choice: exactly-one over {x[l][p]} + {e[l]}
        for l in range(T):
            lits = [self.x(l, p) for p in range(P)] + [self.e(l)]
            self.C(lits[:])                      # at-least-one
            for a in range(len(lits)):
                for b in range(a + 1, len(lits)):
                    self.C([-lits[a], -lits[b]])  # at-most-one
        # SB2: first comparator is (1,2) -> pair (0,1)
        if self.sb2 and self.variant in ("sort", "median") and T > 0:
            self.C([self.x(0, 0)])
        # SB3: no consecutive duplicate pair
        if self.sb3:
            for l in range(T - 1):
                for p in range(P):
                    self.C([-self.x(l, p), -self.x(l + 1, p)])
        # SB1: adjacent disjoint comparators in lex-nondecreasing order
        if self.sb1:
            for l in range(T - 1):
                for pi in range(P):
                    for pj in range(P):
                        if pi != pj and _disjoint(pairs[pi], pairs[pj]) and pj < pi:
                            self.C([-self.x(l, pi), -self.x(l + 1, pj)])
        # simulation semantics (total)
        for vi, v in enumerate(self.inputs):
            for l in range(T):
                # initial level constants handled by unit clauses at l==0
                for c in range(n):
                    if l == 0:
                        lit = self.wid(vi, 0, c)
                        self.C([lit if v[c] else -lit])
                    cur = self.wid(vi, l, c)
                    nxt = self.wid(vi, l + 1, c)
                    # comparators covering c: min/max definitions
                    for p, (i, j) in enumerate(pairs):
                        xv = self.x(l, p)
                        if i == c:
                            o = self.wid(vi, l, j)
                            self._min_def(xv, nxt, cur, o)
                        elif j == c:
                            o = self.wid(vi, l, i)
                            self._max_def(xv, nxt, cur, o)
                        else:
                            # comparator not touching c -> copy
                            self.C([-xv, -nxt, cur])
                            self.C([-xv, nxt, -cur])
                    # empty slot -> copy
                    ev = self.e(l)
                    self.C([-ev, -nxt, cur])
                    self.C([-ev, nxt, -cur])
        # property
        if self.variant == "sort":
            # forbid descent out[k]=1, out[k+1]=0
            for vi in range(len(self.inputs)):
                for k in range(n - 1):
                    self.C([-self.wid(vi, T, k), self.wid(vi, T, k + 1)])
        elif self.variant == "median":
            mid = self.n // 2  # ascending median index for odd n
            for vi, v in enumerate(self.inputs):
                mb = sorted(v)[mid]
                self.C([self.wid(vi, T, mid)] if mb else [-self.wid(vi, T, mid)])
        else:  # merge: full ascending output on constant inputs
            for vi in range(len(self.inputs)):
                for k in range(n - 1):
                    self.C([-self.wid(vi, T, k), self.wid(vi, T, k + 1)])

    def _min_def(self, x, na, a, b):
        # na == min(a,b) under activation x
        self.C([-x, -na, a])
        self.C([-x, -na, b])
        self.C([-x, na, -a, -b])

    def _max_def(self, x, nx, a, b):
        # nx == max(a,b) under activation x
        self.C([-x, -a, nx])
        self.C([-x, -b, nx])
        self.C([-x, -nx, a, b])


def _disjoint(p, q):
    return p[0] not in q and p[1] not in q


def write_dimacs(enc, path):
    with open(path, "w") as f:
        f.write("p cnf %d %d\n" % (enc.nvars, len(enc.clauses)))
        buf = []
        for cl in enc.clauses:
            buf.append(" ".join(map(str, cl)) + " 0")
            if len(buf) > 20000:
                f.write("\n".join(buf) + "\n")
                buf = []
        if buf:
            f.write("\n".join(buf) + "\n")


def extract_model(model_lits, n, T, variant):
    """model_lits: list of nonzero ints. Returns comparator list (1-based chans),
    stripping empty slots."""
    P = n * (n - 1) // 2
    s = set(l for l in model_lits if l > 0)
    comps = []
    for l in range(T):
        picked = None
        for p in range(P):
            if l * P + p + 1 in s:
                i, j = None, None
                # recover pair from index
                idx = 0
                for i in range(n):
                    for j in range(i + 1, n):
                        if idx == p:
                            picked = (i + 1, j + 1)
                        idx += 1
        if picked:
            comps.append(list(picked))
    return {"n": n, "variant": variant, "comparators": comps}


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c1 = sub.add_parser("cnf")
    c1.add_argument("--n", type=int, required=True)
    c1.add_argument("--slots", type=int, required=True)
    c1.add_argument("--variant", default="sort", choices=["sort", "median", "merge"])
    c1.add_argument("--m", type=int, default=None)
    c1.add_argument("--out", required=True)
    c1.add_argument("--no-sb1", action="store_true")
    c1.add_argument("--no-sb2", action="store_true")
    c1.add_argument("--no-sb3", action="store_true")
    c2 = sub.add_parser("model")
    c2.add_argument("--dimacs", required=True)
    c2.add_argument("--n", type=int, required=True)
    c2.add_argument("--slots", type=int, required=True)
    c2.add_argument("--variant", default="sort")
    args = ap.parse_args()

    if args.cmd == "cnf":
        enc = Enc(args.n, args.slots, args.variant, args.m,
                  sb1=not args.no_sb1, sb2=not args.no_sb2, sb3=not args.no_sb3)
        write_dimacs(enc, args.out)
        print("n=%d T=%d variant=%s inputs=%d vars=%d clauses=%d -> %s"
              % (args.n, args.slots, args.variant, len(enc.inputs),
                 enc.nvars, len(enc.clauses), args.out))
    else:
        best = None
        with open(args.dimacs) as f:
            for line in f:
                if line.startswith("v ") or line.startswith("SAT") or line[0] in "-123456789":
                    toks = line.split()
                    vals = [int(t) for t in toks if t not in ("v", "SAT")]
                    if vals:
                        best = vals
        comps = extract_model(best, args.n, args.slots, args.variant)
        print(json.dumps(comps))


if __name__ == "__main__":
    main()
