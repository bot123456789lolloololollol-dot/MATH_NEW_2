#!/usr/bin/env python3
"""
Independent brute-force verifier for comparator networks. Shares NO logic with
the SAT encoder: it reads a JSON network and simulates every relevant input.

Properties:
  sort   : all 2^n binary inputs emerge ascending (0-1 principle => full sorter)
  median : channel n//2 carries the ascending-median bit for every binary input
  merge  : given first m channels initially ascending and the rest ascending,
           output is fully sorted, over ALL (m+1)*(nm-m+1) constant pairs
           (exhaustive over the actual input space of a merger)

Usage:
  python checker.py --net net.json                 # property from net["variant"]
  python checker.py --net net.json --property sort
  python checker.py --net net.json --m 2 --property merge
Exit code 0 = PASS, 1 = FAIL. Prints details.

Self-test (validates the checker itself on known-good and known-bad nets):
  python checker.py --selftest
"""
import argparse
import itertools
import json
import sys


def simulate(n, comparators, inp):
    w = list(inp)
    for (i, j) in comparators:
        if not (1 <= i < j <= n):
            raise ValueError("bad comparator (%d,%d)" % (i, j))
        a, b = w[i - 1], w[j - 1]
        w[i - 1], w[j - 1] = min(a, b), max(a, b)
    return w


def check(net, property_override=None):
    n = net["n"]
    comps = [tuple(c) for c in net["comparators"]]
    prop = property_override or net.get("variant", "sort")
    m = net.get("m")
    tested = 0
    if prop == "sort":
        for v in itertools.product([0, 1], repeat=n):
            out = simulate(n, comps, v)
            if list(out) != sorted(out):
                return False, "input %s -> output %s not sorted" % (v, out), tested + 1
            tested += 1
    elif prop == "median":
        mid = n // 2
        for v in itertools.product([0, 1], repeat=n):
            out = simulate(n, comps, v)
            want = sorted(v)[mid]
            if out[mid] != want:
                return False, ("input %s -> ch%d=%d but median bit is %d"
                               % (v, mid + 1, out[mid], want)), tested + 1
            tested += 1
    elif prop == "merge":
        if not m:
            raise ValueError("merge needs --m")
        for a in range(m + 1):
            left = tuple(1 if k >= m - a else 0 for k in range(m))
            for b in range(n - m + 1):
                right = tuple(1 if k >= (n - m) - b else 0 for k in range(n - m))
                v = left + right
                out = simulate(n, comps, v)
                if list(out) != sorted(out):
                    return False, "input %s -> output %s not merged" % (v, out), tested + 1
                tested += 1
    else:
        raise ValueError("unknown property " + prop)
    return True, "all %d inputs pass" % tested, tested


KNOWN = {
    # validated anchors: (variant, n[, m]) -> minimal size (see S5-novelty.md)
}


def selftest():
    ok = True

    def expect(net, prop, should_pass, label):
        nonlocal ok
        passed, msg, _ = check(net, prop)
        good = (passed == should_pass) or (not should_pass and not passed)
        if passed != should_pass:
            print("SELFTEST FAIL: %s (%s)" % (label, msg))
            ok = False
        else:
            print("SELFTEST OK  : %s" % label)

    # known-good optimal sorting networks
    expect({"n": 4, "comparators": [[1, 2], [3, 4], [1, 3], [2, 4], [2, 3]]},
           "sort", True, "optimal s(4)=5 network sorts")
    # broken: drop last comparator of the above -> must fail on input 0011? test
    expect({"n": 4, "comparators": [[1, 2], [3, 4], [1, 3], [2, 4]]},
           "sort", False, "truncated s(4) network rejected")
    # identity network fails to sort
    expect({"n": 2, "comparators": []}, "sort", False, "empty n=2 network rejected")
    # wrong-shape network: comparator (2,1) style invalid handled by simulate? we
    # forbid i<j violation via ValueError; emulate with swapped roles instead
    expect({"n": 3, "comparators": [[1, 2], [1, 2], [1, 2]]},
           "sort", False, "redundant-only network rejected")
    # median: 3-sorter's middle output is the median of 3
    expect({"n": 3, "comparators": [[1, 2], [2, 3], [1, 2]]},
           "median", True, "3-comparator median-of-3 works")
    expect({"n": 3, "comparators": [[1, 2]]},
           "median", False, "1-comparator median-of-3 rejected")
    # merge: Batcher odd-even merge of two sorted pairs, size 4
    expect({"n": 4, "m": 2,
            "comparators": [[1, 3], [2, 4], [2, 3]]},
           "merge", True, "odd-even merge(2,2) size 4 works")
    expect({"n": 4, "m": 2,
            "comparators": [[1, 3], [2, 4], [1, 2]]},
           "merge", False, "truncated merge(2,2) rejected")
    print("SELFTEST %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net")
    ap.add_argument("--property", dest="prop")
    ap.add_argument("--m", type=int)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        sys.exit(selftest())
    with open(args.net) as f:
        net = json.load(f)
    if args.m:
        net["m"] = args.m
    passed, msg, cnt = check(net, args.prop)
    print("%s: %s | n=%d comps=%d inputs_tested=%d" %
          ("PASS" if passed else "FAIL", msg, net["n"], len(net["comparators"]), cnt))
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
