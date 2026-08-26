"""Core library for comparator (sorting) networks.

Representation: a network is a list of compare-exchange elements (i, j) with
i < j; semantics: after the element, wire i carries min(old_i, old_j) and wire
j carries max(old_i, old_j).

Exhaustive verification uses the classic bitset parallelism trick: for n wires,
each wire's value across ALL 2^n zero-one inputs is stored as one Python int of
2^n bits (bit p = wire value on input pattern p, where p itself encodes the
input: bit i of p is the value on wire i).  Under this encoding a compare-
exchange on wires (i, j) is exactly (w_i & w_j, w_i | w_j), because min/max of
monotone functions are their meet/join.  Correctness of this representation
follows from the fact that every wire function computed by a comparator network
is a monotone Boolean function of the input bits, and on {0,1} values
min = AND, max = OR.
"""
from itertools import combinations

# ---------------------------------------------------------------- wire masks


def wire_masks(n):
    """masks[i] has bit p set iff bit i of p is 1 (initial wire contents)."""
    N = 1 << n
    masks = []
    for i in range(n):
        period = 1 << (i + 1)
        unit = ((1 << (1 << i)) - 1) << (1 << i)          # ones in upper half
        reps = N // period
        # magic number = 0b0...1...1 with reps ones spaced `period` apart
        magic = ((1 << (reps * period)) - 1) // ((1 << period) - 1)
        masks.append(unit * magic)
    return masks


def expected_masks(n):
    """exp[k] has bit p set iff wire k of the ascending-sorted output on
    input p is 1, i.e. iff popcount(p) >= n-k (the k-th smallest of |p| ones
    is 1 exactly when fewer than k+1 zeros exist)."""
    N = 1 << n
    exp = [0] * n
    for ones in range(n + 1):
        m = 0
        for combo in combinations(range(n), ones):
            p = 0
            for c in combo:
                p |= 1 << c
            m |= 1 << p
        # wires k >= n-ones carry 1 (the last `ones` output positions)
        for k in range(n - ones, n):
            exp[k] |= m
    return exp


_MASK_CACHE = {}


def masks_for(n):
    if n not in _MASK_CACHE:
        _MASK_CACHE[n] = (wire_masks(n), expected_masks(n))
    return _MASK_CACHE[n]


# ---------------------------------------------------------------- evaluation


def evaluate(net, n):
    """Return list of final wire masks, or None if any intermediate check
    fails (never fails here; evaluation always succeeds)."""
    wm, _ = masks_for(n)
    w = list(wm)
    for (i, j) in net:
        a, b = w[i], w[j]
        w[i] = a & b
        w[j] = a | b
    return w


def is_sorting(net, n):
    """Exhaustive check over all 2^n zero-one inputs.  By the zero-one
    principle (Knuth, TAOCP 5.3.4, Theorem Z), this is equivalent to sorting
    correctly on all inputs from an arbitrary total order."""
    _, em = masks_for(n)
    w = evaluate(net, n)
    return all(w[k] == em[k] for k in range(n))


def error_bits(net, n):
    """Total number of (input pattern, wire) mismatches -- a gradient signal:
    0 iff the network sorts."""
    _, em = masks_for(n)
    w = evaluate(net, n)
    return sum((w[k] ^ em[k]).bit_count() for k in range(n))


def failing_patterns(net, n):
    """Number of distinct input patterns on which the output is not sorted."""
    _, em = masks_for(n)
    w = evaluate(net, n)
    bad = 0
    for k in range(n):
        bad |= w[k] ^ em[k]
    return bad.bit_count()


# ---------------------------------------------------------------- metrics


def depth(net):
    lvl = [0] * (max((max(i, j) for (i, j) in net), default=-1) + 1)
    d = 0
    for (i, j) in net:
        l = max(lvl[i], lvl[j]) + 1
        lvl[i] = lvl[j] = l
        d = max(d, l)
    return d


def width(net):
    return max((j for (_, j) in net), default=-1) + 1


# ---------------------------------------------------------------- transforms


def prune(net, n):
    """Remove redundant compare-exchange elements until fixpoint."""
    net = list(net)
    changed = True
    while changed:
        changed = False
        r = 0
        while r < len(net):
            cand = net[:r] + net[r + 1:]
            if is_sorting(cand, n):
                net = cand
                changed = True
            else:
                r += 1
    return net


def to_layers(net):
    """Greedy layer assignment -> list of layers, each a list of (i, j)."""
    lvl = {}
    layers = []
    for (i, j) in net:
        l = max(lvl.get(i, 0), lvl.get(j, 0))
        layers.append([])
        li = layers[l]
        li.append((i, j))
        lvl[i] = lvl[j] = l + 1
    while layers and not layers[-1]:
        layers.pop()
    return layers


# ---------------------------------------------------------------- io


def save(path, name, net, n, extra=None):
    import json
    obj = {"name": name, "n": n, "size": len(net),
           "depth": depth(net), "network": [list(e) for e in net]}
    if extra:
        obj.update(extra)
    with open(path, "w") as f:
        json.dump(obj, f)


def load(path):
    import json
    with open(path) as f:
        obj = json.load(f)
    return obj["n"], [tuple(e) for e in obj["network"]]
