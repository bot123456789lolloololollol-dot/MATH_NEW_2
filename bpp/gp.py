"""Genetic-programming representation of *placement rules* for bin packing.

A program is a postfix instruction list evaluated once per placement step,
vectorized over all candidate bins (open bins + one VIRTUAL new bin).  Items
are processed largest-first (classic decreasing order); the program maps a
feature matrix to a score vector and we place the item into the highest-
scoring feasible candidate (ties -> lowest index; NaN scores ignored).

Subsumes Best Fit (score ~ -(resid - s)) and Worst Fit; evolution can discover
lookahead terms via remaining-size-histogram features.

Instruction encoding: int32 pairs flattened [opcode, operand].
"""
import numpy as np
from numba import njit

# ------------------------------------------------------------------ opcodes
PUSH_FEAT = 0   # operand = feature index
PUSH_CONST = 1  # operand = const table index
ADD = 2
SUB = 3
MUL = 4
DIV = 5         # protected
MIN = 6
MAX = 7
ABS = 8
NEG = 9
INV = 10        # protected reciprocal
BINOPS = (ADD, SUB, MUL, DIV, MIN, MAX)
UNOPS = (ABS, NEG, INV)
STACK_DEPTH = 12
N_CONSTS = 8

# ------------------------------------------------------------------ features
FEAT_NAMES = [
    "resid",           # 0  remaining capacity of candidate bin (cap if virtual)
    "fit",             # 1  resid - s
    "perfect",         # 2  1 if fit == 0 else 0
    "load_frac",       # 3  current load / cap
    "resid_frac",      # 4  resid / cap
    "size_frac",       # 5  item size / cap
    "progress",        # 6  fraction of items already placed
    "cnt_exact",       # 7  #remaining items whose size == fit
    "cnt_fit",         # 8  #remaining items with size <= fit
    "sum_fit_frac",    # 9  total size of remaining items <= fit, /cap
    "min_rem_frac",    # 10 smallest remaining size / cap
    "max_rem_frac",    # 11 largest remaining size / cap
    "mean_rem_frac",   # 12 mean remaining size / cap
    "nbins",           # 13 nbins / n
    "waste_frac",      # 14 (used*cap - packed) / cap
    "rem_count_frac",  # 15 #remaining / n
    "cnt_half_fit",    # 16 #remaining items with fit/2 < size <= fit
    "is_new_bin",      # 17 1 for the virtual bin
    "bin_idx_frac",    # 18 bin index / m
    "age_frac",        # 19 (step - opened_at) / n
    "dead_resid",      # 20 1 if fit>0 but NO remaining item fits the leftover
]
NF = len(FEAT_NAMES)


# ------------------------------------------------------------------ interpreter
@njit(cache=True, boundscheck=False)
def pack_gp(sizes, cap, ops, consts):
    """Pack `sizes` (sorted DESCENDING) under the program.

    Returns (nbins, assign).
    """
    n = sizes.shape[0]
    resids = np.empty(n + 1, dtype=np.int64)
    loads = np.empty(n + 1, dtype=np.int64)
    opened = np.empty(n + 1, dtype=np.int64)
    assign = np.empty(n, dtype=np.int64)

    hist = np.zeros(cap + 2, dtype=np.int64)
    pref_cnt = np.zeros(cap + 2, dtype=np.int64)
    pref_sum = np.zeros(cap + 2, dtype=np.int64)
    total_items = n
    rem_max_track = 0
    for i in range(n):
        hist[sizes[i]] += 1
        if sizes[i] > rem_max_track:
            rem_max_track = int(sizes[i])

    nb = 0
    packed = 0
    nrem = n

    # scratch buffers
    X = np.empty((NF, n + 1), dtype=np.float64)
    stack = np.empty((STACK_DEPTH, n + 1), dtype=np.float64)

    for step in range(n):
        s = int(sizes[step])
        hist[s] -= 1
        nrem -= 1
        while rem_max_track > 0 and hist[rem_max_track] == 0:
            rem_max_track -= 1

        # prefix counts/sums up to current max remaining size + min/max/sum
        rem_min = 0
        acc = 0
        cc = 0
        ss = 0
        found_min = False
        for v in range(1, rem_max_track + 1):
            c = hist[v]
            if c != 0:
                cc += c
                ss += c * v
                pref_cnt[v] = cc
                pref_sum[v] = ss
                acc += c * v
                if not found_min:
                    rem_min = v
                    found_min = True
            else:
                pref_cnt[v] = cc
                pref_sum[v] = ss
        rem_mean = acc / nrem if nrem > 0 else 0.0

        m = nb + 1
        prog_frac = step / n
        waste_frac = (nb * cap - packed) / cap

        # ---- feature matrix (col b = candidate bin; col nb = virtual)
        for b in range(m):
            if b < nb:
                r = resids[b]
                is_new = 0.0
                bidx = b / m
                age = step - opened[b]
            else:
                r = cap
                is_new = 1.0
                bidx = 1.0
                age = 0
            fit = r - s
            perfect = 1.0 if fit == 0 else 0.0
            ce = cf = sf = chf = 0.0
            if 0 <= fit <= cap:
                fi = fit if fit < rem_max_track else rem_max_track
                ce = float(hist[fit]) if fit <= cap else 0.0
                cf = float(pref_cnt[fi])
                sf = float(pref_sum[fi]) / cap
                half = (fit // 2) if (fit // 2) < rem_max_track else rem_max_track
                chf = float(pref_cnt[fi] - pref_cnt[half])
            X[0, b] = float(r)
            X[1, b] = float(fit)
            X[2, b] = perfect
            X[3, b] = loads[b] / cap
            X[4, b] = r / cap
            X[5, b] = s / cap
            X[6, b] = prog_frac
            X[7, b] = ce
            X[8, b] = cf
            X[9, b] = sf
            X[10, b] = rem_min / cap
            X[11, b] = rem_max_track / cap
            X[12, b] = rem_mean / cap
            X[13, b] = nb / max(1, total_items)
            X[14, b] = waste_frac
            X[15, b] = nrem / total_items
            X[16, b] = chf
            X[17, b] = is_new
            X[18, b] = bidx
            X[19, b] = age / max(1, total_items)
            X[20, b] = 1.0 if (fit > 0 and cf == 0.0 and nrem > 0) else 0.0

        # ---- run program once, vectorized over candidates
        sp = 0
        valid = True
        p = 0
        nops = ops.shape[0]
        while p < nops:
            op = ops[p, 0]
            arg = ops[p, 1]
            p += 1
            if sp >= STACK_DEPTH:
                valid = False
                break
            if op == PUSH_FEAT:
                for b in range(m):
                    stack[sp, b] = X[arg, b]
                sp += 1
            elif op == PUSH_CONST:
                cv = consts[arg % N_CONSTS]
                for b in range(m):
                    stack[sp, b] = cv
                sp += 1
            elif op == ADD or op == SUB or op == MUL or op == DIV \
                    or op == MIN or op == MAX:
                if sp < 2:
                    valid = False
                    break
                for b in range(m):
                    y = stack[sp - 1, b]
                    x0 = stack[sp - 2, b]
                    if op == ADD:
                        z = x0 + y
                    elif op == SUB:
                        z = x0 - y
                    elif op == MUL:
                        z = x0 * y
                    elif op == DIV:
                        d = y
                        if -1e-12 < d < 1e-12:
                            d = 1e-12
                        z = x0 / d
                    elif op == MIN:
                        z = x0 if x0 < y else y
                    else:
                        z = x0 if x0 > y else y
                    stack[sp - 2, b] = z
                sp -= 1
            else:  # unary ABS/NEG/INV
                if sp < 1:
                    valid = False
                    break
                for b in range(m):
                    y = stack[sp - 1, b]
                    if op == ABS:
                        z = abs(y)
                    elif op == NEG:
                        z = -y
                    else:
                        d = y
                        if -1e-12 < d < 1e-12:
                            d = 1e-12
                        z = 1.0 / d
                    stack[sp - 1, b] = z

        if not valid or sp < 1:
            # degenerate program: fall back to Best Fit behavior this step
            best_b = nb
            best_f = cap + 1
            for b in range(nb):
                f = resids[b] - s
                if f >= 0 and f < best_f:
                    best_f = f
                    best_b = b
        else:
            # mask infeasible real bins so they can never win the argmax
            for b in range(nb):
                if resids[b] < s:
                    stack[sp - 1, b] = np.nan
            best_b = -1
            best_score = 0.0
            have = False
            for b in range(m):
                sc = stack[sp - 1, b]
                if sc != sc:      # NaN skip
                    continue
                if not have or sc > best_score:
                    have = True
                    best_score = sc
                    best_b = b
            if best_b < 0:
                best_b = nb

        if best_b == nb:
            resids[nb] = cap - s
            loads[nb] = s
            opened[nb] = step
            nb += 1
        else:
            resids[best_b] -= s
            loads[best_b] += s
        packed += s
        assign[step] = best_b

    return nb, assign


# ------------------------------------------------------------------ program generation
class ProgSpec:
    __slots__ = ("ops", "consts")

    def __init__(self, ops, consts):
        self.ops = validate_ops(ops)
        self.consts = consts


def validate_ops(ops):
    """Drop instructions that would underflow the stack."""
    seq = []
    sp = 0
    for k in range(ops.shape[0]):
        op, arg = int(ops[k, 0]), int(ops[k, 1])
        if op in (PUSH_FEAT, PUSH_CONST):
            sp += 1
            seq.append((op, arg))
        elif op in BINOPS:
            if sp >= 2:
                sp -= 1
                seq.append((op, arg))
        elif op in UNOPS:
            if sp >= 1:
                seq.append((op, arg))
    if not seq:
        seq = [(PUSH_FEAT, 1)]
    return np.array(seq, dtype=np.int32).reshape(-1, 2)


def rand_prog(rng, max_len=28, depth=(2, 5)):
    body = []

    def leaf():
        if rng.random() < 0.85:
            body.extend([PUSH_FEAT, int(rng.integers(0, NF))])
        else:
            body.extend([PUSH_CONST, int(rng.integers(0, N_CONSTS))])

    def grow(d):
        if d <= 0 or rng.random() < 0.35:
            leaf()
            return
        if rng.random() < 0.75:
            grow(d - 1)
            grow(d - 1)
            body.extend([BINOPS[int(rng.integers(0, len(BINOPS)))], 0])
        else:
            grow(d - 1)
            body.extend([UNOPS[int(rng.integers(0, len(UNOPS)))], 0])

    grow(int(rng.integers(depth[0], depth[1] + 1)))
    if not body:
        body = [PUSH_FEAT, 1]
    return np.array(body[: 2 * max_len], dtype=np.int32).reshape(-1, 2)


CONST_POOL = np.array([0.0, 0.25, 0.5, 1.0, -1.0, 2.0, 0.1, 10.0])


def rand_consts(rng):
    idx = rng.integers(0, len(CONST_POOL), size=N_CONSTS)
    return CONST_POOL[idx].astype(np.float64)


def mutate(prog: ProgSpec, rng, rate=0.15):
    ops = prog.ops.copy()
    k = ops.shape[0]
    n_mut = max(1, int(round(rate * k)))
    for _ in range(n_mut):
        i = int(rng.integers(0, k))
        kind = rng.random()
        if kind < 0.45:                      # tweak operand
            if ops[i, 0] == PUSH_FEAT:
                ops[i, 1] = int(rng.integers(0, NF))
            elif ops[i, 0] == PUSH_CONST:
                ops[i, 1] = int(rng.integers(0, N_CONSTS))
            # operator operand unused; leave
        elif kind < 0.8:                     # swap operator
            r = rng.random()
            if r < 0.6:
                ops[i, 0] = BINOPS[int(rng.integers(0, len(BINOPS)))]
            else:
                ops[i, 0] = UNOPS[int(rng.integers(0, len(UNOPS)))]
        else:                                # insert random instr
            ins = np.array([[PUSH_FEAT, int(rng.integers(0, NF))]], dtype=np.int32)
            ops = np.insert(ops, i, ins, axis=0)
    consts = prog.consts.copy()
    if rng.random() < 0.3:
        j = int(rng.integers(0, N_CONSTS))
        consts[j] = CONST_POOL[int(rng.integers(0, len(CONST_POOL)))] + \
            float(rng.normal() * 0.05)
    return ProgSpec(ops, consts)


def crossover(a: ProgSpec, b: ProgSpec, rng):
    ka, kb = a.ops.shape[0], b.ops.shape[0]
    ia = int(rng.integers(0, ka))
    ib = int(rng.integers(0, kb)) if kb > 0 else 0
    ja = int(rng.integers(ia + 1, ka + 1))
    jb = int(rng.integers(ib + 1 if kb > 0 else 0, kb + 1))
    child = np.vstack([a.ops[:ia], b.ops[ib:jb], a.ops[ja:]])
    consts = a.consts.copy()
    if rng.random() < 0.5:
        take = int(rng.integers(0, N_CONSTS))
        consts[take:] = b.consts[take:]
    return ProgSpec(child, consts)


def program_to_str(prog: ProgSpec):
    """Pretty-print postfix to infix-ish expression."""
    st = []
    for k in range(prog.ops.shape[0]):
        op, arg = int(prog.ops[k, 0]), int(prog.ops[k, 1])
        if op == PUSH_FEAT:
            st.append(FEAT_NAMES[arg])
        elif op == PUSH_CONST:
            st.append(repr(float(prog.consts[arg])))
        elif op in BINOPS:
            b_ = st.pop()
            a_ = st.pop()
            sym = {ADD: "+", SUB: "-", MUL: "*", DIV: "/", MIN: "min", MAX: "max"}[op]
            if sym in ("min", "max"):
                st.append(f"{sym}({a_},{b_})")
            else:
                st.append(f"({a_}{sym}{b_})")
        else:
            a_ = st.pop()
            sym = {ABS: "abs", NEG: "-", INV: "inv"}[op]
            st.append(f"{sym}({a_})")
    return st[-1] if st else "<empty>"
