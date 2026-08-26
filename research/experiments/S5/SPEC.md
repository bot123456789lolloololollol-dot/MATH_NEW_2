# S5 SPEC — SAT Encoding of Size-Minimal Comparator Networks

This document specifies the encoding precisely enough for independent reimplementation.
Code: `sat_net.py` (encoder), `solve_run.py` (solver drivers), `checker.py` (independent verifier).

## 1. Problem

Given channel count `n` and slot count `T`, decide whether there exists a comparator
network on channels `1..n` consisting of at most `T` comparators (a comparator is an
ordered pair `(i,j)`, `i<j`, that replaces the values on channels `i,j` by their min/max)
satisfying a target property, and to certify minimality by proving UNSAT at `T-1`.

Variants implemented:
- `sort`: output fully ascending for every input (`s(n)` minimal sorting network size).
- `median`: output channel `mid = n//2` (0-based; odd `n`) carries the median bit.
- `merge --m m`: channels are partitioned into first block `1..m` and second block
  `m+1..n`; each block is initially ascending; output must be fully ascending.

## 2. Variables

Let `P = n(n-1)/2`, pairs indexed `p = 0..P-1` in lexicographic order
(`p=0` is `(1,2)`). Inputs enumerated as bit vectors `v` (see §4).

| variable      | id                              | meaning                                   |
|---------------|---------------------------------|-------------------------------------------|
| `x(l,p)`      | `l*P + p + 1`                   | slot `l` holds comparator pair `p`         |
| `e(l)`        | `T*P + l + 1`                   | slot `l` is empty (no-op)                  |
| `w(v,l,c)`    | `TP+T + (vi*(T+1)+l)*n + c + 1` | value on channel `c` after slot `l` under input `v` |

`l` ranges over `0..T`; `l=0` is the input level.

## 3. Clauses

### 3.1 Slot choice (exactly one active per slot)
For each slot `l`: at-least-one over `{x(l,p)} ∪ {e(l)}`; pairwise at-most-one.

### 3.2 Simulation semantics (TOTAL — this matters for soundness)
For every input `v`, slot `l`, and channel `c`, let `cur = w(v,l,c)`, `nxt = w(v,l+1,c)`:
- For each pair `p=(i,j)` with `i == c` (comparator's min side):
  `x(l,p) → nxt = min(cur, w(v,l,j))`, encoded as
  `[¬x, ¬nxt, cur]`, `[¬x, ¬nxt, w(v,l,j)]`, `[¬x, nxt, ¬cur, ¬w(v,l,j)]`.
- For each pair `p=(i,j)` with `j == c` (max side): analogous with roles swapped:
  `[¬x, -cur, nxt]`, `[¬x, ¬w(v,l,i), nxt]`, `[¬x, cur, w(v,l,i), ¬nxt]`.
- For each pair not touching `c`: `x(l,p) → copy`: `[¬x, ¬nxt, cur]`, `[¬x, nxt, ¬cur]`.
- Empty: `e(l) → copy`: `[¬e, ¬nxt, cur]`, `[¬e, nxt, ¬cur]`.

Because exactly one of `{x(l,p)} ∪ {e(l)}` is true and the clauses above cover every
possible active choice, `w(v,l+1,c)` is always exactly constrained. A previous draft
constrained only the two compared channels and left others free — that encoding is
UNSOUND (spurious SAT); do not use it.

Input level: unit clauses `w(v,0,c) ↔ v[c]`.

### 3.3 Property
- `sort` / `merge`: forbid descents on the final level:
  `[¬w(v,T,k), w(v,T,k+1)]` for all `k` (forbids `out[k]=1 ∧ out[k+1]=0`).
- `median`: unit clause `w(v,T,mid) ↔ med(v)` where `med(v)` is the median bit of `v`.

## 4. Input enumeration (soundness-critical)

- `sort`/`median`: ALL binary inputs except the ascending-canonical ones
  `0^(n-k) 1^k`, `k=0..n`. Exclusion is sound for both variants because comparators
  preserve ascending-sortedness: a canonical input produces an ascending output, which
  is trivially sorted (`sort`) and whose middle channel is automatically its median
  (`median`), so these inputs can never violate either property.
- `merge`: only the constant inputs arising from ascending blocks: first block
  `0^(m-a) 1^a` for `a=0..m`, second block likewise — total `(m+1)(n-m+1)` inputs.
  This is EXHAUSTIVE over the merger's actual input space (inputs are sorted lists).

## 5. Symmetry breaking (each preserves satisfiability of the existence question)

- **SB1** (default ON): adjacent slots may not hold disjoint comparators in decreasing
  lex order: clause `[¬x(l,pi), ¬x(l+1,pj)]` whenever pairs `pi`, `pj` are disjoint and
  `pj < pi`. Sound because adjacent disjoint comparators commute; any network can be
  bubble-sorted into lex-normal form without changing its function or size.
- **SB2** (ON for `sort`/`merge`... see below, OFF for `median`): unit `[x(0,0)]` — the
  first comparator may be assumed `(1,2)` by channel relabeling. Sound for `sort`
  because relabeling bijects networks preserving size and the sorting property.
  **NOT sound for `median`** (the designated output channel must be fixed by the
  relabeling, which is impossible when the first comparator touches it) — we learned
  this empirically: with SB2 ON the encoder wrongly reports UNSAT at `med(5)=7`.
  For `merge`, channel blocks have roles; SB2 is disabled (not offered).
- **SB3** (default ON): no identical pair in consecutive slots: `[¬x(l,p), ¬x(l+1,p)]`.
  Sound: a repeated consecutive comparator is a no-op; removing duplicates then
  re-bubbling terminates and preserves the function with fewer comparators.
- Compatibility (SB1∧SB2∧SB3 simultaneously achievable): given any network, relabel
  (SB2), remove consecutive duplicates (SB3), bubble adjacent lex-inversions among
  disjoint neighbors (SB1), repeat. Each step preserves the computed function; SB steps
  cannot displace the lex-minimal `(1,2)` from slot 0; termination is immediate from a
  lexicographic potential (comparator count, inversion count).

## 6. Soundness of UNSAT-based lower bounds

UNSAT of the CNF at slot count `T` implies no network with ≤ T comparators satisfies
the property, PROVIDED the encoding admits every such network padded with empty slots:
(3.1)+(3.2) reproduce the exact min/max/copy semantics of any chosen comparator list;
§5 transformations map any network to a model of the SB-restricted formula with equal
or smaller size; §4 enumerates enough inputs by the arguments given (0-1 principle /
threshold argument for median / exhaustiveness of constant inputs for merge).
Hence `UNSAT(T-1) ∧ SAT(T)` proves the minimum equals `T`.

## 7. Independent verification

`checker.py` shares no code with the encoder: it reads a JSON comparator list and
brute-force simulates all `2^n` binary inputs (`sort`, `median`) or all `(m+1)(n-m+1)`
constant inputs (`merge`). Self-test mode validates the checker itself against
known-good and known-bad networks (including a case where the author's hand-written
anchor was wrong and the checker rejected it).

## 8. Known pitfalls encountered (do NOT re-introduce)

1. Inverted sortedness direction (forbidding ascents instead of descents) — produced
   global UNSAT including known-satisfiable instances.
2. Max-gate definitional clauses with wrong literal signs — caught by forcing a
   checker-validated optimal network into the encoder and diffusing every clause
   against the true trajectory (see `REPORT.md`, debugging methodology).
3. SB2 applied to the median variant — unsound, see §5.
4. Partial semantics encoding (only compared channels constrained) — unsound, see §3.2.
