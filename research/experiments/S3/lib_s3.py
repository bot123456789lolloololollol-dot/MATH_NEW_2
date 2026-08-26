"""Shared infrastructure for session S3 experiments.

Deterministic, seeded-free (no randomness needed): exhaustive generators for
unlabeled rooted trees, unlabeled free trees, and unlabeled simple graphs,
with self-validation against OEIS A000081 / A000055 / A000088 counts.

Representations
---------------
- Rooted tree: nested tuple. tree = (child1, child2, ...) with children sorted
  in increasing canonical (tuple) order. The empty tuple () is NOT a tree;
  a single vertex is ().
  We represent a rooted tree on k>=1 vertices as a tuple of its rooted
  subtree branches (each branch itself a rooted tree); sorting the branch
  tuples gives a canonical form. Size = 1 + sum(branch sizes).
- Free tree: (n, edges) with edges as frozenset of frozensets {u,v}, or an
  adjacency bitmask adjacency list. Here we keep free trees as edge lists
  built deterministically from the generator.
- Graph: (n, edges frozenset of (u,v) with u<v) plus helper conversions.
"""

from itertools import combinations_with_replacement

# ---------------------------------------------------------------- rooted trees

def gen_rooted_trees(max_n):
    """Generate all unlabeled rooted trees with sizes 1..max_n.

    Returns list `by_size` where by_size[k] is a list of canonical rooted
    trees of size k (k >= 1). A rooted tree of size k is a tuple of its
    k-1 branches (each a rooted tree tuple), sorted ascending.
    """
    # A rooted tree of size k is represented as the sorted tuple of its
    # branches (each branch itself a rooted tree). Single vertex = ():
    by_size = [None] * (max_n + 1)          # index 1..max_n used
    by_size[0] = []
    by_size[1] = [()]
    all_smaller = [[] for _ in range(max_n + 1)]   # all trees grouped by size
    all_smaller[1] = [()]
    for k in range(2, max_n + 1):
        # branches are trees of size 1..k-1, multiset, total size k-1
        # enumerate multisets in canonical order: choose branches in
        # non-decreasing (tree_tuple) order to hit each multiset once.
        pool = []                            # (size, tree) sorted by tree tuple
        for sz in range(1, k):
            pool.extend((sz, t) for t in by_size[sz])
        pool.sort(key=lambda x: x[1])
        trees = []
        stack = []                           # chosen branches (tuples)

        def rec(start_idx, remaining):
            if remaining == 0:
                trees.append(tuple(sorted(stack)))
                return
            for i in range(start_idx, len(pool)):
                sz, t = pool[i]
                # NOTE: pool is sorted by tree tuple, which does NOT respect
                # size order -> must not break early on sz > remaining.
                if sz > remaining:
                    continue
                # enforce non-decreasing tree order within the multiset
                if stack and t < stack[-1]:
                    continue
                stack.append(t)
                rec(i, remaining - sz)       # i (not i+1): repeats allowed
                stack.pop()

        rec(0, k - 1)
        trees.sort()
        by_size[k] = trees
        all_smaller[k] = by_size[k]
    return by_size


# ------------------------------------------------------------------ free trees

# OEIS A000055 (unlabeled free trees on n vertices), n=0..21
A000055 = [1, 1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159,
           7741, 19320, 48629, 123867, 317955, 823065, 2144505]


def _rtree_to_edges(tree, base=0):
    """Convert canonical rooted tree tuple to edge list with vertices
    numbered depth-first from `base`. Returns (edges, next_free_index)."""
    edges = []
    cur = base + 1
    for br in tree:
        edges.append((base, cur))
        sub, cur = _rtree_to_edges(br, cur)
        edges.extend(sub)
    return edges, cur


def rtree_size(tree):
    return 1 + sum(rtree_size(b) for b in tree)


def gen_free_trees(max_n):
    """Generate all unlabeled free trees on n=1..max_n vertices.

    Canonical construction via centroids:
      - unique-centroid trees: multisets of rooted trees (branches) with
        total size n-1 and each branch of size <= floor((n-1)/2); the root
        is the (unique) centroid.
      - bi-centroid trees (even n): two rooted trees of size n/2 joined by
        an edge between their roots, taken as an unordered pair.
    Returns dict n -> list of edge-lists (sorted tuple of (u,v)).
    """
    hmax = max(1, max_n // 2)
    rooted = gen_rooted_trees(hmax)
    out = {}
    out[1] = [()]
    if max_n >= 2:
        out[2] = [((0, 1),)]
    for n in range(3, max_n + 1):
        seen = set()
        trees = []
        half_floor = (n - 1) // 2          # max allowed branch size (case A)
        pool = []
        for sz in range(1, min(half_floor, hmax) + 1):
            pool.extend((sz, t) for t in rooted[sz])
        pool.sort(key=lambda x: x[1])
        stack = []

        def rec(start_idx, remaining):
            if remaining == 0:
                edges = []
                cur = 1
                for br in stack:
                    edges.append((0, cur))
                    sub, cur = _rtree_to_edges(br, cur)
                    edges.extend(sub)
                key = tuple(sorted(edges))
                if key not in seen:
                    seen.add(key)
                    trees.append(key)
                return
            for i in range(start_idx, len(pool)):
                sz, t = pool[i]
                # pool sorted by tuple order which ignores size: no break
                if sz > remaining or stack and t < stack[-1]:
                    continue
                stack.append(t)
                rec(i, remaining - sz)
                stack.pop()

        rec(0, n - 1)
        if n % 2 == 0 and n >= 4:
            h = n // 2
            rh = rooted[h]
            for i in range(len(rh)):
                for j in range(i + 1):
                    e1, _top1 = _rtree_to_edges(rh[i], 0)
                    e2, _ = _rtree_to_edges(rh[j], h)
                    edges = e1 + e2 + [(0, h)]
                    key = tuple(sorted(edges))
                    if key not in seen:
                        seen.add(key)
                        trees.append(key)
        trees.sort()
        out[n] = trees
    return out


def tree_adjacency(n, edges):
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


# --------------------------------------------------------------------- graphs

def _refine(n, adj, colors):
    """1-D Weisfeiler-Lehman refinement until stable. Colors are ints."""
    while True:
        sig = [(colors[v], tuple(sorted(colors[w] for w in adj[v])))
               for v in range(n)]
        uniq = {}
        new = []
        for s in sig:
            if s not in uniq:
                uniq[s] = len(uniq)
            new.append(uniq[s])
        if len(set(new)) == len(set(colors)) and new == _norm_colors(colors):
            return colors
        colors = _norm_colors(new)
        if len(set(colors)) == n:
            return colors


def _norm_colors(cols):
    m = {}
    out = []
    for c in cols:
        if c not in m:
            m[c] = len(m)
        out.append(m[c])
    return out


def canon_form(n, edges):
    """Canonical relabeling-invariant form of a simple graph on n vertices.

    Returns a string: 'n|sorted-edge-list-under-canonical-labeling'.
    Uses individualization-refinement with exhaustive branching over the
    first non-singleton cell and takes the lexicographically minimal edge
    word. Exact (no pruning that could lose the true minimum except
    standard early termination on complete orders)."""
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    deg = tuple(len(adj[v]) for v in range(n))
    best = [None]

    def edges_word(mapping):
        # mapping[v] = new label of old vertex v
        w = []
        for u, v in edges:
            a, b = mapping[u], mapping[v]
            if a > b:
                a, b = b, a
            w.append((a, b))
        w.sort()
        return w

    def rec(colors, mapping):
        colors = _refine(n, adj, list(colors))
        cells = {}
        for v, c in enumerate(colors):
            cells.setdefault(c, []).append(v)
        nonsing = [c for c in sorted(cells) if len(cells[c]) > 1]
        if not nonsing or len(set(colors)) == n:
            order = sorted(range(n), key=lambda v: colors[v])
            mp = {v: i for i, v in enumerate(order)}
            w = edges_word(mp)
            s = repr(w)
            if best[0] is None or s < best[0]:
                best[0] = s
            return
        c0 = nonsing[0]
        for v in cells[c0]:
            newc = list(colors)
            mx = max(newc) + 1
            newc[v] = mx
            rec(newc, mapping)

    rec(_norm_colors(list(deg)), {})
    return f"{n}|{best[0]}"


def _inv_bucket(n, edges):
    """Isomorphism-invariant bucket key: sorted list over vertices of
    (degree, sorted neighbor degrees, triangles-through-vertex)."""
    adjm = [[0]*n for _ in range(n)]
    for u, v in edges:
        adjm[u][v] = adjm[v][u] = 1
    sig = []
    for v in range(n):
        nb = [w for w in range(n) if adjm[v][w]]
        ndeg = tuple(sorted(len([w for w in range(n) if adjm[x][w]])
                            for x in nb))
        tri = sum(1 for i in range(len(nb)) for j in range(i+1, len(nb))
                  if adjm[nb[i]][nb[j]])
        sig.append((len(nb), ndeg, tri))
    return tuple(sorted(sig))


def gen_unlabeled_graphs(max_n, atlas=None, verbose=False):
    """All unlabeled simple graphs on n=1..max_n vertices (max_n <= 8).

    Extension scheme: every (k+1)-vertex unlabeled graph arises from some
    k-vertex unlabeled graph plus a new vertex with arbitrary neighborhood
    (delete the new vertex). Deduplicate exactly with networkx VF2 inside
    invariant buckets.
    Returns dict n -> sorted list of edge-list tuples."""
    import networkx as nx
    from networkx.algorithms.isomorphism import GraphMatcher
    if atlas is None:
        atlas = nx.graph_atlas_g()
    reps_by_n = {n: [] for n in range(1, max_n + 1)}   # list of (edges, NXGraph)
    buckets = {}                                       # (n, key) -> [indices]

    def try_add(nn, elist):
        edges = frozenset(elist)
        key = _inv_bucket(nn, elist)
        b = buckets.setdefault((nn, key), [])
        g = None
        for idx in b:
            rep_edges, rep_g = reps_by_n[nn][idx]
            if g is None:
                g = nx.Graph(elist)
                g.add_nodes_from(range(nn))
            if GraphMatcher(g, rep_g).is_isomorphic():
                return False
        if g is None:
            g = nx.Graph(elist)
            g.add_nodes_from(range(nn))
        b.append(len(reps_by_n[nn]))
        reps_by_n[nn].append((edges, g))
        return True

    # atlas contains every unlabeled graph on <= 7 vertices: seed directly
    for g in atlas:
        nn = g.number_of_nodes()
        if nn == 0 or nn > max_n or nn > 7:
            continue
        idx = {v: i for i, v in enumerate(g.nodes())}
        elist = []
        for u, v in g.edges():
            a, b = idx[u], idx[v]
            if a > b:
                a, b = b, a
            elist.append((a, b))
        try_add(nn, elist)

    # single extension step(s) beyond the atlas
    for k in range(7, max_n):
        base_graphs = [g for (_e, g) in reps_by_n[k]]
        for g in base_graphs:
            nodes = list(range(k))
            for mask in range(1 << k):     # neighborhood of new vertex k
                elist = [(u, v) for u, v in g.edges()]
                for i in range(k):
                    if mask >> i & 1:
                        elist.append((i, k))
                try_add(k + 1, elist)
        if verbose:
            print(f"  n={k+1}: {len(reps_by_n[k+1])}")
    out = {}
    for n in range(1, max_n + 1):
        out[n] = sorted(e for (e, _g) in reps_by_n[n])
    return out


# ---------------------------------------------------------------- validation

def self_test():
    A000081 = [0, 1, 1, 2, 4, 9, 20, 48, 115, 286, 719, 1842, 4766, 12486,
               32973, 87811, 235381, 634847, 1721159]
    rt = gen_rooted_trees(13)
    for n in range(1, 14):
        assert len(rt[n]) == A000081[n], \
            f"rooted trees n={n}: got {len(rt[n])}, want {A000081[n]}"

    ft = gen_free_trees(16)
    for n in range(1, 17):
        want = [1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159,
                7741, 19320][n - 1]
        assert len(ft[n]) == want, \
            f"free trees n={n}: got {len(ft[n])}, want {want}"

    # graph canon sanity: known counts up to 6 via brute force extension
    gs = gen_unlabeled_graphs(8)
    A000088 = [1, 1, 2, 4, 11, 34, 156, 1044, 12346]
    for n in range(1, 9):
        assert len(gs[n]) == A000088[n], \
            f"graphs n={n}: got {len(gs[n])}, want {A000088[n]}"
    print("lib_s3 self_test passed")


if __name__ == "__main__":
    self_test()
