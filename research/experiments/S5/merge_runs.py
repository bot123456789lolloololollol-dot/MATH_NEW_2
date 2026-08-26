"""Merge worker logs, verify every candidate exhaustively, and emit the
Pareto frontier vs published records.

Usage: python merge_runs.py runs_n13 runs_n14 ...
Output: verified/<tag>_frontier.json  (only networks passing ALL verifiers)
"""
import glob
import json
import os
import sys

import netlib
import verify_selector as vs
import sat_check

PUBLISHED = {   # (size, depth) points listed by Dobbelaere, fetched 2026-08-26
    13: [(45, 10), (46, 9)],
    14: [(51, 10), (52, 9)],
    15: [(56, 10), (57, 9)],
    16: [(60, 10), (61, 9)],
    17: [(71, 12), (72, 11), (74, 10)],
}


def collect(dirs):
    cands = {}
    for d in dirs:
        tag = os.path.basename(d)
        for f in glob.glob(os.path.join(d, "worker_*.jsonl")):
            for line in open(f):
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                net = tuple(map(tuple, r["net"]))
                key = (tag, net)
                m = (r["size"], r["depth"])
                if key not in cands or m < cands[key]:
                    cands[key] = m
    return cands


def main():
    dirs = sys.argv[1:] or sorted(glob.glob("runs_*"))
    cands = collect(dirs)
    print(f"{len(cands)} distinct candidates")
    by_tag = {}
    for (tag, net), m in cands.items():
        n = int(tag.replace("runs_n", "").replace("runs_sel", "")
                .split("_")[0]) if any(
                    c.isdigit() for c in tag) else None
        by_tag.setdefault(tag, []).append((m, net))

    os.makedirs("verified", exist_ok=True)
    summary = {}
    for tag, items in sorted(by_tag.items()):
        # infer n and mode from tag: runs_nNN (sort) / runs_selNN_kK (select)
        parts = tag.split("_")
        mode = "select" if tag.startswith("runs_sel") else "sort"
        n = int(parts[1][1:]) if mode == "sort" else int(parts[1][2:])
        k = int(parts[2][1:]) if mode == "select" else None
        best = {}
        for m, net in items:
            if len(net) != m[0]:
                continue                      # stale log line
            ok = False
            if mode == "sort":
                ok = (len(net) <= n * 3 and netlib.is_sorting(net, n))
            else:
                ok = vs.sel_f1(net, n, k) and vs.sel_f2(net, n, k)
            if not ok:
                print(f"REJECTED invalid candidate in {tag}: size={len(net)}")
                continue
            key = m
            if key not in best or len(net) < len(best[key]):
                best[key] = list(net)
        # Pareto frontier over (size, depth)
        pts = []
        for (s, d), net in sorted(best.items()):
            if any(s2 <= s and d2 <= d and (s2, d2) != (s, d)
                   for (s2, d2) in best):
                continue
            pts.append((s, d))
        pts.sort()
        out_name = f"verified/{tag}_frontier.json"
        payload = {"tag": tag, "mode": mode, "n": n, "k": k,
                   "frontier": [{"size": s, "depth": d, "net": best[(s, d)]}
                                for (s, d) in pts]}
        with open(out_name, "w") as f:
            json.dump(payload, f)
        pub = PUBLISHED.get(n, [])
        novel = [p for p in pts if p not in pub]
        dominated_pub = [p for p in pub if any(
            q[0] <= p[0] and q[1] <= p[1] and q != p for q in pts)]
        summary[tag] = {"points": pts, "novel_vs_published": novel,
                        "dominates_published": dominated_pub}
        print(f"{tag}: frontier={pts} novel={novel} dominates={dominated_pub}")
    with open("verified/summary.json", "w") as f:
        json.dump(summary, f, indent=1)


if __name__ == "__main__":
    main()
