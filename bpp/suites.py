"""Loaders for standard benchmark suites (downloaded into benchmarks/)."""
import os
import re
import zipfile
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BENCH = os.path.join(ROOT, "benchmarks")


def load_falkenauer(classes=("u120", "u250", "u500", "u1000", "t60", "t120")):
    """Returns list of dicts: name, sizes (int32), cap, best_known.

    t-class files have one-decimal sizes -> scaled x10 to integers.
    """
    file_for = {"u120": 1, "u250": 2, "u500": 3, "u1000": 4,
                "t60": 5, "t120": 6, "t249": 7, "t501": 8}
    out = []
    for cls in classes:
        path = os.path.join(BENCH, "falkenauer", f"binpack{file_for[cls]}.txt")
        with open(path) as fh:
            lines = [ln.rstrip() for ln in fh]
        p = int(lines[0].split()[0])
        i = 1
        for _ in range(p):
            name = lines[i].split()[0]
            i += 1
            cap_s, n_s, best_s = lines[i].split()
            i += 1
            cap, n, bk = float(cap_s), int(n_s), int(best_s)
            scale = 10 if "." in cap_s else 1
            sizes = []
            for k in range(n):
                sizes.append(int(round(float(lines[i + k]) * scale)))
            i += n
            cap_i = int(round(cap * scale))
            assert max(sizes) <= cap_i and min(sizes) > 0
            out.append({"name": name, "sizes": np.array(sizes, dtype=np.int32),
                        "cap": cap_i, "best_known": bk})
    return out


def load_scholl_set1(max_n=None):
    """720 instances NxCyWz_v.BPP from bin1data.zip."""
    zf = zipfile.ZipFile(os.path.join(BENCH, "scholl", "bin1data.zip"))
    out = []
    for nm in sorted(zf.namelist()):
        m = re.match(r"N(\d)C(\d)W(\d)_([A-T])\.BPP", nm)
        if not m:
            continue
        data = zf.read(nm).decode().split()
        n, cap = int(data[0]), int(data[1])
        if max_n and n > max_n:
            continue
        sizes = np.array([int(x) for x in data[2:2 + n]], dtype=np.int32)
        out.append({"name": nm[:-4], "sizes": sizes, "cap": cap})
    return out


def load_scholl_hard():
    zf = zipfile.ZipFile(os.path.join(BENCH, "scholl", "bin3data.zip"))
    out = []
    for nm in sorted(zf.namelist(), key=lambda s: (len(s), s)):
        data = zf.read(nm).decode().split()
        n, cap = int(data[0]), int(data[1])
        sizes = np.array([int(x) for x in data[2:2 + n]], dtype=np.int32)
        out.append({"name": nm[:-4], "sizes": sizes, "cap": cap})
    return out


def load_waescher(which="both"):
    files = {"gau1": "WAE_GAU1.BPP", "gau2": "WAE_GAU2.BPP"}
    keys = list(files.keys()) if which == "both" else [which]
    out = []
    for k in keys:
        path = os.path.join(BENCH, "waescher", files[k])
        txt = open(path).read()
        blocks = [b for b in txt.split("'") if b.strip()]
        # split on quoted names: 'NAME' followed by body until next quote
        entries = re.split(r"'([^']*)'", txt)
        # entries: ['', name1, body1, name2, body2, ...]
        for j in range(1, len(entries), 2):
            name = entries[j].strip()
            nums = [int(x) for x in entries[j + 1].split()]
            distinct, cap = nums[0], nums[1]
            pairs = nums[2:]
            sizes = []
            for a in range(0, len(pairs), 2):
                w, c = pairs[a], pairs[a + 1]
                sizes.extend([w] * c)
            assert sum(1 for _ in sizes) >= 1
            out.append({"name": f"WAE_{k}_{name}", "sizes": np.array(sizes, dtype=np.int32),
                        "cap": cap})
    return out


def load_hard28():
    txt = open(os.path.join(BENCH, "hard28", "hard28")).read()
    entries = re.split(r"'([^']*)'", txt)
    out = []
    for j in range(1, len(entries), 2):
        name = entries[j].strip()
        nums = [int(x) for x in entries[j + 1].split()]
        distinct, cap = nums[0], nums[1]   # format: distinct count? verify below
        pairs = nums[2:]
        sizes = []
        ok = False
        # try interpretation A: first two numbers are (distinct, cap)
        w, c = pairs[0], pairs[1]
        if w <= cap:
            for a in range(0, len(pairs), 2):
                ww, cc = pairs[a], pairs[a + 1]
                sizes.extend([ww] * cc)
            ok = True
        assert ok and sizes
        out.append({"name": f"HARD28_{name.replace(' ', '_')}",
                    "sizes": np.array(sizes, dtype=np.int32), "cap": cap})
    return out


def suite_summary():
    rows = []
    for inst in (load_falkenauer() + load_scholl_set1() + load_scholl_hard()
                 + load_waescher() + load_hard28()):
        s = inst["sizes"]
        rows.append((inst["name"], len(s), inst["cap"],
                     round(float(np.mean(s)) / inst["cap"], 3)))
    return rows


if __name__ == "__main__":
    rows = suite_summary()
    print(f"total instances: {len(rows)}")
    for r in rows[:5]:
        print(r)
    for r in rows[-5:]:
        print(r)
