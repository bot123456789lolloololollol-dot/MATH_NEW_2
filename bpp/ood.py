"""Evaluate a program vs baselines across: fresh-seed training families,
held-out families, standard suites, adversarial instances."""
import sys
import json
import numpy as np

from .core import ffd, bfd, wfd, l2_lower_bound
from .baselines import djd_full
from .instances import GENERATORS, HELDOUT_GENERATORS
from .evaluate import prep_dataset
from .suites import load_falkenauer, load_scholl_set1, load_scholl_hard, \
    load_waescher, load_hard28
from .champion import load_champion


def eval_vs_baselines(instances, bins_fn, name="H"):
    """instances: list of dict(name,sizes,cap) or (sizes,cap). Paired stats."""
    rows = []
    for inst in instances:
        if isinstance(inst, dict):
            s, cap = inst["sizes"], int(inst["cap"])
        else:
            s, cap = inst
        sd = np.sort(s)[::-1].copy()
        h = bins_fn(sd, cap)
        f = ffd(sd, cap)
        b = bfd(sd, cap)
        base = min(f, b)
        lb = l2_lower_bound(s, cap)
        rows.append({
            "name": inst.get("name", "") if isinstance(inst, dict) else "",
            "H": h, "ffd": f, "bfd": b, "best_classical": base,
            "l2": int(lb),
            "H_excess": (h - lb) / lb,
            "base_excess": (base - lb) / lb,
        })
    H = np.array([r["H"] for r in rows], dtype=np.int64)
    B = np.array([r["best_classical"] for r in rows], dtype=np.int64)
    L2 = np.array([r["l2"] for r in rows], dtype=np.int64)
    wins = int(np.sum(H < B))
    ties = int(np.sum(H == B))
    losses = int(np.sum(H > B))
    tot_H, tot_B = int(H.sum()), int(B.sum())
    summary = {
        "n": len(rows),
        "mean_excess_H": float(np.mean((H - L2) / L2)),
        "mean_excess_base": float(np.mean((B - L2) / L2)),
        "wins": wins, "ties": ties, "losses": losses,
        "total_bins_H": tot_H, "total_bins_base": tot_B,
        "bins_saved": tot_B - tot_H,
    }
    return rows, summary


def gen_instances(gen_dict, names, rng, per_family):
    out = []
    for nm in names:
        for _ in range(per_family):
            s, c = gen_dict[nm](rng)
            out.append({"name": nm, "sizes": s, "cap": c})
    return out


def main(champion_path, out_json="experiments/ood_eval.json", per_family=25,
         scholl_n_max=100):
    H = load_champion(champion_path)
    print("champion program:", H.str)
    results = {}

    rng = np.random.default_rng(999)
    train_fresh = gen_instances(GENERATORS, list(GENERATORS.keys()), rng,
                                per_family)
    _, r1 = eval_vs_baselines(train_fresh, H)
    results["train_fresh"] = r1
    print("train_fresh:", json.dumps(r1))

    heldout = gen_instances(HELDOUT_GENERATORS, list(HELDOUT_GENERATORS.keys()),
                            rng, per_family)
    _, r2 = eval_vs_baselines(heldout, H)
    results["heldout"] = r2
    print("heldout:", json.dumps(r2))

    # standard suites
    suites = {
        "falkenauer_u": [i for i in load_falkenauer(classes=("u120", "u250", "u500", "u1000"))],
        "falkenauer_t": [i for i in load_falkenauer(classes=("t60", "t120"))],
        "scholl_hard": load_scholl_hard(),
        "waescher": load_waescher(),
        "hard28": load_hard28(),
    }
    scholl = load_scholl_set1()
    suites["scholl_small"] = [i for i in scholl if len(i["sizes"]) <= scholl_n_max]
    suites["scholl_big"] = [i for i in scholl if len(i["sizes"]) > scholl_n_max]
    for k, insts in suites.items():
        _, rk = eval_vs_baselines(insts, H)
        results[k] = rk
        print(f"{k}: {json.dumps(rk)}")

    with open(out_json, "w") as fh:
        json.dump({"champion": H.str, "results": results}, fh, indent=2)
    return results


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "experiments/run1_champion.npz"
    main(path)
