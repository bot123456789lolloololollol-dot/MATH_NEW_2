"""Out-of-sample evaluation of an online policy champion vs baselines.

Baselines: NF/FF/BF/WF (online), Harmonic-Greedy(m*) and TBF(tau*) tuned on
TRAINING streams only, plus PFB (perfect-fit-first BF).  Metrics: paired bins,
excess vs L1, win/loss, Wilcoxon + Cohen's dz + bootstrap CI.
"""
import sys
import json
import numpy as np

from .online import online_nf, online_ff, online_bf, online_wf, \
    harmonic_count
from .tbf import tbf, pfb_bf, tune_tbf, tune_harmonic
from .streams import STREAM_GENS, HELDOUT_STREAMS, make_stream_set
from .champion import load_champion  # reuse loader? it uses offline pack_gp
from .gp import pack_gp  # noqa: F401 (keep import path stable)
from .stats_tools import paired_report


def load_online_champion(path):
    z = np.load(path)
    ops = z["ops"].astype(np.int32)
    consts = z["consts"].astype(np.float64)
    return ops, consts


def policy_bins(ops, consts, s, cap):
    from .online import pack_online
    return int(pack_online(np.asarray(s, dtype=np.int64), int(cap), ops, consts))


def eval_streams(ops, consts, streams):
    """Return dict of per-stream arrays for every policy."""
    res = {"H": [], "nf": [], "ff": [], "bf": [], "wf": [],
           "hg": [], "tbf": [], "pfb": [], "l1": []}
    meta = []
    for nm, s, cap in streams:
        sa = np.asarray(s, dtype=np.int64)
        l1 = max(1, -(-int(sa.sum()) // int(cap)))
        res["H"].append(policy_bins(ops, consts, sa, cap))
        res["nf"].append(int(online_nf(sa, cap)))
        res["ff"].append(int(online_ff(sa, cap)))
        res["bf"].append(int(online_bf(sa, cap)))
        res["wf"].append(int(online_wf(sa, cap)))
        res["hg"].append(int(harmonic_count(sa, cap, EVAL_M)))
        res["tbf"].append(int(tbf(sa, cap, EVAL_TAU)))
        res["pfb"].append(int(pfb_bf(sa, cap)))
        res["l1"].append(l1)
        meta.append(nm)
    return {k: np.array(v) for k, v in res.items()}, meta


EVAL_M = 8      # replaced by tuned value in main()
EVAL_TAU = 0.05


def main(champ_path, out_json="experiments/ood_online.json", n_train_tune=120):
    global EVAL_M, EVAL_TAU
    ops, consts = load_online_champion(champ_path)

    rng = np.random.default_rng(4242)
    tune_streams = [STREAM_GENS[nm](rng)[0:2] for nm in STREAM_GENS
                    for _ in range(n_train_tune // len(STREAM_GENS))]
    m_star, _ = tune_harmonic(tune_streams)
    tau_star, _ = tune_tbf([(np.asarray(s, dtype=np.int64), c)
                            for s, c in tune_streams])
    EVAL_M, EVAL_TAU = m_star, tau_star
    print(f"tuned on train: Harmonic-Greedy m*={m_star}, TBF tau*={tau_star}")

    results = {}
    sets = {}
    sets["train_fresh"] = make_stream_set(np.random.default_rng(9001), 200)
    hold_rng = np.random.default_rng(9002)
    sets["heldout"] = [(nm,) + HELDOUT_STREAMS[nm](hold_rng)
                       for nm in HELDOUT_STREAMS for _ in range(40)]
    big_rng = np.random.default_rng(9003)
    sets["large_n"] = [("uniform10k",) +
                       STREAM_GENS["uniform"](big_rng, n=10000, cap=1000),
                       ("triplet_like10k",) +
                       STREAM_GENS["triplet_like"](big_rng, n=10000, cap=1000),
                       ("discrete10k",) +
                       STREAM_GENS["discrete_3"](big_rng, n=10000, cap=750)]

    for set_name, streams in sets.items():
        arrs, meta = eval_streams(ops, consts, streams)
        block = {}
        for base in ["bf", "ff", "wf", "nf", "hg", "tbf", "pfb"]:
            rep = paired_report(arrs["H"], arrs[base], "H", base)
            rep["mean_excess_H_vs_L1"] = float(
                np.mean((arrs["H"] - arrs["l1"]) / arrs["l1"]))
            rep["mean_excess_base_vs_L1"] = float(
                np.mean((arrs[base] - arrs["l1"]) / arrs["l1"]))
            block[base] = rep
        results[set_name] = block
        b = block["bf"]
        print(f"[{set_name}] H vs BF: mean {b['mean_a']:.1f} vs "
              f"{b['mean_b']:.1f}  wins/ties/losses "
              f"{b['wins_a']}/{b['ties']}/{b['wins_b']}  "
              f"p={b['wilcoxon_p'] if b['wilcoxon_p'] is not None else 'NA'}")

    with open(out_json, "w") as fh:
        json.dump({"champion_path": champ_path, "m_star": m_star,
                   "tau_star": tau_star, "results": results}, fh, indent=2,
                  default=float)
    return results


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "experiments/onrun1_champion.npz"
    main(path)
