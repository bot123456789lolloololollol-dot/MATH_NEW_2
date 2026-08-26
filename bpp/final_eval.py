"""Final comprehensive evaluation: RSBF vs all baselines, all families.

Produces experiments/RSBF_results.json + a markdown table.
Protocol:
  - tuning seed stream draws (555) are used ONLY for choosing
    Harmonic m*, TBF tau*, RSBF hyper-params;
  - evaluation uses fresh seed 20260826 draws, N_EVAL streams per family;
  - every comparison is paired per stream; Wilcoxon signed-rank on nonzero
    diffs, Cohen's dz, bootstrap 95% CI of the mean difference.
"""
import json
import numpy as np

from bpp.online import online_nf, online_ff, online_bf, online_wf, \
    harmonic_count
from bpp.tbf import tbf, pfb_bf, tune_tbf, tune_harmonic
from bpp.rsbf import rsbf, tune_rsbf
from bpp.streams import STREAM_GENS, HELDOUT_STREAMS
from bpp.stats_tools import paired_report

TUNE_SEED = 555
EVAL_SEED = 20260826
N_TUNE = 12
N_EVAL = 300


def main():
    fams = list(STREAM_GENS.keys())
    heldout = list(HELDOUT_STREAMS.keys())

    # ---------------- tuning (never reused for eval draws)
    trng = np.random.default_rng(TUNE_SEED)
    train = []
    for nm in fams:
        for _ in range(N_TUNE):
            s, c = STREAM_GENS[nm](trng)
            train.append((np.asarray(s, dtype=np.int64), int(c)))
    m_star, _ = tune_harmonic(train)
    tau_star, _ = tune_tbf(train)
    (w, th, e, t), _ = tune_rsbf(
        [(s, c) for s, c in train],
        warmups=(20, 50, 100), thetas=(0.80, 0.90, 0.97),
        epss=(0.05, 0.15), taus=(0.02, 0.05))
    print(f"tuned: m*={m_star} tau*={tau_star} rsbf=(warmup {w}, theta {th},"
          f" eps {e}, tau {t})", flush=True)

    # ---------------- evaluation
    erng = np.random.default_rng(EVAL_SEED)

    def draw(gen):
        s, c = gen(erng)
        return np.asarray(s, dtype=np.int64), int(c)

    table = {}
    for nm in fams + heldout:
        gen = HELDOUT_STREAMS[nm] if nm in heldout else STREAM_GENS[nm]
        rows = {k: [] for k in
                ("rsbf", "bf", "ff", "wf", "nf", "hg", "tbf", "pfb", "l1")}
        for _ in range(N_EVAL):
            sa, cap = draw(gen)
            l1 = max(1, -(-int(sa.sum()) // cap))
            rows["rsbf"].append(rsbf(sa, cap, w, th, e,
                                     int(round(t * 1000)), 1000))
            rows["bf"].append(int(online_bf(sa, cap)))
            rows["ff"].append(int(online_ff(sa, cap)))
            rows["wf"].append(int(online_wf(sa, cap)))
            rows["nf"].append(int(online_nf(sa, cap)))
            rows["hg"].append(int(harmonic_count(sa, cap, m_star)))
            rows["tbf"].append(int(tbf(sa, cap, tau_star)))
            rows["pfb"].append(int(pfb_bf(sa, cap)))
            rows["l1"].append(l1)
        arrs = {k: np.array(v, dtype=float) for k, v in rows.items()}
        entry = {"totals": {k: float(arrs[k].sum()) for k in arrs},
                 "excess_vs_L1": {k: float(np.mean((arrs[k] - arrs["l1"])
                                                   / arrs["l1"]))
                                  for k in arrs if k != "l1"}}
        comp = {}
        for base in ("bf", "tbf", "ff", "hg"):
            rep = paired_report(arrs["rsbf"], arrs[base], "rsbf", base)
            rep.pop("boot_ci95", None)
            comp[base] = rep
        entry["paired"] = comp
        table[nm] = entry

    out = {"config": {"m_star": int(m_star), "tau_star": float(tau_star),
                      "rsbf": {"warmup": int(w), "theta_mid": float(th),
                               "eps_small": float(e), "tau": float(t)},
                      "n_eval_per_family": N_EVAL,
                      "eval_seed": EVAL_SEED, "tune_seed": TUNE_SEED},
           "results": table}
    with open("experiments/RSBF_results.json", "w") as fh:
        json.dump(out, fh, indent=2, default=float)

    # markdown summary
    lines = ["| family | RSBF | BF | TBF | FF | HGr(m*) | RSBF vs BF "
             "(W/T/L, p) |",
             "|---|---|---|---|---|---|---|"]
    for nm, e in table.items():
        t = e["totals"]
        p = e["paired"]["bf"]["wilcoxon_p"]
        ps = f"{p:.1e}" if p is not None else "ties"
        lines.append(
            f"| {nm} | {t['rsbf']:.0f} | {t['bf']:.0f} | {t['tbf']:.0f} | "
            f"{t['ff']:.0f} | {t['hg']:.0f} | "
            f"{e['paired']['bf']['wins_a']}/{e['paired']['bf']['ties']}/"
            f"{e['paired']['bf']['wins_b']}, {ps} |")
    md = "\n".join(lines)
    with open("experiments/RSBF_table.md", "w") as fh:
        fh.write(md + "\n")
    print(md)


if __name__ == "__main__":
    main()
