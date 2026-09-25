#!/usr/bin/env python3
"""16-channel per-mutation-class null for the sub-profile 'rescue' fits.

For a signature whose reconstruction is argued within a single mutation class,
test whether the 16-channel (one-class) fit is unusual for a shuffled-eSS library.
The reference is shuffled by permuting channels within each mutation class, so a
class block keeps its 16 values but loses their arrangement.

Tested: SBS24 (C>A), SBS29 (C>A), SBS98 (C>G).

For each (signature, class):
    real_16    = cosine of the COSMIC class block vs the eSS decomposition
                 (De_Novo_map weights), restricted to that class.
    null_16(b) = best achievable 16-channel cosine (greedy k<=3 NNLS) from the
                 b-th within-class-shuffled reference.
    p          = (1 + #{null_16 >= real_16}) / (B + 1),  BH across the tested set.
"""
import os
import re
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import nnls
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
np.seterr(all="ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
ESS_REF = os.path.join(HERE, "input_eSS_reference.txt")
COSMIC = os.path.join(HERE, "input_COSMIC_Human_SBS-96_GRCh38_v3.6.tsv")
DECOMP_CSV = os.path.join(HERE, "input_De_Novo_map_to_COSMIC_SBS96.csv")

TESTS = [("SBS24", "C>A"), ("SBS29", "C>A"), ("SBS98", "C>G")]
B = 1000
SEED = 20260902

eSS = pd.read_csv(ESS_REF, sep="\t").set_index("MutationType")
cosmic = pd.read_csv(COSMIC, sep="\t").set_index("Type").loc[eSS.index]
CHANNELS = list(eSS.index)
E = (eSS / eSS.sum()).to_numpy(float)                 # 96 x 49
C = (cosmic / cosmic.sum())                            # 96 x 101 (normalised)
ESS_NAMES = list(eSS.columns)

CLASS_IDX = {}
for i, mt in enumerate(CHANNELS):
    CLASS_IDX.setdefault(mt[2:5], []).append(i)
CLASS_IDX = {k: np.array(v) for k, v in CLASS_IDX.items()}


def cosine(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else 0.0


def best_k3(y, L):
    """Best cosine of a non-negative <=3-column reconstruction of y from L (greedy)."""
    k = L.shape[1]
    chosen, best = [], 0.0
    for _ in range(3):
        gain, pick = -1.0, -1
        for i in range(k):
            if i in chosen:
                continue
            cols = chosen + [i]
            c = cosine(L[:, cols] @ nnls(L[:, cols], y)[0], y)
            if c > gain:
                gain, pick = c, i
        if pick < 0:
            break
        chosen.append(pick)
        best = max(best, gain)
    return best


def within_class_permutation(rng):
    perm = np.arange(len(CHANNELS))
    for idx in CLASS_IDX.values():
        perm[idx] = rng.permutation(idx)
    return perm


def decomposition_weights():
    """{SBS: {eSS: fraction}} from the real decomposition."""
    out = {}
    for line in open(DECOMP_CSV):
        line = line.strip()
        if not line or line.startswith("De novo"):
            continue
        fields = [f.strip() for f in line.split(",")]
        m = re.search(r"SBS([0-9]+[a-d]?)", fields[0])
        comps = re.findall(r"eSS(\d+)\s*\(([\d.]+)%\)", fields[1])
        if m and comps:
            out["SBS" + m.group(1)] = {"eSS" + n: float(p) / 100.0 for n, p in comps}
    return out


def benjamini_hochberg(pvalues):
    p = np.asarray(pvalues, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    running = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        running = min(running, p[i] * n / (rank + 1))
        q[i] = running
    return np.clip(q, 0, 1)


def main():
    weights = decomposition_weights()
    rows = []
    null_by_test = {}
    for sig, cls in TESTS:
        idx = CLASS_IDX[cls]
        y = C[sig].to_numpy(float)
        y_class = y[idx]

        w = np.zeros(len(ESS_NAMES))
        for ess, frac in weights[sig].items():
            w[ESS_NAMES.index(ess)] = frac
        reconstruction = E @ w
        real_16 = cosine(y_class, reconstruction[idx])
        real_best = best_k3(y_class, E[idx, :])

        null = np.empty(B)
        for b in range(B):
            perm = within_class_permutation(np.random.default_rng(SEED + b))
            null[b] = best_k3(y_class, E[perm, :][idx, :])
        null_by_test[(sig, cls)] = null
        p = (1 + np.sum(null >= real_16)) / (B + 1)

        rows.append(dict(cosmic_signature=sig, mutation_class=cls,
                         whole_profile_cosine=round(cosine(y, reconstruction), 3),
                         perclass_cosine_16=round(real_16, 3),
                         perclass_best_achievable_16=round(real_best, 3),
                         null_p50=round(float(np.percentile(null, 50)), 3),
                         null_p95=round(float(np.percentile(null, 95)), 3),
                         emp_p=round(float(p), 4)))

    table = pd.DataFrame(rows)
    table["emp_q"] = benjamini_hochberg(table["emp_p"]).round(4)
    table.to_csv(os.path.join(HERE, "results_A2_perclass_null.csv"), index=False)
    print(table.to_string(index=False))

    # figure: one null distribution per test with the real per-class cosine marked
    fig, axes = plt.subplots(1, len(TESTS), figsize=(4.2 * len(TESTS), 4.0), sharey=True)
    for ax, (sig, cls) in zip(axes, TESTS):
        null = null_by_test[(sig, cls)]
        real_16 = float(table.loc[(table.cosmic_signature == sig) &
                                  (table.mutation_class == cls), "perclass_cosine_16"].iloc[0])
        p = float(table.loc[(table.cosmic_signature == sig) &
                            (table.mutation_class == cls), "emp_p"].iloc[0])
        ax.hist(null, bins=np.linspace(0, 1, 40), color="#8c8cc0", alpha=0.7)
        ax.axvline(real_16, color="#b3202c", lw=2.4)
        ax.set_title(f"{sig}  ({cls})\nreal {real_16:.2f}, p={p:.3f}", fontsize=11)
        ax.set_xlabel("16-channel cosine", fontsize=10)
        ax.set_xlim(0, 1)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel(f"Shuffled eSS libraries (n={B})", fontsize=10)
    fig.suptitle("Per-mutation-class (16-channel) null", fontsize=12)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, f"fig_A2_perclass_null.{ext}"), dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
