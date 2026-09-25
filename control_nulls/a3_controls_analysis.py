#!/usr/bin/env python3
"""A3 negative control: control samples through the eSS clustering pipeline.

Compares the control-derived clusters (consensus profiles) against the 49 real eSS
and against COSMIC, to test whether untreated samples yield eSS-like clusters that
would pass the same criteria.
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
SPA = os.path.join(os.path.dirname(HERE), "A2_SPA_decompose_null")
CTRL_DIR = os.path.join(HERE, "results_controls_filtered")   # >=307-count filter (matches eSS criteria); unfiltered gives identical flags

ess = pd.read_csv(os.path.join(SPA, "input_eSS_reference.txt"), sep="\t").set_index("MutationType")
ctrl = pd.read_csv(os.path.join(CTRL_DIR, "denovo_clusters_0.85.tsv"), sep="\t").set_index("MutationType").loc[ess.index]
counts = pd.read_csv(os.path.join(CTRL_DIR, "main_clusters_mutation_counts.tsv"), sep="\t")

E = (ess / ess.sum()).to_numpy(float)
C = (ctrl / ctrl.sum()).to_numpy(float)
ess_names, ctrl_names = list(ess.columns), list(ctrl.columns)


def cos(a, b):
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))


def source_label(samples_field):
    first = samples_field.split("|")[0]
    return "_".join(first.split("_")[:2])


M = np.array([[cos(C[:, i], E[:, j]) for j in range(E.shape[1])] for i in range(C.shape[1])])

# the pipeline labels control clusters "eSSn"; rename to control_cluster_n to avoid
# confusion with the real eSS (which stay as eSS in best_eSS / closest columns).
def ctrl_label(name):
    return "control_cluster_" + name[3:]

# control-cluster summary
rows = []
for i, cn in enumerate(ctrl_names):
    row = counts[counts.Cluster_ID == cn]
    n = int(row.N_Samples.iloc[0]) if len(row) else np.nan
    tot = int(row.Total_Mutations.iloc[0]) if len(row) else np.nan
    src = source_label(row.Samples.iloc[0]) if len(row) else ""
    j = M[i].argmax()
    rows.append(dict(control_cluster=ctrl_label(cn), n_samples=n, total_mutations=tot,
                     dominant_source=src, best_eSS=ess_names[j], cosine_to_eSS=round(M[i, j], 3)))
ctrl_summary = pd.DataFrame(rows)
ctrl_summary.to_csv(os.path.join(HERE, "A3_control_cluster_summary.csv"), index=False)

# per-eSS closest control cluster
best = M.max(0)
argbest = M.argmax(0)
ess_tab = pd.DataFrame({
    "eSS": ess_names,
    "max_cosine_to_control": np.round(best, 3),
    "closest_control_cluster": [ctrl_label(ctrl_names[argbest[j]]) for j in range(len(ess_names))],
    "flag": ["background-adjacent (>=0.90)" if b >= 0.90 else
             "similar (0.85-0.90)" if b >= 0.85 else "distinct (<0.85)" for b in best],
}).sort_values("max_cosine_to_control", ascending=False)
ess_tab.to_csv(os.path.join(HERE, "A3_eSS_vs_control_cosine.csv"), index=False)

print(f"control clusters: {len(ctrl_names)} (all >=3 samples: {(ctrl_summary.n_samples >= 3).all()})")
print(f"eSS with control cosine >=0.90: {(best>=0.90).sum()};  >=0.85: {(best>=0.85).sum()};  <0.85 (distinct): {(best<0.85).sum()}")

# figure: max control cosine per eSS
order = np.argsort(best)
colors = ["#b3202c" if best[j] >= 0.90 else "#e08214" if best[j] >= 0.85 else "#7fb069" for j in order]
fig, ax = plt.subplots(figsize=(7.5, 10))
y = np.arange(len(ess_names))
ax.barh(y, best[order], color=colors)
ax.set_yticks(y)
ax.set_yticklabels([ess_names[j] for j in order], fontsize=7)
ax.axvline(0.85, color="0.5", ls="--", lw=1)
ax.axvline(0.90, color="0.3", ls="-.", lw=1)
ax.set_xlabel("Max cosine to any control-derived cluster", fontsize=11)
ax.set_xlim(0.4, 1.0)
ax.set_title("A3: how close is each eSS to a control background cluster?", fontsize=11)
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ["#b3202c", "#e08214", "#7fb069"]]
ax.legend(handles, ["background-adjacent (>=0.90)", "similar (0.85-0.90)", "distinct (<0.85)"],
          fontsize=8, loc="lower right", frameon=False)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, f"fig_A3_eSS_vs_control.{ext}"), dpi=300, bbox_inches="tight")
print("wrote A3_control_cluster_summary.csv, A3_eSS_vs_control_cosine.csv, fig_A3_eSS_vs_control")
