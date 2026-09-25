#!/usr/bin/env python3
"""COSMIC v3.6 -> eSS decomposition bubble plot.

Bubble colour encodes the reconstruction cosine (0.80-1.00, with a break at 0.90);
bubble size encodes each eSS's percent contribution. Pie charts below each column
show the sample composition of that eSS. Rows are the COSMIC signatures reconstructed
at cosine >= 0.90 or significantly above the shuffled-eSS null (q < 0.05).

Row annotation is set by ANNOT:
    "cosstar"  cosine value + significance asterisk  (e.g. "0.99 **")
    "star"     significance asterisk only
    "q"        exact BH q-value

Inputs:
    DECOMP_CSV  De_Novo_map_to_COSMIC_SBS96.csv   (eSS combo + reconstruction cosine)
    Q_CSV       results_A2_SPAnull_per_signature.csv   (empirical p, BH q vs SPA null)
    PIE_TXT     cluster summary for the pie charts
"""
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
DECOMP_CSV = "/Users/mzhivagui/Documents/Papers/eSignatures_MS/MS_Split/Decompositon_COSMICv3.6_eSS/Decompose_Solution/De_Novo_map_to_COSMIC_SBS96.csv"
Q_CSV = os.path.join(HERE, "results_A2_SPAnull_per_signature.csv")
PIE_TXT = "/Users/mzhivagui/Documents/eSigantures_new_work/RESULTS/2025/SigRescueR_filter/Clustering/Poisson_resampling/eSig_update_101725_final_Cardiff_resampling_final/eSig_update_111025_FINAL/0.9_clustering_threshold/for_decomposition_cluster_summary_for_pie_charts.txt"

ANNOT = "cosstar"          # "cosstar" | "star" | "q"
OUTNAME = "decomp_bubble_cosine"
Q_SIG = 0.05
COS_RANGE = (0.80, 1.00)
NS_GREY = "#b6b6b6"

# Warm below 0.90, cool-dark at/above 0.90 (0.90 maps to the colormap midpoint).
COSINE_CMAP = LinearSegmentedColormap.from_list("cosine", [
    (0.00, "#fff7bc"), (0.30, "#fec44f"), (0.499, "#d7301f"),
    (0.50, "#5e3c99"), (0.75, "#3f007d"), (1.00, "#000000")])

# ColorBrewer Set3 order for pie slices, dominant slice first.
PIE_PALETTE = ["#8dd3c7", "#f0c419", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
               "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f"]


def significance(q):
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < Q_SIG:
        return "*"
    return "ns"


def sbs_order(name):
    m = re.match(r"SBS(\d+)([a-d]?)", name)
    return int(m.group(1)), m.group(2)


def read_decomposition(path):
    """Return {SBS: (cosine, [(eSS, percent), ...])} for eSS-decomposed signatures."""
    out = {}
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("De novo"):
            continue
        fields = [f.strip() for f in line.split(",")]
        m = re.search(r"SBS([0-9]+[a-d]?)", fields[0])
        comps = re.findall(r"eSS(\d+)\s*\(([\d.]+)%\)", fields[1])
        if not m or not comps:
            continue
        out["SBS" + m.group(1)] = (float(fields[5]),
                                   [("eSS" + n, float(p)) for n, p in comps])
    return out


def read_pies(path):
    """Return {eSS: [(sample, count), ...]}."""
    out = {}
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        fields = line.split(",")
        comp = []
        for item in fields[1:]:
            if ":" not in item:
                continue
            label, count = item.rsplit(":", 1)
            if label.strip().lower() != "total":
                comp.append((label.strip(), int(count)))
        out[fields[0].strip()] = comp
    return out


def wrap(text, width):
    parts = re.split(r"(_| )", text)
    lines, cur = [], ""
    for p in parts:
        if cur and len(cur) + len(p) > width:
            lines.append(cur)
            cur = p.lstrip()
        else:
            cur += p
    if cur:
        lines.append(cur)
    return lines or [text]


def annotation(cosine, q):
    if ANNOT == "star":
        return significance(q)
    if ANNOT == "q":
        return "q<0.001" if q < 0.001 else f"q={q:.3f}"
    return f"{cosine:.2f} {significance(q)}"


def build(signatures, cosine, pies, qmap):
    ess_list = sorted({e for s in signatures for e, _ in cosine[s][1]},
                      key=lambda e: int(e[3:]))
    xpos = {e: i for i, e in enumerate(ess_list)}
    ypos = {s: i for i, s in enumerate(signatures)}
    n_x, n_y = len(ess_list), len(signatures)
    cnorm = Normalize(*COS_RANGE)

    pie_labels, label_h = {}, 0.0
    line_h, entry_gap, wrap_w, label_fs = 0.24, 0.17, 15, 12.0
    for e in ess_list:
        rows = []
        total = sum(c for _, c in pies.get(e, [])) or 1
        for i, (sample, count) in enumerate(pies.get(e, [])):
            colour = PIE_PALETTE[i % len(PIE_PALETTE)]
            rows.append((colour, wrap(f"{sample} ({100 * count / total:.1f}%)", wrap_w)))
        pie_labels[e] = rows
        height = sum(len(w) * line_h for _, w in rows) + max(0, len(rows) - 1) * entry_gap
        label_h = max(label_h, height)
    label_h += 0.35

    bubble_h, pie_h, legend_h = max(9.0, n_y * 1.05), 5.0, 3.8
    fig_h = bubble_h + pie_h + label_h + legend_h
    strip = {"star": 0.9, "q": 2.4, "cosstar": 1.6}[ANNOT]
    fig_w = max(22, n_x * 2.45)

    fig = plt.figure(figsize=(fig_w * (n_x + strip) / n_x, fig_h))
    gs = fig.add_gridspec(4, 2, width_ratios=[n_x, strip],
                          height_ratios=[bubble_h, pie_h, label_h, legend_h],
                          hspace=0.05, wspace=0.01)
    ax = fig.add_subplot(gs[0, 0])
    sax = fig.add_subplot(gs[0, 1], sharey=ax)

    def area(pct):
        return 500 + (pct / 100.0) * 3200

    for s in signatures:
        cos = cosine[s][0]
        for e, pct in cosine[s][1]:
            ax.scatter(xpos[e], ypos[s], s=area(pct), c=[COSINE_CMAP(cnorm(cos))],
                       edgecolors="white", linewidths=1.0, alpha=0.97, zorder=3)

    ax.set_xticks(range(n_x))
    ax.set_xticklabels([])
    ax.set_yticks(range(n_y))
    ax.set_yticklabels(signatures, fontsize=26)
    ax.set_ylabel("COSMICv3.6 signature", fontsize=32, labelpad=14)
    ax.set_xlim(-0.7, n_x - 0.3)
    ax.set_ylim(-0.7, n_y - 0.3)
    ax.invert_yaxis()
    ax.tick_params(axis="y", length=0)
    ax.grid(True, color="0.88", linewidth=0.8, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for i in range(n_x):
        ax.axvline(i, color="0.94", linewidth=0.6, zorder=0)

    sax.set_xlim(0, 1)
    sax.axis("off")
    for s in signatures:
        q = qmap[s]
        centred = ANNOT == "star"
        sax.text(0.5 if centred else 0.05, ypos[s], annotation(cosine[s][0], q),
                 ha="center" if centred else "left", va="center", fontsize=22,
                 fontweight="bold", color=NS_GREY if q >= Q_SIG else "0.12")

    legend = fig.add_subplot(gs[3, :])
    legend.axis("off")
    bar = legend.inset_axes([0.07, 0.40, 0.11, 0.26])
    cb = fig.colorbar(ScalarMappable(norm=cnorm, cmap=COSINE_CMAP), cax=bar,
                      orientation="horizontal")
    cb.set_ticks([0.80, 0.85, 0.90, 0.95, 1.00])
    cb.ax.set_title("Reconstruction cosine similarity", fontsize=34, pad=16)
    cb.ax.tick_params(labelsize=28, pad=8)
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=np.sqrt(area(p)),
                      markerfacecolor="0.6", markeredgecolor="white", label=f"{p}%")
               for p in (20, 40, 60, 80, 100)]
    legend.legend(handles=handles, ncol=5, loc="center", bbox_to_anchor=(0.66, 0.42),
                  columnspacing=2.3, handletextpad=1.1, labelspacing=1.5, borderpad=1.2,
                  frameon=False, fontsize=34)
    legend.text(0.66, 0.58, "Percent contribution", ha="center", va="bottom",
                fontsize=36, transform=legend.transAxes)

    pax = fig.add_subplot(gs[1, 0])
    pax.axis("off")
    pax.set_xlim(-0.7, n_x - 0.3)
    pax.set_ylim(0, 1)
    inv = fig.transFigure.inverted()
    col_w = inv.transform(pax.transData.transform((1, 0.55)))[0] - \
            inv.transform(pax.transData.transform((0, 0.55)))[0]
    diameter = 0.9 * col_w * fig.get_size_inches()[0] / fig_h
    for e in ess_list:
        cx = inv.transform(pax.transData.transform((xpos[e], 0.55)))
        pie = fig.add_axes([cx[0] - diameter / 2, cx[1] - diameter / 2, diameter, diameter])
        pie.set_aspect("equal")
        pie.set_facecolor("none")
        comp = pies.get(e, [])
        if comp:
            pie.pie([c for _, c in comp],
                    colors=[PIE_PALETTE[i % len(PIE_PALETTE)] for i in range(len(comp))],
                    startangle=90, wedgeprops=dict(edgecolor="white", linewidth=0.8))
        pie.set_title(e, fontsize=34, fontweight="bold", pad=6)

    tax = fig.add_subplot(gs[2, 0])
    tax.axis("off")
    tax.set_xlim(-0.7, n_x - 0.3)
    tax.set_ylim(0, label_h)
    for e in ess_list:
        y = label_h - 0.18
        for colour, lines in pie_labels[e]:
            tax.scatter(xpos[e] - 0.46, y - 0.02, s=70, marker="o", facecolor=colour,
                        edgecolor="none", clip_on=False, zorder=3)
            for text in lines:
                tax.text(xpos[e] - 0.30, y, text, ha="left", va="top",
                         fontsize=label_fs, clip_on=False)
                y -= line_h
            y -= entry_gap

    fig.patch.set_facecolor("white")
    fig.savefig(os.path.join(HERE, OUTNAME + ".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(os.path.join(HERE, OUTNAME + ".png"), dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return n_x, n_y


def main():
    cosine = read_decomposition(DECOMP_CSV)
    pies = read_pies(PIE_TXT)
    pies.setdefault("eSS50", [("HFF_UVA", 1)])
    pies.setdefault("eSS51", [("Human_Intestinal_organoids_5-Fluorouracil", 1)])
    q = pd.read_csv(Q_CSV).set_index("cosmic_signature")["emp_q_SPA"]
    qmap = {s: float(q.get(s, 1.0)) for s in cosine}

    signatures = sorted([s for s, (c, _) in cosine.items()
                         if s != "SBS8" and (c >= 0.90 or qmap[s] < Q_SIG)],
                        key=sbs_order)
    n_x, n_y = build(signatures, cosine, pies, qmap)
    print(f"{OUTNAME}: {n_x} eSS x {n_y} COSMIC signatures (ANNOT={ANNOT})")


if __name__ == "__main__":
    main()
