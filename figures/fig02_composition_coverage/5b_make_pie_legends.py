#!/usr/bin/env python3
"""
Aligned legend strip for the pie charts in the decomposition bubble plots.

Recreates the exact "Sample_Name (pct%)" labelling (with the Set3 colour dot)
that appears next to each pie in the cluster Summary figure, laid out in columns
that line up with the eSS pie columns of the bubble plot. Drop this underneath
the (label-less) pies of the bubble plot.

Outputs:
  - pie_legends_with_SBS5_SBS40.pdf / .png    (28 eSS columns)
  - pie_legends_without_SBS5_SBS40.pdf / .png (26 eSS columns)
"""
import os, re, textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DECOMP_CSV = os.path.join(os.path.dirname(HERE), "De_Novo_map_to_COSMIC_SBS96.csv")
PIE_TXT = os.environ.get(
    "ESS_PIE_TXT",
    os.path.join(os.path.dirname(HERE), "cluster_summary_for_pie_charts.txt"))

SET3 = ["#8dd3c7", "#f0c419", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
        "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f"]
def pie_color(i):
    return SET3[i % len(SET3)]

EXCLUDE = {"SBS5", "SBS40", "SBS40a", "SBS40b", "SBS40c"}

def parse_decomp():
    rows = []
    with open(DECOMP_CSV) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("De novo"):
                continue
            parts = [p.strip() for p in line.split(",")]
            m = re.search(r"SBS([0-9]+[a-d]?)", parts[0])
            if not m:
                continue
            sbs = "SBS" + m.group(1)
            comps = re.findall(r"eSS(\d+)\s*\(([\d.]+)%\)", parts[1])
            if not comps:
                continue
            cos = float(parts[5])
            for num, pct in comps:
                rows.append((sbs, "eSS" + num, float(pct), cos))
    return rows

def parse_pies():
    pies = {}
    with open(PIE_TXT) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            f = line.split(",")
            ess = f[0].strip()
            comp = []
            for x in f[1:]:
                if ":" not in x:
                    continue
                lab, cnt = x.rsplit(":", 1)
                lab = lab.strip()
                if lab.lower() == "total":
                    continue
                comp.append((lab, int(cnt)))
            pies[ess] = comp
    return pies

ROWS_ALL = parse_decomp()
ROWS = [r for r in ROWS_ALL if 0.90 < r[3] < 1.0]
PIES = parse_pies()

def ess_num(e):
    return int(e[3:])

def wrap_label(txt, width):
    """Wrap, allowing breaks at underscores and spaces (sample names lack spaces)."""
    parts = re.split(r"(_| )", txt)        # keep the delimiters
    lines, cur = [], ""
    for p in parts:
        if cur and len(cur) + len(p) > width:
            lines.append(cur)
            cur = p.lstrip()
        else:
            cur += p
    if cur:
        lines.append(cur)
    return lines or [txt]

def entries(ess):
    """Exact legend entries for an eSS: (color, 'Name (pct%)')."""
    comp = PIES.get(ess, [])
    total = sum(c for _, c in comp)
    out = []
    for i, (lab, cnt) in enumerate(comp):
        pct = 100.0 * cnt / total if total else 0.0
        out.append((pie_color(i), f"{lab} ({pct:.1f}%)"))
    return out

def build(rows, outname):
    # same column ordering / geometry as the bubble plot pies
    ess_list = sorted({r[1] for r in rows}, key=ess_num)
    n_x = len(ess_list)
    xpos = {e: i for i, e in enumerate(ess_list)}

    col_w = 2.45                       # inches per column (matches bubble plot)
    fig_w = max(22, n_x * col_w)

    # wrap long sample names to roughly the column width
    WRAP = 22
    wrapped = {}
    max_lines = 0
    for e in ess_list:
        ents = entries(e)
        rows_w = []
        for col, txt in ents:
            lines = wrap_label(txt, WRAP)
            rows_w.append((col, lines))
        wrapped[e] = rows_w
        n_lines = sum(len(l) for _, l in rows_w)
        max_lines = max(max_lines, n_lines)

    line_h = 0.30                      # inches per text line
    top_pad = 0.35
    fig_h = top_pad + max_lines * line_h + 0.3

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-0.7, n_x - 0.3)
    ax.set_ylim(0, fig_h)
    ax.axis("off")

    # keep the dot + text block inside the column [x-0.5, x+0.5] so it sits
    # directly under the pie that is centred on x.
    dot_dx = -0.46                     # dot near the left edge of the column
    fs = 12.5

    for e in ess_list:
        x = xpos[e]
        y = fig_h - top_pad
        for col, lines in wrapped[e]:
            # colour dot on the first line of each entry
            ax.scatter(x + dot_dx, y - 0.02, s=130, marker="o",
                       facecolor=col, edgecolor="none", clip_on=False, zorder=3)
            for j, ln in enumerate(lines):
                ax.text(x + dot_dx + 0.18, y, ln, ha="left", va="top",
                        fontsize=fs, color="black", clip_on=False)
                y -= line_h
    fig.patch.set_alpha(0.0)
    fig.savefig(os.path.join(HERE, outname + ".pdf"), bbox_inches="tight",
                transparent=True)
    fig.savefig(os.path.join(HERE, outname + ".png"), dpi=170,
                bbox_inches="tight", transparent=True)
    print("wrote", outname, "cols:", n_x, "max_lines:", max_lines,
          "fig:", round(fig_w, 1), "x", round(fig_h, 1))

build(ROWS, "pie_legends_with_SBS5_SBS40")
build([r for r in ROWS if r[0] not in EXCLUDE], "pie_legends_without_SBS5_SBS40")
