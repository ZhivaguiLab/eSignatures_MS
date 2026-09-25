#!/usr/bin/env python3
"""
Organoid eSS decomposition bubble plot with per-eSS pie charts.

Each organoid sample is decomposed onto the 49 eSS references. Bubble size is the
eSS percent contribution and bubble colour is the per-sample reconstruction cosine.
Confident decompositions (cosine >= 0.90) are drawn, with a permissive cosine
>= 0.85 variant.

Inputs (in inputs/):
  Assignment_Solution_Activities.txt      per-sample eSS activities
  Assignment_Solution_Samples_Stats.txt   per-sample reconstruction cosine
  cluster_summary_for_pie_charts.txt      eSS cluster composition for the pies

Outputs (in this folder):
  organoid_decomposition_bubble_wavg.pdf / .png          cosine >= 0.90
  organoid_decomposition_bubble_wavg_cos085.pdf / .png   cosine >= 0.85
"""
import os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
INPUTS = os.path.join(HERE, "inputs")
ACT_TXT = os.path.join(INPUTS, "Assignment_Solution_Activities.txt")
STATS_TXT = os.path.join(INPUTS, "Assignment_Solution_Samples_Stats.txt")
PIE_TXT = os.path.join(INPUTS, "cluster_summary_for_pie_charts.txt")

# ------------------------------------------------------- parse sample decomp
# Weighted-average per-compound reconstruction from eSS signatures. Bubble row =
# compound, column = eSS, size = eSS percent contribution, colour = the
# compound's reconstruction cosine (from the Samples_Stats file).
def parse_decomp():
    # per-compound reconstruction cosine (col 2 = "Cosine Similarity")
    cos = {}
    with open(STATS_TXT) as fh:
        fh.readline()                      # header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                cos[f[0]] = float(f[2])
            except ValueError:
                continue
    rows = []  # (sample, ess, pct, cosine)  -- 'sample' takes the y-axis slot
    with open(ACT_TXT) as fh:
        header = fh.readline().rstrip("\n").split("\t")[1:]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sample = f[0]
            vals = [float(x) for x in f[1:]]
            tot = sum(vals)
            if tot <= 0 or sample not in cos:
                continue
            c = cos[sample]
            for ess, v in zip(header, vals):
                if v > 0:
                    rows.append((sample, ess, 100.0 * v / tot, c))
    return rows

# ---------------------------------------------------------------- parse pies
def parse_pies():
    pies = {}  # eSS -> list of (label, count)
    with open(PIE_TXT) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            fields = line.split(",")
            ess = fields[0].strip()
            comp = []
            for f in fields[1:]:
                if ":" not in f:
                    continue
                lab, cnt = f.rsplit(":", 1)
                lab = lab.strip()
                if lab.lower() == "total":
                    continue
                comp.append((lab, int(cnt)))
            pies[ess] = comp
    return pies

ROWS_ALL = parse_decomp()
# Keep only readable, confident decompositions: 0.90 <= cosine < 1.00
ROWS = [r for r in ROWS_ALL if (0.90 - 1e-9) <= r[3] < 1.0]
PIES = parse_pies()

# Pie palette: ColorBrewer "Set3" (the palette used for the cluster pie charts).
# Slices are coloured in contributor order, so the dominant slice is teal (#8dd3c7).
# Set3 with the washed-out pale yellow (#ffffb3) replaced by a clearer, less
# glaring gold so the second slice is actually visible.
SET3 = ["#8dd3c7", "#f0c419", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
        "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f"]

def pie_color(i):
    return SET3[i % len(SET3)]

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

def pie_entries(ess):
    """Exact legend entries for an eSS pie: (color, 'Sample_Name (pct.%)')."""
    comp = PIES.get(ess, [])
    total = sum(c for _, c in comp)
    out = []
    for i, (lab, cnt) in enumerate(comp):
        pct = 100.0 * cnt / total if total else 0.0
        out.append((pie_color(i), f"{lab} ({pct:.1f}%)"))
    return out

# Cosine similarity is bounded high (0.90-1.00) with "higher = better" and no
# meaningful midpoint, so use a SEQUENTIAL map (monotonic, darker = higher)
# instead of a diverging one. Truncate the light end so ~0.90 is still visible.
from matplotlib.colors import LinearSegmentedColormap
_base = plt.get_cmap("Blues")
COS_CMAP = LinearSegmentedColormap.from_list(
    "cos_seq", _base(np.linspace(0.28, 1.0, 256)))
# distinct (non-blue) sequential map for the cos>=0.90 organoid plot
_pbase = plt.get_cmap("Purples")
COS_CMAP_PURPLE = LinearSegmentedColormap.from_list(
    "cos_seq_purple", _pbase(np.linspace(0.30, 1.0, 256)))

def ess_num(e):
    return int(e[3:])

# --- sample-name parsing: "Tissue-Line_Compound_Dose" -----------------------
def sample_tissue(s):
    return s.split("-")[0].split("_")[0]          # Colon / Kidney / Liver ...

def sample_compound(s):
    parts = s.split("_")
    comp = parts[1] if len(parts) > 1 else parts[0]
    comp = comp.split("-")[0]                       # handle "AAI-25" dose style
    return "AAI" if comp.upper() in ("AA1", "AAI") else comp

def sample_sort_key(s):
    # arrange by compound first, then tissue, then the full sample name
    return (sample_compound(s).upper(), sample_tissue(s).upper(), s)

def sbs_key(s):
    m = re.match(r"SBS(\d+)([a-d]?)", s)
    return (int(m.group(1)), m.group(2))

# ---------------------------------------------------------------- bubble figure
def build_figure(rows, outname, title, vmin=0.90, cmap=COS_CMAP,
                 col_w=2.45, fs=1.0, pie_h=5.0, pie_frac=0.9, label_fs=11.5,
                 ysort=None, group_key=None):
    ess_list = sorted({r[1] for r in rows}, key=ess_num)
    # y-axis = organoid samples; order by the supplied key (compound, then tissue)
    sbs_list = sorted({r[0] for r in rows}, key=ysort)
    xpos = {e: i for i, e in enumerate(ess_list)}
    ypos = {s: i for i, s in enumerate(reversed(sbs_list))}

    n_x = len(ess_list)
    n_y = len(sbs_list)
    norm = Normalize(vmin=vmin, vmax=1.0)

    # Geometry: keep columns reasonably narrow (labels are wrapped narrow) so the
    # whole canvas isn't so wide that the fonts look tiny; give rows plenty of
    # vertical room so bubbles never touch.
    row_h = 1.05          # inches per SBS row (bubbles can't overlap at this pitch)
    fig_w = max(22, n_x * col_w)
    bubble_h = max(9.0, n_y * row_h)
    # per-pie composition labels: size the band to the busiest column
    LABEL_WRAP = 15       # wrap narrow so words push to a new line (avoids overlap)
    LABEL_FS = label_fs
    line_h = 0.23 * (label_fs / 11.5)   # inches per text line, scales with font
    entry_gap = 0.16 * (label_fs / 11.5)  # extra inches between sample entries
    label_rows = {}
    max_col_h = 0.0
    for e in ess_list:
        rws = [(c, wrap_label(t, LABEL_WRAP)) for c, t in pie_entries(e)]
        label_rows[e] = rws
        col_h = sum(len(l) * line_h for _, l in rws) + max(0, len(rws) - 1) * entry_gap
        max_col_h = max(max_col_h, col_h)
    label_h = 0.35 + max_col_h
    legend_w = 9.0        # vertical legend panel on the right
    fig_h = bubble_h + pie_h + label_h
    fig = plt.figure(figsize=(fig_w + legend_w, fig_h))
    outer = fig.add_gridspec(1, 2, width_ratios=[fig_w, legend_w], wspace=0.02)
    gs = outer[0, 0].subgridspec(3, 1,
                                 height_ratios=[bubble_h, pie_h, label_h],
                                 hspace=0.05)
    ax = fig.add_subplot(gs[0])

    # bubble area (points^2). Smaller so neighbouring rows/cols never touch.
    def area(pct):
        return 320 + (pct / 100.0) * 1900

    for sbs, ess, pct, cos in rows:
        ax.scatter(xpos[ess], ypos[sbs], s=area(pct),
                   c=[cmap(norm(cos))], edgecolors="white", linewidths=1.0,
                   alpha=0.95, zorder=3)

    ax.set_xticks(range(n_x))
    ax.set_xticklabels([])               # eSS labels appear as pie titles below
    ax.set_yticks(range(n_y))
    ax.set_yticklabels(list(reversed(sbs_list)), fontsize=18 * fs)
    ax.set_ylabel("Compound", fontsize=32 * fs, labelpad=14)
    ax.set_xlim(-0.7, n_x - 0.3)
    ax.set_ylim(-0.7, n_y - 0.3)
    ax.tick_params(axis="y", length=0)
    ax.grid(True, color="0.88", linewidth=0.8, zorder=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_title(title, fontsize=27 * fs, pad=22)
    # light vertical guides aligning bubbles to pies
    for i in range(n_x):
        ax.axvline(i, color="0.94", linewidth=0.6, zorder=0)

    # ---- group rows into boxes (e.g. by compound), with alternating shading ----
    if group_key is not None:
        from matplotlib.patches import Rectangle
        x0, x1 = -0.7, n_x - 0.3
        # contiguous runs of the same group in display order (top -> bottom)
        order = list(reversed(sbs_list))          # row index 0 is the top label
        groups = []
        start = 0
        for i in range(1, len(order) + 1):
            if i == len(order) or group_key(order[i]) != group_key(order[start]):
                groups.append((group_key(order[start]), start, i - 1))
                start = i
        for gi, (gname, i_lo, i_hi) in enumerate(groups):
            # row indices increase downward; data-y is the same as row index
            lo = min(i_lo, i_hi) - 0.5
            hi = max(i_lo, i_hi) + 0.5
            if gi % 2 == 0:
                ax.axhspan(lo, hi, xmin=0, xmax=1, color="0.92",
                           alpha=0.55, zorder=0)
            ax.add_patch(Rectangle((x0, lo), x1 - x0, hi - lo, fill=False,
                                   edgecolor="0.45", linewidth=1.6, zorder=4,
                                   clip_on=False))
            ax.text(x1 + 0.12, (lo + hi) / 2, gname, ha="left", va="center",
                    rotation=90, fontsize=20 * fs, fontweight="bold",
                    color="0.25", clip_on=False)

    # ---- vertical legend panel on the right (colorbar + size legend stacked) ----
    lax = fig.add_subplot(outer[0, 1]); lax.axis("off")
    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    # vertical colorbar near the top of the panel (sized in inches so it isn't
    # squished on a short figure)
    cb_h = min(0.34, 4.8 / fig_h)              # ~4.8in tall, capped
    cax = lax.inset_axes([0.12, 0.94 - cb_h, 0.10, cb_h])
    cb = fig.colorbar(sm, cax=cax, orientation="vertical")
    cb.ax.set_title("Cosine\nsimilarity", fontsize=24 * fs, pad=30)
    cb.ax.tick_params(labelsize=20 * fs, pad=6)

    # size legend stacked vertically below the colorbar; marker DIAMETER (points)
    # = sqrt(scatter area) so it matches the bubbles drawn in the plot.
    leg_handles = []
    for p in (20, 40, 60, 80, 100):
        leg_handles.append(Line2D([0], [0], marker="o", linestyle="",
                                  markersize=np.sqrt(area(p)),
                                  markerfacecolor="0.6", markeredgecolor="white",
                                  label=f"{p}%"))
    lax.legend(handles=leg_handles, title="Percent\ncontribution", ncol=1,
               loc="upper left", bbox_to_anchor=(0.07, 0.86 - cb_h),
               labelspacing=1.0, handletextpad=1.0, borderpad=0.8,
               frameon=False, fontsize=23 * fs, title_fontsize=25 * fs)

    # ---- pie charts band (annotations intentionally omitted for now) ----
    pax = fig.add_subplot(gs[1]); pax.axis("off")
    pax.set_xlim(-0.7, n_x - 0.3); pax.set_ylim(0, 1)
    pie_cy = 0.55                  # vertical centre of pies (data frac of pax)
    # Size pies to nearly fill the column spacing. The circle diameter equals
    # pw*fig_h (the smaller axes dimension); set that just under the per-column
    # spacing in inches so neighbouring pies almost touch without overlapping.
    spacing_in = fig_w * 0.86 / n_x
    pw = pie_frac * spacing_in / fig_h
    inv = fig.transFigure.inverted()
    for e in ess_list:
        x = xpos[e]
        comp = PIES.get(e, [])
        cx_fig = inv.transform(pax.transData.transform((x, pie_cy)))
        ia = fig.add_axes([cx_fig[0] - pw / 2, cx_fig[1] - pw / 2, pw, pw])
        ia.set_aspect("equal")
        ia.set_facecolor("none")
        cols = [pie_color(i) for i in range(len(comp))]
        if comp:
            vals = [c for _, c in comp]
            ia.pie(vals, colors=cols, startangle=90,
                   wedgeprops=dict(edgecolor="white", linewidth=0.8))
        ia.set_title(e, fontsize=24 * fs, fontweight="bold", pad=6)

    # ---- per-pie composition labels (exact Sample_Name (pct%) + colour dot) ----
    tax = fig.add_subplot(gs[2]); tax.axis("off")
    tax.set_xlim(-0.7, n_x - 0.3)
    tax.set_ylim(0, label_h)
    for e in ess_list:
        x = xpos[e]
        y = label_h - 0.18
        for col, lines in label_rows[e]:
            tax.scatter(x - 0.46, y - 0.02, s=70, marker="o",
                        facecolor=col, edgecolor="none", clip_on=False, zorder=3)
            for ln in lines:
                tax.text(x - 0.30, y, ln, ha="left", va="top",
                         fontsize=LABEL_FS, color="black", clip_on=False)
                y -= line_h            # ylim spans label_h inches, so this is 1 line
            y -= entry_gap             # extra space before the next sample entry

    fig.patch.set_alpha(0.0)
    fig.savefig(os.path.join(HERE, outname + ".pdf"), bbox_inches="tight",
                transparent=True)
    fig.savefig(os.path.join(HERE, outname + ".png"), dpi=170,
                bbox_inches="tight", transparent=True)
    plt.close(fig)
    print("wrote", outname, "  eSS:", n_x, " SBS:", n_y)

# ---------------------------------------------------------- faceted bubble figure
def build_faceted(rows, outname, title, vmin=0.90, cmap=COS_CMAP,
                  col_w=3.8, fs=1.5, pie_h=7.0, pie_frac=1.12, label_fs=16,
                  group_key=sample_compound):
    """One facet (box) per compound, sharing the same fixed eSS x-axis. The eSS
    columns line up across every facet and with the pie/label bands below."""
    ess_list = sorted({r[1] for r in rows}, key=ess_num)
    n_x = len(ess_list)
    xpos = {e: i for i, e in enumerate(ess_list)}
    norm = Normalize(vmin=vmin, vmax=1.0)

    # group samples by compound (in compound, then tissue order)
    all_samples = sorted({r[0] for r in rows}, key=sample_sort_key)
    groups = []                      # [(compound, [samples top->bottom]), ...]
    for s in all_samples:
        g = group_key(s)
        if not groups or groups[-1][0] != g:
            groups.append((g, []))
        groups[-1][1].append(s)

    # rows of bubbles per sample
    by_sample = {}
    for samp, ess, pct, cos in rows:
        by_sample.setdefault(samp, []).append((ess, pct, cos))

    row_h = 1.05
    gap_h = 0.55                     # blank spacer between facet boxes
    fig_w = max(22, n_x * col_w)
    legend_w = 9.0

    # per-pie composition labels band height (same logic as build_figure)
    LABEL_WRAP, LABEL_FS = 15, label_fs
    line_h = 0.23 * (label_fs / 11.5)
    entry_gap = 0.16 * (label_fs / 11.5)
    label_rows = {}
    max_col_h = 0.0
    for e in ess_list:
        rws = [(c, wrap_label(t, LABEL_WRAP)) for c, t in pie_entries(e)]
        label_rows[e] = rws
        max_col_h = max(max_col_h, sum(len(l) * line_h for _, l in rws)
                        + max(0, len(rws) - 1) * entry_gap)
    label_h = 0.35 + max_col_h

    # build the row height list: facet, gap, facet, gap, ..., pie, label
    facet_hs = [max(1, len(s)) * row_h for _, s in groups]
    heights, kinds = [], []          # kinds: ('facet', i) | 'gap' | 'pie' | 'label'
    for i, fh in enumerate(facet_hs):
        if i > 0:
            heights.append(gap_h); kinds.append("gap")
        heights.append(fh); kinds.append(("facet", i))
    heights += [pie_h, label_h]; kinds += ["pie", "label"]

    fig_h = sum(heights)
    fig = plt.figure(figsize=(fig_w + legend_w, fig_h))
    outer = fig.add_gridspec(1, 2, width_ratios=[fig_w, legend_w], wspace=0.02)
    gs = outer[0, 0].subgridspec(len(heights), 1, height_ratios=heights, hspace=0.0)

    def area(pct):
        return 500 + (pct / 100.0) * 3200

    # ---- facets ----
    facet_axes = {}
    for ri, kind in enumerate(kinds):
        if not (isinstance(kind, tuple) and kind[0] == "facet"):
            continue
        gi = kind[1]
        gname, samples = groups[gi]
        ax = fig.add_subplot(gs[ri]); facet_axes[gi] = ax
        ypos = {s: j for j, s in enumerate(samples)}
        for s in samples:
            for ess, pct, cos in by_sample.get(s, []):
                ax.scatter(xpos[ess], ypos[s], s=area(pct),
                           c=[cmap(norm(cos))], edgecolors="white",
                           linewidths=1.0, alpha=0.95, zorder=3)
        ax.set_xlim(-0.7, n_x - 0.3)
        ax.set_ylim(len(samples) - 0.5, -0.5)      # first sample at top
        ax.set_xticks(range(n_x)); ax.set_xticklabels([])
        ax.set_yticks(range(len(samples)))
        ax.set_yticklabels(samples, fontsize=18 * fs)
        ax.tick_params(axis="both", length=0)
        for i in range(n_x):
            ax.axvline(i, color="0.94", linewidth=0.6, zorder=0)
        ax.grid(True, axis="y", color="0.90", linewidth=0.6, zorder=0)
        # the facet itself is the "box": keep all four spines
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_edgecolor("0.45"); sp.set_linewidth(1.6)
        if gi % 2 == 1:
            ax.set_facecolor((0.0, 0.0, 0.0, 0.04))
        # compound label on the right strip
        ax.text(1.006, 0.5, gname, transform=ax.transAxes, ha="left",
                va="center", rotation=90, fontsize=20 * fs, fontweight="bold",
                color="0.25")
        if gi == 0:
            ax.set_title(title, fontsize=36 * fs, pad=22)

    # shared y-axis label
    fig.text(0.012, (label_h + pie_h + sum(facet_hs)) / 2 / fig_h,
             "Organoid sample", rotation=90, va="center", ha="center",
             fontsize=32 * fs)

    # ---- pie band ----
    pie_ri = kinds.index("pie")
    pax = fig.add_subplot(gs[pie_ri]); pax.axis("off")
    pax.set_xlim(-0.7, n_x - 0.3); pax.set_ylim(0, 1)
    spacing_in = fig_w * 0.86 / n_x
    pw = pie_frac * spacing_in / fig_h
    inv = fig.transFigure.inverted()
    for e in ess_list:
        cx_fig = inv.transform(pax.transData.transform((xpos[e], 0.55)))
        ia = fig.add_axes([cx_fig[0] - pw / 2, cx_fig[1] - pw / 2, pw, pw])
        ia.set_aspect("equal"); ia.set_facecolor("none")
        comp = PIES.get(e, [])
        if comp:
            ia.pie([c for _, c in comp],
                   colors=[pie_color(i) for i in range(len(comp))],
                   startangle=90, wedgeprops=dict(edgecolor="white", linewidth=0.8))
        ia.set_title(e, fontsize=32 * fs, fontweight="bold", pad=6)

    # ---- per-pie composition labels ----
    lab_ri = kinds.index("label")
    tax = fig.add_subplot(gs[lab_ri]); tax.axis("off")
    tax.set_xlim(-0.7, n_x - 0.3); tax.set_ylim(0, label_h)
    for e in ess_list:
        x = xpos[e]; y = label_h - 0.18
        for col, lines in label_rows[e]:
            tax.scatter(x - 0.46, y - 0.02, s=70, marker="o", facecolor=col,
                        edgecolor="none", clip_on=False, zorder=3)
            for ln in lines:
                tax.text(x - 0.30, y, ln, ha="left", va="top",
                         fontsize=LABEL_FS, color="black", clip_on=False)
                y -= line_h
            y -= entry_gap

    # ---- vertical legend panel on the right ----
    lax = fig.add_subplot(outer[0, 1]); lax.axis("off")
    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cax = lax.inset_axes([0.10, 0.80, 0.16, 0.16])
    cb = fig.colorbar(sm, cax=cax, orientation="vertical")
    cb.ax.set_title("Cosine\nsimilarity", fontsize=32 * fs, pad=16)
    cb.ax.tick_params(labelsize=27 * fs, pad=8)
    leg_handles = [Line2D([0], [0], marker="o", linestyle="",
                          markersize=np.sqrt(area(p)), markerfacecolor="0.6",
                          markeredgecolor="white", label=f"{p}%")
                   for p in (20, 40, 60, 80, 100)]
    lax.legend(handles=leg_handles, title="Percent\ncontribution", ncol=1,
               loc="upper left", bbox_to_anchor=(0.05, 0.72),
               labelspacing=2.2, handletextpad=1.4, borderpad=1.0,
               frameon=False, fontsize=30 * fs, title_fontsize=32 * fs)

    fig.patch.set_alpha(0.0)
    fig.savefig(os.path.join(HERE, outname + ".pdf"), bbox_inches="tight",
                transparent=True)
    fig.savefig(os.path.join(HERE, outname + ".png"), dpi=170,
                bbox_inches="tight", transparent=True)
    plt.close(fig)
    print("wrote", outname, "  facets:", len(groups), " eSS:", n_x)

# weighted-average per-compound decomposition into eSS (original, non-faceted plot)
build_figure(ROWS,
             "organoid_decomposition_bubble_wavg",
             "Organoid eSS decomposition",
             vmin=0.90, cmap=COS_CMAP_PURPLE, col_w=3.8, fs=1.5,
             pie_h=7.0, pie_frac=0.95, label_fs=16)

# less strict version: cos >= 0.85
ROWS_085 = [r for r in ROWS_ALL if 0.85 <= r[3] < 1.0]
build_figure(ROWS_085,
             "organoid_decomposition_bubble_wavg_cos085",
             "Organoid eSS decomposition (weighted average per compound, cos >= 0.85)",
             vmin=0.85, cmap=COS_CMAP, col_w=3.8, fs=1.5,
             pie_h=7.0, pie_frac=1.12, label_fs=16)
print("done")
