#!/usr/bin/env python3
"""
Recreate the COSMIC -> eSS decomposition bubble plot (with pie charts) using the
new COSMIC v3.6 decomposition.

Outputs (in this folder):
  - decomposition_bubble_with_SBS5_SBS40.pdf / .png   (full)
  - decomposition_bubble_without_SBS5_SBS40.pdf / .png
  - eSS37_focus.pdf / .png
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
DECOMP_CSV = os.path.join(os.path.dirname(HERE), "De_Novo_map_to_COSMIC_SBS96.csv")
PIE_TXT = os.environ.get(
    "ESS_PIE_TXT",
    os.path.join(os.path.dirname(HERE), "cluster_summary_for_pie_charts.txt"))

# ---------------------------------------------------------------- parse decomp
def parse_decomp():
    rows = []  # (sbs, ess, pct, cosine)
    with open(DECOMP_CSV) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("De novo"):
                continue
            parts = [p.strip() for p in line.split(",")]
            sbs_m = re.search(r"SBS([0-9]+[a-d]?)", parts[0])
            if not sbs_m:
                continue
            sbs = "SBS" + sbs_m.group(1)
            decomp = parts[1]
            # only keep rows that actually decompose into eSS signatures
            comps = re.findall(r"eSS(\d+)\s*\(([\d.]+)%\)", decomp)
            if not comps:
                continue
            cosine = float(parts[5])
            for ess_num, pct in comps:
                rows.append((sbs, "eSS" + ess_num, float(pct), cosine))
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
# Keep only readable, confident decompositions: cosine >= 0.90.
# Upper bound is inclusive (<= 1.0) so perfect reconstructions are kept, e.g.
# SBS90 -> eSS35 (100%) at cosine 1.00, which a strict "< 1.0" cutoff silently dropped.
# SBS8 is explicitly excluded (its cosine sits below the 0.90 threshold).
ROWS = [r for r in ROWS_ALL
        if (0.90 - 1e-9) <= r[3] <= 1.0 and r[0] != "SBS8"]
PIES = parse_pies()

# eSS50 & eSS51 are singleton clusters added after the 49-eSS pie-summary file was
# written (see eSS_VCFs/_FINAL_REPORT.md UPDATE 8), so they are absent from PIE_TXT
# and would otherwise render as an empty axes box. Each is one contributing sample
# (100%): eSS50 = HFF UVA, eSS51 = human intestinal-organoid 5-Fluorouracil.
PIES.setdefault("eSS50", [("HFF_UVA", 1)])
PIES.setdefault("eSS51", [("Human_Intestinal_organoids_5-Fluorouracil", 1)])

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

def ess_num(e):
    return int(e[3:])

def sbs_key(s):
    m = re.match(r"SBS(\d+)([a-d]?)", s)
    return (int(m.group(1)), m.group(2))

# ---------------------------------------------------------------- bubble figure
def build_figure(rows, outname, title):
    ess_list = sorted({r[1] for r in rows}, key=ess_num)
    sbs_list = sorted({r[0] for r in rows}, key=sbs_key)
    xpos = {e: i for i, e in enumerate(ess_list)}
    ypos = {s: i for i, s in enumerate(sbs_list)}

    n_x = len(ess_list)
    n_y = len(sbs_list)
    norm = Normalize(vmin=0.90, vmax=1.0)
    cmap = COS_CMAP

    # Geometry: keep columns reasonably narrow (labels are wrapped narrow) so the
    # whole canvas isn't so wide that the fonts look tiny; give rows plenty of
    # vertical room so bubbles never touch.
    col_w = 2.45          # inches per eSS column
    row_h = 1.05          # inches per SBS row (bubbles can't overlap at this pitch)
    fig_w = max(22, n_x * col_w)
    bubble_h = max(9.0, n_y * row_h)
    pie_h = 5.0           # pie band
    # per-pie composition labels: size the band to the busiest column
    LABEL_WRAP = 15       # wrap narrow so words push to a new line (avoids overlap)
    LABEL_FS = 11.5
    line_h = 0.23         # inches per text line within one label
    entry_gap = 0.16      # extra inches between separate sample entries
    label_rows = {}
    max_col_h = 0.0
    for e in ess_list:
        rws = [(c, wrap_label(t, LABEL_WRAP)) for c, t in pie_entries(e)]
        label_rows[e] = rws
        col_h = sum(len(l) * line_h for _, l in rws) + max(0, len(rws) - 1) * entry_gap
        max_col_h = max(max_col_h, col_h)
    label_h = 0.35 + max_col_h
    legend_h = 2.6        # horizontal legend band at the very bottom
    fig_h = bubble_h + pie_h + label_h + legend_h
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(4, 1,
                          height_ratios=[bubble_h, pie_h, label_h, legend_h],
                          hspace=0.05)
    ax = fig.add_subplot(gs[0])

    # bubble area (points^2). Max diameter ~0.85in < row pitch 1.05in => no overlap.
    def area(pct):
        return 500 + (pct / 100.0) * 3200

    for sbs, ess, pct, cos in rows:
        ax.scatter(xpos[ess], ypos[sbs], s=area(pct),
                   c=[cmap(norm(cos))], edgecolors="white", linewidths=1.0,
                   alpha=0.95, zorder=3)

    ax.set_xticks(range(n_x))
    ax.set_xticklabels([])               # eSS labels appear as pie titles below
    ax.set_yticks(range(n_y))
    ax.set_yticklabels(sbs_list, fontsize=28)
    ax.set_ylabel("COSMIC v3.6 signature", fontsize=32, labelpad=14)
    ax.set_xlim(-0.7, n_x - 0.3)
    ax.set_ylim(-0.7, n_y - 0.3)
    ax.tick_params(axis="y", length=0)
    ax.grid(True, color="0.88", linewidth=0.8, zorder=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_title(title, fontsize=36, pad=22)
    # light vertical guides aligning bubbles to pies
    for i in range(n_x):
        ax.axvline(i, color="0.94", linewidth=0.6, zorder=0)

    # ---- horizontal legend band at the bottom (colorbar + size legend) ----
    lax = fig.add_subplot(gs[3]); lax.axis("off")
    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    # horizontal colorbar on the left
    cax = lax.inset_axes([0.08, 0.42, 0.40, 0.17])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.ax.set_title("Cosine similarity", fontsize=34, pad=18) #cb.set_label("Cosine similarity", fontsize=34, labelpad=12)
    cb.ax.tick_params(labelsize=29, pad=15)

    # size legend as a horizontal row on the right; marker DIAMETER (points) =
    # sqrt(scatter area) so it matches the bubbles drawn in the plot.
    leg_handles = []
    for p in (20, 40, 60, 80, 100):
        leg_handles.append(Line2D([0], [0], marker="o", linestyle="",
                                  markersize=np.sqrt(area(p)),
                                  markerfacecolor="0.6", markeredgecolor="white",
                                  label=f"{p}%"))
    lax.legend(handles=leg_handles, title="Percent contribution", ncol=5,
               loc="center", bbox_to_anchor=(0.76, 0.5),
               columnspacing=3.6, handletextpad=1.4, borderpad=1.0,
               frameon=False, fontsize=31, title_fontsize=34)

    # ---- pie charts band (annotations intentionally omitted for now) ----
    pax = fig.add_subplot(gs[1]); pax.axis("off")
    pax.set_xlim(-0.7, n_x - 0.3); pax.set_ylim(0, 1)
    pie_cy = 0.55                  # vertical centre of pies (data frac of pax)
    # Size pies to nearly fill the column spacing. The circle diameter equals
    # pw*fig_h (the smaller axes dimension); set that just under the per-column
    # spacing in inches so neighbouring pies almost touch without overlapping.
    spacing_in = fig_w * 0.86 / n_x
    pw = 0.9 * spacing_in / fig_h
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
        ia.set_title(e, fontsize=32, fontweight="bold", pad=6)

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

# full version
_TITLE = "COSMIC v3.6 SBS decomposition into eSS signatures"
build_figure(ROWS,
             "decomposition_bubble_with_SBS5_SBS40",
             _TITLE)

# without SBS5 and SBS40 (all SBS40 variants)
EXCLUDE = {"SBS5", "SBS40", "SBS40a", "SBS40b", "SBS40c"}
rows_wo = [r for r in ROWS if r[0] not in EXCLUDE]
build_figure(rows_wo,
             "decomposition_bubble_without_SBS5_SBS40",
             _TITLE)

# ------------------------------------------------------ single-eSS focus figure
def primary_name(ess):
    """Largest-contributing sample name for an eSS (used as a descriptive title)."""
    comp = PIES.get(ess, [])
    return max(comp, key=lambda lc: lc[1])[0] if comp else ess

def build_ess_focus(target, outname):
    sub = [r for r in ROWS_ALL if r[1] == target]
    sub.sort(key=lambda r: sbs_key(r[0]))
    sbs_labels = [r[0] for r in sub]
    pcts = [r[2] for r in sub]
    coss = [r[3] for r in sub]

    norm = Normalize(vmin=min(min(coss), 0.80), vmax=1.0)
    cmap = COS_CMAP

    fig, (ax, axp) = plt.subplots(1, 2, figsize=(16, 7.5),
                                  gridspec_kw={"width_ratios": [2.4, 1]})
    y = np.arange(len(sub))
    ax.barh(y, pcts, color=[cmap(norm(c)) for c in coss],
            edgecolor="0.3", linewidth=0.6)
    ax.set_yticks(y); ax.set_yticklabels(sbs_labels, fontsize=16)
    ax.invert_yaxis()
    ax.set_xlabel(f"{target} percent contribution to COSMIC signature (%)", fontsize=16)
    ax.tick_params(axis="x", labelsize=13)
    ax.set_title(f"{target} ({primary_name(target)})\nacross COSMIC v3.6 signatures",
                 fontsize=18)
    for yi, p, c in zip(y, pcts, coss):
        ax.text(p + 1, yi, f"{p:.0f}%  (cos {c:.2f})", va="center", fontsize=13)
    ax.set_xlim(0, max(pcts) * 1.30)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("Cosine similarity of decomposition", fontsize=14)
    cb.ax.tick_params(labelsize=12)

    comp = PIES.get(target, [])
    axp.axis("off")
    if comp:
        vals = [c for _, c in comp]; labs = [f"{l} ({c})" for l, c in comp]
        axp.pie(vals, labels=labs, colors=[pie_color(i) for i in range(len(comp))],
                startangle=90, autopct="%1.0f%%",
                wedgeprops=dict(edgecolor="white", linewidth=0.6),
                textprops=dict(fontsize=13))
    axp.set_title(f"{target} sample composition", fontsize=16)

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, outname + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(HERE, outname + ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", outname, " n SBS:", len(sub))

build_ess_focus("eSS38", "eSS38_focus")

# --------------------------------- grid: every other decomposing eSS (one panel each)
def build_ess_grid(exclude, outname):
    # eSS that confidently decomposed COSMIC (same set as the bubble plots)
    targets = sorted({r[1] for r in ROWS} - set(exclude), key=ess_num)
    norm = Normalize(vmin=0.90, vmax=1.0)
    cmap = COS_CMAP

    ncol = 4
    nrow = int(np.ceil(len(targets) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 5.2, nrow * 2.5))
    axes = np.atleast_2d(axes).ravel()

    for ax, t in zip(axes, targets):
        sub = sorted([r for r in ROWS if r[1] == t], key=lambda r: sbs_key(r[0]))
        sbs_labels = [r[0] for r in sub]
        pcts = [r[2] for r in sub]
        coss = [r[3] for r in sub]
        y = np.arange(len(sub))
        ax.barh(y, pcts, color=[cmap(norm(c)) for c in coss],
                edgecolor="0.3", linewidth=0.5)
        ax.set_yticks(y); ax.set_yticklabels(sbs_labels, fontsize=11)
        ax.invert_yaxis()
        ax.tick_params(axis="x", labelsize=9)
        ax.set_title(f"{t}  ({primary_name(t)})", fontsize=11, fontweight="bold")
        for yi, p, c in zip(y, pcts, coss):
            ax.text(p + 1, yi, f"{p:.0f}%", va="center", fontsize=9)
        ax.set_xlim(0, max(pcts) * 1.32 if pcts else 1)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    for ax in axes[len(targets):]:
        ax.axis("off")

    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=axes.tolist(), fraction=0.015, pad=0.01)
    cb.set_label("Cosine similarity of decomposition", fontsize=14)
    cb.ax.tick_params(labelsize=12)

    fig.suptitle("eSS percent contribution to COSMIC v3.6 signatures",
                 fontsize=18, y=0.997)
    fig.savefig(os.path.join(HERE, outname + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(HERE, outname + ".png"), dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("wrote", outname, " panels:", len(targets))

build_ess_grid(exclude={"eSS38"}, outname="eSS_other_decomposers")
print("done")
