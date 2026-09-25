#!/usr/bin/env python3
"""Bubble heatmap: mean mutations per Mb per compound (Exposure) across all models.
Compounds sorted by chemical group (shown as an x-axis color strip), and within
each group by mutagenic potential (max mean mut/Mb across models)."""
import csv, collections, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec

import os
DATA = os.environ.get("ESS_BURDEN_DATA", "../data/all_species_mut_burden.txt")
OUT  = os.environ.get("ESS_HEATMAP_OUT", "heatmap_mut_per_Mb_all_compounds.pdf")

# ---------- load + reproduce R transforms ----------
rows = list(csv.DictReader(open(DATA), delimiter="\t"))
rows = [r for r in rows if not (r["compounds"] == "xenon" and r["model"] == "MEF")]
for r in rows:
    if r["species"] == "Mouse" and r["type"] == "in-vivo":
        r["model"] = r["model"] + "_cancer"
    if "brain" in r["model"].lower():
        r["model"] = "Brain_cancer"

def fnum(x):
    try: return float(x)
    except: return np.nan

# ---------- summarise: mean mut/Mb + n samples per (model, Exposure) ----------
agg = collections.defaultdict(list)
for r in rows:
    v = fnum(r["sum_Mb"])
    if not np.isnan(v):
        agg[(r["model"], r["Exposure"])].append(v)
mean_mb = {k: np.mean(v) for k, v in agg.items()}
n_samp  = {k: len(v) for k, v in agg.items()}

models = sorted({m for m, _ in mean_mb}, key=str.lower)[::-1]   # alphabetical (case-insensitive), A at top
y_of = {m: i for i, m in enumerate(models)}

# ---------- compound ordering: chemical group, then mutagenic potential ----------
chem_of = {}
for r in rows:
    chem_of.setdefault(r["Exposure"], r["Chemical.group"])
potency = collections.defaultdict(float)               # max mean mut/Mb over models
for (m, e), v in mean_mb.items():
    potency[e] = max(potency[e], v)
exposures = sorted({e for _, e in mean_mb},
                   key=lambda e: (chem_of.get(e, "zzz"), -potency[e]))
x_of = {e: i for i, e in enumerate(exposures)}

# group -> ordered list & contiguous spans (sorted gives groups alphabetically)
group_order = []
for e in exposures:
    g = chem_of.get(e, "Unknown")
    if g not in group_order:
        group_order.append(g)
# distinct colors for chemical groups
base = (plt.get_cmap("tab20").colors + plt.get_cmap("tab20b").colors
        + plt.get_cmap("tab20c").colors)
gcolor = {g: base[i % len(base)] for i, g in enumerate(group_order)}

# ---------- scatter arrays ----------
xs, ys, sz, cl = [], [], [], []
for (m, e), mb in mean_mb.items():
    xs.append(x_of[e]); ys.append(y_of[m]); sz.append(n_samp[(m, e)]); cl.append(mb)
xs, ys, sz, cl = map(np.array, (xs, ys, sz, cl))

cval = np.log10(np.clip(cl, 1e-3, None))
norm = Normalize(vmin=np.floor(np.nanpercentile(cval, 2)),
                 vmax=np.ceil(np.nanpercentile(cval, 98)))
cmap = plt.get_cmap("magma_r")
smin, smax = sz.min(), sz.max()
def area(n):
    return 25 + 320 * (np.sqrt(n) - math.sqrt(smin)) / (math.sqrt(smax) - math.sqrt(smin) + 1e-9)

# ---------- figure ----------
W = max(16, len(exposures) * 0.17)
H = len(models) * 0.26 + 4.0
fig = plt.figure(figsize=(W, H))
gs = GridSpec(2, 1, height_ratios=[0.5, len(models)], hspace=0.02,
              left=0.13, right=0.86, top=0.95, bottom=0.30)
strip = fig.add_subplot(gs[0])
ax    = fig.add_subplot(gs[1], sharex=strip)

# chemical-group color strip (top)
i = 0
while i < len(exposures):
    g = chem_of.get(exposures[i], "Unknown"); j = i
    while j < len(exposures) and chem_of.get(exposures[j], "Unknown") == g:
        j += 1
    strip.add_patch(Rectangle((i - 0.5, 0), j - i, 1, facecolor=gcolor[g],
                              edgecolor="white", linewidth=0.4))
    i = j
strip.set_ylim(0, 1); strip.axis("off")

# main scatter
ax.set_axisbelow(True)
ax.grid(True, color="#ececec", lw=0.5, zorder=0)
ax.scatter(xs, ys, s=area(sz), c=cval, cmap=cmap, norm=norm,
           edgecolors="#333333", linewidths=0.4, alpha=0.95, zorder=3)

ax.set_xticks(range(len(exposures)))
ax.set_xticklabels(exposures, rotation=90, ha="center", va="top", fontsize=6)
ax.set_yticks(range(len(models)))
ax.set_yticklabels([m.replace("_", " ") for m in models], fontsize=8.5, color="black")
ax.set_xlim(-0.7, len(exposures) - 0.3)
ax.set_ylim(-0.7, len(models) - 0.3)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_color("#cccccc")
strip.set_title("Mean mutations per Mb", fontsize=15, loc="left", pad=8)

# ---------- legends (kept separated, no overlap) ----------
# colorbar (right, top)
cax = fig.add_axes([0.875, 0.62, 0.012, 0.30])
cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
ticks = range(int(np.floor(norm.vmin)), int(np.ceil(norm.vmax)) + 1)
cb.set_ticks(list(ticks)); cb.set_ticklabels([f"{10**t:g}" for t in ticks])
cb.set_label("Mean mutations / Mb", fontsize=11)
cb.ax.tick_params(labelsize=10); cb.outline.set_visible(False)

# size legend (right, below colorbar)
leg_n = [v for v in [5, 20, 40, 60, 80] if v <= smax]
handles = [Line2D([], [], marker="o", ls="", markerfacecolor="#888",
                  markeredgecolor="#333", markeredgewidth=0.4,
                  markersize=np.sqrt(area(v)), label=str(v)) for v in leg_n]
fig.legend(handles=handles, title="N samples", loc="center left",
           bbox_to_anchor=(0.87, 0.40), labelspacing=1.6, frameon=False,
           fontsize=10, title_fontsize=11, borderpad=1)

# chemical-group legend (bottom, multi-column)
g_handles = [Line2D([], [], marker="s", ls="", markerfacecolor=gcolor[g],
                    markeredgecolor="none", markersize=10, label=g) for g in group_order]
fig.legend(handles=g_handles, title="Chemical group", loc="upper center",
           bbox_to_anchor=(0.5, 0.16), ncol=5, frameon=False,
           fontsize=9, title_fontsize=11, columnspacing=1.2,
           handletextpad=0.4, labelspacing=0.4)

# ---------- CSV of the plotted values ----------
CSV = OUT.replace(".pdf", ".csv")
with open(CSV, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["model", "Exposure", "Chemical_group",
                "mean_mut_per_Mb", "n_samples",
                "group_potency_max_mut_per_Mb", "x_order", "y_order"])
    for (m, e), mb in sorted(mean_mb.items(), key=lambda kv: (x_of[kv[0][1]], y_of[kv[0][0]])):
        w.writerow([m, e, chem_of.get(e, "Unknown"),
                    f"{mb:.6g}", n_samp[(m, e)],
                    f"{potency[e]:.6g}", x_of[e], y_of[m]])
print("saved:", CSV)

fig.savefig(OUT, bbox_inches="tight")
fig.savefig(OUT.replace(".pdf", ".png"), dpi=170, bbox_inches="tight")
print("models:", len(models), "compounds:", len(exposures),
      "groups:", len(group_order), "points:", len(xs))
print("saved:", OUT)
