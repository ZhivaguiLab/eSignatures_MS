#!/usr/bin/env python3
"""
Pan-cancer eSS-TMB bubble plot for the SPE de-novo hybrid run2 output.

A sample carries an eSS at >= 5% of its total mutation burden. Bubble colour is the
median mut/Mb of carriers, bubble size is prevalence within the cancer type, and
only eSS x cancer-type cells with prevalence >= 5% are drawn. Two variants are
written to the smoking_act5pct folder:

  TMB_eSS_pancancer_bubble.pdf          all cancer types
  TMB_eSS_pancancer_bubble_compact.pdf  cancer types with no qualifying eSS dropped
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle

BASE = os.environ.get("ESS_BASE", "../data/pancan_eSS_assignment")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling modules
import spa_generate_all_plots as S
import run2_generate as R

PARENT = os.path.join(BASE, "SPE_DeNovo_eSS51",
                      "outputs_denovo_hybrid_with_penalties_run2")
OUTDIR = os.path.join(PARENT, S.ANALYSIS_DIRNAME, "smoking_act5pct")
# also refresh the byte-identical baseline copy
OUTDIR2 = os.path.join(PARENT, S.ANALYSIS_DIRNAME)


def _fmt(x):
    return S._fmt(x)


def _place_labels_smart(ax, labels, num_font, yr, panel_h_in):
    """Label EVERY drawn bubble. Within each cancer-type column, spread the labels
    vertically so none overlap; a label that has to move is placed just to the right
    of its bubble (still inside the column's rectangle) and joined to it by a thin
    leader line. Labels that do not need to move stay centred on their bubble."""
    from collections import defaultdict
    data_per_pt = yr / (panel_h_in * 72.0)        # data-y units per typographic point
    min_gap = 1.3 * num_font * data_per_pt        # vertical spacing so labels can't touch
    cols = defaultdict(list)
    for x, y, l, prev in labels:
        cols[round(x - 0.5)].append((y, l))
    for i, items in cols.items():
        items.sort(key=lambda t: t[0])
        ys = [t[0] for t in items]
        adj = list(ys)
        for k in range(1, len(adj)):               # greedy upward spread
            if adj[k] < adj[k - 1] + min_gap:
                adj[k] = adj[k - 1] + min_gap
        shift = (sum(ys) / len(ys)) - (sum(adj) / len(adj))   # re-centre the cluster
        adj = [a + shift for a in adj]
        xc = i + 0.5
        for (y0, l), ya in zip(items, adj):
            if abs(ya - y0) < 1e-6:                 # no collision: keep centred
                ax.text(xc, y0, l, ha="center", va="center", fontsize=num_font,
                        fontweight="bold", zorder=4)
            else:                                   # moved: side-place + leader line
                lx = xc + 0.30
                ax.plot([xc, lx - 0.04], [y0, ya], color="#8a8a8a", lw=0.4,
                        zorder=3.4, clip_on=False)
                ax.text(lx, ya, l, ha="left", va="center", fontsize=num_font,
                        fontweight="bold", zorder=4)


def bubble(df, present, ess_cols, outfile, drop_empty=False,
           col_in=0.62, top_in=2.0, min_prev=0.05, min_carriers=1,
           label_prev=0.05, smart_labels=False,
           size_scale=11.0, ylabel="Median eSS mutations per\nmegabase (carriers)",
           ct_font=13, ax_font=15, tick_font=14, num_font=7, leg_font=11):
    """Cosmetically-improved copy of S.bubble_tmb. A cell is drawn only if its
    prevalence >= min_prev AND it has >= min_carriers carriers. drop_empty
    removes cancer types that end up with no qualifying eSS bubble."""
    ess_per_mb = df[ess_cols] / S.MB
    tot_mb = df[ess_cols].sum(axis=1) / S.MB
    order = (tot_mb[tot_mb > 0].groupby(df["CancerType"])
             .median().sort_values().index.tolist())

    # pre-compute which cancer types have >=1 qualifying bubble
    def _has_bubble(ct):
        m = df["CancerType"].values == ct
        for c in ess_cols:
            if present[c].values[m].mean() < min_prev:
                continue
            carriers = m & present[c].values
            if (carriers.sum() >= min_carriers
                    and np.median(ess_per_mb[c].values[carriers]) > 0):
                return True
        return False

    if drop_empty:
        order = [ct for ct in order if _has_bubble(ct)]
    N = len(order)

    xs, ys, sizes, cols, labels = [], [], [], [], []
    used_etio = set()
    for i, ct in enumerate(order):
        m = df["CancerType"].values == ct
        for c in ess_cols:
            prev = present[c].values[m].mean()
            if prev < min_prev:
                continue
            carriers = m & present[c].values
            if carriers.sum() < min_carriers:   # too few >=5%-activity carriers
                continue
            med = np.median(ess_per_mb[c].values[carriers])
            if med <= 0:
                continue
            etio = S.ess_etiology(c)
            used_etio.add(etio)
            xs.append(i + 0.5); ys.append(np.log10(med))
            sizes.append(prev * 100 * size_scale)
            cols.append(S.ETIO_COLOR[etio])
            labels.append((i + 0.5, np.log10(med), c.replace("eSS", ""), prev))
    if not ys:
        return set()

    ymin, ymax = np.floor(min(ys)), np.ceil(max(ys))
    ticks = np.arange(ymin, ymax + 1)
    LEFT_IN, RIGHT_IN, PANEL_H, BOT_IN = 2.0, 2.1, 3.4, 0.5
    panel_w_in = N * col_in
    fig_w = LEFT_IN + panel_w_in + RIGHT_IN
    fig_h = top_in + PANEL_H + BOT_IN
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([LEFT_IN / fig_w, BOT_IN / fig_h,
                       panel_w_in / fig_w, PANEL_H / fig_h])
    yr = ymax - ymin
    for k in range(N):
        face = "#E8E8E8" if k % 2 == 0 else "#FFFFFF"
        ax.add_patch(Rectangle((k, ymin), 1, yr, facecolor=face,
                               edgecolor="none", zorder=0))
    for t in ticks:
        ax.axhline(t, color="black", lw=1.0 if t in (ymin, ymax) else 0.6,
                   ls="-" if t in (ymin, ymax) else (0, (2, 3)), zorder=1)
    ax.scatter(xs, ys, s=sizes, c=cols, alpha=0.7, edgecolors="black",
               linewidths=0.5, zorder=3)
    if smart_labels:
        # label every drawn bubble, with per-column de-overlap + leader lines
        _place_labels_smart(ax, labels, num_font, yr, PANEL_H)
    else:
        for x, y, l, prev in labels:
            if prev >= label_prev:
                ax.text(x, y, l, ha="center", va="center", fontsize=num_font,
                        fontweight="bold", zorder=4)
    for i, ct in enumerate(order):
        ax.text(i + 0.5, ymax + 0.03 * yr, ct, ha="left", va="bottom",
                rotation=40, fontsize=ct_font, rotation_mode="anchor")
    ax.set_xlim(0, N); ax.set_ylim(ymin, ymax)
    ax.set_xticks([]); ax.set_yticks(ticks)
    ax.set_yticklabels([_fmt(10 ** t) for t in ticks], fontsize=tick_font)
    ax.tick_params(axis="y", length=0, pad=6)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.spines["left"].set_visible(True); ax.spines["right"].set_visible(True)
    ax.set_ylabel(ylabel, fontsize=ax_font, fontweight="bold", labelpad=8)

    # ---- size key: bubbles + % centred under a centred "prevalence" title ----
    bub_x, lbl_x = N + 1.05, N + 1.55
    title_x = (bub_x + lbl_x) / 2.0                      # centre of the key
    fracs = [0.10, 0.25, 0.50, 1.00]
    # even 0.42 spacing for 10/25/50%, then drop ONLY the 100% bubble an extra 0.1
    # so it barely clears the 50% without touching the rest of the key.
    ylegs = np.array([ymax - 0.35 - 0.42 * i for i in range(4)])
    ylegs[-1] -= 0.10
    for frac, yleg in zip(fracs, ylegs):
        ax.scatter([bub_x], [yleg], s=frac * 100 * size_scale, c="#999999",
                   alpha=0.7, edgecolors="black", linewidths=0.5,
                   clip_on=False, zorder=3)
        ax.text(lbl_x, yleg, f"{int(frac*100)}%", ha="left", va="center",
                fontsize=leg_font, clip_on=False)
    ax.text(title_x, ymax - 0.1, "prevalence", ha="center", va="bottom",
            fontsize=leg_font + 1, fontweight="bold", clip_on=False)

    fig.savefig(outfile, transparent=True)
    plt.close(fig)
    print("wrote", outfile, f"({N} cancer types)")
    return used_etio


def main():
    act = S.load_pooled(PARENT, sub=R.SUB)
    ess_cols = [c for c in act.columns if c.startswith("eSS")]
    ct = S.META[["Samples", "ctype"]].rename(columns={"ctype": "CancerType"})
    df = (act.merge(ct, on="Samples", how="inner")
          .dropna(subset=["CancerType"]).reset_index(drop=True))

    all_sigs = [c for c in df.columns if c not in ("Samples", "CancerType")]
    tot = df[all_sigs].sum(axis=1).replace(0, np.nan)
    present = df[ess_cols].div(tot, axis=0) >= S.ACT_MIN

    for outdir in (OUTDIR, OUTDIR2):
        # ONE legend is shipped for the whole folder, so it must cover every etiology drawn
        # in ANY variant -- not just the default 5%-prevalence bubble. Aflatoxin (eSS10 in
        # liver, eSS11 in lymphoid) sits at ~2% prevalence and is drawn only in the
        # no-floor variants, so a legend built from the default bubble alone omits it.
        used = set()
        used |= bubble(df, present, ess_cols,
                       os.path.join(outdir, "TMB_eSS_pancancer_bubble.pdf"),
                       drop_empty=False)
        used |= bubble(df, present, ess_cols,
                       os.path.join(outdir, "TMB_eSS_pancancer_bubble_compact.pdf"),
                       drop_empty=True)
        # no prevalence floor: the 5% stays on the signature activity (carrier
        # definition), but EVERY eSS x cancer-type cell with >=1 carrier is drawn
        # (min_prev=0), so this bubble matches the dot-matrix membership.
        used |= bubble(df, present, ess_cols,
                       os.path.join(outdir, "TMB_eSS_pancancer_bubble_allcarriers.pdf"),
                       drop_empty=False, min_prev=0.0)
        # de-noised: no prevalence floor, but require the eSS to be a >=5%-activity
        # carrier in >= 3 samples of the cancer type (drops sparse 1-2 sample dots).
        used |= bubble(df, present, ess_cols,
                       os.path.join(outdir, "TMB_eSS_pancancer_bubble_min3carriers.pdf"),
                       drop_empty=True, min_prev=0.0, min_carriers=3)
        # de-noised by PREVALENCE instead of absolute count: draw an eSS x cancer-type
        # cell only if >1% of that cancer type's samples are >=5%-activity carriers
        # (min_prev = 0.01). Scales the floor to cohort size rather than a fixed >=3.
        used |= bubble(df, present, ess_cols,
                       os.path.join(outdir, "TMB_eSS_pancancer_bubble_min1pct.pdf"),
                       drop_empty=True, min_prev=0.01, min_carriers=1)
        # same prevalence rule at a stricter 2% floor.
        used |= bubble(df, present, ess_cols,
                       os.path.join(outdir, "TMB_eSS_pancancer_bubble_min2pct.pdf"),
                       drop_empty=True, min_prev=0.02, min_carriers=1, smart_labels=True,
                       num_font=9, ylabel="Median eSS mutations per Mb")
        S.etiology_legend(used, os.path.join(outdir, "TMB_eSS_etiology_legend.pdf"))


if __name__ == "__main__":
    main()
    print("done.")
