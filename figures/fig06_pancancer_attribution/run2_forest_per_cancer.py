#!/usr/bin/env python3
"""
Per-cancer, cohort-adjusted FOREST plot for the SPE DeNovo eSS51 run2 output,
mirroring Smoking association analysis/adjusted_cancer_type_cohort/
forest_per_cancer_allsig.pdf but built from the DeNovo eSS51 decompose data.

Reads per_cancer_alleSS_results.csv written by run2_smoking_cohort.py (all eSS x
cancer type, within-tissue cohort-adjusted logistic detected ~ smoker + age + sex
+ C(cohort), FDR across cells, separation cells excluded). One forest panel per eSS
with >=1 significant cell (q<0.05); PER_ROW panels per row. Significant cells
coloured (red = up in smokers, blue = down), n.s. grey; OR/q annotated.

Output: forest_per_cancer_allsig.{pdf,png} into the smoking_act5pct folder of the
active analysis dir (default the 4-cohort subset). Honours ESS_COHORTS /
ESS_MUTO_KIDNEY_ONLY / ESS_ANALYSIS_DIR like the other run2 scripts.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE = os.environ.get("ESS_BASE", "../data/pancan_eSS_assignment")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling modules
import spa_generate_all_plots as S

PARENT = os.path.join(BASE, "SPE_DeNovo_eSS51",
                      "outputs_denovo_hybrid_with_penalties_run2")
# match the source folder used to generate this plot (smoking_act5pct by default)
SUBDIR = os.environ.get("ESS_FOREST_SUBDIR", "smoking_act5pct")
OUTDIR = os.path.join(PARENT, S.ANALYSIS_DIRNAME, SUBDIR)

XMAX = 6.0
PER_ROW = 11
col_sm, col_ns, col_g = "#d1495b", "#3a7ca5", "#bdbdbd"


def short_label(ess):
    return S._ess_cosmic_short(ess)


def main():
    res = pd.read_csv(os.path.join(OUTDIR, "per_cancer_alleSS_results.csv"))
    POS = sorted(res.loc[res["sig"] == True, "eSS"].unique(),
                 key=lambda e: int(e.replace("eSS", "")))
    types = sorted(t for t in res["ctype"].unique() if t != "Glandular_Reproductive")
    ypos = {t: i for i, t in enumerate(types[::-1])}
    if not POS:
        print("no significant eSS cells -> no forest drawn"); return

    ncol = min(PER_ROW, len(POS))
    nrow = int(np.ceil(len(POS) / ncol))
    fig, axes = plt.subplots(nrow, ncol,
                             figsize=(2.0 * ncol, 0.34 * len(types) * nrow + 1.6),
                             squeeze=False, sharey=True)

    for k, ess in enumerate(POS):
        ax = axes[k // ncol][k % ncol]
        sub = res[res.eSS == ess]
        for _, r in sub.iterrows():
            if r.ctype not in ypos or not np.isfinite(r.log2_OR):
                continue
            y = ypos[r.ctype]
            color = col_g if not r.sig else (col_sm if r.log2_OR >= 0 else col_ns)
            lo = max(r.ci_lo, -XMAX); hi = min(r.ci_hi, XMAX)
            pt = np.clip(r.log2_OR, -XMAX, XMAX)
            ax.errorbar(pt, y, xerr=[[pt - lo], [hi - pt]], fmt="o", ms=4, color=color,
                        ecolor=color, elinewidth=1, capsize=2, zorder=3)
            if r.sig:
                OR = 2.0 ** r.log2_OR
                if r.qvalue_BH < 1e-300:
                    qt = "q<1e-300"
                elif r.qvalue_BH < 1e-3:
                    qt = f"q={r.qvalue_BH:.0e}"
                else:
                    qt = f"q={r.qvalue_BH:.3f}"
                tx = min(max(pt, -XMAX + 1.8), XMAX - 1.8)
                ax.text(tx, y + 0.34, f"OR {OR:.1f}, {qt}", ha="center", va="bottom",
                        fontsize=4.6, zorder=5)
        ax.axvline(0, color="0.4", lw=0.8, ls="--")
        ax.set_title(f"{ess}\n{short_label(ess)}", fontsize=8, fontweight="bold")
        ax.set_xlim(-XMAX - 0.4, XMAX + 0.4)
        ax.set_ylim(-0.6, len(types) - 0.1)
        ax.tick_params(labelsize=7)
        if k // ncol == nrow - 1:
            ax.set_xlabel("log2(OR ever/never)", fontsize=7.5)

    for k in range(len(POS), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    for r_ in range(nrow):
        axes[r_][0].set_yticks(range(len(types)))
        axes[r_][0].set_yticklabels(types[::-1], fontsize=7)

    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=col_sm, label="↑ in smokers (q<0.05)", ms=6),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=col_ns, label="↓ in smokers (q<0.05)", ms=6),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=col_g, label="n.s.", ms=6)]
    fig.legend(handles=leg, loc="lower center", ncol=3, fontsize=8, frameon=True,
               bbox_to_anchor=(0.5, -0.05))
    fig.suptitle("Per cancer type, cohort-adjusted — significant eSS "
                 "(presence ~ smoking + age + sex + C(cohort); log2 OR, 95% CI)\n"
                 "DeNovo eSS51 decompose · PCAWG+TCGA+Mutographs(RCC)+Sherlock",
                 fontsize=9, y=1.02)
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    fig.savefig(os.path.join(OUTDIR, "forest_per_cancer_allsig.png"), dpi=140, bbox_inches="tight")
    fig.savefig(os.path.join(OUTDIR, "forest_per_cancer_allsig.pdf"), bbox_inches="tight")
    plt.close()
    print(f"Forest eSS ({len(POS)}): {POS}; tissues={types}")
    print(f"Saved forest_per_cancer_allsig.pdf/.png -> {OUTDIR}")


if __name__ == "__main__":
    main()
