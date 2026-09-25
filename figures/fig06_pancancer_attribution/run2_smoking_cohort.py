#!/usr/bin/env python3
"""
Cohort-adjusted smoking plots for the SPE de-novo HYBRID run2 outputs, adding the
two figures missing from the baseline run2 analysis:

  A. volcano_cancertype_cohort_adjusted[_eSS].{pdf,csv}
       pooled within-tissue logistic, presence ~ smoker + age + sex
       + C(cancer type) + C(cohort)   -> adds COHORT to the existing
       cancer-type-adjusted volcano.

  B. volcano_per_cancer_allsig.{pdf,png} + per_cancer_alleSS_results.csv
       one point per (eSS x cancer type) association from the WITHIN-tissue,
       cohort-adjusted model  detected ~ smoker + age + sex + C(cohort);
       non-significant = grey, significant (q<0.05) COLOURED per eSS and
       SHAPED per cancer type (the "eSS x tissue significance" volcano).

Generated into the run2 smoking folders: analysis/, analysis/smoking_act5pct/,
analysis/smoking_act7pct/  (penalized run2 decomposition).  Reuses
spa_generate_all_plots helpers (volcano, labels, refdr) and run2_generate
(load + activity filter).
3 cohorts PCAWG + TCGA + Mutographs; MB = 2800.
"""
import os, sys, warnings
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import PerfectSeparationWarning, ConvergenceWarning
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
warnings.simplefilter("ignore")

BASE = os.environ.get("ESS_BASE", "../data/pancan_eSS_assignment")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling modules
import spa_generate_all_plots as S
import run2_generate as R

MIN_N, Q_THR = 20, 0.05

# --- metadata incl. COHORT (S.META intentionally drops cohort) ---------------
_m = pd.read_csv(S.META_FILE, sep="\t", dtype=str, low_memory=False)
_m = _m[_m["Sample_ID"] != "Sample_ID"]
_m = _m[_m["Cohort"].isin(S.COHORTS)]
_m["smoker"] = _m["Smoking_Status"].apply(S.classify)
_m["age"] = pd.to_numeric(_m["Age"], errors="coerce")
_m["female"] = _m["Sex"].map({"Female": 1, "Male": 0})
_m["cohort"] = _m["Cohort"]; _m["ctype"] = _m["Cancer_Type_HD_Grouping"]
METAC = _m[["Sample_ID", "smoker", "age", "female", "cohort", "ctype"]].rename(
    columns={"Sample_ID": "Samples"})


def build_sdf(act_df):
    sdf = act_df.merge(METAC, on="Samples", how="inner").dropna(
        subset=["smoker", "ctype"]).reset_index(drop=True)
    sdf["smoker"] = sdf["smoker"].astype(int)
    return sdf


def short_label(ess):
    # COSMIC identity at the 0.90-cosine rule, with curated overrides for the beyond-49 /
    # special signatures (eSS50=SBS7d, eSS51=SBS17b). eSS37's old SBS92-like/tobacco call
    # is retired, so it now resolves to "novel" like any other sub-threshold eSS.
    return S._ess_cosmic_short(ess)


# ============================================================ A. cohort-adjusted pooled
def cohort_adjusted_pooled(sdf, sigs, outdir, tag):
    det = (sdf[sigs] > 0).astype(int)
    base = sdf[["smoker", "age", "female", "ctype", "cohort"]].copy()
    rows = []
    for s in sigs:
        dd = base.copy(); dd["y"] = det[s].values
        dd = dd.dropna(subset=["age", "female"])
        keep = [t for t, g in dd.groupby("ctype") if g.y.nunique() == 2]
        dd = dd[dd.ctype.isin(keep)]
        if dd.y.nunique() < 2 or dd.smoker.nunique() < 2 or len(dd) < 50:
            rows.append({"signature": s, "log2_OR": np.nan, "p": np.nan}); continue
        f = "y ~ smoker + age + female + C(ctype)" + (" + C(cohort)" if dd.cohort.nunique() > 1 else "")
        try:
            m = smf.logit(f, data=dd).fit(disp=0, maxiter=300)
            rows.append({"signature": s, "log2_OR": m.params["smoker"] / np.log(2),
                         "p": m.pvalues["smoker"]})
        except Exception:
            rows.append({"signature": s, "log2_OR": np.nan, "p": np.nan})
    adj = pd.DataFrame(rows)
    ok = adj.p.notna(); adj["q"] = np.nan
    adj.loc[ok, "q"] = multipletests(adj.loc[ok, "p"], method="fdr_bh")[1]
    adj.sort_values("q").to_csv(
        os.path.join(outdir, "volcano_cancertype_cohort_adjusted.csv"), index=False)
    adj_e = S.refdr(adj)
    adj_e.sort_values("q").to_csv(
        os.path.join(outdir, "volcano_cancertype_cohort_adjusted_eSS.csv"), index=False)
    XLAB = r"$\log_2$ OR (ever vs never smoker)"
    SIDE = {"eSS33": "center", "eSS17": "topright"}                 # COSMIC / all-sig plots
    SIDE_E = {"eSS33": "center", "eSS17": "topright", "eSS47": "topright"}  # etiology plot
    S.volcano(adj, "log2_OR",
              "Smoking association (type + cohort adjusted)",
              os.path.join(outdir, "volcano_cancertype_cohort_adjusted.pdf"),
              XLAB, xclip=6.0, label_side=SIDE)
    # eSS-only, labelled by COSMIC match
    S.volcano(adj_e, "log2_OR", "Smoking association across eSS",
              os.path.join(outdir, "volcano_cancertype_cohort_adjusted_eSS.pdf"),
              XLAB, xclip=6.0, label_side=SIDE)
    # eSS-only, labelled by etiology (dominant training mutagen / chemical group)
    S.volcano(adj_e, "log2_OR", "Smoking association across eSS",
              os.path.join(outdir, "volcano_cancertype_cohort_adjusted_eSS_etiology.pdf"),
              XLAB, xclip=6.0, label_fn=S._volcano_label_etio, label_side=SIDE_E,
              bold_novel=False, label_dx={"eSS47": -4},   # nudge left off the right border
              transparent=True)
    return int((adj_e.q < Q_THR).sum())


# ============================================================ B. eSS x tissue volcano (shape/colour)
def _fit_logit(g, ess):
    d = g[[ess, "smoker", "age", "female", "cohort"]].dropna(subset=["age", "female"]).copy()
    d = d.rename(columns={ess: "y"})
    if d["y"].nunique() < 2 or d["smoker"].nunique() < 2 or len(d) < 30:
        return (np.nan,) * 4
    f = "y ~ smoker + age + female" + (" + C(cohort)" if d["cohort"].nunique() > 1 else "")
    try:
        m = smf.logit(f, data=d).fit(disp=0, maxiter=300)
        b = m.params["smoker"]; ci = m.conf_int().loc["smoker"]
        l2, lo, hi, p = b / np.log(2), ci[0] / np.log(2), ci[1] / np.log(2), m.pvalues["smoker"]
        if not np.isfinite([l2, lo, hi, p]).all():
            return (np.nan,) * 4
        return l2, lo, hi, p
    except Exception:
        return (np.nan,) * 4


def per_cancer_allsig_volcano(sdf, ess_cols, outdir, tag):
    det = (sdf[ess_cols] > 0).astype(int)
    dfp = det.copy()
    for c in ["smoker", "age", "female", "cohort", "ctype"]:
        dfp[c] = sdf[c].values
    types = sorted([t for t, g in dfp.groupby("ctype")
                    if (g.smoker == 1).sum() >= MIN_N and (g.smoker == 0).sum() >= MIN_N])
    rows = []
    for ess in ess_cols:
        for t in types:
            l2, lo, hi, p = _fit_logit(dfp[dfp.ctype == t], ess)
            rows.append({"eSS": ess, "ctype": t, "log2_OR": l2, "ci_lo": lo,
                         "ci_hi": hi, "logit_p": p})
    res = pd.DataFrame(rows)
    ciw = res["ci_hi"] - res["ci_lo"]
    sep = res["log2_OR"].notna() & ((ciw < 0.50) | (ciw > 15.0) | (res["log2_OR"].abs() > 10.0))
    res.loc[sep, ["log2_OR", "ci_lo", "ci_hi", "logit_p"]] = np.nan
    ok = res["logit_p"].notna()
    res.loc[ok, "qvalue_BH"] = multipletests(res.loc[ok, "logit_p"], method="fdr_bh")[1]
    res = res[res["log2_OR"].notna() & res["qvalue_BH"].notna()].copy()
    if res.empty:
        print(f"    [{tag}] per_cancer_allsig: no estimable cells"); return 0
    res["y"] = -np.log10(res["qvalue_BH"].clip(lower=1e-300))
    res["sig"] = res["qvalue_BH"] < Q_THR
    res.to_csv(os.path.join(outdir, "per_cancer_alleSS_results.csv"), index=False)

    sig = res[res["sig"]]
    sig_eSS = sorted(sig["eSS"].unique(), key=lambda e: int(e.replace("eSS", "")))
    sig_tis = sorted(sig["ctype"].unique())
    cmap = plt.get_cmap("tab20" if len(sig_eSS) > 10 else "tab10")
    ess_color = {e: cmap(i % cmap.N) for i, e in enumerate(sig_eSS)}
    markers = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">", "h"]
    tis_marker = {t: markers[i % len(markers)] for i, t in enumerate(sig_tis)}

    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    fig.subplots_adjust(right=0.70)
    ns = res[~res["sig"]]
    ax.scatter(ns["log2_OR"], ns["y"], s=26, c="#bdbdbd", alpha=0.5, edgecolor="none", zorder=2)
    for _, r in sig.iterrows():
        ax.scatter(r["log2_OR"], r["y"], s=78, c=[ess_color[r["eSS"]]],
                   marker=tis_marker[r["ctype"]], edgecolor="black", linewidth=0.7,
                   alpha=0.95, zorder=4)
    ax.axhline(-np.log10(Q_THR), ls="--", lw=0.8, color="0.6", zorder=1)
    ax.axvline(0, lw=0.6, color="0.3", zorder=1)
    xmin, xmax = res["log2_OR"].min(), res["log2_OR"].max()
    ymax = res["y"].max(); xpad = 0.12 * (xmax - xmin + 1e-9)
    ax.set_xlim(xmin - xpad, xmax + xpad); ax.set_ylim(-0.05 * ymax, ymax * 1.12)
    ax.tick_params(labelsize=8)
    ax.set_xlabel(r"$\log_2$ OR (ever vs never smoker)", fontsize=9)
    ax.set_ylabel(r"$-\log_{10}$(q-value)", fontsize=9)
    ax.set_title("Smoking association per cancer type", fontsize=9)
    ess_handles = [Line2D([0], [0], marker="o", linestyle="", markerfacecolor=ess_color[e],
                          markeredgecolor="black", markeredgewidth=0.4, markersize=7,
                          label=f"{e}  {short_label(e)}") for e in sig_eSS]
    tis_handles = [Line2D([0], [0], marker=tis_marker[t], linestyle="", markerfacecolor="0.5",
                          markeredgecolor="black", markeredgewidth=0.5, markersize=7,
                          label=t) for t in sig_tis]
    leg1 = ax.legend(handles=ess_handles, title="eSS (COSMIC)", loc="upper left",
                     bbox_to_anchor=(1.03, 1.0), fontsize=8, title_fontsize=8.5,
                     frameon=False, handletextpad=0.5, labelspacing=0.45)
    ax.add_artist(leg1)
    ax.legend(handles=tis_handles, title="Cancer type", loc="lower left",
              bbox_to_anchor=(1.03, 0.0), fontsize=8, title_fontsize=8.5,
              frameon=False, handletextpad=0.5, labelspacing=0.45)
    fig.savefig(os.path.join(outdir, "volcano_per_cancer_allsig.png"), dpi=160,
                bbox_inches="tight", pad_inches=0.45)
    fig.savefig(os.path.join(outdir, "volcano_per_cancer_allsig.pdf"),
                bbox_inches="tight", pad_inches=0.45)
    plt.close()
    return int(res.sig.sum())


def process(act_df, outdir, tag):
    os.makedirs(outdir, exist_ok=True)
    sdf = build_sdf(act_df)
    sigs = [c for c in act_df.columns if c != "Samples"]
    ess_cols = [c for c in sigs if c.startswith("eSS")]
    n_e = cohort_adjusted_pooled(sdf, sigs, outdir, tag)
    n_c = per_cancer_allsig_volcano(sdf, ess_cols, outdir, tag)
    print(f"    {tag}: cohort-adj eSS sig(q<.05)={n_e}; per-cancer allsig sig cells={n_c} -> {outdir}")


if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else None
    for name, parent in R.RUNS:
        if only and only not in name:
            continue
        act = S.load_pooled(parent, sub=R.SUB)
        analysis = os.path.join(parent, S.ANALYSIS_DIRNAME)
        process(act, analysis, name)                                   # detection>0
        for thr in R.CUTOFFS:
            act_f = R.filter_activities(act, thr)
            process(act_f, os.path.join(analysis, f"smoking_act{int(thr*100)}pct"),
                    f"{name}_act{int(thr*100)}pct")
    print("\ncohort-adjusted + per-cancer eSS-tissue volcanoes done.")
