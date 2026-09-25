#!/usr/bin/env python3
"""
Regenerate TMB + smoking + Fisher plots for the SPE de-novo HYBRID run2 output,
reusing the canonical pipeline in spa_generate_all_plots.py without modifying its
statistics.  The run2 output is the per-cancer-type DECOMPOSE solution produced
with the NNLS penalties applied (samples x [endogenous SBS + eSS + SBS96N]); the
per-cancer folders are pooled the same way the canonical script pools the
samples-as-input assignment.

Outputs for the run2 folder:
  analysis/                    baseline: TMB (total-TMB 5% carrier) + smoking
                               (detection > 0) + Fisher
  analysis/smoking_act5pct/    smoking + Fisher, eSS/sig kept only if >= 5% of
                               the sample's TOTAL mutation burden (sub-threshold
                               zeroed), i.e. a 5% activity cutoff
  analysis/smoking_act7pct/    same, 7% cutoff

The activity cutoff is a "trusted activity" filter: activities below the cutoff
(relative to total TMB across ALL signature columns) are set to 0 BEFORE the
smoking/Fisher stats, so no change to the canonical functions is needed.
"""
import os, sys
import numpy as np, pandas as pd

# Data root. Point ESS_BASE at the folder holding the pan-cancer assignment data
# (see README); the default is a generic relative path, not a machine-specific one.
BASE = os.environ.get("ESS_BASE", "../data/pancan_eSS_assignment")
# Sibling modules (spa_generate_all_plots) live next to this script.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spa_generate_all_plots as S     # module-level: loads metadata + labels

SPE = os.path.join(BASE, "SPE_DeNovo_eSS51")
SUB = "Decompose_Solution/Activities/Decompose_Solution_Activities.txt"
CUTOFFS = [0.05, 0.07]                 # trusted-activity cutoffs vs TOTAL TMB

# Only the penalized run2 decomposition is used in the manuscript.
RUNS = [
    ("SPE_denovo_hybrid_with_penalties_run2",
     os.path.join(SPE, "outputs_denovo_hybrid_with_penalties_run2")),
]


def filter_activities(act_df, act_min):
    """Zero every signature contribution below act_min of the sample's TOTAL
    mutation burden (sum over ALL signature columns)."""
    sigs = [c for c in act_df.columns if c != "Samples"]
    counts = act_df[sigs]
    tot = counts.sum(axis=1).replace(0, np.nan)
    keep = counts.div(tot, axis=0) >= act_min
    out = act_df.copy()
    out[sigs] = counts.where(keep, 0.0)
    return out


def smoking_only(act_df, outdir, tag):
    """Run just the smoking volcanoes + Fisher (no TMB) into outdir."""
    os.makedirs(outdir, exist_ok=True)
    sigs = [c for c in act_df.columns if c != "Samples"]
    sdf = act_df.merge(S.META, on="Samples", how="inner").dropna(
        subset=["smoker", "ctype"]).reset_index(drop=True)
    sdf["smoker"] = sdf["smoker"].astype(int)
    n_sm, n_ns = int((sdf.smoker == 1).sum()), int((sdf.smoker == 0).sum())
    S.smoking_block(sdf, sigs, outdir, tag)
    nt, npv, nbd = S.fisher_block(sdf, sigs, outdir, tag)
    print(f"    {tag}: smoking df={len(sdf)} (ever={n_sm} never={n_ns}); "
          f"fisher {nt} tissues, {npv} sig prevalence, {nbd} sig burden -> {outdir}")


if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else None
    for name, parent in RUNS:
        if only and only not in name:
            continue
        act = S.load_pooled(parent, sub=SUB)
        analysis = os.path.join(parent, S.ANALYSIS_DIRNAME)
        # baseline: TMB (fixed denom) + smoking (detection>0) + Fisher
        S.run(name, act, analysis)
        # activity-cutoff comparison folders (smoking + Fisher only)
        for thr in CUTOFFS:
            act_f = filter_activities(act, thr)
            outdir = os.path.join(analysis, f"smoking_act{int(thr*100)}pct")
            smoking_only(act_f, outdir, f"{name}_act{int(thr*100)}pct")
    print("\nrun2 done.")
