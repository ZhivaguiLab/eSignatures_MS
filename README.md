# Code availability

Analysis code for the main-figure panels and the pan-cancer association
statistics in the experimental SBS signature (eSS) compendium. Signature
extraction and assignment used the SigProfiler suite; the scripts here cover the
downstream analysis, statistics, and figure generation.

## Layout

```
figures/
  fig02_composition_coverage/    Sankey, model/exposure coverage, pie legends
  fig03_mutational_landscapes/   per-model mutation-burden heatmap
  fig05_cosmic_decomposition/    eSS vs COSMIC v3.6 decomposition
  fig06_pancancer_attribution/   pan-cancer eSS attribution, TMB, smoking, geography
  fig07_organoid_validation/     organoid NanoSeq eSS decomposition
```

## Figures and scripts

| Figure | Panel | Script |
|---|---|---|
| 2 | Sankey | `fig02_composition_coverage/river_plot.py` |
| 2 | pie legends | `fig02_composition_coverage/5b_make_pie_legends.py` |
| 3 | burden heatmap | `fig03_mutational_landscapes/20_heatmap_mean_mutations_model.py` |
| 5 | decomposition bubble | `fig05_cosmic_decomposition/5b_make_plots.py` |
| 6 | attribution, TMB | `fig06_pancancer_attribution/run2_bubble_improved.py` |
| 6 | smoking association | `fig06_pancancer_attribution/run2_smoking_cohort.py`, `run2_forest_per_cancer.py` |
| 6 | geographical association | `fig06_pancancer_attribution/run2_geographical.py` |
| 6 | etiology by cancer type | `fig06_pancancer_attribution/make_figures.py` |
| 7 | organoid decomposition | `fig07_organoid_validation/organoid_bubble_make_plots.py` |

`run2_generate.py` and `spa_generate_all_plots.py` in `fig06_pancancer_attribution/`
hold the shared pooling and plotting routines the other Figure 6 scripts import.

## Software

Python 3.10: SigProfilerMatrixGenerator 1.3.6, SigProfilerAssignment 1.1.4,
SigProfilerPlotting 1.4.3, pandas 2.3.3, numpy 2.2.6, scipy 1.13.1, statsmodels,
matplotlib. Reference genomes GRCh38 (human) and mm10 (mouse); SBS96 context.

R 4.4: readr, dplyr, tidyr, and base stats (glm, wilcox.test, p.adjust).

Pan-cancer decomposition (SigProfilerAssignment, COSMIC v3.5 database) used
`nnls_add_penalty = 0.10`, `nnls_remove_penalty = 0.01`,
`initial_remove_penalty = 0.05`, `de_novo_fit_penalty = 0.02`.
