# `results_wes_to_wgs_normalized` vs `results_default`

**Same for both:** dataset (main 4,282-sample SBS atlas), clustering threshold
(0.9), averaging method (pooled), COSMIC comparison threshold (0.85). The only
difference is that 152 WES samples (9 human, 143 mouse) had their profiles
corrected onto their species' own WGS trinucleotide-opportunity basis before
clustering (see [pipeline/normalize_wes_to_wgs.py](../pipeline/normalize_wes_to_wgs.py));
the other ~750 WGS samples across all 5 species were untouched.

## Headline numbers

| | Default | WES→WGS corrected | Δ |
|---|---|---|---|
| Main clusters | 49 | 50 | +1 |
| Main-cluster samples | 508 | 512 | +4 |
| Small-cluster samples | 32 | 30 | −2 |
| Singletons | 130 | 128 | −2 |
| COSMIC-matched clusters (≥0.85) | 27 | 27 | 0 |
| De novo clusters | 22 | 23 | +1 |

## What actually moved (not just the counts)

- 6 samples changed main/non-main status — **all 6 are WES samples**:
  1 human `HK2_Aristolochic_acid_I` dropped out of main-cluster status; 5
  mouse `Liver_Diethylnitrosamine` samples newly qualified as main. Confirms
  the effect traces directly to the correction, not spillover onto untouched
  WGS samples.
- Among the 507 samples that stayed "main" in both runs, regrouping happened
  for ~4 clusters, most notably **`5-aza-4-thio-2-deoxycytidine`** (human
  CEM/U937 vs mouse thymus) and the mouse **`Diethylnitrosamine`** liver
  group — both compounds with matched WES+WGS representation.
- Per-sample: of the 152 corrected samples, cosine similarity to their own
  pre-correction profile ranges 0.895–0.991 (mean 0.969) — see
  [results_translation_verification/SBS/cosine_similarity_before_vs_after.tsv](../results_translation_verification/SBS/cosine_similarity_before_vs_after.tsv).
  Worst-affected are the same `5-aza-4-thio-2-deoxycytidine` and
  `Diethylnitrosamine` samples, plus `Mouse_MEF_Benzo[a]pyrene`.

## Bottom line

The overall match rate to COSMIC is essentially unchanged (27/49 → 27/50),
but the correction isn't a no-op — it measurably reshuffles cluster
membership for a handful of specific compounds where both sequencing
technologies are represented in the data.
