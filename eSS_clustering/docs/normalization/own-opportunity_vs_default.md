# `own-opportunity` normalization vs the default run

**Runs compared:** `results/min307/SBS/` (default) and
`results/min307_own-opportunity/SBS/`, made with

```bash
bash run_pipeline.sh SBS 0.9 0.85
bash run_pipeline.sh SBS 0.9 0.85 own-opportunity
```

**Same for both:** the 653 samples with ≥307 SBSs (10 mouse MEF samples
excluded), cosine similarity 0.90 with average linkage, no custom thresholds,
equal-replicate consensus, COSMIC v3.6 match at ≥0.85 (artifact signatures
included). The only difference: every human, mouse and rat sample (560
samples, WGS and WES) has its profile divided by its own genome and
technology trinucleotide opportunity (GRCh38, mm10, rn7; genome or exome) and
renormalized. Counts are unchanged. Chicken and C. elegans have no opportunity
tables and keep their plain profiles. See
[NORMALIZATION_APPROACH.md](../../NORMALIZATION_APPROACH.md).

## Headline numbers

| | Default | own-opportunity | Δ |
|---|---|---|---|
| Samples clustered | 653 | 653 | 0 |
| Main clusters (eSS) | 48 | 52 | +4 |
| Main-cluster samples | 498 | 493 | −5 |
| Small clusters | 16 | 15 | −1 |
| Singletons | 123 | 130 | +7 |
| COSMIC matched (≥0.85) | 25 | **14** | **−11** |
| COSMIC unmatched | 23 | 38 | +15 |

## What changed

- **Profiles change a lot.** Corrected profiles have cosine 0.57–0.996 (mean
  0.87) to their uncorrected profiles, much more than wes-to-wgs (0.90–0.99,
  WES only). The most affected are mouse MEF UVA and mouse mammary
  γ-ray/Fe-ion samples.
- **Only 22 of the 48 eSS keep the same members.** 2 have no counterpart
  afterwards, 11 new ones appear, and the rest gain, lose or merge samples. The largest changes are in
  mouse skin DMBA, mouse intestinal organoids (normal and high-fat diet),
  mouse liver diethylnitrosamine, mouse breast γ-ray/Fe-ion and human BEAS-2B
  benzo[a]pyrene.
- 48 samples change group: 22 (13 mouse, 9 human) leave the main clusters and
  17 (12 mouse, 5 human) join them. No rat, chicken or C. elegans sample
  changes group.
- **COSMIC matches fall from 25 to 14.** 12 eSS matched by default lose their
  match (2 go the other way). Even among the 22 eSS with identical members,
  4 flip from matched to unmatched, and their best COSMIC similarity drops by
  a median of 0.03.

## Why COSMIC matches drop

An earlier version of this summary suggested that the drop came from
clusters mixing samples on different opportunity bases (for example human
WGS with mouse WES). That is **not** the main cause: of the 12 eSS that lose
their match, 9 contain a single species and a single sequencing technology.
Among eSS of the corrected species, the single-species/technology ones lose
more similarity (median −0.033) than the mixed ones (median +0.010), and eSS
made only of chicken/C. elegans samples (uncorrected) are unchanged.

The more likely reason is that COSMIC signatures are expressed as
mutations observed on the human genome; they are not opportunity-normalized.
Dividing a profile by its trinucleotide opportunity puts it on a different
basis from COSMIC, so the comparison is no longer like for like, whatever the
cluster's composition.

**Suggested follow-up:** compare on a common basis, e.g. multiply the
own-opportunity consensus profiles by the human GRCh38 WGS opportunity
(COSMIC's basis) before the COSMIC comparison, and check whether the match
rate recovers.

## Bottom line

Per-sample opportunity normalization changes the clustering substantially
(22 of 48 eSS unchanged), mainly for mouse and human WGS samples. Its COSMIC
match rate should not be compared directly with the default run's until the
profiles and COSMIC are put on the same opportunity basis.

*Earlier version:* this comparison was first made on the previous inputs
(671 profiles, AAI/DBP split): 49 → 52 eSS, 26 → 14 COSMIC matches. That
version is in git history.
