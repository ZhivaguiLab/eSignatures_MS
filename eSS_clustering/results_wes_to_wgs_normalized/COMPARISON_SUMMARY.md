# `wes-to-wgs` normalization vs the default run

**Runs compared:** `results/min307_cos0.90/SBS/` (default) and
`results/min307_cos0.90_wes-to-wgs/SBS/`, made with

```bash
bash run_pipeline.sh SBS 0.9 0.85
bash run_pipeline.sh SBS 0.9 0.85 wes-to-wgs
```

**Same for both:** the 653 samples with ≥307 SBSs (10 mouse MEF samples
excluded), cosine similarity 0.90 with average linkage, no custom thresholds,
equal-replicate consensus, COSMIC v3.6 match at ≥0.85 (artifact signatures
included). The only difference: 152 WES samples (9 human, 143 mouse) have
their counts moved onto their own genome's WGS trinucleotide basis (GRCh38,
mm10) before clustering; each sample keeps its total. All WGS samples and all
rat, chicken and C. elegans samples are unchanged. See
[NORMALIZATION_APPROACH.md](../NORMALIZATION_APPROACH.md).

## Headline numbers

| | Default | wes-to-wgs | Δ |
|---|---|---|---|
| Samples clustered | 653 | 653 | 0 |
| Main clusters (eSS) | 48 | 50 | +2 |
| Main-cluster samples | 498 | 502 | +4 |
| Small clusters | 16 | 15 | −1 |
| Singletons | 123 | 121 | −2 |
| COSMIC matched (≥0.85) | 25 | 26 | +1 |
| COSMIC unmatched | 23 | 24 | +1 |

## What changed

- **45 of the 48 eSS keep exactly the same members.** Their profiles move
  only where they contain corrected WES samples (lowest cosine to the default
  profile 0.969).
- **The 5-aza-4-thio-2-deoxycytidine eSS splits by species.** The default
  48-sample eSS1 loses the 8 human CEM/U937 WES samples, which form their own
  eSS (best COSMIC SBS98, 0.64, unmatched); the mouse samples stay together
  (SBS98, 0.73 → 0.76).
- **Mouse liver diethylnitrosamine regroups.** 5 WES samples form a new eSS
  (best SBS96, 0.68, unmatched) and the existing diethylnitrosamine eSS swaps
  4 samples for 4 others.
- **The AAI/DBP eSS loses `Human_HK2_Aristolochic_acid_I`** (a WES sample),
  going from 11 to 10 samples; still matched to SBS22a (0.96).
- **One eSS gains a COSMIC match:** an SBS40a-like eSS with identical members
  goes from 0.815 to 0.859.
- 6 samples change group: 1 human sample leaves the main clusters and 5 mouse
  diethylnitrosamine samples join them.

**Per-sample effect.** Each corrected WES profile has cosine 0.895–0.991
(mean 0.969) to its uncorrected profile. The most affected are human CEM
5-aza-4-thio-2-deoxycytidine (0.895–0.919) and mouse MEF benzo[a]pyrene
(0.918–0.924). `pipeline/build_translation_verification_matrices.py` writes
the before/after profiles for every corrected sample.

## Bottom line

The correction reshapes only WES samples and changes only the eSS that
contain them. The atlas is otherwise identical (45 of 48 eSS), and the COSMIC
match rate is essentially unchanged (25/48 → 26/50). The effect is
concentrated in compounds profiled by both WES and WGS, where it separates
samples by sequencing technology and species.

*Earlier version:* this comparison was first made on the previous inputs
(671 profiles, AAI/DBP split, pooled averaging): 49 → 50 eSS, 27 → 27 COSMIC
matches, with the same compounds affected. That version is in git history.
