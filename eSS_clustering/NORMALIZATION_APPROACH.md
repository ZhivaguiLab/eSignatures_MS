# Trinucleotide opportunity normalization (manuscript revision item A4)

## The problem this addresses

The manuscript compares SBS-96 mutational profiles directly across species
(human, mouse, rat, chicken, C. elegans) and across sequencing technologies
(WGS vs. WES), interpreting similarity or divergence as biological (shared
or divergent repair/metabolism). But raw profiles from different genomes
and different capture territories aren't directly comparable: each of the
96 trinucleotide contexts occurs at a different baseline frequency
depending on genome composition and, for WES, which fraction of the genome
was captured. A raw similarity score conflates that composition confound
with actual biology. This is manuscript revision item **A4**.

Two independent methods are implemented as options of the pipeline's
preprocessing step, applied after the ≥307 mutation filter, and compared
against the default run. Neither is the default; choosing one is a decision
for the manuscript authors, informed by the results below.

```bash
bash run_pipeline.sh SBS 0.9 0.85                   # default → results/min307_cos0.90/SBS/
bash run_pipeline.sh SBS 0.9 0.85 wes-to-wgs        # Method 1 → results/min307_cos0.90_wes-to-wgs/SBS/
bash run_pipeline.sh SBS 0.9 0.85 own-opportunity   # Method 2 → results/min307_cos0.90_own-opportunity/SBS/
```

**Code:** [pipeline/utils/opportunity_normalization.py](pipeline/utils/opportunity_normalization.py)
(both methods), called from preprocessing in `pipeline/perform_clustering.py`.
**Opportunity tables:** `data/references/context_distributions/`
(SigProfilerMatrixGenerator 1.3.6; GRCh38, mm10, rn7, genome and exome).

---

## Method 1: liftover / translation

For every **WES** sample (human or mouse — the two species with confirmed
genome builds and both technologies represented in this dataset), rescale
its raw per-context counts to what they would look like if the same
underlying mutational process had instead been observed under **WGS**
opportunity for that same species/build:

```
ratio[c]     = WGS_relative_frequency[c] / WES_relative_frequency[c]
rescaled[c]  = raw_count[c] * ratio[c]
corrected[c] = rescaled[c] / sum(rescaled) * original_total_mutations
```

WGS samples, and every rat/chicken/celegans sample, are left completely
untouched. This treats WGS as the reference habit that WES samples get
translated onto, so all samples end up comparable to each other on a
single (WGS) basis. It does not address cross-species composition
differences.

This is the direct translation of the manuscript's existing
genome-to-genome/exome signature-crossover method
(`Marcos_New_reference_signature_files/2_crossover_genomes.R`), applied to
sample profiles instead of COSMIC reference signatures.

**Option:** `wes-to-wgs`
**Results:** [results_wes_to_wgs_normalized/COMPARISON_SUMMARY.md](results_wes_to_wgs_normalized/COMPARISON_SUMMARY.md)

**Headline result** (653 samples, ≥307 SBSs, no custom thresholds): main
clusters 48→50, COSMIC matches essentially unchanged (25/48→26/50); 45 of 48
eSS keep identical members. The changes are in compounds profiled by both WES
and WGS: human CEM/U937 `5-aza-4-thio-2-deoxycytidine` separates from the
mouse samples, and mouse liver `Diethylnitrosamine` regroups.

---

## Method 2: per-sample opportunity normalization

Every human, mouse, and rat sample — **WGS and WES alike** — is normalized
by its own species+genome-build+sequencing-technology trinucleotide
opportunity table, independently of every other sample:

```
rate[c]    = raw_count[c] / opportunity[c]   # that sample's own context-count table
profile[c] = rate[c] / sum(rate)              # renormalize to sum 1
```

`opportunity[c]` comes from a per-trinucleotide context-count table specific
to that exact combination: human WGS → GRCh38-WGS; human WES → GRCh38-exome;
mouse WGS → mm10-WGS; mouse WES → mm10-exome; rat WGS → rn7-WGS. No common
reference genome is chosen or needed — each sample only ever consults its
own opportunity table, never another sample's. Chicken and C. elegans have
no confirmed opportunity table yet, so their input files pass through
unchanged (still plain `raw / raw.sum()`) — this is a partial-coverage run,
not a full cross-species correction.

This is mathematically simpler than, and conceptually distinct from,
Method 1: it corrects every sample against its own opportunity, with no
target genome to choose, so it extends the same way to chicken and
C. elegans the moment their opportunity tables exist — no new design
decision needed, just the missing input data.

**Option:** `own-opportunity`
**Results:** [results_opportunity_normalized/COMPARISON_SUMMARY.md](results_opportunity_normalized/COMPARISON_SUMMARY.md)

**Headline result** (653 samples, ≥307 SBSs, no custom thresholds): main
clusters 48→52, only 22 of 48 eSS keep identical members, and COSMIC matches
drop sharply (25→14). The drop is **not** mainly from clusters mixing
opportunity bases (9 of the 12 eSS that lose their match are single species
and technology); the likelier reason is that COSMIC signatures are not
opportunity-normalized, so opportunity-normalized profiles are compared to
COSMIC on a different basis. See the comparison summary for the check and a
suggested follow-up.

---

## What was verified for both methods

- Empirically checked against the exact bug documented in
  `TRANSLATION_BUG_REPORT.md` (not included in this repository; a R→Python port of a 32-trinucleotide→96
  -context expansion once used `np.repeat` instead of `np.tile`, silently
  collapsing distinct trinucleotide ratios together). Both implementations
  do a label-driven dictionary lookup keyed by trinucleotide string, not
  positional array expansion — confirmed on the bug's own example block
  (`A[C>A]A/C/G/T`), which produced 4 distinct values, not 1, in both cases.
- The current implementation (`utils/opportunity_normalization.py`) gives the
  same profiles as the original standalone scripts it replaced (difference
  ≤1e-16 for every species), and the test suite checks both methods against
  an independent calculation from the opportunity tables.
- `wes-to-wgs` changes only WES counts and keeps every sample's total;
  `own-opportunity` leaves all counts unchanged and changes only the
  normalized profiles. So both keep exactly the samples that pass the 307
  filter.
- Every corrected profile sums to 1.
- Both ran through the full pipeline and were compared against the default
  run with the same settings otherwise (`results/min307_cos0.90/SBS/`); each
  run's settings and input hashes are in its `run_parameters.txt`.

## Status and open items

- **Neither method is the pipeline default.** Both are opt-in
  (`wes-to-wgs`, `own-opportunity`); the default run is unchanged.
- Genome builds (GRCh38, mm10, rn7) used by both methods were the best
  available guess, not confirmed against each source study's actual
  alignment — a real per-study build audit is still needed before either
  could be considered final.
- Chicken and C. elegans remain entirely uncorrected in both methods — no
  opportunity tables exist for those genomes yet; building them (from a
  reference genome FASTA) is the natural next step for either approach.
- The mixed-vs-homogeneous check for Method 2's COSMIC-match drop has been
  done (see its comparison summary): basis mixing is not the main cause.
  Comparing the profiles with COSMIC on a common opportunity basis is the
  suggested next step before drawing conclusions from Method 2's match rate.
