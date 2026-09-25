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

Two independent methods were implemented and compared against the same
`results_default` baseline. Both are complete, verified, and kept
available in the codebase side by side — neither has been adopted as *the*
answer; that's a decision for the manuscript authors, informed by the
results below.

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

**Code:** [pipeline/normalize_wes_to_wgs.py](pipeline/normalize_wes_to_wgs.py)
**Results:** [results_wes_to_wgs_normalized/COMPARISON_SUMMARY.md](results_wes_to_wgs_normalized/COMPARISON_SUMMARY.md)

**Headline result:** main clusters 49→50, COSMIC-match rate essentially
unchanged (27/49→27/50), but real regrouping for compounds with matched
WES/WGS representation (`5-aza-4-thio-2-deoxycytidine`,
`Diethylnitrosamine`).

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

**Code:** [pipeline/normalize_by_own_opportunity.py](pipeline/normalize_by_own_opportunity.py)
**Results:** [results_opportunity_normalized/COMPARISON_SUMMARY.md](results_opportunity_normalized/COMPARISON_SUMMARY.md)

**Headline result:** main clusters 49→52, COSMIC-match rate drops sharply
(26→14), 241 of 486 common main-cluster samples regrouped with different
partners. The COSMIC-match drop has a plausible mechanistic explanation
(clusters mixing samples from different opportunity bases average into
something that's on no single consistent basis, COSMIC's included) rather
than necessarily reflecting worse biology — see that comparison summary for
the full reasoning and an open follow-up to verify it directly.

---

## What was verified for both methods

- Empirically checked against the exact bug documented in
  `TRANSLATION_BUG_REPORT.md` (a R→Python port of a 32-trinucleotide→96
  -context expansion once used `np.repeat` instead of `np.tile`, silently
  collapsing distinct trinucleotide ratios together). Both implementations
  do a label-driven dictionary lookup keyed by trinucleotide string, not
  positional array expansion — confirmed on the bug's own example block
  (`A[C>A]A/C/G/T`), which produced 4 distinct values, not 1, in both cases.
- Raw counts confirmed byte-identical to the originals for every species,
  in both methods — only the normalized profile used for clustering and
  consensus averaging is ever affected.
- Every corrected profile confirmed to sum to exactly 1.0.
- Both ran through the full 6-step pipeline and were compared against the
  same `results_default` baseline (same clustering threshold, averaging
  method, and COSMIC comparison threshold as the rest of the atlas).

## Status and open items

- **Neither method is merged as the pipeline default.** Both scripts and
  their comparison summaries live in the codebase for review; `main`'s
  clustering behavior is unchanged unless one of these is explicitly run.
- Genome builds (GRCh38, mm10, rn7) used by both methods were the best
  available guess, not confirmed against each source study's actual
  alignment — a real per-study build audit is still needed before either
  could be considered final.
- Chicken and C. elegans remain entirely uncorrected in both methods — no
  opportunity tables exist for those genomes yet; building them (from a
  reference genome FASTA) is the natural next step for either approach.
- Method 2's COSMIC-match-rate drop needs the mixed-vs-homogeneous-cluster
  check described in its comparison summary before drawing conclusions
  from it.
