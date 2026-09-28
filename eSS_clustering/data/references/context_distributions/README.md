# Trinucleotide opportunity tables

Per-chromosome trinucleotide counts for each genome, whole genome (`_96`) and
exome capture territory (`_96_exome`), used by the optional opportunity
normalization methods (`wes-to-wgs`, `own-opportunity`; see
`../../../NORMALIZATION_APPROACH.md`).

| File | Genome | Used for |
|------|--------|----------|
| `context_counts_GRCh38_96.csv` / `_exome.csv` | human GRCh38 | human WGS / WES samples |
| `context_counts_mm10_96.csv` / `_exome.csv` | mouse mm10 | mouse WGS / WES samples |
| `context_counts_rn7_96.csv` / `_exome.csv` | rat rn7 | rat samples (`own-opportunity` only) |

Copied unchanged from SigProfilerMatrixGenerator 1.3.6
(`SigProfilerMatrixGenerator/references/chromosomes/context_distributions/`),
BSD 2-Clause License, Copyright (c) 2019, Erik Bergstrom [Alexandrov Lab];
see `LICENSE_SigProfilerMatrixGenerator`.

No tables exist here for chicken or C. elegans, so those species are not
corrected by either method. The genome builds (GRCh38, mm10, rn7) are the
best available match and have not been confirmed against each source study's
alignment.
