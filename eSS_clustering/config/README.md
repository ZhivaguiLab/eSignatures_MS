# config/

Contains the four configuration files that are tracked in git.
None of them contains patient data or sensitive information.

| File | Used by | Affects clustering? |
|------|---------|---------------------|
| `preprocessing.yaml` | `perform_clustering.py` | Yes — decides which samples are clustered |
| `sample_mapping.tsv` | `perform_clustering.py` | Names only |
| `compound_grouping.yaml` | `species_compound_matrix.py` | No — display only |
| `abv_table_clusters.txt` | reports, species-compound matrix | No — display only |

---

## `preprocessing.yaml`

Decides which samples are clustered, per mutation type:

- `exclude`: sample-name patterns (case-insensitive substrings). For SBS it
  excludes the 10 mouse MEF samples that are not part of the atlas (Xenon,
  Xenon XPA-/-, Deoxynivalenol, Deoxynivalenol + Patulin).
- `min_mutations` (SBS): the minimum total SBS count per species — 307 for
  every species (the Poisson-resampling stability threshold: 99% of
  simulations stable, over all samples). Every sample is
  tested, and every species in `data/input/SBS` must have a value. To try a
  different cutoff (e.g. per-species values), edit the numbers here; the next
  run rebuilds the preprocessed data. See "Changing the minimum mutation
  cutoff" in the main README.

For SBS, preprocessing also writes the normalized profiles. The cleaned data
is cached in `data/input_cleaned/<TYPE>/` and rebuilt automatically when the
settings or input files change; removed samples are listed in
`preprocessing_removed_samples.csv` there. See the main README's
Preprocessing section.

---

## `sample_mapping.tsv`

Maps original sample names (as they appear in input count files) to
standardised names used throughout the pipeline.

| Column | Description |
|--------|-------------|
| `sample_name` | Original sample name from input file |
| `standardized_name` | Standardised name with species prefix (e.g. `Mouse_Liver_AFB1`) |

---

## `abv_table_clusters.txt`

Tab-separated compound → acronym abbreviation table used in HTML reports
and the species-compound matrix to shorten long compound names.

| Column | Description |
|--------|-------------|
| `compound` | Full compound name (e.g. `Aflatoxin_B1`) |
| `acronym` | Short label used in plots (e.g. `AFB1`) |

---

## `compound_grouping.yaml`

Rules that collapse compound-name variants (replicate numbers, technical
suffixes) into one column of the species-compound matrix, plus explicit
code → name mappings (e.g. `ATC` → `5-aza-4-thio-2-deoxycytidine`,
`AAI` → `Aristolochic_acid_I`). Display only; it does not change clustering.
