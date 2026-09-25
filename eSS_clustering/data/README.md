# data/

This directory is **not tracked in git** (see `.gitignore`).
Copy or symlink your input files here before running the pipeline.

---

## Required files by mutation type

### SBS — raw counts + pre-normalised (two files per species)

```
filtered_mouse_307.txt
normalized_filtered_mouse_307.tsv

filtered_rat_307.txt
normalized_filtered_rat_307.tsv

filtered_chicken_307.txt
normalized_filtered_chicken_307.tsv

filtered_celegans_307.txt
normalized_filtered_celegans_307.tsv

filtered_human_307.txt
normalized_filtered_human_307.tsv
```

### DBS — raw counts only (normalised internally)

```
filtered_mouse_DBS.txt
filtered_rat_DBS.txt
filtered_chicken_DBS.txt
filtered_celegans_DBS.txt
filtered_human_DBS.txt
```

### ID — raw counts only (normalised internally)

```
filtered_mouse_ID.txt
filtered_rat_ID.txt
filtered_chicken_ID.txt
filtered_celegans_ID.txt
filtered_human_ID.txt
```

---

## COSMIC reference profiles

Download from: https://cancer.sanger.ac.uk/signatures/downloads/

```
COSMIC_v3.6_SBS_GRCh38.txt       ← used for SBS
COSMIC_v3.6_DBS_GRCh38.txt       ← used for DBS
COSMIC_v3.6_ID_GRCh37.txt        ← used for ID
```

---

## File format

All input count files are tab-separated with:
- Rows = mutational contexts (96 for SBS, 78 for DBS, 83 for ID)
- Columns = sample names
- Values = raw mutation counts (integers)

Pre-normalised SBS files have the same shape but values sum to 1 per column.
