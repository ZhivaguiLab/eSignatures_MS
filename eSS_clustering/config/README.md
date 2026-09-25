# config/

Contains the two configuration files that are tracked in git.
Neither file contains patient data or sensitive information.

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
