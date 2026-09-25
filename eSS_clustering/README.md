# eSignatures Clustering Pipeline

Cross-species mutational signature clustering pipeline supporting SBS, DBS, and ID
mutation types. Developed in the Alexandrov Lab (UCSD) for experimental mutational
signature clustering.

---

## Repository structure

```
esignatures-pipeline/
  run_pipeline.sh              ← entry point — run this

  pipeline/                   ← all Python scripts
    utils/
      __init__.py
      naming.py               ← shared sample-name parsing
      mutation_type.py        ← per-type constants (contexts, artifacts, file suffixes)
    perform_clustering.py     ← Step 1: hierarchical clustering (with integrated preprocessing)
    generate_summary.py       ← Step 2: HTML reports
    generate_images_heatmap.py       ← Step 3: PDF → PNG
    generate_interactive_heatmap.py  ← Step 4: interactive heatmap
    generate_static_heatmap.py       ← Step 5: static heatmap
    species_compound_matrix.py       ← Step 6: species × compound matrix

  config/                     ← tracked in git
    sample_mapping.tsv        ← original → standardised sample name lookup
    compound_grouping.yaml    ← compound variant collapsing rules (optional)
    abv_table_clusters.txt    ← compound → acronym abbreviations

  data/                       ← NOT tracked (too large; see data/README.md)
    input/                    ← original data files
      SBS/
      DBS/
      ID/
    input_cleaned/            ← auto-generated preprocessed data (NOT tracked)
      DBS/
      ID/
    filtered_mouse_307.txt
    COSMIC_v3.6_SBS_GRCh38.txt
    ... (see data/README.md for full list)

  results/                    ← NOT tracked (generated at runtime)
    SBS/
    DBS/
    ID/

  requirements.txt
  .gitignore
```

---

## Installation

```bash
git clone https://github.com/<org>/eSignatures-clustering-analysis.git
cd eSignatures-clustering-analysis

conda create -n esig python=3.10
conda activate esig

# poppler required by pdf2image
brew install poppler          # macOS
# sudo apt-get install poppler-utils  # Linux

pip install -r requirements.txt
```

Then populate `data/` with your input files — see `data/README.md` for the
full list of required files and where to download COSMIC profiles.

---

## Running the pipeline

```bash
# From the repo root — no arguments needed beyond type and thresholds
bash run_pipeline.sh SBS 0.9 0.85
bash run_pipeline.sh DBS 0.9 0.85
bash run_pipeline.sh ID  0.9 0.85
```

The second argument is the cosine similarity threshold used for clustering
(`--cosine_similarity`); the third is the threshold used when comparing
resulting eSS profiles against COSMIC signatures (`--threshold` in Step 4). The
project's current convention is 0.9 for clustering and 0.85 for the COSMIC
comparison — adjust for your own analysis as needed.

Each run writes to `results/<MUTATION_TYPE>/` automatically.

### Reproducing the published clusters

The canonical SBS clustering uses:

| Setting | Value | Where it's defined |
|---|---|---|
| Cosine similarity threshold | `0.9` (distance `0.1`), average linkage | `--cosine_similarity` default in `perform_clustering.py` |
| Per-cluster custom thresholds | `Aristolochic_acid_I: 0.095`, `Dibenzo[a,l]pyrene: 0.095` (cosine distance) | `default_custom_thresholds` in `pipeline/utils/mutation_type.py` |
| Sample mapping | `config/sample_mapping.tsv` | `--mapping_file` default |

The custom thresholds re-split any main cluster that contains a sample
matching the pattern, using the tighter distance. They're applied by default,
so these two commands give **identical** clusters. With the current
`data/input/SBS` that's **49 main clusters, 16 small clusters and 131
singletons**:

```bash
bash run_pipeline.sh SBS 0.9 0.85
python pipeline/perform_clustering.py --mutation_type SBS --output_dir results
```

To override the custom thresholds, pass `--custom_thresholds 'pattern:value,...'`.
To turn them off, pass `--custom_thresholds none`, which gives 48 main clusters.

The clustering has no random step, so the same input and settings always give
the same clusters. Cluster IDs (`eSS1`, `eSS2`, …) are numbered left to right
along the dendrogram. Adding or removing a sample, or changing a threshold, can
renumber them. Compare runs by which samples are in each cluster, not by ID.

**Preprocessing is automatic** — if `config/preprocessing.yaml` exists and defines
exclusion patterns for the mutation type, samples will be filtered before clustering.
See [Preprocessing](#preprocessing) below.

### Consensus profile averaging method

`perform_clustering.py` supports two ways to compute each cluster's consensus
profile (`--averaging_method`, default `equal_replicate`). `run_pipeline.sh`
does not expose this as a pass-through argument, so to override it, call the
script directly:

```bash
python pipeline/perform_clustering.py \
    --mutation_type SBS \
    --output_dir results \
    --data_dir data/input/SBS \
    --averaging_method pooled
```

- **`equal_replicate`** (default): normalize each sample to sum to 1 first,
  then take the unweighted mean across the cluster's samples. Every
  sample/replicate contributes equally regardless of mutation burden — avoids
  one high-burden replicate dominating a cluster's apparent profile.
- **`pooled`**: sum raw counts across a cluster's samples, then normalize
  once. Samples with a higher mutation burden contribute proportionally more
  to the consensus profile.

Neither option affects cluster membership, only the reported consensus
profile used for plotting and COSMIC comparison.

### Dendrogram color count

`cluster_signatures_with_custom_thresholds()` in `perform_clustering.py` draws
one distinct color per colored cluster (main + small combined — true
singletons are always drawn gray, not from this palette) from a fixed-size
`seaborn` "husl" palette, set by the `num_colors` parameter (currently `100`,
not exposed via CLI — edit the function default directly). It only needs to
exceed however many colored clusters a run actually produces (main + small
clusters -- run `perform_clustering.py` and check the printed "Main clusters"
/ "Small clusters" counts, or count unique values in the `Color` column of
`main_clusters_summary.tsv` + `small_clusters_summary.tsv`); a *smaller*
`num_colors` spaces the palette's hues further apart for better visual
differentiation, as long as it still exceeds that count.

---

## Preprocessing (Optional)

The pipeline includes integrated preprocessing to filter out unwanted samples before
clustering. This is useful for removing experimental artifacts, controls, or other
samples that should not be included in the analysis.

### Quick Start

1. **Edit the preprocessing config:**

`config/preprocessing.yaml` is already tracked in git with the exclusion
patterns currently in use. To change them, edit it directly:

```yaml
# config/preprocessing.yaml
DBS:
  exclude:
    - hTumour
    - experimental
ID:
  exclude:
    - hTumour
```

2. **Run pipeline normally:**

```bash
bash run_pipeline.sh DBS 0.9 0.85
```

The pipeline will:
- Check if `config/preprocessing.yaml` exists
- If exclusion patterns are defined for DBS → preprocess and save to `data/input_cleaned/DBS/`
- If cleaned data already exists → reuse it (cached)
- If no patterns defined → use original data

### How It Works

**First run (preprocessing needed):**
```
data/input/DBS/human_data.txt (450 samples)
  ↓ preprocessing (removes 'hTumour' matches)
data/input_cleaned/DBS/human_data.txt (427 samples)
  ↓ clustering uses cleaned data
results/DBS/main_clusters/...
```

**Subsequent runs (cached):**
```
data/input_cleaned/DBS/ already exists
  ↓ skip preprocessing, use cached data
  ↓ clustering
results/DBS/main_clusters/...
```

### Configuration Format

```yaml
# config/preprocessing.yaml

DBS:
  exclude:
    - hTumour          # Substring match (case-insensitive)
    - experimental
    - _test

ID:
  exclude:
    - hTumour

# SBS: omit if no preprocessing needed
# Alternatively, use empty list:
# SBS:
#   exclude: []

case_sensitive: false  # Default: false
```

**Pattern Matching:**
- Patterns are **substrings**: `hTumour` matches `Mouse_hTumour_sample1`
- **Case-insensitive** by default: `hTumour`, `htumour`, `HTUMOUR` all match
- Samples matching **any** pattern are excluded

### Advanced Options

```bash
# Force reprocessing even if cleaned data exists
python pipeline/perform_clustering.py \
    --mutation_type DBS \
    --output_dir results \
    --force_preprocess

# Skip preprocessing even if config exists (use original data)
python pipeline/perform_clustering.py \
    --mutation_type DBS \
    --output_dir results \
    --skip_preprocessing

# Use custom preprocessing config location
python pipeline/perform_clustering.py \
    --mutation_type DBS \
    --output_dir results \
    --preprocessing_config my_exclusions.yaml
```

### What Gets Tracked in Git

- ✅ `config/preprocessing.yaml` — tracked (documents filtering decisions)
- ❌ `data/input_cleaned/` — NOT tracked (auto-generated, in `.gitignore`)

This means:
1. Collaborators pull your `config/preprocessing.yaml`
2. Pipeline auto-generates `data/input_cleaned/` on their machine
3. Everyone gets the same filtering

### Example Output

```
==========================================================
DATA DIRECTORY SELECTION
==========================================================
✓ Found preprocessing config: config/preprocessing.yaml
✓ Exclusion patterns for DBS: ['hTumour', 'experimental']

Checking for cached cleaned data...
✗ Not found: data/input_cleaned/DBS/

Running preprocessing...

==========================================================
PREPROCESSING
==========================================================
Input:  data/input/DBS
Output: data/input_cleaned/DBS
Exclusion patterns: ['hTumour', 'experimental']
==========================================================

  Processing: human_DBS_data.txt
    Original: 450 samples
    Removed:  23 samples
    Kept:     427 samples
    Matched patterns (first 3):
      - Human_hTumour_sample1 ('hTumour')
      - Human_hTumour_sample2 ('hTumour')
      - Mouse_experimental_1 ('experimental')

==========================================================
PREPROCESSING COMPLETE
  Total samples: 450
  Removed: 23
  Kept: 427
  Cleaned data saved to: data/input_cleaned/DBS
==========================================================

✓ Using cleaned data: data/input_cleaned/DBS

Loading DBS data from: data/input_cleaned/DBS
  human (human_DBS_data.txt): 427 samples loaded
  ...
```

---

## Species-Compound Matrix Visualization

The pipeline generates species-compound matrices for all mutation types (SBS, DBS, ID), providing a visual summary of how different experimental models cluster across various compound exposures.

### What It Shows

**Matrix Structure:**
- **Rows**: Species-model combinations (e.g., Mouse_Liver, Human_iPSC)
- **Columns**: Compound exposures (e.g., Cisplatin, Benzo[a]pyrene)
- **Cells**: Pie charts showing clustering outcomes for each combination

**Pie Chart Colors:**
- 🔵 **Blue** = Main cluster (>2 samples)
- 🟡 **Yellow** = Small cluster (n=2)
- 🔴 **Red** = Singleton (did not cluster)

### Output Files

Three separate PDFs are generated per mutation type:
```
results/<TYPE>/species_compound_matrix_human.pdf
results/<TYPE>/species_compound_matrix_mouse.pdf
results/<TYPE>/species_compound_matrix_other_species.pdf
```

### Compound Grouping (Optional)

By default, compound names may include replicate numbers or technical suffixes that create redundant columns. Use compound grouping to collapse these variants:

**Example problem:**
```
1_human_tumor_connective_SNV_hg38  ← 22 separate columns
2_human_tumor_connective_SNV_hg38
...
22_human_tumor_connective_SNV_hg38

B2D12_CX1, B2D12_CX2, B2D12_CX3, B2D12_CX4  ← 4 separate columns
```

**Solution:** Create `config/compound_grouping.yaml`:

```yaml
# Regex-based pattern matching
regex_rules:
  # Collapse numbered tumor samples
  - pattern: '^\d+_human_tumor_connective_SNV_hg38$'
    replace_with: 'Human_tumor_connective'
  
  # Remove replicate numbers from genotype variants
  - pattern: '^(B2D12|BRCA1|WT)_(CX|ETO|PDS|U)\d+$'
    replace_with: '\1_\2'
  
  # Remove replicate suffix
  - pattern: '^(.+)_replicate\d+_final$'
    replace_with: '\1'

# Explicit overrides for specific cases
explicit_mappings:
  CHF: 'Cyclophosphamide'
  BaP: 'Benzo[a]pyrene'
  AAI: 'Aristolochic_acid_I'
```

**Usage:**
```bash
python pipeline/species_compound_matrix.py \
    --main_output_dir results/DBS \
    --output_path results/DBS/species_compound_matrix.pdf \
    --mutation_type DBS \
    --compound_grouping config/compound_grouping.yaml
```

**Impact:** Can reduce 195 compound columns → 45 columns (75% reduction!)

### Compound Abbreviations

For short compound name abbreviations (without grouping variants), use `config/abv_table_clusters.txt`:

```tsv
compound	acronym
Benzo[a]pyrene	BaP
Cisplatin	Cis
4-Nitroquinoline-1-oxide	4NQO
```

**Note:** Compound grouping (for collapsing variants) and abbreviations (for shortening display names) serve different purposes and can be used together.

### Example Use Cases

**1. Identify which models cluster consistently:**
- Models with mostly blue cells → reliable clustering
- Models with mostly red cells → high variability

**2. Compare compound behavior across species:**
- Same compound, different species → species-specific effects
- Same compound, all species → universal signature

**3. Spot experimental outliers:**
- Single red cell in sea of blue → potential experimental issue

**4. Quality control for technical replicates:**
- Use compound grouping to verify replicates cluster together
- Red cells after grouping → genuine technical issues


---

## Pipeline steps and outputs

| Step | Script | Key outputs in `results/<TYPE>/` |
|------|--------|----------------------------------|
| 0 (optional) | *integrated preprocessing* | `data/input_cleaned/<TYPE>/` (cached, auto-generated) |
| 1 | `perform_clustering.py` | `main_clusters/`, `small_clusters/`, `singletons/`, `interactive_clustering/` |
| 2 | `generate_summary.py` | `reports/cluster_report.html`, `reports/cluster_summary.html`, `reports/main_clusters_with_dendrogram.html`, `reports/small_clusters_cluster_report.html`, `reports/small_clusters_summary.html` |
| 3 | `generate_images_heatmap.py` | `interactive_heatmap/cluster_plots/*.png`, `interactive_heatmap/cosmic_plots/*.png` |
| 4 | `generate_interactive_heatmap.py` | `interactive_heatmap.html`, `cosmic_similar_clusters_<t>.tsv`, `denovo_clusters_<t>.tsv` |
| 5 | `generate_static_heatmap.py` | `cosine_similarity_heatmap.pdf` (3 versions) |
| 6 | `species_compound_matrix.py` | `species_compound_matrix_{human,mouse,other_species}.pdf` |

### Generating a print-ready PDF of a grid report

`cluster_summary.html`, `main_clusters_with_dendrogram.html`, and
`small_clusters_summary.html` (all produced by Step 2, and all sharing the same
grid of mutational profile + pie chart + sample legend, at 3 columns) are meant
to be printed or exported to PDF for use as a manuscript figure. Use headless
Chrome rather than relying on a manual browser print — each report's CSS is
written specifically so that headless-Chrome print pagination lays the grid out
correctly (flexbox reflows row by row; CSS Grid does not paginate reliably and
will scatter columns across separate pages).

```bash
HTML_PATH="results/SBS/reports/cluster_summary.html"   # input — swap for the report/run you're rendering
PDF_PATH="results/SBS/reports/cluster_summary.pdf"      # output

"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
  --headless --disable-gpu --no-sandbox \
  --print-to-pdf="${PDF_PATH}" \
  --print-to-pdf-no-header \
  "file://$(cd "$(dirname "${HTML_PATH}")" && pwd)/$(basename "${HTML_PATH}")"
```

To convert the PDF to PNG (e.g. for embedding in a slide deck or manuscript figure),
one option is `pdftoppm` (Poppler — already a dependency for `pdf2image`, see
[Installation](#installation)):

```bash
pdftoppm -png -r 300 "${PDF_PATH}" "${PDF_PATH%.pdf}"
```

`pdftoppm` appends a page number to the given prefix, so a multi-page
`cluster_summary.pdf` produces `cluster_summary-1.png`, `cluster_summary-2.png`,
etc. — one PNG per page. `-r 300` sets the resolution in DPI; raise it for a
higher-resolution export. To combine multiple pages into a single figure,
assemble them manually in a layout tool (e.g. Adobe Illustrator) — this step
isn't scripted here.

On Linux, use `google-chrome --headless ...` instead of the macOS app path above.

---

## Mutation type differences

| | SBS | DBS | ID |
|--|-----|-----|----|
| Contexts | 96 | 78 | 83 |
| Cluster prefix | `eSS` | `eDS` | `eIS` |
| Normalisation | pre-computed file | internal (col sum → 1) | internal (col sum → 1) |
| Artifact signatures | SBS27, 43, 45–60, 95 | *(none)* | ID9 |
| Typical preprocessing | *(none)* | hTumour exclusion | hTumour exclusion |
| Matrix generation | ✓ | ✓ | ✓ |

All type-specific constants are in `pipeline/utils/mutation_type.py`.

---

## Trinucleotide Opportunity Normalization (Experimental)

Two independent methods for correcting cross-species and cross-technology
(WGS/WES) trinucleotide-composition bias — a manuscript-revision item — are
implemented and verified, but **not wired into the default pipeline**:
`pipeline/normalize_wes_to_wgs.py` (rescales WES samples onto their species'
WGS opportunity basis) and `pipeline/normalize_by_own_opportunity.py`
(normalizes every sample, WGS included, by its own opportunity table
independently, with no common reference genome). Both were run through the
full 6-step pipeline and compared against the default clustering; see
[NORMALIZATION_APPROACH.md](NORMALIZATION_APPROACH.md) for the method
details, rationale, results, and current limitations (unconfirmed genome
builds, no chicken/C. elegans opportunity tables yet).

---

## Adding a new mutation type

1. Add a `MutationTypeConfig` entry in `pipeline/utils/mutation_type.py`
2. Add input files to `data/input/<TYPE>/` following the naming convention in `data/README.md`
3. Add the COSMIC file path to the `COSMIC_FILES` array in `run_pipeline.sh`
4. (Optional) Add preprocessing patterns to `config/preprocessing.yaml`

No other scripts need changes.

--- 
## Input Data

The pipeline accepts raw mutation count matrices (and, for SBS, pre-normalised profiles) as tab-separated `.txt` or `.tsv` files placed in `data/input/<MUTATION_TYPE>/`. Files can follow any naming convention provided the species is identifiable somewhere in the filename — supported species are `human`, `mouse`, `rat`, `chicken`, and `c.elegans` (the latter is matched flexibly to account for variations such as `celegans`, `celegan`, `c_elegans`, etc.). Rather than requiring a fixed file naming scheme, the pipeline automatically discovers all `.txt`/`.tsv` files in the input directory, detects the species from each filename, and skips any files it cannot recognise with a warning. For SBS runs, pre-normalised files are distinguished from raw counts by the presence of `normaliz` in the filename; all other mutation types are normalised internally by scaling each sample's counts to sum to one. This means input directories do not need to contain the same set of species across runs — the pipeline will work with whatever files are present.

**Preprocessing:** If `config/preprocessing.yaml` defines exclusion patterns for a mutation type, samples matching those patterns will be automatically filtered out before clustering. The cleaned data is cached in `data/input_cleaned/<TYPE>/` for reuse across runs. See [Preprocessing](#preprocessing) for details.

---

## Troubleshooting

### Preprocessing Issues

**Cleaned data not being used:**
- Check if `config/preprocessing.yaml` exists
- Verify patterns are defined for your mutation type (DBS/ID/SBS)
- Look for "DATA DIRECTORY SELECTION" output to see what was detected

**Want to regenerate cleaned data:**
```bash
python pipeline/perform_clustering.py --mutation_type DBS --output_dir results --force_preprocess
```

**Want to use original data (skip preprocessing):**
```bash
python pipeline/perform_clustering.py --mutation_type DBS --output_dir results --skip_preprocessing
```

**Check what was filtered:**
- Preprocessing logs show exactly which samples were removed and why
- Compare file sizes: `ls -lh data/input/DBS/` vs `ls -lh data/input_cleaned/DBS/`

### General Issues

**No files found:**
- Verify files are in `data/input/<MUTATION_TYPE>/`
- Check filenames contain species keywords (human, mouse, rat, chicken, celegans)
- For SBS: ensure normalized files contain 'normaliz' in filename

**Clustering errors:**
- Check data dimensions match expected contexts (96/78/83)
- Verify no all-zero samples (auto-removed with warning)
