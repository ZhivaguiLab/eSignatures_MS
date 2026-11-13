# eSignatures Clustering Pipeline

Cross-species mutational signature clustering pipeline supporting SBS, DBS, and ID
mutation types. Developed in the Alexandrov Lab (UCSD) for experimental mutational
signature clustering.

---

## Quick start

```bash
cd eSS_clustering
conda activate esignatures                 # see Installation
python pipeline/html_to_pdf.py --check     # Poppler (required) and Google Chrome (for PDFs)
python -m unittest discover tests          # 50 tests; should end with OK

bash run_background.sh SBS 0.9 0.85          # default (equal_replicate) → results/min307/SBS/
bash run_background.sh SBS 0.9 0.85 pooled   # pooled consensus      → results/min307_pooled/SBS/
tail -f results/logs/<log printed above>     # ends with "Pipeline complete!"
```

Start the second run after the first finishes (they share a preprocessing
cache). `run_pipeline.sh` takes the same arguments and runs in the foreground.

All figures of a run (cluster summaries, species × compound matrices, COSMIC
heatmaps, dendrogram, signature profiles, decomposition) end up as PDF +
300 DPI PNG in **`results/<run name>/SBS/final_output/`**, with
`run_parameters.txt`.

### Changes in November 2026

- **New inputs (data version 2026-11-13):** the collaborator's human file adds
  70 samples (694 profiles clustered, 50 eSS, 26/24 COSMIC matched/unmatched).
  New sample names are mapped in `config/sample_mapping.tsv`, genotypes
  collapsed into cell line + compound (the raw names keep the genotype).
  See [Input data versions](#input-data-versions).
- **`pooled` option:** pooled consensus profiles alongside the default
  `equal_replicate`, in their own run folder. See
  [Consensus profile averaging method](#consensus-profile-averaging-method).
- **Pie charts:** slices and legends ordered from the largest contributor down.
- **Figures made automatically:** one-page PDFs of the cluster summaries
  (Step 2b), the combined species × compound figure (Step 6b), a 300 DPI PNG of
  each, and all figures collected in `final_output/` (Step 9).
- **`run_background.sh`:** runs the pipeline in the background with a log.
- **Tests by data version:** expected results are stored per dated input
  release and the tests pass on either.

---

## Repository structure

```
eSS_clustering/
  run_pipeline.sh              ← entry point — run this
  run_background.sh            ← same, in the background with a log in results/logs/

  pipeline/                   ← all Python scripts
    utils/
      __init__.py
      naming.py               ← shared sample-name parsing
      mutation_type.py        ← per-type constants (contexts, artifacts, custom-threshold presets)
    perform_clustering.py     ← Step 1: hierarchical clustering (with integrated preprocessing)
    generate_summary.py       ← Step 2: HTML reports
    html_to_pdf.py            ← Step 2b: one-page PDFs of the summary reports
    generate_images_heatmap.py       ← Step 3: PDF → PNG
    generate_interactive_heatmap.py  ← Step 4: interactive heatmap
    generate_static_heatmap.py       ← Step 5: static heatmap
    species_compound_matrix.py       ← Step 6: species × compound matrix
    combined_matrix_figure.py        ← Step 6b: the three matrices as one figure (HTML + PDF)
    run_info.py                      ← run folder naming + run_parameters.txt (Step 7)
    collect_final_outputs.py         ← copies the figure PDFs/PNGs into final_output/ (Step 9)
    run_decomposition.py             ← COSMIC decomposition into eSS, SigProfilerAssignment (Step 8)
    utils/opportunity_normalization.py ← optional wes-to-wgs / own-opportunity normalization
    normalize_wes_to_wgs.py          ← standalone: wes-to-wgs copy of a preprocessed folder
    normalize_by_own_opportunity.py  ← standalone: own-opportunity copy of a preprocessed folder
    build_translation_verification_matrices.py ← before/after check for wes-to-wgs
    archive/                  ← superseded scripts, kept for reference

  config/                     ← tracked in git (see config/README.md)
    preprocessing.yaml        ← samples excluded before clustering
    sample_mapping.tsv        ← original → standardised sample name lookup
    compound_grouping.yaml    ← compound variant collapsing rules
    abv_table_clusters.txt    ← compound → acronym abbreviations

  data/                       ← see data/README.md
    input/SBS/                ← unfiltered SBS count matrices, one per species (tracked)
    input/DBS/, input/ID/     ← not included; add your own
    references/               ← COSMIC v3.6 SBS/DBS/ID profiles (tracked)
      context_distributions/  ← genome/exome trinucleotide counts for normalization (tracked)
    input_cleaned/            ← auto-generated preprocessed data (NOT tracked)

  results/<run name>/<TYPE>/  ← NOT tracked; one folder per run, named by its settings

  tests/
    test_clustering_reproducibility.py  ← checks a run reproduces the expected clusters
    data_versions.py                    ← which input data version is in data/input/SBS
    make_expected.py                    ← registers a new input data version + its expected results
    expected/data_versions.json         ← each input data version: date, what changed, fingerprints, expected numbers
    expected/<YYYY-MM-DD>/              ← expected sample → cluster assignments for that version
                                          (default, aai-split, per-species)
    expected/published_307_input_samples.tsv  ← samples in the earlier ≥307 input files
    expected/main_SBS_cluster_membership.tsv  ← earlier clustering on main (old ≥307 inputs)
    expected/main_SBS_eSS_profiles.tsv        ← earlier eSS consensus profiles on main

  NORMALIZATION_APPROACH.md   ← experimental normalization methods and results
  docs/normalization/         ← normalization runs compared with the default run
  requirements.txt            ← direct dependencies (pinned)
  requirements-lock.txt       ← full locked environment for exact reproduction
```

---

## Installation

```bash
git clone https://github.com/ZhivaguiLab/eSignatures_MS.git
cd eSignatures_MS/eSS_clustering

conda create -n esignatures python=3.11
conda activate esignatures

# Poppler (required: pdf2image, the PNGs and the combined figure use it)
brew install poppler          # macOS
# sudo apt-get install poppler-utils  # Linux
# conda install -c conda-forge poppler  # or inside the conda environment

# Exact environment the results were verified in (recommended)
pip install -r requirements-lock.txt
```

Two programs outside Python are used:

- **Poppler** (required). Not a Python package, so `pip install` does not
  install it; use the commands above. Without it the run stops at start-up.
- **Google Chrome** (optional), for the one-page PDFs and PNGs of the cluster
  summaries (Step 2b) and the combined species × compound figure (Step 6b).
  Without it the run continues and those files are skipped. Set `CHROME_PATH`
  if Chrome is not in the usual place.

Every run checks for both first and says so at the top of its log. To check
after creating the environment:

```bash
python pipeline/html_to_pdf.py --check
```

The pipeline requires Python 3.11 and stops at start-up on any other version.
There are two requirements files:

- `requirements-lock.txt` pins every package, including indirect
  dependencies, to the versions that reproduce the expected clusters. Use
  this to reproduce the manuscript results.
- `requirements.txt` pins only the packages the pipeline imports directly.
  Edit this one when updating a dependency, then regenerate the lock file
  (instructions at the top of `requirements-lock.txt`).

Both include SigProfilerAssignment 1.1.4 for the COSMIC decomposition
(Step 8). Other versions give different decomposition output (0.2.0, for
example, renames the COSMIC signatures), so `run_decomposition.py` stops if a
different version is installed.

After installing, check the install reproduces the expected clusters:

```bash
python -m unittest discover tests
```

The SBS input matrices and COSMIC v3.6 references used for the manuscript are
included in `data/`, so the SBS pipeline runs as soon as the environment is
installed. DBS and ID inputs are not included; see `data/README.md` for the
expected files.

---

## Running the pipeline

> **Default analysis (SBS)** — `bash run_pipeline.sh SBS 0.9 0.85`
>
> - Minimum **307 SBSs per sample, for every species** (applied to every sample, whatever its name)
> - 10 mouse MEF samples excluded by name (Xenon, Deoxynivalenol)
> - Clustering: cosine similarity **0.90**, average linkage; no manual split
> - COSMIC match: max cosine similarity **≥ 0.85** (artifact signatures included)
> - COSMIC decomposition: SigProfilerAssignment 1.1.4, novelty threshold **0.8**
> - No opportunity normalization
>
> Output: `results/min307/SBS/` — 694 profiles, 50 eSS, 26/24 COSMIC
> matched/unmatched (input data version 2026-11-13; see
> [Input data versions](#input-data-versions)). The tests check this default.

```bash
# From eSS_clustering/ (run_pipeline.sh also works from any other directory)
bash run_pipeline.sh SBS 0.9 0.85
bash run_pipeline.sh DBS 0.9 0.85
bash run_pipeline.sh ID  0.9 0.85

# Optional extras (SBS), in any combination:
bash run_pipeline.sh SBS 0.9 0.85 min=250            # another cutoff, every species
bash run_pipeline.sh SBS 0.9 0.85 per-species        # per-species cutoffs
bash run_pipeline.sh SBS 0.9 0.85 aai-split          # testing: manual AAI/DBP split
bash run_pipeline.sh SBS 0.9 0.85 wes-to-wgs         # opportunity normalization
bash run_pipeline.sh SBS 0.9 0.85 own-opportunity    # opportunity normalization
bash run_pipeline.sh SBS 0.9 0.85 pooled             # pooled consensus profiles
```

See [Changing the minimum mutation cutoff](#changing-the-minimum-mutation-cutoff)
for `min=<N>` and `per-species`, and
[Consensus profile averaging method](#consensus-profile-averaging-method) for
`pooled`.

### Running in the background

`run_background.sh` takes the same arguments as `run_pipeline.sh` and runs it
in the background (so the terminal can be closed), in the `esignatures` conda
environment, with a timestamped log in `results/logs/`:

```bash
bash run_background.sh SBS 0.9 0.85                 # default
bash run_background.sh SBS 0.9 0.85 pooled          # pooled consensus profiles
tail -f results/logs/SBS_0.9_0.85_pooled_*.log      # follow a run
```

The log ends with `Pipeline complete!` when the run succeeds. To use another
conda environment, set `ESIG_ENV`, e.g. `ESIG_ENV=myenv bash run_background.sh ...`.
Runs with different options write to different folders, but they share a
preprocessing cache when their cutoff and normalization are the same (e.g.
the default and `pooled`), so start the second one after the first finishes.

The second argument is the cosine similarity threshold used for clustering
(`--cosine_similarity`); the third is the threshold used when comparing
resulting eSS profiles against COSMIC signatures (`--threshold` in Step 4). The
project's current convention is 0.9 for clustering and 0.85 for the COSMIC
comparison — adjust for your own analysis as needed.

### Final figures: `final_output/`

The last step copies every figure of the run into one folder,
`results/<run name>/<TYPE>/final_output/` (e.g. `results/min307/SBS/final_output/`
and `results/min307_pooled/SBS/final_output/`), to share or to compare runs:

| Folder | Files |
|---|---|
| `cluster_summaries/` | `cluster_summary`, `main_clusters_with_dendrogram`, `small_clusters_summary` (.pdf + 300 DPI .png) |
| `species_compound_matrix/` | `species_compound_matrix_combined` (.pdf + .png) and the `_human`, `_mouse`, `_other_species` matrices |
| `cosine_similarity/` | the three eSS vs COSMIC heatmaps |
| `dendrogram/` | `all_esignature_models_dendrogram` (.pdf + .png) |
| `signature_profiles/` | eSS and small-cluster signature plots |
| `decomposition/` | `SBS96_Decomposition_Plots.pdf` (SBS) |
| `run_parameters.txt` | the run's settings and results |
| `MANIFEST.txt` | each file and where in the run folder it was copied from |

These are copies; the originals stay where the other steps wrote them. The
folder is replaced on every run. An expected file that was not made (e.g. the
summary PDFs without Google Chrome) is listed in `MANIFEST.txt` and in the log
as a warning. To collect an existing run folder:

```bash
python pipeline/collect_final_outputs.py --run_dir results/min307/SBS
```

### Output folders and run parameters

Each run writes to `results/<run name>/<MUTATION_TYPE>/`, where the run name
records the minimum mutation count and the optional extras, so those runs
never overwrite each other:

| Command | Output folder |
|---|---|
| `bash run_pipeline.sh SBS 0.9 0.85` | `results/min307/SBS/` |
| `... min=250` | `results/min250/SBS/` |
| `... per-species` | `results/min-per-species/SBS/` |
| `... aai-split` | `results/min307_aai-split/SBS/` |
| `... wes-to-wgs` | `results/min307_wes-to-wgs/SBS/` |
| `... own-opportunity` | `results/min307_own-opportunity/SBS/` |
| `... pooled` | `results/min307_pooled/SBS/` |

`min307` is the minimum mutation count per sample (`min<N>` for `min=<N>`,
`min-per-species` for `per-species`, `nomin` if none is set). Options
combine, e.g. `per-species wes-to-wgs` → `results/min-per-species_wes-to-wgs/`
and `per-species pooled` → `results/min-per-species_pooled/`. The
clustering cosine threshold (0.9) is not in the name, to avoid confusing it
with the COSMIC match threshold; it is recorded in `run_parameters.txt`. Runs
that differ only in the clustering threshold write to the same folder, so
move or rename a run folder before rerunning with another threshold.

Every run folder has **`run_parameters.txt`** (and the same in
`run_parameters.json`) listing:

- preprocessing: minimum mutations per species, samples excluded by name,
  opportunity normalization
- clustering: cosine threshold, linkage, custom thresholds, main-cluster size,
  consensus method
- COSMIC comparison: reference file, match threshold, artifacts included
- results: samples clustered, main / small clusters, singletons, COSMIC
  matched / unmatched
- SHA-256 of every input file, opportunity table and COSMIC reference
- the git commit the code ran from (and whether tracked files had uncommitted
  changes; untracked files are ignored),
  the date, and Python/package versions

For SBS, the run folder also has `decomposition/` (Step 8, below), with its
own `decomposition_parameters.txt`.

### COSMIC decomposition (Step 8)

Step 8 of every SBS run decomposes each COSMIC v3.6 SBS signature into
that run's eSS with SigProfilerAssignment 1.1.4 (`decompose_fit`), writing to
`results/<run name>/SBS/decomposition/`:

- `Decompose_Solution/De_Novo_map_to_COSMIC_SBS96.csv`: for each COSMIC
  signature, the eSS combination that reconstructs it and the reconstruction
  cosine similarity (the input for the Figure 5 decomposition plots)
- `Decompose_Solution/SBS96_Decomposition_Plots.pdf`, `Signatures/`,
  `Activities/`, `Solution_Stats/`: the rest of SigProfilerAssignment's output
- `decomposition_parameters.txt` / `.json`: settings, input hashes, package
  versions, and how many COSMIC signatures the eSS reconstruct

Settings (`pipeline/run_decomposition.py`):

- **Novelty threshold 0.8** (SigProfilerAssignment's default): a COSMIC
  signature whose best eSS reconstruction has cosine similarity below 0.8
  stays as itself, i.e. is not explained by the eSS. This is a different
  question from the ≥0.85 COSMIC match of a single eSS, so the numbers
  need not be the same. With the default run (2026-11-13 inputs), 50
  of the 101 COSMIC signatures are reconstructed from eSS.
- **Samples:** the COSMIC signatures themselves (as in
  `SPA_1000_shuffled_libraries/null_library_fit.py`). The COSMIC-to-eSS map
  does not depend on the samples; only the activities do.
- Other settings are SigProfilerAssignment's defaults, stated explicitly
  (`collapse_to_SBS96=False`, `connected_sigs=False`, NNLS penalties
  0.05 / 0.01 / 0.05).

To rerun it on its own, or with other choices:

```bash
python pipeline/run_decomposition.py --run_dir results/min307/SBS
python pipeline/run_decomposition.py --run_dir results/min307/SBS \
    --threshold 0.85 --samples clustered --output results/min307/SBS/decomposition_085
```

`--samples clustered` assigns the samples clustered in that run instead of
the COSMIC signatures.

### Reproducing the clusters

The SBS clustering uses:

| Setting | Value | Where it's defined |
|---|---|---|
| Minimum mutations per sample | 307 SBSs for every species (Poisson-resampling stability threshold: 99% of simulations stable, over all samples), applied to every sample | `SBS: min_mutations` in `config/preprocessing.yaml`; options `min=<N>` / `per-species` |
| Excluded samples | 10 mouse MEF samples (Xenon, Deoxynivalenol) | `SBS: exclude` in `config/preprocessing.yaml` |
| Cosine similarity threshold | `0.9` (distance `0.1`), average linkage | `--cosine_similarity` default in `perform_clustering.py` |
| Per-cluster custom thresholds | none (the AAI/DBP split is an optional testing preset, see below) | `default_custom_thresholds` / `custom_threshold_presets` in `pipeline/utils/mutation_type.py` |
| Main vs small clusters | main: ≥3 samples; small: 2 samples (`MEF_AID` and `MCF10_cisplatin` 2-sample clusters are kept as main) | `special_patterns` in `perform_clustering.py` |
| Sample mapping | `config/sample_mapping.tsv` | `--mapping_file` default |
| Consensus profile | `equal_replicate` (option `pooled`, see below) | `--averaging_method` default in `perform_clustering.py`; `pooled` option in `run_pipeline.sh` |
| COSMIC match | max cosine similarity ≥ `0.85` to any COSMIC v3.6 SBS signature, artifact signatures included | `MATCH_THRESHOLD` in `generate_static_heatmap.py`; `--threshold` in `generate_interactive_heatmap.py` |

These two commands give **identical** clusters. With the current
`data/input/SBS` (version 2026-11-13), 694 of the 1,552 bundled profiles pass
preprocessing and are clustered into **50 main clusters, 18 small clusters and
130 singletons**:

```bash
bash run_pipeline.sh SBS 0.9 0.85
python pipeline/perform_clustering.py --mutation_type SBS --output_dir results
```

**25 of the 48 main clusters match COSMIC (≥0.85); 23 don't.** Both heatmap
scripts report the same split. Four clusters sit just under the cutoff (eSS7,
eSS20, eSS22, eSS25 at 0.847–0.848) and are counted as not matched, even though
the static heatmap shows their value rounded to "0.85".

### Testing option: the manual AAI/DBP split

The `aai-split` preset re-splits the cluster that contains the Aristolochic
acid I and Dibenzo[a,l]pyrene samples at a tighter cosine distance (0.095). It
is not part of the default analysis; use it to test how that split changes the
results:

```bash
bash run_pipeline.sh SBS 0.9 0.85 aai-split        # → results/min307_aai-split/SBS/
python pipeline/perform_clustering.py --mutation_type SBS --output_dir results/min307_aai-split \
    --custom_thresholds aai-split
```

| | Default | `aai-split` |
|---|---|---|
| Main clusters (eSS) | 50 | 51 |
| Small clusters / singletons | 18 / 130 | 18 / 130 |
| COSMIC matched / unmatched (≥0.85) | 26 / 24 | 27 / 24 |

With the split, the 11-sample AAI/DBP cluster becomes two: 8 samples (6
Aristolochic acid I + 2 Dibenzo[a,l]pyrene/DBPDE, so still mixed) and 3 DBPDE
samples. Every other cluster is identical. The preset's thresholds are defined
once, in `custom_threshold_presets` in `pipeline/utils/mutation_type.py`. Other
thresholds can be passed as `--custom_thresholds 'pattern:value,...'`.

### Changing the minimum mutation cutoff

**Default: 307 SBSs per sample for every species.** Two options change it
without editing any file:

| Cutoff | Command | Output folder | Clustered profiles |
|---|---|---|---|
| **307 for every species (default)** | `bash run_pipeline.sh SBS 0.9 0.85` | `results/min307/` | 694 |
| One number for every species | `bash run_pipeline.sh SBS 0.9 0.85 min=250` | `results/min250/` | depends on N |
| Per species | `bash run_pipeline.sh SBS 0.9 0.85 per-species` | `results/min-per-species/` | 715 |

With `perform_clustering.py` directly, use `--min_mutations 250` or
`--min_mutations per-species`.

The numbers live in `config/preprocessing.yaml`:

```yaml
SBS:
  min_mutations: 307                # DEFAULT: one cutoff for every species
  min_mutations_per_species:        # used only with the 'per-species' option
    mouse:    235
    human:    295
    celegans: 451
    chicken:  300
    rat:      6354    # only 5 rat samples; equals the smallest, so all 5 are kept
```

307 is the Poisson-resampling stability threshold ("SBS for 99% of
simulations") over all samples; the per-species values are the same threshold
computed from each species' own samples. `per-species` gives 715 profiles,
50 eSS, 17 small clusters, 136 singletons and 26/24 COSMIC matched/unmatched
(2026-11-13 inputs). The tests check its sample set, counts and cluster
membership (`tests/expected/<version>/SBS_cluster_membership_per_species.tsv`).

Notes:

- Each cutoff gets its own preprocessing cache (`data/input_cleaned/SBS/` for
  the default, `SBS_min250/`, `SBS_min-per-species/`, …), so runs don't
  overwrite each other's data. `min=307` is the default and uses the default
  cache and folder.
- The cutoffs used are listed in each run's `run_parameters.txt`
  ("Minimum mutations per sample", with "config default" or the option).
- To try other per-species values, edit `min_mutations_per_species` (every
  species in `data/input/SBS` must be listed) and run with `per-species`.
  Editing `min_mutations` itself changes the default; the tests (written for
  307) will then fail and list the differences. Prefer `min=<N>`.

### Cluster numbering

The clustering has no random step, so the same input and settings always give
the same clusters. Cluster IDs (`eSS1`, `eSS2`, …) are numbered left to right
along the dendrogram. When a custom threshold splits a cluster, the pieces keep
that cluster's place and are numbered largest first. Adding or removing a
sample, or changing a threshold, can renumber them. Compare runs by which
samples are in each cluster, not by ID.

### How clusters are defined

1. **Distances.** Cosine distance (1 − cosine similarity) between every pair
   of normalized SBS96 profiles.
2. **Tree.** Average-linkage hierarchical clustering of those distances
   (`scipy.cluster.hierarchy.linkage`).
3. **Cut.** The tree is cut at cosine distance 0.1 with
   `scipy.cluster.hierarchy.fcluster`. Two samples are in the same cluster if
   the tree joins them below that height: every merge inside a cluster joins
   two groups whose average pairwise cosine distance is below 0.1 (similarity
   above 0.9). Membership comes straight from the tree; nobody picks clusters
   by hand.
4. **Custom thresholds (optional, off by default).** With `aai-split`, the
   cluster containing the AAI or DBP samples is cut again, on its own, at
   distance 0.095.
5. **Group.** Clusters with ≥3 samples are main clusters (eSS), 2-sample
   clusters are small clusters, 1-sample clusters are singletons.

The dendrogram colours are only for the figure. The code checks that every
cluster from the tree cut is drawn in exactly one colour and every singleton
is drawn gray, and stops with an error if the figure and the clusters ever
disagree.

The log reports the closest merges on either side of the cutoff. For the
current data they are 0.09851 (joined) and 0.10039 (not joined), so a cutoff
anywhere between cosine similarity 0.8996 and 0.9015 gives the same clusters.

### Checking a run reproduces the expected clusters

```bash
python -m unittest discover tests -v
```

The tests find which [input data version](#input-data-versions) is in
`data/input/SBS`, run preprocessing and clustering on it, and check the result
against that version's `tests/expected/<version>/SBS_cluster_membership.tsv`
(every sample's cluster and eSS number). They also check:

- the profile and cluster counts and COSMIC split recorded for the version
  (see the table under [Input data versions](#input-data-versions)), and for the `aai-split` testing
  option its counts, COSMIC split and own expected membership
  (`SBS_cluster_membership_aai_split.tsv`), which differs from the default
  only in the AAI/DBP cluster
- preprocessing: every clustered sample meets its species cutoff; each cleaned
  file equals an independent re-filter of the unfiltered file (and each
  normalized file equals its counts / total); every input sample is either kept
  or listed in the removed-samples log; a cached rerun is identical
- the cutoff options: `min=<N>` applies N to every species (checked against an
  independent re-filter), `min=307` is the default, and `per-species` gives
  the version's profile and cluster counts and the membership in
  `SBS_cluster_membership_per_species.tsv`; other cutoffs use
  their own cache and leave the default's alone
- on the 2026-09-28 inputs (rebuilt from the current files, see below): with
  307 for every species, preprocessing keeps exactly the samples in the
  earlier published `filtered_*_307.txt` inputs
  (`tests/expected/published_307_input_samples.tsv`), except the 18 C. elegans
  CX-5461 samples below 307 that the old filtering script let through
- **agreement with the earlier clustering on `main`** (which used the AAI/DBP
  split, so these tests use `aai-split`), with 307 for every species:
  - adding back the 18 CX-5461 samples the old script never tested gives
    exactly `main`'s inputs, cluster membership, eSS numbering and eSS
    profiles (671 profiles, 49 / 16 / 131)
  - without them, the result equals `main` with only those 18 samples removed:
    48 eSS with identical members and profiles, the CX-5461 eSS reduced to its 3
    samples ≥307, identical small clusters, and 26 / 23 COSMIC matches
- the optional opportunity normalizations: `wes-to-wgs` and
  `own-opportunity` profiles equal an independent calculation from the
  opportunity tables, keep the same samples, and leave uncorrected species
  unchanged
- run folder naming
- the mouse exclusions, and that the safety checks below stop the run on bad
  input

If you change the settings on purpose, regenerate the expected results
(`python tests/make_expected.py --version <date> --description "..." --force`)
and review the diff. For new input data, register a new version (below).

### Input data versions

Each release of the SBS input matrices in `data/input/SBS` is a **data
version**, named by the date the inputs changed and registered in
`tests/expected/data_versions.json` with what changed, a fingerprint of each
input file's contents and the expected results. Its expected cluster
membership files are in `tests/expected/<date>/`.

| Version | What changed | Profiles (≥307) | eSS / small / singletons | COSMIC matched / unmatched | Per-species: profiles, eSS / small / singletons |
|---|---|---|---|---|---|
| 2026-09-28 | original inputs (manuscript submission), 486 human samples | 653 | 48 / 16 / 123 | 25 / 23 | 674, 48 / 15 / 129 |
| **2026-11-13** (current) | collaborator's new human file: +70 human samples (TK6 4NQO/cisplatin ± POLH/I/K, XP30RO UVA/N-acetylcysteine, hTERT-RPE1 cisplatin P53-/-), 556 human samples; other species unchanged | 694 | 50 / 18 / 130 | 26 / 24 | 715, 50 / 17 / 136 |

The tests detect the version from the fingerprints, so they pass with either
set of inputs. If `data/input/SBS` matches no registered version, the
`DataVersionTest` fails with instructions and the version-specific checks are
skipped. The fingerprint is of the matrices themselves (sample names,
contexts, counts), not the file bytes, so line endings or number formatting
don't matter.

The checks against earlier published results (the old `filtered_*_307.txt`
inputs and the clustering on `main`) always use the 2026-09-28 inputs. New
releases so far only add samples after the existing ones (the 2026-11-13
human file is the 2026-09-28 file with 70 samples appended, counts
unchanged), so the tests rebuild the 2026-09-28 inputs from the current
files and confirm the rebuild with its fingerprint.

When the collaborator sends new inputs:

1. Copy them into `data/input/SBS/` and add any new sample names to
   `config/sample_mapping.tsv`.
2. Register them, dated the day they changed:
   ```bash
   python tests/make_expected.py --version 2027-01-15 --description "what changed"
   ```
   This writes `tests/expected/2027-01-15/` and adds the version, with its
   numbers, to `data_versions.json`.
3. Check the printed numbers against a full run on the same inputs, diff the
   new membership files against the previous version's, add a row to the
   table above, and commit the inputs and the new version together.

### Safety checks

The run stops with an error, rather than continuing, if:

- a sample is in the raw counts but not the normalized profiles (or the
  reverse), a sample name is duplicated, or a normalized profile doesn't equal
  its own counts divided by their total
- a cluster member can't be found when building the consensus profile
- the tree cut and the dendrogram colours disagree, or a sample ends up in
  more than one cluster or in none
- preprocessing fails (it used to fall back to the unfiltered data)
- `min_mutations` is set but a species in the input has no cutoff, or a loaded
  sample is below its species' cutoff

**Preprocessing is automatic** — if `config/preprocessing.yaml` exists and defines
exclusion patterns for the mutation type, samples will be filtered before clustering.
See [Preprocessing](#preprocessing) below.

### Consensus profile averaging method

Each eSS's consensus profile can be computed in two ways:

- **`equal_replicate`** (default): normalize each sample to sum to 1 first,
  then take the unweighted mean across the cluster's samples. Every
  sample/replicate contributes equally regardless of mutation burden — avoids
  one high-burden replicate dominating a cluster's apparent profile.
- **`pooled`**: sum raw counts across a cluster's samples, then normalize
  once. Samples with a higher mutation burden contribute proportionally more
  to the consensus profile.

To compare the two, run the pipeline once with each. The `pooled` run writes
to its own folder, so it does not overwrite the default run:

| Method | Command | Output folder |
|---|---|---|
| `equal_replicate` (default) | `bash run_pipeline.sh SBS 0.9 0.85` | `results/min307/SBS/` |
| `pooled` | `bash run_pipeline.sh SBS 0.9 0.85 pooled` | `results/min307_pooled/SBS/` |

`pooled` combines with the other options (e.g. `per-species pooled` →
`results/min-per-species_pooled/SBS/`), and `run_background.sh` takes it the
same way. The method used is listed in each run's `run_parameters.txt`
("Consensus profile").

The method does not change cluster membership. It changes the consensus
profiles (`main_clusters/main_clusters_NORMALIZED_weighted_avg_profiles.tsv`)
and everything computed from them: the eSS signature plots, the COSMIC
heatmaps and matched/unmatched counts, and the COSMIC decomposition (Step 8).

With `perform_clustering.py` directly, use `--averaging_method pooled`.

### Dendrogram color count

`cluster_signatures_with_custom_thresholds()` in `perform_clustering.py` draws
one distinct color per colored cluster (main + small combined — true
singletons are always drawn gray, not from this palette) from a
`seaborn` "husl" palette of `num_colors` colours (default `100`, not exposed
via CLI). Colours only affect the figure and the `Color` column of the
summary TSVs, never cluster membership. If a run produces 100 or more
clusters, the palette is enlarged automatically so no two clusters share a
colour.

---

## Preprocessing

The pipeline includes integrated preprocessing to filter out unwanted samples before
clustering. This is useful for removing experimental artifacts, controls, or other
samples that should not be included in the analysis. For SBS it is part of the
analysis: `data/input/SBS` holds the **unfiltered** profiles, and preprocessing

1. removes the 10 excluded mouse MEF samples (by name), and
2. removes every sample whose total SBS count is below the cutoff
   (default `min_mutations`: 307 SBSs for every species, the Poisson-resampling
   stability threshold; see
   [Changing the minimum mutation cutoff](#changing-the-minimum-mutation-cutoff)), and
3. writes the normalized profiles (each sample's counts divided by its total).

Every sample is tested against its cutoff, whatever its name. (The notebook
that made the earlier `filtered_*_307.txt` inputs only tested samples whose
name contained "exome" or "genome", so 21 C. elegans CX-5461 samples were never
filtered.) Skipping preprocessing does not reproduce the results, and SBS won't
run without it because the normalized files are only made here.

### Quick Start

1. **Edit the preprocessing config:**

`config/preprocessing.yaml` is already tracked in git with the exclusion
patterns currently in use. For SBS these remove 10 mouse MEF samples (Xenon,
Deoxynivalenol) that are not part of the atlas. To change them, edit it directly:

```yaml
# config/preprocessing.yaml
SBS:
  exclude:
    - Xenon
    - Deoxynivalenol
  min_mutations: 307          # one number for every species, or {species: number}
  min_mutations_per_species:  # only with the 'per-species' option
    mouse: 235
    ...
DBS:
  exclude:
    - hTumor
ID:
  exclude:
    - hTumor
```

2. **Run pipeline normally:**

```bash
bash run_pipeline.sh DBS 0.9 0.85
```

The pipeline will:
- Check if `config/preprocessing.yaml` exists
- If exclusion patterns are defined for DBS → preprocess and save to `data/input_cleaned/DBS/`
- If cleaned data already exists and was built from the same input files and
  patterns → reuse it (cached); otherwise rebuild it
- If no patterns or cutoffs defined → use original data

Removed samples, with the reason and their total mutation count, are listed in
`data/input_cleaned/<TYPE>/preprocessing_removed_samples.csv`. If the
cutoffs are given per species, every species in the input folder must have
one; a missing one stops the run.

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

`data/input_cleaned/<TYPE>/preprocessing_stamp.json` records the patterns
and a SHA-256 hash of each input file. The cache is reused only if both still
match; otherwise it is deleted and rebuilt.

```
data/input_cleaned/DBS/ exists, stamp matches
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

SBS:
  exclude:
    - Xenon

# To exclude nothing for a type, omit it or use an empty list:
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

From an SBS run with the current config (paths shortened):

```
======================================================================
DATA DIRECTORY SELECTION
======================================================================
✓ Found preprocessing config: config/preprocessing.yaml
✓ Exclusion patterns for SBS: ['Xenon', 'Deoxynivalenol']
✓ Minimum mutations for SBS: {'mouse': 307, 'human': 307, 'celegans': 307, 'chicken': 307, 'rat': 307}

Running preprocessing...
  ...
  Processing: unfiltered_celegans_SBS96.txt  (cutoff: 307 mutations)
    Original: 245 samples
    Removed:  155 samples
    Kept:     90 samples
    Removed (first 3):
      - CX-5461_NO_UVA19 (below 307 mutations; 48.0 mutations)
      - CX-5461_NO_UVA1 (below 307 mutations; 301.0 mutations)
      - CX-5461_NO_UVA20 (below 307 mutations; 45.0 mutations)
      ... and 152 more
    Wrote normalized profiles: normalized_celegans_SBS96.tsv
  ...
PREPROCESSING COMPLETE
  Total samples: 1552
  Removed: 858
  Kept: 694
  Cleaned data saved to: data/input_cleaned/SBS
  Removed samples listed in: preprocessing_removed_samples.csv
======================================================================

✓ Using cleaned data: data/input_cleaned/SBS
```

On later runs with unchanged inputs and settings, the log shows
`✓ Using cached cleaned data` instead.

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

**Solution:** Add rules to `config/compound_grouping.yaml`:

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

| Step | Script | Key outputs in `results/<run name>/<TYPE>/` |
|------|--------|----------------------------------|
| 0 (optional) | *integrated preprocessing* | `data/input_cleaned/<TYPE>/` (cached, auto-generated) |
| 1 | `perform_clustering.py` | `main_clusters/`, `small_clusters/`, `singletons/`, `interactive_clustering/` |
| 2 | `generate_summary.py` (pie slices and legends largest first) | `reports/cluster_report.html`, `reports/cluster_summary.html`, `reports/main_clusters_with_dendrogram.html`, `reports/small_clusters_cluster_report.html`, `reports/small_clusters_summary.html` |
| 2b | `html_to_pdf.py` | `reports/cluster_summary.pdf`, `reports/main_clusters_with_dendrogram.pdf`, `reports/small_clusters_summary.pdf` (one page each) and a 300 DPI `.png` of each; needs Google Chrome |
| 3 | `generate_images_heatmap.py` | `interactive_heatmap/cluster_plots/*.png`, `interactive_heatmap/cosmic_plots/*.png` |
| 4 | `generate_interactive_heatmap.py` | `interactive_heatmap.html`, `cosmic_similar_clusters_<t>.tsv`, `denovo_clusters_<t>.tsv` |
| 5 | `generate_static_heatmap.py` | `cosine_similarity_heatmap.pdf` (3 versions) |
| 6 | `species_compound_matrix.py` | `species_compound_matrix_*.pdf` |
| 6b | `combined_matrix_figure.py` | `species_compound_matrix_combined.html`, `.pdf` (one page) and `.png` (300 DPI); PDF/PNG need Google Chrome |
| 7 | `run_info.py` | `run_parameters.txt`, `run_parameters.json` |
| 8 (SBS) | `run_decomposition.py` | `decomposition/` (COSMIC decomposition into eSS) |
| 9 | `collect_final_outputs.py` | `final_output/` (copies of all figure PDFs/PNGs, see [Final figures](#final-figures-final_output)) |

### One-page PDFs of the cluster summaries (Step 2b)

Every run prints the three grid reports (mutational profile + pie chart +
sample legend, 3 columns) to PDFs, plus a 300 DPI PNG of each, next to their
HTML, for use as manuscript figures and slides:

- `reports/cluster_summary.pdf` / `.png`
- `reports/main_clusters_with_dendrogram.pdf` / `.png`
- `reports/small_clusters_summary.pdf` / `.png`

Each PDF is **one page**: A4-landscape width (the reports' print layout) and
as tall as the content, so no rows are split onto a second page and nothing
has to be put back together by hand. `pipeline/html_to_pdf.py` prints them with headless Google
Chrome (the same rendering as printing from the browser), on a page tall
enough for the whole report, then trims the empty space under the last row.

Chrome is found at its usual macOS or Linux location, or set `CHROME_PATH`.
Without Chrome the step prints a warning and the run continues without the
PDFs.

To make the PDFs for an existing run folder, or after editing a report:

```bash
python pipeline/html_to_pdf.py --reports_dir results/min307/SBS/reports
python pipeline/html_to_pdf.py --reports_dir results/min307/SBS/reports --html cluster_report.html
```

The PNG is the whole page as one image at 300 DPI, made with `pdftoppm`
(Poppler). For another resolution:

```bash
pdftoppm -png -r 600 -singlefile results/min307/SBS/reports/cluster_summary.pdf results/min307/SBS/reports/cluster_summary_600dpi
```

### Combined species × compound figure (Step 6b)

Step 6 draws three species × compound matrices (`species_compound_matrix_human.pdf`,
`_mouse.pdf`, `_other_species.pdf`). Step 6b puts them together as one figure,
in the run folder:

- `species_compound_matrix_combined.html`
- `species_compound_matrix_combined.pdf` (one page, A4-landscape width)
- `species_compound_matrix_combined.png` (the same page at 300 DPI)

Panels **a** Human, **b** Mouse and **c** Other experimental models, each
titled with its number of models, exposures (columns) and profiles. The
panels are the Step 6 PDFs themselves, converted to vector SVG with
`pdftocairo` (Poppler), so the figure always matches Step 6. Each PDF's own
legend is cropped off, one shared legend is placed next to panel c, and the
three panels are scaled so their grid cells are the same size. The PDF is
printed with headless Google Chrome, as in Step 2b (skipped with a warning if
Chrome is not found).

To make it for an existing run folder:

```bash
python pipeline/combined_matrix_figure.py --run_dir results/min307/SBS
```

---

## Mutation type differences

| | SBS | DBS | ID |
|--|-----|-----|----|
| Contexts | 96 | 78 | 83 |
| Cluster prefix | `eSS` | `eDS` | `eIS` |
| Normalisation | pre-computed file | internal (col sum → 1) | internal (col sum → 1) |
| Artifact signatures | SBS27, 43, 45–60, 95 | *(none)* | ID9 |
| Preprocessing | min_mutations 307 + 10 mouse MEF samples excluded | hTumor exclusion | hTumor exclusion |
| Matrix generation | ✓ | ✓ | ✓ |

All type-specific constants are in `pipeline/utils/mutation_type.py`.

---

## Trinucleotide Opportunity Normalization (Optional)

Two methods correct for trinucleotide-composition differences between genomes
and between whole-genome and exome sequencing (manuscript revision item A4).
Neither is part of the default run; each is an option that is applied in
preprocessing, after the 307 filter:

```bash
bash run_pipeline.sh SBS 0.9 0.85 wes-to-wgs        # → results/min307_wes-to-wgs/SBS/
bash run_pipeline.sh SBS 0.9 0.85 own-opportunity   # → results/min307_own-opportunity/SBS/
```

- **`wes-to-wgs`**: human and mouse WES samples are rescaled onto their own
  genome's WGS opportunity, keeping each sample's total. Everything else is
  unchanged.
- **`own-opportunity`**: every human, mouse and rat sample (WGS and WES) is
  normalized by its own genome and technology opportunity. Counts are
  unchanged; only the profiles used for clustering change.

Both keep the same samples as the default run (the totals used for the 307
filter are not changed), use the opportunity tables in
`data/references/context_distributions/` (from SigProfilerMatrixGenerator
1.3.6), and leave chicken and C. elegans uncorrected (no tables). The code is
`pipeline/utils/opportunity_normalization.py`; the preprocessed data goes to
its own cache (`data/input_cleaned/SBS_<method>/`). The tests check both
methods against an independent calculation.

See [NORMALIZATION_APPROACH.md](NORMALIZATION_APPROACH.md) for the methods,
the comparison with the default run, and open items (the genome builds used
have not been confirmed against each source study).

---

## Adding a new mutation type

1. Add a `MutationTypeConfig` entry in `pipeline/utils/mutation_type.py`
2. Add input files to `data/input/<TYPE>/` following the naming convention in `data/README.md`
3. Add the COSMIC file path to the `COSMIC_PROFILE` `case` statement in `run_pipeline.sh`
4. (Optional) Add preprocessing patterns to `config/preprocessing.yaml`

No other scripts need changes. If the new type should be protected like SBS,
add an expected-membership file and tests under `tests/`.

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

Not for SBS: the SBS inputs are unfiltered, and the normalized SBS files are
only made by preprocessing, so SBS won't run without it.
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

**"requires Python 3.11":** the pipeline only runs on Python 3.11, the version
the pinned packages are verified on. Create the environment as in
[Installation](#installation).

**`AssertionError` during clustering:** one of the [safety checks](#safety-checks)
found inconsistent data, such as a sample present in the raw counts but not the
normalized profiles. The message names the samples involved; fix the input
files rather than bypassing the check.

**Tests fail after changing inputs or settings:** expected if the change was
intentional. For new inputs, `DataVersionTest` says the data isn't a
registered version: register it as in [Input data versions](#input-data-versions).
For changed settings, check the reported differences, then regenerate the
current version's expected results with `make_expected.py --force`.
