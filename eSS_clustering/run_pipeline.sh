#!/bin/bash
# run_pipeline.sh
#
# Runs the full eSignatures clustering and reporting pipeline for one mutation type.
#
# Usage:
#   bash run_pipeline.sh <mutation-type> <cosine-clust> <cosine-heatmap> [options...]
#
# Default: minimum 307 SBSs per sample for every species, no manual split,
# no opportunity normalization (config/preprocessing.yaml).
#
# Options (any combination, in any order):
#   min=<N>          minimum <N> mutations per sample for every species,
#                    instead of the default 307
#   per-species      per-species minimums (SBS: min_mutations_per_species in
#                    config/preprocessing.yaml: mouse 235, human 295,
#                    celegans 451, chicken 300, rat 6354)
#   aai-split        SBS, for testing: add the manual AAI/DBP split (the cluster
#                    containing the Aristolochic acid I / Dibenzo[a,l]pyrene
#                    samples is re-split at cosine distance 0.095)
#   wes-to-wgs       SBS: opportunity normalization, WES samples moved onto the
#                    WGS basis (see NORMALIZATION_APPROACH.md)
#   own-opportunity  SBS: opportunity normalization, every sample normalized by
#                    its own genome/technology opportunity
#   pooled           consensus profile per cluster from summed raw counts
#                    (perform_clustering.py --averaging_method pooled) instead
#                    of the default equal_replicate (mean of each sample's
#                    normalized profile); writes to results/<run name>_pooled/
#
# Each run writes to results/<run name>/<MUTATION_TYPE>/, where the run name
# records the settings, e.g. results/min307/SBS/, results/min250/SBS/,
# results/min-per-species/SBS/ or results/min307_wes-to-wgs/SBS/, and a
# run_parameters.txt listing
# every setting, input file hash, software version and the result counts.
#
# Examples:
#   bash run_pipeline.sh SBS 0.9 0.85                    # default
#   bash run_pipeline.sh SBS 0.9 0.85 min=250            # 250 for every species
#   bash run_pipeline.sh SBS 0.9 0.85 per-species        # per-species cutoffs
#   bash run_pipeline.sh SBS 0.9 0.85 aai-split          # testing: AAI/DBP split
#   bash run_pipeline.sh SBS 0.9 0.85 wes-to-wgs
#   bash run_pipeline.sh SBS 0.9 0.85 own-opportunity
#   bash run_pipeline.sh DBS 0.9 0.85
#   bash run_pipeline.sh SBS 0.9 0.85 per-species wes-to-wgs
#   bash run_pipeline.sh SBS 0.9 0.85 pooled             # pooled consensus profiles
#
# Step 2b prints the cluster summary reports, and Step 6b the combined species
# × compound figure, to one-page PDFs and 300 DPI PNGs (needs Google Chrome;
# skipped with a warning otherwise). The run first checks for Poppler (needed
# throughout) and Chrome. For SBS, the last step decomposes the COSMIC
# signatures into the run's eSS with SigProfilerAssignment (Step 8).

set -euo pipefail

# ---------------------------------------------------------------------------
# Argument validation
# ---------------------------------------------------------------------------
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <mutation-type> <cosine-clust> <cosine-heatmap> [min=<N> | per-species] [aai-split] [wes-to-wgs | own-opportunity] [pooled]"
    echo "  mutation-type   : SBS | DBS | ID"
    echo "  min=<N>         : optional; minimum <N> mutations per sample for every species (default 307)"
    echo "  per-species     : optional; per-species minimums from config/preprocessing.yaml"
    echo "  aai-split       : optional, SBS testing only; add the manual AAI/DBP split"
    echo "  wes-to-wgs      : optional, SBS; WES -> WGS opportunity normalization"
    echo "  own-opportunity : optional, SBS; per-sample opportunity normalization"
    echo "  pooled          : optional; pooled consensus profiles instead of equal_replicate"
    echo ""
    echo "Example: $0 SBS 0.9 0.85"
    exit 1
fi

CUSTOM_THRESHOLDS="none"
NORMALIZATION="none"
MIN_MUTATIONS=""        # empty: the default from config/preprocessing.yaml
AVERAGING_METHOD="equal_replicate"
for opt in "${@:4}"; do
    case "$opt" in
        min=*|per-species)
            if [ -n "${MIN_MUTATIONS}" ]; then
                echo "Error: choose only one of min=<N> / per-species"
                exit 1
            fi
            if [ "$opt" == "per-species" ]; then
                MIN_MUTATIONS="per-species"
            else
                MIN_MUTATIONS="${opt#min=}"
                if [[ ! "${MIN_MUTATIONS}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
                    echo "Error: min=<N> needs a number (got '$opt')"
                    exit 1
                fi
            fi ;;
        aai-split) CUSTOM_THRESHOLDS="aai-split" ;;
        wes-to-wgs|own-opportunity)
            if [ "${NORMALIZATION}" != "none" ]; then
                echo "Error: choose only one of wes-to-wgs / own-opportunity"
                exit 1
            fi
            NORMALIZATION="$opt" ;;
        pooled) AVERAGING_METHOD="pooled" ;;
        *) echo "Error: unknown option '$opt' (expected min=<N>, per-species, aai-split, wes-to-wgs, own-opportunity or pooled)"
           exit 1 ;;
    esac
done
if [ "${NORMALIZATION}" = "own-opportunity" ] && [ "${AVERAGING_METHOD}" = "pooled" ]; then
    # pooled sums raw counts, but own-opportunity changes only the profiles,
    # not the counts, so the consensus would lose the normalization.
    echo "Error: own-opportunity cannot be combined with pooled (pooled sums raw counts,"
    echo "       which own-opportunity leaves unnormalized). Use wes-to-wgs pooled, or own-opportunity alone."
    exit 1
fi

MUTATION_TYPE=$(echo "$1" | tr '[:lower:]' '[:upper:]')   # uppercase
COSINE_THRES_CLUST="$2"
COSINE_THRES_HEATMAP="$3"

# Validate mutation type
if [[ ! "$MUTATION_TYPE" =~ ^(SBS|DBS|ID)$ ]]; then
    echo "Error: mutation-type must be SBS, DBS, or ID (got '$MUTATION_TYPE')"
    exit 1
fi

# ---------------------------------------------------------------------------
# Paths — all relative to repo root
# ---------------------------------------------------------------------------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Input data: data/input/<MUTATION_TYPE>/
DATA_DIR="${REPO_ROOT}/data/input/${MUTATION_TYPE}"

# COSMIC reference profiles: data/references/
REFERENCES_DIR="${REPO_ROOT}/data/references"

CONFIG_DIR="${REPO_ROOT}/config"
PIPELINE_DIR="${REPO_ROOT}/pipeline"
MIN_MUTATIONS_ARGS=()
if [ -n "${MIN_MUTATIONS}" ]; then
    MIN_MUTATIONS_ARGS=(--min_mutations "${MIN_MUTATIONS}")
fi
# Output folder named by the run's settings, e.g. results/min307/
RUN_NAME=$(python "${REPO_ROOT}/pipeline/run_info.py" name \
    --mutation_type "${MUTATION_TYPE}" \
    --preprocessing_config "${REPO_ROOT}/config/preprocessing.yaml" \
    --custom_thresholds "${CUSTOM_THRESHOLDS}" \
    --normalization "${NORMALIZATION}" \
    --averaging_method "${AVERAGING_METHOD}" \
    ${MIN_MUTATIONS_ARGS[@]+"${MIN_MUTATIONS_ARGS[@]}"})
RESULTS_ROOT="${REPO_ROOT}/results/${RUN_NAME}"

# Validate that input directory exists
if [ ! -d "${DATA_DIR}" ]; then
    echo "Error: Input data directory not found: ${DATA_DIR}"
    echo "Expected structure: data/input/${MUTATION_TYPE}/"
    exit 1
fi

case "$MUTATION_TYPE" in
  SBS) COSMIC_PROFILE="${REFERENCES_DIR}/COSMIC_v3.6_SBS_GRCh38.txt" ;;
  DBS) COSMIC_PROFILE="${REFERENCES_DIR}/COSMIC_v3.6_DBS_GRCh38.txt" ;;
  ID)  COSMIC_PROFILE="${REFERENCES_DIR}/COSMIC_v3.6_ID_GRCh37.txt"  ;;
esac

# Validate that COSMIC file exists
if [ ! -f "${COSMIC_PROFILE}" ]; then
    echo "Error: COSMIC reference file not found: ${COSMIC_PROFILE}"
    echo "Expected location: data/references/"
    exit 1
fi

ABBREVIATION_FILE="${CONFIG_DIR}/abv_table_clusters.txt"
MAPPING_FILE="${CONFIG_DIR}/sample_mapping.tsv"

PERFORM_CLUSTERING_PY="${PIPELINE_DIR}/perform_clustering.py"
GENERATE_SUMMARY_PY="${PIPELINE_DIR}/generate_summary.py"
GENERATE_HEATMAP_IMAGES_PY="${PIPELINE_DIR}/generate_images_heatmap.py"
GENERATE_HEATMAP_PY="${PIPELINE_DIR}/generate_interactive_heatmap.py"
GENERATE_STATIC_HEATMAP_PY="${PIPELINE_DIR}/generate_static_heatmap.py"
GENERATE_MATRIX_PY="${PIPELINE_DIR}/species_compound_matrix.py"

# ---------------------------------------------------------------------------
# Derived paths (mutation-type namespaced)
# ---------------------------------------------------------------------------
FILTER_DIR="${RESULTS_ROOT}/${MUTATION_TYPE}"
MAIN_CLUSTERS_DIR="${FILTER_DIR}/main_clusters"
INTERACTIVE_DIR="${FILTER_DIR}/interactive_clustering"
DENDROGRAM_SVG="${INTERACTIVE_DIR}/all_esignature_models_dendrogram.svg"
ESIG_PROFILES="${MAIN_CLUSTERS_DIR}/main_clusters_NORMALIZED_weighted_avg_profiles.tsv"
ESIG_PDF_DIR="${MAIN_CLUSTERS_DIR}/main_clusters_hierarchical_plots"
INTERACTIVE_HEATMAP_DIR="${FILTER_DIR}/interactive_heatmap"

case "$MUTATION_TYPE" in
  SBS) ESIG_PDF="${ESIG_PDF_DIR}/SBS_96_plots_main_clusters_hierarchical.pdf" ;;
  DBS) ESIG_PDF="${ESIG_PDF_DIR}/DBS_78_plots_main_clusters_hierarchical.pdf" ;;
  ID)  ESIG_PDF="${ESIG_PDF_DIR}/ID_83_plots_main_clusters_hierarchical.pdf"  ;;
esac

mkdir -p "${RESULTS_ROOT}"

echo ""
echo "========================================"
echo "eSignatures Pipeline — ${MUTATION_TYPE}"
echo "Input data:   ${DATA_DIR}"
echo "COSMIC ref:   ${COSMIC_PROFILE}"
echo "Typed output: ${FILTER_DIR}"
echo "Cosine (clust):   ${COSINE_THRES_CLUST}"
echo "Cosine (heatmap): ${COSINE_THRES_HEATMAP}"
echo "Run name:         ${RUN_NAME}"
echo "Min mutations:    ${MIN_MUTATIONS:-default (config/preprocessing.yaml)}"
echo "Custom thresholds: ${CUSTOM_THRESHOLDS}"
echo "Normalization:    ${NORMALIZATION}"
echo "Consensus:        ${AVERAGING_METHOD}"
echo "========================================"

# ---------------------------------------------------------------------------
# Tool check
# ---------------------------------------------------------------------------
# Stops here if Poppler is missing (pdf2image needs it from Step 1 on);
# warns if Google Chrome is missing (only the PDFs/PNGs of Steps 2b and 6b
# need it, and those steps are then skipped).
echo ""
echo "Checking tools..."
python "${PIPELINE_DIR}/html_to_pdf.py" --check

# ---------------------------------------------------------------------------
# Step 1: Clustering
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Step 1: Clustering"
echo "========================================"

CLUSTERING_ARGS=(
    --output_dir "${RESULTS_ROOT}"
    --mutation_type "${MUTATION_TYPE}"
    --data_dir "${REPO_ROOT}/data/input"
    --cosine_similarity "${COSINE_THRES_CLUST}"
    --use_original_names
    --mapping_file "${MAPPING_FILE}"
)
# The AAI/DBP split is passed by preset name; its thresholds are defined
# once, in custom_threshold_presets in pipeline/utils/mutation_type.py.
if [ "${CUSTOM_THRESHOLDS}" != "none" ]; then
    CLUSTERING_ARGS+=(--custom_thresholds "${CUSTOM_THRESHOLDS}")
fi
CLUSTERING_ARGS+=(--normalization "${NORMALIZATION}")
CLUSTERING_ARGS+=(--averaging_method "${AVERAGING_METHOD}")
if [ -n "${MIN_MUTATIONS}" ]; then
    CLUSTERING_ARGS+=(--min_mutations "${MIN_MUTATIONS}")
fi

python "${PERFORM_CLUSTERING_PY}" "${CLUSTERING_ARGS[@]}"

# ---------------------------------------------------------------------------
# Step 2: Clustering summary reports
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Step 2: Clustering Summary Reports"
echo "========================================"
python "${GENERATE_SUMMARY_PY}" \
    --input_df_path "${MAIN_CLUSTERS_DIR}/main_clusters_summary.tsv" \
    --plots_dir "${ESIG_PDF_DIR}" \
    --dendrogram_svg "${DENDROGRAM_SVG}" \
    --main_output_dir "${FILTER_DIR}" \
    --main_working_dir "${FILTER_DIR}" \
    --create_dendrogram_summaries \
    --mutation_type "${MUTATION_TYPE}" \
    --abbreviation_file "${ABBREVIATION_FILE}"

# ---------------------------------------------------------------------------
# Step 2b: One-page PDFs of the cluster summaries
# ---------------------------------------------------------------------------
# cluster_summary.pdf, main_clusters_with_dendrogram.pdf and
# small_clusters_summary.pdf next to their HTML in reports/, each on a single
# page (A4-landscape width, as tall as the content). Printed with headless
# Google Chrome; skipped with a warning if Chrome is not installed.
echo ""
echo "========================================"
echo "Step 2b: Summary PDFs"
echo "========================================"
python "${PIPELINE_DIR}/html_to_pdf.py" \
    --reports_dir "${FILTER_DIR}/reports"

# ---------------------------------------------------------------------------
# Step 3: Heatmap images (PDF → PNG)
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Step 3: Heatmap Images"
echo "========================================"
python "${GENERATE_HEATMAP_IMAGES_PY}" \
    --esignature_pdf "${ESIG_PDF}" \
    --cosmic_txt "${COSMIC_PROFILE}" \
    --output_dir "${INTERACTIVE_HEATMAP_DIR}" \
    --mutation_type "${MUTATION_TYPE}"

# ---------------------------------------------------------------------------
# Step 4: Interactive cosine similarity heatmap
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Step 4: Interactive Heatmap"
echo "========================================"
python "${GENERATE_HEATMAP_PY}" \
    --esignature_profiles "${ESIG_PROFILES}" \
    --cosmic_profiles "${COSMIC_PROFILE}" \
    --image_directory "${INTERACTIVE_HEATMAP_DIR}" \
    --output_file "${FILTER_DIR}/interactive_heatmap.html" \
    --threshold "${COSINE_THRES_HEATMAP}" \
    --output_high_similarity "${FILTER_DIR}/cosmic_similar_clusters_${COSINE_THRES_HEATMAP}.tsv" \
    --output_low_similarity "${FILTER_DIR}/denovo_clusters_${COSINE_THRES_HEATMAP}.tsv"

# ---------------------------------------------------------------------------
# Step 5: Static cosine similarity heatmap
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Step 5: Static Heatmap"
echo "========================================"
python "${GENERATE_STATIC_HEATMAP_PY}" \
    --esignature_profiles "${ESIG_PROFILES}" \
    --cosmic_profiles "${COSMIC_PROFILE}" \
    --output_file "${FILTER_DIR}/cosine_similarity_heatmap.pdf" \
    --mutation_type "${MUTATION_TYPE}"

# ---------------------------------------------------------------------------
# Step 6: Species-model vs compound matrix (all mutation types)
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Step 6: Species-Compound Matrix"
echo "========================================"
python "${GENERATE_MATRIX_PY}" \
    --main_output_dir "${FILTER_DIR}" \
    --output_path "${FILTER_DIR}/species_compound_matrix.pdf" \
    --mutation_type "${MUTATION_TYPE}" \
    --abbreviation_file "${ABBREVIATION_FILE}" \
    --compound_grouping "${CONFIG_DIR}/compound_grouping.yaml"

# ---------------------------------------------------------------------------
# Step 6b: Combined species × compound figure
# ---------------------------------------------------------------------------
# The three Step 6 matrix PDFs as one figure (panels a Human, b Mouse,
# c Other, one shared legend, equal cell size):
# species_compound_matrix_combined.html and a one-page .pdf (needs Google
# Chrome for the PDF; skipped with a warning otherwise).
echo ""
echo "========================================"
echo "Step 6b: Combined Species-Compound Figure"
echo "========================================"
python "${PIPELINE_DIR}/combined_matrix_figure.py" \
    --run_dir "${FILTER_DIR}" \
    --abbreviation_file "${ABBREVIATION_FILE}" \
    --compound_grouping "${CONFIG_DIR}/compound_grouping.yaml"

# ---------------------------------------------------------------------------
# Step 7: Run parameter summary
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Step 7: Run Parameters"
echo "========================================"
python "${REPO_ROOT}/pipeline/run_info.py" summary \
    --output_dir "${FILTER_DIR}" \
    --run_name "${RUN_NAME}" \
    --cosmic_profiles "${COSMIC_PROFILE}" \
    --cosmic_threshold "${COSINE_THRES_HEATMAP}"

# ---------------------------------------------------------------------------
# Step 8: COSMIC decomposition (SBS only)
# ---------------------------------------------------------------------------
# Decompose each COSMIC v3.6 SBS signature into this run's eSS with
# SigProfilerAssignment 1.1.4 (decompose_fit). A COSMIC signature whose best
# eSS reconstruction has cosine similarity below 0.8 (SigProfilerAssignment's
# default novelty threshold) stays as itself. Output, with its own
# decomposition_parameters.txt, goes to <run folder>/<TYPE>/decomposition/.
# See pipeline/run_decomposition.py for the settings.
if [ "${MUTATION_TYPE}" == "SBS" ]; then
    echo ""
    echo "========================================"
    echo "Step 8: COSMIC Decomposition"
    echo "========================================"
    python "${PIPELINE_DIR}/run_decomposition.py" \
        --run_dir "${FILTER_DIR}" \
        --cosmic "${COSMIC_PROFILE}"
fi

# ---------------------------------------------------------------------------
# Step 9: Collect the final figures
# ---------------------------------------------------------------------------
# Copies of every figure PDF/PNG (cluster summaries, species x compound
# matrices, COSMIC heatmaps, dendrogram, signature profiles, decomposition
# plots) and run_parameters.txt into <run folder>/<TYPE>/final_output/,
# with a MANIFEST.txt of where each came from. Replaced on every run.
echo ""
echo "========================================"
echo "Step 9: Final Outputs"
echo "========================================"
python "${PIPELINE_DIR}/collect_final_outputs.py" \
    --run_dir "${FILTER_DIR}" \
    --mutation_type "${MUTATION_TYPE}"

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Pipeline complete! (${MUTATION_TYPE})"
echo "Results in: ${FILTER_DIR}"
echo "Figures in: ${FILTER_DIR}/final_output"
echo "========================================"