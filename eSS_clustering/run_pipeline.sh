#!/bin/bash
# run_pipeline.sh
#
# Runs the full eSignatures clustering and reporting pipeline for one mutation type.
# Outputs are written to results/<MUTATION_TYPE>/ so separate runs are isolated.
#
# Usage:
#   bash run_pipeline.sh <mutation-type> <cosine-clust> <cosine-heatmap>
#
# Examples:
#   bash run_pipeline.sh SBS 0.9 0.85
#   bash run_pipeline.sh DBS 0.9 0.85
#   bash run_pipeline.sh ID  0.9 0.85
#
# The canonical clustering threshold is 0.9 — this is also the default of
# perform_clustering.py, so running that script directly gives the same
# clusters. See "Reproducing the published clusters" in README.md.

set -euo pipefail

# ---------------------------------------------------------------------------
# Argument validation
# ---------------------------------------------------------------------------
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <mutation-type> <cosine-clust> <cosine-heatmap>"
    echo "  mutation-type : SBS | DBS | ID"
    echo ""
    echo "Example: $0 SBS 0.9 0.85"
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
RESULTS_ROOT="${REPO_ROOT}/results"

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
echo "========================================"

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
    --data_dir "${DATA_DIR}"
    --cosine_similarity "${COSINE_THRES_CLUST}"
    --use_original_names
    --mapping_file "${MAPPING_FILE}"
)

# Per-cluster custom thresholds (SBS: Aristolochic_acid_I and
# Dibenzo[a,l]pyrene at 0.095) are not passed here: perform_clustering.py
# applies them by default from pipeline/utils/mutation_type.py, so the
# script and this wrapper cannot drift apart.

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
    --compound_grouping "config/compound_grouping.yaml" 

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Pipeline complete! (${MUTATION_TYPE})"
echo "Results in: ${FILTER_DIR}"
echo "========================================"