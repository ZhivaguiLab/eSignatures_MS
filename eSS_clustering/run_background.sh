#!/bin/bash
# run_background.sh
#
# Runs run_pipeline.sh in the background, so the terminal can be closed,
# with a timestamped log in results/logs/.
#
# Usage: the same arguments as run_pipeline.sh, e.g.
#   bash run_background.sh SBS 0.9 0.85
#   bash run_background.sh SBS 0.9 0.85 pooled
#   bash run_background.sh SBS 0.9 0.85 per-species pooled
#
# Uses the conda environment in ESIG_ENV (default: esignatures), built with
#   conda create -n esignatures python=3.11
#   conda run -n esignatures pip install -r requirements-lock.txt
#
# Prints the log path and the process ID. Follow it with
#   tail -f <log path>
# The log ends with "Pipeline complete!" when the run succeeds.

set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <mutation-type> <cosine-clust> <cosine-heatmap> [options...]"
    echo "  (same arguments as run_pipeline.sh)"
    exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${ESIG_ENV:-esignatures}"
LOG_DIR="${REPO_ROOT}/results/logs"
mkdir -p "${LOG_DIR}"

# Log named after the run's settings, e.g. SBS_0.9_0.85_pooled_20261113_1224.log
LABEL=$(IFS=_; echo "$*" | tr '=/ ' '-__')
LOG="${LOG_DIR}/${LABEL}_$(date +%Y%m%d_%H%M%S).log"

nohup conda run --no-capture-output -n "${ENV_NAME}" \
    bash "${REPO_ROOT}/run_pipeline.sh" "$@" > "${LOG}" 2>&1 &
PID=$!

echo "Started run_pipeline.sh $* (PID ${PID}, env ${ENV_NAME})"
echo "Log: ${LOG}"
echo "Follow it with:  tail -f \"${LOG}\""
