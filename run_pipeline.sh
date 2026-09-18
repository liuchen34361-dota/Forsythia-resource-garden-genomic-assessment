#!/usr/bin/env bash
set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_INPUT="${1:-${CODE_DIR}/config/analysis_config.json}"
PROJECT_ROOT_INPUT="${2:-.}"

CONFIG="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${CONFIG_INPUT}")"
PROJECT_ROOT="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${PROJECT_ROOT_INPUT}")"

if [[ ! -f "${CONFIG}" ]]; then
    echo "[ERROR] Configuration not found: ${CONFIG}" >&2
    exit 1
fi
if [[ ! -d "${PROJECT_ROOT}" ]]; then
    echo "[ERROR] Project root not found: ${PROJECT_ROOT}" >&2
    exit 1
fi

run_step() {
    local script="$1"
    echo "[RUN] ${script}"
    python3 "${CODE_DIR}/scripts/${script}" --config "${CONFIG}"
}

cd "${PROJECT_ROOT}"

run_step 01_wild_minor_allele_representation.py
run_step 02_AC_threshold_and_regional_gap.py
run_step 03_sample_size_standardized_retention.py
run_step 04_callability_and_carrier_quality.py
run_step 05_PIHAT_components_and_core.py
run_step 06_core_SNP_ALT_coverage.py
run_step 07_greedy_plan_selection.py
run_step 08_region_matched_random_control.py
run_step 09_selected_sample_QC.py
run_step 10_final_wild_portfolios.py
run_step validate_outputs.py

echo "[DONE] Cleaned genomic assessment workflow completed successfully."
