#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ ! -d analysis_output/step06_patient_plot_tables ]]; then
    echo "ERROR: analysis_output/step06_patient_plot_tables not found." >&2
    echo "Run ./run_pipeline.sh first." >&2
    exit 1
fi

if [[ ! -d analysis_output/step08_figure2_plot_data ]]; then
    echo "ERROR: analysis_output/step08_figure2_plot_data not found." >&2
    echo "Run ./run_pipeline.sh first." >&2
    exit 1
fi

python figures/plot_all_SNP_stacked_summary.py
python figures/plot_figure1_patient_mutation_landscape.py
python figures/plot_figure2_cross_patient_variable_sites.py

echo "FIGURE REPRODUCTION COMPLETE"
