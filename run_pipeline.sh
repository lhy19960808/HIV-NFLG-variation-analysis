#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# HIV-1 NFLG within-host variation analysis pipeline
# ============================================================

SCRIPT_DIR="$(
    cd "$(dirname "${BASH_SOURCE[0]}")"
    pwd
)"

PROJECT_ROOT="$SCRIPT_DIR"

cd "$PROJECT_ROOT"

OUT="${1:-$SCRIPT_DIR/analysis_output}"

THREADS="${THREADS:--1}"

INPUT_FASTA="${INPUT_FASTA:-$PROJECT_ROOT/data/All_185_clones_no_vector_no_primer_plus_consensus_renamed.fasta}"

HXB2_SOURCE="${HXB2_SOURCE:-$PROJECT_ROOT/data/all_patients_NFLG_with_HXB2_original.fasta}"

CODE="$SCRIPT_DIR/code_clean"
CONFIG="$SCRIPT_DIR/config"

mkdir -p "$OUT/logs"


# ============================================================
# Basic input checks
# ============================================================

for f in \
"$INPUT_FASTA" \
"$HXB2_SOURCE" \
"$CONFIG/patient_clone_counts.tsv" \
"$CONFIG/HXB2_gene_coordinates.tsv" \
"$CONFIG/special_regions.tsv"
do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: required file not found:"
        echo "  $f"
        exit 1
    fi
done


if ! command -v python >/dev/null 2>&1; then
    echo "ERROR: python not found in PATH"
    exit 1
fi

if ! command -v mafft >/dev/null 2>&1; then
    echo "ERROR: mafft not found in PATH"
    exit 1
fi


run_step() {
    local label="$1"
    local log="$2"
    shift 2

    echo
    echo "============================================================"
    echo "$label"
    echo "============================================================"

    if "$@" > "$log" 2>&1
    then
        echo "$label: PASS"
    else
        echo "$label: FAIL"
        echo
        echo "Last 50 log lines:"
        tail -50 "$log"
        exit 1
    fi
}


# ============================================================
# Step 01
# ============================================================

run_step \
"STEP01 patient-specific alignments" \
"$OUT/logs/step01.log" \
python "$CODE/01_prepare_patient_alignments.py" \
    --input-fasta "$INPUT_FASTA" \
    --hxb2-source "$HXB2_SOURCE" \
    --output-dir "$OUT/step01_patient_alignments" \
    --threads "$THREADS" \
    --allow-missing-consensus CN2023AH1-17


# ============================================================
# Step 02
# ============================================================

run_step \
"STEP02 patient-specific analysis consensus" \
"$OUT/logs/step02.log" \
python "$CODE/02_build_patient_consensus.py" \
    --alignment-dir "$OUT/step01_patient_alignments/aligned" \
    --output-dir "$OUT/step02_patient_consensus"


# ============================================================
# Step 03
# ============================================================

run_step \
"STEP03 within-patient mutation calling" \
"$OUT/logs/step03.log" \
python "$CODE/03_call_within_patient_mutations.py" \
    --alignment-dir "$OUT/step01_patient_alignments/aligned" \
    --consensus-dir "$OUT/step02_patient_consensus" \
    --output-dir "$OUT/step03_mutation_calls"


# ============================================================
# Step 04A
# ============================================================

run_step \
"STEP04A mapping-consensus bridge" \
"$OUT/logs/step04a.log" \
python "$CODE/04a_build_mapping_consensus_bridge.py" \
    --alignment-dir "$OUT/step01_patient_alignments/aligned" \
    --output-dir "$OUT/step04a_mapping_bridge"


# ============================================================
# Step 04B
# ============================================================

run_step \
"STEP04B HXB2 coordinate mapping" \
"$OUT/logs/step04b.log" \
python "$CODE/04b_map_patient_consensus_to_HXB2.py" \
    --bridge-dir "$OUT/step04a_mapping_bridge" \
    --hxb2-source "$HXB2_SOURCE" \
    --output-dir "$OUT/step04b_hxb2_mapping"


# ============================================================
# Step 05
# ============================================================

run_step \
"STEP05 mutation projection to HXB2" \
"$OUT/logs/step05.log" \
python "$CODE/05_project_mutations_to_HXB2.py" \
    --mutation-dir "$OUT/step03_mutation_calls" \
    --bridge-dir "$OUT/step04a_mapping_bridge" \
    --hxb2-map-dir "$OUT/step04b_hxb2_mapping" \
    --output-dir "$OUT/step05_final_mutation_projection"


# ============================================================
# Step 06
# ============================================================

run_step \
"STEP06 patient plotting tables" \
"$OUT/logs/step06.log" \
python "$CODE/06_prepare_patient_plot_tables.py" \
    --projection-dir "$OUT/step05_final_mutation_projection" \
    --output-dir "$OUT/step06_patient_plot_tables"


# ============================================================
# Step 07
# ============================================================

run_step \
"STEP07 cross-patient dataset" \
"$OUT/logs/step07.log" \
python "$CODE/07_prepare_cross_patient_dataset.py" \
    --projection-dir "$OUT/step05_final_mutation_projection" \
    --plot-table-dir "$OUT/step06_patient_plot_tables" \
    --clone-counts "$CONFIG/patient_clone_counts.tsv" \
    --gene-coordinates "$CONFIG/HXB2_gene_coordinates.tsv" \
    --special-regions "$CONFIG/special_regions.tsv" \
    --output-dir "$OUT/step07_cross_patient_dataset"


# ============================================================
# Step 08
# ============================================================

run_step \
"STEP08 Figure 2 plotting data" \
"$OUT/logs/step08.log" \
python "$CODE/08_prepare_figure2_plot_data.py" \
    --cross-patient-dir "$OUT/step07_cross_patient_dataset" \
    --output-dir "$OUT/step08_figure2_plot_data"


# ============================================================
# Release-level validation
# ============================================================

if [[ -f "$SCRIPT_DIR/tests/check_release_output.py" ]]; then
    echo
    echo "============================================================"
    echo "FINAL RELEASE QC"
    echo "============================================================"

    python \
        "$SCRIPT_DIR/tests/check_release_output.py" \
        "$OUT"
fi


echo
echo "============================================================"
echo "PIPELINE COMPLETE"
echo "============================================================"
echo "Output directory:"
echo "$OUT"
