#!/usr/bin/env python3

from pathlib import Path
import csv
import sys


if len(sys.argv) != 2:
    raise SystemExit(
        "Usage: check_release_output.py <analysis_output>"
    )

ROOT = Path(sys.argv[1])

failures = []


def read_qc(path):
    if not path.exists():
        raise FileNotFoundError(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as f:
        return list(
            csv.DictReader(
                f,
                delimiter="\t",
            )
        )


def check_table(label, path, expected):
    rows = read_qc(path)

    by_check = {
        row["Check"]: row
        for row in rows
    }

    bad_status = [
        row["Check"]
        for row in rows
        if row.get("Status") != "PASS"
    ]

    if bad_status:
        print(
            f"{label:<12} FAIL "
            f"(non-PASS rows: {bad_status})"
        )
        failures.append(label)
        return

    mismatches = []

    for key, expected_value in expected.items():

        if key not in by_check:
            mismatches.append(
                f"{key}=MISSING"
            )
            continue

        observed = by_check[key]["Observed"]

        if observed != str(expected_value):
            mismatches.append(
                f"{key}: {observed} != {expected_value}"
            )

    if mismatches:
        print(f"{label:<12} FAIL")

        for x in mismatches:
            print(f"    {x}")

        failures.append(label)
        return

    print(f"{label:<12} PASS")


print("=" * 72)
print("HIV NFLG RELEASE OUTPUT VALIDATION")
print("=" * 72)


check_table(
    "Step05",
    ROOT /
    "step05_final_mutation_projection/"
    "03_projection_QC.tsv",
    {
        "Directional_mutation_calls": 3287,
        "SNP_call_n": 2477,
        "DEL_base_call_n": 592,
        "INS_base_call_n": 218,
        "DEL_event_n": 44,
        "INS_event_n": 31,
        "Total_indel_event_n": 75,
        "NonPASS_call_n": 0,
        "Exact_anchor_INS_call_n": 218,
        "Interval_anchor_INS_call_n": 0,
        "Patient_range_n": 11,
        "READY_FOR_DOWNSTREAM_ANALYSIS": "YES",
    },
)


check_table(
    "Step06",
    ROOT /
    "step06_patient_plot_tables/"
    "04_plot_data_QC.tsv",
    {
        "SNP_directional_call_n": 2477,
        "DEL_event_n": 44,
        "INS_event_n": 31,
        "Patient_plot_total_INDEL_event_count": 75,
        "Panel_patient_n": 11,
        "READY_FOR_FIGURE": "YES",
    },
)


check_table(
    "Step07",
    ROOT /
    "step07_cross_patient_dataset/"
    "11_final_data_QC.tsv",
    {
        "SNP_call_n": 2477,
        "DEL_event_n": 44,
        "INS_event_n": 31,
        "Patient_n": 11,
        "P18_missing_region_length_matches_MSA_gap": 2417,
        "P18_missing_region_HXB2_interval": "3487-5903",
        "Figure2_primary_metric": "Patient_n",
        "READY_FOR_CROSS_PATIENT_ANALYSIS": "YES",
    },
)


check_table(
    "Step08",
    ROOT /
    "step08_figure2_plot_data/"
    "05_figure2_plot_QC.tsv",
    {
        "SNP_locus_or_anchor_n": 973,
        "Recurrent_SNP_locus_n": 52,
        "Same_substitution_recurrent_locus_n": 21,
        "Mixed_substitution_recurrent_locus_n": 31,
        "DEL_locus_n": 32,
        "INS_anchor_n": 13,
        "READY_FOR_FIGURE2_PLOTTING": "YES",
    },
)


print("=" * 72)

if failures:
    print(
        "RELEASE OUTPUT VALIDATION: FAIL"
    )
    print("=" * 72)
    sys.exit(1)

print(
    "RELEASE OUTPUT VALIDATION: PASS"
)
print("=" * 72)
