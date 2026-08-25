#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FINAL plot-table preparation
================================

Purpose
-------
Prepare the final plotting tables from the validated mutation coordinates.

IMPORTANT
---------
This script does NOT change:
- mutation calls;
- HXB2 coordinates;
- validated HXB2 coordinate projection;
- indel event definitions.

It fixes only the plotting aggregation rule.

Locked plotting counts
----------------------
SNP:
    Count_for_plot = number of UNIQUE variant clones at that exact
    HXB2 / extended-HXB2 coordinate and substitution type.

DEL / INS:
    Count_for_plot = number of contiguous clone-specific INDEL EVENT
    OCCURRENCES at that locus.

This distinction is important because the same clone can, in principle,
contain >1 distinct indel event that projects to the same HXB2 anchor.
Using only unique clone count would collapse such events.

Cross-patient figure
--------------------
Primary metric:
    Patient_n = number of UNIQUE patients carrying a variant at the locus.

Auxiliary:
    Total_variant_clone_n
    Total_event_n

Inputs
------
Clean Step05 final mutation calls
Clean Step05 final INDEL events
Clean Step05 final patient HXB2 ranges
Clean Step05 no-unique-consensus regions

Outputs
-------
Final plotting-table outputs:
  01_patient_plot_ready.tsv
  02_all_patient_plot_ready.tsv
  03_patient_mutation_summary.tsv
  04_plot_data_QC.tsv
  05_patient_panel_metadata.tsv
  06_no_unique_consensus_regions.tsv

Final figure construction should proceed only if:
    READY_FOR_FIGURE = YES
"""

from pathlib import Path
from collections import defaultdict, Counter
import argparse
import csv
import re


parser = argparse.ArgumentParser(
    description=(
        "Prepare patient-level and combined plotting tables "
        "from the validated clean mutation projection."
    )
)

parser.add_argument(
    "--projection-dir",
    required=True,
    type=Path,
    help="Clean Step 05 final mutation-projection directory."
)

parser.add_argument(
    "--output-dir",
    required=True,
    type=Path,
    help="Output directory for final plotting tables."
)

args = parser.parse_args()

PROJECTION = args.projection_dir

CALLS = (
    PROJECTION /
    "01_final_clone_mutation_calls.tsv"
)

EVENTS = (
    PROJECTION /
    "02_final_indel_events.tsv"
)

RANGES = (
    PROJECTION /
    "06_final_patient_HXB2_ranges.tsv"
)

UNRESOLVED = (
    PROJECTION /
    "07_no_unique_consensus_regions.tsv"
)

OUT = args.output_dir
OUT.mkdir(parents=True, exist_ok=True)

PATIENT_PLOT = OUT / "01_patient_plot_ready.tsv"
ALL_PLOT = OUT / "02_all_patient_plot_ready.tsv"
PATIENT_SUMMARY = OUT / "03_patient_mutation_summary.tsv"
QC_OUT = OUT / "04_plot_data_QC.tsv"
PANEL_META = OUT / "05_patient_panel_metadata.tsv"
UNRESOLVED_OUT = OUT / "06_no_unique_consensus_regions.tsv"


def read_tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, fields, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def require(path):
    if not path.exists():
        raise FileNotFoundError(path)


def patient_sort_key(patient):
    m = re.fullmatch(r"CN(\d{4})AH(\d+)-(\d+)", str(patient))
    if m:
        return tuple(map(int, m.groups()))
    return (9999, 9999, str(patient))


def clone_sort_key(name):
    m = re.search(r"\.clone(\d+)", str(name), re.I)
    if m:
        return int(m.group(1))
    return 10**9


for path in [CALLS, EVENTS, RANGES, UNRESOLVED]:
    require(path)

calls = read_tsv(CALLS)
events = read_tsv(EVENTS)
ranges = read_tsv(RANGES)
unresolved = read_tsv(UNRESOLVED)

# ============================================================
# Sanity checks
# ============================================================

for r in calls:
    if not str(r["Projection_status"]).startswith("PASS_"):
        raise RuntimeError(
            f"Non-PASS mutation call remains: "
            f"{r['Patient']} {r['Clone']} "
            f"{r['Original_patient_MSA_column']} "
            f"{r['Projection_status']}"
        )

for r in events:
    if r["Projection_status"] != "PASS":
        raise RuntimeError(
            f"Non-PASS indel event remains: {r}"
        )

# ============================================================
# Patient plot table
# ============================================================

patient_rows = []

# SNP aggregation: unique clones at exact coord + exact substitution.
snp_groups = defaultdict(lambda: {
    "clones": set(),
    "genes": set(),
    "classes": set(),
})

for r in calls:
    if r["Mutation_type"] != "SNP":
        continue

    key = (
        r["Patient"],
        r["Final_HXB2_coordinate"],
        r["Final_HXB2_anchor"],
        r["Final_plot_coordinate"],
        r["SNP_type"],
    )

    g = snp_groups[key]
    g["clones"].add(r["Clone"])

    for gene in str(r["Final_genes"]).split(","):
        if gene:
            g["genes"].add(gene)

    if r.get("Final_coordinate_class"):
        g["classes"].add(r["Final_coordinate_class"])


for key, g in snp_groups.items():
    patient, coord, anchor, plot_x, snp_type = key

    patient_rows.append({
        "Patient": patient,
        "Record_type": "SNP",
        "Mutation_type": "SNP",
        "SNP_type": snp_type,
        "Exact_HXB2_coordinate": coord,
        "HXB2_anchor": anchor,
        "Plot_coordinate": plot_x,
        "Coordinate_class":
            ",".join(sorted(g["classes"])),
        "Count_for_plot": len(g["clones"]),
        "Variant_clone_n": len(g["clones"]),
        "Event_n": "",
        "Genes": ",".join(sorted(g["genes"])),
        "Display_symbol": "stacked_bar",
        "Projection_status": "PASS",
    })


# INDEL aggregation: EVENT OCCURRENCES, not unique clones.
event_groups = defaultdict(lambda: {
    "events": [],
    "clones": set(),
    "genes": set(),
    "lengths": [],
})

for r in events:
    key = (
        r["Patient"],
        r["Event_type"],
        r["Event_HXB2_locus"],
        r["Plot_coordinate"],
    )

    g = event_groups[key]
    g["events"].append(r["Event_ID"])
    g["clones"].add(r["Clone"])
    g["lengths"].append(str(r["Event_length_bp"]))

    for gene in str(r["Genes"]).split(","):
        if gene:
            g["genes"].add(gene)


for key, g in event_groups.items():
    patient, mut, locus, plot_x = key

    patient_rows.append({
        "Patient": patient,
        "Record_type": "INDEL_EVENT",
        "Mutation_type": mut,
        "SNP_type": "",
        "Exact_HXB2_coordinate": locus,
        "HXB2_anchor": locus,
        "Plot_coordinate": plot_x,
        "Coordinate_class": "HXB2_event_locus",

        # Locked rule:
        "Count_for_plot": len(g["events"]),
        "Variant_clone_n": len(g["clones"]),
        "Event_n": len(g["events"]),

        "Genes": ",".join(sorted(g["genes"])),
        "Display_symbol":
            "diamond" if mut == "DEL" else "inverted_triangle",
        "Event_length_bp_values":
            ",".join(sorted(g["lengths"], key=lambda x: int(float(x)))),
        "Event_IDs":
            ",".join(g["events"]),
        "Projection_status": "PASS",
    })


patient_rows.sort(
    key=lambda r: (
        patient_sort_key(r["Patient"]),
        float(r["Plot_coordinate"])
        if str(r["Plot_coordinate"]).strip()
        else 1e20,
        r["Mutation_type"],
        r["SNP_type"],
    )
)

patient_fields = [
    "Patient",
    "Record_type",
    "Mutation_type",
    "SNP_type",
    "Exact_HXB2_coordinate",
    "HXB2_anchor",
    "Plot_coordinate",
    "Coordinate_class",
    "Count_for_plot",
    "Variant_clone_n",
    "Event_n",
    "Genes",
    "Display_symbol",
    "Event_length_bp_values",
    "Event_IDs",
    "Projection_status",
]

write_tsv(PATIENT_PLOT, patient_fields, patient_rows)


# ============================================================
# Patient summary
# ============================================================

summary = defaultdict(lambda: {
    "SNP": 0,
    "DEL": 0,
    "INS": 0,
})

for r in patient_rows:
    summary[r["Patient"]][r["Mutation_type"]] += int(r["Count_for_plot"])

panel_patients = sorted(
    {r["Patient"] for r in ranges},
    key=patient_sort_key,
)

summary_rows = []

for patient in panel_patients:
    summary_rows.append({
        "Patient": patient,
        "SNP_call_n": summary[patient]["SNP"],
        "DEL_event_n": summary[patient]["DEL"],
        "INS_event_n": summary[patient]["INS"],
        "Total_indel_event_n":
            summary[patient]["DEL"] + summary[patient]["INS"],
    })

write_tsv(
    PATIENT_SUMMARY,
    [
        "Patient",
        "SNP_call_n",
        "DEL_event_n",
        "INS_event_n",
        "Total_indel_event_n",
    ],
    summary_rows,
)


# ============================================================
# All-patient table
# ============================================================

# For SNP:
#   common locus = HXB2_anchor
#   preserve substitution type
#
# For INDEL:
#   common locus = event HXB2 locus/anchor
#
# Patient_n is the PRIMARY across-patient metric.

all_groups = defaultdict(lambda: {
    "patients": set(),
    "clone_sum": 0,
    "event_sum": 0,
    "exact_coords": set(),
    "genes": set(),
})

for r in patient_rows:
    if r["Mutation_type"] == "SNP":
        locus = r["HXB2_anchor"] or r["Exact_HXB2_coordinate"]
        snp_type = r["SNP_type"]
    else:
        locus = r["HXB2_anchor"]
        snp_type = ""

    key = (
        r["Mutation_type"],
        snp_type,
        locus,
        r["Plot_coordinate"],
        r["Display_symbol"],
    )

    g = all_groups[key]
    g["patients"].add(r["Patient"])
    g["clone_sum"] += int(r["Variant_clone_n"])

    if str(r["Event_n"]).strip():
        g["event_sum"] += int(r["Event_n"])

    g["exact_coords"].add(r["Exact_HXB2_coordinate"])

    for gene in str(r["Genes"]).split(","):
        if gene:
            g["genes"].add(gene)


all_rows = []

for key, g in all_groups.items():
    mut, snp_type, locus, plot_x, symbol = key

    all_rows.append({
        "Mutation_type": mut,
        "SNP_type": snp_type,
        "HXB2_locus_or_anchor": locus,
        "Plot_coordinate": plot_x,

        "Patient_n": len(g["patients"]),
        "Patients": ",".join(
            sorted(g["patients"], key=patient_sort_key)
        ),

        "Total_variant_clone_n": g["clone_sum"],
        "Total_event_n":
            g["event_sum"] if mut in {"DEL", "INS"} else "",

        "Exact_coordinate_n": len(g["exact_coords"]),
        "Exact_coordinates": ",".join(sorted(g["exact_coords"])),

        "Genes": ",".join(sorted(g["genes"])),
        "Display_symbol": symbol,
        "Projection_status": "PASS",
    })


all_rows.sort(
    key=lambda r: (
        float(r["Plot_coordinate"])
        if str(r["Plot_coordinate"]).strip()
        else 1e20,
        r["Mutation_type"],
        r["SNP_type"],
    )
)

write_tsv(
    ALL_PLOT,
    [
        "Mutation_type",
        "SNP_type",
        "HXB2_locus_or_anchor",
        "Plot_coordinate",
        "Patient_n",
        "Patients",
        "Total_variant_clone_n",
        "Total_event_n",
        "Exact_coordinate_n",
        "Exact_coordinates",
        "Genes",
        "Display_symbol",
        "Projection_status",
    ],
    all_rows,
)


# ============================================================
# Panel metadata: guarantees all 11 patients remain in Figure 1
# ============================================================

unresolved_by_patient = Counter(
    r["Patient"] for r in unresolved
)

panel_rows = []

for r in sorted(ranges, key=lambda x: patient_sort_key(x["Patient"])):
    patient = r["Patient"]

    panel_rows.append({
        "Patient": patient,
        "HXB2_start": r["Final_HXB2_start"],
        "HXB2_end": r["Final_HXB2_end"],
        "HXB2_span": r["Final_HXB2_span"],
        "Has_directional_mutation_records":
            "YES" if patient in summary else "NO",
        "No_unique_region_n": unresolved_by_patient.get(patient, 0),
    })

write_tsv(
    PANEL_META,
    [
        "Patient",
        "HXB2_start",
        "HXB2_end",
        "HXB2_span",
        "Has_directional_mutation_records",
        "No_unique_region_n",
    ],
    panel_rows,
)

# Copy unresolved table without modifying coordinate interpretation.
if unresolved:
    unresolved_fields = list(unresolved[0].keys())
    write_tsv(UNRESOLVED_OUT, unresolved_fields, unresolved)
else:
    write_tsv(
        UNRESOLVED_OUT,
        [
            "Patient",
            "No_unique_class",
            "Start_original_MSA_column",
            "End_original_MSA_column",
            "Length_alignment_columns",
            "Left_anchor_HXB2_coordinate",
            "Right_anchor_HXB2_coordinate",
            "Display_interval",
        ],
        [],
    )


# ============================================================
# QC
# ============================================================

call_snp_n = sum(r["Mutation_type"] == "SNP" for r in calls)
event_del_n = sum(r["Event_type"] == "DEL" for r in events)
event_ins_n = sum(r["Event_type"] == "INS" for r in events)

plot_snp_n = sum(
    int(r["Count_for_plot"])
    for r in patient_rows
    if r["Mutation_type"] == "SNP"
)

plot_del_n = sum(
    int(r["Count_for_plot"])
    for r in patient_rows
    if r["Mutation_type"] == "DEL"
)

plot_ins_n = sum(
    int(r["Count_for_plot"])
    for r in patient_rows
    if r["Mutation_type"] == "INS"
)

qc = [
    {
        "Check": "SNP_directional_call_n",
        "Observed": call_snp_n,
        "Expected": call_snp_n,
        "Status": "PASS",
    },
    {
        "Check": "Patient_plot_SNP_Count_for_plot_sum",
        "Observed": plot_snp_n,
        "Expected": call_snp_n,
        "Status": "PASS" if plot_snp_n == call_snp_n else "FAIL",
    },
    {
        "Check": "DEL_event_n",
        "Observed": event_del_n,
        "Expected": event_del_n,
        "Status": "PASS",
    },
    {
        "Check": "Patient_plot_DEL_Count_for_plot_sum",
        "Observed": plot_del_n,
        "Expected": event_del_n,
        "Status": "PASS" if plot_del_n == event_del_n else "FAIL",
    },
    {
        "Check": "INS_event_n",
        "Observed": event_ins_n,
        "Expected": event_ins_n,
        "Status": "PASS",
    },
    {
        "Check": "Patient_plot_INS_Count_for_plot_sum",
        "Observed": plot_ins_n,
        "Expected": event_ins_n,
        "Status": "PASS" if plot_ins_n == event_ins_n else "FAIL",
    },
    {
        "Check": "Patient_plot_total_INDEL_event_count",
        "Observed": plot_del_n + plot_ins_n,
        "Expected": len(events),
        "Status":
            "PASS"
            if plot_del_n + plot_ins_n == len(events)
            else "FAIL",
    },
    {
        "Check": "Panel_patient_n",
        "Observed": len(panel_rows),
        "Expected": len(ranges),
        "Status": "PASS" if len(panel_rows) == len(ranges) else "FAIL",
    },
]

ready = all(r["Status"] == "PASS" for r in qc)

qc.append({
    "Check": "READY_FOR_FIGURE",
    "Observed": "YES" if ready else "NO",
    "Expected": "YES",
    "Status": "PASS" if ready else "FAIL",
})

write_tsv(
    QC_OUT,
    ["Check", "Observed", "Expected", "Status"],
    qc,
)


print("=" * 78)
print("FINAL PLOT TABLE PREPARATION COMPLETE")
print("=" * 78)
print()
print("No mutation or coordinate was changed.")
print("SNP plot counts = unique variant clone counts.")
print("DEL/INS plot counts = contiguous event occurrence counts.")
print()
print(f"SNP calls       : {call_snp_n}")
print(f"DEL events      : {event_del_n}")
print(f"INS events      : {event_ins_n}")
print(f"Total INDEL     : {len(events)}")
print(f"Patients        : {len(panel_rows)}")
print()
print("READY_FOR_FIGURE:", "YES" if ready else "NO")
print()
print("Outputs:")
for p in [
    PATIENT_PLOT,
    ALL_PLOT,
    PATIENT_SUMMARY,
    QC_OUT,
    PANEL_META,
    UNRESOLVED_OUT,
]:
    print(p)
