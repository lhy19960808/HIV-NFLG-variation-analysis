#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare final Figure 2 plotting-summary tables
====================================================

This script only summarizes the validated cross-patient mutation dataset.
It DOES NOT modify mutation calls, HXB2 coordinates, event definitions,
or patient-level structural annotations.

Figure 2 includes ONLY:
    SNP
    DEL
    INS

Explicitly excluded from Figure 2:
    No unique consensus
    P17 unresolved structural region
    P18 large missing region
    terminal no-callable regions

Primary cross-patient metric:
    Patient_n

Important:
- SNP recurrence refers to recurrence of a variable HXB2 locus/anchor across
  patients, not necessarily recurrence of the same nucleotide substitution.
- DEL/INS are shown as location markers; all current DEL/INS loci have
  Patient_n = 1.
- For insertions, cross-patient localization is by HXB2 anchor only.
"""

from pathlib import Path
import argparse
import csv
from collections import Counter

parser = argparse.ArgumentParser(
    description=(
        "Prepare final Figure 2 plotting-summary tables "
        "from the validated cross-patient dataset."
    )
)

parser.add_argument(
    "--cross-patient-dir",
    required=True,
    type=Path,
    help="Clean Step07 cross-patient dataset directory."
)

parser.add_argument(
    "--output-dir",
    required=True,
    type=Path,
    help="Output directory."
)

args = parser.parse_args()

ROOT = args.cross_patient_dir

SNP_FILE = ROOT / "01_SNP_locus_recurrence.tsv"
INDEL_FILE = ROOT / "02_INDEL_locus_recurrence.tsv"

OUT = args.output_dir
OUT.mkdir(parents=True, exist_ok=True)

OUT_ALL = OUT / "01_figure2_plotting_summary.tsv"
OUT_RECURRENT = OUT / "02_recurrent_SNP_sites.tsv"
OUT_SNP_SUMMARY = OUT / "03_SNP_recurrence_summary.tsv"
OUT_INDEL = OUT / "04_INDEL_plotting_summary.tsv"
OUT_QC = OUT / "05_figure2_plot_QC.tsv"
OUT_README = OUT / "README_figure2_plot_summary.md"

SNP_TYPES = [
    "A>C", "A>G", "A>T",
    "C>A", "C>G", "C>T",
    "G>A", "G>C", "G>T",
    "T>A", "T>C", "T>G",
]


def read_tsv(path):
    if not path.exists():
        raise FileNotFoundError(path)
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, fields, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})


snps = read_tsv(SNP_FILE)
indels = read_tsv(INDEL_FILE)

# ============================================================
# SNP final plotting rows
# ============================================================

snp_plot_rows = []
recurrent_rows = []

for r in snps:
    patient_n = int(r["Patient_n"])
    clone_n = int(r["Total_variant_clone_n"])

    substitutions = []
    substitution_patient_counts = []

    for s in SNP_TYPES:
        n = int(r[f"{s}_patient_n"])
        if n > 0:
            substitutions.append(s)
            substitution_patient_counts.append((s, n))

    substitution_text = ";".join(
        f"{s}:{n}"
        for s, n in substitution_patient_counts
    )

    if patient_n == 1:
        pattern = "SINGLE_PATIENT_SITE"
        recurrent_class = "Patient_n=1"
        same_substitution = "NA"
        display_priority = "BACKGROUND_SNP"
    else:
        same = any(n == patient_n for s, n in substitution_patient_counts)

        if same:
            pattern = "SAME_SUBSTITUTION"
            same_substitution = "YES"
        else:
            pattern = "MIXED_SUBSTITUTIONS"
            same_substitution = "NO"

        recurrent_class = f"Patient_n={patient_n}"
        display_priority = "HIGHLIGHT_RECURRENT_SNP"

    row = {
        "Record_type": "SNP",
        "Mutation_type": "SNP",
        "HXB2_locus_or_anchor": r["HXB2_locus_or_anchor"],
        "Locus_class": r["Locus_class"],
        "Plot_coordinate": r["Plot_coordinate"],
        "Patient_n": patient_n,
        "Patients": r["Patients"],
        "Total_variant_clone_n": clone_n,
        "Total_event_n": "",
        "Genes": r["Genes"],
        "Recurrent_class": recurrent_class,
        "Substitution_pattern": pattern,
        "Same_substitution_across_all_recurrent_patients": same_substitution,
        "Substitution_patient_counts": substitution_text,
        "Exact_coordinates_observed": r["Exact_coordinates_observed"],
        "Display_priority": display_priority,
        "Figure2_role":
            "SNP bar/line height = Patient_n",
        "Interpretation":
            "Cross-patient recurrence of within-host variability at this HXB2 locus/anchor."
    }

    snp_plot_rows.append(row)

    if patient_n >= 2:
        recurrent_rows.append(row.copy())


# ============================================================
# INDEL final plotting rows
# ============================================================

indel_plot_rows = []

for r in indels:
    patient_n = int(r["Patient_n"])

    if r["Mutation_type"] == "DEL":
        display_priority = "DEL_MARKER"
        role = "DEL marker at exact projected event locus"
    else:
        display_priority = "INS_MARKER"
        role = "INS marker at HXB2 anchor"

    row = {
        "Record_type": "INDEL_EVENT",
        "Mutation_type": r["Mutation_type"],
        "HXB2_locus_or_anchor": r["HXB2_locus_or_anchor"],
        "Locus_class": r["Aggregation_class"],
        "Plot_coordinate": r["Plot_coordinate"],
        "Patient_n": patient_n,
        "Patients": r["Patients"],
        "Total_variant_clone_n": r["Total_variant_clone_n"],
        "Total_event_n": r["Total_event_n"],
        "Genes": r["Genes"],
        "Recurrent_class": f"Patient_n={patient_n}",
        "Substitution_pattern": "NA",
        "Same_substitution_across_all_recurrent_patients": "NA",
        "Substitution_patient_counts": "NA",
        "Exact_coordinates_observed": r["HXB2_locus_or_anchor"],
        "Display_priority": display_priority,
        "Figure2_role": role,
        "Interpretation":
            "Location of clone-vs-patient-consensus indel variation; "
            "Patient_n is the cross-patient metric."
    }

    indel_plot_rows.append(row)


# ============================================================
# Combined final plotting table
# ============================================================

all_rows = snp_plot_rows + indel_plot_rows

all_rows.sort(
    key=lambda x: (
        float(x["Plot_coordinate"]),
        {"SNP": 0, "DEL": 1, "INS": 2}.get(x["Mutation_type"], 9)
    )
)

fields = [
    "Record_type",
    "Mutation_type",
    "HXB2_locus_or_anchor",
    "Locus_class",
    "Plot_coordinate",
    "Patient_n",
    "Patients",
    "Total_variant_clone_n",
    "Total_event_n",
    "Genes",
    "Recurrent_class",
    "Substitution_pattern",
    "Same_substitution_across_all_recurrent_patients",
    "Substitution_patient_counts",
    "Exact_coordinates_observed",
    "Display_priority",
    "Figure2_role",
    "Interpretation",
]

write_tsv(OUT_ALL, fields, all_rows)
write_tsv(OUT_RECURRENT, fields, recurrent_rows)
write_tsv(OUT_INDEL, fields, indel_plot_rows)


# ============================================================
# SNP recurrence summary
# ============================================================

patient_n_dist = Counter(int(r["Patient_n"]) for r in snps)

same_n = sum(
    r["Substitution_pattern"] == "SAME_SUBSTITUTION"
    for r in recurrent_rows
)
mixed_n = sum(
    r["Substitution_pattern"] == "MIXED_SUBSTITUTIONS"
    for r in recurrent_rows
)

summary_rows = []

for n in sorted(patient_n_dist):
    summary_rows.append({
        "Category": f"SNP_loci_with_Patient_n_{n}",
        "Count": patient_n_dist[n],
        "Definition":
            f"Number of SNP loci/anchors observed as variable in exactly {n} unique patient(s)."
    })

summary_rows += [
    {
        "Category": "Total_SNP_loci_or_anchors",
        "Count": len(snps),
        "Definition": "All Figure 2 SNP loci/anchors."
    },
    {
        "Category": "Recurrent_SNP_loci_Patient_n_ge_2",
        "Count": len(recurrent_rows),
        "Definition": "SNP loci/anchors variable in at least two patients."
    },
    {
        "Category": "Recurrent_same_substitution_loci",
        "Count": same_n,
        "Definition":
            "Recurrent SNP loci where at least one substitution direction is present in all recurrent patients."
    },
    {
        "Category": "Recurrent_mixed_substitution_loci",
        "Count": mixed_n,
        "Definition":
            "Recurrent SNP loci where no single substitution direction is shared by all recurrent patients."
    },
]

write_tsv(
    OUT_SNP_SUMMARY,
    ["Category", "Count", "Definition"],
    summary_rows
)


# ============================================================
# QC
# ============================================================

del_rows = [r for r in indel_plot_rows if r["Mutation_type"] == "DEL"]
ins_rows = [r for r in indel_plot_rows if r["Mutation_type"] == "INS"]

bad_types = [
    r for r in all_rows
    if r["Mutation_type"] not in {"SNP", "DEL", "INS"}
]

patient_n_mismatch = []

for r in all_rows:
    patients = [
        x.strip()
        for x in r["Patients"].split(",")
        if x.strip()
    ]
    if len(set(patients)) != int(r["Patient_n"]):
        patient_n_mismatch.append(r)

qc = [
    {
        "Check": "SNP_locus_or_anchor_n",
        "Observed": len(snp_plot_rows),
        "Expected": 973,
        "Status": "PASS" if len(snp_plot_rows) == 973 else "FAIL",
    },
    {
        "Check": "Recurrent_SNP_locus_n",
        "Observed": len(recurrent_rows),
        "Expected": 52,
        "Status": "PASS" if len(recurrent_rows) == 52 else "FAIL",
    },
    {
        "Check": "Same_substitution_recurrent_locus_n",
        "Observed": same_n,
        "Expected": 21,
        "Status": "PASS" if same_n == 21 else "FAIL",
    },
    {
        "Check": "Mixed_substitution_recurrent_locus_n",
        "Observed": mixed_n,
        "Expected": 31,
        "Status": "PASS" if mixed_n == 31 else "FAIL",
    },
    {
        "Check": "DEL_locus_n",
        "Observed": len(del_rows),
        "Expected": 32,
        "Status": "PASS" if len(del_rows) == 32 else "FAIL",
    },
    {
        "Check": "INS_anchor_n",
        "Observed": len(ins_rows),
        "Expected": 13,
        "Status": "PASS" if len(ins_rows) == 13 else "FAIL",
    },
    {
        "Check": "DEL_all_Patient_n_1",
        "Observed": max(int(r["Patient_n"]) for r in del_rows) if del_rows else 0,
        "Expected": 1,
        "Status": "PASS" if del_rows and all(int(r["Patient_n"]) == 1 for r in del_rows) else "FAIL",
    },
    {
        "Check": "INS_all_Patient_n_1",
        "Observed": max(int(r["Patient_n"]) for r in ins_rows) if ins_rows else 0,
        "Expected": 1,
        "Status": "PASS" if ins_rows and all(int(r["Patient_n"]) == 1 for r in ins_rows) else "FAIL",
    },
    {
        "Check": "Only_SNP_DEL_INS_in_Figure2",
        "Observed": len(bad_types),
        "Expected": 0,
        "Status": "PASS" if len(bad_types) == 0 else "FAIL",
    },
    {
        "Check": "Patient_n_matches_patient_list",
        "Observed": len(patient_n_mismatch),
        "Expected": 0,
        "Status": "PASS" if len(patient_n_mismatch) == 0 else "FAIL",
    },
]

ready = all(r["Status"] == "PASS" for r in qc)

qc.append({
    "Check": "READY_FOR_FIGURE2_PLOTTING",
    "Observed": "YES" if ready else "NO",
    "Expected": "YES",
    "Status": "PASS" if ready else "FAIL",
})

write_tsv(
    OUT_QC,
    ["Check", "Observed", "Expected", "Status"],
    qc
)


# ============================================================
# README
# ============================================================

readme = f"""# Figure 2 final plotting summary

## Figure 2 scope

Figure 2 contains only:
- SNP
- DEL
- INS

The following are excluded from Figure 2:
- No unique consensus
- P17 unresolved structural region
- P18 HXB2 3487–5903 large missing region
- terminal no-callable regions

These exclusions do not modify the underlying mutation calls.

## Primary metric

`Patient_n`

This is the number of unique patients showing clone-vs-patient-consensus
variation at a given HXB2 locus/anchor.

## SNP recurrence

Total SNP loci/anchors: {len(snps)}

- Patient_n = 1: {patient_n_dist.get(1, 0)}
- Patient_n = 2: {patient_n_dist.get(2, 0)}
- Patient_n = 3: {patient_n_dist.get(3, 0)}
- Recurrent SNP loci (Patient_n >= 2): {len(recurrent_rows)}

Among recurrent SNP loci:
- same-substitution loci: {same_n}
- mixed-substitution loci: {mixed_n}

A recurrent SNP locus therefore means a recurrently variable genomic position,
not necessarily an identical recurrent nucleotide substitution.

## INDEL

DEL loci: {len(del_rows)}
INS anchors: {len(ins_rows)}

All current DEL and INS loci have `Patient_n = 1`.

For Figure 2:
- DEL should be represented as a location marker at its projected event locus.
- INS should be represented as a location marker at its HXB2 anchor.
- Do not use `Total_event_n` as the cross-patient height.

## Special case: INS 7289|7290

CN2024AH2-6 has:
- Patient_n = 1
- Total_variant_clone_n = 3
- Total_event_n = 6

Three clones each contain two separate patient-consensus-relative insertion
blocks, both mapping to the same HXB2 anchor 7289|7290. Therefore the Figure 2
representation remains a single INS anchor marker with Patient_n = 1.

## Plotting recommendation

SNP:
- x = Plot_coordinate
- y = Patient_n
- Patient_n=1 = background/single-patient variable site
- Patient_n>=2 = recurrent variable site highlight

DEL:
- marker only, exact projected event locus

INS:
- marker only, HXB2 anchor

No no-unique-consensus or patient-level missing-region track should be added
to Figure 2.
"""

OUT_README.write_text(readme, encoding="utf-8")

print("=" * 84)
print("FIGURE 2 PLOTTING SUMMARY COMPLETE")
print("=" * 84)
print(f"SNP loci/anchors             : {len(snp_plot_rows)}")
print(f"Recurrent SNP loci >=2       : {len(recurrent_rows)}")
print(f"Same-substitution recurrent  : {same_n}")
print(f"Mixed-substitution recurrent : {mixed_n}")
print(f"DEL loci                     : {len(del_rows)}")
print(f"INS anchors                  : {len(ins_rows)}")
print()
print("READY_FOR_FIGURE2_PLOTTING:", "YES" if ready else "NO")
print()
print("Output directory:")
print(OUT)
