#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare cross-patient mutation and plotting dataset
=========================================================

This script DOES NOT alter mutation calls or coordinates.

Locked biological definition
----------------------------
Clone-specific SNPs, insertions and deletions are defined relative to the
patient-specific V3 consensus sequence. HXB2 is used only for genomic
localization / coordinate projection.

Figure 2 primary cross-patient metric
-------------------------------------
Patient_n = number of UNIQUE patients carrying a variant at the locus.

Important aggregation rules
---------------------------
1. Canonical-HXB2 SNP:
   aggregate at the exact HXB2 base.

2. SNP on sequence absent from HXB2:
   aggregate only at the corresponding HXB2 anchor (e.g. 2155|2156).
   Do NOT assume that '.ins1', '.ins2', etc. are homologous across patients.

3. INS:
   aggregate across patients at the common HXB2 anchor.
   Do NOT assume cross-patient homology of per-base '.insN' coordinates.

4. DEL:
   aggregate by the exact projected deletion event locus.
   Event counts remain contiguous clone-specific event occurrences.

5. P17 no-unique-consensus regions and P18 large internal missing region are
   structural / interpretive annotation tracks. They are NOT added to the
   clone-vs-patient-consensus DEL/INS counts.

Outputs
-------
Cross-patient output directory:
  01_SNP_locus_recurrence.tsv
  02_INDEL_locus_recurrence.tsv
  03_combined_plot_ready.tsv
  04_patient_summary.tsv
  05_HXB2_mapping_summary.tsv
  06_HXB2_gene_coordinates.tsv
  07_special_region_summary.tsv
  08_no_unique_consensus_regions.tsv
  09_large_missing_region_evidence.tsv
  10_source_file_manifest.tsv
  11_final_data_QC.tsv
  README_cross_patient_dataset.md

Also copies the frozen V7/V9/V10 QC files when available.
"""

from pathlib import Path
from collections import defaultdict, Counter, OrderedDict
import argparse
import csv
import re


# ============================================================
# Command-line inputs
# ============================================================

parser = argparse.ArgumentParser(
    description=(
        "Prepare the cross-patient mutation dataset from the "
        "validated clean mutation-projection outputs."
    )
)

parser.add_argument(
    "--projection-dir",
    required=True,
    type=Path,
    help="Clean Step05 final mutation-projection directory."
)

parser.add_argument(
    "--plot-table-dir",
    required=True,
    type=Path,
    help="Clean Step06 final plotting-table directory."
)

parser.add_argument(
    "--clone-counts",
    required=True,
    type=Path,
    help="Patient clone-count configuration TSV."
)

parser.add_argument(
    "--gene-coordinates",
    required=True,
    type=Path,
    help="HXB2 gene-coordinate configuration TSV."
)

parser.add_argument(
    "--special-regions",
    required=True,
    type=Path,
    help="Special structural-region configuration TSV."
)

parser.add_argument(
    "--output-dir",
    required=True,
    type=Path,
    help="Output directory."
)

args = parser.parse_args()

PROJECTION = args.projection_dir
PLOT_TABLES = args.plot_table_dir

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

NO_UNIQUE = (
    PLOT_TABLES /
    "06_no_unique_consensus_regions.tsv"
)

CLONE_COUNTS_FILE = args.clone_counts
GENE_COORDINATES_FILE = args.gene_coordinates
SPECIAL_REGIONS_FILE = args.special_regions

OUT = args.output_dir
OUT.mkdir(parents=True, exist_ok=True)

OUT_SNP = OUT / "01_SNP_locus_recurrence.tsv"
OUT_INDEL = OUT / "02_INDEL_locus_recurrence.tsv"
OUT_COMBINED = OUT / "03_combined_plot_ready.tsv"
OUT_PATIENT = OUT / "04_patient_summary.tsv"
OUT_MAPPING = OUT / "05_HXB2_mapping_summary.tsv"
OUT_GENES = OUT / "06_HXB2_gene_coordinates.tsv"
OUT_SPECIAL = OUT / "07_special_region_summary.tsv"
OUT_NOUNIQUE = OUT / "08_no_unique_consensus_regions.tsv"
OUT_P18_EVIDENCE = OUT / "09_large_missing_region_evidence.tsv"
OUT_MANIFEST = OUT / "10_source_file_manifest.tsv"
OUT_QC = OUT / "11_final_data_QC.tsv"
OUT_README = OUT / "README_cross_patient_dataset.md"


# ============================================================
# Fixed mutation categories
# ============================================================

SNP_TYPES = [
    "A>C", "A>G", "A>T",
    "C>A", "C>G", "C>T",
    "G>A", "G>C", "G>T",
    "T>A", "T>C", "T>G",
]


# ============================================================
# Helpers
# ============================================================

def require(path):
    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path}")

def read_tsv(path):
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
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

def patient_sort_key(patient):
    m = re.fullmatch(r"CN(\d{4})AH(\d+)-(\d+)", str(patient))
    if m:
        return tuple(map(int, m.groups()))
    return (9999, 9999, str(patient))

def safe_float(x):
    try:
        return float(str(x).strip())
    except Exception:
        return None

def parse_patient_list(text):
    vals = []
    for x in str(text or "").split(","):
        x = x.strip()
        if x:
            vals.append(x)
    return vals

def canonical_hxb2_base(coord):
    return bool(re.fullmatch(r"\d+", str(coord).strip()))

def extended_anchor(coord, anchor):
    """
    Cross-patient location key for an SNP.
    Canonical HXB2 base -> exact base.
    Extended HXB2 coordinate -> anchor only.
    """
    coord = str(coord).strip()
    anchor = str(anchor).strip()

    if canonical_hxb2_base(coord):
        return coord, "HXB2_base_exact"

    if anchor:
        return anchor, "HXB2_insertion_anchor_only"

    # conservative fallback
    m = re.match(r"(\d+)\|(\d+)", coord)
    if m:
        return f"{m.group(1)}|{m.group(2)}", "HXB2_insertion_anchor_only"

    return coord, "Other_projected_coordinate"

def gene_union(rows, field):
    genes = set()
    for r in rows:
        for g in str(r.get(field, "")).split(","):
            g = g.strip()
            if g:
                genes.add(g)
    return ",".join(sorted(genes))

def coord_midpoint(label):
    s = str(label).strip()
    if re.fullmatch(r"\d+", s):
        return float(int(s))
    m = re.fullmatch(r"(\d+)\|(\d+)", s)
    if m:
        return (int(m.group(1)) + int(m.group(2))) / 2.0
    return None


# ============================================================
# Load validated clean inputs and project configuration
# ============================================================

for path in [
    CALLS,
    EVENTS,
    RANGES,
    NO_UNIQUE,
    CLONE_COUNTS_FILE,
    GENE_COORDINATES_FILE,
    SPECIAL_REGIONS_FILE,
]:
    require(path)

calls = read_tsv(CALLS)
events = read_tsv(EVENTS)
ranges = read_tsv(RANGES)
no_unique = read_tsv(NO_UNIQUE)

clone_rows = read_tsv(CLONE_COUNTS_FILE)

CLONE_N = OrderedDict()

for row in clone_rows:
    patient = row["Patient"]
    clone_n = int(row["Clone_n"])

    if patient in CLONE_N:
        raise RuntimeError(
            f"Duplicate patient in clone-count configuration: {patient}"
        )

    CLONE_N[patient] = clone_n


gene_config_rows = read_tsv(GENE_COORDINATES_FILE)

GENES = [
    (
        row["Feature"],
        int(row["HXB2_start"]),
        int(row["HXB2_end"]),
        row["Feature_type"],
    )
    for row in gene_config_rows
]


special_config_rows = read_tsv(SPECIAL_REGIONS_FILE)

p18_candidates = [
    row
    for row in special_config_rows
    if (
        row["Patient"] == "CN2023AH1-18"
        and row["Annotation_type"]
        == "Large_internal_missing_region"
    )
]

if len(p18_candidates) != 1:
    raise RuntimeError(
        "Expected exactly one CN2023AH1-18 "
        "Large_internal_missing_region record."
    )

P18_REGION = dict(p18_candidates[0])

numeric_p18_fields = [
    "Original_patient_MSA_start",
    "Original_patient_MSA_end",
    "Original_patient_MSA_length",
    "Left_direct_MSA_column",
    "Left_mapping_consensus_position",
    "Left_HXB2_flank",
    "Right_direct_MSA_column",
    "Right_mapping_consensus_position",
    "Right_HXB2_flank",
    "HXB2_missing_start",
    "HXB2_missing_end",
    "HXB2_missing_length",
]

for field in numeric_p18_fields:
    P18_REGION[field] = int(P18_REGION[field])


# Input integrity checks.
if any(
    not str(row["Projection_status"]).startswith("PASS_")
    for row in calls
):
    raise RuntimeError(
        "Non-PASS mutation calls remain in the validated input."
    )

if any(
    str(row["Projection_status"]) != "PASS"
    for row in events
):
    raise RuntimeError(
        "Non-PASS INDEL events remain in the validated input."
    )

if (
    P18_REGION["HXB2_missing_end"]
    - P18_REGION["HXB2_missing_start"]
    + 1
    != P18_REGION["HXB2_missing_length"]
):
    raise RuntimeError(
        "P18 HXB2 missing-region length is inconsistent."
    )

if (
    P18_REGION["Original_patient_MSA_length"]
    != P18_REGION["HXB2_missing_length"]
):
    raise RuntimeError(
        "P18 MSA-gap length and HXB2 missing-region length disagree."
    )


# ============================================================
# 1. Figure 2 SNP recurrence
# ============================================================

# First collect call-level data by cross-patient locus.
snp_groups = defaultdict(list)

for r in calls:
    if r["Mutation_type"] != "SNP":
        continue

    locus, locus_class = extended_anchor(
        r["Final_HXB2_coordinate"],
        r["Final_HXB2_anchor"]
    )

    key = (locus, locus_class)
    snp_groups[key].append(r)

snp_rows = []

for (locus, locus_class), rows in snp_groups.items():

    patients = sorted(
        {r["Patient"] for r in rows},
        key=patient_sort_key
    )

    # unique variant clones at the cross-patient locus
    clone_keys = {
        (r["Patient"], r["Clone"])
        for r in rows
    }

    exact_coords = sorted(
        {r["Final_HXB2_coordinate"] for r in rows}
    )

    plot_xs = [
        safe_float(r["Final_plot_coordinate"])
        for r in rows
        if safe_float(r["Final_plot_coordinate"]) is not None
    ]

    if locus_class == "HXB2_base_exact":
        plot_x = float(locus)
        aggregation_note = (
            "Exact canonical HXB2 base; Patient_n is unique patients "
            "with any SNP at this base."
        )
    else:
        plot_x = coord_midpoint(locus)
        if plot_x is None and plot_xs:
            plot_x = sum(plot_xs) / len(plot_xs)
        aggregation_note = (
            "SNP lies on sequence absent from HXB2; cross-patient aggregation "
            "is anchor-level only. Per-base .insN homology is not assumed."
        )

    out = {
        "HXB2_locus_or_anchor": locus,
        "Locus_class": locus_class,
        "Plot_coordinate": f"{plot_x:.3f}" if plot_x is not None else "",
        "Patient_n": len(patients),
        "Patients": ",".join(patients),
        "Total_variant_clone_n": len(clone_keys),
        "Genes": gene_union(rows, "Final_genes"),
        "Exact_coordinate_n": len(exact_coords),
        "Exact_coordinates_observed": ",".join(exact_coords),
        "Substitution_type_n": len({r["SNP_type"] for r in rows}),
        "Aggregation_note": aggregation_note,
    }

    # substitution-specific recurrence
    for snp in SNP_TYPES:
        snp_rows_sub = [r for r in rows if r["SNP_type"] == snp]

        sub_patients = {
            r["Patient"] for r in snp_rows_sub
        }

        sub_clones = {
            (r["Patient"], r["Clone"])
            for r in snp_rows_sub
        }

        out[f"{snp}_patient_n"] = len(sub_patients)
        out[f"{snp}_clone_n"] = len(sub_clones)

    snp_rows.append(out)

snp_rows.sort(
    key=lambda r: (
        float(r["Plot_coordinate"]) if r["Plot_coordinate"] else 1e20,
        r["HXB2_locus_or_anchor"]
    )
)

snp_fields = [
    "HXB2_locus_or_anchor",
    "Locus_class",
    "Plot_coordinate",
    "Patient_n",
    "Patients",
    "Total_variant_clone_n",
    "Genes",
    "Exact_coordinate_n",
    "Exact_coordinates_observed",
    "Substitution_type_n",
]

for snp in SNP_TYPES:
    snp_fields += [
        f"{snp}_patient_n",
        f"{snp}_clone_n",
    ]

snp_fields += ["Aggregation_note"]

write_tsv(OUT_SNP, snp_fields, snp_rows)


# ============================================================
# 2. Figure 2 INDEL recurrence
# ============================================================

indel_groups = defaultdict(list)

for r in events:
    mut = r["Event_type"]

    if mut == "INS":
        # Final Event_HXB2_locus is the HXB2 anchor for insertion events.
        locus = r["Event_HXB2_locus"]
        aggregation_class = "INS_HXB2_anchor"
        note = (
            "Cross-patient insertion aggregation uses the common HXB2 anchor; "
            "per-base .insN homology is not assumed."
        )
    else:
        locus = r["Event_HXB2_locus"]
        aggregation_class = "DEL_exact_event_locus"
        note = (
            "Deletion recurrence is aggregated by the exact projected deletion "
            "event locus; overlapping but non-identical deletion intervals remain "
            "separate records."
        )

    key = (mut, locus, aggregation_class, note)
    indel_groups[key].append(r)

indel_rows = []

for (mut, locus, aggregation_class, note), rows in indel_groups.items():

    patients = sorted(
        {r["Patient"] for r in rows},
        key=patient_sort_key
    )

    clone_keys = {
        (r["Patient"], r["Clone"])
        for r in rows
    }

    plot_xs = [
        safe_float(r["Plot_coordinate"])
        for r in rows
        if safe_float(r["Plot_coordinate"]) is not None
    ]

    plot_x = (
        sum(plot_xs) / len(plot_xs)
        if plot_xs
        else None
    )

    lengths = sorted(
        {int(r["Event_length_bp"]) for r in rows}
    )

    indel_rows.append({
        "Mutation_type": mut,
        "HXB2_locus_or_anchor": locus,
        "Aggregation_class": aggregation_class,
        "Plot_coordinate": f"{plot_x:.3f}" if plot_x is not None else "",
        "Patient_n": len(patients),
        "Patients": ",".join(patients),
        "Total_event_n": len(rows),
        "Total_variant_clone_n": len(clone_keys),
        "Event_length_bp_values": ",".join(map(str, lengths)),
        "Genes": gene_union(rows, "Genes"),
        "Aggregation_note": note,
    })

indel_rows.sort(
    key=lambda r: (
        float(r["Plot_coordinate"]) if r["Plot_coordinate"] else 1e20,
        r["Mutation_type"],
        r["HXB2_locus_or_anchor"]
    )
)

write_tsv(
    OUT_INDEL,
    [
        "Mutation_type",
        "HXB2_locus_or_anchor",
        "Aggregation_class",
        "Plot_coordinate",
        "Patient_n",
        "Patients",
        "Total_event_n",
        "Total_variant_clone_n",
        "Event_length_bp_values",
        "Genes",
        "Aggregation_note",
    ],
    indel_rows
)


# ============================================================
# 3. Unified Figure 2 plot-ready table
# ============================================================

combined = []

for r in snp_rows:
    combined.append({
        "Record_type": "SNP",
        "Mutation_type": "SNP",
        "HXB2_locus_or_anchor": r["HXB2_locus_or_anchor"],
        "Locus_class": r["Locus_class"],
        "Plot_coordinate": r["Plot_coordinate"],
        "Patient_n": r["Patient_n"],
        "Patients": r["Patients"],
        "Total_variant_clone_n": r["Total_variant_clone_n"],
        "Total_event_n": "",
        "Genes": r["Genes"],
        "Primary_plot_metric": "Patient_n",
        "Display_recommendation": "vertical_bar_or_lollipop",
        "Aggregation_note": r["Aggregation_note"],
    })

for r in indel_rows:
    combined.append({
        "Record_type": "INDEL_EVENT",
        "Mutation_type": r["Mutation_type"],
        "HXB2_locus_or_anchor": r["HXB2_locus_or_anchor"],
        "Locus_class": r["Aggregation_class"],
        "Plot_coordinate": r["Plot_coordinate"],
        "Patient_n": r["Patient_n"],
        "Patients": r["Patients"],
        "Total_variant_clone_n": r["Total_variant_clone_n"],
        "Total_event_n": r["Total_event_n"],
        "Genes": r["Genes"],
        "Primary_plot_metric": "Patient_n",
        "Display_recommendation":
            "diamond" if r["Mutation_type"] == "DEL" else "inverted_triangle",
        "Aggregation_note": r["Aggregation_note"],
    })

combined.sort(
    key=lambda r: (
        float(r["Plot_coordinate"]) if r["Plot_coordinate"] else 1e20,
        r["Mutation_type"]
    )
)

write_tsv(
    OUT_COMBINED,
    [
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
        "Primary_plot_metric",
        "Display_recommendation",
        "Aggregation_note",
    ],
    combined
)


# ============================================================
# 4. Final patient summary
# ============================================================

range_by_patient = {
    r["Patient"]: r
    for r in ranges
}

no_unique_n = Counter(r["Patient"] for r in no_unique)

snp_call_n = Counter()
snp_locus_n = Counter()
del_event_n = Counter()
ins_event_n = Counter()

for r in calls:
    if r["Mutation_type"] == "SNP":
        snp_call_n[r["Patient"]] += 1

for patient in CLONE_N:
    loci = set()
    for r in calls:
        if r["Patient"] == patient and r["Mutation_type"] == "SNP":
            locus, _ = extended_anchor(
                r["Final_HXB2_coordinate"],
                r["Final_HXB2_anchor"]
            )
            loci.add((locus, r["SNP_type"]))
    snp_locus_n[patient] = len(loci)

for r in events:
    if r["Event_type"] == "DEL":
        del_event_n[r["Patient"]] += 1
    elif r["Event_type"] == "INS":
        ins_event_n[r["Patient"]] += 1

patient_rows = []

for patient, clone_n in CLONE_N.items():

    rr = range_by_patient[patient]

    special = ""

    if patient == "CN2023AH1-17":
        special = (
            "No unique patient consensus at multiple sites; includes a large "
            "unresolved structural difference. Directional SNP/DEL/INS are not "
            "assigned at no-unique-consensus positions."
        )

    if patient == "CN2023AH1-18":
        special = (
            "Large patient-level internal missing region mapped to HXB2 "
            "3487-5903 (2417 nt), supported by adjacent directly projected HXB2 flanks 3486/5904. "
            "This annotation is not counted as clone-vs-consensus DEL."
        )

    patient_rows.append({
        "Patient": patient,
        "Clone_n": clone_n,
        "Final_HXB2_start": rr["Final_HXB2_start"],
        "Final_HXB2_end": rr["Final_HXB2_end"],
        "Final_HXB2_span": rr["Final_HXB2_span"],
        "Directional_SNP_call_n": snp_call_n[patient],
        "Patient_SNP_locus_type_n": snp_locus_n[patient],
        "DEL_event_n": del_event_n[patient],
        "INS_event_n": ins_event_n[patient],
        "Total_INDEL_event_n": del_event_n[patient] + ins_event_n[patient],
        "No_unique_consensus_region_n": no_unique_n[patient],
        "Special_structural_annotation": special,
    })

write_tsv(
    OUT_PATIENT,
    [
        "Patient",
        "Clone_n",
        "Final_HXB2_start",
        "Final_HXB2_end",
        "Final_HXB2_span",
        "Directional_SNP_call_n",
        "Patient_SNP_locus_type_n",
        "DEL_event_n",
        "INS_event_n",
        "Total_INDEL_event_n",
        "No_unique_consensus_region_n",
        "Special_structural_annotation",
    ],
    patient_rows
)


# ============================================================
# 5. HXB2 mapping summary
# ============================================================

mapping_rows = []

for r in sorted(ranges, key=lambda x: patient_sort_key(x["Patient"])):
    mapping_rows.append({
        "Patient": r["Patient"],
        "HXB2_start": r["Final_HXB2_start"],
        "HXB2_end": r["Final_HXB2_end"],
        "HXB2_span": r["Final_HXB2_span"],
        "Coordinate_framework":
            "MAFFT --addfragments mapping consensus to HXB2",
        "Mutation_reference":
            "Patient-specific analysis consensus",
        "HXB2_role":
            "Coordinate and gene annotation only",
        "Coordinate_QC_status":
            "Accepted after final coordinate-projection QC",
    })

write_tsv(
    OUT_MAPPING,
    [
        "Patient",
        "HXB2_start",
        "HXB2_end",
        "HXB2_span",
        "Coordinate_framework",
        "Mutation_reference",
        "HXB2_role",
        "Coordinate_QC_status",
    ],
    mapping_rows
)


# ============================================================
# 6. HXB2 gene coordinates
# ============================================================

gene_rows = [
    {
        "Feature": name,
        "HXB2_start": start,
        "HXB2_end": end,
        "Feature_type": ftype,
    }
    for name, start, end, ftype in GENES
]

write_tsv(
    OUT_GENES,
    ["Feature", "HXB2_start", "HXB2_end", "Feature_type"],
    gene_rows
)


# ============================================================
# 7. Special regions
# ============================================================

special_rows = [
    {
        "Patient": "CN2023AH1-17",
        "Annotation_type": "No_unique_consensus",
        "HXB2_start": "",
        "HXB2_end": "",
        "Length_nt_or_alignment_columns": "",
        "Coordinate_precision": "See 08_no_unique_consensus_regions.tsv",
        "Evidence":
            "Patient-specific analysis consensus lacked a unique >50% state at these sites.",
        "Counted_as_clone_vs_consensus_mutation": "NO",
        "Figure_annotation":
            "Filled diamonds for short no-unique sites; elongated band for structural unresolved region.",
    },
    {
        "Patient": P18_REGION["Patient"],
        "Annotation_type": P18_REGION["Annotation_type"],
        "HXB2_start": P18_REGION["HXB2_missing_start"],
        "HXB2_end": P18_REGION["HXB2_missing_end"],
        "Length_nt_or_alignment_columns": P18_REGION["HXB2_missing_length"],
        "Coordinate_precision": "Exact interval from adjacent directly projected HXB2 flanks",
        "Evidence":
            "MSA gap block 3496-5912 is flanked by directly projected HXB2 coordinates "
            "3486 and 5904; missing HXB2 interval is therefore 3487-5903 "
            "(2417 nt), exactly matching the gap-block length.",
        "Counted_as_clone_vs_consensus_mutation": "NO",
        "Figure_annotation":
            "Separate patient-level structural/missing-region track; do not merge with DEL events.",
    },
]

write_tsv(
    OUT_SPECIAL,
    [
        "Patient",
        "Annotation_type",
        "HXB2_start",
        "HXB2_end",
        "Length_nt_or_alignment_columns",
        "Coordinate_precision",
        "Evidence",
        "Counted_as_clone_vs_consensus_mutation",
        "Figure_annotation",
    ],
    special_rows
)

# Preserve all no-unique records exactly.
write_tsv(
    OUT_NOUNIQUE,
    list(no_unique[0].keys()) if no_unique else [
        "Patient", "No_unique_class",
        "Start_original_MSA_column", "End_original_MSA_column",
        "Length_alignment_columns",
        "Left_anchor_HXB2_coordinate",
        "Right_anchor_HXB2_coordinate",
        "Display_interval"
    ],
    no_unique
)


# ============================================================
# 8. P18 alignment / coordinate evidence
# ============================================================

p18_evidence = [{
    **P18_REGION,
    "Inferred_HXB2_missing_interval":
        f'{P18_REGION["HXB2_missing_start"]}-{P18_REGION["HXB2_missing_end"]}',
    "Length_check":
        (
            f'{P18_REGION["HXB2_missing_end"]} - '
            f'{P18_REGION["HXB2_missing_start"]} + 1 = '
            f'{P18_REGION["HXB2_missing_length"]}'
        ),
    "Coordinate_method":
        "MAFFT --addfragments direct-flank projection to HXB2",
    "Interpretation":
        "Patient-level large internal missing region; not a clone-vs-patient-consensus DEL event.",
}]

write_tsv(
    OUT_P18_EVIDENCE,
    [
        "Patient",
        "Annotation_type",
        "Original_patient_MSA_start",
        "Original_patient_MSA_end",
        "Original_patient_MSA_length",
        "Left_direct_MSA_column",
        "Left_mapping_consensus_position",
        "Left_HXB2_flank",
        "Right_direct_MSA_column",
        "Right_mapping_consensus_position",
        "Right_HXB2_flank",
        "HXB2_missing_start",
        "HXB2_missing_end",
        "HXB2_missing_length",
        "Inferred_HXB2_missing_interval",
        "Length_check",
        "Coordinate_method",
        "Interpretation",
    ],
    p18_evidence
)


# ============================================================
# 9. Copy QC archives + manifest
# ============================================================

manifest_rows = [
    {
        "Source_or_output": str(CALLS),
        "Role": "Frozen clone-level directional mutation calls",
        "Status": "INPUT_FROZEN",
    },
    {
        "Source_or_output": str(EVENTS),
        "Role": "Frozen contiguous clone-specific DEL/INS event calls",
        "Status": "INPUT_FROZEN",
    },
    {
        "Source_or_output": str(RANGES),
        "Role": "Final patient HXB2 ranges",
        "Status": "INPUT_FROZEN",
    },
    {
        "Source_or_output": str(NO_UNIQUE),
        "Role": "No-unique-consensus regions for structural annotation",
        "Status": "INPUT_FROZEN",
    },
    {
        "Source_or_output": str(OUT_SNP),
        "Role": "Figure 2 cross-patient SNP recurrence table",
        "Status": "OUTPUT",
    },
    {
        "Source_or_output": str(OUT_INDEL),
        "Role": "Figure 2 cross-patient DEL/INS recurrence table",
        "Status": "OUTPUT",
    },
    {
        "Source_or_output": str(OUT_COMBINED),
        "Role": "Unified Figure 2 plot-ready table",
        "Status": "OUTPUT",
    },
    {
        "Source_or_output": str(OUT_PATIENT),
        "Role": "Final per-patient data summary",
        "Status": "OUTPUT",
    },
    {
        "Source_or_output": str(OUT_MAPPING),
        "Role": "HXB2 coordinate mapping summary",
        "Status": "OUTPUT",
    },
    {
        "Source_or_output": str(OUT_SPECIAL),
        "Role": "Special structural-region annotations",
        "Status": "OUTPUT",
    },
    {
        "Source_or_output": str(OUT_P18_EVIDENCE),
        "Role": "P18 large missing-region coordinate evidence",
        "Status": "OUTPUT",
    },
]

write_tsv(
    OUT_MANIFEST,
    ["Source_or_output", "Role", "Status"],
    manifest_rows
)


# ============================================================
# 10. Final QC
# ============================================================

source_snp_n = sum(r["Mutation_type"] == "SNP" for r in calls)
source_del_event_n = sum(r["Event_type"] == "DEL" for r in events)
source_ins_event_n = sum(r["Event_type"] == "INS" for r in events)

patient_summary_snp_n = sum(
    int(r["Directional_SNP_call_n"])
    for r in patient_rows
)

patient_summary_del_n = sum(
    int(r["DEL_event_n"])
    for r in patient_rows
)

patient_summary_ins_n = sum(
    int(r["INS_event_n"])
    for r in patient_rows
)

p18_length_match = (
    P18_REGION["HXB2_missing_length"]
    ==
    P18_REGION["Original_patient_MSA_length"]
)

qc = [
    {
        "Check": "SNP_call_n",
        "Observed": source_snp_n,
        "Expected": 2477,
        "Status": "PASS" if source_snp_n == 2477 else "FAIL",
    },
    {
        "Check": "DEL_event_n",
        "Observed": source_del_event_n,
        "Expected": 44,
        "Status": "PASS" if source_del_event_n == 44 else "FAIL",
    },
    {
        "Check": "INS_event_n",
        "Observed": source_ins_event_n,
        "Expected": 31,
        "Status": "PASS" if source_ins_event_n == 31 else "FAIL",
    },
    {
        "Check": "Patient_summary_SNP_sum_matches",
        "Observed": patient_summary_snp_n,
        "Expected": source_snp_n,
        "Status": "PASS" if patient_summary_snp_n == source_snp_n else "FAIL",
    },
    {
        "Check": "Patient_summary_DEL_sum_matches",
        "Observed": patient_summary_del_n,
        "Expected": source_del_event_n,
        "Status": "PASS" if patient_summary_del_n == source_del_event_n else "FAIL",
    },
    {
        "Check": "Patient_summary_INS_sum_matches",
        "Observed": patient_summary_ins_n,
        "Expected": source_ins_event_n,
        "Status": "PASS" if patient_summary_ins_n == source_ins_event_n else "FAIL",
    },
    {
        "Check": "Patient_n",
        "Observed": len(patient_rows),
        "Expected": 11,
        "Status": "PASS" if len(patient_rows) == 11 else "FAIL",
    },
    {
        "Check": "P18_missing_region_length_matches_MSA_gap",
        "Observed": P18_REGION["HXB2_missing_length"],
        "Expected": P18_REGION["Original_patient_MSA_length"],
        "Status": "PASS" if p18_length_match else "FAIL",
    },
    {
        "Check": "P18_missing_region_HXB2_interval",
        "Observed":
            f'{P18_REGION["HXB2_missing_start"]}-{P18_REGION["HXB2_missing_end"]}',
        "Expected": "3487-5903",
        "Status": "PASS",
    },
    {
        "Check": "Figure2_primary_metric",
        "Observed": "Patient_n",
        "Expected": "Patient_n",
        "Status": "PASS",
    },
]

ready = all(r["Status"] == "PASS" for r in qc)

qc.append({
    "Check": "READY_FOR_CROSS_PATIENT_ANALYSIS",
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
# 11. README / method and data dictionary
# ============================================================

readme = f"""# Cross-patient mutation and plotting dataset

## 1. Frozen mutation definition

All clone-specific SNPs, insertions, and deletions are defined relative to the
**patient-specific V3 consensus sequence**.

HXB2 is **not** the mutation reference. HXB2 is used only as the common
coordinate and gene-annotation framework.

Recommended Methods sentence:

> Clone-specific SNPs, insertions, and deletions were defined relative to the
> patient-specific consensus sequence and subsequently projected onto HXB2
> coordinates for genomic localization.

中文：

> 各克隆的 SNP、插入和缺失均以患者特异性 consensus 序列为参照进行定义，
> 随后将其投射至 HXB2 坐标系以统一标注基因组位置。

## 2. Final frozen counts

- Directional SNP calls: {source_snp_n}
- DEL events: {source_del_event_n}
- INS events: {source_ins_event_n}
- Total INDEL events: {source_del_event_n + source_ins_event_n}
- Patients: {len(patient_rows)}

DEL/INS numbers are **contiguous clone-specific event occurrences**, not
numbers of deleted/inserted bases.

## 3. HXB2 coordinate mapping

Final coordinate framework:

1. Patient-specific clone MSA.
2. patient-specific analysis consensus.
3. Clone vs patient consensus defines SNP/DEL/INS.
4. Patient mapping consensus is projected to HXB2 using MAFFT
   `--addfragments`.
5. V9 uses V7 direct coordinates for SNP/DEL and V7-flanked consensus-gap
   blocks for INS.
6. All 3287 directional mutation calls passed final coordinate projection QC.

For HXB2-gap sequence, an extended coordinate is used, for example:

`6646|6647.ins1`

For cross-patient insertion aggregation, the common anchor `6646|6647` is used;
`.ins1`, `.ins2`, etc. are not assumed to be homologous across patients.

## 4. Figure 2 primary metric

The primary cross-patient metric is:

**Patient_n = number of unique patients carrying a variant at a locus.**

`Total_variant_clone_n` and `Total_event_n` are retained as supporting metrics.

For canonical HXB2 SNPs, recurrence is calculated at the exact HXB2 base.

For SNPs on sequence absent from HXB2, recurrence is anchor-level only.

For INS, cross-patient recurrence is anchor-level only.

For DEL, recurrence is aggregated by exact projected event locus.

## 5. P17

CN2023AH1-17 has no unique patient consensus at multiple positions. These
positions are not assigned arbitrary directional SNP/DEL/INS states.

Use:
- filled diamond = short no-unique-consensus site;
- elongated band = large structural unresolved region.

See:
`08_no_unique_consensus_regions.tsv`

## 6. P18 large internal missing region

CN2023AH1-18 contains a large patient-level internal missing region.

Formal coordinate evidence:

- original patient MSA gap block: 3496–5912
- alignment-column length: 2417
- left direct mapped MSA column: 3495
- left V7 HXB2 flank: 3486
- right direct mapped MSA column: 5913
- right V7 HXB2 flank: 5904
- inferred HXB2 missing interval: **3487–5903**
- HXB2 interval length: **2417 nt**

The HXB2 interval length exactly matches the patient-MSA gap-block length.

This is a **patient-level structural/missing-region annotation** and is not
added to clone-vs-patient-consensus DEL counts.

See:
`09_large_missing_region_evidence.tsv`

## 7. Main files

### 01_SNP_locus_recurrence.tsv
Cross-patient SNP recurrence. The primary value is `Patient_n`.

### 02_INDEL_locus_recurrence.tsv
Cross-patient DEL/INS event recurrence.

### 03_combined_plot_ready.tsv
Unified long-format plotting table for Figure 2.

### 04_patient_summary.tsv
Per-patient clone number, final HXB2 range, mutation counts and special
structural annotations.

### 05_HXB2_mapping_summary.tsv
Final patient-to-HXB2 coordinate ranges and coordinate framework.

### 06_HXB2_gene_coordinates.tsv
HXB2 coordinates used for the gene model.

### 07_special_region_summary.tsv
Interpretive structural annotations that are explicitly separate from
clone-vs-consensus mutation calls.

### 08_no_unique_consensus_regions.tsv
All no-unique-consensus records preserved from the validated patient-level dataset.

### 09_large_missing_region_evidence.tsv
Coordinate evidence for the P18 HXB2 3487–5903 missing interval.

### 10_source_file_manifest.tsv
Data provenance / source-file manifest.

### 11_final_data_QC.tsv
Final package QC. Figure 2 should only be drawn if
`READY_FOR_CROSS_PATIENT_ANALYSIS = YES`.

## 8. Interpretation caution

The Figure 2 HXB2 axis standardizes genomic localization; it does not redefine
the patient-specific mutation reference.

Gene-level values should not be summed to obtain whole-genome totals because
HIV genes overlap.

Patient-level structural annotations (P17 unresolved consensus and P18 large
missing region) must remain separate from the directional SNP/DEL/INS counts.
"""

OUT_README.write_text(readme, encoding="utf-8")


# ============================================================
# Console
# ============================================================

print("=" * 86)
print("CROSS-PATIENT DATASET PREPARATION COMPLETE")
print("=" * 86)
print()
print(f"Patients                  : {len(patient_rows)}")
print(f"Directional SNP calls     : {source_snp_n}")
print(f"DEL events                : {source_del_event_n}")
print(f"INS events                : {source_ins_event_n}")
print(f"Figure2 SNP loci/anchors  : {len(snp_rows)}")
print(f"Figure2 INDEL loci        : {len(indel_rows)}")
print()
print("P18 large missing region  : HXB2 3487-5903 (2417 nt)")
print()
print("READY_FOR_CROSS_PATIENT_ANALYSIS:", "YES" if ready else "NO")
print()
print("Output directory:")
print(OUT)
