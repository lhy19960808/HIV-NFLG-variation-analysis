#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Project within-patient mutation calls onto the HXB2 coordinate system.

Biological definition
---------------------
Mutations are defined exclusively relative to the patient-specific
analysis consensus.

HXB2 is used only for genomic localization / coordinate annotation.

Coordinate rules
----------------
1. SNP and DEL:
   original patient MSA column
       -> mapping-consensus bridge
       -> MAFFT --addfragments HXB2 coordinate.

2. INS:
   patient-consensus gap columns are absent from the ungapped mapping
   consensus. Therefore the complete contiguous consensus-gap block is
   identified in the original patient MSA and anchored by the nearest
   directly HXB2-mapped columns on both sides.

   If the two flanks are adjacent canonical HXB2 bases, coordinates are
   represented as:

       6646|6647.ins1
       6646|6647.ins2
       ...

   If flanks are not adjacent, no false single-base precision is invented;
   an interval anchor is retained.

This implementation depends only on the publication-release upstream outputs.
"""

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


HXB2_LEN = 9719

GENES = [
    ("gag", 790, 2292),
    ("pol", 2085, 5096),
    ("vif", 5041, 5619),
    ("vpr", 5559, 5850),
    ("tat_exon1", 5831, 6045),
    ("rev_exon1", 5970, 6045),
    ("vpu", 6062, 6310),
    ("gp120", 6225, 7757),
    ("gp41", 7758, 8795),
    ("tat_exon2", 8379, 8469),
    ("rev_exon2", 8379, 8653),
    ("nef", 8797, 9417),
]


# ============================================================
# CLI
# ============================================================

parser = argparse.ArgumentParser(
    description=(
        "Project clone-vs-patient-consensus SNP/DEL/INS calls "
        "onto the final HXB2 coordinate framework."
    )
)

parser.add_argument(
    "--mutation-dir",
    required=True,
    type=Path,
    help="Clean Step 03 mutation-call output directory."
)

parser.add_argument(
    "--bridge-dir",
    required=True,
    type=Path,
    help="Clean Step 04A mapping-consensus / MSA-column bridge directory."
)

parser.add_argument(
    "--hxb2-map-dir",
    required=True,
    type=Path,
    help="Clean Step 04B final MAFFT --addfragments HXB2 mapping directory."
)

parser.add_argument(
    "--output-dir",
    required=True,
    type=Path,
    help="Output directory."
)

args = parser.parse_args()


CALL_FILE = (
    args.mutation_dir /
    "01_clone_mutation_calls.tsv"
)

PHASE2_EVENTS = (
    args.mutation_dir /
    "02_indel_events.tsv"
)

BRIDGE_FILE = (
    args.bridge_dir /
    "02_original_MSA_column_consensus_index.tsv"
)

HXB2_MAP_FILE = (
    args.hxb2_map_dir /
    "01_mapping_consensus_to_HXB2.tsv"
)

HXB2_RANGE_FILE = (
    args.hxb2_map_dir /
    "02_patient_HXB2_range_summary.tsv"
)

OUTDIR = args.output_dir
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_CALLS = OUTDIR / "01_final_clone_mutation_calls.tsv"
OUT_EVENTS = OUTDIR / "02_final_indel_events.tsv"
OUT_QC = OUTDIR / "03_projection_QC.tsv"
OUT_BLOCK_QC = OUTDIR / "04_insertion_block_QC.tsv"
OUT_REVIEW = OUTDIR / "05_review_or_unmapped_calls.tsv"
OUT_RANGES = OUTDIR / "06_final_patient_HXB2_ranges.tsv"
OUT_UNRESOLVED = OUTDIR / "07_no_unique_consensus_regions.tsv"


# ============================================================
# Helpers
# ============================================================

def read_tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, fields, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for row in rows:
            writer.writerow({
                key: row.get(key, "")
                for key in fields
            })


def require(path):
    if not path.exists():
        raise FileNotFoundError(
            f"Required input not found: {path}"
        )


def as_int(value):
    return int(float(str(value).strip()))


def patient_sort_key(patient):
    m = re.fullmatch(
        r"CN(\d{4})AH(\d+)-(\d+)",
        str(patient)
    )

    if m:
        return tuple(map(int, m.groups()))

    return (9999, 9999, str(patient))


def clone_sort_key(name):
    m = re.search(
        r"\.clone(\d+)",
        str(name),
        re.I
    )

    return int(m.group(1)) if m else 10**9


def detect_column(rows, candidates, label):
    if not rows:
        raise RuntimeError(
            f"{label}: empty table"
        )

    available = list(rows[0].keys())

    for candidate in candidates:
        if candidate in available:
            return candidate

    raise RuntimeError(
        f"Cannot detect {label}. "
        f"Available columns: {available}"
    )


def parse_hxb2_coord(coord):
    s = str(coord).strip()

    m = re.fullmatch(r"(\d+)", s)

    if m:
        p = int(m.group(1))

        return {
            "kind": "HXB2_base",
            "left": p,
            "right": p,
            "plot_x": float(p),
            "anchor": s,
        }

    m = re.fullmatch(
        r"(\d+)\|(\d+)\.ins(\d+)",
        s
    )

    if m:
        left = int(m.group(1))
        right = int(m.group(2))

        return {
            "kind": "HXB2_insertion_segment",
            "left": left,
            "right": right,
            "plot_x": (left + right) / 2.0,
            "anchor": f"{left}|{right}",
        }

    m = re.fullmatch(
        r"(\d+)\.\.(\d+)",
        s
    )

    if m:
        left = int(m.group(1))
        right = int(m.group(2))

        return {
            "kind": "HXB2_interval",
            "left": left,
            "right": right,
            "plot_x": (left + right) / 2.0,
            "anchor": s,
        }

    raise ValueError(
        f"Unrecognized HXB2 coordinate: {coord!r}"
    )


def gene_names_for_interval(left, right):
    lo = max(1, min(left, right))
    hi = min(HXB2_LEN, max(left, right))

    hits = []

    for gene, start, end in GENES:
        if not (hi < start or lo > end):
            hits.append(gene)

    return (
        ",".join(hits)
        if hits
        else "intergenic_or_unassigned"
    )


def gene_names_for_coord(coord):
    try:
        parsed = parse_hxb2_coord(coord)

        return gene_names_for_interval(
            parsed["left"],
            parsed["right"]
        )

    except Exception:
        return "intergenic_or_unassigned"


# ============================================================
# Inputs
# ============================================================

for path in [
    CALL_FILE,
    PHASE2_EVENTS,
    BRIDGE_FILE,
    HXB2_MAP_FILE,
    HXB2_RANGE_FILE,
]:
    require(path)


calls = read_tsv(CALL_FILE)
phase2_events = read_tsv(PHASE2_EVENTS)
bridge = read_tsv(BRIDGE_FILE)
hxb2_map = read_tsv(HXB2_MAP_FILE)
hxb2_ranges = read_tsv(HXB2_RANGE_FILE)


PATIENT_COL = detect_column(
    calls,
    ["Patient", "patient"],
    "patient column"
)

CLONE_COL = detect_column(
    calls,
    ["Clone", "clone"],
    "clone column"
)

MSA_COL = detect_column(
    calls,
    [
        "Original_patient_MSA_column",
        "MSA_col",
        "MSA_column",
        "Patient_MSA_column",
    ],
    "patient MSA column"
)

MUT_COL = detect_column(
    calls,
    ["Mutation_type", "Mutation"],
    "mutation column"
)


# ============================================================
# V6 bridge + direct V7 lookup
# ============================================================

bridge_by_patient = defaultdict(dict)

for row in bridge:
    patient = row["Patient"]
    col = as_int(
        row["Original_patient_MSA_column"]
    )

    bridge_by_patient[patient][col] = row


hxb2_by_mapping_pos = {}

for row in hxb2_map:
    key = (
        row["Patient"],
        as_int(
            row["Mapping_consensus_position"]
        )
    )

    hxb2_by_mapping_pos[key] = row


direct_hxb2_by_msa = {}
direct_cols_by_patient = defaultdict(list)

for patient, colmap in bridge_by_patient.items():

    for msa_col, bridge_row in colmap.items():

        mapping_pos = str(
            bridge_row[
                "Mapping_consensus_position"
            ]
        ).strip()

        if not mapping_pos:
            continue

        key = (
            patient,
            as_int(mapping_pos)
        )

        if key not in hxb2_by_mapping_pos:
            raise RuntimeError(
                f"Missing V7 mapping for "
                f"{patient}, mapping position "
                f"{mapping_pos}"
            )

        direct_hxb2_by_msa[
            (patient, msa_col)
        ] = hxb2_by_mapping_pos[key]

        direct_cols_by_patient[
            patient
        ].append(msa_col)


for patient in direct_cols_by_patient:
    direct_cols_by_patient[
        patient
    ].sort()


# ============================================================
# Build contiguous patient-consensus gap blocks
# ============================================================

gap_block_by_msa = {}
gap_block_rows = []

for patient, colmap in bridge_by_patient.items():

    gap_cols = sorted(
        col
        for col, row in colmap.items()
        if (
            str(
                row["Analysis_consensus_state"]
            ).strip() == "-"
            or
            str(
                row[
                    "Mapping_inclusion_class"
                ]
            ).strip()
            == "Analysis_consensus_gap_excluded"
        )
    )

    if not gap_cols:
        continue

    runs = []
    current = []

    for col in gap_cols:

        if not current:
            current = [col]

        elif col == current[-1] + 1:
            current.append(col)

        else:
            runs.append(current)
            current = [col]

    if current:
        runs.append(current)

    direct_cols = (
        direct_cols_by_patient.get(
            patient,
            []
        )
    )

    for block_index, run in enumerate(
        runs,
        start=1
    ):

        start = run[0]
        end = run[-1]

        left_candidates = [
            c for c in direct_cols
            if c < start
        ]

        right_candidates = [
            c for c in direct_cols
            if c > end
        ]

        left_col = (
            left_candidates[-1]
            if left_candidates
            else None
        )

        right_col = (
            right_candidates[0]
            if right_candidates
            else None
        )

        left_coord = (
            direct_hxb2_by_msa[
                (patient, left_col)
            ]["HXB2_coordinate"]
            if left_col is not None
            else ""
        )

        right_coord = (
            direct_hxb2_by_msa[
                (patient, right_col)
            ]["HXB2_coordinate"]
            if right_col is not None
            else ""
        )

        block_id = (
            f"{patient}.gapblock"
            f"{block_index:04d}"
        )

        status = "REVIEW"
        precision = "UNRESOLVED"
        anchor = ""
        plot_x = ""

        try:
            left_parsed = (
                parse_hxb2_coord(
                    left_coord
                )
            )

            right_parsed = (
                parse_hxb2_coord(
                    right_coord
                )
            )

            if (
                left_parsed["kind"]
                == "HXB2_base"
                and
                right_parsed["kind"]
                == "HXB2_base"
                and
                right_parsed["left"]
                == left_parsed["left"] + 1
            ):
                left_base = (
                    left_parsed["left"]
                )

                right_base = (
                    right_parsed["left"]
                )

                anchor = (
                    f"{left_base}|"
                    f"{right_base}"
                )

                plot_x = (
                    left_base + right_base
                ) / 2.0

                precision = (
                    "EXACT_ADJACENT_HXB2_ANCHOR"
                )

                status = "PASS"

            else:
                anchor = (
                    f"{left_coord}.."
                    f"{right_coord}"
                )

                plot_x = (
                    left_parsed["plot_x"]
                    + right_parsed["plot_x"]
                ) / 2.0

                precision = (
                    "INTERVAL_BETWEEN_HXB2_FLANKS"
                )

                status = "PASS_INTERVAL"

        except Exception:
            pass

        block_info = {
            "Patient": patient,
            "Gap_block_ID": block_id,
            "Gap_block_start_MSA_column":
                start,
            "Gap_block_end_MSA_column":
                end,
            "Gap_block_length_alignment_columns":
                len(run),
            "Left_flank_MSA_column":
                left_col if left_col is not None else "",
            "Left_flank_HXB2_coordinate":
                left_coord,
            "Right_flank_MSA_column":
                right_col if right_col is not None else "",
            "Right_flank_HXB2_coordinate":
                right_coord,
            "Insertion_HXB2_anchor":
                anchor,
            "Insertion_plot_coordinate":
                (
                    f"{plot_x:.3f}"
                    if isinstance(
                        plot_x,
                        (float, int)
                    )
                    else ""
                ),
            "Coordinate_precision":
                precision,
            "Block_projection_status":
                status,
        }

        gap_block_rows.append(
            block_info
        )

        for rank, msa_col in enumerate(
            run,
            start=1
        ):
            gap_block_by_msa[
                (patient, msa_col)
            ] = {
                **block_info,
                "Insertion_alignment_rank_in_gap_block":
                    rank,
            }


# ============================================================
# Project mutation calls
# ============================================================

final_calls = []

for row in calls:

    patient = str(
        row[PATIENT_COL]
    ).strip()

    clone = str(
        row[CLONE_COL]
    ).strip()

    msa_col = as_int(
        row[MSA_COL]
    )

    mut = str(
        row[MUT_COL]
    ).strip().upper()

    if mut not in {
        "SNP",
        "DEL",
        "INS"
    }:
        continue

    bridge_row = (
        bridge_by_patient[
            patient
        ][msa_col]
    )

    mapping_pos = str(
        bridge_row[
            "Mapping_consensus_position"
        ]
    ).strip()

    consensus_allele = str(
        row.get(
            "Analysis_consensus_state",
            ""
        )
    ).strip()

    clone_allele = str(
        row.get(
            "Clone_interpreted_state",
            ""
        )
    ).strip()

    snp_type = (
        str(
            row.get(
                "Change",
                ""
            )
        ).strip()
        if mut == "SNP"
        else ""
    )

    out = {
        "Patient":
            patient,
        "Clone":
            clone,
        "Original_patient_MSA_column":
            msa_col,
        "Mutation_type":
            mut,
        "SNP_type":
            snp_type,
        "Consensus_allele":
            consensus_allele,
        "Clone_allele":
            clone_allele,
        "Final_HXB2_coordinate":
            "",
        "Final_coordinate_class":
            "",
        "Final_HXB2_anchor":
            "",
        "Final_plot_coordinate":
            "",
        "Final_genes":
            "",
        "Coordinate_source":
            "",
        "Projection_status":
            "",
        "Analysis_consensus_state_at_MSA_column":
            bridge_row[
                "Analysis_consensus_state"
            ],
        "Analysis_consensus_reason":
            bridge_row[
                "Consensus_reason"
            ],
        "Mapping_inclusion_class":
            bridge_row[
                "Mapping_inclusion_class"
            ],
        "Gap_block_ID":
            "",
        "Gap_block_start_MSA_column":
            "",
        "Gap_block_end_MSA_column":
            "",
        "Gap_block_length_alignment_columns":
            "",
        "Insertion_alignment_rank_in_gap_block":
            "",
        "Left_flank_MSA_column":
            "",
        "Left_flank_HXB2_coordinate":
            "",
        "Right_flank_MSA_column":
            "",
        "Right_flank_HXB2_coordinate":
            "",
        "Coordinate_precision":
            "",
    }

    # --------------------------------------------------------
    # SNP / DEL -> direct V7
    # --------------------------------------------------------

    if mut in {"SNP", "DEL"}:

        if not mapping_pos:
            out[
                "Projection_status"
            ] = (
                "REVIEW_NO_DIRECT_HXB2_MAPPING"
            )

        else:
            hxb2row = (
                direct_hxb2_by_msa[
                    (patient, msa_col)
                ]
            )

            coord = hxb2row[
                "HXB2_coordinate"
            ]

            out[
                "Final_HXB2_coordinate"
            ] = coord

            out[
                "Final_coordinate_class"
            ] = hxb2row[
                "Coordinate_class"
            ]

            out[
                "Final_HXB2_anchor"
            ] = hxb2row[
                "HXB2_anchor"
            ]

            out[
                "Final_plot_coordinate"
            ] = hxb2row[
                "Plot_coordinate"
            ]

            out[
                "Final_genes"
            ] = gene_names_for_coord(
                coord
            )

            out[
                "Coordinate_source"
            ] = "DIRECT_HXB2"

            out[
                "Coordinate_precision"
            ] = "DIRECT_HXB2"

            out[
                "Projection_status"
            ] = "PASS_DIRECT_HXB2"

    # --------------------------------------------------------
    # INS -> V7-flanked consensus gap block
    # --------------------------------------------------------

    else:

        info = gap_block_by_msa.get(
            (patient, msa_col)
        )

        if info is None:

            out[
                "Projection_status"
            ] = (
                "REVIEW_INS_NOT_IN_"
                "CONSENSUS_GAP_BLOCK"
            )

            out[
                "Coordinate_source"
            ] = (
                "HXB2_GAP_BLOCK"
            )

        else:

            for key in [
                "Gap_block_ID",
                "Gap_block_start_MSA_column",
                "Gap_block_end_MSA_column",
                "Gap_block_length_alignment_columns",
                "Insertion_alignment_rank_in_gap_block",
                "Left_flank_MSA_column",
                "Left_flank_HXB2_coordinate",
                "Right_flank_MSA_column",
                "Right_flank_HXB2_coordinate",
                "Coordinate_precision",
            ]:
                out[key] = info[key]

            anchor = info[
                "Insertion_HXB2_anchor"
            ]

            rank = info[
                "Insertion_alignment_rank_in_gap_block"
            ]

            out[
                "Final_HXB2_anchor"
            ] = anchor

            out[
                "Final_plot_coordinate"
            ] = info[
                "Insertion_plot_coordinate"
            ]

            out[
                "Coordinate_source"
            ] = "HXB2_GAP_BLOCK"

            if (
                info[
                    "Block_projection_status"
                ]
                == "PASS"
            ):

                out[
                    "Final_HXB2_coordinate"
                ] = (
                    f"{anchor}.ins{rank}"
                )

                out[
                    "Final_coordinate_class"
                ] = (
                    "HXB2_insertion_segment"
                )

                out[
                    "Projection_status"
                ] = (
                    "PASS_HXB2_INSERTION_BLOCK"
                )

            elif (
                info[
                    "Block_projection_status"
                ]
                == "PASS_INTERVAL"
            ):

                out[
                    "Final_HXB2_coordinate"
                ] = anchor

                out[
                    "Final_coordinate_class"
                ] = (
                    "HXB2_interval_"
                    "insertion_anchor"
                )

                out[
                    "Projection_status"
                ] = (
                    "PASS_HXB2_INTERVAL_"
                    "INSERTION_BLOCK"
                )

            else:

                out[
                    "Projection_status"
                ] = (
                    "REVIEW_UNRESOLVED_"
                    "HXB2_INSERTION_BLOCK"
                )

            try:

                left = parse_hxb2_coord(
                    info[
                        "Left_flank_HXB2_coordinate"
                    ]
                )

                right = parse_hxb2_coord(
                    info[
                        "Right_flank_HXB2_coordinate"
                    ]
                )

                out[
                    "Final_genes"
                ] = gene_names_for_interval(
                    left["left"],
                    right["right"]
                )

            except Exception:

                out[
                    "Final_genes"
                ] = (
                    "intergenic_or_unassigned"
                )

    final_calls.append(out)


final_calls.sort(
    key=lambda row: (
        patient_sort_key(
            row["Patient"]
        ),
        clone_sort_key(
            row["Clone"]
        ),
        int(
            row[
                "Original_patient_MSA_column"
            ]
        ),
        row[
            "Mutation_type"
        ],
    )
)


CALL_FIELDS = [
    "Patient",
    "Clone",
    "Original_patient_MSA_column",
    "Mutation_type",
    "SNP_type",
    "Consensus_allele",
    "Clone_allele",
    "Final_HXB2_coordinate",
    "Final_coordinate_class",
    "Final_HXB2_anchor",
    "Final_plot_coordinate",
    "Final_genes",
    "Coordinate_source",
    "Projection_status",
    "Analysis_consensus_state_at_MSA_column",
    "Analysis_consensus_reason",
    "Mapping_inclusion_class",
    "Gap_block_ID",
    "Gap_block_start_MSA_column",
    "Gap_block_end_MSA_column",
    "Gap_block_length_alignment_columns",
    "Insertion_alignment_rank_in_gap_block",
    "Left_flank_MSA_column",
    "Left_flank_HXB2_coordinate",
    "Right_flank_MSA_column",
    "Right_flank_HXB2_coordinate",
    "Coordinate_precision",
]

write_tsv(
    OUT_CALLS,
    CALL_FIELDS,
    final_calls
)


# ============================================================
# Reconstruct contiguous clone-specific INDEL events
# ============================================================

indel_calls = [
    row
    for row in final_calls
    if row["Mutation_type"]
    in {"DEL", "INS"}
]

groups = defaultdict(list)

for row in indel_calls:

    groups[
        (
            row["Patient"],
            row["Clone"],
            row["Mutation_type"],
        )
    ].append(row)


event_rows = []
event_counter = 0

for (
    patient,
    clone,
    mutation_type
), rows in groups.items():

    rows = sorted(
        rows,
        key=lambda row:
            int(
                row[
                    "Original_patient_MSA_column"
                ]
            )
    )

    runs = []
    current = []

    for row in rows:

        col = int(
            row[
                "Original_patient_MSA_column"
            ]
        )

        if not current:
            current = [row]

        else:

            previous_col = int(
                current[-1][
                    "Original_patient_MSA_column"
                ]
            )

            if col == previous_col + 1:
                current.append(row)

            else:
                runs.append(current)
                current = [row]

    if current:
        runs.append(current)

    for run in runs:

        event_counter += 1

        coords = [
            row[
                "Final_HXB2_coordinate"
            ]
            for row in run
            if row[
                "Final_HXB2_coordinate"
            ]
        ]

        anchors = [
            row[
                "Final_HXB2_anchor"
            ]
            for row in run
            if row[
                "Final_HXB2_anchor"
            ]
        ]

        xs = [
            float(
                row[
                    "Final_plot_coordinate"
                ]
            )
            for row in run
            if row[
                "Final_plot_coordinate"
            ] not in ("", None)
        ]

        statuses = sorted(
            set(
                row[
                    "Projection_status"
                ]
                for row in run
            )
        )

        genes = sorted(
            set(
                gene
                for row in run
                for gene in row[
                    "Final_genes"
                ].split(",")
                if gene
            )
        )

        if mutation_type == "INS":

            unique_anchors = sorted(
                set(anchors)
            )

            if len(unique_anchors) == 1:
                locus = unique_anchors[0]

            else:
                locus = "..".join(
                    unique_anchors
                )

        else:

            if coords:

                locus = (
                    coords[0]
                    if coords[0] == coords[-1]
                    else (
                        f"{coords[0]}.."
                        f"{coords[-1]}"
                    )
                )

            else:
                locus = ""

        plot_x = (
            (min(xs) + max(xs)) / 2.0
            if xs
            else None
        )

        status = (
            "PASS"
            if all(
                s.startswith("PASS_")
                for s in statuses
            )
            else "REVIEW"
        )

        event_rows.append({
            "Patient":
                patient,
            "Clone":
                clone,
            "Event_ID":
                (
                    f"{patient}."
                    f"{mutation_type}."
                    f"event{event_counter:04d}"
                ),
            "Event_type":
                mutation_type,
            "Start_original_MSA_column":
                run[0][
                    "Original_patient_MSA_column"
                ],
            "End_original_MSA_column":
                run[-1][
                    "Original_patient_MSA_column"
                ],
            "Event_length_bp":
                len(run),
            "Start_HXB2_coordinate":
                coords[0] if coords else "",
            "End_HXB2_coordinate":
                coords[-1] if coords else "",
            "Event_HXB2_locus":
                locus,
            "Plot_coordinate":
                (
                    f"{plot_x:.3f}"
                    if plot_x is not None
                    else ""
                ),
            "Genes":
                ",".join(genes),
            "Projection_status":
                status,
            "Constituent_projection_statuses":
                ",".join(statuses),
        })


event_rows.sort(
    key=lambda row: (
        patient_sort_key(
            row["Patient"]
        ),
        clone_sort_key(
            row["Clone"]
        ),
        int(
            row[
                "Start_original_MSA_column"
            ]
        ),
        row[
            "Event_type"
        ],
    )
)


EVENT_FIELDS = [
    "Patient",
    "Clone",
    "Event_ID",
    "Event_type",
    "Start_original_MSA_column",
    "End_original_MSA_column",
    "Event_length_bp",
    "Start_HXB2_coordinate",
    "End_HXB2_coordinate",
    "Event_HXB2_locus",
    "Plot_coordinate",
    "Genes",
    "Projection_status",
    "Constituent_projection_statuses",
]

write_tsv(
    OUT_EVENTS,
    EVENT_FIELDS,
    event_rows
)


# ============================================================
# Final patient HXB2 ranges
# ============================================================

range_rows = []

for row in sorted(
    hxb2_ranges,
    key=lambda row:
        patient_sort_key(
            row["Patient"]
        )
):

    range_rows.append({
        "Patient":
            row["Patient"],
        "Final_HXB2_start":
            row["HXB2_start"],
        "Final_HXB2_end":
            row["HXB2_end"],
        "Final_HXB2_span":
            row["HXB2_span"],
        "Coordinate_framework":
            (
                "MAFFT_addfragments_"
                "plus_HXB2_flanked_"
                "consensus_gap_blocks"
            ),
    })


write_tsv(
    OUT_RANGES,
    [
        "Patient",
        "Final_HXB2_start",
        "Final_HXB2_end",
        "Final_HXB2_span",
        "Coordinate_framework",
    ],
    range_rows
)


# ============================================================
# No-unique-consensus regions
# ============================================================

unresolved_rows = []

for patient, colmap in bridge_by_patient.items():

    cols = sorted(
        col
        for col, row in colmap.items()
        if row[
            "Consensus_reason"
        ] == "No_unique_consensus"
    )

    if not cols:
        continue

    runs = []
    current = []

    for col in cols:

        if not current:
            current = [col]
            continue

        same_class = (
            colmap[col][
                "Mapping_inclusion_class"
            ]
            ==
            colmap[
                current[-1]
            ][
                "Mapping_inclusion_class"
            ]
        )

        if (
            col == current[-1] + 1
            and same_class
        ):
            current.append(col)

        else:
            runs.append(current)
            current = [col]

    if current:
        runs.append(current)

    direct_cols = (
        direct_cols_by_patient.get(
            patient,
            []
        )
    )

    for run in runs:

        start = run[0]
        end = run[-1]

        left_candidates = [
            c
            for c in direct_cols
            if c < start
        ]

        right_candidates = [
            c
            for c in direct_cols
            if c > end
        ]

        left_col = (
            left_candidates[-1]
            if left_candidates
            else None
        )

        right_col = (
            right_candidates[0]
            if right_candidates
            else None
        )

        left_coord = (
            direct_hxb2_by_msa[
                (patient, left_col)
            ][
                "HXB2_coordinate"
            ]
            if left_col is not None
            else ""
        )

        right_coord = (
            direct_hxb2_by_msa[
                (patient, right_col)
            ][
                "HXB2_coordinate"
            ]
            if right_col is not None
            else ""
        )

        unresolved_rows.append({
            "Patient":
                patient,
            "No_unique_class":
                colmap[
                    start
                ][
                    "Mapping_inclusion_class"
                ],
            "Start_original_MSA_column":
                start,
            "End_original_MSA_column":
                end,
            "Length_alignment_columns":
                len(run),
            "Left_anchor_HXB2_coordinate":
                left_coord,
            "Right_anchor_HXB2_coordinate":
                right_coord,
            "Display_interval":
                (
                    f"{left_coord}.."
                    f"{right_coord}"
                    if left_coord
                    and right_coord
                    else ""
                ),
        })


write_tsv(
    OUT_UNRESOLVED,
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
    unresolved_rows
)


# ============================================================
# QC
# ============================================================

mutation_counts = Counter(
    row["Mutation_type"]
    for row in final_calls
)

event_counts = Counter(
    row["Event_type"]
    for row in event_rows
)

review_calls = [
    row
    for row in final_calls
    if not row[
        "Projection_status"
    ].startswith("PASS_")
]

exact_ins_calls = sum(
    1
    for row in final_calls
    if (
        row["Mutation_type"] == "INS"
        and
        row[
            "Projection_status"
        ]
        == "PASS_HXB2_INSERTION_BLOCK"
    )
)

interval_ins_calls = sum(
    1
    for row in final_calls
    if (
        row["Mutation_type"] == "INS"
        and
        row[
            "Projection_status"
        ]
        == (
            "PASS_HXB2_INTERVAL_"
            "INSERTION_BLOCK"
        )
    )
)

qc = [
    (
        "Directional_mutation_calls",
        len(final_calls),
        3287
    ),
    (
        "SNP_call_n",
        mutation_counts["SNP"],
        2477
    ),
    (
        "DEL_base_call_n",
        mutation_counts["DEL"],
        592
    ),
    (
        "INS_base_call_n",
        mutation_counts["INS"],
        218
    ),
    (
        "DEL_event_n",
        event_counts["DEL"],
        44
    ),
    (
        "INS_event_n",
        event_counts["INS"],
        31
    ),
    (
        "Total_indel_event_n",
        len(event_rows),
        75
    ),
    (
        "NonPASS_call_n",
        len(review_calls),
        0
    ),
    (
        "Exact_anchor_INS_call_n",
        exact_ins_calls,
        218
    ),
    (
        "Interval_anchor_INS_call_n",
        interval_ins_calls,
        0
    ),
    (
        "Patient_range_n",
        len(range_rows),
        11
    ),
]


qc_rows = []

for check, observed, expected in qc:

    qc_rows.append({
        "Check":
            check,
        "Observed":
            observed,
        "Expected":
            expected,
        "Status":
            (
                "PASS"
                if observed == expected
                else "FAIL"
            ),
    })


ready = all(
    row["Status"] == "PASS"
    for row in qc_rows
)

qc_rows.append({
    "Check":
        "READY_FOR_DOWNSTREAM_ANALYSIS",
    "Observed":
        "YES" if ready else "NO",
    "Expected":
        "YES",
    "Status":
        "PASS" if ready else "FAIL",
})


write_tsv(
    OUT_QC,
    [
        "Check",
        "Observed",
        "Expected",
        "Status"
    ],
    qc_rows
)

write_tsv(
    OUT_BLOCK_QC,
    [
        "Patient",
        "Gap_block_ID",
        "Gap_block_start_MSA_column",
        "Gap_block_end_MSA_column",
        "Gap_block_length_alignment_columns",
        "Left_flank_MSA_column",
        "Left_flank_HXB2_coordinate",
        "Right_flank_MSA_column",
        "Right_flank_HXB2_coordinate",
        "Insertion_HXB2_anchor",
        "Insertion_plot_coordinate",
        "Coordinate_precision",
        "Block_projection_status",
    ],
    gap_block_rows
)

write_tsv(
    OUT_REVIEW,
    CALL_FIELDS,
    review_calls
)


print("=" * 78)
print("CLEAN FINAL MUTATION COORDINATE PROJECTION COMPLETE")
print("=" * 78)
print()
print("Mutation reference: patient-specific analysis consensus")
print("HXB2 role: coordinate / gene annotation only")
print()
print(f"Directional mutation calls : {len(final_calls)}")
print(f"SNP calls                  : {mutation_counts['SNP']}")
print(f"DEL base calls             : {mutation_counts['DEL']}")
print(f"INS base calls             : {mutation_counts['INS']}")
print(f"DEL events                 : {event_counts['DEL']}")
print(f"INS events                 : {event_counts['INS']}")
print(f"Exact-anchor INS calls     : {exact_ins_calls}")
print(f"Interval-anchor INS calls  : {interval_ins_calls}")
print(f"Non-PASS calls             : {len(review_calls)}")
print()
print(
    "READY_FOR_DOWNSTREAM_ANALYSIS:",
    "YES" if ready else "NO"
)

if not ready:
    raise SystemExit(1)
