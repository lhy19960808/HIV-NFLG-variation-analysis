#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Patient-specific HXB2 coordinate mapping with MAFFT --addfragments
=======================================================================

Goal
----
Use each patient's Step04A mapping consensus as a FRAGMENT and
add it independently to the HXB2 reference sequence.

This is COORDINATE VALIDATION ONLY.

It does NOT:
- rebuild the within-patient clone alignment;
- redefine the patient-specific consensus rule;
- re-call SNP / DEL / INS;
- create one global all-patient MSA.

Locked logic
------------
Mutation definition:
    clone vs patient-specific analysis consensus

Coordinate reference:
    HXB2

Workflow:
    Step04A patient mapping consensus
            +
          HXB2
            |
            | MAFFT --addfragments
            | NO --keeplength
            v
    patient-specific HXB2 coordinate map

Why no --keeplength?
--------------------
Patient insertions relative to HXB2 must be retained. Therefore HXB2 is allowed
to contain alignment gaps, which are represented as extended HXB2 coordinates,
for example:

    6646|6647.ins1
    6646|6647.ins2

HXB2 nucleotide numbering itself never changes. Only alignment-column numbers
change when gaps are inserted into the HXB2 row.

Inputs
------
21_pairwise_consensus_HXB2_v6/per_patient/*.mapping_consensus.fasta
00_raw/all_patients_NFLG_with_HXB2_original.fasta

Optional comparison inputs
--------------------------
21_pairwise_consensus_HXB2_v6/01_pairwise_consensus_position_to_HXB2.tsv
21_pairwise_consensus_HXB2_v6/03_patient_pairwise_HXB2_range_summary.tsv
20_patient_hxb2_coordinate_map_v5/03_patient_HXB2_range_summary_v5.tsv
18_hxb2_mapping_qc/04_clone_HXB2_mapping_summary.tsv

Outputs
-------
22_addfragments_HXB2_v7/
  01_mapping_consensus_to_HXB2.tsv
  02_patient_HXB2_range_summary.tsv
  07_addfragments_alignment_QC.tsv
  per_patient/
      <patient>.HXB2_reference.fasta
      <patient>.addfragments.input_mapping_consensus.fasta
      <patient>.HXB2_plus_mapping_consensus.addfragments.fasta
"""

from pathlib import Path
from collections import defaultdict
import argparse
import csv
import math
import re
import shutil
import subprocess


# ============================================================
# Paths
# ============================================================

parser = argparse.ArgumentParser(
    description="Map patient mapping consensuses independently to HXB2 using MAFFT --addfragments."
)

parser.add_argument(
    "--bridge-dir",
    required=True,
    type=Path,
    help="Step 04A directory containing mapping consensuses and bridge tables."
)

parser.add_argument(
    "--hxb2-source",
    required=True,
    type=Path,
    help="FASTA containing the HXB2 reference sequence."
)

parser.add_argument(
    "--output-dir",
    required=True,
    type=Path,
    help="Output directory for final patient-specific HXB2 coordinate maps."
)



args = parser.parse_args()

BRIDGE_DIR = args.bridge_dir
BRIDGE_PER_PATIENT = BRIDGE_DIR / "per_patient"

HXB2_SOURCE = args.hxb2_source




OUTDIR = args.output_dir
PER_PATIENT = OUTDIR / "per_patient"

OUTDIR.mkdir(parents=True, exist_ok=True)
PER_PATIENT.mkdir(parents=True, exist_ok=True)

MAP_OUT = OUTDIR / "01_mapping_consensus_to_HXB2.tsv"
RANGE_OUT = OUTDIR / "02_patient_HXB2_range_summary.tsv"
QC_OUT = OUTDIR / "07_addfragments_alignment_QC.tsv"

HXB2_LEN = 9719
BASES = set("ACGT")
QUERY_BASES = set("ACGTN")


# ============================================================
# General helpers
# ============================================================

def read_fasta(path):
    seqs = {}
    name = None
    buf = []

    with open(path, encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()

                name = line[1:].strip()
                buf = []
            else:
                buf.append(line)

        if name is not None:
            seqs[name] = "".join(buf).upper()

    return seqs


def write_fasta(path, seqs, width=80):
    with open(path, "w", encoding="utf-8") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n")

            for i in range(0, len(seq), width):
                f.write(seq[i:i+width] + "\n")


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

        for row in rows:
            w.writerow({
                k: row.get(k, "")
                for k in fields
            })


def patient_from_mapping_fasta(path):
    suffix = ".mapping_consensus.fasta"

    if not path.name.endswith(suffix):
        raise ValueError(path.name)

    return path.name[:-len(suffix)]


def patient_sort_key(patient):
    m = re.fullmatch(
        r"CN(\d{4})AH(\d+)-(\d+)",
        patient,
    )

    if m:
        return (
            int(m.group(1)),
            int(m.group(2)),
            int(m.group(3)),
        )

    return (9999, 9999, patient)


def median(values):
    values = sorted(values)

    if not values:
        return None

    n = len(values)

    if n % 2:
        return values[n // 2]

    return (
        values[n // 2 - 1]
        + values[n // 2]
    ) / 2.0


def exact_coord_to_numeric(coord):
    """
    Numeric helper for coordinate comparisons.

    HXB2 base:
        6646 -> 6646.0

    Extended HXB2:
        6646|6647.ins3 -> 6646.5

    Edge insertions:
        before1.insN -> 0.5
        after9719.insN -> 9719.5
    """

    s = str(coord).strip()

    if re.fullmatch(r"\d+", s):
        return float(int(s))

    m = re.fullmatch(
        r"(\d+)\|(\d+)\.ins(\d+)",
        s,
    )

    if m:
        left = int(m.group(1))
        right = int(m.group(2))

        return (left + right) / 2.0

    if re.fullmatch(r"before1\.ins\d+", s):
        return 0.5

    if re.fullmatch(
        rf"after{HXB2_LEN}\.ins\d+",
        s,
    ):
        return HXB2_LEN + 0.5

    return None


# ============================================================
# HXB2 extraction
# ============================================================

def extract_hxb2():
    seqs = read_fasta(HXB2_SOURCE)

    hits = [
        (name, seq)
        for name, seq in seqs.items()
        if "HXB2" in name.upper()
    ]

    if len(hits) != 1:
        raise RuntimeError(
            "Expected exactly one HXB2 sequence; "
            f"found {[x[0] for x in hits]}"
        )

    name, seq = hits[0]
    seq = seq.replace("-", "")

    if len(seq) != HXB2_LEN:
        raise RuntimeError(
            f"HXB2 length={len(seq)}, expected {HXB2_LEN}"
        )

    return name, seq


# ============================================================
# MAFFT --addfragments
# ============================================================

def run_addfragments(
    hxb2_ref_fasta,
    mapping_fasta,
    output_fasta,
):
    """
    Important:
      Do NOT add --keeplength.
      Query-specific insertions must be retained.
    """

    cmd = [
        "mafft",
        "--addfragments",
        str(mapping_fasta),
        "--thread", "8",
        "--inputorder",
        str(hxb2_ref_fasta),
    ]

    with open(output_fasta, "w", encoding="utf-8") as out:
        subprocess.run(
            cmd,
            stdout=out,
            check=True,
        )


# ============================================================
# HXB2 coordinate labels
# ============================================================

def build_hxb2_labels(hxb2_aligned):
    """
    Label EVERY alignment column.

    If HXB2 has a nucleotide:
        coordinate = 1..9719

    If HXB2 has a gap:
        coordinate = left|right.insN

    Example:
        HXB2 ...A--T...
        query...ACGT...

        A -> 6646
        C -> 6646|6647.ins1
        G -> 6646|6647.ins2
        T -> 6647
    """

    n = len(hxb2_aligned)

    base_pos = [None] * n
    hpos = 0

    for i, c in enumerate(hxb2_aligned):
        if c in BASES:
            hpos += 1
            base_pos[i] = hpos

    if hpos != HXB2_LEN:
        raise RuntimeError(
            f"Aligned HXB2 contains {hpos} bases; expected {HXB2_LEN}"
        )

    labels = [None] * n
    classes = [None] * n
    anchors = [None] * n
    plot_x = [None] * n

    i = 0

    while i < n:

        if hxb2_aligned[i] in BASES:
            pos = base_pos[i]

            labels[i] = str(pos)
            classes[i] = "HXB2_base"
            anchors[i] = str(pos)
            plot_x[i] = float(pos)

            i += 1
            continue

        # HXB2 gap run.
        j = i

        while (
            j < n
            and hxb2_aligned[j] == "-"
        ):
            j += 1

        # nearest HXB2 base on the left
        left = 0

        k = i - 1

        while k >= 0:
            if base_pos[k] is not None:
                left = base_pos[k]
                break

            k -= 1

        # nearest HXB2 base on the right
        right = HXB2_LEN + 1

        k = j

        while k < n:
            if base_pos[k] is not None:
                right = base_pos[k]
                break

            k += 1

        for rank, col in enumerate(
            range(i, j),
            start=1,
        ):
            if left == 0:
                labels[col] = (
                    f"before1.ins{rank}"
                )
                classes[col] = (
                    "HXB2_edge_insertion"
                )
                anchors[col] = "before1"
                plot_x[col] = 0.5

            elif right == HXB2_LEN + 1:
                labels[col] = (
                    f"after{HXB2_LEN}.ins{rank}"
                )
                classes[col] = (
                    "HXB2_edge_insertion"
                )
                anchors[col] = (
                    f"after{HXB2_LEN}"
                )
                plot_x[col] = (
                    HXB2_LEN + 0.5
                )

            else:
                labels[col] = (
                    f"{left}|{right}.ins{rank}"
                )
                classes[col] = (
                    "HXB2_insertion_segment"
                )
                anchors[col] = (
                    f"{left}|{right}"
                )
                plot_x[col] = (
                    (left + right) / 2.0
                )

        i = j

    return (
        labels,
        classes,
        anchors,
        plot_x,
    )


# ============================================================
# Parse one addfragments alignment
# ============================================================

def parse_alignment(
    patient,
    aligned,
    hxb2_raw,
    query_name,
    query_raw,
):
    if "HXB2_reference" not in aligned:
        raise RuntimeError(
            f"{patient}: HXB2_reference missing from output"
        )

    if query_name not in aligned:
        raise RuntimeError(
            f"{patient}: query sequence missing from output: {query_name}"
        )

    h = aligned["HXB2_reference"]
    q = aligned[query_name]

    if len(h) != len(q):
        raise RuntimeError(
            f"{patient}: aligned sequence lengths differ"
        )

    # Strong sequence-integrity checks.
    h_ungapped = h.replace("-", "")
    q_ungapped = q.replace("-", "")

    if h_ungapped != hxb2_raw:
        raise RuntimeError(
            f"{patient}: HXB2 sequence changed after MAFFT --addfragments"
        )

    if q_ungapped != query_raw:
        raise RuntimeError(
            f"{patient}: mapping consensus sequence changed after MAFFT --addfragments"
        )

    (
        labels,
        classes,
        anchors,
        plot_x,
    ) = build_hxb2_labels(h)

    rows = []

    query_pos = 0

    mapped_hxb2_positions = []

    aligned_acgt_pair_n = 0
    exact_match_n = 0

    query_bases_at_hxb2_gap_n = 0
    hxb2_bases_at_query_gap_n = 0

    for aln_col, (hc, qc) in enumerate(
        zip(h, q),
        start=1,
    ):
        idx = aln_col - 1

        if qc != "-":
            query_pos += 1

            rows.append({
                "Patient":
                    patient,

                "Mapping_consensus_position":
                    query_pos,

                "Pairwise_alignment_column":
                    aln_col,

                "Mapping_consensus_base":
                    qc,

                "HXB2_state":
                    hc,

                "HXB2_coordinate":
                    labels[idx],

                "Coordinate_class":
                    classes[idx],

                "HXB2_anchor":
                    anchors[idx],

                "Plot_coordinate":
                    f"{plot_x[idx]:.3f}",
            })

            if (
                hc in BASES
                and qc in QUERY_BASES
            ):
                mapped_hxb2_positions.append(
                    int(labels[idx])
                )

            if (
                hc in BASES
                and qc in BASES
            ):
                aligned_acgt_pair_n += 1

                if hc == qc:
                    exact_match_n += 1

            if (
                hc == "-"
                and qc in QUERY_BASES
            ):
                query_bases_at_hxb2_gap_n += 1

        elif hc in BASES:
            hxb2_bases_at_query_gap_n += 1

    if query_pos != len(query_raw):
        raise RuntimeError(
            f"{patient}: parsed query length={query_pos}; "
            f"expected={len(query_raw)}"
        )

    if mapped_hxb2_positions:
        start = min(mapped_hxb2_positions)
        end = max(mapped_hxb2_positions)
        span = end - start + 1
    else:
        start = ""
        end = ""
        span = ""

    identity = (
        exact_match_n / aligned_acgt_pair_n
        if aligned_acgt_pair_n
        else float("nan")
    )

    qc_row = {
        "Patient":
            patient,

        "Mapping_consensus_length":
            len(query_raw),

        "Addfragments_alignment_length":
            len(h),

        "HXB2_gap_columns_n":
            h.count("-"),

        "Query_gap_columns_n":
            q.count("-"),

        "HXB2_start":
            start,

        "HXB2_end":
            end,

        "HXB2_span":
            span,

        "Aligned_ACGT_pair_n":
            aligned_acgt_pair_n,

        "Exact_match_n":
            exact_match_n,

        "Identity_ACGT_only":
            (
                f"{identity:.6f}"
                if math.isfinite(identity)
                else ""
            ),

        "Query_bases_at_HXB2_gap_n":
            query_bases_at_hxb2_gap_n,

        "HXB2_bases_at_query_gap_n":
            hxb2_bases_at_query_gap_n,

        "HXB2_sequence_integrity":
            "PASS",

        "Query_sequence_integrity":
            "PASS",
    }

    return rows, qc_row


# ============================================================
# Input checks
# ============================================================

if shutil.which("mafft") is None:
    raise RuntimeError(
        "mafft is not available in PATH"
    )

if not HXB2_SOURCE.exists():
    raise FileNotFoundError(
        HXB2_SOURCE
    )

mapping_files = sorted(
    BRIDGE_PER_PATIENT.glob(
        "*.mapping_consensus.fasta"
    )
)

if not mapping_files:
    raise RuntimeError(
        "No mapping-consensus FASTA files found under:\n"
        f"  {BRIDGE_PER_PATIENT}\n"
        "Run Step04A successfully before this step."
    )


# ============================================================
# HXB2 reference
# ============================================================

hxb2_original_name, hxb2_raw = extract_hxb2()


# ============================================================
# Run all patients
# ============================================================

all_map_rows = []
range_rows = []
qc_rows = []

for mapping_path in mapping_files:

    patient = patient_from_mapping_fasta(
        mapping_path
    )

    mseqs = read_fasta(
        mapping_path
    )

    if len(mseqs) != 1:
        raise RuntimeError(
            f"{patient}: expected exactly one mapping consensus sequence"
        )

    original_query_name = next(
        iter(mseqs)
    )

    query_raw = mseqs[
        original_query_name
    ].replace("-", "")

    query_name = (
        f"{patient}.mapping_consensus"
    )

    # Standardize the query FASTA name.
    standardized_mapping = (
        PER_PATIENT /
        f"{patient}.addfragments.input_mapping_consensus.fasta"
    )

    hxb2_ref_fasta = (
        PER_PATIENT /
        f"{patient}.HXB2_reference.fasta"
    )

    output_fasta = (
        PER_PATIENT /
        f"{patient}.HXB2_plus_mapping_consensus.addfragments.fasta"
    )

    write_fasta(
        standardized_mapping,
        {
            query_name:
                query_raw
        }
    )

    write_fasta(
        hxb2_ref_fasta,
        {
            "HXB2_reference":
                hxb2_raw
        }
    )

    run_addfragments(
        hxb2_ref_fasta,
        standardized_mapping,
        output_fasta,
    )

    aligned = read_fasta(
        output_fasta
    )

    map_rows, qc = parse_alignment(
        patient=patient,
        aligned=aligned,
        hxb2_raw=hxb2_raw,
        query_name=query_name,
        query_raw=query_raw,
    )

    all_map_rows.extend(
        map_rows
    )

    qc_rows.append(
        qc
    )

    range_rows.append({
        "Patient":
            patient,

        "Mapping_consensus_length":
            len(query_raw),

        "HXB2_start":
            qc["HXB2_start"],

        "HXB2_end":
            qc["HXB2_end"],

        "HXB2_span":
            qc["HXB2_span"],

        "HXB2_gap_columns_n":
            qc["HXB2_gap_columns_n"],

        "Query_bases_at_HXB2_gap_n":
            qc["Query_bases_at_HXB2_gap_n"],

        "HXB2_bases_at_query_gap_n":
            qc["HXB2_bases_at_query_gap_n"],

        "Identity_ACGT_only":
            qc["Identity_ACGT_only"],
    })


# ============================================================
# Write V7 primary outputs
# ============================================================

all_map_rows.sort(
    key=lambda r: (
        patient_sort_key(
            r["Patient"]
        ),
        int(
            r["Mapping_consensus_position"]
        ),
    )
)

range_rows.sort(
    key=lambda r: patient_sort_key(
        r["Patient"]
    )
)

qc_rows.sort(
    key=lambda r: patient_sort_key(
        r["Patient"]
    )
)


write_tsv(
    MAP_OUT,
    [
        "Patient",
        "Mapping_consensus_position",
        "Pairwise_alignment_column",
        "Mapping_consensus_base",
        "HXB2_state",
        "HXB2_coordinate",
        "Coordinate_class",
        "HXB2_anchor",
        "Plot_coordinate",
    ],
    all_map_rows,
)


write_tsv(
    RANGE_OUT,
    [
        "Patient",
        "Mapping_consensus_length",
        "HXB2_start",
        "HXB2_end",
        "HXB2_span",
        "HXB2_gap_columns_n",
        "Query_bases_at_HXB2_gap_n",
        "HXB2_bases_at_query_gap_n",
        "Identity_ACGT_only",
    ],
    range_rows,
)


write_tsv(
    QC_OUT,
    [
        "Patient",
        "Mapping_consensus_length",
        "Addfragments_alignment_length",
        "HXB2_gap_columns_n",
        "Query_gap_columns_n",
        "HXB2_start",
        "HXB2_end",
        "HXB2_span",
        "Aligned_ACGT_pair_n",
        "Exact_match_n",
        "Identity_ACGT_only",
        "Query_bases_at_HXB2_gap_n",
        "HXB2_bases_at_query_gap_n",
        "HXB2_sequence_integrity",
        "Query_sequence_integrity",
    ],
    qc_rows,
)


# ============================================================
print("MAFFT --addfragments HXB2 COORDINATE MAPPING COMPLETE")
print("=" * 78)

print()
print("No global all-patient alignment was created.")
print("No SNP / DEL / INS re-calling was performed.")
print("No --keeplength was used.")
print("HXB2 and mapping-consensus sequence integrity were checked.")

print()
print(f"HXB2 source: {hxb2_original_name}")
print(f"HXB2 length: {HXB2_LEN}")

print()
print("Patient HXB2 ranges:")

for r in range_rows:
    print(
        f"{r['Patient']:15s} "
        f"{int(r['HXB2_start']):>5}-"
        f"{int(r['HXB2_end']):<5} "
        f"| HXB2 gaps={int(r['HXB2_gap_columns_n']):>4} "
        f"| identity={r['Identity_ACGT_only']}"
    )

print(MAP_OUT)
print(RANGE_OUT)
print(QC_OUT)

print()
print(
    "Recommended review before mutation re-projection:\n"
    "  02_patient_HXB2_range_summary.tsv\n"
    "  07_addfragments_alignment_QC.tsv\n"
    "and selected per-patient addfragments FASTA files for complex indel cases."
)
