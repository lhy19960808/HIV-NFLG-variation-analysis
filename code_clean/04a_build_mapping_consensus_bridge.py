#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mapping-consensus bridge construction
=============================================================

Purpose
-------
Validate HXB2 coordinates independently for each patient, without creating a
global all-patient alignment and without changing the existing within-patient
clone alignment.

For each patient:

    existing patient clone MSA
            |
            | reconstruct locked patient-specific analysis consensus
            v
    gapped analysis consensus
            |
            | create a coordinate-mapping sequence
            | (remove consensus-gap / no-coverage columns)
            v
    patient mapping consensus  +  HXB2
            |
            | patient-specific mapping-consensus construction
            v
    pairwise alignment
            |
            +--> consensus nucleotide position <-> HXB2 coordinate
            +--> patient HXB2 start/end QC
            +--> comparison with previous V5 profile-based mapping
            +--> comparison with minimap2 endpoint QC

IMPORTANT
---------
This script performs COORDINATE VALIDATION ONLY.

It does NOT:
- re-call SNP / DEL / INS;
- change the patient-specific consensus rule;
- change the original patient clone MSA;
- create one global clone/HXB2 MSA.

Mutation definition remains:
    clone vs patient-specific analysis consensus

HXB2 remains:
    coordinate reference only

patient-specific consensus rule
-----------------
At each original patient-MSA column:
1. >=90% clones must be callable.
2. A/C/G/T and internal gap are callable states.
3. Clone N is treated as gap allele (known DEL padding).
4. Terminal alignment gap is missing.
5. Strict >50% among callable states is required.
6. Otherwise consensus = N.

Coordinate-mapping sequence
---------------------------
For pairwise alignment to HXB2:
- V3 consensus A/C/G/T columns are retained.
- "No unique consensus" columns containing BASE states only are retained as N
  so that known nucleotide presence is not deleted from the mapping sequence.
- Consensus-gap columns are excluded.
- No-callable and <90%-callable N columns are excluded.
- Base-vs-gap no-unique columns are excluded and flagged as structural
  unresolved columns.

This mapping sequence is ONLY a coordinate-validation surrogate; it does not
replace the analytical consensus.

Inputs
------
10_patient_alignment/aligned/*.mafft.fasta
00_raw/all_patients_NFLG_with_HXB2_original.fasta

Optional comparison inputs
--------------------------
20_patient_hxb2_coordinate_map_v5/03_patient_HXB2_range_summary_v5.tsv
18_hxb2_mapping_qc/04_clone_HXB2_mapping_summary.tsv

Outputs
-------
21_pairwise_consensus_HXB2_v6/
  01_pairwise_consensus_position_to_HXB2.tsv
  02_original_MSA_column_consensus_index.tsv
  03_patient_pairwise_HXB2_range_summary.tsv
  06_pairwise_alignment_QC.tsv
  per_patient/
      <patient>.mapping_consensus.fasta
      <patient>.HXB2_plus_mapping_consensus.input.fasta
      <patient>.HXB2_plus_mapping_consensus.mafft.fasta
"""

from pathlib import Path
from collections import Counter, defaultdict
import argparse
import csv
import math
import re


# ============================================================
# Paths
# ============================================================

parser = argparse.ArgumentParser(
    description="Build patient mapping consensuses and the original-MSA-column bridge for HXB2 coordinate projection."
)

parser.add_argument(
    "--alignment-dir",
    required=True,
    type=Path,
    help="Directory containing per-patient *.mafft.fasta alignments."
)


parser.add_argument(
    "--output-dir",
    required=True,
    type=Path,
    help="Output directory for mapping consensuses, bridge tables, and QC."
)



args = parser.parse_args()

ALIGN_DIR = args.alignment_dir



OUTDIR = args.output_dir
PER_PATIENT = OUTDIR / "per_patient"

OUTDIR.mkdir(parents=True, exist_ok=True)
PER_PATIENT.mkdir(parents=True, exist_ok=True)

SUMMARY_OUT = OUTDIR / "01_mapping_consensus_summary.tsv"
MSA_INDEX_OUT = OUTDIR / "02_original_MSA_column_consensus_index.tsv"

MIN_CALLABLE_FRACTION = 0.90
BASES = set("ACGT")


# ============================================================
# Helpers
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


def patient_from_filename(path):
    return path.name.replace(".mafft.fasta", "")


def clone_number(name):
    m = re.search(r"\.clone(\d+)", name, re.I)

    if m:
        return int(m.group(1))

    return 10**9


def get_patient_clones(seqs, patient):
    clones = [
        n
        for n in seqs
        if n.startswith(patient + ".clone")
    ]

    return sorted(
        clones,
        key=clone_number,
    )


def first_last_non_gap(seq):
    idxs = [
        i
        for i, c in enumerate(seq)
        if c != "-"
    ]

    if not idxs:
        return None, None

    return idxs[0], idxs[-1]


def callable_state(char, idx, first, last):
    """
    Locked project semantics.

    Returns:
      A/C/G/T
      "-"
      None  (missing / uncallable)
    """

    c = char.upper()

    if first is None:
        return None

    # Outside actual sequence span = terminal missing.
    if idx < first or idx > last:
        return None

    if c in BASES:
        return c

    # Project-specific known deletion padding.
    if c == "N":
        return "-"

    # Internal alignment gap = structural gap allele.
    if c == "-":
        return "-"

    # Other ambiguity.
    return None


def build_analysis_consensus(clone_seqs):
    names = list(clone_seqs)

    if not names:
        raise RuntimeError("No clones supplied")

    lengths = {
        len(clone_seqs[n])
        for n in names
    }

    if len(lengths) != 1:
        raise RuntimeError(
            f"Clone alignment lengths differ: {sorted(lengths)}"
        )

    aln_len = next(iter(lengths))
    clone_n = len(names)

    spans = {
        n: first_last_non_gap(clone_seqs[n])
        for n in names
    }

    consensus = []
    meta = []

    for idx in range(aln_len):

        states = []

        for name in names:
            first, last = spans[name]

            state = callable_state(
                clone_seqs[name][idx],
                idx,
                first,
                last,
            )

            if state is not None:
                states.append(state)

        callable_n = len(states)
        callable_fraction = (
            callable_n / clone_n
            if clone_n
            else 0.0
        )

        counts = Counter(states)

        if callable_n == 0:
            cons = "N"
            reason = "No_callable_state"

        elif callable_fraction < MIN_CALLABLE_FRACTION:
            cons = "N"
            reason = "Insufficient_callable_coverage"

        else:
            major_state, major_n = counts.most_common(1)[0]

            if major_n / callable_n > 0.50:
                cons = major_state
                reason = "Majority_state"
            else:
                cons = "N"
                reason = "No_unique_consensus"

        consensus.append(cons)

        meta.append({
            "callable_n": callable_n,
            "clone_n": clone_n,
            "callable_fraction": callable_fraction,
            "reason": reason,
            "counts": counts,
        })

    return "".join(consensus), meta


def build_mapping_consensus(analysis_consensus, meta):
    """
    Build ungapped sequence used ONLY for pairwise coordinate validation.

    Retain:
      - analysis-consensus A/C/G/T
      - base-only no-unique columns as N

    Exclude:
      - consensus gap
      - insufficient coverage N
      - no-callable N
      - base-vs-gap no-unique columns

    Returns:
      mapping_sequence
      original_msa_col_to_mapping_pos
      mapping_pos_to_original_msa_col
      inclusion labels
    """

    chars = []

    msa_to_pos = {}
    pos_to_msa = {}
    labels = {}

    mapping_pos = 0

    for idx, (cons, m) in enumerate(
        zip(analysis_consensus, meta)
    ):
        counts = m["counts"]
        reason = m["reason"]

        include = False
        map_char = None
        category = None

        if cons in BASES:
            include = True
            map_char = cons
            category = "Analysis_consensus_base"

        elif (
            cons == "N"
            and reason == "No_unique_consensus"
            and counts
            and set(counts).issubset(BASES)
        ):
            # Presence is unambiguous, nucleotide identity is not.
            include = True
            map_char = "N"
            category = "No_unique_consensus_base_only"

        elif cons == "-":
            category = "Analysis_consensus_gap_excluded"

        elif (
            cons == "N"
            and reason == "No_unique_consensus"
            and "-" in counts
        ):
            category = "No_unique_consensus_structural_excluded"

        elif (
            cons == "N"
            and reason == "Insufficient_callable_coverage"
        ):
            category = "Low_coverage_excluded"

        elif (
            cons == "N"
            and reason == "No_callable_state"
        ):
            category = "No_callable_excluded"

        else:
            category = "Other_excluded"

        if include:
            mapping_pos += 1
            chars.append(map_char)

            msa_to_pos[idx] = mapping_pos
            pos_to_msa[mapping_pos] = idx
        else:
            msa_to_pos[idx] = None

        labels[idx] = category

    return (
        "".join(chars),
        msa_to_pos,
        pos_to_msa,
        labels,
    )


def patient_sort_key(patient):
    m = re.fullmatch(
        r"CN(\d{4})AH(\d+)-(\d+)",
        patient
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


# ============================================================
# Input checks
# ============================================================

patient_files = sorted(
    ALIGN_DIR.glob("*.mafft.fasta")
)

if not patient_files:
    raise RuntimeError(
        f"No patient alignments found in {ALIGN_DIR}"
    )


# ============================================================
# Main
# ============================================================

summary_rows = []
msa_index_rows = []

for aln_path in patient_files:

    patient = patient_from_filename(
        aln_path
    )

    seqs = read_fasta(
        aln_path
    )

    clones = get_patient_clones(
        seqs,
        patient
    )

    if not clones:
        print(
            f"WARNING: no clones found for {patient}; skipping"
        )
        continue

    clone_seqs = {
        c: seqs[c]
        for c in clones
    }

    lengths = {
        len(s)
        for s in clone_seqs.values()
    }

    if len(lengths) != 1:
        raise RuntimeError(
            f"{patient}: clone MSA lengths differ"
        )

    original_msa_len = next(
        iter(lengths)
    )

    # --------------------------------------------------------
    # Reconstruct locked analytical consensus.
    # --------------------------------------------------------

    analysis_consensus, meta = build_analysis_consensus(
        clone_seqs
    )

    if len(analysis_consensus) != original_msa_len:
        raise RuntimeError(
            f"{patient}: analysis consensus length mismatch"
        )

    # --------------------------------------------------------
    # Build mapping-only consensus sequence.
    # --------------------------------------------------------

    (
        mapping_consensus,
        msa_to_pos,
        pos_to_msa,
        inclusion_labels,
    ) = build_mapping_consensus(
        analysis_consensus,
        meta,
    )

    if not mapping_consensus:
        raise RuntimeError(
            f"{patient}: empty mapping consensus"
        )

    # MSA index output before pairwise mapping.
    category_counts = Counter()

    for idx in range(original_msa_len):

        category = inclusion_labels[idx]
        category_counts[category] += 1

        m = meta[idx]

        counts_text = ",".join(
            f"{state}:{count}"
            for state, count in sorted(
                m["counts"].items()
            )
        )

        msa_index_rows.append({
            "Patient": patient,
            "Original_patient_MSA_column": idx + 1,

            "Analysis_consensus_state":
                analysis_consensus[idx],

            "Consensus_reason":
                m["reason"],

            "Callable_clone_n":
                m["callable_n"],

            "Total_clone_n":
                m["clone_n"],

            "Callable_fraction":
                f"{m['callable_fraction']:.6f}",

            "Callable_state_counts":
                counts_text,

            "Mapping_consensus_position":
                (
                    msa_to_pos[idx]
                    if msa_to_pos[idx] is not None
                    else ""
                ),

            "Mapping_inclusion_class":
                category,
        })

    # --------------------------------------------------------
    # Write patient mapping consensus.
    # --------------------------------------------------------

    query_name = f"{patient}.mapping_consensus"

    mapping_consensus_path = (
        PER_PATIENT /
        f"{patient}.mapping_consensus.fasta"
    )

    write_fasta(
        mapping_consensus_path,
        {
            query_name:
                mapping_consensus
        }
    )

    # --------------------------------------------------------
    # Mapping-consensus construction summary.
    # --------------------------------------------------------

    summary_rows.append({
        "Patient":
            patient,

        "Clone_n":
            len(clones),

        "Original_patient_MSA_length":
            original_msa_len,

        "Mapping_consensus_length":
            len(mapping_consensus),

        "Analysis_consensus_base_columns":
            category_counts[
                "Analysis_consensus_base"
            ],

        "No_unique_base_only_columns_retained_as_N":
            category_counts[
                "No_unique_consensus_base_only"
            ],

        "Consensus_gap_columns_excluded":
            category_counts[
                "Analysis_consensus_gap_excluded"
            ],

        "No_unique_structural_columns_excluded":
            category_counts[
                "No_unique_consensus_structural_excluded"
            ],

        "Low_coverage_columns_excluded":
            category_counts[
                "Low_coverage_excluded"
            ],

        "No_callable_columns_excluded":
            category_counts[
                "No_callable_excluded"
            ],
    })


# ============================================================
# Write primary outputs
# ============================================================

summary_rows.sort(
    key=lambda r: patient_sort_key(
        r["Patient"]
    )
)

msa_index_rows.sort(
    key=lambda r: (
        patient_sort_key(r["Patient"]),
        int(r["Original_patient_MSA_column"]),
    )
)


write_tsv(
    SUMMARY_OUT,
    [
        "Patient",
        "Clone_n",
        "Original_patient_MSA_length",
        "Mapping_consensus_length",
        "Analysis_consensus_base_columns",
        "No_unique_base_only_columns_retained_as_N",
        "Consensus_gap_columns_excluded",
        "No_unique_structural_columns_excluded",
        "Low_coverage_columns_excluded",
        "No_callable_columns_excluded",
    ],
    summary_rows,
)


write_tsv(
    MSA_INDEX_OUT,
    [
        "Patient",
        "Original_patient_MSA_column",
        "Analysis_consensus_state",
        "Consensus_reason",
        "Callable_clone_n",
        "Total_clone_n",
        "Callable_fraction",
        "Callable_state_counts",
        "Mapping_consensus_position",
        "Mapping_inclusion_class",
    ],
    msa_index_rows,
)


# ============================================================
print("=" * 78)
print("MAPPING-CONSENSUS BRIDGE CONSTRUCTION COMPLETE")
print("=" * 78)

print()
print("No HXB2 coordinate assignment was performed in this step.")
print("No global all-patient alignment was created.")
print("No SNP / DEL / INS re-calling was performed.")
print("Original patient clone MSAs were not modified.")

print()
print(f"Patients processed: {len(summary_rows)}")

for r in summary_rows:
    print(
        f"{r['Patient']:15s} "
        f"| mapping consensus="
        f"{r['Mapping_consensus_length']} nt "
        f"| structural no-unique excluded="
        f"{r['No_unique_structural_columns_excluded']}"
    )

print()
print("Outputs:")
print(SUMMARY_OUT)
print(MSA_INDEX_OUT)
print(PER_PATIENT)

print()
print(
    "Next step: map each patient mapping consensus to HXB2 "
    "using Step04B MAFFT --addfragments."
)
