#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build patient-specific ANALYSIS_CONSENSUS from existing MAFFT alignments.

Final conservative rules
------------------------
1. Clone A/C/G/T = valid nucleotide allele.
2. Clone N = known DEL-padding in this dataset, therefore treated as gap (-).
3. Internal MAFFT gap (-) = valid gap allele.
4. Terminal MAFFT gap (-) = uncovered/missing, not an allele.
5. Other ambiguity codes = missing.
6. Minimum site-level callable coverage is 90% of ALL patient clones.
7. If callable_fraction < 0.90:
      Analysis consensus = N
      Reason = Insufficient_callable_coverage
8. If callable_fraction >= 0.90, a consensus state is accepted only when one
   A/C/G/T/- state is STRICTLY >50% among callable clones.
9. If coverage passes but no state is >50% among callable clones:
      Analysis consensus = N
      Reason = No_unique_consensus_...
10. SnapGene/manual consensus is retained only for QC comparison.

Consensus decision rules
-------
v2 prevented very low-coverage regions from generating artificial consensus,
but required the major state to exceed 50% of ALL clones. That was too strict
when a small number of clones were legitimately missing. the final workflow separates:
  (a) site-level callable coverage threshold = 90%, and
  (b) majority threshold = strictly >50% among callable clones.
"""

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


# ============================================================
# Paths
# ============================================================

parser = argparse.ArgumentParser(
    description="Build patient-specific analysis consensuses from patient MAFFT alignments."
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
    help="Output directory for consensus FASTA and QC tables."
)

args = parser.parse_args()

ALIGN_DIR = args.alignment_dir
OUTDIR = args.output_dir
FASTA_DIR = OUTDIR / "analysis_consensus_fasta"

OUTDIR.mkdir(parents=True, exist_ok=True)
FASTA_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_OUT = OUTDIR / "01_analysis_consensus_summary.tsv"
SNAP_N_OUT = OUTDIR / "02_snapgene_N_resolution.tsv"
DIFF_OUT = OUTDIR / "03_snapgene_vs_analysis_consensus_differences.tsv"
NOUNIQUE_OUT = OUTDIR / "04_no_unique_consensus_sites.tsv"
LOWCOV_OUT = OUTDIR / "05_insufficient_callable_coverage_sites.tsv"
ALL_CONS_OUT = OUTDIR / "06_all_analysis_consensus_aligned.fasta"


BASES = set("ACGT")
VALID_STATES = ("A", "C", "G", "T", "-")

# Final site-level minimum callable coverage threshold
MIN_CALLABLE_FRACTION = 0.90


# ============================================================
# HXB2 gene annotation
# ============================================================

GENES = {
    "gag": [(790, 2292)],
    "pol": [(2085, 5096)],
    "vif": [(5041, 5619)],
    "vpr": [(5559, 5850)],
    "tat": [(5831, 6045), (8379, 8469)],
    "rev": [(5970, 6045), (8379, 8653)],
    "vpu": [(6062, 6310)],
    "gp120": [(6225, 7757)],
    "gp41": [(7758, 8795)],
    "nef": [(8797, 9417)],
}

GENE_ORDER = [
    "gag", "pol", "vif", "vpr", "tat",
    "rev", "vpu", "gp120", "gp41", "nef"
]


def genes_for_pos(pos):
    if pos is None:
        return "intergenic_or_unassigned"

    hits = []

    for gene in GENE_ORDER:
        for start, end in GENES[gene]:
            if start <= pos <= end:
                hits.append(gene)
                break

    return ",".join(hits) if hits else "intergenic_or_unassigned"


# ============================================================
# FASTA helpers
# ============================================================

def read_fasta(path):
    seqs = {}
    name = None
    parts = []

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts).upper()

                name = line[1:].strip()
                parts = []

            else:
                parts.append(line)

        if name is not None:
            seqs[name] = "".join(parts).upper()

    return seqs


def write_fasta_record(handle, name, seq, width=80):
    handle.write(f">{name}\n")

    for i in range(0, len(seq), width):
        handle.write(seq[i:i + width] + "\n")


# ============================================================
# Header helpers
# ============================================================

def natural_clone_key(name):
    m = re.search(
        r"\.clone(\d+)",
        name,
        flags=re.IGNORECASE
    )

    if m:
        return int(m.group(1))

    return 10**9


def find_clone_names(seqs):
    return sorted(
        [
            name
            for name in seqs
            if re.search(
                r"\.clone\d+",
                name,
                flags=re.IGNORECASE
            )
        ],
        key=natural_clone_key
    )


def find_hxb2_name(seqs):
    if "HXB2_reference" in seqs:
        return "HXB2_reference"

    candidates = [
        name
        for name in seqs
        if "HXB2" in name.upper()
    ]

    if len(candidates) == 1:
        return candidates[0]

    raise RuntimeError(
        "Could not uniquely identify HXB2 reference: "
        f"{candidates}"
    )


def find_manual_consensus_name(seqs, patient):
    preferred = f"{patient}.MANUAL_CONSENSUS"

    if preferred in seqs:
        return preferred

    candidates = [
        name
        for name in seqs
        if "CONSENSUS" in name.upper()
        and ".clone" not in name.lower()
        and "ANALYSIS_CONSENSUS" not in name.upper()
    ]

    if len(candidates) == 1:
        return candidates[0]

    return None


# ============================================================
# MSA -> HXB2 map
# ============================================================

def build_hxb2_map(hxb2):
    mapping = {}

    hpos = 0
    insertion_counter = defaultdict(int)

    for idx, hb in enumerate(hxb2):
        hb = hb.upper()

        if hb in BASES:
            hpos += 1

            mapping[idx] = {
                "column_type": "HXB2_base",
                "coordinate": str(hpos),
                "hxb2_pos": hpos,
            }

        elif hb == "-":
            left = hpos
            right = hpos + 1 if hpos < 9719 else None

            if left == 0:
                anchor = "before_1"
                gene_pos = 1

            elif right is None:
                anchor = "after_9719"
                gene_pos = 9719

            else:
                anchor = f"{left}|{right}"
                gene_pos = right

            insertion_counter[anchor] += 1
            ins_n = insertion_counter[anchor]

            mapping[idx] = {
                "column_type": "HXB2_relative_insertion",
                "coordinate": f"{anchor}.ins{ins_n}",
                "hxb2_pos": gene_pos,
            }

        else:
            raise RuntimeError(
                f"Unexpected HXB2 character at MSA col {idx+1}: {hb!r}"
            )

    if hpos != 9719:
        raise RuntimeError(
            f"HXB2 mapped length = {hpos}; expected 9719"
        )

    return mapping


# ============================================================
# Clone state interpretation
# ============================================================

def clone_span(seq):
    """
    First/last MAFFT non-gap character.
    Terminal '-' outside this span is missing/uncovered.
    """

    indices = [
        i
        for i, c in enumerate(seq)
        if c != "-"
    ]

    if not indices:
        return None, None

    return indices[0], indices[-1]


def interpreted_clone_state(seq, idx, first_non_gap, last_non_gap):
    """
    Returns:
      A/C/G/T : nucleotide allele
      -       : gap/deletion allele
      None    : missing/uncovered

    Project-specific:
      clone N is known DEL-padding and is treated as '-'.
    """

    c = seq[idx].upper()

    if c in BASES:
        return c

    if c == "N":
        return "-"

    if c == "-":
        if (
            first_non_gap is not None
            and first_non_gap <= idx <= last_non_gap
        ):
            return "-"

        return None

    # Other ambiguity symbols are missing.
    return None


# ============================================================
# No-unique class
# ============================================================

def classify_no_unique(counter):
    bases_present = [
        b
        for b in "ACGT"
        if counter[b] > 0
    ]

    gap_present = counter["-"] > 0

    if len(bases_present) >= 2 and not gap_present:
        return "No_unique_consensus_SNP_like"

    if len(bases_present) == 1 and gap_present:
        return "No_unique_consensus_INDEL_like"

    if len(bases_present) >= 2 and gap_present:
        return "No_unique_consensus_mixed"

    return "No_unique_consensus_other"


# ============================================================
# Input discovery
# ============================================================

alignment_files = sorted(
    ALIGN_DIR.glob("*.mafft.fasta")
)

if not alignment_files:
    raise RuntimeError(
        f"No *.mafft.fasta found in {ALIGN_DIR}"
    )


# ============================================================
# Output rows
# ============================================================

summary_rows = []
snap_n_rows = []
diff_rows = []
nounique_rows = []
lowcov_rows = []
analysis_records = []


# ============================================================
# Main
# ============================================================

for path in alignment_files:

    patient = path.name.replace(
        ".mafft.fasta",
        ""
    )

    seqs = read_fasta(path)

    hxb2_name = find_hxb2_name(seqs)
    manual_name = find_manual_consensus_name(
        seqs,
        patient
    )
    clone_names = find_clone_names(seqs)

    if not clone_names:
        raise RuntimeError(
            f"{patient}: no clone sequences found"
        )

    clone_n = len(clone_names)

    hxb2 = seqs[hxb2_name]
    msa_len = len(hxb2)

    lengths = {
        len(seq)
        for seq in seqs.values()
    }

    if len(lengths) != 1:
        raise RuntimeError(
            f"{patient}: inconsistent MSA lengths: "
            f"{sorted(lengths)}"
        )

    mapping = build_hxb2_map(hxb2)

    spans = {
        clone: clone_span(
            seqs[clone]
        )
        for clone in clone_names
    }

    manual = (
        seqs[manual_name]
        if manual_name is not None
        else None
    )

    analysis_chars = []

    analysis_base_n = 0
    analysis_gap_n = 0
    no_unique_n = 0
    lowcov_n = 0
    zero_callable_n = 0

    manual_N_total = 0
    manual_N_to_base = 0
    manual_N_to_gap = 0
    manual_N_to_no_unique = 0
    manual_N_to_lowcov = 0

    for idx in range(msa_len):

        counter = Counter()
        missing_n = 0

        for clone in clone_names:
            first_non_gap, last_non_gap = spans[
                clone
            ]

            state = interpreted_clone_state(
                seqs[clone],
                idx,
                first_non_gap,
                last_non_gap
            )

            if state is None:
                missing_n += 1

            else:
                counter[state] += 1

        callable_n = sum(
            counter[state]
            for state in VALID_STATES
        )

        callable_fraction = (
            callable_n / clone_n
        )

        major_state = max(
            VALID_STATES,
            key=lambda s: counter[s]
        )

        major_n = counter[
            major_state
        ]

        major_fraction_total = (
            major_n / clone_n
        )

        major_fraction_callable = (
            major_n / callable_n
            if callable_n > 0
            else 0.0
        )

        info = mapping[idx]

        # ----------------------------------------------------
        # Final consensus decision
        # ----------------------------------------------------

        if callable_n == 0:

            final_cons = "N"
            reason = "No_callable_state"

            zero_callable_n += 1

            lowcov_rows.append([
                patient,
                idx + 1,
                info["coordinate"],
                info["column_type"],
                genes_for_pos(
                    info["hxb2_pos"]
                ),

                clone_n,

                counter["A"],
                counter["C"],
                counter["G"],
                counter["T"],
                counter["-"],

                missing_n,
                callable_n,
                f"{callable_fraction:.6f}",

                major_state,
                major_n,
                f"{major_fraction_total:.6f}",
                f"{major_fraction_callable:.6f}",

                final_cons,
                reason
            ])

        elif callable_fraction < MIN_CALLABLE_FRACTION:

            final_cons = "N"
            reason = (
                "Insufficient_callable_coverage"
            )

            lowcov_n += 1

            lowcov_rows.append([
                patient,
                idx + 1,
                info["coordinate"],
                info["column_type"],
                genes_for_pos(
                    info["hxb2_pos"]
                ),

                clone_n,

                counter["A"],
                counter["C"],
                counter["G"],
                counter["T"],
                counter["-"],

                missing_n,
                callable_n,
                f"{callable_fraction:.6f}",

                major_state,
                major_n,
                f"{major_fraction_total:.6f}",
                f"{major_fraction_callable:.6f}",

                final_cons,
                reason
            ])

        elif major_fraction_callable > 0.50:

            final_cons = major_state
            reason = (
                "Majority_gt50_of_callable"
            )

            if final_cons == "-":
                analysis_gap_n += 1
            else:
                analysis_base_n += 1

        else:

            final_cons = "N"
            reason = classify_no_unique(
                counter
            )

            no_unique_n += 1

            nounique_rows.append([
                patient,
                idx + 1,
                info["coordinate"],
                info["column_type"],
                genes_for_pos(
                    info["hxb2_pos"]
                ),

                clone_n,

                counter["A"],
                counter["C"],
                counter["G"],
                counter["T"],
                counter["-"],

                missing_n,
                callable_n,
                f"{callable_fraction:.6f}",

                major_state,
                major_n,
                f"{major_fraction_total:.6f}",
                f"{major_fraction_callable:.6f}",

                reason
            ])

        analysis_chars.append(
            final_cons
        )

        # ----------------------------------------------------
        # SnapGene N QC
        # ----------------------------------------------------

        if manual is not None:

            snap = manual[idx].upper()

            if snap == "N":

                manual_N_total += 1

                if final_cons in BASES:
                    manual_N_to_base += 1

                elif final_cons == "-":
                    manual_N_to_gap += 1

                elif reason.startswith(
                    "No_unique_consensus"
                ):
                    manual_N_to_no_unique += 1

                else:
                    manual_N_to_lowcov += 1

                snap_n_rows.append([
                    patient,
                    idx + 1,
                    info["coordinate"],
                    info["column_type"],
                    genes_for_pos(
                        info["hxb2_pos"]
                    ),

                    snap,

                    clone_n,

                    counter["A"],
                    counter["C"],
                    counter["G"],
                    counter["T"],
                    counter["-"],

                    missing_n,
                    callable_n,
                    f"{callable_fraction:.6f}",

                    major_state,
                    major_n,
                    f"{major_fraction_total:.6f}",
                    f"{major_fraction_callable:.6f}",

                    final_cons,
                    reason
                ])

            # all SnapGene-vs-analysis differences
            if snap != final_cons:

                diff_rows.append([
                    patient,
                    idx + 1,
                    info["coordinate"],
                    info["column_type"],
                    genes_for_pos(
                        info["hxb2_pos"]
                    ),

                    snap,
                    final_cons,
                    reason,

                    clone_n,

                    counter["A"],
                    counter["C"],
                    counter["G"],
                    counter["T"],
                    counter["-"],

                    missing_n,
                    callable_n,
                    f"{callable_fraction:.6f}",

                    major_state,
                    major_n,
                    f"{major_fraction_total:.6f}",
                    f"{major_fraction_callable:.6f}"
                ])

    analysis_consensus = "".join(
        analysis_chars
    )

    # --------------------------------------------------------
    # Per-patient aligned analysis consensus FASTA
    # --------------------------------------------------------

    patient_out = (
        FASTA_DIR /
        f"{patient}.analysis_consensus.aligned.fasta"
    )

    with open(
        patient_out,
        "w",
        encoding="utf-8"
    ) as f:

        write_fasta_record(
            f,
            f"{patient}.ANALYSIS_CONSENSUS",
            analysis_consensus
        )

    analysis_records.append(
        (
            f"{patient}.ANALYSIS_CONSENSUS",
            analysis_consensus
        )
    )

    summary_rows.append([
        patient,
        clone_n,
        msa_len,

        manual_name
        if manual_name is not None
        else "NO_MANUAL_CONSENSUS",

        manual_N_total
        if manual is not None
        else "N/A",

        manual_N_to_base
        if manual is not None
        else "N/A",

        manual_N_to_gap
        if manual is not None
        else "N/A",

        manual_N_to_no_unique
        if manual is not None
        else "N/A",

        manual_N_to_lowcov
        if manual is not None
        else "N/A",

        analysis_base_n,
        analysis_gap_n,
        no_unique_n,
        lowcov_n,
        zero_callable_n
    ])


# ============================================================
# Write tables
# ============================================================

with open(
    SUMMARY_OUT,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow([
        "Patient",
        "Clone_n",
        "MSA_length",
        "SnapGene_consensus_name",

        "SnapGene_N_total_n",
        "SnapGene_N_resolved_to_base_n",
        "SnapGene_N_resolved_to_gap_n",
        "SnapGene_N_remained_no_unique_n",
        "SnapGene_N_remained_low_coverage_n",

        "Analysis_consensus_base_column_n",
        "Analysis_consensus_gap_column_n",
        "Analysis_consensus_no_unique_N_column_n",
        "Analysis_consensus_low_coverage_N_column_n",
        "Analysis_consensus_zero_callable_N_column_n"
    ])

    writer.writerows(
        summary_rows
    )


common_header = [
    "Patient",
    "MSA_col",
    "HXB2_coordinate",
    "Column_type",
    "Gene",

    "Clone_n",

    "A_n",
    "C_n",
    "G_n",
    "T_n",
    "Gap_n",

    "Missing_or_terminal_n",
    "Callable_n",
    "Callable_fraction",

    "Major_state",
    "Major_n",
    "Major_fraction_of_all_clones",
    "Major_fraction_of_callable"
]


with open(
    SNAP_N_OUT,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow([
        "Patient",
        "MSA_col",
        "HXB2_coordinate",
        "Column_type",
        "Gene",

        "SnapGene_consensus",

        "Clone_n",

        "A_n",
        "C_n",
        "G_n",
        "T_n",
        "Gap_n",

        "Missing_or_terminal_n",
        "Callable_n",
        "Callable_fraction",

        "Major_state",
        "Major_n",
        "Major_fraction_of_all_clones",
        "Major_fraction_of_callable",

        "Analysis_consensus",
        "Reason"
    ])

    writer.writerows(
        snap_n_rows
    )


with open(
    DIFF_OUT,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow([
        "Patient",
        "MSA_col",
        "HXB2_coordinate",
        "Column_type",
        "Gene",

        "SnapGene_consensus",
        "Analysis_consensus",
        "Reason",

        "Clone_n",

        "A_n",
        "C_n",
        "G_n",
        "T_n",
        "Gap_n",

        "Missing_or_terminal_n",
        "Callable_n",
        "Callable_fraction",

        "Major_state",
        "Major_n",
        "Major_fraction_of_all_clones",
        "Major_fraction_of_callable"
    ])

    writer.writerows(
        diff_rows
    )


with open(
    NOUNIQUE_OUT,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow(
        common_header
        + ["No_unique_type"]
    )

    writer.writerows(
        nounique_rows
    )


with open(
    LOWCOV_OUT,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow(
        common_header
        + [
            "Analysis_consensus",
            "Reason"
        ]
    )

    writer.writerows(
        lowcov_rows
    )


with open(
    ALL_CONS_OUT,
    "w",
    encoding="utf-8"
) as f:

    for name, seq in analysis_records:
        write_fasta_record(
            f,
            name,
            seq
        )


# ============================================================
# Report
# ============================================================

print()
print("=" * 78)
print("PATIENT-SPECIFIC ANALYSIS CONSENSUS COMPLETE")
print("=" * 78)

print()
print("FINAL CONSENSUS RULE:")
print("  clone N          -> DEL-padding -> gap allele (-)")
print("  internal '-'     -> gap allele")
print("  terminal '-'     -> missing")
print(f"  minimum callable -> {MIN_CALLABLE_FRACTION:.0%} of ALL patient clones")
print("  callable <90%    -> N / Insufficient_callable_coverage")
print("  callable >=90%   -> major state must be >50% among callable clones")
print("  no state >50%    -> N / No_unique_consensus")

print()
print("Outputs:")
print(SUMMARY_OUT)
print(SNAP_N_OUT)
print(DIFF_OUT)
print(NOUNIQUE_OUT)
print(LOWCOV_OUT)
print(ALL_CONS_OUT)

print()
print(
    "Existing mutation calls are unchanged. "
    "Do NOT start Phase 2 until these consensus QC tables are reviewed."
)
print()
