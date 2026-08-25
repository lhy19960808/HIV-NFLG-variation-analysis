#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Call clone-level SNP / DEL / INS using the locked
patient-specific ANALYSIS_CONSENSUS.

This script does NOT modify any previous result directory.

Locked interpretation
---------------------
Clone state:
  A/C/G/T        -> nucleotide
  N              -> known DEL-padding -> treated as gap (-)
  internal '-'   -> gap (-)
  terminal '-'   -> missing / uncovered
  other ambiguity codes -> missing

Mutation calls:
  consensus A/C/G/T + clone different A/C/G/T -> SNP
  consensus A/C/G/T + clone gap                -> DEL
  consensus gap       + clone A/C/G/T          -> INS
  consensus gap       + clone gap              -> no mutation
  consensus N                                   -> not directionally callable
  clone missing                                  -> not callable

Indel event:
  consecutive MSA columns with the same DEL or INS state in the same clone
  are grouped into one event.

Important:
- Detailed mutation calls include HXB2-base and HXB2-relative insertion columns.
- The HXB2-base position table contains SNP/DEL/INS calls at actual HXB2 bases.
- HXB2-relative insertion events are summarized separately by anchor.
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
    description="Call within-patient SNPs and indels relative to patient-specific analysis consensuses."
)

parser.add_argument(
    "--alignment-dir",
    required=True,
    type=Path,
    help="Directory containing per-patient *.mafft.fasta alignments."
)

parser.add_argument(
    "--consensus-dir",
    required=True,
    type=Path,
    help="Step 02 output directory containing analysis_consensus_fasta and consensus QC tables."
)

parser.add_argument(
    "--output-dir",
    required=True,
    type=Path,
    help="Output directory for clone-level mutation calls and QC tables."
)

args = parser.parse_args()

ALIGN_DIR = args.alignment_dir

CONS_ROOT = args.consensus_dir
CONS_DIR = CONS_ROOT / "analysis_consensus_fasta"

PHASE1_NOUNIQUE = CONS_ROOT / "04_no_unique_consensus_sites.tsv"
PHASE1_LOWCOV = CONS_ROOT / "05_insufficient_callable_coverage_sites.tsv"

OUTDIR = args.output_dir
OUTDIR.mkdir(parents=True, exist_ok=True)

CALLS_OUT = OUTDIR / "01_clone_mutation_calls.tsv"
EVENTS_OUT = OUTDIR / "02_indel_events.tsv"
PATIENT_SUMMARY_OUT = OUTDIR / "03_patient_mutation_summary.tsv"
HXB2_BASE_OUT = OUTDIR / "04_HXB2_base_position_mutation_counts.tsv"
HXB2_INS_ANCHOR_OUT = OUTDIR / "05_HXB2_relative_insertion_event_counts.tsv"
CLONE_N_QC_OUT = OUTDIR / "06_clone_N_as_DEL_QC.tsv"
QC_OUT = OUTDIR / "07_mutation_calling_QC.tsv"


BASES = set("ACGT")
VALID_CONS = set("ACGT-N")


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


def natural_clone_key(name):
    m = re.search(
        r"\.clone(\d+)",
        name,
        flags=re.IGNORECASE
    )

    return int(m.group(1)) if m else 10**9


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
        f"Cannot uniquely identify HXB2: {candidates}"
    )


def find_analysis_consensus(seqs, patient):
    preferred = f"{patient}.ANALYSIS_CONSENSUS"

    if preferred in seqs:
        return preferred

    candidates = [
        name
        for name in seqs
        if "ANALYSIS_CONSENSUS" in name.upper()
    ]

    if len(candidates) == 1:
        return candidates[0]

    raise RuntimeError(
        f"{patient}: cannot identify ANALYSIS_CONSENSUS "
        f"in {list(seqs)}"
    )


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
                "anchor": "",
                "plot_x": float(hpos),
            }

        elif hb == "-":
            left = hpos
            right = hpos + 1 if hpos < 9719 else None

            if left == 0:
                anchor = "before_1"
                gene_pos = 1
                plot_x = 0.5

            elif right is None:
                anchor = "after_9719"
                gene_pos = 9719
                plot_x = 9719.5

            else:
                anchor = f"{left}|{right}"
                gene_pos = right
                plot_x = left + 0.5

            insertion_counter[anchor] += 1
            ins_n = insertion_counter[anchor]

            mapping[idx] = {
                "column_type": "HXB2_relative_insertion",
                "coordinate": f"{anchor}.ins{ins_n}",
                "hxb2_pos": gene_pos,
                "anchor": anchor,
                "plot_x": plot_x,
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
    First and last MAFFT non-gap characters.
    N is deliberately non-gap here because raw clone N is a real
    sequence placeholder for a known deletion block in this project.
    """

    idxs = [
        i
        for i, c in enumerate(seq)
        if c != "-"
    ]

    if not idxs:
        return None, None

    return idxs[0], idxs[-1]


def interpret_clone_state(seq, idx, first_non_gap, last_non_gap):
    """
    Returns:
      interpreted_state: A/C/G/T/-/None
      state_source: nucleotide / N_as_DEL / internal_gap / terminal_missing /
                    ambiguity_missing
    """

    raw = seq[idx].upper()

    if raw in BASES:
        return raw, "nucleotide"

    if raw == "N":
        return "-", "N_as_DEL"

    if raw == "-":
        if (
            first_non_gap is not None
            and first_non_gap <= idx <= last_non_gap
        ):
            return "-", "internal_gap"

        return None, "terminal_missing"

    return None, "ambiguity_missing"


# ============================================================
# Phase-1 status maps for consensus=N
# ============================================================

def load_phase1_status(path, default_reason):
    out = {}

    if not path.exists():
        return out

    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for row in reader:
            patient = row["Patient"]
            msa_col = int(row["MSA_col"])
            reason = (
                row.get("Reason")
                or row.get("No_unique_type")
                or default_reason
            )

            out[(patient, msa_col)] = reason

    return out


nounique_status = load_phase1_status(
    PHASE1_NOUNIQUE,
    "No_unique_consensus"
)

lowcov_status = load_phase1_status(
    PHASE1_LOWCOV,
    "Insufficient_callable_coverage"
)


# ============================================================
# Output containers
# ============================================================

call_rows = []
event_rows = []
clone_n_qc_rows = []

patient_stats = defaultdict(
    lambda: {
        "clone_n": 0,

        "SNP_call_n_all_columns": 0,
        "SNP_unique_site_all": set(),

        "SNP_call_n_HXB2_base": 0,
        "SNP_unique_site_HXB2_base": set(),

        "DEL_base_call_n": 0,
        "DEL_event_n": 0,
        "DEL_unique_event_locus": set(),

        "INS_base_call_n": 0,
        "INS_event_n": 0,
        "INS_unique_event_locus": set(),

        "clone_N_as_DEL_base_call_n": 0,
        "clone_N_as_DEL_event_n": 0,

        "consensus_no_unique_col_n": 0,
        "consensus_lowcov_col_n": 0,
        "consensus_zero_callable_col_n": 0,
    }
)

hxb2_pos_stats = defaultdict(
    lambda: {
        "genes": set(),

        "SNP_clone_calls": 0,
        "SNP_clones": set(),
        "SNP_patients": set(),

        "DEL_clone_calls": 0,
        "DEL_clones": set(),
        "DEL_patients": set(),

        "INS_clone_calls": 0,
        "INS_clones": set(),
        "INS_patients": set(),
    }
)

ins_anchor_stats = defaultdict(
    lambda: {
        "genes": set(),
        "event_n": 0,
        "clones": set(),
        "patients": set(),
        "total_inserted_base_n": 0,
        "event_lengths": [],
        "plot_x": None,
    }
)

qc_rows = []


# ============================================================
# Indel event helpers
# ============================================================

def genes_union_for_indices(indices, mapping):
    genes = []

    seen = set()

    for idx in indices:
        g = genes_for_pos(
            mapping[idx]["hxb2_pos"]
        )

        for item in g.split(","):
            if item not in seen:
                seen.add(item)
                genes.append(item)

    return ",".join(genes)


def event_locus(event_type, indices, mapping):
    """
    Stable event locus identifier for across-clone unique event counting.
    Uses exact aligned coordinates, preserving HXB2-base vs insertion columns.
    """

    start = mapping[indices[0]]["coordinate"]
    end = mapping[indices[-1]]["coordinate"]

    return f"{event_type}:{start}..{end}"


def summarize_event(patient, clone, event_type, indices, raw_sources, mapping, event_idx):
    start_i = indices[0]
    end_i = indices[-1]

    start_coord = mapping[start_i]["coordinate"]
    end_coord = mapping[end_i]["coordinate"]

    hxb2_base_n = sum(
        1
        for i in indices
        if mapping[i]["column_type"] == "HXB2_base"
    )

    hxb2_ins_col_n = len(indices) - hxb2_base_n

    n_as_del_n = sum(
        1
        for s in raw_sources
        if s == "N_as_DEL"
    )

    genes = genes_union_for_indices(
        indices,
        mapping
    )

    locus = event_locus(
        event_type,
        indices,
        mapping
    )

    return {
        "Patient": patient,
        "Clone": clone,
        "Event_type": event_type,
        "Event_index_within_clone": event_idx,

        "MSA_start": start_i + 1,
        "MSA_end": end_i + 1,

        "HXB2_start_coordinate": start_coord,
        "HXB2_end_coordinate": end_coord,

        "Genes": genes,

        "Event_length_alignment_columns": len(indices),
        "HXB2_base_column_n": hxb2_base_n,
        "HXB2_relative_insertion_column_n": hxb2_ins_col_n,

        "Raw_clone_N_as_DEL_column_n": n_as_del_n,
        "Event_locus_key": locus,
    }


# ============================================================
# Main per-patient loop
# ============================================================

alignment_files = sorted(
    ALIGN_DIR.glob("*.mafft.fasta")
)

if not alignment_files:
    raise RuntimeError(
        f"No *.mafft.fasta files found in {ALIGN_DIR}"
    )


for aln_path in alignment_files:

    patient = aln_path.name.replace(
        ".mafft.fasta",
        ""
    )

    cons_path = (
        CONS_DIR /
        f"{patient}.analysis_consensus.aligned.fasta"
    )

    if not cons_path.exists():
        raise FileNotFoundError(
            f"Missing analysis consensus: {cons_path}"
        )

    aln = read_fasta(
        aln_path
    )

    cons_fa = read_fasta(
        cons_path
    )

    hxb2_name = find_hxb2_name(
        aln
    )

    clone_names = find_clone_names(
        aln
    )

    cons_name = find_analysis_consensus(
        cons_fa,
        patient
    )

    hxb2 = aln[hxb2_name]
    consensus = cons_fa[cons_name]

    msa_len = len(hxb2)

    if len(consensus) != msa_len:
        raise RuntimeError(
            f"{patient}: consensus MSA length {len(consensus)} "
            f"!= patient MSA length {msa_len}"
        )

    for clone in clone_names:
        if len(aln[clone]) != msa_len:
            raise RuntimeError(
                f"{patient} {clone}: clone MSA length mismatch"
            )

    mapping = build_hxb2_map(
        hxb2
    )

    patient_stats[patient]["clone_n"] = len(
        clone_names
    )

    # Consensus-N status counts
    for idx, c in enumerate(consensus):
        if c != "N":
            continue

        key = (patient, idx + 1)

        if key in nounique_status:
            patient_stats[patient][
                "consensus_no_unique_col_n"
            ] += 1

        elif key in lowcov_status:
            reason = lowcov_status[key]

            if reason == "No_callable_state":
                patient_stats[patient][
                    "consensus_zero_callable_col_n"
                ] += 1
            else:
                patient_stats[patient][
                    "consensus_lowcov_col_n"
                ] += 1

    # --------------------------------------------------------
    # Clone-level mutation calling
    # --------------------------------------------------------

    for clone in clone_names:

        seq = aln[clone]

        first_non_gap, last_non_gap = clone_span(
            seq
        )

        n_total_raw = seq.count("N")
        n_called_del = 0

        event_state = {
            "type": None,
            "indices": [],
            "sources": [],
            "event_index": 0,
        }

        def flush_event():
            if not event_state["indices"]:
                return

            event_state["event_index"] += 1

            event = summarize_event(
                patient,
                clone,
                event_state["type"],
                event_state["indices"],
                event_state["sources"],
                mapping,
                event_state["event_index"]
            )

            event_rows.append(event)

            locus = event["Event_locus_key"]

            if event_state["type"] == "DEL":
                patient_stats[patient]["DEL_event_n"] += 1
                patient_stats[patient]["DEL_unique_event_locus"].add(locus)

                if event["Raw_clone_N_as_DEL_column_n"] > 0:
                    patient_stats[patient]["clone_N_as_DEL_event_n"] += 1

            elif event_state["type"] == "INS":
                patient_stats[patient]["INS_event_n"] += 1
                patient_stats[patient]["INS_unique_event_locus"].add(locus)

                col_types = {
                    mapping[i]["column_type"]
                    for i in event_state["indices"]
                }

                anchors = {
                    mapping[i]["anchor"]
                    for i in event_state["indices"]
                    if mapping[i]["column_type"] == "HXB2_relative_insertion"
                }

                if (
                    col_types == {"HXB2_relative_insertion"}
                    and len(anchors) == 1
                ):
                    anchor = next(iter(anchors))
                    stat = ins_anchor_stats[anchor]

                    stat["event_n"] += 1
                    stat["clones"].add(clone)
                    stat["patients"].add(patient)
                    stat["total_inserted_base_n"] += len(event_state["indices"])
                    stat["event_lengths"].append(len(event_state["indices"]))

                    for i in event_state["indices"]:
                        stat["genes"].update(
                            genes_for_pos(mapping[i]["hxb2_pos"]).split(",")
                        )

                    stat["plot_x"] = mapping[
                        event_state["indices"][0]
                    ]["plot_x"]

            event_state["type"] = None
            event_state["indices"] = []
            event_state["sources"] = []

        for idx in range(msa_len):

            cons = consensus[idx]

            if cons not in VALID_CONS:
                raise RuntimeError(
                    f"{patient}: unexpected analysis consensus state "
                    f"{cons!r} at MSA col {idx+1}"
                )

            clone_state, source = interpret_clone_state(
                seq,
                idx,
                first_non_gap,
                last_non_gap
            )

            mutation = None

            if cons == "N":
                mutation = None

            elif clone_state is None:
                mutation = None

            elif cons in BASES:

                if clone_state in BASES:
                    if clone_state != cons:
                        mutation = "SNP"

                elif clone_state == "-":
                    mutation = "DEL"

            elif cons == "-":

                if clone_state in BASES:
                    mutation = "INS"

            # ---------------- mutation call row ----------------

            if mutation is not None:

                info = mapping[idx]
                genes = genes_for_pos(
                    info["hxb2_pos"]
                )

                if mutation == "SNP":
                    change = f"{cons}>{clone_state}"

                elif mutation == "DEL":
                    change = f"{cons}>-"

                else:
                    change = f"->{clone_state}"

                call_rows.append([
                    patient,
                    clone,
                    idx + 1,

                    info["coordinate"],
                    info["column_type"],
                    info["hxb2_pos"],
                    info["anchor"],
                    genes,

                    cons,
                    seq[idx],
                    clone_state,
                    source,

                    mutation,
                    change,
                ])

                # patient stats
                if mutation == "SNP":

                    patient_stats[patient][
                        "SNP_call_n_all_columns"
                    ] += 1

                    patient_stats[patient][
                        "SNP_unique_site_all"
                    ].add(idx + 1)

                    if info[
                        "column_type"
                    ] == "HXB2_base":

                        patient_stats[patient][
                            "SNP_call_n_HXB2_base"
                        ] += 1

                        patient_stats[patient][
                            "SNP_unique_site_HXB2_base"
                        ].add(
                            info["hxb2_pos"]
                        )

                elif mutation == "DEL":

                    patient_stats[patient][
                        "DEL_base_call_n"
                    ] += 1

                    if source == "N_as_DEL":

                        patient_stats[patient][
                            "clone_N_as_DEL_base_call_n"
                        ] += 1

                        n_called_del += 1

                elif mutation == "INS":

                    patient_stats[patient][
                        "INS_base_call_n"
                    ] += 1

                # HXB2-base position stats
                if info[
                    "column_type"
                ] == "HXB2_base":

                    p = info["hxb2_pos"]
                    st = hxb2_pos_stats[p]

                    st["genes"].update(
                        genes.split(",")
                    )

                    if mutation == "SNP":
                        st[
                            "SNP_clone_calls"
                        ] += 1
                        st[
                            "SNP_clones"
                        ].add(
                            (patient, clone)
                        )
                        st[
                            "SNP_patients"
                        ].add(patient)

                    elif mutation == "DEL":
                        st[
                            "DEL_clone_calls"
                        ] += 1
                        st[
                            "DEL_clones"
                        ].add(
                            (patient, clone)
                        )
                        st[
                            "DEL_patients"
                        ].add(patient)

                    elif mutation == "INS":
                        st[
                            "INS_clone_calls"
                        ] += 1
                        st[
                            "INS_clones"
                        ].add(
                            (patient, clone)
                        )
                        st[
                            "INS_patients"
                        ].add(patient)

            # ---------------- event grouping ----------------

            if mutation in {"DEL", "INS"}:

                if (
                    event_state["type"] == mutation
                    and event_state["indices"]
                    and idx == event_state["indices"][-1] + 1
                ):
                    event_state["indices"].append(idx)
                    event_state["sources"].append(source)

                else:
                    flush_event()

                    event_state["type"] = mutation
                    event_state["indices"] = [idx]
                    event_state["sources"] = [source]

            else:
                flush_event()

        flush_event()

        clone_n_qc_rows.append([
            patient,
            clone,
            n_total_raw,
            n_called_del,
            (
                "PASS"
                if n_called_del <= n_total_raw
                else "FAIL"
            )
        ])

    qc_rows.append([
        patient,
        len(clone_names),
        msa_len,
        "PASS"
    ])


# ============================================================
# Write detailed mutation calls
# ============================================================

with open(
    CALLS_OUT,
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
        "Clone",
        "MSA_col",

        "HXB2_coordinate",
        "Column_type",
        "HXB2_gene_mapping_position",
        "HXB2_insertion_anchor",
        "Gene",

        "Analysis_consensus_state",
        "Clone_raw_state",
        "Clone_interpreted_state",
        "Clone_state_source",

        "Mutation_type",
        "Change"
    ])

    writer.writerows(
        call_rows
    )


# ============================================================
# Write indel events
# ============================================================

with open(
    EVENTS_OUT,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    fieldnames = [
        "Patient",
        "Clone",
        "Event_type",
        "Event_index_within_clone",

        "MSA_start",
        "MSA_end",

        "HXB2_start_coordinate",
        "HXB2_end_coordinate",

        "Genes",

        "Event_length_alignment_columns",
        "HXB2_base_column_n",
        "HXB2_relative_insertion_column_n",

        "Raw_clone_N_as_DEL_column_n",
        "Event_locus_key"
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        event_rows
    )


# ============================================================
# Patient summary
# ============================================================

with open(
    PATIENT_SUMMARY_OUT,
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

        "SNP_clone_call_n_all_alignment_columns",
        "SNP_unique_site_n_all_alignment_columns",

        "SNP_clone_call_n_HXB2_base_only",
        "SNP_unique_site_n_HXB2_base_only",

        "DEL_base_call_n",
        "DEL_event_occurrence_n",
        "DEL_unique_event_locus_n",

        "INS_base_call_n",
        "INS_event_occurrence_n",
        "INS_unique_event_locus_n",

        "Clone_N_as_DEL_base_call_n",
        "Clone_N_as_DEL_event_occurrence_n",

        "Consensus_no_unique_column_n",
        "Consensus_low_coverage_column_n",
        "Consensus_zero_callable_column_n"
    ])

    for patient in sorted(
        patient_stats
    ):
        s = patient_stats[patient]

        writer.writerow([
            patient,
            s["clone_n"],

            s[
                "SNP_call_n_all_columns"
            ],
            len(
                s["SNP_unique_site_all"]
            ),

            s[
                "SNP_call_n_HXB2_base"
            ],
            len(
                s["SNP_unique_site_HXB2_base"]
            ),

            s["DEL_base_call_n"],
            s["DEL_event_n"],
            len(
                s["DEL_unique_event_locus"]
            ),

            s["INS_base_call_n"],
            s["INS_event_n"],
            len(
                s["INS_unique_event_locus"]
            ),

            s[
                "clone_N_as_DEL_base_call_n"
            ],
            s[
                "clone_N_as_DEL_event_n"
            ],

            s[
                "consensus_no_unique_col_n"
            ],
            s[
                "consensus_lowcov_col_n"
            ],
            s[
                "consensus_zero_callable_col_n"
            ],
        ])


# ============================================================
# HXB2-base position counts
# ============================================================

with open(
    HXB2_BASE_OUT,
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
        "HXB2_position",
        "Gene",

        "SNP_clone_call_n",
        "SNP_variant_clone_n",
        "SNP_variant_patient_n",

        "DEL_clone_call_n",
        "DEL_variant_clone_n",
        "DEL_variant_patient_n",

        "INS_clone_call_n",
        "INS_variant_clone_n",
        "INS_variant_patient_n"
    ])

    for pos in sorted(
        hxb2_pos_stats
    ):
        s = hxb2_pos_stats[pos]

        writer.writerow([
            pos,
            ",".join(
                sorted(
                    s["genes"]
                )
            ),

            s[
                "SNP_clone_calls"
            ],
            len(
                s["SNP_clones"]
            ),
            len(
                s["SNP_patients"]
            ),

            s[
                "DEL_clone_calls"
            ],
            len(
                s["DEL_clones"]
            ),
            len(
                s["DEL_patients"]
            ),

            s[
                "INS_clone_calls"
            ],
            len(
                s["INS_clones"]
            ),
            len(
                s["INS_patients"]
            ),
        ])


# ============================================================
# HXB2-relative insertion event counts
# ============================================================

def anchor_sort_key(anchor):
    if anchor == "before_1":
        return (-1, 0)

    if anchor == "after_9719":
        return (9720, 0)

    left, right = anchor.split("|")

    return (
        int(left),
        int(right)
    )


with open(
    HXB2_INS_ANCHOR_OUT,
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
        "HXB2_insertion_anchor",
        "Plot_x",
        "Gene",

        "INS_event_occurrence_n",
        "INS_variant_clone_n",
        "INS_variant_patient_n",

        "INS_total_inserted_base_n",
        "INS_min_event_length",
        "INS_max_event_length"
    ])

    for anchor in sorted(
        ins_anchor_stats,
        key=anchor_sort_key
    ):
        s = ins_anchor_stats[
            anchor
        ]

        writer.writerow([
            anchor,
            s["plot_x"],
            ",".join(
                sorted(
                    s["genes"]
                )
            ),

            s["event_n"],
            len(
                s["clones"]
            ),
            len(
                s["patients"]
            ),

            s[
                "total_inserted_base_n"
            ],
            min(
                s["event_lengths"]
            ),
            max(
                s["event_lengths"]
            ),
        ])


# ============================================================
# Clone N-as-DEL QC
# ============================================================

with open(
    CLONE_N_QC_OUT,
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
        "Clone",
        "Raw_clone_N_total_n",
        "Raw_clone_N_called_as_DEL_n",
        "QC"
    ])

    writer.writerows(
        clone_n_qc_rows
    )


# ============================================================
# Phase2 QC
# ============================================================

with open(
    QC_OUT,
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
        "Alignment_and_consensus_length_QC"
    ])

    writer.writerows(
        qc_rows
    )


# ============================================================
# Console summary
# ============================================================

print()
print("=" * 78)
print("WITHIN-PATIENT MUTATION CALLING COMPLETE")
print("=" * 78)

print()
print("Outputs:")
print(CALLS_OUT)
print(EVENTS_OUT)
print(PATIENT_SUMMARY_OUT)
print(HXB2_BASE_OUT)
print(HXB2_INS_ANCHOR_OUT)
print(CLONE_N_QC_OUT)
print(QC_OUT)

print()
print(f"Directional mutation calls: {len(call_rows)}")
print(f"Indel events: {len(event_rows)}")

print()
print(
    "Next: inspect patient summary, P19 INS events, "
    "P6 clone-N deletion QC, and HXB2 position tables "
    "before making final gene/figure tables."
)
print()
