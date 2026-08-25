#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare patient-specific HIV-1 clone alignments.

Purpose
-------
1. Read the combined clone FASTA.
2. Group clone sequences by patient.
3. Add HXB2 and, where available, the supplied manual consensus as
   alignment guide sequences.
4. Align each patient independently using MAFFT.
5. Verify that MAFFT does not alter any ungapped input sequence.

Important
---------
The supplied manual consensus and HXB2 are included only as alignment
guide sequences. They are NOT used as mutation references.

Downstream mutation calling uses the patient-specific analytical consensus
constructed from clone sequences only.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path


HXB2_LENGTH = 9719


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare independent patient-specific HIV-1 clone alignments."
    )

    parser.add_argument(
        "--input-fasta",
        required=True,
        type=Path,
        help="Combined FASTA containing patient clones and optional manual consensuses.",
    )

    parser.add_argument(
        "--hxb2-source",
        required=True,
        type=Path,
        help="FASTA containing the HXB2 reference sequence.",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Output directory.",
    )

    parser.add_argument(
        "--threads",
        default="-1",
        help="Number of MAFFT threads. Default: -1 (automatic/all available).",
    )

    parser.add_argument(
        "--mafft",
        default="mafft",
        help="MAFFT executable. Default: mafft",
    )

    parser.add_argument(
        "--allow-missing-consensus",
        action="append",
        default=[],
        help=(
            "Patient allowed to lack a supplied manual consensus. "
            "May be specified multiple times."
        ),
    )

    return parser.parse_args()


def read_fasta(path: Path):
    seqs = {}
    header = None
    parts = []

    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(parts).upper()

                header = line[1:].strip()
                parts = []
            else:
                parts.append(re.sub(r"\s+", "", line))

        if header is not None:
            seqs[header] = "".join(parts).upper()

    return seqs


def write_fasta(records, path: Path, width=80):
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")

            for i in range(0, len(seq), width):
                handle.write(seq[i:i + width] + "\n")


def ungap(seq: str):
    return seq.replace("-", "").replace(".", "")


def patient_from_header(header: str):
    match = re.search(r"(CN\d{4}AH\d+-\d+)", header)

    if not match:
        return None

    return match.group(1)


def clone_number(header: str):
    match = re.search(r"\.clone(\d+)", header, re.I)

    if not match:
        return None

    return int(match.group(1))


def is_consensus(header: str):
    return "consensus" in header.lower()


def extract_hxb2(path: Path):
    records = read_fasta(path)

    candidates = []

    for header, seq in records.items():
        clean = ungap(seq)

        if "HXB2" in header.upper() and len(clean) == HXB2_LENGTH:
            candidates.append((header, clean))

    unique = {}

    for header, seq in candidates:
        unique[seq] = header

    if len(unique) != 1:
        raise RuntimeError(
            "Could not uniquely identify one HXB2 sequence "
            f"of length {HXB2_LENGTH} in {path}. "
            f"Found {len(unique)} unique candidate sequence(s)."
        )

    sequence = next(iter(unique))
    original_header = unique[sequence]

    return original_header, sequence


def main():
    args = parse_args()

    if shutil.which(args.mafft) is None:
        raise RuntimeError(
            f"MAFFT executable not found in PATH: {args.mafft}"
        )

    if not args.input_fasta.exists():
        raise FileNotFoundError(args.input_fasta)

    if not args.hxb2_source.exists():
        raise FileNotFoundError(args.hxb2_source)

    input_dir = args.output_dir / "input"
    align_dir = args.output_dir / "aligned"
    log_dir = args.output_dir / "logs"
    qc_dir = args.output_dir / "qc"

    for directory in [
        input_dir,
        align_dir,
        log_dir,
        qc_dir,
    ]:
        directory.mkdir(parents=True, exist_ok=True)

    all_records = read_fasta(args.input_fasta)

    patient_clones = defaultdict(list)
    patient_consensus = {}

    for header, seq in all_records.items():
        patient = patient_from_header(header)

        if patient is None:
            continue

        if is_consensus(header):
            if patient in patient_consensus:
                raise RuntimeError(
                    f"More than one supplied consensus detected for {patient}"
                )

            patient_consensus[patient] = (header, seq)

        else:
            number = clone_number(header)

            if number is None:
                continue

            patient_clones[patient].append(
                (number, header, seq)
            )

    for patient in patient_clones:
        patient_clones[patient].sort(
            key=lambda x: x[0]
        )

    patients = sorted(patient_clones)

    if not patients:
        raise RuntimeError(
            "No patient clone sequences were detected."
        )

    allowed_missing = set(args.allow_missing_consensus)

    unexpected_missing = [
        patient
        for patient in patients
        if (
            patient not in patient_consensus
            and patient not in allowed_missing
        )
    ]

    if unexpected_missing:
        raise RuntimeError(
            "The following patients lack a supplied manual consensus "
            "but were not explicitly allowed:\n  "
            + "\n  ".join(unexpected_missing)
        )

    hxb2_original_header, hxb2_sequence = extract_hxb2(
        args.hxb2_source
    )

    print("=" * 72)
    print("PATIENT ALIGNMENT PREPARATION")
    print("=" * 72)
    print(f"Patients detected : {len(patients)}")
    print(
        "Clone sequences   : "
        f"{sum(len(x) for x in patient_clones.values())}"
    )
    print(f"Manual consensuses: {len(patient_consensus)}")
    print(f"HXB2 source       : {hxb2_original_header}")
    print(f"HXB2 length       : {len(hxb2_sequence)}")
    print()

    alignment_files = {}

    summary_rows = []

    for patient in patients:
        records = [
            ("HXB2_reference", hxb2_sequence)
        ]

        if patient in patient_consensus:
            _, seq = patient_consensus[patient]

            records.append(
                (
                    f"{patient}.MANUAL_CONSENSUS",
                    ungap(seq),
                )
            )

        for _, header, seq in patient_clones[patient]:
            records.append(
                (
                    header,
                    ungap(seq),
                )
            )

        input_path = input_dir / f"{patient}.fasta"
        aligned_path = align_dir / f"{patient}.mafft.fasta"
        log_path = log_dir / f"{patient}.mafft.log"

        write_fasta(records, input_path)

        print(
            f"Running MAFFT: {patient} "
            f"({len(patient_clones[patient])} clones)"
        )

        cmd = [
            args.mafft,
            "--auto",
            "--thread",
            str(args.threads),
            str(input_path),
        ]

        with (
            aligned_path.open("w", encoding="utf-8") as out_handle,
            log_path.open("w", encoding="utf-8") as log_handle,
        ):
            subprocess.run(
                cmd,
                stdout=out_handle,
                stderr=log_handle,
                check=True,
            )

        alignment_files[patient] = aligned_path

        summary_rows.append({
            "Patient": patient,
            "Clone_n": len(patient_clones[patient]),
            "Manual_consensus_in_alignment":
                "YES" if patient in patient_consensus else "NO",
            "Input_sequence_n": len(records),
            "Alignment_file": str(aligned_path),
        })

    # ---------------------------------------------------------
    # Alignment sequence-integrity QC
    # ---------------------------------------------------------

    qc_rows = []
    failure_n = 0

    for patient in patients:
        input_path = input_dir / f"{patient}.fasta"
        aligned_path = alignment_files[patient]

        original = read_fasta(input_path)
        aligned = read_fasta(aligned_path)

        for header, seq in original.items():
            if header not in aligned:
                status = "MISSING"
                aligned_ungapped = ""
            else:
                aligned_ungapped = ungap(aligned[header])

                status = (
                    "PASS"
                    if aligned_ungapped == seq
                    else "FAIL"
                )

            if status != "PASS":
                failure_n += 1

            qc_rows.append({
                "Patient": patient,
                "Sequence": header,
                "Input_ungapped_length": len(seq),
                "Aligned_ungapped_length": len(aligned_ungapped),
                "Sequence_integrity": status,
            })

    qc_file = qc_dir / "patient_MAFFT_integrity.tsv"

    with qc_file.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Patient",
                "Sequence",
                "Input_ungapped_length",
                "Aligned_ungapped_length",
                "Sequence_integrity",
            ],
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(qc_rows)

    summary_file = qc_dir / "patient_alignment_summary.tsv"

    with summary_file.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Patient",
                "Clone_n",
                "Manual_consensus_in_alignment",
                "Input_sequence_n",
                "Alignment_file",
            ],
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(summary_rows)

    print()
    print("=" * 72)

    if failure_n:
        print(
            f"MAFFT sequence-integrity QC: FAIL "
            f"({failure_n} sequence(s))"
        )
        print(f"Review: {qc_file}")

        raise SystemExit(1)

    print("MAFFT sequence-integrity QC: PASS")
    print(f"QC table       : {qc_file}")
    print(f"Patient summary: {summary_file}")
    print("=" * 72)


if __name__ == "__main__":
    main()
