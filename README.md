# HIV-1 NFLG Within-Host Variation Analysis

## Overview

This repository contains the reproducible analysis pipeline used to characterize within-host nucleotide variation among HIV-1 near-full-length genome (NFLG) clone sequences.

Clone-specific SNPs, insertions, and deletions are defined relative to the patient-specific consensus sequence and subsequently projected onto HXB2 coordinates for genomic localization.

HXB2 is used only as a standardized coordinate and gene-annotation framework and is not used as the mutation reference.

## Workflow

The production pipeline consists of nine scripts:

1. `01_prepare_patient_alignments.py`
2. `02_build_patient_consensus.py`
3. `03_call_within_patient_mutations.py`
4. `04a_build_mapping_consensus_bridge.py`
5. `04b_map_patient_consensus_to_HXB2.py`
6. `05_project_mutations_to_HXB2.py`
7. `06_prepare_patient_plot_tables.py`
8. `07_prepare_cross_patient_dataset.py`
9. `08_prepare_figure2_plot_data.py`

Patient sequences are processed independently. No global all-patient multiple-sequence alignment is used for mutation definition.

## Software

The minimal environment contains:

- Python 3.14.6
- MAFFT 7.526

The production Python scripts use only the Python standard library.

Create the environment with:

```bash
conda env create -f environment.yml
conda activate hiv_nflg
```

## Input files

The default runner expects:

```text
data/All_185_clones_no_vector_no_primer_plus_consensus_renamed.fasta
data/all_patients_NFLG_with_HXB2_original.fasta
```

Configuration files are located in `config/`.

## Run the pipeline

From the project root:

```bash
conda activate hiv_nflg
bash ./run_pipeline.sh
```

A custom output directory can be supplied:

```bash
bash ./run_pipeline.sh /path/to/output
```

The default output directory is:

```text
analysis_output/
```

## Main output stages

```text
step01_patient_alignments/
step02_patient_consensus/
step03_mutation_calls/
step04a_mapping_bridge/
step04b_hxb2_mapping/
step05_final_mutation_projection/
step06_patient_plot_tables/
step07_cross_patient_dataset/
step08_figure2_plot_data/
logs/
```

Step 05 contains the final clone-level HXB2-projected mutation calls.
Step 06 contains patient-level plotting tables.
Step 07 contains the cross-patient recurrence dataset.
Step 08 contains final Figure 2 plotting tables.

## Frozen reference results

The finalized dataset contains:

```text
Patients                         11
Directional mutation calls      3287
SNP calls                        2477
DEL base calls                    592
INS base calls                    218
DEL events                         44
INS events                         31
Total INDEL events                 75

SNP loci/anchors                  973
Recurrent SNP loci                 52
Same-substitution recurrent        21
Mixed-substitution recurrent       31
DEL loci                            32
INS anchors                         13
```

The large internal missing region in `CN2023AH1-18` corresponds to HXB2 positions 3487-5903 (2417 nt) and is treated as a missing-region annotation, not as a clone-versus-consensus deletion event.

## Interpretation

For example, `A>G` means:

```text
patient-specific consensus A -> clone G
```

It does not represent an ancestral or evolutionary substitution direction.

Terminal alignment gaps caused by incomplete sequence coverage are treated as missing/unsequenced data rather than deletions.

DEL and INS summary counts represent contiguous clone-specific event occurrences rather than affected nucleotide counts.

## HXB2 projection

Mutation calling and coordinate assignment are separate operations.

Step 03 defines mutations relative to the patient-specific consensus.
Step 04A constructs the mapping-consensus bridge.
Step 04B maps each patient mapping consensus independently to HXB2 using MAFFT `--addfragments`.
Step 05 projects the already-defined mutations onto the HXB2 coordinate framework.

## Cross-patient recurrence

The primary cross-patient metric is `Patient_n`, defined as the number of unique patients showing within-host variation at the same HXB2 locus or anchor.

Insertion recurrence is compared using the common HXB2 anchor rather than assuming that `.ins1`, `.ins2`, and related extended coordinates are homologous among patients.

## Validation

After a complete run:

```bash
python tests/check_release_output.py \
    analysis_output
```

A successful frozen-data reproduction ends with:

```text
RELEASE OUTPUT VALIDATION: PASS
```

The frozen production code can be checked with:

```bash
sha256sum -c tests/core_code.sha256
```

## Citation

Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution. 2013;30(4):772-780. doi:10.1093/molbev/mst010.

## Figure reproduction

After completing the analysis pipeline, regenerate the publication figures with:

```bash
./run_figures.sh
```

Figure outputs are written to `figure_output/`.

The plotting scripts require Matplotlib v3.11.1.
