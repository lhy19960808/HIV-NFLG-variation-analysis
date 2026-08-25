# Reproducibility Guide

## Environment

Create the minimal environment:

```bash
conda env create -f environment.yml
conda activate hiv_nflg
```

The frozen analysis used:

```text
Python 3.14.6
MAFFT 7.526
```

The exact analysis environment is additionally recorded in:

- `docs/software_versions.txt`

## Verify frozen source code

Before running the analysis:

```bash
sha256sum -c tests/core_code.sha256
```

All nine production scripts should return `OK`.

## Input files

The default runner expects:

```text
data/All_185_clones_no_vector_no_primer_plus_consensus_renamed.fasta
data/all_patients_NFLG_with_HXB2_original.fasta
```

Alternative input locations can be provided through the `INPUT_FASTA` and `HXB2_SOURCE` environment variables.

## Complete reproduction

From the project root:

```bash
conda activate hiv_nflg
bash ./run_pipeline.sh
```

To use a custom output directory:

```bash
bash ./run_pipeline.sh /path/to/analysis_output
```

To set the number of MAFFT threads:

```bash
THREADS=8 \
bash ./run_pipeline.sh /path/to/analysis_output
```

The pipeline uses `set -euo pipefail` and stops when an analytical step returns a non-zero exit status.

## Expected frozen results

A successful reproduction of the frozen dataset should produce:

```text
Directional mutation calls       3287
SNP calls                         2477
DEL base calls                     592
INS base calls                     218
DEL events                          44
INS events                          31
Total INDEL events                  75

SNP loci/anchors                   973
Recurrent SNP loci                  52
Same-substitution recurrent         21
Mixed-substitution recurrent        31
DEL loci                             32
INS anchors                          13
```

The expected large internal missing region for `CN2023AH1-18` is:

```text
HXB2 interval: 3487-5903
Length: 2417 nt
```

## Release-level validation

After a complete run:

```bash
python tests/check_release_output.py \
    analysis_output
```

Expected result:

```text
Step05       PASS
Step06       PASS
Step07       PASS
Step08       PASS

RELEASE OUTPUT VALIDATION: PASS
```

## Interpretation boundary

The mutation reference and coordinate reference are deliberately distinct:

```text
Mutation identity:
clone versus patient-specific analytical consensus

Genomic localization:
HXB2 coordinate framework
```

Therefore HXB2 alleles must not be interpreted as the reference alleles used to define within-patient mutation direction.

## Development regression records



The public release-level validation script is `tests/check_release_output.py`.

## Frozen-code policy

The scientific algorithms and public schema of `code_clean/01-08` are frozen. Documentation, packaging, and release metadata may be updated without changing the frozen production scripts.
