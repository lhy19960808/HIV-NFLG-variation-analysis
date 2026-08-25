# Data Dictionary

## Scope and interpretation

This document describes the tabular outputs of the finalized HIV-1 near-full-length genome (NFLG) within-host variation pipeline. Mutation identity and coordinate localization are deliberately separated:

- **Mutation reference:** patient-specific analytical consensus.
- **Coordinate reference:** HXB2, used only for standardized genomic localization and gene annotation.
- A label such as `A>G` means **patient consensus A → clone G** and does not imply ancestral or evolutionary direction.
- Terminal alignment gaps from incomplete coverage are missing/unsequenced, not deletions.
- Clone `N` states designated as deletion padding are interpreted as gap/DEL alleles; other genuine IUPAC ambiguity is missing.
- DEL/INS summary counts represent contiguous clone-specific event occurrences, not numbers of affected bases.
- Cross-patient recurrence uses `Patient_n`, not clone count or event count, as the primary recurrence metric.
- Cross-patient insertions are aggregated by common HXB2 anchor; patient-specific `.insN` ranks are not assumed homologous across patients.
- No-unique-consensus regions and large missing regions are annotation tracks and are not added to directional DEL/INS counts.

## Common field glossary

Repeated fields use the definitions below unless a table description states otherwise.

- `Analysis_consensus` — Analytical consensus state reconstructed from patient clone sequences.
- `Analysis_consensus_state` — Analytical consensus state at the MSA column.
- `Callable_fraction` — Callable_n divided by total clone count.
- `Callable_n` — Number of clones with callable A/C/G/T/gap state at the column.
- `Check` — QC check identifier.
- `Clone` — Clone identifier within a patient.
- `Clone_allele` — Clone allele interpreted for the mutation call.
- `Clone_n` — Number of analyzed clones for the patient.
- `Consensus_allele` — Patient-specific analytical-consensus allele.
- `Coordinate_QC_status` — QC status of HXB2 coordinate mapping.
- `Coordinate_class` — Class of HXB2 coordinate (canonical base versus HXB2-gap/extended representation).
- `Coordinate_precision` — Whether coordinate assignment is exact or interval/other precision class.
- `Coordinate_source` — Source/projection logic used to assign the HXB2 coordinate.
- `Count_for_plot` — Patient-level plotted count: variant-clone count for SNPs and contiguous event-occurrence count for DEL/INS.
- `Display_interval` — HXB2 interval used to display a no-unique region as an annotation track.
- `Event_length_bp` — Length of the reconstructed contiguous indel event.
- `Event_length_bp_values` — Observed event lengths represented in the aggregate record.
- `Event_n` — Number of contiguous event occurrences contributing to a patient-level record.
- `Exact_coordinate_n` — Number of distinct exact coordinates represented.
- `Exact_coordinates_observed` — Exact projected coordinates represented within an anchor-level aggregate.
- `Expected` — Frozen/reference expected value.
- `Feature` — HXB2 feature name.
- `Final_HXB2_anchor` — Final canonical HXB2 anchor associated with the call.
- `Final_HXB2_coordinate` — Final HXB2 coordinate assigned after patient-specific Step 04B mapping.
- `Final_coordinate_class` — Final coordinate class after Step 05 projection.
- `Final_genes` — HXB2 gene/feature annotation(s) overlapping the projected coordinate.
- `Final_plot_coordinate` — Final numeric plotting coordinate.
- `Gap_n` — Number of callable gap states.
- `Gene` — HXB2 gene/feature annotation.
- `Genes` — One or more overlapping HXB2 gene/feature annotations.
- `HXB2_anchor` — Canonical HXB2 anchor used for HXB2-gap/inserted sequence.
- `HXB2_coordinate` — HXB2 coordinate assigned to the alignment/mapping position; extended coordinates may be used for sequence aligned to HXB2 gaps.
- `HXB2_end` — Last canonical HXB2 base in a mapped range/feature.
- `HXB2_locus_or_anchor` — Cross-patient aggregation key: exact HXB2 base for canonical SNPs, common HXB2 anchor for insertion-associated sequence, or projected deletion locus.
- `HXB2_role` — Explicit statement that HXB2 is used for coordinate/gene annotation rather than mutation definition.
- `HXB2_span` — Inclusive HXB2 span.
- `HXB2_start` — First canonical HXB2 base in a mapped range/feature.
- `MSA_col` — 1-based column in the patient-specific multiple-sequence alignment.
- `MSA_end` — Last MSA column of an event.
- `MSA_length` — Length of the patient-specific MSA in alignment columns.
- `MSA_start` — First MSA column of an event.
- `Major_fraction_of_callable` — Major_n divided by Callable_n.
- `Major_n` — Number of callable clones carrying Major_state.
- `Major_state` — Most frequent callable state at the column.
- `Missing_or_terminal_n` — Number of clone states treated as missing/unsequenced, including terminal alignment gaps.
- `Mutation_reference` — Explicit statement of the mutation reference; should indicate patient-specific consensus.
- `Mutation_type` — Mutation class, typically SNP, DEL, or INS.
- `No_unique_class` — Classification of a region without a unique analytical consensus.
- `Observed` — Observed value for the QC check.
- `Original_patient_MSA_column` — 1-based column in the original patient-specific MSA.
- `Patient` — Patient identifier.
- `Patient_n` — Number of unique patients showing the relevant within-host variable state at the locus/anchor; primary cross-patient recurrence metric.
- `Patients` — Delimited list of contributing patient identifiers.
- `Plot_coordinate` — Numeric x-coordinate used for plotting. Extended insertion coordinates may be placed between adjacent canonical HXB2 bases.
- `Projection_status` — Status assigned during HXB2 projection.
- `Recurrent_class` — Recurrence category based on Patient_n and mutation/substitution pattern.
- `Role` — Role of the listed source/output file.
- `SNP_type` — Patient-consensus-relative substitution label (for example A>G).
- `Same_substitution_across_all_recurrent_patients` — Whether all recurrent patients share the same patient-consensus-relative substitution.
- `Source_or_output` — Path recorded in the provenance manifest.
- `Status` — QC/provenance status.
- `Substitution_patient_counts` — Per-substitution patient counts at the locus/anchor.
- `Substitution_pattern` — Summary of substitution types observed among contributing patients.
- `Total_event_n` — Total contiguous DEL/INS event occurrences contributing to an aggregated record.
- `Total_variant_clone_n` — Total number of distinct variant clones contributing to the aggregated record.
- `Variant_clone_n` — Number of distinct variant clones contributing to a patient-level record.

## Output tables

### Step 02 — Patient-specific analysis consensus

#### `01_analysis_consensus_summary.tsv`

Per-patient summary of clone count, MSA length, supplied consensus N-resolution, and analytical-consensus state counts.

**Columns:** `Patient`, `Clone_n`, `MSA_length`, `SnapGene_consensus_name`, `SnapGene_N_total_n`, `SnapGene_N_resolved_to_base_n`, `SnapGene_N_resolved_to_gap_n`, `SnapGene_N_remained_no_unique_n`, `SnapGene_N_remained_low_coverage_n`, `Analysis_consensus_base_column_n`, `Analysis_consensus_gap_column_n`, `Analysis_consensus_no_unique_N_column_n`, `Analysis_consensus_low_coverage_N_column_n`, `Analysis_consensus_zero_callable_N_column_n`.

#### `02_snapgene_N_resolution.tsv`

Column-level audit of positions where the supplied manual consensus contains N and how clone evidence resolves each position.

**Columns:** `Patient`, `MSA_col`, `HXB2_coordinate`, `Column_type`, `Gene`, `SnapGene_consensus`, `Clone_n`, `A_n`, `C_n`, `G_n`, `T_n`, `Gap_n`, `Missing_or_terminal_n`, `Callable_n`, `Callable_fraction`, `Major_state`, `Major_n`, `Major_fraction_of_all_clones`, `Major_fraction_of_callable`, `Analysis_consensus`, `Reason`.

#### `03_snapgene_vs_analysis_consensus_differences.tsv`

Column-level differences between the supplied manual consensus and the reconstructed analytical consensus.

**Columns:** `Patient`, `MSA_col`, `HXB2_coordinate`, `Column_type`, `Gene`, `SnapGene_consensus`, `Analysis_consensus`, `Reason`, `Clone_n`, `A_n`, `C_n`, `G_n`, `T_n`, `Gap_n`, `Missing_or_terminal_n`, `Callable_n`, `Callable_fraction`, `Major_state`, `Major_n`, `Major_fraction_of_all_clones`, `Major_fraction_of_callable`.

#### `04_no_unique_consensus_sites.tsv`

Alignment columns that pass callable coverage but lack a unique >50% consensus state.

**Columns:** `Patient`, `MSA_col`, `HXB2_coordinate`, `Column_type`, `Gene`, `Clone_n`, `A_n`, `C_n`, `G_n`, `T_n`, `Gap_n`, `Missing_or_terminal_n`, `Callable_n`, `Callable_fraction`, `Major_state`, `Major_n`, `Major_fraction_of_all_clones`, `Major_fraction_of_callable`, `No_unique_type`.

#### `05_insufficient_callable_coverage_sites.tsv`

Alignment columns that fail the minimum callable-coverage requirement.

**Columns:** `Patient`, `MSA_col`, `HXB2_coordinate`, `Column_type`, `Gene`, `Clone_n`, `A_n`, `C_n`, `G_n`, `T_n`, `Gap_n`, `Missing_or_terminal_n`, `Callable_n`, `Callable_fraction`, `Major_state`, `Major_n`, `Major_fraction_of_all_clones`, `Major_fraction_of_callable`, `Analysis_consensus`, `Reason`.


### Step 03 — Within-patient mutation calling

#### `01_clone_mutation_calls.tsv`

Initial clone-versus-patient-consensus SNP/DEL/INS calls in patient MSA coordinates before final HXB2 projection.

**Columns:** `Patient`, `Clone`, `MSA_col`, `HXB2_coordinate`, `Column_type`, `HXB2_gene_mapping_position`, `HXB2_insertion_anchor`, `Gene`, `Analysis_consensus_state`, `Clone_raw_state`, `Clone_interpreted_state`, `Clone_state_source`, `Mutation_type`, `Change`.

#### `02_indel_events.tsv`

Contiguous clone-specific insertion and deletion events reconstructed from base/column-level calls.

**Columns:** `Patient`, `Clone`, `Event_type`, `Event_index_within_clone`, `MSA_start`, `MSA_end`, `HXB2_start_coordinate`, `HXB2_end_coordinate`, `Genes`, `Event_length_alignment_columns`, `HXB2_base_column_n`, `HXB2_relative_insertion_column_n`, `Raw_clone_N_as_DEL_column_n`, `Event_locus_key`.

#### `03_patient_mutation_summary.tsv`

Compact patient-level SNP and indel event summary used in plotting.

**Columns:** `Patient`, `Clone_n`, `SNP_clone_call_n_all_alignment_columns`, `SNP_unique_site_n_all_alignment_columns`, `SNP_clone_call_n_HXB2_base_only`, `SNP_unique_site_n_HXB2_base_only`, `DEL_base_call_n`, `DEL_event_occurrence_n`, `DEL_unique_event_locus_n`, `INS_base_call_n`, `INS_event_occurrence_n`, `INS_unique_event_locus_n`, `Clone_N_as_DEL_base_call_n`, `Clone_N_as_DEL_event_occurrence_n`, `Consensus_no_unique_column_n`, `Consensus_low_coverage_column_n`, `Consensus_zero_callable_column_n`.

#### `04_HXB2_base_position_mutation_counts.tsv`

Exploratory mutation counts at HXB2 base positions from the Step 03 coordinate representation; final HXB2 projection occurs later in Step 05.

**Columns:** `HXB2_position`, `Gene`, `SNP_clone_call_n`, `SNP_variant_clone_n`, `SNP_variant_patient_n`, `DEL_clone_call_n`, `DEL_variant_clone_n`, `DEL_variant_patient_n`, `INS_clone_call_n`, `INS_variant_clone_n`, `INS_variant_patient_n`.

#### `05_HXB2_relative_insertion_event_counts.tsv`

Exploratory insertion-event counts by HXB2-relative anchor from Step 03; final cross-patient recurrence is prepared later.

**Columns:** `HXB2_insertion_anchor`, `Plot_x`, `Gene`, `INS_event_occurrence_n`, `INS_variant_clone_n`, `INS_variant_patient_n`, `INS_total_inserted_base_n`, `INS_min_event_length`, `INS_max_event_length`.

#### `06_clone_N_as_DEL_QC.tsv`

QC audit confirming how clone N states designated as deletion padding were interpreted.

**Columns:** `Patient`, `Clone`, `Raw_clone_N_total_n`, `Raw_clone_N_called_as_DEL_n`, `QC`.

#### `07_mutation_calling_QC.tsv`

Per-patient alignment/consensus length consistency checks for mutation calling.

**Columns:** `Patient`, `Clone_n`, `MSA_length`, `Alignment_and_consensus_length_QC`.


### Step 04A — Mapping-consensus bridge

#### `01_mapping_consensus_summary.tsv`

Per-patient summary of mapping-consensus construction and excluded/retained MSA-column classes.

**Columns:** `Patient`, `Clone_n`, `Original_patient_MSA_length`, `Mapping_consensus_length`, `Analysis_consensus_base_columns`, `No_unique_base_only_columns_retained_as_N`, `Consensus_gap_columns_excluded`, `No_unique_structural_columns_excluded`, `Low_coverage_columns_excluded`, `No_callable_columns_excluded`.

#### `02_original_MSA_column_consensus_index.tsv`

Bridge from original patient MSA columns to analytical-consensus states and mapping-consensus positions.

**Columns:** `Patient`, `Original_patient_MSA_column`, `Analysis_consensus_state`, `Consensus_reason`, `Callable_clone_n`, `Total_clone_n`, `Callable_fraction`, `Callable_state_counts`, `Mapping_consensus_position`, `Mapping_inclusion_class`.


### Step 04B — Patient-specific HXB2 mapping

#### `01_mapping_consensus_to_HXB2.tsv`

Position-level map from patient mapping consensus to the independently aligned HXB2 coordinate framework.

**Columns:** `Patient`, `Mapping_consensus_position`, `Pairwise_alignment_column`, `Mapping_consensus_base`, `HXB2_state`, `HXB2_coordinate`, `Coordinate_class`, `HXB2_anchor`, `Plot_coordinate`.

#### `02_patient_HXB2_range_summary.tsv`

Per-patient HXB2 mapped range and alignment identity summary.

**Columns:** `Patient`, `Mapping_consensus_length`, `HXB2_start`, `HXB2_end`, `HXB2_span`, `HXB2_gap_columns_n`, `Query_bases_at_HXB2_gap_n`, `HXB2_bases_at_query_gap_n`, `Identity_ACGT_only`.

#### `07_addfragments_alignment_QC.tsv`

QC metrics for MAFFT --addfragments alignment and sequence integrity.

**Columns:** `Patient`, `Mapping_consensus_length`, `Addfragments_alignment_length`, `HXB2_gap_columns_n`, `Query_gap_columns_n`, `HXB2_start`, `HXB2_end`, `HXB2_span`, `Aligned_ACGT_pair_n`, `Exact_match_n`, `Identity_ACGT_only`, `Query_bases_at_HXB2_gap_n`, `HXB2_bases_at_query_gap_n`, `HXB2_sequence_integrity`, `Query_sequence_integrity`.


### Step 05 — Final HXB2 mutation projection

#### `01_final_clone_mutation_calls.tsv`

Authoritative clone-level SNP/DEL/INS calls after projection onto the final HXB2 coordinate framework.

**Columns:** `Patient`, `Clone`, `Original_patient_MSA_column`, `Mutation_type`, `SNP_type`, `Consensus_allele`, `Clone_allele`, `Final_HXB2_coordinate`, `Final_coordinate_class`, `Final_HXB2_anchor`, `Final_plot_coordinate`, `Final_genes`, `Coordinate_source`, `Projection_status`, `Analysis_consensus_state_at_MSA_column`, `Analysis_consensus_reason`, `Mapping_inclusion_class`, `Gap_block_ID`, `Gap_block_start_MSA_column`, `Gap_block_end_MSA_column`, `Gap_block_length_alignment_columns`, `Insertion_alignment_rank_in_gap_block`, `Left_flank_MSA_column`, `Left_flank_HXB2_coordinate`, `Right_flank_MSA_column`, `Right_flank_HXB2_coordinate`, `Coordinate_precision`.

#### `02_final_indel_events.tsv`

Authoritative contiguous clone-specific DEL/INS events after HXB2 projection.

**Columns:** `Patient`, `Clone`, `Event_ID`, `Event_type`, `Start_original_MSA_column`, `End_original_MSA_column`, `Event_length_bp`, `Start_HXB2_coordinate`, `End_HXB2_coordinate`, `Event_HXB2_locus`, `Plot_coordinate`, `Genes`, `Projection_status`, `Constituent_projection_statuses`.

#### `03_projection_QC.tsv`

Global projection checks comparing observed values with frozen expected values.

**Columns:** `Check`, `Observed`, `Expected`, `Status`.

#### `04_insertion_block_QC.tsv`

QC for consensus-gap blocks used to assign exact HXB2 insertion anchors.

**Columns:** `Patient`, `Gap_block_ID`, `Gap_block_start_MSA_column`, `Gap_block_end_MSA_column`, `Gap_block_length_alignment_columns`, `Left_flank_MSA_column`, `Left_flank_HXB2_coordinate`, `Right_flank_MSA_column`, `Right_flank_HXB2_coordinate`, `Insertion_HXB2_anchor`, `Insertion_plot_coordinate`, `Coordinate_precision`, `Block_projection_status`.

#### `05_review_or_unmapped_calls.tsv`

Calls not passing final projection; expected to be empty in the frozen release.

**Columns:** `Patient`, `Clone`, `Original_patient_MSA_column`, `Mutation_type`, `SNP_type`, `Consensus_allele`, `Clone_allele`, `Final_HXB2_coordinate`, `Final_coordinate_class`, `Final_HXB2_anchor`, `Final_plot_coordinate`, `Final_genes`, `Coordinate_source`, `Projection_status`, `Analysis_consensus_state_at_MSA_column`, `Analysis_consensus_reason`, `Mapping_inclusion_class`, `Gap_block_ID`, `Gap_block_start_MSA_column`, `Gap_block_end_MSA_column`, `Gap_block_length_alignment_columns`, `Insertion_alignment_rank_in_gap_block`, `Left_flank_MSA_column`, `Left_flank_HXB2_coordinate`, `Right_flank_MSA_column`, `Right_flank_HXB2_coordinate`, `Coordinate_precision`.

#### `06_final_patient_HXB2_ranges.tsv`

Final patient-specific HXB2 coordinate ranges used downstream.

**Columns:** `Patient`, `Final_HXB2_start`, `Final_HXB2_end`, `Final_HXB2_span`, `Coordinate_framework`.

#### `07_no_unique_consensus_regions.tsv`

Collapsed no-unique-consensus regions for structural/figure annotation, not directional mutation counts.

**Columns:** `Patient`, `No_unique_class`, `Start_original_MSA_column`, `End_original_MSA_column`, `Length_alignment_columns`, `Left_anchor_HXB2_coordinate`, `Right_anchor_HXB2_coordinate`, `Display_interval`.


### Step 06 — Patient-level plotting tables

#### `01_patient_plot_ready.tsv`

Patient-level plot-ready SNP/DEL/INS records with plotting counts and display metadata.

**Columns:** `Patient`, `Record_type`, `Mutation_type`, `SNP_type`, `Exact_HXB2_coordinate`, `HXB2_anchor`, `Plot_coordinate`, `Coordinate_class`, `Count_for_plot`, `Variant_clone_n`, `Event_n`, `Genes`, `Display_symbol`, `Event_length_bp_values`, `Event_IDs`, `Projection_status`.

#### `02_all_patient_plot_ready.tsv`

Combined across-patient plot-ready aggregation with patient, clone, event, and coordinate counts.

**Columns:** `Mutation_type`, `SNP_type`, `HXB2_locus_or_anchor`, `Plot_coordinate`, `Patient_n`, `Patients`, `Total_variant_clone_n`, `Total_event_n`, `Exact_coordinate_n`, `Exact_coordinates`, `Genes`, `Display_symbol`, `Projection_status`.

#### `03_patient_mutation_summary.tsv`

Compact patient-level SNP and indel event summary used in plotting.

**Columns:** `Patient`, `SNP_call_n`, `DEL_event_n`, `INS_event_n`, `Total_indel_event_n`.

#### `04_plot_data_QC.tsv`

QC checks for patient-level plotting tables.

**Columns:** `Check`, `Observed`, `Expected`, `Status`.

#### `05_patient_panel_metadata.tsv`

Per-patient panel metadata including HXB2 span and presence of directional mutation records.

**Columns:** `Patient`, `HXB2_start`, `HXB2_end`, `HXB2_span`, `Has_directional_mutation_records`, `No_unique_region_n`.

#### `06_no_unique_consensus_regions.tsv`

No-unique-consensus annotation regions propagated for patient plotting.

**Columns:** `Patient`, `No_unique_class`, `Start_original_MSA_column`, `End_original_MSA_column`, `Length_alignment_columns`, `Left_anchor_HXB2_coordinate`, `Right_anchor_HXB2_coordinate`, `Display_interval`.


### Step 07 — Cross-patient dataset

#### `01_SNP_locus_recurrence.tsv`

Cross-patient SNP recurrence by exact HXB2 base or HXB2 anchor, including substitution-specific patient/clone counts.

**Columns:** `HXB2_locus_or_anchor`, `Locus_class`, `Plot_coordinate`, `Patient_n`, `Patients`, `Total_variant_clone_n`, `Genes`, `Exact_coordinate_n`, `Exact_coordinates_observed`, `Substitution_type_n`, `A>C_patient_n`, `A>C_clone_n`, `A>G_patient_n`, `A>G_clone_n`, `A>T_patient_n`, `A>T_clone_n`, `C>A_patient_n`, `C>A_clone_n`, `C>G_patient_n`, `C>G_clone_n`, `C>T_patient_n`, `C>T_clone_n`, `G>A_patient_n`, `G>A_clone_n`, `G>C_patient_n`, `G>C_clone_n`, `G>T_patient_n`, `G>T_clone_n`, `T>A_patient_n`, `T>A_clone_n`, `T>C_patient_n`, `T>C_clone_n`, `T>G_patient_n`, `T>G_clone_n`, `Aggregation_note`.

#### `02_INDEL_locus_recurrence.tsv`

Cross-patient DEL/INS recurrence by projected deletion locus or common HXB2 insertion anchor.

**Columns:** `Mutation_type`, `HXB2_locus_or_anchor`, `Aggregation_class`, `Plot_coordinate`, `Patient_n`, `Patients`, `Total_event_n`, `Total_variant_clone_n`, `Event_length_bp_values`, `Genes`, `Aggregation_note`.

#### `03_combined_plot_ready.tsv`

Unified SNP/DEL/INS cross-patient table prepared for Figure 2-style visualization.

**Columns:** `Record_type`, `Mutation_type`, `HXB2_locus_or_anchor`, `Locus_class`, `Plot_coordinate`, `Patient_n`, `Patients`, `Total_variant_clone_n`, `Total_event_n`, `Genes`, `Primary_plot_metric`, `Display_recommendation`, `Aggregation_note`.

#### `04_patient_summary.tsv`

Final per-patient summary integrating clone count, mapped range, mutations, no-unique regions, and special annotations.

**Columns:** `Patient`, `Clone_n`, `Final_HXB2_start`, `Final_HXB2_end`, `Final_HXB2_span`, `Directional_SNP_call_n`, `Patient_SNP_locus_type_n`, `DEL_event_n`, `INS_event_n`, `Total_INDEL_event_n`, `No_unique_consensus_region_n`, `Special_structural_annotation`.

#### `05_HXB2_mapping_summary.tsv`

Publication-facing summary of coordinate framework, mutation reference, HXB2 role, and QC status.

**Columns:** `Patient`, `HXB2_start`, `HXB2_end`, `HXB2_span`, `Coordinate_framework`, `Mutation_reference`, `HXB2_role`, `Coordinate_QC_status`.

#### `06_HXB2_gene_coordinates.tsv`

HXB2 feature coordinates copied from the release configuration.

**Columns:** `Feature`, `HXB2_start`, `HXB2_end`, `Feature_type`.

#### `07_special_region_summary.tsv`

Special structural/missing-region annotations kept separate from directional mutation counts.

**Columns:** `Patient`, `Annotation_type`, `HXB2_start`, `HXB2_end`, `Length_nt_or_alignment_columns`, `Coordinate_precision`, `Evidence`, `Counted_as_clone_vs_consensus_mutation`, `Figure_annotation`.

#### `08_no_unique_consensus_regions.tsv`

No-unique-consensus regions propagated into the cross-patient dataset.

**Columns:** `Patient`, `No_unique_class`, `Start_original_MSA_column`, `End_original_MSA_column`, `Length_alignment_columns`, `Left_anchor_HXB2_coordinate`, `Right_anchor_HXB2_coordinate`, `Display_interval`.

#### `09_large_missing_region_evidence.tsv`

Detailed coordinate evidence supporting large internal missing-region annotations.

**Columns:** `Patient`, `Annotation_type`, `Original_patient_MSA_start`, `Original_patient_MSA_end`, `Original_patient_MSA_length`, `Left_direct_MSA_column`, `Left_mapping_consensus_position`, `Left_HXB2_flank`, `Right_direct_MSA_column`, `Right_mapping_consensus_position`, `Right_HXB2_flank`, `HXB2_missing_start`, `HXB2_missing_end`, `HXB2_missing_length`, `Inferred_HXB2_missing_interval`, `Length_check`, `Coordinate_method`, `Interpretation`.

#### `10_source_file_manifest.tsv`

Provenance manifest of source and output paths. Run-directory prefixes can differ between reproductions without changing scientific results.

**Columns:** `Source_or_output`, `Role`, `Status`.

#### `11_final_data_QC.tsv`

Final cross-patient dataset QC and frozen expected-value checks.

**Columns:** `Check`, `Observed`, `Expected`, `Status`.


### Step 08 — Figure 2 plotting data

#### `01_figure2_plotting_summary.tsv`

Primary final plotting table combining SNP, DEL, and INS locus/anchor summaries and display annotations.

**Columns:** `Record_type`, `Mutation_type`, `HXB2_locus_or_anchor`, `Locus_class`, `Plot_coordinate`, `Patient_n`, `Patients`, `Total_variant_clone_n`, `Total_event_n`, `Genes`, `Recurrent_class`, `Substitution_pattern`, `Same_substitution_across_all_recurrent_patients`, `Substitution_patient_counts`, `Exact_coordinates_observed`, `Display_priority`, `Figure2_role`, `Interpretation`.

#### `02_recurrent_SNP_sites.tsv`

Subset of recurrent SNP loci used for recurrence-focused Figure 2 interpretation.

**Columns:** `Record_type`, `Mutation_type`, `HXB2_locus_or_anchor`, `Locus_class`, `Plot_coordinate`, `Patient_n`, `Patients`, `Total_variant_clone_n`, `Total_event_n`, `Genes`, `Recurrent_class`, `Substitution_pattern`, `Same_substitution_across_all_recurrent_patients`, `Substitution_patient_counts`, `Exact_coordinates_observed`, `Display_priority`, `Figure2_role`, `Interpretation`.

#### `03_SNP_recurrence_summary.tsv`

Compact category/count summary of SNP recurrence classes.

**Columns:** `Category`, `Count`, `Definition`.

#### `04_INDEL_plotting_summary.tsv`

Final DEL/INS plotting summary with patient/event/clone counts and interpretation metadata.

**Columns:** `Record_type`, `Mutation_type`, `HXB2_locus_or_anchor`, `Locus_class`, `Plot_coordinate`, `Patient_n`, `Patients`, `Total_variant_clone_n`, `Total_event_n`, `Genes`, `Recurrent_class`, `Substitution_pattern`, `Same_substitution_across_all_recurrent_patients`, `Substitution_patient_counts`, `Exact_coordinates_observed`, `Display_priority`, `Figure2_role`, `Interpretation`.

#### `05_figure2_plot_QC.tsv`

Final Figure 2 plotting-data QC against frozen expected values.

**Columns:** `Check`, `Observed`, `Expected`, `Status`.

## Counting conventions

- **SNP call counts** count patient-consensus-relative clone calls.
- **Variant clone counts** count distinct clones carrying the relevant variant.
- **DEL/INS event counts** count contiguous reconstructed event occurrences.
- **Patient_n** counts unique patients contributing the locus/anchor and is the primary Figure 2 recurrence metric.
- Overlapping HXB2 genes/features can cause one coordinate or event to be annotated to more than one feature; gene-level counts therefore need not sum to whole-genome totals.

## Coordinate conventions

Canonical HXB2 numbering advances only across non-gap HXB2 bases. Sequence aligned to an HXB2 gap can receive an extended coordinate between adjacent canonical HXB2 bases. Exact extended coordinates are retained in detailed tables, while anchor-level cross-patient analyses use the corresponding common HXB2 anchor.

## Special structural annotations

The release retains no-unique-consensus regions and large internal missing regions separately from directional mutation calls. In the frozen dataset, the major internal missing region for `CN2023AH1-18` is localized to HXB2 3487–5903 (2417 nt) and is not counted as a clone-versus-consensus deletion event.
