# Methods: HIV-1 NFLG Within-Host Sequence Variation Analysis

## Patient-specific sequence alignment

Near-full-length HIV-1 clone sequences were analyzed independently for each patient. Multiple-sequence alignments were performed using MAFFT v7.526. Patient grouping was determined from sequence identifiers, and no global all-patient multiple-sequence alignment was used for mutation definition.

HXB2 and supplied manual consensus sequences could be included as alignment guide sequences during alignment preparation, but neither was used as the mutation reference.

## Patient-specific consensus construction

An analytical consensus was reconstructed independently from the clone sequences of each patient.

At each alignment column, at least 90% of patient clones were required to contain a callable state. Among callable A, C, G, T, or gap states, a state was accepted as the consensus only when its frequency was strictly greater than 50%.

Sites failing the 90% callable-coverage criterion were classified as `Insufficient_callable_coverage`. Sites meeting the coverage criterion but lacking a state with >50% support were classified as `No_unique_consensus`.

The analytical consensus used for mutation calling was reconstructed from patient clone sequences rather than taken directly from a supplied manual consensus.

## Within-patient mutation calling

Clone-specific SNPs, insertions, and deletions were defined relative to the patient-specific analytical consensus.

Accordingly, a notation such as `A>G` represents a patient-consensus allele A and a clone allele G and does not imply an ancestral or evolutionary substitution direction.

Terminal alignment gaps attributable to incomplete sequence coverage were treated as missing/unsequenced sequence and were not called as deletions.

Clone `N` states representing known deletion padding were treated as gap alleles for deletion calling, whereas other genuine ambiguous nucleotide states were treated as missing.

Contiguous insertion or deletion calls belonging to the same clone were reconstructed into indel events. Consequently, summary DEL and INS counts represent event occurrences rather than numbers of affected nucleotides.

## Mapping-consensus construction

For coordinate projection, a patient-specific mapping consensus was constructed separately from mutation calling.

A bridge table retained the relationship between original patient multiple-sequence-alignment columns and positions in the mapping consensus. This step did not assign final HXB2 coordinates and did not recall SNPs, insertions, or deletions.

## HXB2 coordinate projection

Each patient mapping consensus was independently aligned to the HXB2 reference using the MAFFT `--addfragments` procedure.

HXB2 was used only as a standardized genomic coordinate and gene-annotation framework. Mutation identities remained defined relative to the patient-specific analytical consensus.

Clone-specific SNPs, insertions, and deletions were subsequently projected through the original-MSA-to-mapping-consensus bridge onto the final patient-specific HXB2 coordinate map.

Only non-gap HXB2 nucleotides incremented the canonical HXB2 coordinate. Positions aligned against gaps in HXB2 were represented using extended coordinates without renumbering downstream HXB2 bases.

Insertions were assigned using their corresponding HXB2 anchor or flanking coordinate interval. Exact extended insertion coordinates were retained in the detailed tables, whereas the common HXB2 anchor was used for cross-patient recurrence analyses.

Because gap placement can be locally uncertain in highly variable or complex-indel regions, HXB2 projection was interpreted as coordinate localization rather than as a redefinition of the underlying patient-consensus-relative mutation.

## Special structural regions

Regions lacking a unique patient consensus were retained as annotation tracks rather than converted into directional clone-versus-consensus mutations.

Similarly, large internal sequence regions identified as missing coverage were treated as missing-region annotations and were excluded from directional deletion counts.

For CN2023AH1-18, the major internal missing region corresponded to HXB2 positions 3487-5903 (2417 nt).

## Patient-level summaries

For patient-specific mutation profiles, SNPs were summarized at each HXB2-projected position according to substitution type, and the number of clones carrying each patient-consensus-relative substitution was used as the plotted count.

Deletion and insertion counts used for patient-level plotting represented contiguous event occurrences.

## Cross-patient recurrence analysis

Cross-patient recurrence was calculated only after independent patient-specific HXB2 projection.

The primary recurrence metric, `Patient_n`, was defined as the number of unique patients showing within-host variability at the same HXB2 locus or anchor.

For canonical HXB2 SNPs, recurrence was evaluated at the exact HXB2 position. SNPs corresponding to sequence aligned against an HXB2 gap were aggregated by the associated HXB2 anchor across patients rather than by assuming homology of patient-specific extended insertion coordinates.

Insertions were likewise compared using the common HXB2 anchor. Deletion events were compared using their projected HXB2 deletion loci.

No formal inferential statistical tests were applied to these mutation recurrence summaries.

## Software and quality control

Multiple-sequence alignment was performed with MAFFT v7.526. Data processing, mutation calling, HXB2 coordinate projection, aggregation, and summary-table generation were performed using custom Python scripts developed for this study with Python v3.14.6.

Consistency checks were implemented at each processing step to verify sequence integrity, mutation counts, indel-event reconstruction, and HXB2 coordinate projection.

## Reference

Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution. 2013;30(4):772-780. doi:10.1093/molbev/mst010.
