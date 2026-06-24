# Detailed Methods Draft for DSM123 Metabolic Model Reconstruction and Curation

## 1. Input datasets and reference resources

We reconstructed the DSM123 genome-scale metabolic model using a hybrid workflow that combines orthology-guided transfer, flux-informed gap recovery, and iterative manual curation.

Primary inputs:
- Reference model: `iDT1294Photo.xml` (converted to `iDT1294Photo.json` when needed).
- Target genome: `genomes/DSM123.gb`.
- Reference/target protein and nucleotide FASTA files generated from GenBank CDS features.
- Bidirectional BLAST outputs in `bbh/` and `nucl/`.

All reconstruction scripts were run from the project root (`workflow1`).

## 2. Automated draft reconstruction pipeline

The automated pipeline is implemented in `src/dsm123_pipeline/pipeline.py` and can be invoked through `scripts/run_pipeline.py`.

### 2.1 Genome parsing

CDS entries were parsed from target GenBank, and both protein and nucleotide FASTA files were generated.
- Protein FASTA: `prots/DSM123.fa`
- Nucleotide FASTA: `nucl/DSM123.fa`

Gene identifiers were assigned by priority:
1. `locus_tag`
2. `gene`
3. fallback synthetic ID when both are absent.

### 2.2 Orthology inference (bidirectional BLAST)

Bidirectional BLASTP was performed between reference and target proteomes.
Reciprocal best-hit logic was used to construct a BBH table.

Thresholds used in pipeline config (`src/dsm123_pipeline/config.py`):
- BLAST E-value cutoff: `1e-3`
- BBH coverage threshold: `0.2`
- Ortholog PID threshold for binary orthology matrix: `65.0`

Outputs:
- Parsed BBH table: `bbh/Rpal_BisA53_vs_DSM123_parsed.csv`
- Orthology matrix: `ortho_matrix.csv`
- Gene ID mapping matrix: `geneIDs_matrix.csv`

### 2.3 Unannotated ORF rescue using nucleotide alignment

Reference-vs-target BLASTN was run using target genomic FASTA built from GenBank.
Candidate missing genes were screened with:
- BLASTN PID > `70.0`
- Alignment length ratio > `0.8`
- Unannotated rescue threshold PID >= `80.0`

Candidates satisfying criteria were added as ortholog-supported entries.

### 2.4 Draft model assembly and gene renaming

The reference JSON model was copied, non-homologous genes were removed, and corresponding unsupported reactions were deleted (`cobra.manipulation.delete.remove_genes(..., remove_reactions=True)`).
Remaining genes were renamed to target IDs using `geneIDs_matrix.csv` mapping.

Draft output:
- `Models/DSM123.json` (or target-specific draft path in config)

### 2.5 Flux-informed reaction recovery (gap recovery)

To reduce false negatives introduced by strict orthology filtering, reactions with non-zero reference flux/reduced-cost support were scanned and missing reactions were reintroduced into the target draft.

Key outputs:
- `missing_reactions.csv`
- `reaction_gene_relationships.csv`
- `updated_consensus.json`

### 2.6 FBA checkpoints

FBA was computed under default and M9-like medium settings.
Outputs:
- `fba_final_default_fluxes.csv`
- `fba_final_m9_fluxes.csv`
- `fba_summary.csv`

## 3. Step1 manual curation: gene-reaction consistency audit

Step1 artifacts were regenerated using:
- `python scripts/regenerate_step1_step2_reviews.py`

The script reconstructs a pre-step2 model and compares model-assigned gene-reaction links against BBH-supported reference mapping.

Decision-support tables:
- `manual_curation_outputs_merged/step1_review/model_gpr_conflicts_report.csv`
- `manual_curation_outputs_merged/step1_review/model_missing_genes_candidates.csv`

Step1 scoring thresholds for candidate filtering:
- Identity >= `50%`
- E-value <= `1e-3`

These tables are intended for curator action annotation before subsequent rewrites.

## 4. Step2 manual curation: metabolite namespace + structure cleanup

Step2 uses the same script and writes into:
- `manual_curation_outputs_merged/step2_review/`

### 4.1 `_u` suffix normalization

Metabolites ending with `_u` were normalized:
- If base ID existed: merge `_u` metabolite into existing target and rewire stoichiometry.
- If base ID did not exist: rename by dropping `_u` and refresh compartment suffix.

### 4.2 Quinolinate ID unification

`quln_c` was set as the canonical quinolinate ID.
Alias IDs (if present) were merged into `quln_c`, and NNDPR-like reaction stoichiometry was rewritten to canonical participants when all required metabolites were available.

### 4.3 Duplicate metabolite candidate detection (manual review table)

Duplicate candidates were grouped strictly by:
- same compartment
- same molecular formula

The candidate table:
- `step2_table_duplicate_metabolite_candidates.csv`

Fields include IDs and names, and curator-provided `new` instructions.

Manual syntax for merge plan supports:
- `a | b = target`
- `a | b | c = target`
- Single target token (interpreted as merge all listed IDs to that target).

Blank `new` means no rewrite.

### 4.4 Applying duplicate merge plan

For each reviewed row with non-empty `new`, source metabolite occurrences in all reactions were transferred to target metabolite, then source nodes were removed.

### 4.5 Compartment bridge audit and insertion

Bridge reactions were generated when metabolite base tokens existed in compartment pairs but lacked shared reaction connectivity.
Audited pairs:
- `c <-> p` (`cp`)
- `e <-> p` (`ep`)
- `c <-> e` (`ce`)

Bridge output table:
- `step2_table_compartment_bridge_reactions.csv`

### 4.6 Orphan metabolite table and action-driven deletion

A post-merge orphan scan generated:
- `step2_table_orphan_metabolites.csv`

Rows with `action = no` were removed from model; other rows were retained.
Empty reactions generated by metabolite removal were dropped.

### 4.7 Step2 outputs retained by design

Per curation request, only core review tables are preserved:
- `step2_table_duplicate_metabolite_candidates.csv`
- `step2_table_compartment_bridge_reactions.csv`
- `step2_table_orphan_metabolites.csv`

Final step2-updated working model:
- `Models/DSM123_manual_working.json`

## 5. Final model used for downstream summary

For reporting and figure-data extraction, the curated final model analyzed here is:
- `Models/purple_bacteriav_DSM123.json`

### 5.1 Operational definition of manual curation additions/changes

In this study, all differences between `updated_consensus.json` and `Models/purple_bacteriav_DSM123.json` were treated as manually reviewed and manually applied curation outcomes.

Specifically, the following events were all classified as manual curation edits:
- reaction additions
- reaction removals
- reaction-level GPR rewrites
- gene additions/removals associated with those curated reaction/GPR updates

This definition was used consistently in downstream figure generation and change accounting.

## 6. Stage-wise quantitative tracking for figure generation

To support bar plots and delta analysis, we exported standardized tables using:
- `python scripts/export_purple_stage_stats.py`

Stages analyzed:
1. `iDT1294Photo.json`
2. `Models/DSM123.json`
3. `updated_consensus.json`
4. `Models/purple_bacteriav_DSM123.json`

For each stage we report:
- number of genes
- number of metabolites
- number of reactions
- number of exogenous genes

Exogenous genes are defined as model genes not found among CDS locus tags in `genomes/DSM123.gb`.

Output:
- `manual_curation_outputs_merged/figure_data/plot_table_stage_4d_counts.csv`

## 7. Final-model enrichment summaries for pie/bar charts

For the final model (`purple_bacteriav_DSM123.json`) we exported:

1. Reaction type enrichment (counts and proportions)
- `plot_table_final_reaction_type_enrichment.csv`

2. Reaction subsystem enrichment (counts and proportions)
- `plot_table_final_reaction_subsystem_enrichment.csv`

3. Gene-linked enrichment by reaction type
- `plot_table_final_gene_enrichment_by_reaction_type.csv`

4. Gene-linked enrichment by subsystem
- `plot_table_final_gene_enrichment_by_subsystem.csv`

Gene-linked enrichment is based on parsed gene-reaction-rule (GPR) associations (gene-reaction pairs and unique-gene coverage per category).

### 7.1 Subsystem enrichment of reactions (what was actually calculated)

The reaction-level subsystem analysis is a composition-based enrichment summary (category frequency), not a hypothesis-test enrichment against an external background.

Procedure:
1. Parse every reaction in `purple_bacteriav_DSM123.json`.
2. Read `reaction["subsystem"]`.
3. If subsystem is missing/empty, assign category `Unknown`.
4. Count reactions per subsystem:
   - `count(subsystem_i) = number of reactions annotated as subsystem_i`
5. Compute proportions:
   - `fraction(subsystem_i) = count(subsystem_i) / total_reactions`
   - `percent(subsystem_i) = 100 * fraction(subsystem_i)`
6. Sort by descending count and export:
   - `plot_table_final_reaction_subsystem_enrichment.csv`

Therefore, this table should be interpreted as subsystem composition of the final model, suitable for pie/bar plotting.

### 7.2 Gene-linked subsystem enrichment (GPR-based)

For each reaction with non-empty GPR:
1. Parse gene tokens from `gene_reaction_rule` (logical words `and/or/not` removed).
2. Keep only parsed genes present in the model gene list.
3. For each valid gene associated with a reaction, record one gene-reaction pair under that reaction's subsystem.
4. For each subsystem, export:
   - `gene_reaction_pair_count`: total number of gene-reaction pairs in that subsystem
   - `pair_fraction`, `pair_percent`: proportion of pairs among all gene-reaction pairs
   - `unique_gene_count`: unique genes linked to that subsystem
   - `unique_gene_fraction_of_model`, `unique_gene_percent_of_model`: coverage in all model genes
   - `reaction_count_with_gene_rule`: number of reactions in that subsystem with non-empty GPR

Output file:
- `plot_table_final_gene_enrichment_by_subsystem.csv`

### 7.3 Important interpretation note

In the current final model, most reactions have missing subsystem annotation and are therefore labeled `Unknown`.
As a result, subsystem composition is dominated by `Unknown`, which reflects annotation completeness rather than pathway absence.

If statistical over-representation analysis is needed later (e.g., hypergeometric/Fisher against a defined reaction universe), that should be implemented as a separate analysis with an explicit background set.

## 8. Delta reporting between `updated_consensus` and final curated model

To document reconstruction-to-curation changes, we exported:

1. Added genes with linked reactions and GPR rules
- `plot_table_updated_to_final_added_genes_with_gpr.csv`

2. Reaction-level GPR deltas
- `plot_table_updated_to_final_gpr_changes.csv`
- Change types: `added_reaction`, `removed_reaction`, `gpr_changed`

3. One-row overall delta summary
- `plot_table_updated_to_final_overall_delta.csv`

Under the operational definition in Section 5.1, these delta tables represent the manually curated change set from `updated_consensus` to `purple_bacteriav_DSM123`.

## 9. Reproducibility and environment notes

The workflow relies on:
- COBRApy for model IO/manipulation/FBA
- Biopython (`SeqIO`) for genome parsing
- pandas for tabular audit outputs
- BLAST+ (`blastp`, `blastn`, `makeblastdb`) for homology inference

Main reproducibility scripts:
- `scripts/run_pipeline.py`
- `scripts/regenerate_step1_step2_reviews.py`
- `scripts/run_step3_checkpoint.py` (optional standalone step3 audit)
- `scripts/export_purple_stage_stats.py`

This Methods draft is intentionally verbose so sections can be trimmed during manuscript preparation.
