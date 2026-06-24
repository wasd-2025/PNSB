# Workflow 1: DSM123 Metabolic Model Reconstruction

This directory contains the complete reconstruction, manual curation, validation, and downstream modeling workflow for the *Rhodopseudomonas palustris* DSM123 genome-scale metabolic model.

The workflow starts from the photoheterotrophic reference model `iDT1294Photo`, transfers gene-reaction content using sequence homology, recovers flux-supported reactions, applies manual metabolite and GPR curation, and produces the final DSM123 model used by the analyses in `plot/`.

## Workflow overview

```text
Reference model + DSM123 genome
        |
        v
Bidirectional BLAST and orthology mapping
        |
        v
DSM123 draft model
        |
        v
Flux-supported reaction recovery
        |
        v
Step 1: GPR consistency review
        |
        v
Step 2: metabolite and compartment curation
        |
        v
Final model: Models/purple_bacteriav_DSM123.json
        |
        +--> model-stage and enrichment summaries
        `--> plot/: dFBA, FVA, kinetic sensitivity, and heatmaps
```

## Directory structure

```text
workflow_1_model_reconstruction/
|-- src/dsm123_pipeline/          reusable reconstruction package
|-- scripts/                      pipeline, audit, export, and plotting scripts
|-- genomes/                      reference and DSM123 GenBank/genome files
|-- prots/                        protein FASTA files and BLAST databases
|-- nucl/                         nucleotide FASTA files and BLAST results
|-- bbh/                          bidirectional best-hit results
|-- bigg_resources/               local BiGG models and sequence resources
|-- Models/                       draft, working, and final DSM123 models
|-- manual_curation_outputs_merged/
|   |-- step1_review/             GPR conflicts and missing-gene candidates
|   |-- step2_review/             metabolite, bridge, and orphan audits
|   `-- figure_data/              reconstruction statistics and figure inputs
|-- plot/                         downstream model applications and figures
|-- manual_curation.ipynb         phase-organized manual curation notebook
|-- metabolic_network_reconstruction_project.py
|-- pyproject.toml
`-- README.md
```

## Software requirements

- Python 3.10 or later
- BLAST+ (`makeblastdb`, `blastp`, and `blastn` available on `PATH`)
- A COBRApy-compatible optimization solver

Install the Python package and core dependencies from this directory:

```bash
python -m pip install -e .
```

The declared dependencies are `appdirs`, `biopython`, `cobra`, and `pandas`. The downstream notebooks in `plot/` additionally require Jupyter, NumPy, Matplotlib, seaborn, SciPy, and openpyxl.

## Automated reconstruction

Run from `workflow_1_model_reconstruction/`:

```bash
python scripts/run_pipeline.py \
  --target-id DSM123 \
  --target-gb genomes/DSM123.gb
```

The compatibility entry point is equivalent:

```bash
python metabolic_network_reconstruction_project.py \
  --target-id DSM123 \
  --target-gb genomes/DSM123.gb
```

Existing generated outputs are replaced by default. Add `--no-overwrite` only when previous outputs must be retained.

### Reconstruction thresholds

| Parameter | Value |
|---|---:|
| BLAST E-value cutoff | `1e-3` |
| BBH coverage threshold | `0.2` |
| Orthology PID threshold | `65%` |
| BLASTN PID threshold | `70%` |
| BLASTN alignment/query-length ratio | `0.8` |
| Unannotated ORF rescue PID | `80%` |

The thresholds are defined in `src/dsm123_pipeline/config.py` and recorded in `pipeline_run_report.json`.

## Reconstruction stages

### 1. Genome parsing and homology mapping

CDS features are extracted from GenBank files to protein and nucleotide FASTA files. Bidirectional BLASTP results are filtered by coverage and reduced to reciprocal best hits. The resulting PID matrix is converted into a binary orthology matrix using the 65% PID threshold.

Key outputs:

- `bbh/Rpal_BisA53_vs_DSM123_parsed.csv`
- `ortho_matrix.csv`
- `geneIDs_matrix.csv`
- `na_matrix.csv`

### 2. Draft model and flux-supported recovery

Unsupported reference genes and their reactions are removed, retained genes are renamed to DSM123 identifiers, and nucleotide evidence is used to rescue candidate unannotated ORFs. Reactions supported by reference-model flux or reduced-cost evidence are then restored.

Key outputs:

- `Models/DSM123.json`: orthology-derived draft
- `missing_reactions.csv`: reactions identified for recovery
- `reaction_gene_relationships.csv`: recovered reaction-GPR relationships
- `updated_consensus.json`: flux-recovered model
- `pipeline_run_report.json`: machine-readable run summary

The recorded DSM123 run identified and restored 210 flux-supported reactions. The reported biomass objectives were 0.709771 under the default medium and 0.216185 under the M9-like medium.

### 3. Manual curation

The phase-organized notebook is `manual_curation.ipynb`. Reproducible Step 1 and Step 2 review artifacts can be regenerated with:

```bash
python scripts/regenerate_step1_step2_reviews.py
```

Step 1 compares model GPR assignments with filtered BBH evidence (`identity >= 50%`, `E-value <= 1e-3`). The stored review contains 280 GPR conflicts, 154 missing-gene candidates, and 3,262 filtered BBH hits.

Step 2 performs `_u` metabolite normalization, quinolinate unification, formula-and-compartment duplicate review, reviewed metabolite merges, compartment-bridge auditing, and action-based orphan removal. The stored run applied 11 reviewed merge rows, added 82 cytosol-periplasm, 36 extracellular-periplasm, and 144 cytosol-extracellular bridge reactions, and removed 12 reviewed orphan metabolites plus 19 reactions made empty by those removals.

Core review tables are retained in `manual_curation_outputs_merged/step1_review/` and `manual_curation_outputs_merged/step2_review/`.

## Final model and recorded model sizes

The final model used for downstream analyses is:

```text
Models/purple_bacteriav_DSM123.json
```

| Stage | Genes | Metabolites | Reactions |
|---|---:|---:|---:|
| `iDT1294Photo` | 1,294 | 2,038 | 2,588 |
| `Models/DSM123.json` | 976 | 2,038 | 2,259 |
| `updated_consensus.json` | 1,168 | 2,038 | 2,469 |
| `Models/purple_bacteriav_DSM123.json` | 1,003 | 2,030 | 2,355 |

Relative to `updated_consensus.json`, the final curated model contains 32 added and 197 removed genes, 64 added and 178 removed reactions, and 114 reactions with changed GPR rules. These values describe the stored model transition; they are not a statistical enrichment test.

## Reconstruction summaries and figures

Generate stage statistics and final-model composition tables with:

```bash
python scripts/export_purple_stage_stats.py
python scripts/export_model_transition_appendix_tables.py
python scripts/plot_stage_and_gene_enrichment.py
```

The principal tabular outputs are under `manual_curation_outputs_merged/figure_data/`. Subsystem tables summarize model annotation composition. Entries labeled `Unknown` reflect missing subsystem annotation and must not be interpreted as pathway absence.

The KEGG export script may require internet access and local BLAST+:

```bash
python scripts/export_final_kegg_subsystem_enrichment.py
```

## Downstream applications in `plot/`

The `plot/` directory is the application layer of Workflow 1. It uses `../Models/purple_bacteriav_DSM123.json` and contains:

- `enzyme_constrained_dfba.ipynb`: RNA-seq-scaled enzyme constraints, active transport versus diffusion, ForT sensitivity, loopless FVA, and module-capacity scaling.
- `formate_light_scan.ipynb`: steady-state formate-light FBA scan.
- kinetic, substrate, RNA-seq allocation, and reaction-parameter tables.
- source data, SVG/PNG figures, and `Supplementary_Methods.docx`.

Run the notebooks from `plot/` so their relative paths resolve correctly:

```bash
cd plot
jupyter nbconvert --to notebook --execute --inplace enzyme_constrained_dfba.ipynb
jupyter nbconvert --to notebook --execute --inplace formate_light_scan.ipynb
```

See `plot/README.md` for equations, parameter choices, output tables, and interpretation limits.

## Reproducibility notes

1. Several stored JSON and CSV reports contain absolute paths from the original analysis workstation. These paths are provenance strings; executable scripts resolve inputs relative to the current project root.
2. BLAST database index files are retained to reproduce the archived run but can be regenerated from the corresponding FASTA files.
3. `updated_consensus.json` is the automated flux-recovered model, whereas `Models/purple_bacteriav_DSM123.json` is the manually curated final model used by `plot/`.
4. The downstream enzyme allocation uses RNA-seq as a relative proxy and a model-derived biomass protein fraction. It is not a direct absolute proteomics measurement.

