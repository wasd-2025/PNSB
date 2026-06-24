# DSM123 Metabolic Reconstruction and Tn-seq Analysis

This repository contains the metabolic model reconstruction and Tn-seq analyses for *Rhodopseudomonas palustris* DSM123.

## Repository structure

```text
.
|-- workflow_1_model_reconstruction/  metabolic reconstruction and model applications
`-- Tnseq/                             Tn-seq processing and visualization
```

## Workflow 1: metabolic reconstruction

`workflow_1_model_reconstruction/` contains:

- automated homology-based model reconstruction;
- flux-supported reaction recovery;
- GPR, metabolite, compartment, and orphan-metabolite curation;
- the final model `Models/purple_bacteriav_DSM123.json`;
- reconstruction statistics and figure tables;
- downstream dFBA, FVA, and kinetic analyses in `plot/`.

Install and run the reconstruction pipeline:

```bash
cd workflow_1_model_reconstruction
python -m pip install -e .
python scripts/run_pipeline.py --target-id DSM123 --target-gb genomes/DSM123.gb
```

The main manual-curation notebook is `manual_curation.ipynb`. See `workflow_1_model_reconstruction/README.md` for thresholds, reconstruction stages, recorded results, and additional commands.

## Downstream model analyses

The `workflow_1_model_reconstruction/plot/` directory contains enzyme-constrained dFBA, formate-transport sensitivity, loopless FVA, module-capacity scaling, and formate-light scan analyses.

Run the notebooks from the `plot/` directory:

```bash
cd workflow_1_model_reconstruction/plot
jupyter nbconvert --to notebook --execute --inplace enzyme_constrained_dfba.ipynb
jupyter nbconvert --to notebook --execute --inplace formate_light_scan.ipynb
```

## Tn-seq analysis

`Tnseq/` contains the R scripts, input tables, and figures used for pathway-level fitness heatmaps and genome-level visualization.

```r
setwd("Tnseq")
source("Tnseq_heatmap.R")
source("heatmap.R")
source("circos.R")
```

Generated Tn-seq figures and tables are written to `Tnseq/output/`. See `Tnseq/README.md` for R package requirements and output descriptions.

## Main outputs

- Final metabolic model: `workflow_1_model_reconstruction/Models/purple_bacteriav_DSM123.json`
- Reconstruction summaries: `workflow_1_model_reconstruction/manual_curation_outputs_merged/figure_data/`
- dFBA/FVA figures: `workflow_1_model_reconstruction/plot/`
- Tn-seq figures: `Tnseq/output/`

