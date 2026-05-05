# Tn-seq Git Upload (Minimal)

This folder is the merged minimal package from the two previous locations:
- main analysis folder
- annotation/circos inputs

## Structure

- `Tnseq_heatmap.R`
- `heatmap.R`
- `circos.R`
- `input/` (all required input files)
- `output/` (generated files)
- `.gitignore`

## Input Files in `input/`

- `130.csv`
- `cga009_vs_consensus_parsed.csv`
- `consensus.gff`
- `fitness.Rdata`
- `fitness_gene.Rdata`
- `reaction_data(dsm123).csv`
- `chr_length.txt`
- `circos_log2FC_track.txt`
- `gc_content_result.tsv`
- `gene_density.tsv`

## Run

```r
source("Tnseq_heatmap.R")
source("heatmap.R")
source("circos.R")
```

`heatmap.R` now auto-sources `Tnseq_heatmap.R` when required objects are missing.

All generated figures/tables are written to `output/`.
