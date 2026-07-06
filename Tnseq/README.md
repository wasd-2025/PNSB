# PNSB Tn-seq Analysis

Minimal R scripts and input files for DSM123 Tn-seq visualization.

## Structure

- `Tnseq_heatmap.R`: builds annotated fitness tables and overview figures.
- `heatmap.R`: draws pathway heatmaps together with genome cluster plots.
- `circos.R`: draws the genome-level circos plot.
- `input/`: required input data.
- `output/`: generated figures and tables.

## Requirements

Install the CRAN packages used by the scripts:

```r
install.packages(c(
  "dplyr", "tidyr", "forcats", "colorspace", "lattice",
  "ggplot2", "ggrepel", "cluster", "RColorBrewer", "tibble",
  "readr", "circlize"
))
```

Install the Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ComplexHeatmap")
```

Optional packages:

```r
install.packages("Rtsne")
```

## Run

Run from the project root:

```r
setwd("path/to/Tnseq")
source("Tnseq_heatmap.R")
source("heatmap.R")
source("circos.R")
```

`heatmap.R` will automatically source `Tnseq_heatmap.R` if the required
objects are not already in the R environment.

## Main Outputs

All outputs are written to `output/`.

- `fdh_heatmap_genome.svg`
- `cbb_heatmap_genome.svg`
- `etc_heatmap_genome.svg`
- `rc_heatmap_genome.svg`
- `log2fc_130_heatmap.svg`
- `combined_plots.svg`
- `grouped_heatmap.svg`
- `circos_genome.pdf`
