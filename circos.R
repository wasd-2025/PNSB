# =============================================================================
# circos.R
# =============================================================================
# Purpose:
#   Draw a circular genome plot using files in ./input and write outputs to ./output.
#
# Required input:
#   - input/chr_length.txt
#
# Optional inputs (drawn when present):
#   - input/circos_log2FC_track.txt
#   - input/gc_content_result.tsv
#   - input/gene_density.tsv
#   - input/gene_density.txt
#   - input/repeat_density.txt
#   - input/LTR_density.txt
#   - input/GC_content.txt
#   - input/gene_pair_link.xls
#
# Main output:
#   - output/circos_genome.pdf
# =============================================================================

load_packages <- function(pkgs, section_name) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      sprintf(
        "Missing package(s) for %s: %s",
        section_name,
        paste(missing_pkgs, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

load_packages(c("gtools", "circlize", "RColorBrewer", "grid"), "circos plotting")

# Resolve base directory from current working directory.
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

input_dir <- file.path(base_dir, "input")
output_dir <- file.path(base_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_file <- file.path(input_dir, "chr_length.txt")
if (!file.exists(required_file)) {
  stop(sprintf("Required input not found: %s", required_file), call. = FALSE)
}

read_chr_info <- function(path) {
  x <- utils::read.delim(path, stringsAsFactors = FALSE, header = FALSE)
  if (ncol(x) < 2) {
    stop("chr_length.txt must contain at least 2 columns: contig and length.", call. = FALSE)
  }
  x <- x[, 1:2]
  colnames(x) <- c("chr", "end")
  x$start <- 1
  x <- x[, c("chr", "start", "end")]
  x
}

read_track_with_fallback <- function(primary_path, fallback_path = NULL, colnames_expected = NULL) {
  path <- NULL
  if (!is.null(primary_path) && file.exists(primary_path)) {
    path <- primary_path
  } else if (!is.null(fallback_path) && file.exists(fallback_path)) {
    path <- fallback_path
  }

  if (is.null(path)) {
    return(NULL)
  }

  dat <- tryCatch(
    utils::read.delim(path, stringsAsFactors = FALSE, header = TRUE),
    error = function(e) utils::read.delim(path, stringsAsFactors = FALSE, header = FALSE)
  )

  if (!is.null(colnames_expected) && ncol(dat) >= length(colnames_expected)) {
    colnames(dat)[seq_along(colnames_expected)] <- colnames_expected
  }

  dat
}

chr_info <- read_chr_info(required_file)
n_chr <- nrow(chr_info)

palette_fun <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))

output_pdf <- file.path(output_dir, "circos_genome.pdf")
grDevices::pdf(output_pdf, width = 14, height = 10)
on.exit(grDevices::dev.off(), add = TRUE)

circlize::circos.clear()
circlize::circos.par(start.degree = 90, gap.degree = c(rep(2, n_chr - 1), 10))
circlize::circos.genomicInitialize(chr_info, plotType = "axis")

# Inner chromosome label ring
circlize::circos.track(
  ylim = c(0, 1),
  track.height = 0.08,
  bg.border = NA,
  bg.col = palette_fun(n_chr),
  panel.fun = function(x, y) {
    xlim <- circlize::CELL_META$xlim
    chr_id <- circlize::CELL_META$sector.index
    circlize::circos.text(
      x = mean(xlim),
      y = 0.5,
      labels = chr_id,
      cex = 0.7,
      col = "black",
      facing = "inside",
      niceFacing = FALSE
    )
  }
)

# Optional track: RNA-seq style log2FC blocks
rna_track <- read_track_with_fallback(
  primary_path = file.path(input_dir, "circos_log2FC_track.txt"),
  colnames_expected = c("chr", "start", "end", "value")
)
if (!is.null(rna_track) && ncol(rna_track) >= 4) {
  rna_track <- rna_track[, 1:4]
  colnames(rna_track) <- c("chr", "start", "end", "value")
  col_fun_rna <- circlize::colorRamp2(c(-5, 0, 6), c("#4575b4", "#ffffbf", "#d73027"))

  circlize::circos.genomicTrackPlotRegion(
    rna_track,
    stack = TRUE,
    track.height = 0.10,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      circlize::circos.genomicRect(region, value, col = col_fun_rna(value[[1]]), border = NA, ...)
    }
  )
}

# Optional track: gene density (tsv format)
gene_density_tsv <- read_track_with_fallback(
  primary_path = file.path(input_dir, "gene_density.tsv")
)
if (!is.null(gene_density_tsv) && all(c("Contig", "Start", "End", "Gene_Count") %in% colnames(gene_density_tsv))) {
  gene_track <- gene_density_tsv[, c("Contig", "Start", "End", "Gene_Count")]
  colnames(gene_track) <- c("chr", "start", "end", "value")
  col_fun_gene <- circlize::colorRamp2(
    c(min(gene_track$value, na.rm = TRUE), mean(gene_track$value, na.rm = TRUE), max(gene_track$value, na.rm = TRUE)),
    c("#FFE5E0", "#E35036", "#71281B")
  )

  circlize::circos.genomicTrackPlotRegion(
    gene_track,
    stack = TRUE,
    track.height = 0.10,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      circlize::circos.genomicRect(region, value, col = col_fun_gene(value[[1]]), border = NA, ...)
    }
  )
}

# Optional track: gene density (legacy txt format)
gene_density_txt <- read_track_with_fallback(
  primary_path = file.path(input_dir, "gene_density.txt")
)
if (!is.null(gene_density_txt) && all(c("number") %in% colnames(gene_density_txt)) && ncol(gene_density_txt) >= 4) {
  gd <- gene_density_txt[, 1:4]
  colnames(gd) <- c("chr", "start", "end", "number")
  col_fun <- circlize::colorRamp2(
    c(min(gd$number, na.rm = TRUE), mean(gd$number, na.rm = TRUE), max(gd$number, na.rm = TRUE)),
    c("white", "#E35036", "#71281B")
  )
  circlize::circos.genomicTrackPlotRegion(
    gd,
    stack = TRUE,
    track.height = 0.08,
    panel.fun = function(region, value, ...) {
      circlize::circos.genomicRect(region, value, col = col_fun(value[[1]]), border = NA, ...)
    },
    bg.border = NA
  )
}

# Optional track: GC content (tsv format)
gc_tsv <- read_track_with_fallback(
  primary_path = file.path(input_dir, "gc_content_result.tsv")
)
if (!is.null(gc_tsv)) {
  gc_col <- grep("GC_Content|GC", colnames(gc_tsv), value = TRUE)[1]
  if (!is.na(gc_col) && all(c("Contig", "Start", "End") %in% colnames(gc_tsv))) {
    gc_track <- gc_tsv[, c("Contig", "Start", "End", gc_col)]
    colnames(gc_track) <- c("chr", "start", "end", "value")
    mean_gc <- mean(gc_track$value, na.rm = TRUE)

    circlize::circos.genomicTrack(
      gc_track,
      track.height = 0.12,
      bg.col = "#EEEEEE6E",
      bg.border = NA,
      panel.fun = function(region, value, ...) {
        circlize::circos.genomicLines(region, value, col = "#7126D1", lwd = 0.4, ...)
        xlim <- circlize::CELL_META$xlim
        circlize::circos.lines(xlim, c(mean_gc, mean_gc), col = "blue2", lwd = 0.2, lty = 2)
      }
    )
  }
}

# Optional track: GC content (legacy txt format)
gc_txt <- read_track_with_fallback(
  primary_path = file.path(input_dir, "GC_content.txt")
)
if (!is.null(gc_txt) && all(c("number") %in% colnames(gc_txt)) && ncol(gc_txt) >= 4) {
  gc <- gc_txt[, 1:4]
  colnames(gc) <- c("chr", "start", "end", "number")
  mean_gc <- mean(gc$number, na.rm = TRUE)

  circlize::circos.genomicTrack(
    gc,
    track.height = 0.08,
    bg.col = "#EEEEEE6E",
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      circlize::circos.genomicLines(region, value, col = "#7126D1", lwd = 0.35, ...)
      xlim <- circlize::CELL_META$xlim
      circlize::circos.lines(xlim, c(mean_gc, mean_gc), col = "blue2", lwd = 0.15, lty = 2)
    }
  )
}

# Optional legacy tracks: repeat/LTR density
repeat_density <- read_track_with_fallback(file.path(input_dir, "repeat_density.txt"))
if (!is.null(repeat_density) && all(c("number") %in% colnames(repeat_density)) && ncol(repeat_density) >= 4) {
  rp <- repeat_density[, 1:4]
  colnames(rp) <- c("chr", "start", "end", "number")
  col_fun <- circlize::colorRamp2(
    c(min(rp$number, na.rm = TRUE), mean(rp$number, na.rm = TRUE), max(rp$number, na.rm = TRUE)),
    c("white", "#6FC5FF", "#0073BF")
  )
  circlize::circos.genomicTrackPlotRegion(
    rp,
    stack = TRUE,
    track.height = 0.08,
    panel.fun = function(region, value, ...) {
      circlize::circos.genomicRect(region, value, col = col_fun(value[[1]]), border = NA, ...)
    },
    bg.border = NA
  )
}

ltr_density <- read_track_with_fallback(file.path(input_dir, "LTR_density.txt"))
if (!is.null(ltr_density) && all(c("number") %in% colnames(ltr_density)) && ncol(ltr_density) >= 4) {
  lt <- ltr_density[, 1:4]
  colnames(lt) <- c("chr", "start", "end", "number")
  col_fun <- circlize::colorRamp2(
    c(min(lt$number, na.rm = TRUE), mean(lt$number, na.rm = TRUE), max(lt$number, na.rm = TRUE)),
    c("white", "#59EC6B", "#19421E")
  )
  circlize::circos.genomicTrackPlotRegion(
    lt,
    stack = TRUE,
    track.height = 0.08,
    panel.fun = function(region, value, ...) {
      circlize::circos.genomicRect(region, value, col = col_fun(value[[1]]), border = NA, ...)
    },
    bg.border = NA
  )
}

# Optional links
link_file <- file.path(input_dir, "gene_pair_link.xls")
if (file.exists(link_file)) {
  link_tab <- utils::read.table(link_file, sep = "\t", header = FALSE, check.names = FALSE)
  if (ncol(link_tab) >= 6) {
    chr_split <- split(link_tab[, 1:6], link_tab$V1)
    chr_split <- chr_split[gtools::mixedorder(names(chr_split))]

    for (n in seq_along(chr_split)) {
      tmp_link <- do.call(rbind, chr_split[n])
      chr_name <- unique(chr_split[[n]]$V1)
      row_idx <- which(chr_info$chr == chr_name)
      if (length(row_idx) > 0) {
        circlize::circos.genomicLink(
          tmp_link[c(1, 2, 3)],
          tmp_link[c(4, 5, 6)],
          col = palette_fun(n_chr)[row_idx[1]],
          lwd = 0.05
        )
      }
    }
  }
}

circlize::circos.clear()
message(sprintf("Circos plot written: %s", output_pdf))
dev.off()
