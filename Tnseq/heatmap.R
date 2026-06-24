library(dplyr)
library(forcats)
library(tidyr)
library(colorspace)
library(lattice)

has_geneplotter <- requireNamespace("geneplotter", quietly = TRUE)
has_copynumber <- requireNamespace("copynumber", quietly = TRUE)
if (has_geneplotter) {
  suppressPackageStartupMessages(library(geneplotter))
}
if (has_copynumber) {
  suppressPackageStartupMessages(library(copynumber))
}

# Keep compatibility with older scripts that expected panel.geneplot to exist
# after loading geneplotter/copynumber.
if (!exists("panel.geneplot", mode = "function")) {
  if (has_geneplotter &&
      exists("panel.geneplot", envir = asNamespace("geneplotter"), mode = "function", inherits = FALSE)) {
    panel.geneplot <- get("panel.geneplot", envir = asNamespace("geneplotter"))
  } else if (has_copynumber &&
             exists("panel.geneplot", envir = asNamespace("copynumber"), mode = "function", inherits = FALSE)) {
    panel.geneplot <- get("panel.geneplot", envir = asNamespace("copynumber"))
  } else {
    panel.geneplot <- function(x, y, arrows = TRUE, tip = 0.1, rot_labels = 0,
                               gene_strand = NULL, gene_name = NULL,
                               subscripts = seq_along(x), ...) {
      if (length(x) == 0) return(invisible(NULL))

      strand <- rep("+", length(x))
      if (!is.null(gene_strand)) {
        strand <- if (length(subscripts) > 0 && length(gene_strand) >= max(subscripts)) {
          gene_strand[subscripts]
        } else {
          rep_len(gene_strand, length(x))
        }
      }

      gene_labels <- rep("", length(x))
      if (!is.null(gene_name)) {
        gene_labels <- if (length(subscripts) > 0 && length(gene_name) >= max(subscripts)) {
          gene_name[subscripts]
        } else {
          rep_len(gene_name, length(x))
        }
      }

      x_left <- pmin(x, y)
      x_right <- pmax(x, y)
      y_mid <- 0
      half_height <- 0.24

      if (isTRUE(arrows)) {
        strand_code <- toupper(trimws(as.character(strand)))
        points_right <- strand_code %in% c("+", "1", "F", "FORWARD", "PLUS")
        points_left <- strand_code %in% c("-", "-1", "R", "REVERSE", "MINUS")
        points_right[!(points_right | points_left)] <- x[!(points_right | points_left)] <= y[!(points_right | points_left)]

        span <- diff(range(c(x_left, x_right), na.rm = TRUE))
        if (!is.finite(span) || span <= 0) span <- 1

        gene_width <- x_right - x_left
        head_width <- pmin(gene_width * 0.35, span * 0.04)
        fill_col <- ifelse(points_right, "#E69F00", "#0072B2")

        for (i in seq_along(x_left)) {
          if (!is.finite(gene_width[i]) || gene_width[i] <= 0) next

          if (points_right[i]) {
            body_end <- max(x_left[i], x_right[i] - head_width[i])
            px <- c(x_left[i], body_end, x_right[i], body_end, x_left[i])
          } else {
            body_start <- min(x_right[i], x_left[i] + head_width[i])
            px <- c(x_right[i], body_start, x_left[i], body_start, x_right[i])
          }

          py <- c(
            y_mid - half_height,
            y_mid - half_height,
            y_mid,
            y_mid + half_height,
            y_mid + half_height
          )

          lattice::panel.polygon(px, py, col = fill_col[i], border = "black", lwd = 1)

          tail_x <- if (points_right[i]) x_left[i] else x_right[i]
          lattice::panel.segments(
            tail_x,
            y_mid - half_height,
            tail_x,
            y_mid + half_height,
            col = "black",
            lwd = 1
          )
        }
      } else {
        lattice::panel.rect(
          xleft = x_left,
          ybottom = y_mid - half_height,
          xright = x_right,
          ytop = y_mid + half_height,
          col = "#E69F00",
          border = "black"
        )
      }

      label_x <- (x_left + x_right) / 2
      label_y <- rep(-0.55, length(x))
      span <- diff(range(c(x_left, x_right), na.rm = TRUE))
      if (!is.finite(span) || span <= 0) span <- 1
      close_cutoff <- span * 0.06
      label_order <- order(label_x)
      for (j in seq_along(label_order)[-1]) {
        current <- label_order[j]
        previous <- label_order[j - 1]
        if (is.finite(label_x[current] - label_x[previous]) &&
            abs(label_x[current] - label_x[previous]) < close_cutoff) {
          label_y[current] <- -0.85
        }
      }
      lattice::panel.text(label_x, label_y, labels = gene_labels, cex = 0.72)
    }
  }
}

# Load prerequisite objects when running heatmap.R standalone.
if (!exists("fitness123_annotated", inherits = TRUE) ||
    !exists("group1_locus", inherits = TRUE) ||
    !exists("group2_locus", inherits = TRUE) ||
    !exists("group3_locus", inherits = TRUE) ||
    !exists("group4_locus", inherits = TRUE)) {
  source("Tnseq_heatmap.R")
}

if (!exists("output_dir", inherits = TRUE)) {
  output_dir <- "output"
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("outfile", mode = "function", inherits = TRUE)) {
  outfile <- function(name) {
    file.path(output_dir, name)
  }
}

fitness123_annotated1 <- fitness123_annotated %>%
  mutate(time = recode(time, `0` = "T0-YE", `1` = "T1-Formate", `2` = "T2-Formate"))

# Use locusId as the heatmap label, so no manual gene-name input is needed.
fill_gene_name <- function(df) {
  df$gene_name <- as.character(df$locusId)
  df
}

heatmap_fitness <- function(data, key = TRUE, max_value = 6) {
  heat_cols <- diverging_hcl(n = 7, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3))
  levelplot(
    norm_gene_fitness_mean ~ factor(gene_name) * fct_rev(factor(time)),
    data = data,
    col.regions = colorRampPalette(heat_cols)(25),
    at = seq(-6, 6, 0.5),
    aspect = "iso",
    xlab = "",
    ylab = "",
    scales = list(cex = 0.7, x = list(rot = 90)),
    panel = function(x, y, z, ...) {
      panel.levelplot(x, y, z, ...)
      panel.abline(v = seq_along(unique(x)) + 0.5,
                   h = seq_along(unique(y)) + 0.5,
                   col = "white", lwd = 2)
    }
  )
}

genome_plot <- function(df, xlim = NULL, title = "", rot_labels = 0) {
  df <- replace_na(df, list(strains_per_gene = 0))
  plot_theme <- list(
    axis.line = list(col = grey(0.3)),
    axis.text = list(col = grey(0.3), cex = 0.7),
    background = list(col = grey(0.7))
  )

  if (!is.null(xlim)) {
    xscale <- list(limits = xlim)
  } else {
    xscale <- list()
  }

  xyplot(
    end / 1000 ~ start / 1000, df,
    groups = strand, cex = 0.7, lwd = 1,
    par.settings = plot_theme, strains = df$strains_per_gene,
    scales = list(y = list(draw = FALSE), x = xscale),
    ylim = c(-3, 2), xlab = "", ylab = "",
    gene_strand = df[["strand"]],
    gene_name = df[["gene_name"]],
    panel = function(x, y, strains = NULL, gene_strand = NULL, gene_name = NULL,
                     subscripts = seq_along(x), ...) {
      if (length(x) == 0) return(invisible(NULL))
      if (is.null(strains)) strains <- rep("", length(x))
      strain_labels <- if (length(subscripts) > 0 && length(strains) >= max(subscripts)) strains[subscripts] else strains
      gene_strand_panel <- if (length(subscripts) > 0 && length(gene_strand) >= max(subscripts)) gene_strand[subscripts] else gene_strand
      gene_name_panel <- if (length(subscripts) > 0 && length(gene_name) >= max(subscripts)) gene_name[subscripts] else gene_name
      panel_limits <- lattice::current.panel.limits()
      panel.segments(panel_limits$xlim[1], 0, panel_limits$xlim[2], 0, col = "black", lwd = 1)
      panel.geneplot(
        x,
        y,
        subscripts = seq_along(x),
        gene_strand = gene_strand_panel,
        gene_name = gene_name_panel,
        arrows = TRUE,
        tip = 0.1,
        rot_labels = rot_labels,
        ...
      )
      panel.text((x + y) / 2, rep(0, length(x)), labels = strain_labels, cex = 0.7)
    }
  )
}

prepare_genome_df <- function(df) {
  df %>%
    mutate(gene_name = ifelse(is.na(gene) | gene == "", as.character(locusId), as.character(gene))) %>%
    group_by(locusId, gene, scaffold, strand, start, end, gene_name) %>%
    summarize(strains_per_gene = min(unique(strains_per_gene)), .groups = "drop") %>%
    arrange(start, end)
}

save_group_combined <- function(filename, heatmap_plot, genome_plots, genome_positions,
                                heatmap_position = c(0, 0.30, 1, 1.00),
                                width = 10, height = 8) {
  svg(outfile(filename), width = width, height = height)
  print(heatmap_plot, position = heatmap_position, more = TRUE)
  for (i in seq_along(genome_plots)) {
    more_flag <- i < length(genome_plots)
    print(genome_plots[[i]], position = genome_positions[[i]], more = more_flag)
  }
  dev.off()
}

# FDH ------------------------------------------------------------------------
fdh_df <- fitness123_annotated1 %>%
  filter(locusId %in% group1_locus) %>%
  fill_gene_name()

plot_fdh_fit <- heatmap_fitness(fdh_df)
fdh_genome_df <- prepare_genome_df(fdh_df)
plot_moco_g1 <- genome_plot(fdh_genome_df, xlim = c(4948.7, 4957.2), title = "chr 1")

save_group_combined(
  filename = "fdh_heatmap_genome.svg",
  heatmap_plot = plot_fdh_fit,
  genome_plots = list(plot_moco_g1),
  genome_positions = list(c(0.08, 0.00, 0.92, 0.36)),
  heatmap_position = c(0, 0.38, 1, 1.00)
)

# CBB ------------------------------------------------------------------------
preset_genes <- c("rbcS", "rlp2")
cbb_df <- fitness123_annotated1 %>%
  filter(locusId %in% group2_locus) %>%
  fill_gene_name() %>%
  mutate(order = case_when(
    gene %in% preset_genes ~ 1,
    TRUE ~ 2
  )) %>%
  arrange(order) %>%
  select(-order)

plot_fdh_fit <- heatmap_fitness(cbb_df)
cbb_genome_df <- prepare_genome_df(cbb_df)
plot_moco_g1 <- genome_plot(cbb_genome_df, xlim = c(1247, 1250.1), title = "chr 1")
plot_moco_g2 <- genome_plot(cbb_genome_df, xlim = c(4462.1, 4469.3), title = "chr 1")
plot_moco_g3 <- genome_plot(cbb_genome_df, xlim = c(3757.7, 3762.6), title = "chr 1")

save_group_combined(
  filename = "cbb_heatmap_genome.svg",
  heatmap_plot = plot_fdh_fit,
  genome_plots = list(plot_moco_g1, plot_moco_g2, plot_moco_g3),
  genome_positions = list(
    c(0.00, 0.00, 0.33, 0.34),
    c(0.335, 0.00, 0.665, 0.34),
    c(0.67, 0.00, 1.00, 0.34)
  )
)

# ETC ------------------------------------------------------------------------
etc_df <- fitness123_annotated1 %>%
  filter(locusId %in% group3_locus) %>%
  fill_gene_name()

plot_fdh_fit <- heatmap_fitness(etc_df)
etc_genome_df <- prepare_genome_df(etc_df)
plot_moco_g1 <- genome_plot(etc_genome_df, xlim = c(4982.4, 4987.2), title = "chr 1")
plot_moco_g2 <- genome_plot(etc_genome_df, xlim = c(4832.4, 4838), title = "chr 1")

save_group_combined(
  filename = "etc_heatmap_genome.svg",
  heatmap_plot = plot_fdh_fit,
  genome_plots = list(plot_moco_g1, plot_moco_g2),
  genome_positions = list(
    c(0.05, 0.00, 0.48, 0.34),
    c(0.52, 0.00, 0.95, 0.34)
  )
)

# RC -------------------------------------------------------------------------
rc_df <- fitness123_annotated1 %>%
  filter(locusId %in% group4_locus) %>%
  fill_gene_name()

plot_fdh_fit <- heatmap_fitness(rc_df)
rc_genome_df <- prepare_genome_df(rc_df)
plot_moco_g1 <- genome_plot(rc_genome_df, xlim = c(4006.1, 4009), title = "chr 1")
plot_moco_g2 <- genome_plot(rc_genome_df, xlim = c(1230.1, 1237.8), title = "chr 1")

save_group_combined(
  filename = "rc_heatmap_genome.svg",
  heatmap_plot = plot_fdh_fit,
  genome_plots = list(plot_moco_g1, plot_moco_g2),
  genome_positions = list(
    c(0.05, 0.00, 0.48, 0.34),
    c(0.52, 0.00, 0.95, 0.34)
  )
)

# 130.csv heatmap ------------------------------------------------------------
file_130 <- if (file.exists("130.csv")) "130.csv" else file.path("input", "130.csv")
if (file.exists(file_130)) {
  df130 <- read.csv(file_130, stringsAsFactors = FALSE)

  fitness123_annotated2 <- fitness123_annotated %>%
    mutate(time = recode(time, `0` = "0-6", `1` = "0-48", `2` = "6-48"))

  heatmap_fitness_130 <- function(data, key = TRUE, max_value = 8) {
    heat_cols <- diverging_hcl(n = 8, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3))
    levelplot(
      Log2FC_mean ~ factor(locusId) * fct_rev(factor(time)),
      data = data,
      col.regions = colorRampPalette(heat_cols)(25),
      at = seq(-5, 8, 0.6),
      aspect = "iso",
      xlab = "",
      ylab = "",
      scales = list(cex = 0.7, x = list(rot = 90)),
      panel = function(x, y, z, ...) {
        panel.levelplot(x, y, z, ...)
        panel.abline(v = seq_along(unique(x)) + 0.5,
                     h = seq_along(unique(y)) + 0.5,
                     col = "white", lwd = 2)
      }
    )
  }

  plot_fdh_fit <- heatmap_fitness_130(df130)
  svg(outfile("log2fc_130_heatmap.svg"), height = 6)
  print(plot_fdh_fit, position = c(0, 0.25, 1, 1.05), more = FALSE)
  dev.off()
}

