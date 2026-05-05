# =============================================================================
# Tnseq_heatmap.R
# =============================================================================
# Purpose:
#   Build annotated fitness tables and generate overview plots for DSM123 Tn-seq.
#
# Main outputs:
#   - combined_plots.svg
#   - figure_silhouette.svg
#   - heatmap.svg
#   - tSNE_plot_with_labels.pdf
#   - kegg_enrichment_plot.svg
#   - kegg_bubble_plot.svg
#   - kegg_enrichment_bubble.png
#   - fitness123_with_category.csv
#   - gene_cluster_dendrogram.png
#   - heatmap(22).svg
#
# Notes:
#   1) This script also creates `fitness123_annotated` and locus groups used by
#      heatmap.R.
#   2) Each figure section loads only the packages it needs.
# =============================================================================

# ---- Utility: package loader ------------------------------------------------
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
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }
}

# ---- Utility: basic guards --------------------------------------------------
require_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required file not found: %s", path), call. = FALSE)
  }
}

save_placeholder_plot <- function(filename, message, width = 8, height = 6) {
  ext <- tolower(tools::file_ext(filename))

  if (ext == "svg") {
    grDevices::svg(filename, width = width, height = height)
  } else if (ext == "pdf") {
    grDevices::pdf(filename, width = width, height = height)
  } else if (ext == "png") {
    grDevices::png(filename, width = width, height = height, units = "in", res = 150)
  } else {
    stop(sprintf("Unsupported placeholder format: %s", filename), call. = FALSE)
  }

  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot.new()
  graphics::text(0.5, 0.5, message, cex = 1)
}

save_ggplot_svg <- function(plot_obj, filename, width = 10, height = 8) {
  grDevices::svg(filename, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot_obj)
}

save_ggplot_png <- function(plot_obj, filename, width = 10, height = 6, dpi = 300) {
  grDevices::png(filename, width = width, height = height, units = "in", res = dpi)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot_obj)
}

# ---- I/O directories --------------------------------------------------------
input_dir <- "input"
output_dir <- "output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

infile <- function(name) {
  file.path(input_dir, name)
}

outfile <- function(name) {
  file.path(output_dir, name)
}

# ---- 1) Build annotated dataset --------------------------------------------
load_packages(
  c("dplyr", "tidyr", "tibble", "stringr", "rtracklayer"),
  "data preparation"
)

require_file(infile("fitness_gene.Rdata"))
require_file(infile("consensus.gff"))

load(infile("fitness_gene.Rdata"))
if (file.exists(infile("fitness.Rdata"))) {
  load(infile("fitness.Rdata"))
}

if (!exists("fitness_gene")) {
  stop("Object `fitness_gene` was not found after loading fitness_gene.Rdata.", call. = FALSE)
}

extract_gff_attr <- function(attr_vector, key) {
  if (is.null(attr_vector)) {
    return(rep(NA_character_, 0))
  }
  pattern <- paste0(key, "=([^;]+)")
  stringr::str_match(attr_vector, pattern)[, 2]
}

fitness123 <- fitness_gene %>%
  dplyr::group_by(locusId, scaffold, time, strains_per_gene) %>%
  dplyr::summarise(
    norm_gene_fitness_mean = mean(norm_gene_fitness, na.rm = TRUE),
    log2FC_mean = mean(log2FC, na.rm = TRUE),
    tstat_mean = mean(t, na.rm = TRUE),
    .groups = "drop"
  )

gff_df <- as.data.frame(rtracklayer::import(infile("consensus.gff")))
attr_col <- if ("attributes" %in% names(gff_df)) as.character(gff_df$attributes) else rep(NA_character_, nrow(gff_df))

gff_df$locus_col <- if ("ID" %in% names(gff_df)) as.character(gff_df$ID) else extract_gff_attr(attr_col, "ID")
gff_df$name_col <- if ("Name" %in% names(gff_df)) as.character(gff_df$Name) else extract_gff_attr(attr_col, "Name")
gff_df$gene_col <- if ("gene" %in% names(gff_df)) as.character(gff_df$gene) else extract_gff_attr(attr_col, "gene")

genes123 <- gff_df %>%
  dplyr::filter(type == "CDS") %>%
  dplyr::transmute(
    locusId = locus_col,
    start = as.numeric(start),
    end = as.numeric(end),
    strand = as.character(strand),
    Name = name_col,
    gene = gene_col
  ) %>%
  dplyr::filter(!is.na(locusId), locusId != "") %>%
  dplyr::distinct(locusId, .keep_all = TRUE)

fitness123_annotated <- fitness123 %>%
  dplyr::left_join(genes123, by = "locusId") %>%
  dplyr::mutate(time = as.numeric(as.character(time)))

fitness123_annotated_com <- fitness123_annotated %>%
  dplyr::select(locusId, time, gene, Name, norm_gene_fitness_mean) %>%
  tidyr::pivot_wider(names_from = time, values_from = norm_gene_fitness_mean)

needed_time_cols <- c("0", "1", "2")
missing_time_cols <- setdiff(needed_time_cols, colnames(fitness123_annotated_com))
if (length(missing_time_cols) > 0) {
  stop(
    sprintf(
      "The wide fitness table is missing expected time column(s): %s",
      paste(missing_time_cols, collapse = ", ")
    ),
    call. = FALSE
  )
}

# ---- 2) Figure: combined time-course line plots -----------------------------
load_packages(c("lattice", "gridExtra", "RColorBrewer"), "combined line plots")

stdcol <- RColorBrewer::brewer.pal(8, "Dark2")

plot_hist_log2FC <- lattice::xyplot(
  log2FC_mean ~ time,
  data = fitness123_annotated,
  groups = locusId,
  as.table = TRUE,
  col = stdcol[1],
  alpha = 0.8,
  layout = c(1, 1),
  xlab = "Timepoints",
  ylab = expression("Log"[2] * " FC"),
  xlim = c(0, 2),
  type = "l",
  between = list(x = 0.5, y = 0.5),
  scales = list(
    alternating = FALSE,
    x = list(at = c(0, 1, 2), labels = c("T0-YE", "T1-Formate", "T2-Formate"), rot = 45, cex = 0.8)
  ),
  lwd = 2,
  panel = function(x, y, ...) {
    lattice::panel.grid(h = -1, v = -1, col = grey(0.9))
    lattice::panel.xyplot(x, y, ...)
  }
)

plot_hist_normfg <- lattice::xyplot(
  norm_gene_fitness_mean ~ time,
  data = fitness123_annotated,
  groups = locusId,
  as.table = TRUE,
  col = stdcol[2],
  alpha = 0.2,
  layout = c(1, 1),
  xlab = "Timepoints",
  ylab = "Fitness",
  type = "l",
  xlim = c(0, 2),
  between = list(x = 0.5, y = 0.5),
  scales = list(
    alternating = FALSE,
    x = list(at = c(0, 1, 2), labels = c("T0-YE", "T1-Formate", "T2-Formate"), rot = 45, cex = 0.8)
  ),
  lwd = 2,
  par.settings = list(superpose.line = list(col = stdcol)),
  panel = function(x, y, ...) {
    lattice::panel.grid(h = -1, v = -1, col = grey(0.9))
    lattice::panel.xyplot(x, y, ...)
  }
)

grDevices::svg(outfile("combined_plots.svg"), width = 12, height = 6)
gridExtra::grid.arrange(plot_hist_log2FC, plot_hist_normfg, ncol = 2)
grDevices::dev.off()

# ---- 3) Matrix preparation for clustering and heatmap -----------------------
load_packages(
  c("cluster", "colorspace", "ComplexHeatmap", "circlize", "RColorBrewer", "dplyr", "tibble"),
  "clustering and heatmap"
)

mat_heatmap <- fitness123_annotated_com %>%
  dplyr::filter(
    dplyr::if_any(
      dplyr::all_of(needed_time_cols),
      ~ !dplyr::between(., -1.5, 0.5)
    )
  ) %>%
  dplyr::select(locusId, dplyr::all_of(needed_time_cols)) %>%
  tibble::column_to_rownames("locusId") %>%
  as.matrix()

mat_heatmap <- mat_heatmap[stats::complete.cases(mat_heatmap), , drop = FALSE]

if (nrow(mat_heatmap) < 2) {
  stop("Not enough rows in `mat_heatmap` to run clustering.", call. = FALSE)
}

mat_cluster <- stats::hclust(stats::dist(mat_heatmap), method = "ward.D2")
mat_heatmap_ordered <- mat_heatmap[
  order.dendrogram(as.dendrogram(mat_cluster)),
  ,
  drop = FALSE
]

n_clusters <- min(6, max(2, nrow(mat_heatmap) - 1))

# Cluster assignment for the heatmap / t-SNE branch
cluster_assignment_main <- stats::cutree(mat_cluster, k = n_clusters)

cluster_df_main <- tibble::tibble(
  locusId = rownames(mat_heatmap),
  cluster = as.factor(cluster_assignment_main[rownames(mat_heatmap)])
)

# ---- 4) Figure: silhouette plot ---------------------------------------------
if (nrow(mat_heatmap) >= 3) {
  set.seed(123)
  kmeans_result <- stats::kmeans(mat_heatmap, centers = n_clusters)
  silhouette_result <- cluster::silhouette(
    kmeans_result$cluster,
    stats::dist(mat_heatmap)
  )
  
  grDevices::svg(outfile("figure_silhouette.svg"), width = 4, height = 4)
  graphics::plot(
    silhouette_result,
    col = RColorBrewer::brewer.pal(8, "Dark2")
  )
  grDevices::dev.off()
} else {
  save_placeholder_plot(
    filename = outfile("figure_silhouette.svg"),
    message = "Not enough points for silhouette analysis.",
    width = 4,
    height = 4
  )
}


# ---- 5) Figure: overview heatmap --------------------------------------------
transposed_heatmap_matrix <- t(mat_heatmap_ordered)

main_heatmap <- ComplexHeatmap::Heatmap(
  matrix = transposed_heatmap_matrix,
  name = "Fitness",
  col = circlize::colorRamp2(
    c(-6, -3, 0, 1, 2),
    c("#002F70", "#829EDA", "#E2E2E2", "#E6BCC3", "#E495A5")
  ),
  cluster_rows = FALSE,
  cluster_columns = as.dendrogram(mat_cluster),
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 45,
  row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 6),
  heatmap_height = grid::unit(4, "cm"),
  width = grid::unit(20, "cm"),
  show_parent_dend_line = FALSE
)

grDevices::svg(outfile("heatmap.svg"), width = 20, height = 5)
ComplexHeatmap::draw(main_heatmap)
grDevices::dev.off()


# ---- 6) Figure: t-SNE with gene labels (old workflow reproduction) ----------
load_packages(
  c("lattice", "RColorBrewer", "tibble", "dplyr", "magrittr"),
  "t-SNE figure"
)

if (!requireNamespace("tsne", quietly = TRUE)) {
  save_placeholder_plot(
    filename = outfile("tSNE_plot_with_labels.pdf"),
    message = "Package tsne is not installed.",
    width = 10,
    height = 8
  )
} else if (nrow(mat_heatmap) < 5) {
  save_placeholder_plot(
    filename = outfile("tSNE_plot_with_labels.pdf"),
    message = "Not enough points for t-SNE.",
    width = 10,
    height = 8
  )
} else {
  set.seed(123)
  
  # Old workflow: t-SNE is performed on the distance matrix
  tsne_result <- tsne::tsne(
    stats::dist(mat_heatmap),
    max_iter = 500
  )
  
  df_tsne <- tsne_result %>%
    stats::setNames(c("V1", "V2")) %>%
    as.data.frame() %>%
    dplyr::mutate(
      locusId = rownames(mat_heatmap),
      cluster = as.factor(cluster_assignment_main[rownames(mat_heatmap)])
    ) %>%
    dplyr::left_join(
      fitness123_annotated_com %>%
        dplyr::distinct(locusId, .keep_all = TRUE) %>%
        dplyr::select(locusId, gene, Name),
      by = "locusId"
    )
  
  utils::write.csv(
    df_tsne,
    outfile("tSNE_coordinates.csv"),
    row.names = FALSE
  )
  
  custom_group_colors <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999"
  )
  
  cluster_colors <- custom_group_colors[seq_len(nlevels(df_tsne$cluster))]
  
  # Keep the old panel style and old axis orientation: V1 ~ V2
  plot_tsne <- lattice::xyplot(
    V1 ~ V2,
    data = df_tsne,
    groups = cluster,
    col = cluster_colors,
    pch = 16,
    xlab = "tSNE 1",
    ylab = "tSNE 2",
    par.settings = list(
      superpose.symbol = list(
        col = cluster_colors,
        pch = 16,
        alpha = 0.8
      ),
      strip.background = list(col = "transparent")
    ),
    auto.key = list(
      space = "right",
      title = "Cluster",
      cex.title = 1,
      columns = min(4, length(unique(df_tsne$cluster)))
    ),
    panel = function(x, y, ..., subscripts) {
      lattice::panel.grid(h = -1, v = -1, col = grDevices::grey(0.92), lty = 2)
      lattice::panel.xyplot(x, y, ..., subscripts = subscripts)
      
      gene_labels <- df_tsne$gene[subscripts]
      non_na <- !is.na(gene_labels) & gene_labels != ""
      
      if (any(non_na)) {
        lattice::panel.text(
          x = x[non_na],
          y = y[non_na],
          labels = gene_labels[non_na],
          pos = 3,
          cex = 0.7,
          offset = 0.5,
          col = "black"
        )
      }
    }
  )
  
  grDevices::pdf(outfile("tSNE_plot_with_labels.pdf"), width = 10, height = 8)
  print(plot_tsne)
  grDevices::dev.off()
}


# ---- 7) Figure: KEGG enrichment plots (old workflow reproduction) -----------
load_packages(
  c("limma", "ggplot2", "scales", "forcats", "stringr", "dplyr", "tibble"),
  "KEGG enrichment"
)

if (!file.exists(infile("cga009_vs_consensus_parsed.csv"))) {
  save_placeholder_plot(
    filename = outfile("kegg_enrichment_plot.svg"),
    message = "Mapping file cga009_vs_consensus_parsed.csv not found.",
    width = 10,
    height = 8
  )
  save_placeholder_plot(
    filename = outfile("kegg_bubble_plot.svg"),
    message = "Mapping file cga009_vs_consensus_parsed.csv not found.",
    width = 10,
    height = 8
  )
  save_placeholder_plot(
    filename = outfile("kegg_enrichment_bubble.png"),
    message = "Mapping file cga009_vs_consensus_parsed.csv not found.",
    width = 10,
    height = 6
  )
} else {
  mapping_table <- utils::read.csv(
    infile("cga009_vs_consensus_parsed.csv"),
    stringsAsFactors = FALSE
  )
  
  if (!all(c("subject", "rpa") %in% names(mapping_table))) {
    save_placeholder_plot(
      filename = outfile("kegg_enrichment_plot.svg"),
      message = "Mapping table must include columns: subject and rpa.",
      width = 10,
      height = 8
    )
    save_placeholder_plot(
      filename = outfile("kegg_bubble_plot.svg"),
      message = "Mapping table must include columns: subject and rpa.",
      width = 10,
      height = 8
    )
    save_placeholder_plot(
      filename = outfile("kegg_enrichment_bubble.png"),
      message = "Mapping table must include columns: subject and rpa.",
      width = 10,
      height = 6
    )
  } else {
    # Reproduce the old mapping-based KEGG workflow
    fitness123_rich <- fitness123_annotated_com %>%
      dplyr::left_join(mapping_table, by = c("locusId" = "subject")) %>%
      dplyr::filter(!is.na(rpa), rpa != "")
    
    data_rich <- fitness123_rich %>%
      dplyr::filter(
        dplyr::if_any(
          dplyr::all_of(needed_time_cols),
          ~ !dplyr::between(., -1.5, 0.5)
        )
      ) %>%
      dplyr::select(rpa, dplyr::all_of(needed_time_cols)) %>%
      dplyr::distinct(rpa, .keep_all = TRUE) %>%
      tibble::column_to_rownames("rpa") %>%
      as.matrix()
    
    data_rich <- data_rich[stats::complete.cases(data_rich), , drop = FALSE]
    
    if (nrow(data_rich) < 2) {
      save_placeholder_plot(
        filename = outfile("kegg_enrichment_plot.svg"),
        message = "Not enough mapped genes for KEGG enrichment.",
        width = 10,
        height = 8
      )
      save_placeholder_plot(
        filename = outfile("kegg_bubble_plot.svg"),
        message = "Not enough mapped genes for KEGG enrichment.",
        width = 10,
        height = 8
      )
      save_placeholder_plot(
        filename = outfile("kegg_enrichment_bubble.png"),
        message = "Not enough mapped genes for KEGG enrichment.",
        width = 10,
        height = 6
      )
    } else {
      mat_cluster_rich <- stats::hclust(stats::dist(data_rich), method = "ward.D2")
      k_rich <- min(6, max(2, nrow(data_rich) - 1))
      
      # Old workflow compatibility:
      # order the cutree labels by dendrogram leaf order
      cutreeord <- function(hc, k) {
        cl <- stats::cutree(hc, k = k)
        ord <- hc$order
        cl_ord <- cl[ord]
        
        relabeled <- integer(length(cl_ord))
        seen <- character(0)
        next_id <- 0L
        
        for (i in seq_along(cl_ord)) {
          lab <- as.character(cl_ord[i])
          if (!lab %in% seen) {
            seen <- c(seen, lab)
            next_id <- next_id + 1L
          }
          relabeled[i] <- match(lab, seen)
        }
        
        out <- relabeled
        names(out) <- names(cl_ord)
        
        out[names(cl)]
      }
      
      cluster_assignment_rich <- cutreeord(mat_cluster_rich, k = k_rich)
      
      df_cluster_rich <- tibble::enframe(
        cluster_assignment_rich,
        name = "locus_tag",
        value = "cluster"
      ) %>%
        dplyr::arrange(cluster)
      
      enrichment_list <- lapply(seq_len(k_rich), function(clust) {
        cluster_genes <- df_cluster_rich %>%
          dplyr::filter(cluster == clust) %>%
          dplyr::pull(locus_tag)
        
        if (length(cluster_genes) == 0) {
          return(NULL)
        }
        
        out <- tryCatch(
          limma::kegga(cluster_genes, species.KEGG = "rpa"),
          error = function(e) NULL
        )
        
        if (is.null(out) || nrow(out) == 0) {
          return(NULL)
        }
        
        out_df <- as.data.frame(out)
        
        if ("Pathway" %in% names(out_df)) {
          out_df$Pathway <- NULL
        }
        
        tibble::rownames_to_column(out_df, var = "Pathway") %>%
          dplyr::mutate(cluster = clust)
      })
      
      df_kegg_enrichment <- dplyr::bind_rows(enrichment_list)
      
      if (
        nrow(df_kegg_enrichment) == 0 ||
        !all(c("P.DE", "N", "DE", "Pathway", "cluster") %in% names(df_kegg_enrichment))
      ) {
        save_placeholder_plot(
          filename = outfile("kegg_enrichment_plot.svg"),
          message = "No KEGG enrichment result was produced.",
          width = 10,
          height = 8
        )
        save_placeholder_plot(
          filename = outfile("kegg_bubble_plot.svg"),
          message = "No KEGG enrichment result was produced.",
          width = 10,
          height = 8
        )
        save_placeholder_plot(
          filename = outfile("kegg_enrichment_bubble.png"),
          message = "No KEGG enrichment result was produced.",
          width = 10,
          height = 6
        )
      } else {
        plot_keggenrich <- df_kegg_enrichment %>%
          dplyr::filter(is.finite(P.DE), P.DE > 0, N >= 5, N <= 200) %>%
          dplyr::mutate(
            log10_p_value = log10(P.DE),
            pathway_short = Pathway %>%
              tolower() %>%
              stringr::str_remove_all(" - cupriavidus necator h16|path:") %>%
              stringr::str_remove_all(".?biosynthesis| of| metabolism") %>%
              stringr::str_sub(1, 40)
          ) %>%
          dplyr::filter(log10_p_value < 0) %>%
          dplyr::group_by(cluster) %>%
          dplyr::slice_min(
            order_by = log10_p_value,
            n = 3,
            with_ties = FALSE
          ) %>%
          dplyr::ungroup()
        
        if (nrow(plot_keggenrich) == 0) {
          save_placeholder_plot(
            filename = outfile("kegg_enrichment_plot.svg"),
            message = "No significant KEGG pathway passed the filters.",
            width = 10,
            height = 8
          )
          save_placeholder_plot(
            filename = outfile("kegg_bubble_plot.svg"),
            message = "No significant KEGG pathway passed the filters.",
            width = 10,
            height = 8
          )
          save_placeholder_plot(
            filename = outfile("kegg_enrichment_bubble.png"),
            message = "No significant KEGG pathway passed the filters.",
            width = 10,
            height = 6
          )
        } else {
          p_enrichment <- ggplot2::ggplot(
            plot_keggenrich,
            ggplot2::aes(
              x = log10_p_value,
              y = stats::reorder(pathway_short, log10_p_value),
              color = as.factor(cluster),
              size = DE
            )
          ) +
            ggplot2::geom_point(alpha = 0.85) +
            ggplot2::theme_bw(base_size = 11) +
            ggplot2::labs(
              x = expression("log"[10] * " p-value"),
              y = "Pathway",
              color = "Cluster",
              size = "Gene count",
              title = "KEGG enrichment (top pathways per cluster)"
            )
          
          save_ggplot_svg(
            p_enrichment,
            outfile("kegg_enrichment_plot.svg"),
            width = 10,
            height = 8
          )
          
          bubble_data <- plot_keggenrich %>%
            dplyr::mutate(
              EnrichmentScore = -log10(P.DE),
              cluster = paste("Cluster", cluster)
            )

          p_bubble <- ggplot2::ggplot(
            bubble_data,
            ggplot2::aes(
              x = EnrichmentScore,
              y = stats::reorder(pathway_short, EnrichmentScore)
            )
          ) +
            ggplot2::geom_point(
              ggplot2::aes(size = DE, color = P.DE),
              alpha = 0.85
            ) +
            ggplot2::scale_color_gradient(
              name = "P-value",
              low = "red",
              high = "blue",
              trans = "log10",
              breaks = scales::trans_breaks("log10", function(x) 10^x),
              labels = scales::trans_format("log10", scales::math_format(10^.x))
            ) +
            ggplot2::theme_bw(base_size = 12) +
            ggplot2::facet_grid(cluster ~ ., scales = "free_y", space = "free") +
            ggplot2::labs(
              x = "Enrichment Score (-log10 p-value)",
              y = "Pathway",
              size = "Gene count",
              title = "KEGG pathway enrichment bubble plot"
            )
          
          save_ggplot_svg(
            p_bubble,
            outfile("kegg_bubble_plot.svg"),
            width = 10,
            height = 8
          )
          save_ggplot_png(
            p_bubble,
            outfile("kegg_enrichment_bubble.png"),
            width = 10,
            height = 6,
            dpi = 300
          )
        }
      }
    }
  }
}

# ---- 8) Define focused locus groups (used by heatmap.R) --------------------
group1_locus <- c(
  "PGIDNB_23480", "PGIDNB_23485", "PGIDNB_23475", "PGIDNB_23470", "PGIDNB_23465", "PGIDNB_23490", "PGIDNB_23460"
)

group2_locus <- c(
  "PGIDNB_03300", "PGIDNB_03315", "PGIDNB_03325", "PGIDNB_03330", "PGIDNB_06415",
  "PGIDNB_06420", "PGIDNB_06425", "PGIDNB_11355", "PGIDNB_13255", "PGIDNB_14840",
  "PGIDNB_14930", "PGIDNB_15395", "PGIDNB_16830", "PGIDNB_18030", "PGIDNB_18035",
  "PGIDNB_18040", "PGIDNB_18045", "PGIDNB_21305", "PGIDNB_21310", "PGIDNB_21315",
  "PGIDNB_21320", "PGIDNB_21325", "PGIDNB_22485", "PGIDNB_22490", "PGIDNB_22495",
  "PGIDNB_22500", "PGIDNB_22505", "PGIDNB_23825", "PGIDNB_23830", "PGIDNB_23835",
  "PGIDNB_23840", "PGIDNB_23845"
)

group3_locus <- c(
  "PGIDNB_02970", "PGIDNB_02975", "PGIDNB_03840", "PGIDNB_04125",
  "PGIDNB_06315", "PGIDNB_09720", "PGIDNB_11980", "PGIDNB_12075", "PGIDNB_19470",
  "PGIDNB_20005", "PGIDNB_20010", "PGIDNB_21725", "PGIDNB_22540", "PGIDNB_22545",
  "PGIDNB_22560", "PGIDNB_22565", "PGIDNB_22800", "PGIDNB_22880", "PGIDNB_22895",
  "PGIDNB_22900", "PGIDNB_22905", "PGIDNB_22910", "PGIDNB_22915", "PGIDNB_23600",
  "PGIDNB_23605", "PGIDNB_23610", "PGIDNB_23615"
)

group4_locus <- c(
  "PGIDNB_01310", "PGIDNB_03365", "PGIDNB_05575", "PGIDNB_05580", "PGIDNB_06090",
  "PGIDNB_06095", "PGIDNB_06165", "PGIDNB_06170", "PGIDNB_06175", "PGIDNB_06180", "PGIDNB_06195",
  "PGIDNB_06200", "PGIDNB_06230", "PGIDNB_06235", "PGIDNB_06240", "PGIDNB_06245",
  "PGIDNB_06250", "PGIDNB_06255", "PGIDNB_06260", "PGIDNB_06265", "PGIDNB_06270",
  "PGIDNB_06275", "PGIDNB_06280", "PGIDNB_06290", "PGIDNB_06300", "PGIDNB_06340",
  "PGIDNB_06345", "PGIDNB_06350", "PGIDNB_06355", "PGIDNB_06360", "PGIDNB_06365",
  "PGIDNB_06375", "PGIDNB_06380", "PGIDNB_06385", "PGIDNB_06395", "PGIDNB_06400",
  "PGIDNB_06405", "PGIDNB_11620", "PGIDNB_12445", "PGIDNB_12450", "PGIDNB_13875",
  "PGIDNB_13880", "PGIDNB_13910", "PGIDNB_13915", "PGIDNB_19090", "PGIDNB_19095"
)

# ---- 9) Category table for grouped heatmap ---------------------------------
load_packages(c("dplyr", "tidyr", "ComplexHeatmap", "circlize", "grid"), "category grouped heatmap")

data_filtered <- fitness123_annotated_com %>%
  dplyr::filter(dplyr::if_any(dplyr::all_of(needed_time_cols), ~ !dplyr::between(., -1.5, 0.5)))

fitness123_name <- fitness123_annotated_com %>%
  dplyr::mutate(
    category = dplyr::case_when(
      locusId %in% group1_locus ~ "FDH",
      locusId %in% group2_locus ~ "CBB",
      locusId %in% group3_locus ~ "ETC",
      locusId %in% group4_locus ~ "RC",
      locusId %in% data_filtered$locusId ~ "Other",
      TRUE ~ NA_character_
    )
  )

utils::write.csv(fitness123_name, outfile("fitness123_with_category.csv"), row.names = FALSE)

# Optional subsystem parsing block (kept for traceability)
if (file.exists(infile("reaction_data(dsm123).csv"))) {
  gene_subsystem_df <- utils::read.csv(infile("reaction_data(dsm123).csv"), stringsAsFactors = FALSE) %>%
    dplyr::select(gene_reaction_rule, SUBSYSTEM)

  gene_mapping <- gene_subsystem_df %>%
    tidyr::separate_rows(gene_reaction_rule, sep = "or|and") %>%
    dplyr::mutate(gene_reaction_rule = trimws(gene_reaction_rule)) %>%
    dplyr::mutate(gene_reaction_rule = stringr::str_remove_all(gene_reaction_rule, "\\(|\\)")) %>%
    dplyr::distinct(gene_reaction_rule, SUBSYSTEM) %>%
    dplyr::rename(gene = gene_reaction_rule)

  result_data <- fitness123_name %>%
    dplyr::left_join(gene_mapping, by = c("locusId" = "gene")) %>%
    dplyr::group_by(locusId) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}

filtered_fitness123_name <- fitness123_name %>%
  dplyr::filter(!is.na(category))

if (nrow(filtered_fitness123_name) >= 2) {
  mat_heatmap_category <- filtered_fitness123_name %>%
    dplyr::select(locusId, dplyr::all_of(needed_time_cols)) %>%
    tibble::column_to_rownames("locusId") %>%
    as.matrix()

  mat_heatmap_category <- mat_heatmap_category[stats::complete.cases(mat_heatmap_category), , drop = FALSE]

  transposed_category_matrix <- t(mat_heatmap_category)

  if (ncol(transposed_category_matrix) >= 2) {
    gene_cluster_columns <- stats::hclust(stats::dist(t(transposed_category_matrix)), method = "ward.D2")
    gene_dendrogram <- as.dendrogram(gene_cluster_columns)

    grDevices::png(outfile("gene_cluster_dendrogram.png"), width = 10, height = 6, units = "in", res = 150)
    graphics::plot(gene_dendrogram)
    grDevices::dev.off()
  } else {
    save_placeholder_plot(
      filename = outfile("gene_cluster_dendrogram.png"),
      message = "Not enough genes to build a dendrogram.",
      width = 10,
      height = 6
    )
    gene_dendrogram <- FALSE
  }

  column_split <- filtered_fitness123_name$category[
    match(colnames(transposed_category_matrix), filtered_fitness123_name$locusId)
  ]

  column_split <- factor(
    column_split,
    levels = c("FDH", "ETC", "CBB", "RC", "Other")
  )

  grouped_heatmap <- ComplexHeatmap::Heatmap(
    matrix = transposed_category_matrix,
    name = "Fitness",
    col = circlize::colorRamp2(
      c(-6, -3, 0, 1, 2),
      c("#002F70", "#829EDA", "#E2E2E2", "#E6BCC3", "#E495A5")
    ),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    column_split = column_split,
    show_parent_dend_line = FALSE,
    column_names_side = "bottom",
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 6, fontface = "italic"),
    row_names_side = "right",
    row_names_gp = grid::gpar(fontsize = 10),
    heatmap_height = grid::unit(4, "cm"),
    width = grid::unit(20, "cm"),
    show_heatmap_legend = TRUE
  )

  grDevices::svg(outfile("heatmap(22).svg"), width = 20, height = 5)
  ComplexHeatmap::draw(grouped_heatmap)
  grDevices::dev.off()
} else {
  save_placeholder_plot(
    filename = outfile("gene_cluster_dendrogram.png"),
    message = "Not enough rows in category table.",
    width = 10,
    height = 6
  )
  save_placeholder_plot(
    filename = outfile("heatmap(22).svg"),
    message = "Not enough rows in category table.",
    width = 20,
    height = 5
  )
}

message("Tnseq_heatmap.R completed.")

