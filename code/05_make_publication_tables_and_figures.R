## Generate publication-style tables and figures from analysis result tables.

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
} else {
  getwd()
}
source(file.path(script_dir, "00_setup_paths_and_functions.R"))
load_required_packages(c("data.table", "dplyr", "readr", "stringr", "tidyr", "ggplot2", "patchwork"))

pal_disease <- c("No-MASLD" = "#D2D2D2", "MASL" = "#FFAF5A", "eMASH" = "#87EBFA", "aMASH" = "#FFCDCD")
pal_pair <- c("A0-1" = "#87EBFA", "A2-4" = "#FFCDCD")
pal_saf_full <- c("A0" = "#C7F6FB", "A1" = "#87EBFA", "A2" = "#FFE3E3", "A3" = "#FFCDCD", "A4" = "#F5A7B2")
pal_state <- c("MA-low" = "#87EBFA", "MA-high" = "#FFCDCD")
pal_bulk <- c("Control" = "#D2D2D2", "F0-1" = "#FFAF5A", "F2" = "#87EBFA", "F3-4" = "#FFCDCD")
ink <- "#202020"
mid <- "#666666"
group_levels <- names(pal_disease)

read_table <- function(name, required = TRUE) {
  path <- file.path(paths$tables, name)
  if (!file.exists(path)) {
    if (required) stop("Missing result table: ", path)
    warning("Skipping missing optional table: ", path)
    return(NULL)
  }
  readr::read_tsv(path, show_col_types = FALSE)
}

add_simple_bracket <- function(plot, x1, x2, y, label, span, show_ns = FALSE) {
  if (length(label) == 0 || is.na(label) || (label == "ns" && !show_ns)) return(plot)
  plot +
    ggplot2::annotate("segment", x = x1, xend = x1, y = y - span * 0.04, yend = y, linewidth = 0.42) +
    ggplot2::annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.42) +
    ggplot2::annotate("segment", x = x2, xend = x2, y = y - span * 0.04, yend = y, linewidth = 0.42) +
    ggplot2::annotate("text", x = (x1 + x2) / 2, y = y + span * 0.03, label = label, size = 3.0, fontface = "bold", family = "Helvetica")
}

make_box <- function(df, x_col, y_col, fill_values, x_lab, y_lab, title) {
  plot_df <- df |>
    dplyr::filter(is.finite(.data[[y_col]]), !is.na(.data[[x_col]]))
  ggplot2::ggplot(plot_df, ggplot2::aes(.data[[x_col]], .data[[y_col]], fill = .data[[x_col]])) +
    ggplot2::geom_boxplot(width = 0.52, linewidth = 0.42, outlier.shape = NA, colour = ink, alpha = 0.72) +
    ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.28, linewidth = 0.44, colour = "black") +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.09, height = 0, seed = 20260705),
      shape = 16, size = 1.35, colour = mid, alpha = 0.50
    ) +
    ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
    ggplot2::labs(x = x_lab, y = y_lab, title = title) +
    theme_publication(8)
}

donor_summary <- read_table("GSE289173_donor_metabolic_aging_summary.tsv", required = FALSE)
if (!is.null(donor_summary)) {
  donor_summary <- donor_summary |>
    dplyr::mutate(disease_group = factor(disease_group, levels = group_levels))
  p_fig1c <- make_box(
    donor_summary,
    "disease_group",
    "ma_high_pct",
    pal_disease,
    NULL,
    "MA-high hepatocytes (%)",
    "MA-high hepatocytes"
  )
  save_plot(p_fig1c, file.path(paths$figures, "Fig1c_GSE289173_MA_high_percent_disease_grade"), 2.9, 2.9)

  p_sup1e <- make_box(
    donor_summary,
    "disease_group",
    "metabolic_aging_score",
    pal_disease,
    NULL,
    "Metabolic-aging score",
    "Metabolic-aging score"
  )
  save_plot(p_sup1e, file.path(paths$figures, "SupFig1e_GSE289173_metabolic_aging_score_disease_grade"), 2.9, 2.9)
}

fig1d_df <- read_table("Fig1d_GSE202379_SAF_A01_A24_MA_high_percent_source_data.tsv", required = FALSE)
fig1d_stats <- read_table("Fig1d_GSE202379_SAF_A01_A24_Wilcoxon.tsv", required = FALSE)
if (!is.null(fig1d_df)) {
  fig1d_df <- fig1d_df |>
    dplyr::mutate(activity_pair = factor(activity_pair, levels = c("A0-1", "A2-4")))
  ymax <- max(fig1d_df$ma_high_pct, na.rm = TRUE)
  ymin <- min(fig1d_df$ma_high_pct, na.rm = TRUE)
  span <- max(ymax - ymin, 1)
  p_fig1d <- make_box(fig1d_df, "activity_pair", "ma_high_pct", pal_pair, "SAF activity score", "MA-high hepatocytes (%)", "SAF activity") +
    ggplot2::coord_cartesian(ylim = c(ymin - 0.05 * span, ymax + 0.25 * span), clip = "off")
  if (!is.null(fig1d_stats)) {
    p_fig1d <- add_simple_bracket(p_fig1d, 1, 2, ymax + span * 0.12, fig1d_stats$significance[[1]], span, show_ns = TRUE)
  }
  save_plot(p_fig1d, file.path(paths$figures, "Fig1d_GSE202379_SAF_A01_A24_MA_high_percent"), 2.55, 2.9)
}

gsea_selected <- read_table("Fig1e_GSE289173_selected_Hallmark_Reactome_GSEA_terms.tsv", required = FALSE)
if (!is.null(gsea_selected) && nrow(gsea_selected) > 0) {
  gsea_selected <- gsea_selected |>
    dplyr::mutate(
      term_label = factor(term_label, levels = rev(unique(term_label))),
      neg_log10_fdr = -log10(pmax(p.adjust, 1e-300))
    )
  p_fig1e <- ggplot2::ggplot(gsea_selected, ggplot2::aes(NES, term_label, size = neg_log10_fdr, fill = NES)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey55") +
    ggplot2::geom_point(shape = 21, color = ink, stroke = 0.20, alpha = 0.95) +
    ggplot2::scale_fill_gradient2(low = "#87EBFA", mid = "#F7F7F7", high = "#FFCDCD", midpoint = 0, name = "NES") +
    ggplot2::scale_size_continuous(range = c(1.6, 4.2), name = "-log10 FDR") +
    ggplot2::labs(x = "Normalized enrichment score", y = NULL, title = "MA-high vs MA-low pathways") +
    theme_publication(7.4, legend_position = "right") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7.0), legend.title = ggplot2::element_text(size = 7))
  save_plot(p_fig1e, file.path(paths$figures, "Fig1e_GSE289173_donor_internal_MA_high_vs_low_Hallmark_Reactome_GSEA"), 4.7, 3.25)
}

rps_donor <- read_table("GSE289173_hepatocyte_RPS27L_pseudobulk_by_donor.tsv", required = FALSE)
if (!is.null(rps_donor)) {
  rps_donor <- rps_donor |>
    dplyr::mutate(disease_group = factor(disease_group, levels = group_levels))
  p_fig1f <- make_box(
    rps_donor,
    "disease_group",
    "RPS27L_log1p_cpm",
    pal_disease,
    NULL,
    "Hepatocyte pseudobulk\nRPS27L (log1p CPM)",
    "RPS27L by disease grade"
  )
  save_plot(p_fig1f, file.path(paths$figures, "Fig1f_GSE289173_hepatocyte_RPS27L_disease_grade"), 2.9, 2.9)
}

rps_state <- read_table("Fig1g_GSE289173_donor_internal_MA_high_low_RPS27L_source_data.tsv", required = FALSE)
paired_stats <- read_table("Fig1g_GSE289173_donor_internal_MA_high_low_RPS27L_paired_Wilcoxon.tsv", required = FALSE)
if (!is.null(rps_state)) {
  paired_long <- rps_state |>
    dplyr::mutate(state = factor(state, levels = c("MA-low", "MA-high")))
  ymax <- max(paired_long$RPS27L_log1p_cpm, na.rm = TRUE)
  ymin <- min(paired_long$RPS27L_log1p_cpm, na.rm = TRUE)
  span <- max(ymax - ymin, 1)
  p_fig1g <- ggplot2::ggplot(paired_long, ggplot2::aes(state, RPS27L_log1p_cpm, group = orig.ident)) +
    ggplot2::geom_line(linewidth = 0.25, color = "grey65", alpha = 0.65) +
    ggplot2::geom_point(ggplot2::aes(fill = state), shape = 21, size = 1.8, stroke = 0.25, color = ink, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = pal_state, drop = FALSE) +
    ggplot2::labs(x = NULL, y = "Hepatocyte pseudobulk\nRPS27L (log1p CPM)", title = "RPS27L in MA-low vs MA-high") +
    ggplot2::coord_cartesian(ylim = c(ymin - 0.05 * span, ymax + 0.22 * span), clip = "off") +
    theme_publication(8)
  if (!is.null(paired_stats)) {
    p_fig1g <- add_simple_bracket(p_fig1g, 1, 2, ymax + span * 0.10, paired_stats$significance[[1]], span, show_ns = TRUE)
  }
  save_plot(p_fig1g, file.path(paths$figures, "Fig1g_GSE289173_RPS27L_donor_internal_MA_high_low_paired"), 2.55, 2.9)
}

fig1j_df <- read_table("Fig1j_GSE213621_external_bulk_RPS27L_source_data.tsv", required = FALSE)
if (!is.null(fig1j_df)) {
  fig1j_df <- fig1j_df |>
    dplyr::mutate(fibrosis_group_external = factor(fibrosis_group_external, levels = names(pal_bulk)))
  p_fig1j <- make_box(
    fig1j_df,
    "fibrosis_group_external",
    "RPS27L_log2_FPKM",
    pal_bulk,
    NULL,
    "Bulk liver RPS27L\n(log2 FPKM)",
    "External bulk RPS27L"
  )
  save_plot(p_fig1j, file.path(paths$figures, "Fig1j_GSE213621_external_bulk_RPS27L_by_fibrosis"), 3.05, 2.9)
}

sup1f_df <- read_table("SupFig1f_GSE202379_full_SAF_MA_high_percent_source_data.tsv", required = FALSE)
if (!is.null(sup1f_df)) {
  sup1f_df <- sup1f_df |>
    dplyr::mutate(activity_full = factor(activity_full, levels = names(pal_saf_full)))
  p_sup1f <- make_box(sup1f_df, "activity_full", "ma_high_pct", pal_saf_full, "SAF activity score", "MA-high hepatocytes (%)", "Full SAF activity")
  save_plot(p_sup1f, file.path(paths$figures, "SupFig1f_GSE202379_full_SAF_MA_high_percent"), 3.1, 2.9)
}

umap_path <- file.path(paths$processed, "GSE289173_hepatocyte_umap_coordinates.tsv.gz")
cell_scores_path <- file.path(paths$processed, "GSE289173_hepatocyte_metabolic_aging_cell_scores.tsv.gz")
if (file.exists(umap_path) && file.exists(cell_scores_path)) {
  umap_df <- data.table::fread(cmd = sprintf("gzip -dc %s", shQuote(umap_path))) |>
    tibble::as_tibble() |>
    dplyr::inner_join(
      data.table::fread(cmd = sprintf("gzip -dc %s", shQuote(cell_scores_path))) |>
        tibble::as_tibble() |>
        dplyr::select(cell_barcode, disease_group, ma_state_global),
      by = "cell_barcode"
    ) |>
    dplyr::mutate(
      disease_group = factor(disease_group, levels = group_levels),
      ma_state_umap = factor(dplyr::if_else(ma_state_global == "MA_high", "MA-high", "Other hepatocytes"),
                             levels = c("Other hepatocytes", "MA-high"))
    )
  p_fig1b <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP_1, UMAP_2, color = ma_state_umap)) +
    ggplot2::geom_point(data = dplyr::filter(umap_df, ma_state_umap == "Other hepatocytes"), size = 0.14, alpha = 0.72, stroke = 0) +
    ggplot2::geom_point(data = dplyr::filter(umap_df, ma_state_umap == "MA-high"), size = 0.12, alpha = 0.82, stroke = 0) +
    ggplot2::facet_wrap(~disease_group, nrow = 1) +
    ggplot2::scale_color_manual(values = c("Other hepatocytes" = "#35B7C8", "MA-high" = "#E60012"), drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "snRNA-seq of human hepatocytes") +
    ggplot2::theme_void(base_size = 8, base_family = "Helvetica") +
    ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold"))
  save_plot(p_fig1b, file.path(paths$figures, "Fig1b_GSE289173_hepatocyte_UMAP_MA_high_split"), 7.1, 2.3)
}

message("Publication figures written to results/figures.")
