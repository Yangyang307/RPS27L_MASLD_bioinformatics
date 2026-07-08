## GSE202379 SAF activity validation using donor-level hepatocyte MA-high summaries.

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
} else {
  getwd()
}
source(file.path(script_dir, "00_setup_paths_and_functions.R"))
load_required_packages(c("dplyr", "readr", "stringr", "tidyr"))

candidate_paths <- c(
  file.path(paths$processed, "GSE202379_donor_level_summary.tsv"),
  file.path(paths$raw, "GSE202379", "donor_level_summary.tsv")
)
donor_path <- candidate_paths[file.exists(candidate_paths)][1]
if (is.na(donor_path)) {
  stop(
    "Missing GSE202379 donor-level summary. Expected one of:\n",
    paste(candidate_paths, collapse = "\n")
  )
}

gse202379_donor <- readr::read_tsv(donor_path, show_col_types = FALSE)
required_cols <- c("patient_id", "disease_clean", "saf_score_raw", "ma_high_pct", "metabolic_aging_score")
missing_cols <- setdiff(required_cols, names(gse202379_donor))
if (length(missing_cols) > 0) {
  stop("GSE202379 donor table is missing required columns: ", paste(missing_cols, collapse = ", "))
}
if (!"fibrosis_stage" %in% names(gse202379_donor)) {
  gse202379_donor$fibrosis_stage <- NA_character_
}

gse202379_donor <- gse202379_donor |>
  dplyr::mutate(
    activity_score = suppressWarnings(as.integer(stringr::str_match(saf_score_raw, "A([0-4])")[, 2])),
    activity_pair = dplyr::case_when(
      activity_score %in% c(0, 1) ~ "A0-1",
      activity_score %in% c(2, 3, 4) ~ "A2-4",
      TRUE ~ NA_character_
    ),
    is_healthy_control = stringr::str_detect(disease_clean, stringr::regex("healthy|control", ignore_case = TRUE)),
    activity_full = dplyr::case_when(
      activity_score %in% 0:4 ~ paste0("A", activity_score),
      TRUE ~ NA_character_
    )
  )

fig1d_excluded_healthy <- gse202379_donor |>
  dplyr::filter(activity_pair %in% c("A0-1", "A2-4"), is_healthy_control, is.finite(ma_high_pct))

fig1d_df <- gse202379_donor |>
  dplyr::filter(activity_pair %in% c("A0-1", "A2-4"), !is_healthy_control, is.finite(ma_high_pct)) |>
  dplyr::mutate(activity_pair = factor(activity_pair, levels = c("A0-1", "A2-4")))

fig1d_w <- suppressWarnings(stats::wilcox.test(ma_high_pct ~ activity_pair, data = fig1d_df, exact = FALSE))
fig1d_stats <- tibble::tibble(
  metric = "MA-high hepatocytes (%)",
  comparison = "SAF activity A0-1 vs A2-4, healthy/control excluded",
  n_A0_1 = sum(fig1d_df$activity_pair == "A0-1"),
  n_A2_4 = sum(fig1d_df$activity_pair == "A2-4"),
  statistic = unname(fig1d_w$statistic),
  p_value = fig1d_w$p.value,
  significance = p_to_stars(fig1d_w$p.value)
)

sup1f_df <- gse202379_donor |>
  dplyr::filter(activity_score %in% 0:4, is.finite(ma_high_pct)) |>
  dplyr::mutate(
    activity_full = factor(paste0("A", activity_score), levels = c("A0", "A1", "A2", "A3", "A4")),
    sample_origin = dplyr::if_else(is_healthy_control, "Healthy/control assigned by A score", "NAFLD/NASH")
  )
sup1f_kw <- suppressWarnings(stats::kruskal.test(ma_high_pct ~ activity_full, data = sup1f_df))
sup1f_dunn <- dunn_bh(sup1f_df, "ma_high_pct", "activity_full", "MA-high hepatocytes (%)")
sup1f_group_summary <- sup1f_df |>
  dplyr::group_by(activity_full) |>
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_healthy_control_assigned = sum(sample_origin == "Healthy/control assigned by A score"),
    median_ma_high_pct = median(ma_high_pct, na.rm = TRUE),
    mean_ma_high_pct = mean(ma_high_pct, na.rm = TRUE),
    .groups = "drop"
  )
sup1f_kw_tbl <- tibble::tibble(
  metric = "MA-high hepatocytes (%)",
  test = "Kruskal-Wallis with healthy/control samples assigned by parsed A score",
  statistic = unname(sup1f_kw$statistic),
  p_value = sup1f_kw$p.value,
  significance = p_to_stars(sup1f_kw$p.value)
)

readr::write_tsv(
  fig1d_excluded_healthy |>
    dplyr::select(patient_id, activity_pair, activity_full, disease_clean, saf_score_raw, ma_high_pct, metabolic_aging_score, fibrosis_stage),
  file.path(paths$tables, "Fig1d_GSE202379_excluded_Healthy_control_with_SAF_scores.tsv")
)
readr::write_tsv(
  fig1d_df |>
    dplyr::select(patient_id, activity_pair, activity_full, disease_clean, saf_score_raw, ma_high_pct, metabolic_aging_score, fibrosis_stage),
  file.path(paths$tables, "Fig1d_GSE202379_SAF_A01_A24_MA_high_percent_source_data.tsv")
)
readr::write_tsv(fig1d_stats, file.path(paths$tables, "Fig1d_GSE202379_SAF_A01_A24_Wilcoxon.tsv"))
readr::write_tsv(
  sup1f_df |>
    dplyr::select(patient_id, activity_full, sample_origin, disease_clean, saf_score_raw, ma_high_pct, metabolic_aging_score, fibrosis_stage),
  file.path(paths$tables, "SupFig1f_GSE202379_full_SAF_MA_high_percent_source_data.tsv")
)
readr::write_tsv(sup1f_group_summary, file.path(paths$tables, "SupFig1f_GSE202379_full_SAF_group_summary.tsv"))
readr::write_tsv(
  sup1f_df |>
    dplyr::filter(sample_origin == "Healthy/control assigned by A score") |>
    dplyr::select(patient_id, activity_full, disease_clean, saf_score_raw, ma_high_pct, metabolic_aging_score, fibrosis_stage),
  file.path(paths$tables, "SupFig1f_GSE202379_Healthy_control_assigned_to_A_groups.tsv")
)
readr::write_tsv(sup1f_kw_tbl, file.path(paths$tables, "SupFig1f_GSE202379_full_SAF_Kruskal.tsv"))
readr::write_tsv(sup1f_dunn, file.path(paths$tables, "SupFig1f_GSE202379_full_SAF_Dunn_BH.tsv"))

message("GSE202379 SAF validation tables written to results/tables.")
