## GSE213621 bulk liver RPS27L validation across fibrosis groups.

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
} else {
  getwd()
}
source(file.path(script_dir, "00_setup_paths_and_functions.R"))
load_required_packages(c("data.table", "dplyr", "readr", "stringr", "tidyr"))

gse_dir <- file.path(paths$raw, "GSE213621")
expr_path <- file.path(gse_dir, "GSE213621_FPKMs_allsamples.txt.gz")
series_path <- file.path(gse_dir, "GSE213621_series_matrix.txt.gz")
stopifnot(file.exists(expr_path), file.exists(series_path))

series_lines <- readLines(gzfile(series_path), warn = FALSE)
sample_title <- extract_quoted(get_geo_row(series_lines, "!Sample_title"))
sample_gsm <- extract_quoted(get_geo_row(series_lines, "!Sample_geo_accession"))
cell_type_line <- series_lines[grep('^!Sample_characteristics_ch1.*cell type:', series_lines)[1]]
fibrosis_line <- series_lines[grep('^!Sample_characteristics_ch1.*fibrotic stage:', series_lines)[1]]
cell_type <- sub("^cell type: ", "", extract_quoted(cell_type_line))
fibrotic_stage <- sub("^fibrotic stage: ", "", extract_quoted(fibrosis_line))

sample_meta <- tibble::tibble(
  geo_accession = sample_gsm,
  sample_title = sample_title,
  sample_id = stringr::str_replace(sample_title, "^human liver sample ", ""),
  sample_col = paste0("Sample_", sample_id),
  cell_type = cell_type,
  fibrotic_stage_raw = fibrotic_stage,
  fibrosis_group_external = factor(
    dplyr::recode(fibrotic_stage_raw, F0F1 = "F0-1", F3F4 = "F3-4"),
    levels = c("Control", "F0-1", "F2", "F3-4")
  ),
  fibrosis_midpoint = dplyr::case_when(
    fibrotic_stage_raw == "Control" ~ 0,
    fibrotic_stage_raw == "F0F1" ~ 0.5,
    fibrotic_stage_raw == "F2" ~ 2,
    fibrotic_stage_raw == "F3F4" ~ 3.5,
    TRUE ~ NA_real_
  )
)

expr <- data.table::fread(expr_path)
if (!"V1" %in% names(expr)) {
  names(expr)[1] <- "V1"
}
rps <- expr[expr$V1 == "RPS27L", , drop = FALSE]
if (nrow(rps) != 1) stop("Expected exactly one RPS27L row in ", expr_path)

rps_df <- tibble::tibble(
  sample_col = names(rps)[-1],
  RPS27L_FPKM = as.numeric(rps[1, -1]),
  sample_id = stringr::str_replace(sample_col, "^Sample_", "")
) |>
  dplyr::mutate(RPS27L_log2_FPKM = log2(RPS27L_FPKM + 1)) |>
  dplyr::left_join(sample_meta, by = c("sample_col", "sample_id")) |>
  dplyr::filter(!is.na(fibrosis_group_external), is.finite(RPS27L_log2_FPKM))

group_summary <- rps_df |>
  dplyr::group_by(fibrosis_group_external) |>
  dplyr::summarise(
    n_samples = dplyr::n(),
    median_RPS27L_log2_FPKM = median(RPS27L_log2_FPKM, na.rm = TRUE),
    mean_RPS27L_log2_FPKM = mean(RPS27L_log2_FPKM, na.rm = TRUE),
    .groups = "drop"
  )

kw <- suppressWarnings(stats::kruskal.test(RPS27L_log2_FPKM ~ fibrosis_group_external, data = rps_df))
sp_all <- suppressWarnings(stats::cor.test(rps_df$RPS27L_log2_FPKM, rps_df$fibrosis_midpoint, method = "spearman", exact = FALSE))
fib_only <- rps_df |>
  dplyr::filter(fibrosis_group_external != "Control")
sp_fib <- suppressWarnings(stats::cor.test(fib_only$RPS27L_log2_FPKM, fib_only$fibrosis_midpoint, method = "spearman", exact = FALSE))
overall_stats <- tibble::tibble(
  test = c("Kruskal-Wallis across external fibrosis groups", "Spearman trend including Control", "Spearman trend within NAFLD fibrosis groups"),
  statistic = c(unname(kw$statistic), unname(sp_all$statistic), unname(sp_fib$statistic)),
  estimate = c(NA_real_, unname(sp_all$estimate), unname(sp_fib$estimate)),
  p_value = c(kw$p.value, sp_all$p.value, sp_fib$p.value),
  n = c(nrow(rps_df), nrow(rps_df), nrow(fib_only)),
  significance = p_to_stars(p_value)
)
dunn_tbl <- dunn_bh(rps_df, "RPS27L_log2_FPKM", "fibrosis_group_external", "Bulk liver RPS27L (log2 FPKM)")

readr::write_tsv(rps_df, file.path(paths$tables, "Fig1j_GSE213621_external_bulk_RPS27L_source_data.tsv"))
readr::write_tsv(group_summary, file.path(paths$tables, "Fig1j_GSE213621_external_bulk_RPS27L_group_summary.tsv"))
readr::write_tsv(overall_stats, file.path(paths$tables, "Fig1j_GSE213621_external_bulk_RPS27L_overall_stats.tsv"))
readr::write_tsv(dunn_tbl, file.path(paths$tables, "Fig1j_GSE213621_external_bulk_RPS27L_Dunn_BH.tsv"))
readr::write_tsv(sample_meta, file.path(paths$tables, "GSE213621_GEO_sample_metadata.tsv"))

message("GSE213621 external bulk RPS27L validation tables written to results/tables.")
