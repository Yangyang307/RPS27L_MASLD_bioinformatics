## GSE289173 hepatocyte metabolic-aging score and RPS27L pseudobulk analysis.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
} else {
  getwd()
}
source(file.path(script_dir, "00_setup_paths_and_functions.R"))
load_required_packages(c("Matrix", "data.table", "dplyr", "readr", "stringr", "tidyr"))

group_levels <- c("No-MASLD", "MASL", "eMASH", "aMASH")

gse_dir <- file.path(paths$raw, "GSE289173")
mtx_dir <- file.path(gse_dir, "matrix")
barcode_dir <- file.path(gse_dir, "barcodes")
feature_dir <- file.path(gse_dir, "features")

series_path <- file.path(gse_dir, "GSE289173_series_matrix.txt.gz")
cov_path <- file.path(gse_dir, "covariates.hepatocyte.txt.gz")
gene_set_path <- file.path(paths$gene_sets, "metabolic_aging_gene_sets.tsv")

stopifnot(file.exists(series_path), file.exists(cov_path), file.exists(gene_set_path))

series_lines <- readLines(gzfile(series_path), warn = FALSE)
sample_title <- extract_quoted(get_geo_row(series_lines, "!Sample_title"))
sample_gsm <- extract_quoted(get_geo_row(series_lines, "!Sample_geo_accession"))
disease_line <- series_lines[grep('^!Sample_characteristics_ch1.*disease group:', series_lines)[1]]
disease_raw <- sub("^disease group: ", "", extract_quoted(disease_line))
barcode_urls <- extract_quoted(get_geo_row(series_lines, "!Sample_supplementary_file_1"))
barcode_urls <- sub("^ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov", barcode_urls)
matrix_urls <- extract_quoted(get_geo_row(series_lines, "!Sample_supplementary_file_3"))
matrix_urls <- sub("^ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov", matrix_urls)

sample_meta <- tibble::tibble(
  orig.ident = sample_title,
  gsm = sample_gsm,
  disease_group_raw = disease_raw,
  disease_group = dplyr::recode(disease_raw, ctrl = "No-MASLD", MASL = "MASL", eMASH = "eMASH", aMASH = "aMASH"),
  severity_score = match(disease_group, group_levels) - 1,
  barcode_url = barcode_urls,
  matrix_url = matrix_urls,
  barcode_file = file.path(barcode_dir, basename(barcode_urls)),
  matrix_file = file.path(mtx_dir, basename(matrix_urls))
) |>
  dplyr::mutate(disease_group = factor(disease_group, levels = group_levels))

readr::write_tsv(sample_meta, file.path(paths$tables, "GSE289173_GEO_sample_metadata.tsv"))

feature_files <- list.files(feature_dir, pattern = "features[.]tsv[.]gz$", full.names = TRUE)
if (length(feature_files) == 0) {
  stop("No feature file found under ", feature_dir, ". Provide one GEO feature TSV file, for example GSM*_features.tsv.gz.")
}
feature_path <- feature_files[[1]]
features <- data.table::fread(cmd = sprintf("gzip -dc %s", shQuote(feature_path)), header = FALSE)
if (ncol(features) < 2) stop("Feature table must contain gene ID and gene symbol columns.")
data.table::setnames(features, seq_len(ncol(features)), paste0("V", seq_len(ncol(features))))
features <- features |>
  as_tibble() |>
  dplyr::transmute(gene_id = V1, gene = V2)

gene_sets_raw <- readr::read_tsv(gene_set_path, show_col_types = FALSE) |>
  dplyr::mutate(module = as.character(module), gene = as.character(gene))
score_excluded_genes <- "RPS27L"
score_gene_sets <- gene_sets_raw |>
  dplyr::filter(!gene %in% score_excluded_genes)
score_genes <- unique(score_gene_sets$gene)
score_rows <- which(features$gene %in% score_genes)
score_gene_names <- features$gene[score_rows]

availability <- gene_sets_raw |>
  dplyr::group_by(module) |>
  dplyr::summarise(
    requested_genes = dplyr::n_distinct(gene),
    excluded_from_score = paste(intersect(unique(gene), score_excluded_genes), collapse = ";"),
    n_excluded_from_score = length(intersect(unique(gene), score_excluded_genes)),
    available_genes = paste(intersect(unique(gene), unique(features$gene)), collapse = ";"),
    n_available = length(intersect(unique(gene), unique(features$gene))),
    available_scoring_genes = paste(intersect(setdiff(unique(gene), score_excluded_genes), unique(features$gene)), collapse = ";"),
    n_available_scoring_genes = length(intersect(setdiff(unique(gene), score_excluded_genes), unique(features$gene))),
    .groups = "drop"
  )
readr::write_tsv(availability, file.path(paths$tables, "GSE289173_metabolic_aging_gene_set_availability.tsv"))
readr::write_tsv(
  tibble::tibble(excluded_gene = score_excluded_genes, reason = "Excluded before MA score calculation to avoid circular association testing."),
  file.path(paths$tables, "GSE289173_metabolic_aging_score_excluded_genes.tsv")
)

row_by_gene <- split(seq_along(score_gene_names), score_gene_names)
module_rows <- lapply(split(score_gene_sets$gene, score_gene_sets$module), function(genes) {
  rows <- unlist(row_by_gene[intersect(unique(genes), names(row_by_gene))], use.names = FALSE)
  sort(unique(rows))
})
module_rows <- module_rows[vapply(module_rows, length, integer(1)) >= 2]

cov <- data.table::fread(cmd = sprintf("gzip -dc %s", shQuote(cov_path)))
prefix <- paste0(cov$orig.ident, "_", cov$orig.ident, "_")
cov[, raw_barcode := data.table::fifelse(startsWith(cell_barcode, prefix),
                                         substring(cell_barcode, nchar(prefix) + 1),
                                         NA_character_)]
cov <- cov[!is.na(raw_barcode) & orig.ident %in% sample_meta$orig.ident]

cell_score_list <- vector("list", nrow(sample_meta))
sample_qc <- vector("list", nrow(sample_meta))
rps_rows <- list()

for (i in seq_len(nrow(sample_meta))) {
  smp <- sample_meta[i, ]
  message(sprintf("[%02d/%02d] %s", i, nrow(sample_meta), smp$orig.ident))
  download_if_missing(smp$barcode_url, smp$barcode_file)
  download_if_missing(smp$matrix_url, smp$matrix_file)

  barcodes <- read_barcodes(smp$barcode_file)
  cov_s <- cov[orig.ident == smp$orig.ident]
  col_idx <- match(cov_s$raw_barcode, barcodes)
  valid <- !is.na(col_idx)
  cov_s <- cov_s[valid]
  col_idx <- col_idx[valid]

  mat <- read_mtx_gz(smp$matrix_file)
  if (nrow(mat) != nrow(features)) stop("Feature count does not match matrix rows for ", smp$orig.ident)
  mat_s <- mat[, col_idx, drop = FALSE]
  total_counts <- Matrix::colSums(mat_s)
  keep <- is.finite(total_counts) & total_counts > 0
  mat_s <- mat_s[, keep, drop = FALSE]
  cov_s <- cov_s[keep]
  total_counts <- total_counts[keep]

  score_counts <- mat_s[score_rows, , drop = FALSE]
  norm_factor <- 10000 / total_counts
  score_norm <- score_counts %*% Matrix::Diagonal(x = norm_factor)
  score_norm@x <- log1p(score_norm@x)

  scores <- lapply(names(module_rows), function(module) {
    Matrix::colMeans(score_norm[module_rows[[module]], , drop = FALSE])
  })
  names(scores) <- names(module_rows)
  score_dt <- data.table::as.data.table(scores)
  score_dt[, `:=`(
    cell_barcode = cov_s$cell_barcode,
    raw_barcode = cov_s$raw_barcode,
    orig.ident = smp$orig.ident,
    gsm = smp$gsm,
    disease_group_raw = smp$disease_group_raw,
    disease_group = as.character(smp$disease_group),
    severity_score = smp$severity_score
  )]
  data.table::setcolorder(score_dt, c("cell_barcode", "raw_barcode", "orig.ident", "gsm",
                                      "disease_group_raw", "disease_group", "severity_score",
                                      names(module_rows)))
  cell_score_list[[i]] <- score_dt

  rps_row <- which(features$gene == "RPS27L")
  if (length(rps_row) == 1) {
    total_umis <- sum(total_counts)
    rps_counts <- sum(mat_s[rps_row, , drop = TRUE])
    rps_rows[[as.character(smp$orig.ident)]] <- tibble::tibble(
      orig.ident = smp$orig.ident,
      gsm = as.character(smp$gsm),
      disease_group = as.character(smp$disease_group),
      severity_score = smp$severity_score,
      n_hepatocytes = ncol(mat_s),
      total_umis = total_umis,
      RPS27L_counts = rps_counts,
      RPS27L_detected_cells = sum(mat_s[rps_row, , drop = TRUE] > 0),
      RPS27L_positive_pct = 100 * RPS27L_detected_cells / n_hepatocytes,
      RPS27L_cpm = 1e6 * RPS27L_counts / total_umis,
      RPS27L_log1p_cpm = log1p(RPS27L_cpm)
    )
  }

  sample_qc[[i]] <- tibble::tibble(
    orig.ident = smp$orig.ident,
    gsm = smp$gsm,
    disease_group = as.character(smp$disease_group),
    n_matrix_barcodes = length(barcodes),
    n_covariate_hepatocytes = nrow(cov[orig.ident == smp$orig.ident]),
    n_matched_hepatocytes = nrow(score_dt),
    n_missing_covariate_barcodes = sum(!valid)
  )

  rm(mat, mat_s, score_counts, score_norm, score_dt)
  gc(verbose = FALSE)
}

cell_scores <- data.table::rbindlist(cell_score_list, use.names = TRUE, fill = TRUE)
cell_scores[, disease_group := factor(disease_group, levels = group_levels)]

positive_modules <- intersect(c("senescence", "p53_dna_damage", "oxidative_stress", "er_proteostasis", "sasp_inflammatory"), names(cell_scores))
function_modules <- intersect(c("hepatocyte_function", "fatty_acid_oxidation", "xenobiotic_metabolism", "urea_cycle", "mature_hepatocyte_identity"), names(cell_scores))

for (module in c(positive_modules, function_modules)) {
  cell_scores[[paste0(module, "_z")]] <- zscore(cell_scores[[module]])
}
cell_scores[, stress_senescence_score := rowMeans(.SD, na.rm = TRUE), .SDcols = paste0(positive_modules, "_z")]
cell_scores[, hepatocyte_function_score := rowMeans(.SD, na.rm = TRUE), .SDcols = paste0(function_modules, "_z")]
cell_scores[, metabolic_aging_score := stress_senescence_score - hepatocyte_function_score]

q25 <- quantile(cell_scores$metabolic_aging_score, 0.25, na.rm = TRUE)
q75 <- quantile(cell_scores$metabolic_aging_score, 0.75, na.rm = TRUE)
cell_scores[, ma_state_global := data.table::fcase(
  metabolic_aging_score >= q75, "MA_high",
  metabolic_aging_score <= q25, "MA_low",
  default = "MA_mid"
)]

donor_thresholds <- cell_scores |>
  as_tibble() |>
  dplyr::group_by(orig.ident, disease_group, severity_score) |>
  dplyr::summarise(
    donor_q25 = quantile(metabolic_aging_score, 0.25, na.rm = TRUE),
    donor_q75 = quantile(metabolic_aging_score, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

cell_scores <- cell_scores |>
  as_tibble() |>
  dplyr::left_join(donor_thresholds, by = c("orig.ident", "disease_group", "severity_score")) |>
  dplyr::mutate(
    ma_state_donor = dplyr::case_when(
      metabolic_aging_score >= donor_q75 ~ "MA_high",
      metabolic_aging_score <= donor_q25 ~ "MA_low",
      TRUE ~ "MA_mid"
    )
  )

donor_summary <- cell_scores |>
  dplyr::group_by(orig.ident, gsm, disease_group, disease_group_raw, severity_score) |>
  dplyr::summarise(
    n_hepatocytes = dplyr::n(),
    n_ma_high = sum(ma_state_global == "MA_high", na.rm = TRUE),
    ma_high_pct = 100 * n_ma_high / n_hepatocytes,
    metabolic_aging_score = median(metabolic_aging_score, na.rm = TRUE),
    metabolic_aging_score_mean = mean(metabolic_aging_score, na.rm = TRUE),
    stress_senescence_score = median(stress_senescence_score, na.rm = TRUE),
    hepatocyte_function_score = median(hepatocyte_function_score, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(disease_group = factor(disease_group, levels = group_levels))

rps_donor <- dplyr::bind_rows(rps_rows) |>
  dplyr::mutate(disease_group = factor(disease_group, levels = group_levels))

trend_stats <- dplyr::bind_rows(
  kw_spearman(donor_summary, "ma_high_pct", "disease_group", "severity_score", "MA-high hepatocytes (%)"),
  kw_spearman(donor_summary, "metabolic_aging_score", "disease_group", "severity_score", "Metabolic-aging score"),
  kw_spearman(rps_donor, "RPS27L_log1p_cpm", "disease_group", "severity_score", "Hepatocyte pseudobulk RPS27L (log1p CPM)")
)

readr::write_tsv(dplyr::bind_rows(sample_qc), file.path(paths$tables, "GSE289173_sample_qc.tsv"))
readr::write_tsv(donor_summary, file.path(paths$tables, "GSE289173_donor_metabolic_aging_summary.tsv"))
readr::write_tsv(donor_thresholds, file.path(paths$tables, "GSE289173_donor_internal_metabolic_aging_thresholds.tsv"))
readr::write_tsv(tibble::tibble(threshold = c("global_q25", "global_q75"), metabolic_aging_score = c(q25, q75)),
                 file.path(paths$tables, "GSE289173_global_metabolic_aging_thresholds.tsv"))
readr::write_tsv(rps_donor, file.path(paths$tables, "GSE289173_hepatocyte_RPS27L_pseudobulk_by_donor.tsv"))
readr::write_tsv(trend_stats, file.path(paths$tables, "GSE289173_trend_statistics.tsv"))
data.table::fwrite(cell_scores, file.path(paths$processed, "GSE289173_hepatocyte_metabolic_aging_cell_scores.tsv.gz"), sep = "\t")

message("GSE289173 metabolic-aging and RPS27L outputs written to results/tables and data/processed.")
