## Donor-internal MA-high vs MA-low pseudobulk differential expression and GSEA.

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
} else {
  getwd()
}
source(file.path(script_dir, "00_setup_paths_and_functions.R"))
load_required_packages(c("Matrix", "data.table", "dplyr", "readr", "stringr", "tidyr", "edgeR", "clusterProfiler", "msigdbr"))

group_levels <- c("No-MASLD", "MASL", "eMASH", "aMASH")

get_msig <- function(collection, subcollection = NULL) {
  args_new <- list(species = "Homo sapiens", collection = collection)
  args_old <- list(species = "Homo sapiens", category = collection)
  if (!is.null(subcollection)) {
    args_new$subcollection <- subcollection
    args_old$subcategory <- subcollection
  }
  tryCatch(do.call(msigdbr::msigdbr, args_new), error = function(e) do.call(msigdbr::msigdbr, args_old))
}

clean_term <- function(x) {
  x |>
    stringr::str_replace("^GOBP_", "") |>
    stringr::str_replace("^REACTOME_", "") |>
    stringr::str_replace("^HALLMARK_", "") |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_to_sentence()
}

short_pathway_label <- function(id, description) {
  dplyr::case_when(
    stringr::str_detect(id, stringr::regex("TNFA_SIGNALING_VIA_NFKB|NFKB", ignore_case = TRUE)) ~ "TNFA/NFKB signaling",
    stringr::str_detect(id, stringr::regex("P53", ignore_case = TRUE)) ~ "p53 pathway",
    stringr::str_detect(id, stringr::regex("INFLAMMATORY_RESPONSE", ignore_case = TRUE)) ~ "Inflammatory response",
    stringr::str_detect(id, stringr::regex("ATF6|CHAPERON", ignore_case = TRUE)) ~ "ATF6 chaperone response",
    stringr::str_detect(id, stringr::regex("TRANSLATION|RIBOSOM", ignore_case = TRUE)) ~ "Translation / ribosome",
    stringr::str_detect(id, stringr::regex("APOPTOSIS", ignore_case = TRUE)) ~ "Apoptosis",
    stringr::str_detect(id, stringr::regex("FATTY_ACID", ignore_case = TRUE)) ~ "Fatty acid metabolism",
    stringr::str_detect(id, stringr::regex("PHASE_I_FUNCTIONALIZATION|CYTOCHROME_P450|XENOBIOTIC", ignore_case = TRUE)) ~ "Xenobiotic metabolism",
    TRUE ~ clean_term(description)
  )
}

cell_scores_path <- file.path(paths$processed, "GSE289173_hepatocyte_metabolic_aging_cell_scores.tsv.gz")
sample_meta_path <- file.path(paths$tables, "GSE289173_GEO_sample_metadata.tsv")
donor_summary_path <- file.path(paths$tables, "GSE289173_donor_metabolic_aging_summary.tsv")
feature_dir <- file.path(paths$raw, "GSE289173", "features")

stopifnot(file.exists(cell_scores_path), file.exists(sample_meta_path), file.exists(donor_summary_path))

feature_files <- list.files(feature_dir, pattern = "features[.]tsv[.]gz$", full.names = TRUE)
if (length(feature_files) == 0) {
  stop("No GSE289173 feature TSV found under ", feature_dir)
}
features <- data.table::fread(cmd = sprintf("gzip -dc %s", shQuote(feature_files[[1]])), header = FALSE)
data.table::setnames(features, seq_len(ncol(features)), paste0("V", seq_len(ncol(features))))
features <- features |>
  tibble::as_tibble() |>
  dplyr::transmute(gene_id = V1, gene = V2)
rps_row <- which(features$gene == "RPS27L")
if (length(rps_row) < 1) stop("RPS27L was not found in the feature table.")
rps_row <- rps_row[[1]]

sample_meta <- readr::read_tsv(sample_meta_path, show_col_types = FALSE) |>
  dplyr::mutate(disease_group = factor(disease_group, levels = group_levels))
cell_scores <- data.table::fread(cmd = sprintf("gzip -dc %s", shQuote(cell_scores_path))) |>
  tibble::as_tibble() |>
  dplyr::mutate(disease_group = factor(disease_group, levels = group_levels))
donor_summary <- readr::read_tsv(donor_summary_path, show_col_types = FALSE) |>
  dplyr::mutate(disease_group = factor(disease_group, levels = group_levels))

state_rows <- list()
pb_counts <- list()
pb_meta_rows <- list()

for (i in seq_len(nrow(sample_meta))) {
  donor <- sample_meta$orig.ident[[i]]
  if (!file.exists(sample_meta$barcode_file[[i]]) || !file.exists(sample_meta$matrix_file[[i]])) {
    stop("Missing matrix or barcode file for ", donor)
  }
  donor_cells <- cell_scores |>
    dplyr::filter(orig.ident == donor, ma_state_donor %in% c("MA_low", "MA_high")) |>
    dplyr::select(cell_barcode, raw_barcode, ma_state_donor)
  if (nrow(donor_cells) == 0) next

  message(sprintf("[%02d/%02d] Aggregating %s", i, nrow(sample_meta), donor))
  barcodes <- read_barcodes(sample_meta$barcode_file[[i]])
  mat <- read_mtx_gz(sample_meta$matrix_file[[i]])
  if (nrow(mat) != nrow(features)) stop("Feature count does not match matrix rows for ", donor)

  col_idx <- match(donor_cells$raw_barcode, barcodes)
  valid <- !is.na(col_idx)
  donor_cells <- donor_cells[valid, , drop = FALSE]
  col_idx <- col_idx[valid]

  for (state in c("MA_low", "MA_high")) {
    state_idx <- col_idx[donor_cells$ma_state_donor == state]
    if (length(state_idx) < 10) next
    sample_id <- paste(donor, state, sep = "_")
    total_umis <- sum(Matrix::colSums(mat[, state_idx, drop = FALSE]))
    rps_counts <- sum(mat[rps_row, state_idx, drop = TRUE])
    state_rows[[sample_id]] <- tibble::tibble(
      sample_id = sample_id,
      orig.ident = donor,
      gsm = as.character(sample_meta$gsm[[i]]),
      disease_group = as.character(sample_meta$disease_group[[i]]),
      severity_score = sample_meta$severity_score[[i]],
      state = state,
      n_cells = length(state_idx),
      total_umis = total_umis,
      RPS27L_counts = rps_counts,
      RPS27L_detected_cells = sum(mat[rps_row, state_idx, drop = TRUE] > 0),
      RPS27L_positive_pct = 100 * RPS27L_detected_cells / n_cells,
      RPS27L_cpm = 1e6 * RPS27L_counts / total_umis,
      RPS27L_log1p_cpm = log1p(RPS27L_cpm)
    )
    pb_counts[[sample_id]] <- Matrix::rowSums(mat[, state_idx, drop = FALSE])
    pb_meta_rows[[sample_id]] <- tibble::tibble(
      sample_id = sample_id,
      orig.ident = donor,
      gsm = as.character(sample_meta$gsm[[i]]),
      disease_group = as.character(sample_meta$disease_group[[i]]),
      severity_score = sample_meta$severity_score[[i]],
      state = state,
      n_cells = length(state_idx)
    )
  }
  rm(mat)
  gc(verbose = FALSE)
}

pb_meta <- dplyr::bind_rows(pb_meta_rows)
rps_state <- dplyr::bind_rows(state_rows) |>
  dplyr::mutate(
    disease_group = factor(disease_group, levels = group_levels),
    state = factor(state, levels = c("MA_low", "MA_high"), labels = c("MA-low", "MA-high"))
  )

counts_mat <- do.call(cbind, pb_counts)
colnames(counts_mat) <- names(pb_counts)
valid_gene <- !is.na(features$gene) & features$gene != ""
counts_by_gene <- rowsum(as.matrix(counts_mat[valid_gene, , drop = FALSE]), group = features$gene[valid_gene], reorder = FALSE)

paired_donors <- pb_meta |>
  dplyr::count(orig.ident, state) |>
  tidyr::pivot_wider(names_from = state, values_from = n, values_fill = 0) |>
  dplyr::filter(MA_low > 0, MA_high > 0) |>
  dplyr::pull(orig.ident)
de_meta <- pb_meta |>
  dplyr::filter(orig.ident %in% paired_donors, n_cells >= 10) |>
  dplyr::arrange(orig.ident, state)
de_counts <- counts_by_gene[, de_meta$sample_id, drop = FALSE]

dge <- edgeR::DGEList(counts = round(de_counts), samples = de_meta)
keep <- edgeR::filterByExpr(dge, group = de_meta$state)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge)
de_meta$orig.ident <- factor(de_meta$orig.ident)
de_meta$state <- relevel(factor(de_meta$state), ref = "MA_low")
design <- model.matrix(~ orig.ident + state, data = de_meta)
dge <- edgeR::estimateDisp(dge, design)
fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
coef_name <- grep("^stateMA_high$", colnames(design), value = TRUE)
qlf <- edgeR::glmQLFTest(fit, coef = coef_name)
de_res <- edgeR::topTags(qlf, n = Inf)$table |>
  tibble::rownames_to_column("gene") |>
  tibble::as_tibble() |>
  dplyr::arrange(FDR)

ranks <- de_res |>
  dplyr::filter(is.finite(logFC), !is.na(gene), gene != "") |>
  dplyr::group_by(gene) |>
  dplyr::summarise(rank = logFC, .groups = "drop") |>
  dplyr::arrange(dplyr::desc(rank))
rank_vec <- ranks$rank
names(rank_vec) <- ranks$gene
rank_vec <- sort(rank_vec, decreasing = TRUE)

msig_hallmark <- get_msig("H") |>
  dplyr::select(gs_name, gene_symbol) |>
  dplyr::distinct()
msig_reactome <- get_msig("C2", "CP:REACTOME") |>
  dplyr::select(gs_name, gene_symbol) |>
  dplyr::distinct()
run_gsea <- function(term2gene, source_name) {
  out <- suppressWarnings(clusterProfiler::GSEA(
    geneList = rank_vec,
    TERM2GENE = term2gene,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = TRUE
  ))
  tibble::as_tibble(out@result) |>
    dplyr::mutate(source = source_name)
}
gsea_all <- dplyr::bind_rows(
  run_gsea(msig_hallmark, "Hallmark"),
  run_gsea(msig_reactome, "Reactome")
) |>
  dplyr::arrange(p.adjust)

fig1e_terms <- tibble::tribble(
  ~ID, ~word_group, ~word_priority_note,
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "MA-high enriched", "TNFA/NFKB signaling",
  "HALLMARK_P53_PATHWAY", "MA-high enriched", "p53 signaling",
  "HALLMARK_INFLAMMATORY_RESPONSE", "MA-high enriched", "Inflammatory response",
  "HALLMARK_APOPTOSIS", "MA-high enriched", "Apoptosis",
  "REACTOME_ATF6_ATF6_ALPHA_ACTIVATES_CHAPERONES", "MA-high enriched", "ATF6-mediated chaperone response",
  "REACTOME_PHASE_I_FUNCTIONALIZATION_OF_COMPOUNDS", "MA-low enriched", "Xenobiotic metabolism",
  "REACTOME_FATTY_ACID_METABOLISM", "MA-low enriched", "Fatty acid metabolism",
  "REACTOME_TRANSLATION", "MA-low enriched", "Translation/ribosome"
)
gsea_selected <- gsea_all |>
  dplyr::inner_join(fig1e_terms, by = "ID") |>
  dplyr::mutate(
    term_label = short_pathway_label(ID, Description),
    neg_log10_fdr = -log10(pmax(p.adjust, 1e-300)),
    group_rank = dplyr::if_else(word_group == "MA-high enriched", 1L, 2L),
    display_rank = dplyr::if_else(word_group == "MA-high enriched", -NES, NES)
  ) |>
  dplyr::arrange(group_rank, display_rank)

paired_rps <- rps_state |>
  dplyr::select(orig.ident, disease_group, state, RPS27L_log1p_cpm) |>
  tidyr::pivot_wider(names_from = state, values_from = RPS27L_log1p_cpm)
paired_test <- suppressWarnings(stats::wilcox.test(paired_rps$`MA-high`, paired_rps$`MA-low`, paired = TRUE, exact = FALSE))
paired_stats <- tibble::tibble(
  metric = "Paired RPS27L pseudobulk in donor-internal MA-low vs MA-high hepatocytes",
  comparison = "MA-high vs MA-low",
  n_donors = sum(is.finite(paired_rps$`MA-high`) & is.finite(paired_rps$`MA-low`)),
  statistic = unname(paired_test$statistic),
  p_value = paired_test$p.value,
  significance = p_to_stars(paired_test$p.value)
)

readr::write_tsv(pb_meta, file.path(paths$tables, "Fig1e_GSE289173_donor_internal_MA_high_vs_low_pseudobulk_sample_metadata.tsv"))
readr::write_tsv(de_res, file.path(paths$tables, "Fig1e_GSE289173_donor_internal_MA_high_vs_low_edgeR_DE.tsv"))
readr::write_tsv(gsea_all, file.path(paths$tables, "Fig1e_GSE289173_donor_internal_MA_high_vs_low_Hallmark_Reactome_GSEA.tsv"))
readr::write_tsv(gsea_selected, file.path(paths$tables, "Fig1e_GSE289173_selected_Hallmark_Reactome_GSEA_terms.tsv"))
readr::write_tsv(rps_state, file.path(paths$tables, "Fig1g_GSE289173_donor_internal_MA_high_low_RPS27L_source_data.tsv"))
readr::write_tsv(paired_stats, file.path(paths$tables, "Fig1g_GSE289173_donor_internal_MA_high_low_RPS27L_paired_Wilcoxon.tsv"))
saveRDS(list(counts = counts_by_gene, meta = pb_meta), file.path(paths$processed, "GSE289173_donor_internal_MA_high_low_pseudobulk_counts_edgeR_input.rds"))

message("GSE289173 MA-high/MA-low edgeR and GSEA outputs written to results/tables.")
