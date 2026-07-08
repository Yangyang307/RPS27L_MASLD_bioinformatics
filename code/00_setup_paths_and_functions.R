## Shared setup for RPS27L MASLD bioinformatics analyses.

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE)))
  }
  getwd()
}

repo_root <- normalizePath(
  Sys.getenv("RPS27L_REPO_ROOT", unset = file.path(get_script_dir(), "..")),
  mustWork = FALSE
)

paths <- list(
  repo = repo_root,
  code = file.path(repo_root, "code"),
  raw = file.path(repo_root, "data", "raw"),
  processed = file.path(repo_root, "data", "processed"),
  results = file.path(repo_root, "results"),
  tables = file.path(repo_root, "results", "tables"),
  figures = file.path(repo_root, "results", "figures"),
  source_data = file.path(repo_root, "results", "source_data"),
  logs = file.path(repo_root, "results", "logs"),
  gene_sets = file.path(repo_root, "data_availability", "gene_sets")
)

invisible(lapply(paths[c("processed", "results", "tables", "figures", "source_data", "logs")],
                 dir.create, recursive = TRUE, showWarnings = FALSE))

required_packages <- c(
  "Matrix", "data.table", "dplyr", "readr", "stringr", "tidyr",
  "ggplot2", "patchwork"
)
optional_packages <- c("edgeR", "clusterProfiler", "msigdbr", "Seurat", "uwot", "ragg")

load_required_packages <- function(pkgs = required_packages) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "))
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

check_optional_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required package(s) for this step: ", paste(missing, collapse = ", "))
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

read_mtx_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  as(Matrix::readMM(con), "dgCMatrix")
}

read_barcodes <- function(path) {
  data.table::fread(cmd = sprintf("gzip -dc %s", shQuote(path)), header = FALSE)[[1]]
}

download_if_missing <- function(url, dest, allow_download = Sys.getenv("RPS27L_ALLOW_DOWNLOAD", "0") == "1") {
  if (file.exists(dest) && file.info(dest)$size > 0) return(invisible(dest))
  if (!allow_download) {
    stop("Missing input file: ", dest, "\nSet RPS27L_ALLOW_DOWNLOAD=1 to allow scripted download from: ", url)
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  status <- system2("curl", c("-L", "--fail", "--retry", "3", "--retry-delay", "2", "-o", dest, url))
  if (!identical(status, 0L) || !file.exists(dest) || file.info(dest)$size == 0) {
    stop("Download failed: ", url)
  }
  invisible(dest)
}

extract_quoted <- function(line) {
  x <- regmatches(line, gregexpr('"[^"]*"', line))[[1]]
  gsub('^"|"$', "", x)
}

get_geo_row <- function(lines, prefix) {
  hit <- lines[startsWith(lines, prefix)]
  if (length(hit) == 0) stop("Missing GEO series matrix row: ", prefix)
  hit[[1]]
}

format_p <- function(p) {
  if (length(p) == 0 || is.na(p)) return("P = NA")
  if (p < 0.001) return("P < 0.001")
  sprintf("P = %.3f", p)
}

p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

zscore <- function(x) {
  sx <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / sx)
}

sem <- function(x) {
  stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
}

dunn_bh <- function(df, value_col, group_col, metric_label) {
  load_required_packages(c("dplyr", "tidyr"))
  d <- df |>
    dplyr::transmute(value = .data[[value_col]], group = as.character(.data[[group_col]])) |>
    dplyr::filter(is.finite(value), !is.na(group))
  d$group <- factor(d$group, levels = unique(d$group))
  n_total <- nrow(d)
  if (n_total < 3 || nlevels(d$group) < 2) {
    return(tibble::tibble(metric = metric_label, group_1 = character(), group_2 = character()))
  }
  d$rank <- rank(d$value, ties.method = "average")
  ties <- table(d$value)
  tie_adj <- 1 - sum(ties^3 - ties) / (n_total^3 - n_total)
  tie_adj <- ifelse(is.finite(tie_adj) && tie_adj > 0, tie_adj, 1)
  rank_mean <- tapply(d$rank, d$group, mean)
  n_by_group <- table(d$group)
  pairs <- combn(names(rank_mean), 2, simplify = FALSE)
  out <- dplyr::bind_rows(lapply(pairs, function(pair) {
    g1 <- pair[[1]]
    g2 <- pair[[2]]
    se <- sqrt((n_total * (n_total + 1) / 12) * tie_adj *
                 (1 / n_by_group[[g1]] + 1 / n_by_group[[g2]]))
    z <- as.numeric((rank_mean[[g1]] - rank_mean[[g2]]) / se)
    p <- 2 * pnorm(-abs(z))
    tibble::tibble(
      metric = metric_label,
      group_1 = g1,
      group_2 = g2,
      n_1 = as.integer(n_by_group[[g1]]),
      n_2 = as.integer(n_by_group[[g2]]),
      z = z,
      p_value = p
    )
  }))
  out |>
    dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"),
                  significance = p_to_stars(p_adj_BH))
}

kw_spearman <- function(df, value_col, group_col, score_col, metric_label) {
  d <- df |>
    dplyr::filter(is.finite(.data[[value_col]]), !is.na(.data[[group_col]]),
                  is.finite(.data[[score_col]]))
  kw <- suppressWarnings(stats::kruskal.test(as.formula(paste(value_col, "~", group_col)), data = d))
  sp <- suppressWarnings(stats::cor.test(d[[value_col]], d[[score_col]], method = "spearman", exact = FALSE))
  tibble::tibble(
    metric = metric_label,
    n = nrow(d),
    kruskal_statistic = unname(kw$statistic),
    kruskal_p = kw$p.value,
    spearman_rho = unname(sp$estimate),
    spearman_p = sp$p.value
  )
}

save_plot <- function(plot, filename, width = 6, height = 4.5, dpi = 600) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(filename, ".pdf"), plot, width = width, height = height, bg = "white")
  ggplot2::ggsave(paste0(filename, ".png"), plot, width = width, height = height, dpi = dpi, bg = "white")
}

theme_publication <- function(base_size = 8, legend_position = "none") {
  ggplot2::theme_classic(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.45, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.45, color = "black"),
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(color = "grey25", hjust = 0)
    )
}

message("Repository root: ", repo_root)

