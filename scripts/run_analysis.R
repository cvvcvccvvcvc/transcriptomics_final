#!/usr/bin/env Rscript

.libPaths(c(normalizePath("r_libs"), .libPaths()))
suppressPackageStartupMessages(library(limma))

dir.create("tables", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

series_soft <- "geo/GSE7877_family.soft.gz"
platform_soft <- "geo/GPL5216_family.soft.gz"
raw_dir <- "GSE7877_RAW"

control_patterns <- c(
  "^Blank$",
  "^No Data$",
  "^Unknown$",
  "^Human Cot 1$",
  "^Mouse Cot 1$",
  "^Poly A$",
  "^Salmon Sperm$",
  "^3XSSC",
  "^Spot Report Product",
  "^B-Human Actin$"
)

is_control_name <- function(x) {
  x <- ifelse(is.na(x), "", trimws(x))
  Reduce(`|`, lapply(control_patterns, grepl, x = x, perl = TRUE))
}

trim_field <- function(x) {
  y <- trimws(x)
  y[y %in% c("", "NA", "NULL")] <- NA_character_
  y
}

read_platform_table <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  start <- match("!platform_table_begin", lines)
  end <- match("!platform_table_end", lines)
  if (is.na(start) || is.na(end) || end <= start + 1) {
    stop("Could not locate platform table in ", path)
  }
  table_lines <- lines[(start + 1):(end - 1)]
  con <- textConnection(table_lines)
  on.exit(close(con), add = TRUE)
  tab <- read.delim(con, sep = "\t", quote = "", fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  names(tab) <- gsub(" ", "_", names(tab), fixed = TRUE)
  tab$GB_ACC <- trim_field(tab$GB_ACC)
  tab$GENE_SYMBOL <- trim_field(tab$GENE_SYMBOL)
  tab
}

parse_sample_block <- function(block_lines) {
  get_value <- function(prefix) {
    hit <- grep(prefix, block_lines, value = TRUE)
    if (!length(hit)) {
      return(NA_character_)
    }
    sub(prefix, "", hit[1], perl = TRUE)
  }

  get_processed_table <- function() {
    table_start <- match("!sample_table_begin", block_lines)
    table_end <- match("!sample_table_end", block_lines)
    if (is.na(table_start) || is.na(table_end) || table_end <= table_start + 1) {
      return(NULL)
    }
    con <- textConnection(block_lines[(table_start + 1):(table_end - 1)])
    on.exit(close(con), add = TRUE)
    read.delim(con, sep = "\t", quote = "", fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }

  gsm <- sub("^\\^SAMPLE = ", "", block_lines[1])
  processed <- get_processed_table()
  list(
    sample = data.frame(
      gsm = gsm,
      title = get_value("^!Sample_title = "),
      source_ch1 = get_value("^!Sample_source_name_ch1 = "),
      source_ch2 = get_value("^!Sample_source_name_ch2 = "),
      characteristics_ch1 = get_value("^!Sample_characteristics_ch1 = "),
      characteristics_ch2 = get_value("^!Sample_characteristics_ch2 = "),
      label_ch1 = get_value("^!Sample_label_ch1 = "),
      label_ch2 = get_value("^!Sample_label_ch2 = "),
      description = get_value("^!Sample_description = "),
      supplementary_file = get_value("^!Sample_supplementary_file = "),
      stringsAsFactors = FALSE
    ),
    processed = processed
  )
}

read_series_samples <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  sample_starts <- grep("^\\^SAMPLE = ", lines)
  sample_ends <- c(sample_starts[-1] - 1L, length(lines))

  sample_info <- vector("list", length(sample_starts))
  processed_tables <- vector("list", length(sample_starts))
  names(processed_tables) <- character(length(sample_starts))

  for (i in seq_along(sample_starts)) {
    block <- lines[sample_starts[i]:sample_ends[i]]
    parsed <- parse_sample_block(block)
    sample_info[[i]] <- parsed$sample
    processed_tables[[i]] <- parsed$processed
    names(processed_tables)[i] <- parsed$sample$gsm
  }

  samples <- do.call(rbind, sample_info)
  list(samples = samples, processed = processed_tables)
}

parse_gpr_header <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header_lines <- character()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line) || startsWith(line, "Block\tColumn\tRow\tName")) {
      break
    }
    header_lines <- c(header_lines, line)
  }
  kv <- grep("=", header_lines, value = TRUE)
  keys <- gsub('^"|"$', "", sub("=.*$", "", kv))
  vals <- gsub('^"|"$', "", sub("^[^=]+=", "", kv))
  meta <- stats::setNames(as.list(vals), keys)
  image_files <- strsplit(meta$ImageFiles %||% "", "\t", fixed = TRUE)[[1]]
  cy5_image <- if (length(image_files) >= 1L) image_files[1] else NA_character_
  slide_id <- sub("^.*\\\\", "", cy5_image)
  slide_id <- sub("_Cy5.tif 0\"?$", "", slide_id)
  list(
    raw_file_gsm = sub("\\.gpr\\.gz$", "", basename(path)),
    slide_id = trim_field(slide_id),
    scan_datetime = trim_field(meta$DateTime %||% NA_character_),
    scanner = trim_field(meta$Scanner %||% NA_character_)
  )
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

collapse_rows_mean <- function(mat, group) {
  split_idx <- split(seq_along(group), group)
  collapsed <- vapply(
    split_idx,
    function(idx) {
      vals <- mat[idx, , drop = FALSE]
      out <- colMeans(vals, na.rm = TRUE)
      out[apply(is.na(vals), 2, all)] <- NA_real_
      out
    },
    numeric(ncol(mat))
  )
  out <- t(collapsed)
  rownames(out) <- names(split_idx)
  colnames(out) <- colnames(mat)
  out
}

plot_density_panel <- function(mat, file, main) {
  finite_rows <- rowSums(is.finite(mat)) > 0
  pdf(file, width = 8, height = 6)
  plotDensities(mat[finite_rows, , drop = FALSE], main = main, legend = "topright")
  dev.off()
}

plot_matrix_heatmap <- function(mat, file, main, row_labels = rownames(mat), col_labels = colnames(mat), digits = 2) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  cols <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
  pdf(file, width = 10, height = 9)
  par(mar = c(10, 10, 4, 2))
  image(
    x = seq_len(nc),
    y = seq_len(nr),
    z = t(mat[nr:1, , drop = FALSE]),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main,
    col = cols,
    zlim = c(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE))
  )
  axis(1, at = seq_len(nc), labels = col_labels, las = 2, cex.axis = 0.7)
  axis(2, at = seq_len(nr), labels = rev(row_labels), las = 2, cex.axis = 0.7)
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      text(j, nr - i + 1, labels = format(round(mat[i, j], digits), nsmall = digits), cex = 0.55)
    }
  }
  box()
  dev.off()
}

plot_pca <- function(mat, sample_meta, file, main) {
  keep <- rowSums(is.finite(mat)) >= ncol(mat) / 2
  x <- t(mat[keep, , drop = FALSE])
  x[!is.finite(x)] <- 0
  pca <- prcomp(x, center = TRUE, scale. = FALSE)
  pct <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

  scan_levels <- unique(sample_meta$scan_date)
  cols <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")[seq_along(scan_levels)], scan_levels)
  scan_col <- cols[sample_meta$scan_date]
  scanner_levels <- unique(sample_meta$scanner)
  pchs <- setNames(c(16, 17, 15, 18)[seq_along(scanner_levels)], scanner_levels)
  scan_pch <- pchs[sample_meta$scanner]

  pdf(file, width = 8, height = 6)
  plot(
    pca$x[, 1], pca$x[, 2],
    col = scan_col,
    pch = scan_pch,
    xlab = sprintf("PC1 (%.1f%%)", pct[1]),
    ylab = sprintf("PC2 (%.1f%%)", pct[2]),
    main = main
  )
  text(pca$x[, 1], pca$x[, 2], labels = sample_meta$subject_label, pos = 3, cex = 0.8)
  legend("topright", legend = names(cols), col = cols, pch = 16, title = "Scan date", bty = "n")
  legend("bottomright", legend = names(pchs), pch = pchs, title = "Scanner", bty = "n")
  dev.off()

  data.frame(
    gsm = sample_meta$gsm,
    subject_label = sample_meta$subject_label,
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
}

volcano_plot <- function(tab, file) {
  y <- -log10(tab$adj.P.Val)
  cols <- ifelse(tab$deg_class == "Up", "#d7301f", ifelse(tab$deg_class == "Down", "#2b8cbe", "#bdbdbd"))
  pdf(file, width = 8, height = 6)
  plot(tab$logFC, y, pch = 16, cex = 0.6, col = cols, xlab = "log2 fold change (CrM / Placebo)", ylab = "-log10(FDR)", main = "DEG Volcano Plot")
  abline(h = -log10(0.05), lty = 2, col = "grey40")
  top_idx <- head(order(tab$adj.P.Val, -abs(tab$logFC)), 15)
  with(tab[top_idx, ], text(logFC, -log10(adj.P.Val), labels = display_label, pos = 3, cex = 0.7))
  dev.off()
}

session_info_path <- "logs/session_info.txt"
writeLines(capture.output(sessionInfo()), session_info_path)

platform <- read_platform_table(platform_soft)
series_data <- read_series_samples(series_soft)

sample_meta <- series_data$samples
sample_meta$soft_gsm <- trim_field(sample_meta$gsm)
sample_meta$soft_subject_id <- as.integer(sub(".*subject", "", sample_meta$title))
sample_meta$soft_subject_label <- trimws(sub(" - CrM vs. PL", "", sample_meta$description, fixed = TRUE))
sample_meta <- sample_meta[order(sample_meta$soft_subject_id), ]
sample_meta$raw_file_gsm <- sample_meta$soft_gsm

gpr_files <- file.path(raw_dir, paste0(sample_meta$raw_file_gsm, ".gpr.gz"))
header_meta <- do.call(rbind, lapply(gpr_files, function(x) as.data.frame(parse_gpr_header(x), stringsAsFactors = FALSE)))
header_meta$scan_date <- sub(" .*", "", header_meta$scan_datetime)

sample_meta <- merge(sample_meta, header_meta, by = "raw_file_gsm", all.x = TRUE, sort = FALSE)
sample_meta <- sample_meta[order(sample_meta$soft_subject_id), ]
sample_meta$contrast_label <- "CrM_vs_Placebo"
sample_meta$channel_635 <- "Placebo"
sample_meta$channel_532 <- "Creatine_supplement"
sample_meta$log_ratio_orientation <- "log2(Creatine_supplement / Placebo)"
sample_meta$matched_gsm <- NA_character_
sample_meta$matched_subject_id <- NA_integer_
sample_meta$matched_subject_label <- NA_character_
sample_meta$matched_description <- NA_character_
sample_meta$cor_with_matched_geo_processed <- NA_real_

write.table(platform, file = "tables/platform_annotation.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

RG <- read.maimages(
  files = gpr_files,
  source = "genepix",
  columns = list(R = "F635 Median", G = "F532 Median", Rb = "B635 Median", Gb = "B532 Median"),
  other.columns = c("Flags")
)

for (slot_name in c("R", "G", "Rb", "Gb")) {
  colnames(RG[[slot_name]]) <- sample_meta$raw_file_gsm
}
colnames(RG$other$Flags) <- sample_meta$raw_file_gsm

platform_key <- platform[!is.na(platform$GB_ACC) & nzchar(platform$GB_ACC), ]
platform_key <- platform_key[!duplicated(platform_key$GB_ACC), ]
annot <- platform_key[match(trim_field(RG$genes$Name), platform_key$GB_ACC), ]
names(annot) <- paste0("platform_", names(annot))
RG$genes <- cbind(RG$genes, annot)
RG$targets <- sample_meta

control_spots <- is_control_name(RG$genes$Name) | is.na(RG$genes$platform_GB_ACC)
bad_flags <- RG$other$Flags < 0

RGb <- backgroundCorrect(RG, method = "normexp", offset = 50)
RGb$weights <- matrix(1, nrow = nrow(RGb$R), ncol = ncol(RGb$R))
RGb$weights[bad_flags] <- 0

MA_raw <- normalizeWithinArrays(RGb, method = "none")
MA_norm <- normalizeWithinArrays(RGb, method = "printtiploess", controlspots = control_spots)

# Reorient from native 635/532 (placebo/CrM) to CrM/placebo for interpretation.
MA_raw$M <- -MA_raw$M
MA_norm$M <- -MA_norm$M

MA_raw$M[bad_flags] <- NA_real_
MA_raw$A[bad_flags] <- NA_real_
MA_norm$M[bad_flags] <- NA_real_
MA_norm$A[bad_flags] <- NA_real_

noncontrol_idx <- which(!control_spots)
raw_m <- MA_raw$M[noncontrol_idx, , drop = FALSE]
norm_m <- MA_norm$M[noncontrol_idx, , drop = FALSE]
norm_a <- MA_norm$A[noncontrol_idx, , drop = FALSE]
gb_acc <- RG$genes$platform_GB_ACC[noncontrol_idx]

raw_gene <- collapse_rows_mean(raw_m, gb_acc)
norm_gene <- collapse_rows_mean(norm_m, gb_acc)
norm_gene_a <- collapse_rows_mean(norm_a, gb_acc)

gene_annot <- platform_key[match(rownames(norm_gene), platform_key$GB_ACC), ]
rownames(gene_annot) <- gene_annot$GB_ACC

median_a <- apply(norm_gene_a, 1, median, na.rm = TRUE)
nonmissing_n <- rowSums(is.finite(norm_gene))
keep_gene <- nonmissing_n >= 8 & median_a > 8

analysis_matrix <- norm_gene[keep_gene, , drop = FALSE]
analysis_a <- norm_gene_a[keep_gene, , drop = FALSE]
analysis_annot <- gene_annot[rownames(analysis_matrix), , drop = FALSE]

write.table(
  cbind(GB_ACC = rownames(analysis_matrix), analysis_annot, analysis_matrix, median_A = median_a[keep_gene], nonmissing_n = nonmissing_n[keep_gene]),
  file = "results/gene_level_matrix.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_density_panel(raw_gene[keep_gene, , drop = FALSE], "figures/density_raw_gene_level.pdf", "Raw log2 ratios before within-array normalization")
plot_density_panel(analysis_matrix, "figures/density_normalized_gene_level.pdf", "Gene-level log2 ratios after print-tip loess normalization")

pdf("figures/ma_plot_subject1_before_after.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plotMA(MA_raw, array = 1, main = "Subject 1 before normalization")
plotMA(MA_norm, array = 1, main = "Subject 1 after print-tip loess")
dev.off()

design <- matrix(1, nrow = ncol(analysis_matrix), ncol = 1)
colnames(design) <- "CrM_vs_Placebo"
rownames(design) <- sample_meta$raw_file_gsm

fit <- lmFit(analysis_matrix, design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
deg_table <- topTable(fit, coef = "CrM_vs_Placebo", number = Inf, sort.by = "P")
deg_table$GB_ACC <- rownames(deg_table)
deg_table <- cbind(deg_table, analysis_annot[deg_table$GB_ACC, c("ID", "UG_CLUSTER", "DESCRIPTION", "GENE_SYMBOL"), drop = FALSE])
deg_table$display_label <- ifelse(is.na(deg_table$GENE_SYMBOL) | deg_table$GENE_SYMBOL == "", deg_table$GB_ACC, deg_table$GENE_SYMBOL)
deg_table$deg_class <- ifelse(
  deg_table$adj.P.Val < 0.05 & deg_table$logFC > 0,
  "Up",
  ifelse(deg_table$adj.P.Val < 0.05 & deg_table$logFC < 0, "Down", "Not_significant")
)

write.table(deg_table, file = "results/deg_all.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(subset(deg_table, deg_class == "Up"), file = "results/deg_up.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(subset(deg_table, deg_class == "Down"), file = "results/deg_down.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

volcano_plot(deg_table, "figures/volcano_deg.pdf")

deg_summary <- data.frame(
  metric = c(
    "genes_after_control_and_mapping_filter",
    "genes_after_intensity_and_missingness_filter",
    "DEGs_FDR_lt_0.05_total",
    "DEGs_FDR_lt_0.05_up",
    "DEGs_FDR_lt_0.05_down"
  ),
  value = c(
    nrow(norm_gene),
    nrow(analysis_matrix),
    sum(deg_table$deg_class %in% c("Up", "Down")),
    sum(deg_table$deg_class == "Up"),
    sum(deg_table$deg_class == "Down")
  ),
  stringsAsFactors = FALSE
)
write.table(deg_summary, file = "results/deg_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

processed_list <- series_data$processed
processed_matrix <- lapply(processed_list[sample_meta$soft_gsm], function(tab) {
  vals <- tab$VALUE
  vals[vals == ""] <- NA
  out <- as.numeric(vals)
  names(out) <- tab$ID_REF
  out
})
processed_ids <- Reduce(intersect, lapply(processed_matrix, names))
processed_mat <- do.call(cbind, lapply(processed_matrix, function(x) x[processed_ids]))
colnames(processed_mat) <- sample_meta$soft_gsm
rownames(processed_mat) <- processed_ids

norm_by_platform_id <- analysis_matrix
platform_ids <- analysis_annot$ID
valid_platform <- !is.na(platform_ids) & nzchar(platform_ids)
rownames(norm_by_platform_id) <- platform_ids
norm_by_platform_id <- norm_by_platform_id[valid_platform, , drop = FALSE]

shared_ids <- intersect(rownames(norm_by_platform_id), rownames(processed_mat))
sample_cor_matrix <- cor(norm_by_platform_id[shared_ids, , drop = FALSE], processed_mat[shared_ids, , drop = FALSE], use = "pairwise.complete.obs")
best_match_idx <- max.col(sample_cor_matrix, ties.method = "first")
matched_gsms <- colnames(sample_cor_matrix)[best_match_idx]
if (anyDuplicated(matched_gsms)) {
  stop("Sample matching from raw files to GEO processed tables was not one-to-one.")
}
sample_meta$matched_gsm <- matched_gsms[match(sample_meta$raw_file_gsm, rownames(sample_cor_matrix))]
sample_meta$cor_with_matched_geo_processed <- sample_cor_matrix[cbind(match(sample_meta$raw_file_gsm, rownames(sample_cor_matrix)), best_match_idx)]
matched_rows <- match(sample_meta$matched_gsm, series_data$samples$gsm)
sample_meta$matched_subject_id <- as.integer(sub(".*subject", "", series_data$samples$title[matched_rows]))
sample_meta$matched_subject_label <- trimws(sub(" - CrM vs. PL", "", series_data$samples$description[matched_rows], fixed = TRUE))
sample_meta$matched_description <- series_data$samples$description[matched_rows]
sample_meta$file_label_matches_processed_profile <- sample_meta$raw_file_gsm == sample_meta$matched_gsm

orientation_check <- data.frame(
  raw_file_gsm = sample_meta$raw_file_gsm,
  matched_gsm = sample_meta$matched_gsm,
  correlation_current_orientation = vapply(
    seq_len(nrow(sample_meta)),
    function(i) cor(norm_by_platform_id[shared_ids, sample_meta$raw_file_gsm[i]], processed_mat[shared_ids, sample_meta$matched_gsm[i]], use = "pairwise.complete.obs"),
    numeric(1)
  ),
  correlation_sign_flipped = vapply(
    seq_len(nrow(sample_meta)),
    function(i) cor(-norm_by_platform_id[shared_ids, sample_meta$raw_file_gsm[i]], processed_mat[shared_ids, sample_meta$matched_gsm[i]], use = "pairwise.complete.obs"),
    numeric(1)
  ),
  stringsAsFactors = FALSE
)
orientation_check$orientation_consistent <- orientation_check$correlation_current_orientation > 0 & orientation_check$correlation_sign_flipped < 0
orientation_check$correlation_gap <- orientation_check$correlation_current_orientation - orientation_check$correlation_sign_flipped

plot_meta <- sample_meta
plot_meta$subject_label <- ifelse(is.na(plot_meta$matched_subject_label), plot_meta$soft_subject_label, plot_meta$matched_subject_label)
plot_meta$gsm <- plot_meta$raw_file_gsm
raw_pca <- plot_pca(raw_gene[keep_gene, , drop = FALSE], plot_meta, "figures/pca_raw_gene_level.pdf", "PCA before within-array normalization")
norm_pca <- plot_pca(analysis_matrix, plot_meta, "figures/pca_normalized_gene_level.pdf", "PCA after print-tip loess normalization")
write.table(raw_pca, file = "tables/pca_raw_scores.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(norm_pca, file = "tables/pca_normalized_scores.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

heatmap_row_labels <- paste0(sample_meta$raw_file_gsm, "\n", sample_meta$slide_id)
heatmap_col_labels <- paste0(sample_meta$soft_gsm, "\n", sample_meta$soft_subject_label)
plot_matrix_heatmap(
  sample_cor_matrix,
  "figures/raw_to_geo_processed_correlation_heatmap.pdf",
  "Raw-to-GEO processed profile correlation",
  row_labels = heatmap_row_labels,
  col_labels = heatmap_col_labels
)

pdf("figures/orientation_sanity_check.pdf", width = 9, height = 5)
par(mar = c(8, 4, 4, 1))
bar_centers <- barplot(
  t(as.matrix(orientation_check[, c("correlation_current_orientation", "correlation_sign_flipped")])),
  beside = TRUE,
  col = c("#1b9e77", "#d95f02"),
  ylim = c(-1, 1),
  las = 2,
  names.arg = orientation_check$raw_file_gsm,
  ylab = "Correlation with matched GEO processed profile",
  main = "Orientation sanity check"
)
abline(h = 0, lty = 2, col = "grey40")
legend("topleft", legend = c("Current orientation", "Sign-flipped"), fill = c("#1b9e77", "#d95f02"), bty = "n")
dev.off()

write.table(sample_cor_matrix, file = "results/raw_to_geo_processed_correlation_matrix.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(
  sample_meta[, c(
    "raw_file_gsm",
    "matched_gsm",
    "soft_subject_label",
    "matched_subject_label",
    "cor_with_matched_geo_processed",
    "file_label_matches_processed_profile",
    "slide_id",
    "scan_date",
    "scanner"
  )],
  file = "results/raw_to_geo_sample_mapping.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  data.frame(raw_file_gsm = sample_meta$raw_file_gsm, matched_gsm = sample_meta$matched_gsm, cor_with_matched_geo_processed = sample_meta$cor_with_matched_geo_processed, stringsAsFactors = FALSE),
  file = "results/correlation_with_geo_processed.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(orientation_check, file = "results/orientation_sanity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(sample_meta, file = "tables/sample_metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("Analysis complete.\n")
