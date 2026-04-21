#!/usr/bin/env Rscript

Sys.setenv(XDG_CACHE_HOME = normalizePath("cache", mustWork = FALSE))
.libPaths(c(normalizePath("r_libs"), .libPaths()))
suppressPackageStartupMessages({
  library(limma)
  library(msigdbr)
})

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

expr_tbl <- read.delim("results/gene_level_matrix.tsv", sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
sample_meta <- read.delim("tables/sample_metadata.tsv", sep = "\t", stringsAsFactors = FALSE)

expr_mat <- as.matrix(expr_tbl[, grep("^GSM", names(expr_tbl)), drop = FALSE])
storage.mode(expr_mat) <- "double"
rownames(expr_mat) <- expr_tbl$GB_ACC
sample_meta <- sample_meta[match(colnames(expr_mat), sample_meta$raw_file_gsm), ]

sample_meta$scan_date_c <- ifelse(sample_meta$scan_date == sort(unique(sample_meta$scan_date))[1], -0.5, 0.5)
sample_meta$scanner_c <- ifelse(sample_meta$scanner == sort(unique(sample_meta$scanner))[1], -0.5, 0.5)

display_label <- ifelse(is.na(expr_tbl$GENE_SYMBOL) | expr_tbl$GENE_SYMBOL == "", expr_tbl$GB_ACC, expr_tbl$GENE_SYMBOL)
names(display_label) <- expr_tbl$GB_ACC

plot_small_matrix <- function(mat, file, main, digits = 2) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  cols <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
  pdf(file, width = 7, height = 6)
  par(mar = c(6, 6, 4, 1))
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
  axis(1, at = seq_len(nc), labels = colnames(mat), las = 2, cex.axis = 0.8)
  axis(2, at = seq_len(nr), labels = rev(rownames(mat)), las = 2, cex.axis = 0.8)
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      text(j, nr - i + 1, labels = format(round(mat[i, j], digits), nsmall = digits), cex = 0.8)
    }
  }
  box()
  dev.off()
}

make_pairwise_overlap <- function(deg_sets) {
  model_names <- names(deg_sets)
  do.call(
    rbind,
    lapply(model_names, function(m1) {
      do.call(
        rbind,
        lapply(model_names, function(m2) {
          inter <- length(intersect(deg_sets[[m1]], deg_sets[[m2]]))
          union_n <- length(union(deg_sets[[m1]], deg_sets[[m2]]))
          data.frame(
            model_1 = m1,
            model_2 = m2,
            overlap_n = inter,
            union_n = union_n,
            jaccard = if (union_n == 0) NA_real_ else inter / union_n,
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )
}

fit_deg_model <- function(design, label) {
  fit <- lmFit(expr_mat, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
  tt$GB_ACC <- rownames(tt)
  tt <- cbind(tt, expr_tbl[match(tt$GB_ACC, expr_tbl$GB_ACC), c("ID", "UG_CLUSTER", "DESCRIPTION", "GENE_SYMBOL"), drop = FALSE])
  tt$display_label <- ifelse(is.na(tt$GENE_SYMBOL) | tt$GENE_SYMBOL == "", tt$GB_ACC, tt$GENE_SYMBOL)
  tt$deg_class <- ifelse(
    tt$adj.P.Val < 0.05 & tt$logFC > 0,
    "Up",
    ifelse(tt$adj.P.Val < 0.05 & tt$logFC < 0, "Down", "Not_significant")
  )
  list(label = label, design = design, fit = fit, table = tt)
}

design_base <- cbind(effect = rep(1, ncol(expr_mat)))
design_scan <- cbind(effect = rep(1, ncol(expr_mat)), scan_date = sample_meta$scan_date_c)
design_scanner <- cbind(effect = rep(1, ncol(expr_mat)), scanner = sample_meta$scanner_c)
design_both <- cbind(effect = rep(1, ncol(expr_mat)), scan_date = sample_meta$scan_date_c, scanner = sample_meta$scanner_c)

models <- list(
  base = fit_deg_model(design_base, "base"),
  scan_adjusted = fit_deg_model(design_scan, "scan_adjusted"),
  scanner_adjusted = fit_deg_model(design_scanner, "scanner_adjusted"),
  both_adjusted = fit_deg_model(design_both, "both_adjusted")
)

model_names <- names(models)
model_logfc <- sapply(models, function(x) x$table$logFC[match(expr_tbl$GB_ACC, x$table$GB_ACC)])
model_t <- sapply(models, function(x) x$table$t[match(expr_tbl$GB_ACC, x$table$GB_ACC)])
rownames(model_logfc) <- expr_tbl$GB_ACC
rownames(model_t) <- expr_tbl$GB_ACC

logfc_cor <- cor(model_logfc, use = "pairwise.complete.obs")
t_cor <- cor(model_t, use = "pairwise.complete.obs")
deg_sets <- lapply(models, function(x) x$table$GB_ACC[x$table$deg_class %in% c("Up", "Down")])
deg_overlap <- make_pairwise_overlap(deg_sets)

write.table(logfc_cor, file = "results/model_logFC_correlations.tsv", sep = "\t", quote = FALSE, col.names = NA)
write.table(t_cor, file = "results/model_t_correlations.tsv", sep = "\t", quote = FALSE, col.names = NA)
write.table(deg_overlap, file = "results/model_deg_overlap.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

plot_small_matrix(logfc_cor, "figures/model_logFC_correlation_heatmap.pdf", "Pairwise correlation of logFC across models")
plot_small_matrix(t_cor, "figures/model_t_correlation_heatmap.pdf", "Pairwise correlation of moderated t across models")
plot_small_matrix(
  xtabs(jaccard ~ model_1 + model_2, deg_overlap),
  "figures/model_deg_jaccard_heatmap.pdf",
  "Pairwise Jaccard index of FDR-significant DEGs"
)

write.table(
  data.frame(
    model = names(models),
    n_deg_fdr_0_05 = vapply(models, function(x) sum(x$table$deg_class %in% c("Up", "Down")), integer(1)),
    n_up = vapply(models, function(x) sum(x$table$deg_class == "Up"), integer(1)),
    n_down = vapply(models, function(x) sum(x$table$deg_class == "Down"), integer(1)),
    stringsAsFactors = FALSE
  ),
  file = "results/model_sensitivity.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

for (nm in names(models)) {
  tt <- models[[nm]]$table
  write.table(tt, file = file.path("results", paste0("deg_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

main_model <- models$both_adjusted
main_deg <- main_model$table

write.table(main_deg, file = "results/deg_all.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(subset(main_deg, deg_class == "Up"), file = "results/deg_up.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(subset(main_deg, deg_class == "Down"), file = "results/deg_down.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(
  data.frame(
    metric = c(
      "main_model",
      "genes_tested",
      "DEGs_FDR_lt_0.05_total",
      "DEGs_FDR_lt_0.05_up",
      "DEGs_FDR_lt_0.05_down"
    ),
    value = c(
      "both_adjusted_scan_date_and_scanner",
      nrow(main_deg),
      sum(main_deg$deg_class %in% c("Up", "Down")),
      sum(main_deg$deg_class == "Up"),
      sum(main_deg$deg_class == "Down")
    ),
    stringsAsFactors = FALSE
  ),
  file = "results/deg_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf("figures/volcano_deg.pdf", width = 8, height = 6)
y <- -log10(main_deg$adj.P.Val)
cols <- ifelse(main_deg$deg_class == "Up", "#d7301f", ifelse(main_deg$deg_class == "Down", "#2b8cbe", "#bdbdbd"))
plot(main_deg$logFC, y, pch = 16, cex = 0.6, col = cols, xlab = "log2 fold change (CrM / Placebo)", ylab = "-log10(FDR)", main = "Volcano Plot: Both-adjusted DEG model")
abline(h = -log10(0.05), lty = 2, col = "grey40")
top_idx <- head(order(main_deg$adj.P.Val, -abs(main_deg$logFC)), 15)
with(main_deg[top_idx, ], text(logFC, -log10(adj.P.Val), labels = display_label, pos = 3, cex = 0.7))
dev.off()

complete_idx <- rowSums(!is.finite(expr_mat)) == 0
enrich_mat <- expr_mat[complete_idx, , drop = FALSE]
enrich_symbols <- display_label[rownames(enrich_mat)]
main_t <- main_model$fit$t[rownames(enrich_mat), 1]
names(main_t) <- enrich_symbols

build_geneset_indices <- function(collection, subcollection = NULL) {
  gs <- msigdbr(species = "human", collection = collection, subcollection = subcollection)
  gs <- gs[!is.na(gs$gene_symbol) & gs$gene_symbol != "", c("gs_name", "gene_symbol")]
  idx <- ids2indices(split(gs$gene_symbol, gs$gs_name), identifiers = enrich_symbols)
  idx[lengths(idx) >= 5 & lengths(idx) <= 500]
}

run_camera_collection <- function(design, idx) {
  out <- camera(enrich_mat, idx, design = design, contrast = 1, sort = TRUE)
  out$geneset <- rownames(out)
  out
}

geneset_indices <- list(
  hallmark = build_geneset_indices("H"),
  reactome = build_geneset_indices("C2", "CP:REACTOME"),
  kegg_legacy = build_geneset_indices("C2", "CP:KEGG_LEGACY"),
  go_bp = build_geneset_indices("C5", "GO:BP")
)

enrichment_by_model <- lapply(
  models,
  function(mdl) {
    lapply(geneset_indices, function(idx) run_camera_collection(mdl$design, idx))
  }
)

hallmark <- enrichment_by_model$both_adjusted$hallmark
reactome <- enrichment_by_model$both_adjusted$reactome
kegg <- enrichment_by_model$both_adjusted$kegg_legacy
go_bp <- enrichment_by_model$both_adjusted$go_bp

write.table(hallmark, file = "results/enrichment_hallmark.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(reactome, file = "results/enrichment_reactome.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(kegg, file = "results/enrichment_kegg_legacy.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(go_bp, file = "results/enrichment_go_bp.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

top_summary <- rbind(
  transform(head(hallmark[, c("geneset", "NGenes", "Direction", "PValue", "FDR")], 10), source = "Hallmark"),
  transform(head(reactome[, c("geneset", "NGenes", "Direction", "PValue", "FDR")], 10), source = "Reactome"),
  transform(head(kegg[, c("geneset", "NGenes", "Direction", "PValue", "FDR")], 10), source = "KEGG_legacy"),
  transform(head(go_bp[, c("geneset", "NGenes", "Direction", "PValue", "FDR")], 10), source = "GO_BP")
)
write.table(top_summary, file = "results/enrichment_top_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

selected_pathways <- data.frame(
  source = c(
    "hallmark", "hallmark", "hallmark",
    "reactome", "reactome", "reactome",
    "kegg_legacy", "kegg_legacy",
    "go_bp", "go_bp"
  ),
  geneset = c(
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MYOGENESIS",
    "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
    "REACTOME_RIBOSOME_ASSOCIATED_QUALITY_CONTROL",
    "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
    "KEGG_RIBOSOME",
    "KEGG_OXIDATIVE_PHOSPHORYLATION",
    "GOBP_CYTOPLASMIC_TRANSLATION",
    "GOBP_OXIDATIVE_PHOSPHORYLATION"
  ),
  stringsAsFactors = FALSE
)

pathway_stability <- do.call(
  rbind,
  lapply(model_names, function(model_name) {
    do.call(
      rbind,
      lapply(seq_len(nrow(selected_pathways)), function(i) {
        source_name <- selected_pathways$source[i]
        gs_name <- selected_pathways$geneset[i]
        tab <- enrichment_by_model[[model_name]][[source_name]]
        hit <- tab[match(gs_name, tab$geneset), c("NGenes", "Direction", "PValue", "FDR"), drop = FALSE]
        data.frame(
          model = model_name,
          source = source_name,
          geneset = gs_name,
          NGenes = hit$NGenes,
          Direction = hit$Direction,
          PValue = hit$PValue,
          FDR = hit$FDR,
          stringsAsFactors = FALSE
        )
      })
    )
  })
)
write.table(pathway_stability, file = "results/pathway_stability.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

pdf("figures/pathway_stability_across_models.pdf", width = 10, height = 7)
plot(
  NA,
  xlim = c(0.5, length(model_names) + 0.5),
  ylim = c(0.5, nrow(selected_pathways) + 0.5),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = "",
  main = "Selected pathway stability across models"
)
axis(1, at = seq_along(model_names), labels = model_names, las = 2)
axis(2, at = seq_len(nrow(selected_pathways)), labels = rev(selected_pathways$geneset), las = 2, cex.axis = 0.7)
for (i in seq_along(model_names)) {
  for (j in seq_len(nrow(selected_pathways))) {
    pathway_name <- selected_pathways$geneset[nrow(selected_pathways) - j + 1]
    row_hit <- pathway_stability[pathway_stability$model == model_names[i] & pathway_stability$geneset == pathway_name, ]
    point_col <- ifelse(row_hit$Direction == "Up", "#d7301f", "#2b8cbe")
    point_cex <- max(0.8, min(3.2, -log10(row_hit$FDR + 1e-12) / 3))
    points(i, j, pch = 16, col = point_col, cex = point_cex)
  }
}
legend("topright", legend = c("Up", "Down"), col = c("#d7301f", "#2b8cbe"), pch = 16, bty = "n")
dev.off()

selected_indices <- lapply(seq_len(nrow(selected_pathways)), function(i) geneset_indices[[selected_pathways$source[i]]][[selected_pathways$geneset[i]]])
names(selected_indices) <- selected_pathways$geneset
selected_indices <- selected_indices[lengths(selected_indices) > 0]
fry_validation <- fry(enrich_mat, index = selected_indices, design = design_both, contrast = 1, sort = "mixed")
fry_validation$geneset <- rownames(fry_validation)
write.table(fry_validation, file = "results/pathway_validation_fry.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

pathway_of_interest <- "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION"
reactome_gs <- msigdbr(species = "human", collection = "C2", subcollection = "CP:REACTOME")
pathway_genes <- unique(reactome_gs$gene_symbol[reactome_gs$gs_name == pathway_of_interest])
pathway_idx <- which(enrich_symbols %in% pathway_genes)

pdf("figures/pathway_barcode_translation_initiation.pdf", width = 9, height = 5)
barcodeplot(main_t, index = pathway_idx, labels = c("Down in CrM", "Up in CrM"), main = "Pathway of Interest: Reactome translation initiation")
dev.off()

heatmap_symbols <- intersect(pathway_genes, enrich_symbols)
heatmap_rows <- which(enrich_symbols %in% heatmap_symbols)
heatmap_rows <- heatmap_rows[order(abs(main_t[heatmap_rows]), decreasing = TRUE)]
heatmap_rows <- heatmap_rows[seq_len(min(25, length(heatmap_rows)))]
heatmap_mat <- enrich_mat[heatmap_rows, , drop = FALSE]
rownames(heatmap_mat) <- enrich_symbols[heatmap_rows]
heatmap_scaled <- t(scale(t(heatmap_mat)))
heatmap_scaled[!is.finite(heatmap_scaled)] <- 0
col_labels <- ifelse(is.na(sample_meta$matched_subject_label), sample_meta$raw_file_gsm, sample_meta$matched_subject_label)

pdf("figures/pathway_heatmap_translation_initiation.pdf", width = 9, height = 10)
heatmap(
  heatmap_scaled,
  Colv = NA,
  scale = "none",
  col = colorRampPalette(c("#2166ac", "white", "#b2182b"))(64),
  margins = c(8, 10),
  labCol = col_labels,
  main = "Top translation-initiation genes\nrow-scaled log2 ratios"
)
dev.off()

write.table(
  data.frame(
    pathway_of_interest = pathway_of_interest,
    direction = reactome$Direction[match(pathway_of_interest, reactome$geneset)],
    p_value = reactome$PValue[match(pathway_of_interest, reactome$geneset)],
    fdr = reactome$FDR[match(pathway_of_interest, reactome$geneset)],
    stringsAsFactors = FALSE
  ),
  file = "results/pathway_of_interest.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(capture.output(sessionInfo()), "logs/session_info_finalize.txt")
cat("Final modeling and enrichment complete.\n")
