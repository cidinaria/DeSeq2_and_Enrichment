#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(data.table)
  library(pheatmap)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggrepel)
  library(RColorBrewer)
  library(ashr)
})
# =========================================================
# GSE156171 - DESeq2 pipeline
# Contrast: Infiltrating vs Other
# Other = Superficial + Solid + Mixed
# =========================================================
counts_file   <- "/home/basal_cell_carcinoma/GSE156171//data/GSE156171_raw-counts.csv.gz"
geo_de_file   <- "/home/basal_cell_carcinoma/GSE156171/data/GSE156171_DE_InfiltVsOther_results.csv.gz"
metadata_file <- "/home/basal_cell_carcinoma/GSE156171/metadata"
outdir <- "/home/basal_cell_carcinoma/GSE156171/results"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# -----------------------------
# Helpers
# -----------------------------
save_plot <- function(filename, plot_obj, width = 8, height = 6) {
  ggsave(
    filename = file.path(outdir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )
}
first_existing <- function(x, choices) {
  hit <- intersect(choices, x)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}
detect_gene_col <- function(df) {
  possible <- c(
    "gene", "Gene", "GENE", "symbol", "SYMBOL", "GeneSymbol",
    "external_gene_name", "hgnc_symbol", "X1", "V1"
  )
  col <- first_existing(colnames(df), possible)
  if (is.na(col)) col <- colnames(df)[1]
  col
}
map_gene_ids <- function(ids, keytype = c("SYMBOL", "ENSEMBL")) {
  keytype <- match.arg(keytype)
  mapped <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(ids),
    keytype = keytype,
    columns = c(keytype, "SYMBOL", "ENTREZID")
  ) %>%
    tibble::as_tibble() %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::distinct(.data[[keytype]], .keep_all = TRUE)
  mapped
}
guess_id_type <- function(x) {
  x <- unique(na.omit(x))
  x <- x[seq_len(min(length(x), 100))]
  if (length(x) == 0) return("SYMBOL")
  ens_like <- grepl("^ENSG[0-9]+(\\.[0-9]+)?$", x)
  if (mean(ens_like) > 0.5) "ENSEMBL" else "SYMBOL"
}
clean_ensembl <- function(x) sub("\\..*$", "", x)
make_group <- function(subtype_vec) {
  case_when(
    subtype_vec == "Infiltrating" ~ "Infiltrating",
    subtype_vec %in% c("Superficial", "Solid", "Mixed") ~ "Other",
    TRUE ~ NA_character_
  )
}
# -----------------------------
# Read metadata first
# -----------------------------
message("Reading metadata...")
meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
stopifnot(all(c("sample_id", "subtype") %in% colnames(meta)))
# -----------------------------
# Read counts
# -----------------------------
message("Reading counts...")
counts_df <- fread(counts_file, data.table = FALSE)
gene_col <- detect_gene_col(counts_df)
message("Detected gene column: ", gene_col)
sample_cols <- intersect(colnames(counts_df), meta$sample_id)
if (length(sample_cols) == 0) {
  stop("No sample columns from metadata were found in the count matrix.")
}
message("Detected sample columns: ", paste(sample_cols, collapse = ", "))
gene_ids <- counts_df[[gene_col]]
count_df <- counts_df[, sample_cols, drop = FALSE]
count_df[] <- lapply(count_df, function(x) as.numeric(as.character(x)))
count_mat <- as.matrix(count_df)
rownames(count_mat) <- make.unique(as.character(gene_ids))
keep_non_na <- !is.na(rownames(count_mat)) & rownames(count_mat) != ""
count_mat <- count_mat[keep_non_na, , drop = FALSE]
message("Count matrix dimensions: ", nrow(count_mat), " genes x ", ncol(count_mat), " samples")
# -----------------------------
# Organize metadata to match counts
# -----------------------------
if (!all(colnames(count_mat) %in% meta$sample_id)) {
  missing_meta <- setdiff(colnames(count_mat), meta$sample_id)
  stop("Samples in count matrix missing in metadata: ", paste(missing_meta, collapse = ", "))
}
meta <- meta[match(colnames(count_mat), meta$sample_id), , drop = FALSE]
rownames(meta) <- meta$sample_id
meta$subtype <- trimws(meta$subtype)
meta$subtype <- recode(
  meta$subtype,
  "Infiltrative" = "Infiltrating"
)
valid_subtypes <- c("Superficial", "Solid", "Mixed", "Infiltrating")
if (!all(meta$subtype %in% valid_subtypes)) {
  bad <- unique(meta$subtype[!meta$subtype %in% valid_subtypes])
  stop("Unexpected subtype(s) in metadata: ", paste(bad, collapse = ", "))
}
meta$group <- make_group(meta$subtype)
if (any(is.na(meta$group))) {
  stop("Some samples could not be assigned to Infiltrating/Other.")
}
meta$group <- factor(meta$group, levels = c("Other", "Infiltrating"))
meta$subtype <- factor(meta$subtype, levels = c("Superficial", "Solid", "Mixed", "Infiltrating"))
write.csv(meta, file.path(outdir, "00_metadata_used.csv"), row.names = FALSE)
message("Group counts:")
print(table(meta$group))
print(table(meta$subtype))
# -----------------------------
# Optional gene-id harmonization
# -----------------------------
id_type <- guess_id_type(rownames(count_mat))
message("Guessed gene ID type: ", id_type)
if (id_type == "ENSEMBL") {
  rownames(count_mat) <- clean_ensembl(rownames(count_mat))
  count_mat <- rowsum(count_mat, group = rownames(count_mat))
}
# -----------------------------
# Filter low-count genes
# -----------------------------
keep <- rowSums(count_mat >= 10) >= 3
count_mat_filt <- count_mat[keep, , drop = FALSE]
message("Genes before filtering: ", nrow(count_mat))
message("Genes after filtering: ", nrow(count_mat_filt))
# -----------------------------
# DESeq2
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat_filt),
  colData = meta,
  design = ~ group
)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
# -----------------------------
# PCA by collapsed group
# -----------------------------
pca_group <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_group, "percentVar"))
p1 <- ggplot(pca_group, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_bw(base_size = 12) +
  ggtitle("PCA - Infiltrating vs Other")
save_plot("01_PCA_group.png", p1, 8, 6)
# PCA by original subtype
pca_subtype <- plotPCA(vsd, intgroup = "subtype", returnData = TRUE)
percentVar2 <- round(100 * attr(pca_subtype, "percentVar"))
p2 <- ggplot(pca_subtype, aes(PC1, PC2, color = subtype, label = name)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  xlab(paste0("PC1: ", percentVar2[1], "%")) +
  ylab(paste0("PC2: ", percentVar2[2], "%")) +
  theme_bw(base_size = 12) +
  ggtitle("PCA - original subtypes")
save_plot("02_PCA_subtype.png", p2, 8, 6)
# -----------------------------
# Sample distance heatmap
# -----------------------------
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
ann_col <- data.frame(
  group = meta$group,
  subtype = meta$subtype
)
rownames(ann_col) <- rownames(meta)
png(file.path(outdir, "03_sample_distance_heatmap.png"), width = 2400, height = 2200, res = 250)
pheatmap(
  sampleDistMatrix,
  annotation_col = ann_col,
  annotation_row = ann_col,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
  main = "Sample-to-sample distances"
)
dev.off()
# -----------------------------
# Differential expression
# -----------------------------
res <- results(dds, contrast = c("group", "Infiltrating", "Other"))
res <- lfcShrink(dds, contrast = c("group", "Infiltrating", "Other"), res = res, type = "ashr")
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)
if (id_type == "ENSEMBL") {
annot <- map_gene_ids(res_df$gene_id, keytype = "ENSEMBL") %>%
  dplyr::rename(gene_id = ENSEMBL)
  res_df <- dplyr::left_join(res_df, annot, by = "gene_id")
  res_df$plot_label <- ifelse(is.na(res_df$SYMBOL), res_df$gene_id, res_df$SYMBOL)
} else {
  annot <- map_gene_ids(res_df$gene_id, keytype = "SYMBOL") %>%
   dplyr::rename(gene_id = SYMBOL, ENTREZID_mapped = ENTREZID)
  rres_df <- dplyr::left_join(res_df, annot, by = "gene_id")
  res_df$SYMBOL <- res_df$gene_id
  res_df$ENTREZID <- res_df$ENTREZID_mapped
  res_df$plot_label <- res_df$gene_id
}
write.csv(res_df, file.path(outdir, "04_DESeq2_results_Infiltrating_vs_Other.csv"), row.names = FALSE)
sig_res <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1)
write.csv(sig_res, file.path(outdir, "05_DESeq2_significant_genes.csv"), row.names = FALSE)
# -----------------------------
# GEO comparison
# -----------------------------
if (file.exists(geo_de_file)) {
  geo_df <- fread(geo_de_file, data.table = FALSE)
  write.csv(geo_df, file.path(outdir, "06_GEO_original_DE_results_copy.csv"), row.names = FALSE)
}
# -----------------------------
# MA plot
# -----------------------------
png(file.path(outdir, "07_MA_plot.png"), width = 2200, height = 1800, res = 300)
plotMA(res, ylim = c(-5, 5), main = "MA plot - Infiltrating vs Other")
dev.off()
# -----------------------------
# Volcano
# -----------------------------
volcano_df <- res_df %>%
  mutate(label_to_show = ifelse(rank(padj, ties.method = "first") <= 20, plot_label, NA))
png(file.path(outdir, "08_volcano_plot.png"), width = 2400, height = 2000, res = 300)
EnhancedVolcano(
  volcano_df,
  lab = volcano_df$plot_label,
  x = "log2FoldChange",
  y = "padj",
  selectLab = na.omit(volcano_df$label_to_show),
  xlab = expression(Log[2]~fold~change),
  ylab = expression(-Log[10]~adjusted~italic(P)),
  title = "Infiltrating vs Other",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 1.8,
  labSize = 3.3,
  colAlpha = 0.85,
  legendPosition = "right",
  drawConnectors = TRUE,
  widthConnectors = 0.4
)
dev.off()
# -----------------------------
# Heatmap top DE genes
# -----------------------------
top_n <- 50
top_genes <- sig_res %>%
  arrange(padj) %>%
  slice_head(n = top_n) %>%
  pull(gene_id)
top_genes <- intersect(top_genes, rownames(vsd))
if (length(top_genes) >= 2) {
  mat <- assay(vsd)[top_genes, , drop = FALSE]
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0
  row_labels <- top_genes
  if ("SYMBOL" %in% colnames(res_df)) {
symbol_map <- res_df %>%
  dplyr::select(gene_id, SYMBOL) %>%
  dplyr::distinct()
    row_labels <- symbol_map$SYMBOL[match(top_genes, symbol_map$gene_id)]
    row_labels[is.na(row_labels) | row_labels == ""] <- top_genes[is.na(row_labels) | row_labels == ""]
  }
  rownames(mat_z) <- make.unique(row_labels)
  png(file.path(outdir, "09_heatmap_top50_DE_genes.png"), width = 2300, height = 2800, res = 250)
  pheatmap(
    mat_z,
    annotation_col = ann_col,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 8,
    main = "Top 50 DE genes (z-score)"
  )
  dev.off()
}
# -----------------------------
# Heatmap of selected article-related genes
# -----------------------------
genes_of_interest <- c("WISP1", "POSTN", "COL1A1", "COL1A2", "FN1", "ITGA5", "ITGB1")
available_goi <- intersect(genes_of_interest, res_df$SYMBOL)
if (length(available_goi) >= 2) {
  gene_map <- res_df %>%
    filter(SYMBOL %in% available_goi) %>%
    distinct(SYMBOL, gene_id)
  mat_goi <- assay(vsd)[gene_map$gene_id, , drop = FALSE]
  rownames(mat_goi) <- gene_map$SYMBOL
  mat_goi_z <- t(scale(t(mat_goi)))
  mat_goi_z[is.na(mat_goi_z)] <- 0
  png(file.path(outdir, "10_heatmap_genes_of_interest.png"), width = 2200, height = 1600, res = 250)
  pheatmap(
    mat_goi_z,
    annotation_col = ann_col,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "Genes of interest (z-score)"
  )
  dev.off()
}
# -----------------------------
# Enrichment helper
# -----------------------------
run_enrichment <- function(sig_table, direction = c("up", "down"), out_prefix = "UP") {
  direction <- match.arg(direction)
  if (direction == "up") {
    genes <- sig_table %>% filter(log2FoldChange > 0)
  } else {
    genes <- sig_table %>% filter(log2FoldChange < 0)
  }
  entrez <- unique(na.omit(genes$ENTREZID))
  bg_entrez <- unique(na.omit(res_df$ENTREZID))
  if (length(entrez) < 5) {
    message("Skipping enrichment for ", out_prefix, ": too few mapped genes.")
    return(NULL)
  }
  ego_bp <- enrichGO(
    gene = entrez,
    universe = bg_entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  write.csv(as.data.frame(ego_bp), file.path(outdir, paste0(out_prefix, "_GO_BP.csv")), row.names = FALSE)
  if (nrow(as.data.frame(ego_bp)) > 0) {
    png(file.path(outdir, paste0(out_prefix, "_GO_BP_dotplot.png")), width = 2400, height = 1800, res = 300)
    print(dotplot(ego_bp, showCategory = 20) + ggtitle(paste0(out_prefix, " GO Biological Process")))
    dev.off()
  }
  ego_mf <- enrichGO(
    gene = entrez,
    universe = bg_entrez,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  write.csv(as.data.frame(ego_mf), file.path(outdir, paste0(out_prefix, "_GO_MF.csv")), row.names = FALSE)
  if (nrow(as.data.frame(ego_mf)) > 0) {
    png(file.path(outdir, paste0(out_prefix, "_GO_MF_dotplot.png")), width = 2400, height = 1800, res = 300)
    print(dotplot(ego_mf, showCategory = 20) + ggtitle(paste0(out_prefix, " GO Molecular Function")))
    dev.off()
  }
  ego_cc <- enrichGO(
    gene = entrez,
    universe = bg_entrez,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  write.csv(as.data.frame(ego_cc), file.path(outdir, paste0(out_prefix, "_GO_CC.csv")), row.names = FALSE)
  if (nrow(as.data.frame(ego_cc)) > 0) {
    png(file.path(outdir, paste0(out_prefix, "_GO_CC_dotplot.png")), width = 2400, height = 1800, res = 300)
    print(dotplot(ego_cc, showCategory = 20) + ggtitle(paste0(out_prefix, " GO Cellular Component")))
    dev.off()
  }
  ekegg <- enrichKEGG(
    gene = entrez,
    universe = bg_entrez,
    organism = "hsa",
    pvalueCutoff = 0.05
  )
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    write.csv(as.data.frame(ekegg), file.path(outdir, paste0(out_prefix, "_KEGG.csv")), row.names = FALSE)
    png(file.path(outdir, paste0(out_prefix, "_KEGG_dotplot.png")), width = 2400, height = 1800, res = 300)
    print(dotplot(ekegg, showCategory = 20) + ggtitle(paste0(out_prefix, " KEGG pathways")))
    dev.off()
  }
}
run_enrichment(sig_res, "up",   "11_UP")
run_enrichment(sig_res, "down", "12_DOWN")
# -----------------------------
# Optional subtype-specific model
# -----------------------------
dds_multi <- DESeqDataSetFromMatrix(
  countData = round(count_mat_filt),
  colData = meta,
  design = ~ subtype
)
dds_multi <- DESeq(dds_multi)
res_infil_vs_superficial <- results(dds_multi, contrast = c("subtype", "Infiltrating", "Superficial")) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")
write.csv(
  res_infil_vs_superficial,
  file.path(outdir, "13_DESeq2_Infiltrating_vs_Superficial.csv"),
  row.names = FALSE
)
# -----------------------------
# Summary
# -----------------------------
summary_lines <- c(
  "GSE156171 DESeq2 summary",
  paste("Samples:", ncol(count_mat)),
  paste("Genes after filtering:", nrow(count_mat_filt)),
  paste("Infiltrating samples:", sum(meta$group == 'Infiltrating')),
  paste("Other samples:", sum(meta$group == 'Other')),
  paste("Significant genes (padj<0.05 & |log2FC|>=1):", nrow(sig_res))
)
writeLines(summary_lines, file.path(outdir, "14_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "15_sessionInfo.txt"))
message("Done. Results written to: ", outdir)
