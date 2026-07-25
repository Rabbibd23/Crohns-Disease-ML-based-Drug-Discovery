############################################################
# Differential Expression with limma on batch-corrected data
# - Input:  batch_corrected_data.csv
# - Output: deg_all.csv, deg_upregulated.csv, deg_downregulated.csv
# - Plots:  PCA, Volcano, Heatmap (PNG, high resolution)
############################################################

########################
# 1. Install packages
########################
packages <- c("limma", "ggplot2", "pheatmap", "RColorBrewer",
              "ggrepel", "data.table", "matrixStats")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(data.table)
library(matrixStats)

########################
# 2. Load expression data
########################
setwd("E:/My Work/CD & T2D/T2D/Datasets")
# Make sure batch_corrected_data.csv is in your working directory
expr <- read.csv("batch_corrected_data.csv",
                 row.names = 1,
                 check.names = FALSE)

# Convert to numeric matrix
expr <- as.matrix(expr)

cat("Expression matrix dimensions: ",
    nrow(expr), "genes x", ncol(expr), "samples\n")

########################
# 3. Define groups (Control vs Disease)
########################
# Assumes column names contain "control" or "disease"
# If your labels differ, edit the patterns below.
sample_names <- colnames(expr)

group <- ifelse(grepl("control", sample_names, ignore.case = TRUE),
                "Control",
                ifelse(grepl("disease", sample_names, ignore.case = TRUE),
                       "Disease", NA))

# Check any NA (unassigned) samples
if (any(is.na(group))) {
  cat("WARNING: Some samples were not classified into Control/Disease:\n")
  print(sample_names[is.na(group)])
  stop("Please fix group assignment patterns.")
}

group <- factor(group, levels = c("Control", "Disease"))
table(group)

########################
# 4. Quick QC: PCA plot
########################
# log2 data is already appropriate for PCA + limma
# PCA on samples (columns)
pca <- prcomp(t(expr), scale. = TRUE)

pca_df <- data.frame(
  Sample = colnames(expr),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = group
)

# Save PCA plot
png("PCA_before_DEG.png", width = 2000, height = 1600, res = 300)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(level = 0.95, linetype = 2) +
  xlab(paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "%)")) +
  theme_bw(base_size = 14) +
  ggtitle("PCA of Batch-Corrected Expression Data") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )
dev.off()

########################
# 5. Design matrix and limma model
########################
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "Disease_vs_Control")
design

# Fit linear model and empirical Bayes
fit <- lmFit(expr, design)
fit <- eBayes(fit)

########################
# 6. Extract DEG results
########################
# Contrast of interest: Disease vs Control
# coef = "Disease_vs_Control" corresponds to logFC(Disease - Control)
topTable_out <- topTable(fit,
                         coef = "Disease_vs_Control",
                         number = Inf,
                         sort.by = "P")

# Add gene column
deg <- as.data.table(topTable_out, keep.rownames = "Gene")

# Write all DEGs (no threshold)
fwrite(deg, file = "deg_all.csv")

########################
# 7. Apply thresholds for up/down-regulated
########################
# You can change these thresholds as needed
logFC_cutoff <- 1       # |log2FC| >= 1
adjP_cutoff  <- 0.05    # FDR < 0.05

deg[, Regulation := "NotSig"]
deg[logFC >=  logFC_cutoff & adj.P.Val < adjP_cutoff, Regulation := "Up"]
deg[logFC <= -logFC_cutoff & adj.P.Val < adjP_cutoff, Regulation := "Down"]

# Subsets
deg_up   <- deg[Regulation == "Up"]
deg_down <- deg[Regulation == "Down"]

# Save filtered tables with only key columns if you want
cols_to_keep <- c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Regulation")
fwrite(deg[, ..cols_to_keep],      file = "deg_filtered.csv")
fwrite(deg_up[, ..cols_to_keep],   file = "deg_upregulated.csv")
fwrite(deg_down[, ..cols_to_keep], file = "deg_downregulated.csv")

cat("Significant upregulated genes:   ", nrow(deg_up), "\n")
cat("Significant downregulated genes: ", nrow(deg_down), "\n")

########################
# 8. Volcano plot
########################
deg$Status <- "NotSig"
deg$Status[deg$Regulation == "Up"]   <- "Upregulated"
deg$Status[deg$Regulation == "Down"] <- "Downregulated"

# For labelling: top 15 by |logFC| among significant
deg$label <- ""
top_label <- deg[Regulation != "NotSig"][order(-abs(logFC))][1:min(15, .N)]
deg[label %in% top_label$Gene, label := Gene]

png("Volcano_plot.png", width = 2200, height = 1600, res = 300)
ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val), color = Status)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c(
    "Upregulated" = "#D73027",
    "Downregulated" = "#4575B4",
    "NotSig" = "grey70"
  )) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(adjP_cutoff),
             linetype = "dashed", color = "black") +
  theme_bw(base_size = 14) +
  xlab("log2 Fold Change (Disease vs Control)") +
  ylab(expression(-log[10]("adj. p-value"))) +
  ggtitle("Volcano Plot of Differential Expression (limma)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  geom_text_repel(
    data = top_label,
    aes(label = Gene),
    size = 3,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size = 0.3
  )
dev.off()

########################
# 9. Heatmap of top DEGs
########################
# Choose top N most significant DEGs for heatmap
top_n <- 50
top_genes <- deg[Regulation != "NotSig"][order(adj.P.Val)][1:min(top_n, .N)]$Gene

expr_top <- expr[top_genes, ]

# Z-score per gene for nicer heatmap
expr_top_scaled <- t(scale(t(expr_top)))

annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(expr_top_scaled)

heat_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)

png("Heatmap_top_DEGs.png", width = 2200, height = 2600, res = 300)
pheatmap(expr_top_scaled,
         color = heat_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize = 8,
         main = paste0("Top ", min(top_n, length(top_genes)),
                       " Differentially Expressed Genes"))
dev.off()

cat("\nAnalysis complete.\n")
cat("Generated files:\n")
cat("  - deg_all.csv\n")
cat("  - deg_filtered.csv (significant + non-significant with flag)\n")
cat("  - deg_upregulated.csv\n")
cat("  - deg_downregulated.csv\n")
cat("  - PCA_before_DEG.png\n")
cat("  - Volcano_plot.png\n")
cat("  - Heatmap_top_DEGs.png\n")
