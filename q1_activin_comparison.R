# Question 1: Compare DGE for 0 ng/ml vs 15 ng/ml Activin - Exp1 vs Exp2
# Straightforward academic comparison with clean visualizations

library(DESeq2)
library(readr)
library(tidyverse)
library(readxl)
library(pheatmap)

# ============================================================================
# LOAD AND PREPARE DATA
# ============================================================================

counts_raw <- read.table("salmon.merged.gene_counts.tsv", header = TRUE, row.names = 1)
counts_raw <- counts_raw[, -1]
counts_int <- round(counts_raw)

samples <- read_csv("samples.csv")
parse_sample_names <- function(sample_name) {
  sub("^\\d{8}_R\\d+_", "", sample_name)
}
clean_treatments <- sapply(samples$sample_description, parse_sample_names, USE.NAMES = FALSE)

metadata <- data.frame(
  sample = paste0("S", samples$requests_sample_sample_id),
  treatment = clean_treatments,
  row.names = paste0("S", samples$requests_sample_sample_id)
)

metadata$experiment <- ifelse(grepl("DMSO|SB50", metadata$treatment), "Exp2", "Exp1")
metadata$concentration <- case_when(
  grepl("^0ngmlActivin", metadata$treatment) ~ "0ngml",
  grepl("^5ngmlActivin", metadata$treatment) ~ "5ngml",
  grepl("^10ngmlActivin", metadata$treatment) ~ "10ngml",
  grepl("^15ngmlActivin", metadata$treatment) ~ "15ngml",
  grepl("DMSO", metadata$treatment) & grepl("^0ngml", metadata$treatment) ~ "0ngml_DMSO",
  grepl("DMSO", metadata$treatment) & grepl("^15ngml", metadata$treatment) ~ "15ngml_DMSO",
  TRUE ~ NA_character_
)
metadata$time_min <- as.numeric(str_extract(metadata$treatment, "\\d+(?=min$)"))

# Filter counts
common_samples <- intersect(colnames(counts_int), rownames(metadata))
counts_int <- counts_int[, common_samples]
metadata <- metadata[common_samples, ]

# Keep genes with at least 10 counts in at least 3 samples
keep <- rowSums(counts_int >= 10) >= 3
counts_filtered <- counts_int[keep, ]

cat("Total samples:", ncol(counts_filtered), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")

# ============================================================================
# ANALYSIS: 0 vs 15 ng/ml Activin - Experiment 1
# ============================================================================

cat("\n========================================\n")
cat("EXPERIMENT 1: 0 vs 15 ng/ml Activin\n")
cat("========================================\n")

# Use 240min timepoint
exp1_240min <- metadata[
  metadata$experiment == "Exp1" &
  metadata$concentration %in% c("0ngml", "15ngml") &
  metadata$time_min == 240,
]

cat("Samples:", nrow(exp1_240min), "\n")

exp1_240min$concentration <- factor(exp1_240min$concentration, levels = c("0ngml", "15ngml"))

dds_exp1 <- DESeqDataSetFromMatrix(
    countData = counts_filtered[, rownames(exp1_240min)],
    colData = exp1_240min,
    design = ~ concentration
  )
dds_exp1 <- DESeq(dds_exp1, quiet = TRUE)
res_exp1 <- results(dds_exp1, contrast = c("concentration", "15ngml", "0ngml"))

res_exp1_df <- as.data.frame(res_exp1) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

n_de <- sum(res_exp1_df$padj < 0.05 & abs(res_exp1_df$log2FoldChange) >= 1, na.rm = TRUE)
n_up <- sum(res_exp1_df$padj < 0.05 & res_exp1_df$log2FoldChange > 1, na.rm = TRUE)
n_down <- sum(res_exp1_df$padj < 0.05 & res_exp1_df$log2FoldChange < -1, na.rm = TRUE)

cat("DE genes (padj < 0.05):", n_de, "\n")
cat("  Upregulated:", n_up, "\n")
cat("  Downregulated:", n_down, "\n")
cat("\nTop 10 DE genes (by padj):\n")
print(head(res_exp1_df[, c("gene", "log2FoldChange", "padj")], 10))

write_csv(res_exp1_df, "q1_exp1_0vs15_activin_240min.csv")

# Volcano plot
 res_exp1_df$sig <- ifelse(res_exp1_df$padj < 0.05 & abs(res_exp1_df$log2FoldChange) > 1.0, "DE", "Not significant")

p <- ggplot(res_exp1_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("DE" = "#B2182B", "Not significant" = "#CCCCCC"), name = "") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#666666", size = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#666666", size = 0.4) +
  labs(x = "log₂(Fold Change)", y = "−log₁₀(p-value)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10)
  )
  print(p)
  cat("\nSaved: q1_volcano_exp1.pdf\n")


# ============================================================================
# ANALYSIS: 0 vs 15 ng/ml DMSO - Experiment 2
# ============================================================================

cat("\n========================================\n")
cat("EXPERIMENT 2: 0 vs 15 ng/ml DMSO\n")
cat("========================================\n")

exp2_dmso <- metadata[
  metadata$experiment == "Exp2" &
  metadata$concentration %in% c("0ngml", "15ngml"),
]

cat("Samples:", nrow(exp2_dmso), "\n")

exp2_dmso$concentration <- factor(exp2_dmso$concentration, levels = c("0ngml", "15ngml"))

dds_exp2 <- DESeqDataSetFromMatrix(
  countData = counts_filtered[, rownames(exp2_dmso)],
  colData = exp2_dmso,
  design = ~ concentration
)
dds_exp2 <- DESeq(dds_exp2, quiet = TRUE)
res_exp2 <- results(dds_exp2, contrast = c("concentration", "15ngml", "0ngml"))

res_exp2_df <- as.data.frame(res_exp2) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

n_de <- sum(res_exp2_df$padj < 0.05 & abs(res_exp2_df$log2FoldChange) >= 1, na.rm = TRUE)
n_up <- sum(res_exp2_df$padj < 0.05 & res_exp2_df$log2FoldChange >= 1, na.rm = TRUE)
n_down <- sum(res_exp2_df$padj < 0.05 & res_exp2_df$log2FoldChange <= -1, na.rm = TRUE)

cat("DE genes (padj < 0.05):", n_de, "\n")
cat("  Upregulated:", n_up, "\n")
cat("  Downregulated:", n_down, "\n")
cat("\nTop 10 DE genes (by padj):\n")
print(head(res_exp2_df[, c("gene", "log2FoldChange", "padj")], 10))

write_csv(res_exp2_df, "q1_exp2_0vs15_dmso.csv")

# Volcano plot
res_exp2_df$sig <- ifelse(res_exp2_df$padj < 0.05 & abs(res_exp2_df$log2FoldChange) > 1, "DE", "Not significant")

p <- ggplot(res_exp2_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("DE" = "#2166AC", "Not significant" = "#CCCCCC"), name = "") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#666666", size = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#666666", size = 0.4) +
  labs(x = "log₂(Fold Change)", y = "−log₁₀(p-value)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10)
  )
print(p)
cat("\nSaved: q1_volcano_exp2.pdf\n")

# ============================================================================
# COMPARISON BETWEEN EXPERIMENTS
# ============================================================================

de_exp1 <- res_exp1_df$gene[res_exp1_df$padj < 0.05 & abs(res_exp1_df$log2FoldChange) >= 1]
de_exp2 <- res_exp2_df$gene[res_exp2_df$padj < 0.05 & abs(res_exp2_df$log2FoldChange) >= 1]
overlap <- intersect(de_exp1, de_exp2)

cat("\n========================================\n")
cat("COMPARISON: Exp1 vs Exp2\n")
cat("========================================\n")
cat("DE genes unique to Exp1:", length(setdiff(de_exp1, de_exp2)), "\n")
cat("DE genes unique to Exp2:", length(setdiff(de_exp2, de_exp1)), "\n")
cat("Shared DE genes:", length(overlap), "\n")

# ============================================================================
# COMBINED VOLCANO PLOT WITH SHARED GENES HIGHLIGHTED
# ============================================================================

library(ggrepel)
library(patchwork)

# Prepare data with shared gene annotation
res_exp1_df$category <- case_when(
  res_exp1_df$gene %in% overlap ~ "Shared",
  res_exp1_df$padj < 0.05 & abs(res_exp1_df$log2FoldChange) > 1.0 ~ "DE (Exp1 only)",
  TRUE ~ "NS"
)
res_exp1_df$category <- factor(res_exp1_df$category, levels = c("NS", "DE (Exp1 only)", "Shared"))

res_exp2_df$category <- case_when(
  res_exp2_df$gene %in% overlap ~ "Shared",
  res_exp2_df$padj < 0.05 & abs(res_exp2_df$log2FoldChange) > 1.0 ~ "DE (Exp2 only)",
  TRUE ~ "NS"
)
res_exp2_df$category <- factor(res_exp2_df$category, levels = c("NS", "DE (Exp2 only)", "Shared"))

# Color palette
colors_exp1 <- c("NS" = "#E0E0E0", "DE (Exp1 only)" = "#B2182B", "Shared" = "#7570B3")
colors_exp2 <- c("NS" = "#E0E0E0", "DE (Exp2 only)" = "#2166AC", "Shared" = "#7570B3")

# Volcano plot Exp1
# Plot shared genes first (behind), then all genes on top
p1 <- ggplot(res_exp1_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = res_exp1_df[res_exp1_df$category == "Shared", ],
             aes(color = category), size = 2, alpha = 0.7) +
  geom_point(data = res_exp1_df[res_exp1_df$category != "Shared", ],
             aes(color = category), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = colors_exp1, name = "") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#666666", linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#666666", linewidth = 0.3) +
  geom_text_repel(
    data = res_exp1_df[res_exp1_df$category == "Shared", ][1:min(15, sum(res_exp1_df$category == "Shared")), ],
    aes(label = gene), size = 2.5, max.overlaps = 20, segment.size = 0.2,
    color = "#333333", fontface = "italic"
  ) +
  labs(x = expression(log[2]~"(Fold Change)"), y = expression(-log[10]~"(adj. p-value)"),
       title = "Experiment 1", subtitle = "15 ng/ml vs 0 ng/ml Activin") +
  coord_cartesian(xlim = c(-10, 10)) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "#666666"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

# Volcano plot Exp2
# Plot shared genes first (behind), then all genes on top
p2 <- ggplot(res_exp2_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = res_exp2_df[res_exp2_df$category == "Shared", ],
             aes(color = category), size = 2, alpha = 0.7) +
  geom_point(data = res_exp2_df[res_exp2_df$category != "Shared", ],
             aes(color = category), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = colors_exp2, name = "") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#666666", linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#666666", linewidth = 0.3) +
  geom_text_repel(
    data = res_exp2_df[res_exp2_df$category == "Shared", ][1:min(15, sum(res_exp2_df$category == "Shared")), ],
    aes(label = gene), size = 2.5, max.overlaps = 20, segment.size = 0.2,
    color = "#333333", fontface = "italic"
  ) +
  labs(x = expression(log[2]~"(Fold Change)"), y = expression(-log[10]~"(adj. p-value)"),
       title = "Experiment 2", subtitle = "15 ng/ml Activin vs 0 ng/ml Activin") +
  coord_cartesian(xlim = c(-10, 10)) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "#666666"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

# Summary counts panel
n_exp1_only <- length(setdiff(de_exp1, de_exp2))
n_exp2_only <- length(setdiff(de_exp2, de_exp1))
n_shared <- length(overlap)

summary_df <- data.frame(
  category = factor(c("Exp1 only", "Shared", "Exp2 only"),
                    levels = c("Exp1 only", "Shared", "Exp2 only")),
  count = c(n_exp1_only, n_shared, n_exp2_only),
  fill = c("#B2182B", "#7570B3", "#2166AC")
)

p3 <- ggplot(summary_df, aes(x = category, y = count, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = count), vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Exp1 only" = "#B2182B", "Shared" = "#7570B3", "Exp2 only" = "#2166AC")) +
  labs(x = "", y = "Number of DE genes", title = "DE Gene Overlap") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

# Log2FC comparison scatter plot for shared genes
if(length(overlap) > 0) {
  shared_fc <- merge(
    res_exp1_df[res_exp1_df$gene %in% overlap, c("gene", "log2FoldChange", "padj")],
    res_exp2_df[res_exp2_df$gene %in% overlap, c("gene", "log2FoldChange", "padj")],
    by = "gene", suffixes = c("_Exp1", "_Exp2")
  )

  p4 <- ggplot(shared_fc, aes(x = log2FoldChange_Exp1, y = log2FoldChange_Exp2)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#999999", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "#CCCCCC", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "#CCCCCC", linewidth = 0.3) +
    geom_point(color = "#7570B3", size = 2.5, alpha = 0.8) +
    geom_text_repel(aes(label = gene), size = 2.5, max.overlaps = 20,
                    segment.size = 0.2, color = "#333333", fontface = "italic") +
    labs(x = expression(log[2]~"FC (Exp1)"), y = expression(log[2]~"FC (Exp2)"),
         title = "Shared DE Genes") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )

  # Combine all plots
  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "Differential Expression: 0 vs 15 ng/ml Activin (240 min)",
      subtitle = "Comparison between Experiment 1 and Experiment 2 ",
      theme = theme(
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#666666")
      )
    )

  # Save combined plot
  pdf("q1_combined_volcano.pdf", width = 10, height = 6, family = "Helvetica")
  print(combined)
  dev.off()
  cat("\nSaved: q1_combined_volcano.pdf\n")

  # Save shared genes CSV with log2FC from both datasets
  shared_output <- shared_fc %>%
    arrange(desc(abs(log2FoldChange_Exp1) + abs(log2FoldChange_Exp2))) %>%
    dplyr::select(gene, log2FoldChange_Exp1, padj_Exp1, log2FoldChange_Exp2, padj_Exp2)

  write_csv(shared_output, "q1_shared_de_genes.csv")
  cat("Saved: q1_shared_de_genes.csv\n")
  cat("  Contains", nrow(shared_output), "shared DE genes\n")
} else {
  cat("\nNo shared DE genes found between experiments.\n")
}

