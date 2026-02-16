# PCA for Experiments 1 & 2 + GO enrichment on Exp2 PC1 drivers
# Combined figure: A (Exp1 PCA) | B (Exp2 PCA biplot + GO inset)

source("preprocess.R")
library(clusterProfiler)
library(org.Dr.eg.db)
library(patchwork)

# NOTE: PCA uses metadata_all / counts_int_all (from preprocess.R) so that the
# outlier S343239 is *visible* in the PCA — the whole point is to show it.
# GO enrichment (section 4) uses the clean counts_filtered (without outlier).

# ============================================================================
# 1. EXPERIMENT 1 PCA
# ============================================================================
cat("\n========== EXPERIMENT 1 PCA ==========\n")

exp1_metadata <- metadata_all %>%
  filter(experiment == "Exp1", !is.na(concentration), !is.na(time_min))

common_exp1 <- intersect(colnames(counts_int_all), rownames(exp1_metadata))
counts_exp1 <- counts_int_all[, common_exp1]
counts_exp1 <- counts_exp1[rowSums(counts_exp1 >= 10) >= 3, ]
exp1_metadata <- exp1_metadata[common_exp1, ]

dds1 <- DESeqDataSetFromMatrix(counts_exp1, exp1_metadata, ~ 1)
dds1 <- estimateSizeFactors(dds1)
vsd1 <- vst(dds1, blind = TRUE)

pca1 <- prcomp(t(assay(vsd1)), center = TRUE, scale. = FALSE)
pct1 <- round(100 * summary(pca1)$importance[2, 1:2], 1)

pca1_df <- data.frame(
  PC1 = pca1$x[, 1], PC2 = pca1$x[, 2],
  concentration = exp1_metadata$concentration,
  time_min = factor(exp1_metadata$time_min),
  sample = rownames(exp1_metadata)
)

exp1_conc_colors <- c(
  "0ngml" = "#BDBDBD", "5ngml" = "#A6D96A", "10ngml" = "#4DAF4A", "15ngml" = "#1B7837"
)
exp1_time_shapes <- c("60" = 16, "120" = 17, "180" = 15, "240" = 18)

p_exp1 <- ggplot(pca1_df, aes(x = PC1, y = PC2, color = concentration, shape = time_min)) +
  geom_point(size = 2.5, alpha = 0.85) +
  stat_ellipse(aes(group = concentration), level = 0.68, linewidth = 0.4,
               linetype = "dashed", show.legend = FALSE) +
  scale_color_manual(values = exp1_conc_colors, name = "Concentration") +
  scale_shape_manual(values = exp1_time_shapes, name = "Time (min)") +
  labs(
    x = paste0("PC1 (", pct1[1], "%)"),
    y = paste0("PC2 (", pct1[2], "%)")
  ) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    axis.title = element_text(size = 10)
  ) +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

# ============================================================================
# 2. EXPERIMENT 2 PCA
# ============================================================================
cat("\n========== EXPERIMENT 2 PCA ==========\n")

metadata_q2 <- metadata_all %>%
  mutate(
    condition = case_when(
      concentration == "0ngml_DMSO"  ~ "0ngml_DMSO",
      concentration == "15ngml_DMSO" ~ "15ngml_DMSO",
      concentration == "SB50" & time_min == 60  ~ "SB50_60min",
      concentration == "SB50" & time_min == 120 ~ "SB50_120min",
      concentration == "SB50" & time_min == 180 ~ "SB50_180min"
    )
  ) %>%
  filter(!is.na(condition))

counts_exp2 <- counts_int_all[, rownames(metadata_q2)]
counts_q2 <- counts_exp2[rowSums(counts_exp2 >= 10) >= 3, ]
metadata_q2$condition <- factor(metadata_q2$condition,
  levels = c("0ngml_DMSO", "15ngml_DMSO", "SB50_60min", "SB50_120min", "SB50_180min"))

dds2 <- DESeqDataSetFromMatrix(counts_q2, metadata_q2, ~ 1)
dds2 <- estimateSizeFactors(dds2)
vsd2 <- vst(dds2, blind = TRUE)

pca2 <- prcomp(t(assay(vsd2)), center = TRUE, scale. = FALSE)
pct2 <- round(100 * summary(pca2)$importance[2, 1:2], 1)

conditions <- c("0ngml_DMSO", "15ngml_DMSO", "SB50_60min", "SB50_120min", "SB50_180min")
cond_labels <- c("0 ng/ml (Control)", "15 ng/ml Activin",
                 "SB50 60 min", "SB50 120 min", "SB50 180 min")
names(cond_labels) <- conditions
cond_colors <- c("0ngml_DMSO" = "#BDBDBD", "15ngml_DMSO" = "#5FB358",
                 "SB50_60min" = "#FDBF6F", "SB50_120min" = "#FF7F00",
                 "SB50_180min" = "#E31A1C")

pca2_df <- data.frame(
  PC1 = pca2$x[, 1], PC2 = pca2$x[, 2],
  condition = metadata_q2$condition,
  sample = rownames(metadata_q2)
)

p_exp2 <- ggplot(pca2_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 2.5, alpha = 0.85) +
  stat_ellipse(level = 0.68, linewidth = 0.4, linetype = "dashed",
               show.legend = FALSE) +
  scale_color_manual(values = cond_colors, labels = cond_labels, name = "") +
  labs(
    x = paste0("PC1 (", pct2[1], "%)"),
    y = paste0("PC2 (", pct2[2], "%)")
  ) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    axis.title = element_text(size = 10)
  ) +
  guides(color = guide_legend(nrow = 2))

# Label S343239 outlier
outlier <- pca2_df %>% filter(sample == "S343239")
if (nrow(outlier) > 0) {
  p_exp2 <- p_exp2 + ggrepel::geom_text_repel(
    data = outlier, aes(label = sample),
    size = 2.5, color = "black", fontface = "bold",
    nudge_x = 2, nudge_y = 2,
    segment.color = "#999999", segment.size = 0.3
  )
}

# ============================================================================
# 3. EXP2 PC1 LOADINGS
# ============================================================================
cat("\n========== EXP2 PC1 DRIVER GENES ==========\n")

loadings <- data.frame(
  gene = rownames(pca2$rotation),
  PC1 = pca2$rotation[, 1],
  PC2 = pca2$rotation[, 2]
) %>%
  mutate(
    abs_PC1 = abs(PC1),
    abs_PC2 = abs(PC2),
    pc1_specificity = PC1^2 / (PC1^2 + PC2^2)
  ) %>%
  arrange(desc(abs_PC1))

cat("=== Top 20 genes by |PC1| loading ===\n")
print(head(loadings[, c("gene", "PC1", "PC2", "pc1_specificity")], 20))

pc1_specific <- loadings %>%
  filter(pc1_specificity > 0.8) %>%
  arrange(desc(abs_PC1))

n_top <- ceiling(nrow(loadings) * 0.05)
top_genes <- pc1_specific$gene[1:min(n_top, nrow(pc1_specific))]
cat(sprintf("PC1-specific genes (>80%%): %d | Top %d for GO\n",
            nrow(pc1_specific), length(top_genes)))

write_csv(loadings, results_path("q2_pca_loadings.csv"))

# ============================================================================
# 4. GO ENRICHMENT (Exp2 PC1 drivers): MF, BP, CC
# ============================================================================
cat("\n========== GO ENRICHMENT (Exp2 PC1 drivers) ==========\n")

gene_map <- bitr(top_genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Dr.eg.db)
bg_map  <- bitr(rownames(counts_q2), fromType = "SYMBOL", toType = "ENTREZID",
                OrgDb = org.Dr.eg.db)
cat(sprintf("Mapped %d/%d genes | Background %d/%d\n",
            nrow(gene_map), length(top_genes), nrow(bg_map), nrow(counts_q2)))

run_go <- function(ont_name) {
  enrichGO(
    gene = gene_map$ENTREZID,
    universe = bg_map$ENTREZID,
    OrgDb = org.Dr.eg.db,
    ont = ont_name,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
}

ego_mf <- run_go("MF")
ego_bp <- run_go("BP")
ego_cc <- run_go("CC")

for (info in list(
  list(ego = ego_mf, name = "Molecular Function"),
  list(ego = ego_bp, name = "Biological Process"),
  list(ego = ego_cc, name = "Cellular Component")
)) {
  cat(sprintf("\n=== GO %s ===\n", info$name))
  df <- as.data.frame(info$ego)
  if (nrow(df) > 0) {
    print(head(df[, c("Description", "GeneRatio", "p.adjust", "Count")], 10))
  } else {
    cat("No significant terms\n")
  }
}

write_csv(as.data.frame(ego_mf), results_path("q2_pca_drivers_go_mf.csv"))
write_csv(as.data.frame(ego_bp), results_path("q2_pca_drivers_go_bp.csv"))
write_csv(as.data.frame(ego_cc), results_path("q2_pca_drivers_go_cc.csv"))

# ============================================================================
# 5. COMBINED FIGURE: A (Exp1 PCA) | B (Exp2 biplot + GO inset)
# ============================================================================
cat("\n========== COMBINED FIGURE ==========\n")

# -- Build GO summary text for inset in panel B --
go_lines <- c()
for (info in list(
  list(ego = ego_bp, tag = "BP"),
  list(ego = ego_mf, tag = "MF"),
  list(ego = ego_cc, tag = "CC")
)) {
  df <- as.data.frame(info$ego)
  if (nrow(df) > 0) {
    top3 <- head(df, 3)
    for (i in seq_len(nrow(top3))) {
      desc <- top3$Description[i]
      # Truncate long descriptions
      if (nchar(desc) > 30) desc <- paste0(substr(desc, 1, 28), "...")
      padj <- formatC(top3$p.adjust[i], format = "e", digits = 1)
      go_lines <- c(go_lines,
        sprintf("%s: %s (p=%s)", info$tag, desc, padj))
    }
  }
}

go_text <- if (length(go_lines) > 0) {
  paste0("GO enrichment (PC1 drivers)\n",
         paste(go_lines, collapse = "\n"))
} else {
  NULL
}

# -- Panel B: Exp2 biplot with loading arrows + GO inset --
top15 <- pc1_specific[1:min(15, nrow(pc1_specific)), ]
arrow_scale <- max(abs(pca2_df$PC1)) / max(sqrt(top15$PC1^2 + top15$PC2^2)) * 0.6

p_biplot <- ggplot(pca2_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_segment(data = top15,
               aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
               arrow = arrow(length = unit(0.12, "cm")),
               color = "#AAAAAA", linewidth = 0.35, inherit.aes = FALSE) +
  ggrepel::geom_text_repel(
    data = top15,
    aes(x = PC1 * arrow_scale, y = PC2 * arrow_scale, label = gene),
    size = 2.2, color = "#555555", fontface = "italic", inherit.aes = FALSE,
    max.overlaps = 20, segment.size = 0.2
  ) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = cond_colors, labels = cond_labels, name = "") +
  labs(
    x = paste0("PC1 (", pct2[1], "%)"),
    y = paste0("PC2 (", pct2[2], "%)")
  ) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    axis.title = element_text(size = 10)
  ) +
  guides(color = guide_legend(nrow = 2))

# Label S343239 outlier on biplot
if (nrow(outlier) > 0) {
  p_biplot <- p_biplot + ggrepel::geom_text_repel(
    data = outlier, aes(label = sample),
    size = 2.5, color = "black", fontface = "bold",
    nudge_x = 2, nudge_y = 2,
    segment.color = "#999999", segment.size = 0.3
  )
}

# Add GO text box inset (bottom-right corner)
if (!is.null(go_text)) {
  p_biplot <- p_biplot +
    annotate("label",
             x = Inf, y = -Inf, hjust = 1, vjust = 0,
             label = go_text,
             size = 2.1, family = "Helvetica",
             lineheight = 0.9,
             fill = alpha("white", 0.88),
             linewidth = 0.3,        # border thickness
             color = "#333333",
             label.padding = unit(0.3, "lines"),
             label.r = unit(0.15, "lines"))
}

# -- Combine A + B side by side --
p_combined <- p_exp1 + p_biplot +
  plot_annotation(
    tag_levels = "A",
    caption = "VST-transformed counts | 68% confidence ellipses | Arrows: top PC1-specific genes (>80% PC1 variance)"
  ) &
  theme(
    plot.tag     = element_text(size = 14, face = "bold"),
    plot.caption = element_text(size = 7, hjust = 0, color = "#666666")
  )

pdf(results_path("pca_combined.pdf"), width = 14, height = 6.5, family = "Helvetica")
print(p_combined)
dev.off()
cat("Saved:", results_path("pca_combined.pdf"), "\n")

# -- Also save standalone biplot --
pdf(results_path("q2_pca_biplot.pdf"), width = 8, height = 6.5, family = "Helvetica")
print(p_biplot)
dev.off()
cat("Saved:", results_path("q2_pca_biplot.pdf"), "\n")

# ============================================================================
# 6. FULL GO DOTPLOT (separate, all terms, vertically faceted)
# ============================================================================
n_show <- 10

go_combined_full <- bind_rows(
  as.data.frame(ego_bp) %>% mutate(Ontology = "Biological Process")  %>% head(n_show),
  as.data.frame(ego_mf) %>% mutate(Ontology = "Molecular Function") %>% head(n_show),
  as.data.frame(ego_cc) %>% mutate(Ontology = "Cellular Component")  %>% head(n_show)
)

if (nrow(go_combined_full) > 0) {
  go_combined_full <- go_combined_full %>%
    mutate(
      GeneRatioNum = sapply(GeneRatio, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]]); parts[1] / parts[2]
      }),
      Ontology = factor(Ontology,
        levels = c("Biological Process", "Molecular Function", "Cellular Component")),
      Description = forcats::fct_reorder(Description, -log10(p.adjust))
    )

  p_go_all <- ggplot(go_combined_full,
                     aes(x = GeneRatioNum, y = Description,
                         size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#B2182B", high = "#2166AC",
                         name = "p.adjust", trans = "log10") +
    scale_size_continuous(range = c(2, 6), name = "Count") +
    facet_wrap(~ Ontology, scales = "free_y", ncol = 1) +
    labs(x = "Gene Ratio", y = NULL,
         caption = "Exp2 PC1-specific genes (>80% PC1 variance) | Top 10 per ontology") +
    theme_bw(base_size = 11, base_family = "Helvetica") +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y  = element_text(size = 8),
      strip.text   = element_text(size = 10, face = "bold"),
      legend.position = "right",
      plot.caption = element_text(size = 7, hjust = 0, color = "#666666")
    )

  pdf_h <- max(6, nrow(go_combined_full) * 0.3 + 3)
  pdf(results_path("q2_pca_drivers_go_all.pdf"), width = 9, height = pdf_h,
      family = "Helvetica")
  print(p_go_all)
  dev.off()
  cat("Saved:", results_path("q2_pca_drivers_go_all.pdf"), "\n")
}

cat("\nDone.\n")
