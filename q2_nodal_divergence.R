# ============================================================================
# Q2 NODAL SCORE DIVERGENCE DRIVERS
# ============================================================================
#
# Which genes drive the difference between the 15ngml Activin (green)
# and SB50 180min (red) lines in the nodal score cumulative plot?
#
# For each nodal gene, compute the per-gene contribution to divergence
# at the final rank, and rank genes by their divergence contribution.
#

source("preprocess.R")
library(readxl)
library(patchwork)

# ============================================================================
# 1. PREPARE DATA (Exp2 without outlier — from preprocess.R)
# ============================================================================

metadata_q2 <- metadata %>%
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

counts_exp2 <- counts_filtered[, rownames(metadata_q2)]
counts_q2 <- counts_exp2[rowSums(counts_exp2 >= 10) >= 3, ]

metadata_q2$condition <- factor(metadata_q2$condition)
dds <- DESeqDataSetFromMatrix(counts_q2, metadata_q2, ~ 1)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

# Nodal genes
nodal_genes <- unique(tolower(na.omit(read_excel("docs/nodal-score-genes_complete.xlsx",
                                                  skip = 1)$`Nodal score`)))
nodal_counts <- norm_counts[tolower(rownames(norm_counts)) %in% nodal_genes, , drop = FALSE]
n_nodal <- nrow(nodal_counts)

conditions <- c("0ngml_DMSO", "15ngml_DMSO", "SB50_60min", "SB50_120min", "SB50_180min")
cond_labels <- c("0 ng/ml (Control)", "15 ng/ml Activin",
                 "SB50 60 min", "SB50 120 min", "SB50 180 min")
names(cond_labels) <- conditions
cond_colors <- c("0ngml_DMSO" = "#BDBDBD", "15ngml_DMSO" = "#5FB358",
                 "SB50_60min" = "#FDBF6F", "SB50_120min" = "#FF7F00",
                 "SB50_180min" = "#E31A1C")

# ============================================================================
# 2. PER-GENE DIVERGENCE BETWEEN CONDITIONS
# ============================================================================

cat("\n========== NODAL SCORE DIVERGENCE DRIVERS ==========\n")

# For each gene, compute mean log2(count+1) per condition
gene_means <- map_dfr(rownames(nodal_counts), function(gene) {
  map_dfr(conditions, function(cond) {
    samples <- rownames(metadata_q2)[metadata_q2$condition == cond]
    vals <- log2(nodal_counts[gene, samples] + 1)
    tibble(
      gene = gene,
      condition = cond,
      mean_log2 = mean(vals),
      sd_log2 = sd(vals)
    )
  })
})

# Divergence: 15ngml_DMSO - SB50_180min (what makes Activin line higher than SB50 180)
divergence_activin_sb180 <- gene_means %>%
  filter(condition %in% c("15ngml_DMSO", "SB50_180min")) %>%
  pivot_wider(id_cols = gene, names_from = condition, values_from = mean_log2) %>%
  mutate(
    divergence = `15ngml_DMSO` - SB50_180min,
    abs_divergence = abs(divergence),
    direction = ifelse(divergence > 0, "Higher in Activin", "Higher in SB50 180")
  ) %>%
  arrange(desc(abs_divergence))

cat("\nTop divergence drivers (Activin vs SB50 180min):\n")
print(head(divergence_activin_sb180 %>%
  dplyr::select(gene, divergence, direction) %>%
  mutate(divergence = round(divergence, 3)), 15))

write_csv(divergence_activin_sb180, results_path("q2_nodal_divergence_activin_vs_sb180.csv"))

# Also compute divergence for all SB50 timepoints
divergence_all <- map_dfr(c("SB50_60min", "SB50_120min", "SB50_180min"), function(sb_cond) {
  gene_means %>%
    filter(condition %in% c("15ngml_DMSO", sb_cond)) %>%
    pivot_wider(id_cols = gene, names_from = condition, values_from = mean_log2) %>%
    mutate(
      sb_condition = sb_cond,
      divergence = `15ngml_DMSO` - .data[[sb_cond]],
      abs_divergence = abs(divergence)
    ) %>%
    dplyr::select(gene, sb_condition, divergence, abs_divergence)
})

write_csv(divergence_all, results_path("q2_nodal_divergence_all_sb50.csv"))

# ============================================================================
# 3. DIVERGENCE BARPLOT
# ============================================================================

# Plot: per-gene divergence between Activin and SB50 180min
top_n <- min(n_nodal, 25)
top_div <- head(divergence_activin_sb180, top_n)

p_divergence <- ggplot(top_div, aes(x = reorder(gene, divergence), y = divergence,
                                     fill = direction)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Higher in Activin" = "#5FB358",
                                "Higher in SB50 180" = "#E31A1C"),
                    name = "") +
  labs(
    x = "", y = "Divergence: log2(Activin) − log2(SB50 180min)",
    title = "Nodal score gene divergence drivers",
    subtitle = "Genes ranked by contribution to Activin vs SB50 180min separation",
    caption = "Positive = gene more expressed with Activin | Negative = gene more expressed with SB50"
  ) +
  theme_minimal(base_size = 11, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "italic", size = 9),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50"),
    plot.caption = element_text(size = 8, hjust = 0, color = "#666666")
  )

ggsave(results_path("q2_nodal_divergence_drivers.pdf"), p_divergence, width = 8, height = 7)
cat("Saved:", results_path("q2_nodal_divergence_drivers.pdf"), "\n")

# ============================================================================
# 4. PER-GENE EXPRESSION ACROSS ALL CONDITIONS (individual gene panels)
# ============================================================================

# Plot top 12 divergence genes individually
top12 <- head(divergence_activin_sb180$gene, 12)

gene_expr_long <- gene_means %>%
  filter(gene %in% top12) %>%
  mutate(
    condition = factor(condition, levels = conditions),
    gene = factor(gene, levels = top12)
  )

p_gene_panels <- ggplot(gene_expr_long, aes(x = condition, y = mean_log2,
                                              fill = condition)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_log2 - sd_log2, ymax = mean_log2 + sd_log2),
                width = 0.3, linewidth = 0.4) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = cond_colors, labels = cond_labels, guide = "none") +
  scale_x_discrete(labels = c("Ctrl", "Act", "60", "120", "180")) +
  labs(x = "", y = "mean log2(counts + 1)",
       title = "Top divergence drivers: per-gene expression",
       subtitle = "Genes driving the separation between Activin and SB50 180min lines") +
  theme_minimal(base_size = 9, base_family = "Helvetica") +
  theme(
    strip.text = element_text(face = "bold.italic", size = 9),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50")
  )

ggsave(results_path("q2_nodal_divergence_gene_panels.pdf"), p_gene_panels, width = 10, height = 8)
cat("Saved:", results_path("q2_nodal_divergence_gene_panels.pdf"), "\n")

cat("\n========== NODAL DIVERGENCE ANALYSIS COMPLETE ==========\n")
