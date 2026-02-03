# Question 2: Nodal Score - Cumulative Expression with Variance Ribbons
# Publication-ready figure for LaTeX (Helvetica)

library(DESeq2)
library(readr)
library(tidyverse)
library(readxl)

# Load data
counts_raw <- read.table("salmon.merged.gene_counts.tsv", header = TRUE, row.names = 1)[, -1]
counts_int <- round(counts_raw)
samples <- read_csv("samples.csv", show_col_types = FALSE)

# Build metadata
metadata <- data.frame(

  sample = paste0("S", samples$requests_sample_sample_id),
  treatment = sub("^\\d{8}_R\\d+_", "", samples$sample_description),
  row.names = paste0("S", samples$requests_sample_sample_id)
) %>%
  mutate(
    experiment = ifelse(grepl("DMSO|SB50", treatment), "Exp2", "Exp1"),
    condition = case_when(
      grepl("^0ngml.*DMSO", treatment) ~ "0ngml_DMSO",
      grepl("^15ngml.*DMSO", treatment) ~ "15ngml_DMSO",
      grepl("50uMSB50.*60min", treatment) ~ "SB50_60min",
      grepl("50uMSB50.*120min", treatment) ~ "SB50_120min",
      grepl("50uMSB50.*180min", treatment) ~ "SB50_180min"
    )
  ) %>%
  filter(experiment == "Exp2", !is.na(condition))

# Filter and normalize
common_samples <- intersect(colnames(counts_int), rownames(metadata))
counts_exp2 <- counts_int[, common_samples]
metadata <- metadata[common_samples, ]
counts_filtered <- counts_exp2[rowSums(counts_exp2 >= 10) >= 3, ]

metadata$condition <- factor(metadata$condition)
dds <- DESeqDataSetFromMatrix(counts_filtered, metadata, ~ 1)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

# Nodal genes
nodal_genes <- unique(tolower(na.omit(read_excel("docs/nodal-score-genes_complete.xlsx", skip = 1)$`Nodal score`)))
nodal_counts <- norm_counts[tolower(rownames(norm_counts)) %in% nodal_genes, , drop = FALSE]
n_nodal <- nrow(nodal_counts)
cat("Nodal genes:", n_nodal, "\n")

# Settings
conditions <- c("0ngml_DMSO", "15ngml_DMSO", "SB50_60min", "SB50_120min", "SB50_180min")
cond_labels <- c("0 ng/ml (Control)", "15 ng/ml Activin", "SB50 60 min", "SB50 120 min", "SB50 180 min")
names(cond_labels) <- conditions
cond_colors <- c("0ngml_DMSO" = "#BDBDBD", "15ngml_DMSO" = "#5FB358",
                 "SB50_60min" = "#FDBF6F", "SB50_120min" = "#FF7F00", "SB50_180min" = "#E31A1C")

# Order genes by control expression
gene_order <- order(rowMeans(nodal_counts[, metadata$condition == "0ngml_DMSO"]))
ordered_genes <- rownames(nodal_counts)[gene_order]

# Cumulative expression per sample
cumsum_mat <- sapply(colnames(nodal_counts), function(s) cumsum(log2(nodal_counts[ordered_genes, s] + 1)))

# Stats per condition at each rank
cumsum_stats <- expand.grid(gene_rank = 1:n_nodal, condition = conditions) %>%
  rowwise() %>%
  mutate(
    vals = list(cumsum_mat[gene_rank, rownames(metadata)[metadata$condition == condition]]),
    mean = mean(unlist(vals)),
    sd = sd(unlist(vals))
  ) %>%
  select(-vals) %>%
  ungroup() %>%
  mutate(condition = factor(condition, levels = conditions))

# Per-sample CUMULATIVE nodal score (final rank) for statistics
# This is the correct metric: total cumulative expression per sample
final_cumsum_per_sample <- data.frame(
  sample = colnames(cumsum_mat),
  condition = factor(metadata[colnames(cumsum_mat), "condition"], levels = conditions),
  cumsum_score = cumsum_mat[n_nodal, ]  # Final cumulative value
)

# ============================================================================
# STATISTICAL TESTS: Permutation test on mean curve difference (AUC-style)
# ============================================================================
# Compare the ENTIRE curve, not just the endpoint
# Test statistic: mean difference in cumulative score across all gene ranks

cat("\n--- Statistical Tests (Permutation on Curve Difference) ---\n")
cat(sprintf("Testing mean difference across all %d gene ranks (AUC-like)\n", n_nodal))
cat("Permutation test: 10,000 iterations\n\n")

set.seed(42)
n_perm <- 10000

# Function to calculate mean cumulative score per condition across all ranks
calc_mean_curve <- function(samples_idx, mat) {
  rowMeans(mat[, samples_idx, drop = FALSE])
}

ctrl_samples <- which(metadata$condition == "0ngml_DMSO")
ctrl_curve <- calc_mean_curve(ctrl_samples, cumsum_mat)

pairwise <- map_dfr(conditions[-1], function(cond) {
  test_samples <- which(metadata$condition == cond)
  test_curve <- calc_mean_curve(test_samples, cumsum_mat)

  # Observed statistic: mean difference across all ranks
  obs_diff <- mean(test_curve - ctrl_curve)

  # Also compute max deviation (KS-like)
  obs_max_dev <- max(abs(test_curve - ctrl_curve))

  # Permutation test
  all_samples <- c(ctrl_samples, test_samples)
  n_ctrl <- length(ctrl_samples)
  n_test <- length(test_samples)

  perm_diffs <- replicate(n_perm, {
    perm_idx <- sample(all_samples)
    perm_ctrl <- calc_mean_curve(perm_idx[1:n_ctrl], cumsum_mat)
    perm_test <- calc_mean_curve(perm_idx[(n_ctrl+1):(n_ctrl+n_test)], cumsum_mat)
    mean(perm_test - perm_ctrl)
  })

  # Two-tailed p-value
  p_perm <- (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1)

  sig <- case_when(p_perm < 0.001 ~ "***", p_perm < 0.01 ~ "**", p_perm < 0.05 ~ "*", TRUE ~ "ns")

  cat(sprintf("  %s:\n", cond_labels[cond]))
  cat(sprintf("    Mean curve diff = %+.2f, Max deviation = %.2f\n", obs_diff, obs_max_dev))
  cat(sprintf("    Permutation p = %.4f %s\n", p_perm, sig))

  tibble(condition = cond, mean_diff = obs_diff, max_dev = obs_max_dev, p_perm = p_perm, sig = sig)
})

write_csv(pairwise, "q2_statistical_tests.csv")

# Kruskal-Wallis on final scores (for caption, less important now)
kw <- kruskal.test(cumsum_score ~ condition, data = final_cumsum_per_sample)

# Build annotation data - place bars at end of lines, vertically stacked
annot <- pairwise %>%
  filter(sig != "ns") %>%
  arrange(desc(abs(mean_diff))) %>%
  mutate(
    y_test = final$mean[match(condition, final$condition)],
    bar_x = n_nodal * 1.02  # All bars at same x, just past the lines
  )

# Plot
p <- ggplot(cumsum_stats, aes(x = gene_rank, y = mean, color = condition, fill = condition)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = cond_colors, labels = cond_labels, name = "") +
  scale_fill_manual(values = cond_colors, labels = cond_labels, name = "") +
  labs(
    x = bquote("Nodal Score Genes (ranked, " * italic(n) * " = " * .(n_nodal) * ")"),
    y = expression(Sigma ~ log[2](counts + 1)),
    caption = "Shading: ±1 SD | Permutation test (10,000 iter) on mean curve difference vs control"
  ) +
  coord_cartesian(xlim = c(1, n_nodal * 1.15)) +
  theme_bw(base_size = 11, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 11),
    plot.caption = element_text(size = 8, hjust = 0, color = "#666666"),
    plot.margin = margin(10, 40, 10, 10)
  )

# Add significance annotations as colored brackets with asterisks
# Each bracket connects control line endpoint to test condition endpoint
bar_width <- 0.4
for(i in seq_len(nrow(annot))) {
  x_offset <- (i - 1) * 1.2  # Small horizontal offset for each bar
  bar_x <- n_nodal + 0.5 + x_offset

  # Use condition color for the bracket
  bar_color <- cond_colors[annot$condition[i]]

  p <- p +
    # Vertical bar
    annotate("segment", x = bar_x, xend = bar_x,
             y = ctrl_final, yend = annot$y_test[i],
             color = bar_color, linewidth = 0.8) +
    # Bottom cap (at control)
    annotate("segment", x = bar_x - bar_width, xend = bar_x,
             y = ctrl_final, yend = ctrl_final,
             color = bar_color, linewidth = 0.8) +
    # Top cap (at test condition)
    annotate("segment", x = bar_x - bar_width, xend = bar_x,
             y = annot$y_test[i], yend = annot$y_test[i],
             color = bar_color, linewidth = 0.8) +
    # Asterisks - larger, bold, black
    annotate("text", x = bar_x + 0.3, y = (ctrl_final + annot$y_test[i]) / 2,
             label = annot$sig[i], size = 5, hjust = 0, fontface = "bold", color = "black")
}

# Save
pdf("q2_nodal_score_cumulative.pdf", width = 10, height = 5, family = "Helvetica")
print(p)
dev.off()

png("q2_nodal_score_cumulative.png", width = 10, height = 5, units = "in", res = 300)
print(p)
dev.off()

write_csv(final_cumsum_per_sample, "q2_nodal_scores_per_sample.csv")
cat("\nSaved: q2_nodal_score_cumulative.pdf/png\n")
