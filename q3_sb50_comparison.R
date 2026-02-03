# Question 3: Does SB50 block Activin-induced gene expression changes?
#
# EXPERIMENTAL DESIGN (Exp2):
#   - All samples collected at 240 min (4h)
#   - SB50 samples = SB50 + Activin (SB50 added at different times: 60, 120, 180 min)
#   - 15ngml_DMSO = Activin only (no SB50) - "what Activin does" control
#   - 0ngml_DMSO = No Activin, no SB50 - baseline control
#
# QUESTION: Which Activin-induced genes are blocked by SB50?
#
# COMPARISONS (all within Exp2):
#   1. Activin effect = 15ngml_DMSO vs 0ngml_DMSO (what does Activin do?)
#   2. SB50+Activin effect = SB50_Xmin vs 0ngml_DMSO (what does adding SB50 do?)
#
# CATEGORIES:
#   BLOCKED:      Activin changes gene, SB50+Activin does NOT
#                 → SB50 prevents the Activin-induced change
#   NOT BLOCKED:  Both change gene in same direction
#                 → Activin effect persists despite SB50
#   REVERSED:     Activin UP but SB50+Activin DOWN (or vice versa)
#                 → SB50 overcorrects
#   SB50-SPECIFIC: Only SB50+Activin changes gene, not Activin alone
#                 → Off-target effect of SB50

library(DESeq2)
library(readr)
library(tidyverse)
library(ggplot2)
library(patchwork)

# ============================================================================
# DATA PREPARATION
# ============================================================================

counts_raw <- read.table("salmon.merged.gene_counts.tsv", header = TRUE, row.names = 1)
counts_raw <- counts_raw[, -1]
counts_int <- round(counts_raw)

samples <- read_csv("samples.csv", show_col_types = FALSE)
metadata <- data.frame(
  sample = paste0("S", samples$requests_sample_sample_id),
  treatment = sub("^\\d{8}_R\\d+_", "", samples$sample_description),
  row.names = paste0("S", samples$requests_sample_sample_id)
)

metadata$experiment <- ifelse(grepl("DMSO|SB50", metadata$treatment), "Exp2", "Exp1")
metadata$concentration <- case_when(
  grepl("^50uMSB50", metadata$treatment) ~ "SB50",
  grepl("DMSO", metadata$treatment) & grepl("^15ngml", metadata$treatment) ~ "15ngml_DMSO",
  grepl("DMSO", metadata$treatment) & grepl("^0ngml", metadata$treatment) ~ "0ngml_DMSO",
  grepl("^15ngmlActivin", metadata$treatment) ~ "15ngml_Activin",
  grepl("^0ngmlActivin", metadata$treatment) & !grepl("DMSO", metadata$treatment) ~ "0ngml_Activin",
  TRUE ~ "other"
)
metadata$time_min <- as.numeric(str_extract(metadata$treatment, "\\d+(?=min$)"))

common_samples <- intersect(colnames(counts_int), rownames(metadata))
counts_int <- counts_int[, common_samples]
metadata <- metadata[common_samples, ]
keep <- rowSums(counts_int >= 10) >= 3
counts_filtered <- counts_int[keep, ]

# ============================================================================
# DEFINE SAMPLE GROUPS (all Exp2 - collected at 240min)
# ============================================================================

# Exp2: SB50 + Activin samples (SB50 added at different times)
sb50_60min  <- metadata[metadata$concentration == "SB50" & metadata$time_min == 60, ]
sb50_120min <- metadata[metadata$concentration == "SB50" & metadata$time_min == 120, ]
sb50_180min <- metadata[metadata$concentration == "SB50" & metadata$time_min == 180, ]

# Exp2: Activin only control (15ngml Activin + DMSO, no SB50)
activin_only <- metadata[metadata$concentration == "15ngml_DMSO", ]

# Exp2: Baseline control (0ngml Activin + DMSO, no SB50)
baseline_ctrl <- metadata[metadata$concentration == "0ngml_DMSO", ]

cat("=== SAMPLE COUNTS (all Exp2, collected at 240min) ===\n")
cat("SB50+Activin: 60min=", nrow(sb50_60min), " 120min=", nrow(sb50_120min), " 180min=", nrow(sb50_180min), "\n")
cat("  (timepoint = when SB50 was added)\n")
cat("Activin only (15ngml_DMSO): ", nrow(activin_only), "\n")
cat("Baseline (0ngml_DMSO): ", nrow(baseline_ctrl), "\n\n")

# ============================================================================
# DESeq2 FUNCTIONS
# ============================================================================

run_deseq <- function(test_samples, ref_samples, counts_mat) {
  combined <- rbind(
    data.frame(test_samples, group = "test"),
    data.frame(ref_samples, group = "ref")
  )
  combined$group <- factor(combined$group, levels = c("ref", "test"))

  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat[, rownames(combined)],
    colData = combined,
    design = ~ group
  )
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", "test", "ref"))

  as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    arrange(padj)
}

print_summary <- function(res_df, name, lfc_thresh = 1.5) {
  sig <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= lfc_thresh
  n_de <- sum(sig)
  n_up <- sum(sig & res_df$log2FoldChange > 0)
  n_down <- sum(sig & res_df$log2FoldChange < 0)
  cat(sprintf("%s: %d DE (↑%d ↓%d)\n", name, n_de, n_up, n_down))
  data.frame(comparison = name, de_genes = n_de, up = n_up, down = n_down)
}

# ============================================================================
# COMPARISONS (all within Exp2, all vs baseline control)
# ============================================================================

# COMPARISON 1: What does Activin do? (Activin only vs Baseline)
cat("=== ACTIVIN EFFECT: 15ngml Activin vs Baseline (0ngml) ===\n")
res_activin_vs_baseline <- run_deseq(activin_only, baseline_ctrl, counts_filtered)
print_summary(res_activin_vs_baseline, "Activin_vs_Baseline")

# COMPARISON 2: What does SB50+Activin do? (SB50+Activin vs Baseline)
cat("\n=== SB50+ACTIVIN EFFECT vs BASELINE ===\n")
res_sb50_60_vs_baseline <- run_deseq(sb50_60min, baseline_ctrl, counts_filtered)
sum1 <- print_summary(res_sb50_60_vs_baseline, "SB50_60min_vs_Baseline")

res_sb50_120_vs_baseline <- run_deseq(sb50_120min, baseline_ctrl, counts_filtered)
sum2 <- print_summary(res_sb50_120_vs_baseline, "SB50_120min_vs_Baseline")

res_sb50_180_vs_baseline <- run_deseq(sb50_180min, baseline_ctrl, counts_filtered)
sum3 <- print_summary(res_sb50_180_vs_baseline, "SB50_180min_vs_Baseline")

# ============================================================================
# SAVE RESULTS
# ============================================================================

write_csv(res_activin_vs_baseline, "q3_Activin_vs_Baseline.csv")
write_csv(res_sb50_60_vs_baseline, "q3_SB50_60min_vs_Baseline.csv")
write_csv(res_sb50_120_vs_baseline, "q3_SB50_120min_vs_Baseline.csv")
write_csv(res_sb50_180_vs_baseline, "q3_SB50_180min_vs_Baseline.csv")

summary_all <- bind_rows(sum1, sum2, sum3)
write_csv(summary_all, "q3_summary.csv")
print(summary_all)

# ============================================================================
# BLOCKING ANALYSIS
# ============================================================================
#
# For each gene, we compare:
#   - Activin effect = 15ngml_DMSO vs 0ngml_DMSO (does Activin change this gene?)
#   - SB50+Activin effect = SB50_Xmin vs 0ngml_DMSO (does SB50+Activin change this gene?)
#
# CATEGORIES:
#   BLOCKED (blue):      Activin changes gene, but SB50+Activin does NOT
#                        → SB50 prevents the Activin-induced change (gene stays at baseline)
#   NOT BLOCKED (red):   Both change gene in SAME direction vs baseline
#                        → The Activin effect persists despite SB50
#   REVERSED (purple):   Activin UP but SB50+Activin DOWN (or vice versa) vs baseline
#                        → SB50 overcorrects
#   SB50-SPECIFIC (green): Only SB50+Activin changes gene, not Activin alone
#                        → Off-target effect of SB50
#   NS (grey):           Neither treatment significantly changes gene vs baseline
#

LFC_THRESH <- 1.5
PADJ_THRESH <- 0.05
is_sig <- function(padj, lfc) !is.na(padj) & padj < PADJ_THRESH & abs(lfc) >= LFC_THRESH

get_blocking_data <- function(res_activin, res_sb50, timepoint) {
  inner_join(
    res_activin %>% select(gene, lfc_activin = log2FoldChange, padj_activin = padj),
    res_sb50 %>% select(gene, lfc_sb50 = log2FoldChange, padj_sb50 = padj),
    by = "gene"
  ) %>%
    mutate(
      activin_sig = is_sig(padj_activin, lfc_activin),
      sb50_sig = is_sig(padj_sb50, lfc_sb50),
      category = case_when(
        activin_sig & !sb50_sig ~ "Blocked",
        activin_sig & sb50_sig & sign(lfc_activin) == sign(lfc_sb50) ~ "Not blocked",
        activin_sig & sb50_sig & sign(lfc_activin) != sign(lfc_sb50) ~ "Reversed",
        !activin_sig & sb50_sig ~ "SB50-specific",
        TRUE ~ "NS"
      ),
      timepoint = timepoint
    )
}

# All comparisons use the SAME Activin effect (single comparison within Exp2)
data_60 <- get_blocking_data(res_activin_vs_baseline, res_sb50_60_vs_baseline, 60)
data_120 <- get_blocking_data(res_activin_vs_baseline, res_sb50_120_vs_baseline, 120)
data_180 <- get_blocking_data(res_activin_vs_baseline, res_sb50_180_vs_baseline, 180)

cat("\n=== BLOCKING SUMMARY ===\n")
cat("(SB50 timepoint = when SB50 was added; longer time = SB50 present longer)\n\n")
for (d in list(list(data_60, 60), list(data_120, 120), list(data_180, 180))) {
  dat <- d[[1]]; tp <- d[[2]]
  n_activin <- sum(dat$activin_sig)
  cat(sprintf("SB50 added at %dmin (present for %dmin):\n", tp, 240 - tp))
  cat(sprintf("  Activin-responsive genes: %d\n", n_activin))
  cat(sprintf("  Blocked: %d (%.0f%%)\n", sum(dat$category == "Blocked"),
              100 * sum(dat$category == "Blocked") / max(1, n_activin)))
  cat(sprintf("  Not blocked: %d\n", sum(dat$category == "Not blocked")))
  cat(sprintf("  Reversed: %d\n", sum(dat$category == "Reversed")))
  cat(sprintf("  SB50-specific: %d\n\n", sum(dat$category == "SB50-specific")))
}

# Color palette
cat_colors <- c("Blocked" = "#2166AC", "Not blocked" = "#B2182B",
                "Reversed" = "#762A83", "SB50-specific" = "#1B7837", "NS" = "grey85")

# ============================================================================
# SCATTER PLOT
# ============================================================================
#
# HOW TO READ:
#   X-axis: Activin effect (Activin only vs Baseline)
#   Y-axis: SB50+Activin effect (SB50+Activin vs Baseline)
#
#   If SB50 BLOCKS Activin:
#     → Genes that Activin changes (high |X|) should NOT change with SB50+Activin (low |Y|)
#     → Points along X-axis = BLOCKED (blue)
#
#   If SB50 does NOT block:
#     → Genes change with both treatments
#     → Points along diagonal = NOT BLOCKED (red)
#

make_scatter_plot <- function(data, timepoint) {
  sb50_duration <- 240 - timepoint
  ggplot(data, aes(x = lfc_activin, y = lfc_sb50, color = category)) +
    geom_point(data = filter(data, category == "NS"), alpha = 0.1, size = 0.5) +
    geom_point(data = filter(data, category != "NS"), alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted", linewidth = 0.6) +
    scale_color_manual(
      values = cat_colors,
      breaks = c("Blocked", "Not blocked", "Reversed", "SB50-specific")
    ) +
    labs(
      title = sprintf("SB50 added at %dmin", timepoint),
      subtitle = sprintf("(SB50 present for %dmin)", sb50_duration),
      x = "Activin effect (log2FC)",
      y = "SB50+Activin effect (log2FC)",
      color = NULL
    ) +
    coord_fixed(ratio = 1, xlim = c(-10, 10), ylim = c(-10, 10)) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey50"),
      legend.position = "none"
    )
}

p_scatter_60 <- make_scatter_plot(data_60, 60)
p_scatter_120 <- make_scatter_plot(data_120, 120)
p_scatter_180 <- make_scatter_plot(data_180, 180)

# ============================================================================
# 100% STACKED BARPLOT
# ============================================================================

all_data <- bind_rows(data_60, data_120, data_180) %>%
  filter(category != "NS") %>%
  mutate(
    timepoint = factor(timepoint, levels = c(60, 120, 180)),
    sb50_duration = factor(240 - as.numeric(as.character(timepoint)))
  )

prop_data <- all_data %>%
  group_by(timepoint, category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(timepoint) %>%
  mutate(total = sum(n), pct = 100 * n / total) %>%
  ungroup() %>%
  mutate(category = factor(category, levels = c("SB50-specific", "Reversed", "Not blocked", "Blocked")))

# Chi-square test
# NOTE: This tests whether the PROPORTIONS of gene categories (Blocked, Not blocked, Reversed)
# differ across SB50 timepoints. Each "observation" is a GENE, not a sample.
# We're asking: does the distribution of blocking outcomes change depending on when SB50 was added?
contingency <- all_data %>%
  filter(category %in% c("Blocked", "Not blocked", "Reversed")) %>%
  dplyr::count(timepoint, category) %>%
  pivot_wider(names_from = category, values_from = n, values_fill = 0)

chi_test <- chisq.test(contingency[, -1])
cat("=== CHI-SQUARE TEST ===\n")
cat("Tests if category proportions differ across SB50 timepoints (unit = genes, not samples)\n")
cat(sprintf("Chi-squared = %.2f, df = %d, p = %.2e\n", chi_test$statistic, chi_test$parameter, chi_test$p.value))

p_bar <- ggplot(prop_data, aes(x = timepoint, y = pct, fill = category)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(
    aes(label = ifelse(pct > 8, sprintf("%d\n(%.0f%%)", n, pct), "")),
    position = position_stack(vjust = 0.5),
    size = 2.5, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(
    values = cat_colors,
    breaks = c("Blocked", "Not blocked", "Reversed", "SB50-specific")
  ) +
  scale_x_discrete(labels = c("60" = "60min", "120" = "120min", "180" = "180min")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  labs(
    title = "Category proportions",
    x = "SB50 added at",
    y = "% of DE genes",
    fill = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 7),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(fill = guide_legend(nrow = 2))

# ============================================================================
# COMBINED 2x2 FIGURE (equal sized panels)
# ============================================================================

# Wrap barplot to match scatter plot dimensions
p_bar_wrapped <- p_bar +
  theme(aspect.ratio = 1)

combined_fig <- (p_scatter_60 + p_scatter_120) / (p_scatter_180 + p_bar_wrapped) +
  plot_layout(widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(
    title = "Analysis 1: Does SB50 block Activin-induced gene expression changes?",
    subtitle = "All comparisons vs baseline (0ngml_DMSO). X = Activin effect, Y = SB50+Activin effect",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
    )
  )

ggsave("q3_analysis1_blocking.pdf", combined_fig, width = 8, height = 8)
cat("\nSaved: q3_analysis1_blocking.pdf\n")

# Save results
write_csv(prop_data %>% select(timepoint, category, n, total, pct), "q3_category_summary.csv")

# ============================================================================
# INTERPRETATION
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("HOW TO INTERPRET\n")
cat(strrep("=", 70), "\n\n")

cat("SCATTER PLOTS:\n")
cat("  X-axis: Activin effect (15ngml Activin vs 0ngml baseline)\n")
cat("  Y-axis: SB50+Activin effect (SB50+Activin vs 0ngml baseline)\n\n")

cat("  BLUE (Blocked): Gene changes with Activin alone, but NOT with SB50+Activin\n")
cat("    → SB50 successfully blocks this Activin-induced change\n")
cat("    → Gene stays at baseline level despite Activin being present\n\n")

cat("  RED (Not blocked): Gene changes with BOTH treatments (same direction)\n")
cat("    → Activin effect persists even with SB50 present\n")
cat("    → SB50 fails to block this gene\n\n")

cat("  PURPLE (Reversed): Opposite effects vs baseline\n")
cat("    → SB50 overcorrects (e.g., Activin UP → SB50+Activin DOWN)\n\n")

cat("  GREEN (SB50-specific): Only SB50+Activin changes gene, not Activin alone\n")
cat("    → Off-target effect of SB50\n\n")

cat("BARPLOT:\n")
cat("  Shows proportion of genes in each category\n")
cat("  More BLUE = better blocking by SB50\n")
cat("  Longer SB50 exposure (earlier addition) may show more blocking\n")

# ============================================================================
# ============================================================================
# ANALYSIS 2: Direct comparisons (your original question)
# ============================================================================
# ============================================================================
#
# Compare each SB50 condition to:
#   1. 15ngml DMSO 4h (Exp2) - Activin-only control at same collection time
#   2. 15ngml Activin at same timepoint (Exp1) - time-matched Activin
#
# DIFFERENT INSIGHT:
#   Analysis 1 asks: "Does SB50+Activin look like baseline?"
#   Analysis 2 asks: "Does SB50+Activin look like Activin-only OR like time-matched Activin?"
#
#   - If SB50 blocks Activin → SB50 samples should be DIFFERENT from Activin-only (many DE genes)
#   - If SB50 blocks Activin → SB50 samples should be SIMILAR to baseline (few DE genes)
#
# This provides a different perspective: instead of asking "what changed vs baseline",
# we ask "how similar/different are SB50 samples from Activin-treated samples?"
#

cat("\n")
cat(strrep("=", 70), "\n")
cat("ANALYSIS 2: Direct comparisons to Activin references\n")
cat(strrep("=", 70), "\n\n")

# Need Exp1 Activin samples (time-matched)
activin_exp1_60min  <- metadata[metadata$concentration == "15ngml_Activin" & metadata$time_min == 60, ]
activin_exp1_120min <- metadata[metadata$concentration == "15ngml_Activin" & metadata$time_min == 120, ]
activin_exp1_180min <- metadata[metadata$concentration == "15ngml_Activin" & metadata$time_min == 180, ]

cat("Exp1 Activin samples: 60min=", nrow(activin_exp1_60min),
    " 120min=", nrow(activin_exp1_120min),
    " 180min=", nrow(activin_exp1_180min), "\n\n")

# ============================================================================
# Comparison A: SB50 vs Activin-only (15ngml_DMSO, Exp2)
# Both collected at 240min, both have Activin - only difference is SB50
# ============================================================================

cat("=== SB50 vs ACTIVIN-ONLY (both Exp2, 240min) ===\n")
cat("If similar (few DE) → SB50 does NOT block\n")
cat("If different (many DE) → SB50 DOES block (changes the outcome)\n\n")

res_sb50_60_vs_activin_only <- run_deseq(sb50_60min, activin_only, counts_filtered)
res_sb50_120_vs_activin_only <- run_deseq(sb50_120min, activin_only, counts_filtered)
res_sb50_180_vs_activin_only <- run_deseq(sb50_180min, activin_only, counts_filtered)

print_summary(res_sb50_60_vs_activin_only, "SB50_60min vs Activin-only")
print_summary(res_sb50_120_vs_activin_only, "SB50_120min vs Activin-only")
print_summary(res_sb50_180_vs_activin_only, "SB50_180min vs Activin-only")

# ============================================================================
# Comparison B: SB50 vs Time-matched Activin (Exp1)
# SB50 collected at 240min, Exp1 Activin at matching timepoint
# ============================================================================

cat("\n=== SB50 (Exp2) vs TIME-MATCHED ACTIVIN (Exp1) ===\n")
cat("Compares SB50+Activin (collected 240min) to Activin-only (collected at SB50 timepoint)\n\n")

res_sb50_60_vs_exp1_60 <- run_deseq(sb50_60min, activin_exp1_60min, counts_filtered)
res_sb50_120_vs_exp1_120 <- run_deseq(sb50_120min, activin_exp1_120min, counts_filtered)
res_sb50_180_vs_exp1_180 <- run_deseq(sb50_180min, activin_exp1_180min, counts_filtered)

print_summary(res_sb50_60_vs_exp1_60, "SB50_60min(Exp2) vs Activin_60min(Exp1)")
print_summary(res_sb50_120_vs_exp1_120, "SB50_120min(Exp2) vs Activin_120min(Exp1)")
print_summary(res_sb50_180_vs_exp1_180, "SB50_180min(Exp2) vs Activin_180min(Exp1)")

# Save these results
write_csv(res_sb50_60_vs_activin_only, "q3_SB50_60min_vs_ActivinOnly.csv")
write_csv(res_sb50_120_vs_activin_only, "q3_SB50_120min_vs_ActivinOnly.csv")
write_csv(res_sb50_180_vs_activin_only, "q3_SB50_180min_vs_ActivinOnly.csv")
write_csv(res_sb50_60_vs_exp1_60, "q3_SB50_60min_vs_Exp1_60min.csv")
write_csv(res_sb50_120_vs_exp1_120, "q3_SB50_120min_vs_Exp1_120min.csv")
write_csv(res_sb50_180_vs_exp1_180, "q3_SB50_180min_vs_Exp1_180min.csv")

# ============================================================================
# VISUALIZATION: Scatter plots of the two comparisons
# ============================================================================
#
# X-axis: log2FC (SB50 vs Activin-only at 240min, Exp2)
# Y-axis: log2FC (SB50 vs time-matched Activin, Exp1)
#
# INTERPRETATION:
#   - If SB50 = Activin-only → X near 0
#   - If SB50 = time-matched Activin → Y near 0
#   - Points in bottom-left or top-right: SB50 differs from BOTH references
#   - Points along X-axis: SB50 differs from Activin-only but similar to Exp1
#   - Points along Y-axis: SB50 similar to Activin-only but differs from Exp1
#

get_comparison_data <- function(res_vs_activin_only, res_vs_exp1, timepoint) {
  inner_join(
    res_vs_activin_only %>% select(gene, lfc_vs_activin_only = log2FoldChange, padj_vs_activin_only = padj),
    res_vs_exp1 %>% select(gene, lfc_vs_exp1 = log2FoldChange, padj_vs_exp1 = padj),
    by = "gene"
  ) %>%
    mutate(
      sig_vs_activin_only = is_sig(padj_vs_activin_only, lfc_vs_activin_only),
      sig_vs_exp1 = is_sig(padj_vs_exp1, lfc_vs_exp1),
      category = case_when(
        sig_vs_activin_only & sig_vs_exp1 ~ "DE in both",
        sig_vs_activin_only & !sig_vs_exp1 ~ "DE vs Activin-only",
        !sig_vs_activin_only & sig_vs_exp1 ~ "DE vs Exp1",
        TRUE ~ "NS"
      ),
      timepoint = timepoint
    )
}

comp_60 <- get_comparison_data(res_sb50_60_vs_activin_only, res_sb50_60_vs_exp1_60, 60)
comp_120 <- get_comparison_data(res_sb50_120_vs_activin_only, res_sb50_120_vs_exp1_120, 120)
comp_180 <- get_comparison_data(res_sb50_180_vs_activin_only, res_sb50_180_vs_exp1_180, 180)

cat("\n=== ANALYSIS 2 SUMMARY ===\n")
for (d in list(list(comp_60, 60), list(comp_120, 120), list(comp_180, 180))) {
  dat <- d[[1]]; tp <- d[[2]]
  cat(sprintf("SB50 added at %dmin:\n", tp))
  cat(sprintf("  DE vs Activin-only (Exp2): %d\n", sum(dat$sig_vs_activin_only)))
  cat(sprintf("  DE vs time-matched Activin (Exp1): %d\n", sum(dat$sig_vs_exp1)))
  cat(sprintf("  DE in both: %d\n\n", sum(dat$category == "DE in both")))
}

# Color palette for Analysis 2
comp_colors <- c("DE in both" = "#7570B3", "DE vs Activin-only" = "#D95F02",
                 "DE vs Exp1" = "#1B9E77", "NS" = "grey85")

make_comparison_scatter <- function(data, timepoint) {
  sb50_duration <- 240 - timepoint
  n_vs_activin <- sum(data$sig_vs_activin_only)
  n_vs_exp1 <- sum(data$sig_vs_exp1)

  ggplot(data, aes(x = lfc_vs_activin_only, y = lfc_vs_exp1, color = category)) +
    geom_point(data = filter(data, category == "NS"), alpha = 0.1, size = 0.5) +
    geom_point(data = filter(data, category != "NS"), alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted", linewidth = 0.6) +
    scale_color_manual(
      values = comp_colors,
      breaks = c("DE vs Activin-only", "DE vs Exp1", "DE in both")
    ) +
    labs(
      title = sprintf("SB50 at %dmin", timepoint),
      subtitle = sprintf("DE: %d vs Act-only, %d vs Exp1", n_vs_activin, n_vs_exp1),
      x = "vs Activin-only 240min (log2FC)",
      y = sprintf("vs Activin %dmin Exp1 (log2FC)", timepoint),
      color = NULL
    ) +
    coord_fixed(ratio = 1, xlim = c(-10, 10), ylim = c(-10, 10)) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey50"),
      legend.position = "none"
    )
}

p2_scatter_60 <- make_comparison_scatter(comp_60, 60)
p2_scatter_120 <- make_comparison_scatter(comp_120, 120)
p2_scatter_180 <- make_comparison_scatter(comp_180, 180)

# Summary barplot for Analysis 2
all_comp_data <- bind_rows(comp_60, comp_120, comp_180) %>%
  filter(category != "NS") %>%
  mutate(timepoint = factor(timepoint, levels = c(60, 120, 180)))

comp_prop_data <- all_comp_data %>%
  group_by(timepoint, category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(timepoint) %>%
  mutate(total = sum(n), pct = 100 * n / total) %>%
  ungroup() %>%
  mutate(category = factor(category, levels = c("DE vs Exp1", "DE in both", "DE vs Activin-only")))

p2_bar <- ggplot(comp_prop_data, aes(x = timepoint, y = pct, fill = category)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(
    aes(label = ifelse(pct > 8, sprintf("%d\n(%.0f%%)", n, pct), "")),
    position = position_stack(vjust = 0.5),
    size = 2.5, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(
    values = comp_colors,
    breaks = c("DE vs Activin-only", "DE vs Exp1", "DE in both")
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  labs(
    title = "DE gene distribution",
    x = "SB50 added at",
    y = "% of DE genes",
    fill = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 7),
    panel.grid.major.x = element_blank(),
    aspect.ratio = 1
  ) +
  guides(fill = guide_legend(nrow = 2))

# Combined figure for Analysis 2
combined_fig2 <- (p2_scatter_60 + p2_scatter_120) / (p2_scatter_180 + p2_bar) +
  plot_layout(widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(
    title = "Analysis 2: How does SB50 compare to Activin references?",
    subtitle = "X = SB50 vs Activin-only (Exp2, 240min) | Y = SB50 vs time-matched Activin (Exp1)",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
    )
  )

ggsave("q3_analysis2_comparisons.pdf", combined_fig2, width = 8, height = 8)
cat("\nSaved: q3_analysis2_comparisons.pdf\n")

# Save comparison summary
write_csv(comp_prop_data %>% select(timepoint, category, n, total, pct), "q3_analysis2_summary.csv")

# ============================================================================
# INTERPRETATION OF ANALYSIS 2
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("HOW TO INTERPRET ANALYSIS 2\n")
cat(strrep("=", 70), "\n\n")

cat("QUESTION: How similar/different is SB50+Activin from Activin-only?\n\n")

cat("X-axis: SB50 vs Activin-only (both Exp2, collected at 240min)\n")
cat("  → Measures: Does adding SB50 change the outcome vs Activin alone?\n")
cat("  → Many DE genes = SB50 has strong effect\n\n")

cat("Y-axis: SB50 vs time-matched Activin (Exp1, collected at SB50 timepoint)\n")
cat("  → Measures: Is SB50+Activin at 240min similar to Activin at earlier time?\n")
cat("  → If similar = SB50 'freezes' gene expression at an earlier state\n\n")

cat("CATEGORIES:\n")
cat("  ORANGE (DE vs Activin-only): SB50 differs from Activin-only but not from Exp1\n")
cat("    → SB50 makes cells look like earlier Activin timepoint\n")
cat("  GREEN (DE vs Exp1): SB50 similar to Activin-only but differs from Exp1\n")
cat("    → Time effect but not SB50 effect\n")
cat("  PURPLE (DE in both): SB50 differs from both references\n")
cat("    → SB50 creates unique state (neither Activin-like nor frozen)\n\n")

cat("INSIGHT COMPARISON:\n")
cat("  Analysis 1: Asks 'does SB50 keep cells at baseline?'\n")
cat("  Analysis 2: Asks 'does SB50 make cells different from Activin-treated?'\n")
cat("  Together: Full picture of SB50's blocking mechanism\n")
