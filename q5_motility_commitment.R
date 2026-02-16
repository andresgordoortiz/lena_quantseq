# ============================================================================
# Q5: CELL MOTILITY GENE COMMITMENT ANALYSIS
# ============================================================================
#
# Hypothesis: Cell motility/migration genes are NOT committed at SB50-60min
# (consistent with no motility phenotype when SB50 is added early), but ARE
# committed at SB50-120min and 180min (consistent with protrusions/motility
# phenotype observed when SB50 is added at those later points).
#
# Approach:
#   1. Get cell motility/migration genes from GO annotations
#   2. Cross-reference with SB50 blocking data at 60, 120, 180 min
#   3. Test whether migration genes are enriched in "Blocked" at 60min
#      and enriched in "Not blocked" at 120/180min
#   4. Compare commitment profile of motility genes vs all other DEGs
#

source("preprocess.R")
library(patchwork)
library(ggrepel)
library(clusterProfiler)
library(org.Dr.eg.db)

# ============================================================================
# 1. RECOMPUTE BLOCKING DATA (same logic as q3_extended_analysis.R)
# ============================================================================

LFC_THRESH <- 1.0
PADJ_THRESH <- 0.05
is_sig <- function(padj, lfc) !is.na(padj) & padj < PADJ_THRESH & abs(lfc) >= LFC_THRESH

# Sample groups
sb50_60min  <- metadata[metadata$concentration == "SB50" & metadata$time_min == 60, ]
sb50_120min <- metadata[metadata$concentration == "SB50" & metadata$time_min == 120, ]
sb50_180min <- metadata[metadata$concentration == "SB50" & metadata$time_min == 180, ]
activin_only <- metadata[metadata$concentration == "15ngml_DMSO", ]
baseline_ctrl <- metadata[metadata$concentration == "0ngml_DMSO", ]

# DE comparisons
res_activin_vs_baseline <- run_deseq(activin_only, baseline_ctrl, counts_filtered)
res_sb50_60_vs_baseline <- run_deseq(sb50_60min, baseline_ctrl, counts_filtered)
res_sb50_120_vs_baseline <- run_deseq(sb50_120min, baseline_ctrl, counts_filtered)
res_sb50_180_vs_baseline <- run_deseq(sb50_180min, baseline_ctrl, counts_filtered)

get_blocking_data <- function(res_activin, res_sb50, timepoint) {
  inner_join(
    res_activin %>% dplyr::select(gene, lfc_activin = log2FoldChange, padj_activin = padj),
    res_sb50 %>% dplyr::select(gene, lfc_sb50 = log2FoldChange, padj_sb50 = padj),
    by = "gene"
  ) %>%
    mutate(
      activin_sig = is_sig(padj_activin, lfc_activin),
      sb50_sig = is_sig(padj_sb50, lfc_sb50),
      activin_direction = case_when(
        lfc_activin > 0 ~ "Up",
        lfc_activin < 0 ~ "Down",
        TRUE ~ "None"
      ),
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

data_60 <- get_blocking_data(res_activin_vs_baseline, res_sb50_60_vs_baseline, 60)
data_120 <- get_blocking_data(res_activin_vs_baseline, res_sb50_120_vs_baseline, 120)
data_180 <- get_blocking_data(res_activin_vs_baseline, res_sb50_180_vs_baseline, 180)

# ============================================================================
# 2. GET CELL MOTILITY/MIGRATION & ADHESION GENES FROM GO
# ============================================================================

cat("\n========== IDENTIFYING MOTILITY/MIGRATION GENES ==========\n")

get_go_genes <- function(go_ids, label) {
  entrez <- tryCatch({
    AnnotationDbi::select(org.Dr.eg.db, keys = go_ids,
                           columns = "ENTREZID", keytype = "GOALL")$ENTREZID
  }, error = function(e) character(0))
  if (length(entrez) == 0) return(character(0))
  symbols <- AnnotationDbi::select(org.Dr.eg.db,
    keys = unique(entrez), columns = "SYMBOL", keytype = "ENTREZID")$SYMBOL
  symbols <- unique(na.omit(symbols))
  # Restrict to genes in our dataset
  in_data <- intersect(symbols, rownames(counts_filtered))
  cat(sprintf("  %s: %d GO-annotated genes, %d in dataset\n",
              label, length(symbols), length(in_data)))
  in_data
}

# Cell motility / migration GO terms
motility_genes <- get_go_genes(
  c("GO:0048870",   # cell motility
    "GO:0016477",   # cell migration
    "GO:0030334",   # regulation of cell migration
    "GO:0040011"),  # locomotion
  "Cell motility/migration"
)

# Cell adhesion (for comparison)
adhesion_genes <- get_go_genes(
  c("GO:0007155"),  # cell adhesion
  "Cell adhesion"
)

# ============================================================================
# 3. CROSS-REFERENCE: MOTILITY GENES × BLOCKING CATEGORIES
# ============================================================================

cat("\n========== MOTILITY GENES × COMMITMENT STATUS ==========\n")

# Tag genes as motility or not in each blocking dataset
tag_motility <- function(blocking_data, motility, adhesion) {
  blocking_data %>%
    mutate(
      is_motility = gene %in% motility,
      is_adhesion = gene %in% adhesion,
      gene_class = case_when(
        is_motility & is_adhesion ~ "Motility + Adhesion",
        is_motility ~ "Motility",
        is_adhesion ~ "Adhesion",
        TRUE ~ "Other"
      )
    )
}

data_60_tagged <- tag_motility(data_60, motility_genes, adhesion_genes)
data_120_tagged <- tag_motility(data_120, motility_genes, adhesion_genes)
data_180_tagged <- tag_motility(data_180, motility_genes, adhesion_genes)

# Focus on Activin-responsive genes only (Blocked, Not blocked, or Reversed)
activin_cats <- c("Blocked", "Not blocked", "Reversed")

# Count motility genes per category at each timepoint
count_by_class <- function(data, tp_label) {
  data %>%
    filter(category %in% activin_cats) %>%
    group_by(gene_class, category) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(timepoint = tp_label)
}

motility_counts <- bind_rows(
  count_by_class(data_60_tagged, "SB50 at 60 min"),
  count_by_class(data_120_tagged, "SB50 at 120 min"),
  count_by_class(data_180_tagged, "SB50 at 180 min")
)

# Detailed breakdown
cat("\n--- Motility gene commitment status ---\n")
for (tp in c("SB50 at 60 min", "SB50 at 120 min", "SB50 at 180 min")) {
  cat(sprintf("\n%s:\n", tp))
  subset <- motility_counts %>%
    filter(timepoint == tp, gene_class %in% c("Motility", "Motility + Adhesion"))
  if (nrow(subset) > 0) print(as.data.frame(subset))
  else cat("  No Activin-responsive motility genes\n")
}

# ============================================================================
# 4. COMPUTE COMMITMENT PROPORTIONS: MOTILITY vs ALL OTHER DEGs
# ============================================================================

cat("\n========== COMMITMENT PROPORTIONS ==========\n")

compute_proportions <- function(data, tp_label) {
  activin_responsive <- data %>% filter(category %in% activin_cats)

  motility_subset <- activin_responsive %>% filter(is_motility)
  other_subset <- activin_responsive %>% filter(!is_motility)

  bind_rows(
    motility_subset %>%
      count(category) %>%
      mutate(total = sum(n), pct = round(100 * n / total, 1),
             gene_set = sprintf("Motility (n=%d)", sum(n))),
    other_subset %>%
      count(category) %>%
      mutate(total = sum(n), pct = round(100 * n / total, 1),
             gene_set = sprintf("All other DEGs (n=%d)", sum(n)))
  ) %>%
    mutate(timepoint = tp_label)
}

proportions <- bind_rows(
  compute_proportions(data_60_tagged, "SB50 at 60 min"),
  compute_proportions(data_120_tagged, "SB50 at 120 min"),
  compute_proportions(data_180_tagged, "SB50 at 180 min")
)

cat("\n--- Commitment proportions ---\n")
print(as.data.frame(proportions %>% arrange(timepoint, gene_set, category)))

# ============================================================================
# 5. FISHER'S EXACT TEST: Are motility genes enriched in Blocked at 60
#    and enriched in Not blocked at 120/180?
# ============================================================================

cat("\n========== STATISTICAL TESTS ==========\n")

run_fisher <- function(data, tp_label) {
  activin_responsive <- data %>% filter(category %in% c("Blocked", "Not blocked"))

  # Contingency table: motility (yes/no) vs blocked (yes/no)
  tbl <- table(
    Motility = activin_responsive$is_motility,
    Blocked = activin_responsive$category == "Blocked"
  )

  cat(sprintf("\n%s — Contingency table (Motility × Blocked):\n", tp_label))
  print(tbl)

  if (all(dim(tbl) == c(2, 2))) {
    ft <- fisher.test(tbl)
    cat(sprintf("  Fisher's exact test: OR = %.2f, p = %.4f\n",
                ft$estimate, ft$p.value))

    # Also compute % blocked in motility vs other
    n_mot_blocked <- sum(activin_responsive$is_motility & activin_responsive$category == "Blocked")
    n_mot_total <- sum(activin_responsive$is_motility)
    n_other_blocked <- sum(!activin_responsive$is_motility & activin_responsive$category == "Blocked")
    n_other_total <- sum(!activin_responsive$is_motility)
    cat(sprintf("  Motility genes blocked: %d/%d (%.1f%%)\n",
                n_mot_blocked, n_mot_total, 100 * n_mot_blocked / n_mot_total))
    cat(sprintf("  Other genes blocked:    %d/%d (%.1f%%)\n",
                n_other_blocked, n_other_total, 100 * n_other_blocked / n_other_total))

    tibble(
      timepoint = tp_label,
      n_motility_blocked = n_mot_blocked,
      n_motility_total = n_mot_total,
      pct_motility_blocked = round(100 * n_mot_blocked / n_mot_total, 1),
      n_other_blocked = n_other_blocked,
      n_other_total = n_other_total,
      pct_other_blocked = round(100 * n_other_blocked / n_other_total, 1),
      odds_ratio = ft$estimate,
      p_value = ft$p.value
    )
  } else {
    cat("  Cannot run Fisher's test (empty cells)\n")
    NULL
  }
}

fisher_results <- bind_rows(
  run_fisher(data_60_tagged, "SB50 at 60 min"),
  run_fisher(data_120_tagged, "SB50 at 120 min"),
  run_fisher(data_180_tagged, "SB50 at 180 min")
)

if (nrow(fisher_results) > 0) {
  cat("\n--- Summary: Motility gene blocking enrichment ---\n")
  print(as.data.frame(fisher_results))
  write_csv(fisher_results, results_path("q5_motility_commitment_fisher.csv"))
}

# ============================================================================
# 6. GENE-LEVEL TRACKING: Individual motility genes across timepoints
# ============================================================================

cat("\n========== INDIVIDUAL MOTILITY GENE TRACKING ==========\n")

# Build wide table: each motility gene's category at 60, 120, 180 min
motility_activin_resp <- data_60_tagged %>%
  filter(is_motility, category %in% activin_cats) %>%
  dplyr::select(gene, cat_60 = category, lfc_activin, padj_activin, activin_direction)

motility_tracking <- motility_activin_resp %>%
  left_join(
    data_120_tagged %>%
      filter(is_motility) %>%
      dplyr::select(gene, cat_120 = category, lfc_sb50_120 = lfc_sb50),
    by = "gene"
  ) %>%
  left_join(
    data_180_tagged %>%
      filter(is_motility) %>%
      dplyr::select(gene, cat_180 = category, lfc_sb50_180 = lfc_sb50),
    by = "gene"
  ) %>%
  mutate(
    commitment_path = paste(cat_60, "→", cat_120, "→", cat_180),
    # Does this gene match the hypothesis? Blocked at 60, Not blocked at 120 or 180
    supports_hypothesis = (cat_60 == "Blocked") &
                          (cat_120 == "Not blocked" | cat_180 == "Not blocked")
  ) %>%
  arrange(desc(supports_hypothesis), commitment_path, padj_activin)

cat(sprintf("\nActivin-responsive motility genes tracked: %d\n", nrow(motility_tracking)))
cat("\nCommitment paths:\n")
path_counts <- motility_tracking %>% count(commitment_path, sort = TRUE)
print(as.data.frame(path_counts))

n_support <- sum(motility_tracking$supports_hypothesis, na.rm = TRUE)
cat(sprintf("\nGenes supporting hypothesis (Blocked@60 → Not blocked@120/180): %d / %d (%.1f%%)\n",
            n_support, nrow(motility_tracking),
            100 * n_support / nrow(motility_tracking)))

write_csv(motility_tracking, results_path("q5_motility_gene_tracking.csv"))

# ============================================================================
# 7. VISUALIZATION
# ============================================================================

cat("\n========== GENERATING PLOTS ==========\n")

cat_colors <- c("Blocked" = "#2166AC", "Not blocked" = "#B2182B",
                "Reversed" = "#762A83")

# ---- Panel A: Motility gene commitment across timepoints (bar chart) ----
# One bar per timepoint showing how motility genes are split across categories

motility_prop <- bind_rows(
  data_60_tagged  %>% filter(is_motility, category %in% activin_cats) %>%
    count(category) %>% mutate(timepoint = "SB50 at\n60 min"),
  data_120_tagged %>% filter(is_motility, category %in% activin_cats) %>%
    count(category) %>% mutate(timepoint = "SB50 at\n120 min"),
  data_180_tagged %>% filter(is_motility, category %in% activin_cats) %>%
    count(category) %>% mutate(timepoint = "SB50 at\n180 min")
) %>%
  group_by(timepoint) %>%
  mutate(total = sum(n), pct = 100 * n / total) %>%
  ungroup() %>%
  mutate(
    category = factor(category, levels = c("Reversed", "Not blocked", "Blocked")),
    timepoint = factor(timepoint, levels = c("SB50 at\n60 min", "SB50 at\n120 min", "SB50 at\n180 min"))
  )

p_bars <- ggplot(motility_prop, aes(x = timepoint, y = pct, fill = category)) +
  geom_col(position = "stack", width = 0.65) +
  geom_text(
    aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
    position = position_stack(vjust = 0.5),
    size = 3.2, color = "white", fontface = "bold", lineheight = 0.85
  ) +
  scale_fill_manual(values = cat_colors,
                    breaks = c("Blocked", "Not blocked", "Reversed")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 102)) +
  labs(
    title = "Motility gene commitment",
    x = "", y = "% of Activin-responsive\nmotility genes", fill = NULL
  ) +
  theme_minimal(base_size = 11, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    aspect.ratio = 1.4,
    plot.margin = margin(t = 5, r = 25, b = 5, l = 5)
  )

# ---- Panel B: Heatmap — ONLY delayed-commitment motility genes ----
# Genes blocked at 60min but not blocked at 120 and/or 180min

delayed_genes <- motility_tracking %>%
  filter(supports_hypothesis) %>%
  arrange(commitment_path, padj_activin)

cat(sprintf("\nDelayed-commitment motility genes (Blocked@60 → Not blocked later): %d\n",
            nrow(delayed_genes)))
cat("Genes:", paste(delayed_genes$gene, collapse = ", "), "\n")

write_csv(delayed_genes, results_path("q5_motility_delayed_commitment_genes.csv"))

if (nrow(delayed_genes) > 0) {
  heatmap_data <- delayed_genes %>%
    dplyr::select(gene, cat_60, cat_120, cat_180) %>%
    pivot_longer(cols = c(cat_60, cat_120, cat_180),
                 names_to = "timepoint", values_to = "category") %>%
    mutate(
      timepoint = recode(timepoint,
                          "cat_60" = "SB50 at 60 min",
                          "cat_120" = "SB50 at 120 min",
                          "cat_180" = "SB50 at 180 min"),
      timepoint = factor(timepoint, levels = c("SB50 at 60 min",
                                                "SB50 at 120 min",
                                                "SB50 at 180 min")),
      gene = factor(gene, levels = rev(delayed_genes$gene))
    )

  p_heatmap <- ggplot(heatmap_data, aes(x = timepoint, y = gene, fill = category)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(
      values = c(cat_colors, "NS" = "grey90", "SB50-specific" = "#1B7837"),
      na.value = "grey95", name = ""
    ) +
    labs(
      title = "Delayed-commitment genes",
      x = "", y = ""
    ) +
    theme_minimal(base_size = 10, base_family = "Helvetica") +
    theme(
      axis.text.y = element_text(face = "italic", size = 7),
      axis.text.x = element_text(size = 9),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "bottom",
      plot.margin = margin(t = 5, r = 5, b = 5, l = 15)
    )
}

# ---- Panel C: Temporal expression of delayed-commitment genes (Exp1) ----

# Normalize Exp1
exp1_meta <- metadata %>% filter(experiment == "Exp1")
exp1_counts <- counts_filtered[, rownames(exp1_meta)]
dds_exp1 <- DESeqDataSetFromMatrix(exp1_counts, exp1_meta, ~ 1)
dds_exp1 <- estimateSizeFactors(dds_exp1)
norm_exp1 <- counts(dds_exp1, normalized = TRUE)

exp1_15ngml <- exp1_meta %>% filter(concentration == "15ngml")
exp1_0ngml  <- exp1_meta %>% filter(concentration == "0ngml")
timepoints_exp1 <- sort(unique(exp1_15ngml$time_min))

compute_temporal_profile <- function(genes, norm_mat, meta_15, meta_0, tp) {
  map_dfr(genes, function(gene) {
    row_idx <- which(tolower(rownames(norm_mat)) == gene)
    if (length(row_idx) == 0) return(NULL)
    map_dfr(tp, function(t) {
      treated <- norm_mat[row_idx, rownames(meta_15)[meta_15$time_min == t]]
      control <- norm_mat[row_idx, rownames(meta_0)[meta_0$time_min == t]]
      if (length(treated) == 0 || length(control) == 0) return(NULL)
      tibble(gene = gene, time_min = t,
             log2fc = log2(mean(treated) + 1) - log2(mean(control) + 1))
    })
  })
}

delayed_profiles <- compute_temporal_profile(
  tolower(delayed_genes$gene), norm_exp1, exp1_15ngml, exp1_0ngml, timepoints_exp1
)

if (nrow(delayed_profiles) > 0) {
  # Add direction from the tracking data
  direction_map <- delayed_genes %>%
    dplyr::select(gene, activin_direction) %>%
    mutate(gene = tolower(gene))
  delayed_profiles <- delayed_profiles %>%
    left_join(direction_map, by = "gene") %>%
    mutate(direction = ifelse(activin_direction == "Up", "Upregulated", "Downregulated"))

  dir_means <- delayed_profiles %>%
    group_by(direction, time_min) %>%
    summarise(mean_fc = mean(log2fc), .groups = "drop")

  endpoints <- delayed_profiles %>% filter(time_min == max(time_min))

  n_up <- sum(delayed_genes$activin_direction == "Up")
  n_down <- sum(delayed_genes$activin_direction == "Down")

  dir_colors <- c("Upregulated" = "#D73027", "Downregulated" = "#4575B4")

  p_expression <- ggplot(delayed_profiles, aes(x = time_min, y = log2fc)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_line(aes(group = gene, color = direction), alpha = 0.3, linewidth = 0.4) +
    geom_line(data = dir_means, aes(y = mean_fc, color = direction),
              linewidth = 1.4) +
    geom_point(data = dir_means, aes(y = mean_fc, color = direction), size = 2.5) +
    geom_text_repel(data = endpoints,
                    aes(label = gene, color = direction), size = 2.3,
                    fontface = "italic", direction = "y", hjust = 0,
                    nudge_x = 10, segment.size = 0.2, max.overlaps = 40,
                    show.legend = FALSE) +
    facet_wrap(~ direction, ncol = 2, scales = "free_y") +
    scale_color_manual(values = dir_colors, guide = "none") +
    scale_x_continuous(breaks = timepoints_exp1,
                       expand = expansion(mult = c(0.05, 0.3))) +
    labs(
      x = "Time (min)", y = "log2FC (15 vs 0 ng/ml Activin)",
      title = "Temporal expression — delayed-commitment motility genes (Exp1)",
      subtitle = sprintf("Genes blocked at SB50-60min but committed later (n=%d: %d up, %d down)",
                          nrow(delayed_genes), n_up, n_down)
    ) +
    theme_minimal(base_size = 10, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5),
      strip.text = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )

  write_csv(delayed_profiles, results_path("q5_motility_delayed_expression.csv"))
}

# ---- Combine panels ----

combined <- (p_bars | p_heatmap) /
  p_expression +
  plot_layout(heights = c(1, 1), widths = c(1, 1.2)) +
  plot_annotation(
    title = "Cell motility genes: commitment to Activin programme",
    subtitle = paste0(
      "Motility genes not committed at SB50-60min (no protrusions) ",
      "become committed at 120-180min (protrusions observed)"
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30")
    )
  )

ggsave(results_path("q5_motility_commitment.pdf"), combined,
       width = 14, height = 12)
cat("Saved:", results_path("q5_motility_commitment.pdf"), "\n")

# ============================================================================
# 8. SUMMARY
# ============================================================================

cat("\n========== SUMMARY ==========\n")
cat("\nKey question: Are motility genes disproportionately 'Blocked' at SB50-60min\n")
cat("(uncommitted → no motility phenotype) but 'Not blocked' at 120/180min\n")
cat("(committed → protrusions observed)?\n\n")

if (nrow(fisher_results) > 0) {
  for (i in seq_len(nrow(fisher_results))) {
    r <- fisher_results[i, ]
    cat(sprintf("%s:\n", r$timepoint))
    cat(sprintf("  Motility genes blocked: %.1f%% (%d/%d)\n",
                r$pct_motility_blocked, r$n_motility_blocked, r$n_motility_total))
    cat(sprintf("  Other genes blocked:    %.1f%% (%d/%d)\n",
                r$pct_other_blocked, r$n_other_blocked, r$n_other_total))
    cat(sprintf("  Fisher OR=%.2f, p=%.4f %s\n\n",
                r$odds_ratio, r$p_value,
                ifelse(r$p_value < 0.05, "(*)", "")))
  }
}

cat(sprintf("\nDelayed-commitment motility genes: %d / %d (%.1f%%)\n",
            nrow(delayed_genes), nrow(motility_tracking),
            100 * nrow(delayed_genes) / nrow(motility_tracking)))
cat("These genes:", paste(delayed_genes$gene, collapse = ", "), "\n")

cat("\n========== Q5 MOTILITY COMMITMENT ANALYSIS COMPLETE ==========\n")
cat("\nOutput files:\n")
cat("  - q5_motility_commitment.pdf (combined: bars + heatmap + expression)\n")
cat("  - q5_motility_delayed_commitment_genes.csv (gene list)\n")
cat("  - q5_motility_delayed_expression.csv (temporal profiles)\n")
cat("  - q5_motility_gene_tracking.csv (all motility genes tracking)\n")
cat("  - q5_motility_commitment_fisher.csv (statistical tests)\n")
