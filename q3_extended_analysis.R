# ============================================================================
# Q3 EXTENDED ANALYSIS
# ============================================================================
#
# Additional analyses requested:
#   1. Blocking figure WITHOUT SB50-specific (green) category
#   2. Gene lists (CSVs) for each category at each timepoint
#   3. Gene transfer tracking: which genes change category across timepoints
#   4. Expression patterns of genes that ARE blocked by SB50
#   5. Blocking analysis restricted to shared DEGs (Exp1 ∩ Exp2)
#   6. GO enrichment on SB50-specific genes
#

source("preprocess.R")
library(patchwork)
library(ggrepel)
library(clusterProfiler)
library(org.Dr.eg.db)

# ============================================================================
# RECOMPUTE BLOCKING DATA (from q3_sb50_comparison.R logic)
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

cat_colors <- c("Blocked" = "#2166AC", "Not blocked" = "#B2182B",
                "Reversed" = "#762A83", "SB50-specific" = "#1B7837", "NS" = "grey85")

# ============================================================================
# 1. GENE LISTS (CSVs) FOR EACH CATEGORY AT EACH TIMEPOINT
# ============================================================================

cat("\n========== GENE LISTS PER CATEGORY ==========\n")

save_category_lists <- function(data, tp) {
  for (cat_name in c("Blocked", "Not blocked", "Reversed", "SB50-specific")) {
    genes_df <- data %>%
      filter(category == cat_name) %>%
      dplyr::select(gene, lfc_activin, padj_activin, lfc_sb50, padj_sb50, activin_direction) %>%
      arrange(padj_activin)
    fname <- results_path(sprintf("q3_genelist_%s_%dmin.csv",
                                   gsub("[- ]", "_", tolower(cat_name)), tp))
    write_csv(genes_df, fname)
    cat(sprintf("  %s @ %dmin: %d genes → %s\n", cat_name, tp, nrow(genes_df), fname))
  }
}

save_category_lists(data_60, 60)
save_category_lists(data_120, 120)
save_category_lists(data_180, 180)

# ============================================================================
# 2. BLOCKING FIGURE WITHOUT SB50-SPECIFIC CATEGORY
# ============================================================================

cat("\n========== BLOCKING FIGURE (NO SB50-SPECIFIC) ==========\n")

make_scatter_clean <- function(data, timepoint) {
  sb50_duration <- 240 - timepoint
  data_clean <- data %>% filter(category != "SB50-specific")
  ggplot(data_clean, aes(x = lfc_activin, y = lfc_sb50, color = category)) +
    geom_point(data = filter(data_clean, category == "NS"), alpha = 0.1, size = 0.5) +
    geom_point(data = filter(data_clean, category != "NS"), alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted", linewidth = 0.6) +
    scale_color_manual(
      values = cat_colors,
      breaks = c("Blocked", "Not blocked", "Reversed")
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

p_clean_60 <- make_scatter_clean(data_60, 60)
p_clean_120 <- make_scatter_clean(data_120, 120)
p_clean_180 <- make_scatter_clean(data_180, 180)

# Stacked bar without SB50-specific
all_data_clean <- bind_rows(data_60, data_120, data_180) %>%
  filter(category %in% c("Blocked", "Not blocked", "Reversed")) %>%
  mutate(timepoint = factor(timepoint, levels = c(60, 120, 180)))

prop_clean <- all_data_clean %>%
  group_by(timepoint, category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(timepoint) %>%
  mutate(total = sum(n), pct = 100 * n / total) %>%
  ungroup() %>%
  mutate(category = factor(category, levels = c("Reversed", "Not blocked", "Blocked")))

p_bar_clean <- ggplot(prop_clean, aes(x = timepoint, y = pct, fill = category)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(
    aes(label = ifelse(pct > 5, sprintf("%d\n(%.0f%%)", n, pct), "")),
    position = position_stack(vjust = 0.5),
    size = 2.5, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(values = cat_colors, breaks = c("Blocked", "Not blocked", "Reversed")) +
  scale_x_discrete(labels = c("60" = "60min", "120" = "120min", "180" = "180min")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  labs(title = "Category proportions (Activin-responsive only)",
       x = "SB50 added at", y = "% of Activin DE genes", fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    panel.grid.major.x = element_blank(),
    aspect.ratio = 1
  ) +
  guides(fill = guide_legend(nrow = 1))

combined_clean <- (p_clean_60 + p_clean_120) / (p_clean_180 + p_bar_clean) +
  plot_annotation(
    title = "SB50 blocking of Activin-induced genes (SB50-specific removed)",
    subtitle = "Only genes significantly DE by Activin vs baseline shown",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
    )
  )

ggsave(results_path("q3_blocking_no_sb50specific.pdf"), combined_clean, width = 8, height = 8)
cat("Saved:", results_path("q3_blocking_no_sb50specific.pdf"), "\n")

# ============================================================================
# 3. GENE TRANSFER ANALYSIS: ONTOLOGY OF GENES SWITCHING CATEGORIES
# ============================================================================
#
# Track genes that switch between "Blocked" and "Not blocked" across
# timepoints (60 → 120 → 180 min), then identify what biological functions
# (GO) are enriched in each transition group.
#

cat("\n========== GENE TRANSFER ANALYSIS ==========\n")

# Build wide table: one row per gene, columns = category at each timepoint
transfer <- data_60 %>%
  dplyr::select(gene, cat_60 = category, lfc_activin, activin_direction) %>%
  left_join(
    data_120 %>% dplyr::select(gene, cat_120 = category),
    by = "gene"
  ) %>%
  left_join(
    data_180 %>% dplyr::select(gene, cat_180 = category),
    by = "gene"
  ) %>%
  mutate(
    full_path = paste(cat_60, "→", cat_120, "→", cat_180)
  )

# Summary of full paths
cat("\nFull transition paths (60 → 120 → 180 min):\n")
path_summary <- transfer %>%
  count(full_path, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 1))
print(as.data.frame(path_summary))

# Save full tracking
write_csv(transfer, results_path("q3_gene_transfer_tracking.csv"))
write_csv(path_summary, results_path("q3_transfer_path_summary.csv"))

# ---- Define key transition groups ----
# Focus on genes that are Activin-responsive (significant in at least one tp)
transfer_sig <- transfer %>%
  filter(cat_60 != "NS" | cat_120 != "NS" | cat_180 != "NS")

# Key groups for GO analysis
# "Blocked" = gene's activin response is reversible (needs longer exposure to commit)
# "Not blocked" = gene already committed to activin-induced expression
transition_groups <- list(
  # Blocked at all 3 cutoffs → gene never commits regardless of exposure
  "Always reversible" = transfer_sig %>%
    filter(cat_60 == "Blocked" & cat_120 == "Blocked" & cat_180 == "Blocked") %>%
    pull(gene),
  # Not blocked at all 3 cutoffs → gene committed within 60 min of activin
  "Early commitment (<60 min)" = transfer_sig %>%
    filter(cat_60 == "Not blocked" & cat_120 == "Not blocked" & cat_180 == "Not blocked") %>%
    pull(gene),
  # Blocked at 60 but not at 180 → gene commits between 60-180 min of activin
  # These genes required prolonged activin to lock in their expression change
  "Delayed commitment (60-180 min)" = transfer_sig %>%
    filter(cat_60 == "Blocked" & cat_180 == "Not blocked") %>%
    pull(gene),
  # Not blocked at 60 but blocked at 180 → paradoxical: committed early but
  # then expression returns toward baseline. Reflects transient activin response:
  # gene peaks early and naturally declines, so SB50 at 60 catches it high
  # (looks committed) but SB50 at 180 finds it already declining (looks blocked)
  "Transient response" = transfer_sig %>%
    filter(cat_60 == "Not blocked" & cat_180 == "Blocked") %>%
    pull(gene)
)

# Also capture the 60→120 and 120→180 transitions
transition_groups[["Blocked→Activated (60→120)"]] <- transfer_sig %>%
  filter(cat_60 == "Blocked" & cat_120 == "Not blocked") %>% pull(gene)
transition_groups[["Blocked→Activated (120→180)"]] <- transfer_sig %>%
  filter(cat_120 == "Blocked" & cat_180 == "Not blocked") %>% pull(gene)
transition_groups[["Activated→Blocked (60→120)"]] <- transfer_sig %>%
  filter(cat_60 == "Not blocked" & cat_120 == "Blocked") %>% pull(gene)
transition_groups[["Activated→Blocked (120→180)"]] <- transfer_sig %>%
  filter(cat_120 == "Not blocked" & cat_180 == "Blocked") %>% pull(gene)

cat("\n--- Transition group sizes ---\n")
for (nm in names(transition_groups)) {
  cat(sprintf("  %s: %d genes\n", nm, length(transition_groups[[nm]])))
}

# Save gene lists for each group
for (nm in names(transition_groups)[1:4]) {
  genes <- transition_groups[[nm]]
  if (length(genes) > 0) {
    df <- transfer_sig %>%
      filter(gene %in% genes) %>%
      dplyr::select(gene, lfc_activin, activin_direction, cat_60, cat_120, cat_180) %>%
      arrange(desc(abs(lfc_activin)))
    fname <- gsub("[→ ()]", "_", tolower(nm))
    fname <- gsub("_+", "_", fname)
    fname <- gsub("_$", "", fname)
    write_csv(df, results_path(paste0("q3_transfer_", fname, ".csv")))
  }
}

# ---- Panel A: Gene commitment bar chart ----
# Show how gene commitment to activin program changes with exposure time.
# SB50 severs activin signaling at administration time; all collected at 240 min.
# "Blocked" = gene's activin effect is reversible after X min → not yet committed
# "Not blocked" = gene committed to activin-induced state within X min
#
# Key transitions as activin exposure increases:
#   Blocked → Blocked       = Still uncommitted (needs even more activin)
#   Blocked → Not blocked   = Newly committed (commits with additional exposure)
#   Not blocked → Blocked   = Commitment reverted (paradoxical, dynamic regulation)
#   Not blocked → Not blocked = Already committed (locked in early)

count_transitions <- function(data, from_col, to_col, interval_label) {
  data %>%
    filter(!!sym(from_col) %in% c("Blocked", "Not blocked", "Reversed"),
           !!sym(to_col) %in% c("Blocked", "Not blocked", "Reversed")) %>%
    count(from = !!sym(from_col), to = !!sym(to_col)) %>%
    mutate(interval = interval_label)
}

trans_60_120 <- count_transitions(transfer_sig, "cat_60", "cat_120", "60 → 120 min")
trans_120_180 <- count_transitions(transfer_sig, "cat_120", "cat_180", "120 → 180 min")

# Focus on the 4 biologically meaningful Blocked <-> Not blocked transitions
key_trans <- bind_rows(trans_60_120, trans_120_180) %>%
  filter(from %in% c("Blocked", "Not blocked"),
         to %in% c("Blocked", "Not blocked")) %>%
  mutate(
    bio_label = case_when(
      from == "Blocked" & to == "Blocked"       ~ "Still\nuncommitted",
      from == "Blocked" & to == "Not blocked"   ~ "Newly\ncommitted",
      from == "Not blocked" & to == "Blocked"   ~ "Commitment\nreverted",
      from == "Not blocked" & to == "Not blocked" ~ "Already\ncommitted"
    ),
    bio_label = factor(bio_label, levels = c(
      "Still\nuncommitted", "Newly\ncommitted",
      "Already\ncommitted", "Commitment\nreverted"
    )),
    interval = factor(interval, levels = c("60 → 120 min", "120 → 180 min"))
  )

# Colors: blue = reversible (uncommitted), red = committed
transition_colors <- c(
  "Still\nuncommitted"   = "#2166AC",
  "Newly\ncommitted"     = "#EF8A62",
  "Already\ncommitted"   = "#B2182B",
  "Commitment\nreverted" = "#67A9CF"
)

p_transit <- ggplot(key_trans, aes(x = bio_label, y = n, fill = bio_label)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5, fontface = "bold") +
  facet_wrap(~ interval) +
  scale_fill_manual(values = transition_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = NULL, y = "Number of genes",
    title = "Gene commitment to activin program by exposure time",
    subtitle = "Blue: activin effect still reversible | Red: gene committed to activin-induced expression"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save full transition matrix for reference (including Reversed)
all_trans <- bind_rows(trans_60_120, trans_120_180)
write_csv(all_trans, results_path("q3_transition_matrix.csv"))

# ---- Panel B: GO enrichment on transition groups (matching barplot) ----
# Run GO on the EXACT same gene sets shown in each bar above,
# so gene counts in facet labels match the bar values directly.
bg_map_transfer <- bitr(rownames(counts_filtered), fromType = "SYMBOL",
                         toType = "ENTREZID", OrgDb = org.Dr.eg.db)

# Build transition gene lists that match the barplot bars exactly
bar_go_groups <- list()
for (interval_info in list(
  list(from_col = "cat_60", to_col = "cat_120", interval = "60\u2192120 min"),
  list(from_col = "cat_120", to_col = "cat_180", interval = "120\u2192180 min")
)) {
  for (trans in list(
    list(from = "Blocked", to = "Blocked", bio = "Still uncommitted"),
    list(from = "Blocked", to = "Not blocked", bio = "Newly committed"),
    list(from = "Not blocked", to = "Not blocked", bio = "Already committed"),
    list(from = "Not blocked", to = "Blocked", bio = "Commitment reverted")
  )) {
    genes <- transfer_sig %>%
      filter(!!sym(interval_info$from_col) == trans$from,
             !!sym(interval_info$to_col) == trans$to) %>%
      pull(gene)
    if (length(genes) > 0) {
      label <- paste0(trans$bio, " (", interval_info$interval, ")")
      bar_go_groups[[label]] <- genes
    }
  }
}

cat("\n--- Transition bar groups for GO ---\n")
for (nm in names(bar_go_groups)) {
  cat(sprintf("  %s: %d genes\n", nm, length(bar_go_groups[[nm]])))
}

# Run GO BP and MF on each bar group
go_results_list <- list()
for (nm in names(bar_go_groups)) {
  genes <- bar_go_groups[[nm]]
  if (length(genes) < 5) {
    cat(sprintf("  Skipping GO for '%s' (%d genes)\n", nm, length(genes)))
    next
  }

  gene_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb = org.Dr.eg.db)
  if (nrow(gene_ids) < 3) next

  for (ont in c("BP", "MF")) {
    ego <- enrichGO(
      gene = gene_ids$ENTREZID,
      universe = bg_map_transfer$ENTREZID,
      OrgDb = org.Dr.eg.db,
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      go_df <- as.data.frame(ego) %>%
        head(6) %>%
        mutate(group = nm, ontology = ont)
      go_results_list[[paste(nm, ont)]] <- go_df
      cat(sprintf("  GO %s for '%s': %d terms (showing top 6)\n",
                  ont, nm, nrow(as.data.frame(ego))))
    }
  }
}

# Combine GO results
if (length(go_results_list) > 0) {
  go_combined <- bind_rows(go_results_list) %>%
    mutate(
      GeneRatioNum = sapply(GeneRatio, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]]); parts[1] / parts[2]
      }),
      neg_log_p = -log10(p.adjust),
      Description = str_wrap(Description, width = 35),
      ontology = factor(ontology, levels = c("BP", "MF"),
                        labels = c("Biological Process", "Molecular Function"))
    )

  # Build facet labels with gene counts matching barplot bars
  # group names are "bio_label (interval)" — add (n = X genes) suffix
  group_n_labels <- setNames(
    sapply(names(bar_go_groups), function(nm) {
      paste0(nm, "  [", length(bar_go_groups[[nm]]), " genes]")
    }),
    names(bar_go_groups)
  )
  # Keep only groups that appear in GO results
  used_groups <- intersect(names(bar_go_groups),
                           unique(go_combined$group))
  go_combined$group <- factor(
    group_n_labels[as.character(go_combined$group)],
    levels = group_n_labels[used_groups]
  )

  # Order terms within each group×ontology facet by gene ratio
  go_combined <- go_combined %>%
    group_by(group, ontology) %>%
    mutate(Description = factor(Description,
                                levels = Description[order(GeneRatioNum)])) %>%
    ungroup()

  n_groups <- length(used_groups)
  go_height <- max(14, 3 + n_groups * 2.5)  # dynamic height

  p_go_transfer <- ggplot(go_combined,
                           aes(x = GeneRatioNum, y = Description,
                               size = Count, fill = neg_log_p)) +
    geom_point(shape = 21, color = "white", stroke = 0.6) +
    scale_fill_gradient(low = "#80CDC1", high = "#01665E",
                        name = expression(-log[10](p[adj]))) +
    scale_size_continuous(range = c(3, 8), name = "Genes",
                          breaks = scales::breaks_pretty(n = 4)) +
    guides(size = guide_legend(override.aes = list(fill = "grey40"))) +
    facet_grid(group ~ ontology, scales = "free_y", space = "free_y") +
    labs(x = "Gene Ratio", y = NULL,
         title = "GO enrichment by commitment transition",
         subtitle = "Biological processes and molecular functions for each pairwise transition") +
    theme_minimal(base_size = 10, base_family = "Helvetica") +
    theme(
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.6),
      axis.text.y = element_text(size = 7.5, color = "grey15"),
      strip.text.y = element_text(size = 8, face = "bold", angle = 0),
      strip.text.x = element_text(size = 10, face = "bold"),
      legend.position = "right",
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey50")
    )

  # Combined figure: commitment chart + GO
  combined_transfer <- p_transit / p_go_transfer +
    plot_layout(heights = c(1, 4)) +
    plot_annotation(
      title = "Activin commitment dynamics and associated biological functions",
      subtitle = "Gene commitment transitions between consecutive SB50 time-points (60, 120, 180 min)",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
      )
    )

  ggsave(results_path("q3_gene_transfer_go.pdf"), combined_transfer,
         width = 13, height = go_height)
  cat("Saved:", results_path("q3_gene_transfer_go.pdf"), "\n")

  # Save GO table
  write_csv(go_combined %>% dplyr::select(group, ontology, ID, Description,
                                           GeneRatio, pvalue, p.adjust, Count, geneID),
            results_path("q3_transfer_go_enrichment.csv"))
} else {
  # No GO results — save transition plot alone
  ggsave(results_path("q3_gene_transfer_go.pdf"), p_transit, width = 8, height = 6)
  cat("Saved:", results_path("q3_gene_transfer_go.pdf"), " (no GO results)\n")
}

# Also save the simpler transition-only plot
ggsave(results_path("q3_gene_transfer_plot.pdf"), p_transit, width = 8, height = 5)
cat("Saved:", results_path("q3_gene_transfer_plot.pdf"), "\n")

# ============================================================================
# 4. NOT-BLOCKED GENES: HIGHLIGHT EXPRESSION DIFFERENCES
# ============================================================================
#
# Genes in the "Not blocked" (red) category — Activin changes them and
# SB50 fails to block them. Show their expression patterns.
#

cat("\n========== NOT-BLOCKED GENE PATTERNS ==========\n")

# Normalize Exp2 for heatmap
exp2_meta <- metadata %>% filter(experiment == "Exp2")
exp2_counts <- counts_filtered[, rownames(exp2_meta)]
dds_exp2 <- DESeqDataSetFromMatrix(exp2_counts, exp2_meta, ~ 1)
dds_exp2 <- estimateSizeFactors(dds_exp2)
norm_exp2 <- counts(dds_exp2, normalized = TRUE)

# Get not-blocked genes at each timepoint
not_blocked_genes <- list(
  "60"  = data_60  %>% filter(category == "Not blocked") %>% pull(gene),
  "120" = data_120 %>% filter(category == "Not blocked") %>% pull(gene),
  "180" = data_180 %>% filter(category == "Not blocked") %>% pull(gene)
)

# Genes that are always not blocked — most interesting
always_nb <- Reduce(intersect, not_blocked_genes)
cat(sprintf("Always not blocked (all 3 timepoints): %d genes\n", length(always_nb)))
write_csv(
  data.frame(gene = always_nb) %>%
    left_join(data_60 %>% dplyr::select(gene, lfc_activin, padj_activin), by = "gene") %>%
    arrange(padj_activin),
  results_path("q3_always_not_blocked_genes.csv")
)

# ============================================================================
# 5. BLOCKED GENES: EXPRESSION PATTERNS
# ============================================================================

cat("\n========== BLOCKED GENE EXPRESSION ==========\n")

# Genes blocked at 60min (most blocked)
blocked_60 <- data_60 %>%
  filter(category == "Blocked") %>%
  arrange(padj_activin) %>%
  pull(gene)

cat(sprintf("Blocked at 60min: %d genes\n", length(blocked_60)))

# Show expression of blocked genes across conditions
blocked_expr <- map_dfr(blocked_60, function(g) {
  row_idx <- which(rownames(norm_exp2) == g)
  if (length(row_idx) == 0) return(NULL)
  data.frame(
    gene = g,
    sample = colnames(norm_exp2),
    norm_count = norm_exp2[row_idx, ],
    condition = exp2_meta$concentration,
    time_min = exp2_meta$time_min
  )
}) %>%
  mutate(
    cond_label = case_when(
      condition == "0ngml_DMSO" ~ "Baseline",
      condition == "15ngml_DMSO" ~ "Activin",
      condition == "SB50" & time_min == 60 ~ "SB50 60min",
      condition == "SB50" & time_min == 120 ~ "SB50 120min",
      condition == "SB50" & time_min == 180 ~ "SB50 180min"
    ),
    cond_label = factor(cond_label, levels = c("Baseline", "Activin",
                                                "SB50 60min", "SB50 120min", "SB50 180min"))
  ) %>%
  filter(!is.na(cond_label))

# Aggregate by condition
blocked_agg <- blocked_expr %>%
  group_by(gene, cond_label) %>%
  summarise(mean_expr = mean(log2(norm_count + 1)), .groups = "drop") %>%
  group_by(gene) %>%
  mutate(z_expr = {
    s <- sd(mean_expr)
    if (is.na(s) || s == 0) 0 else (mean_expr - mean(mean_expr)) / s
  }) %>%
  ungroup()

# Select top 30 blocked genes by significance for heatmap
top_blocked <- head(blocked_60, 30)

blocked_heatmap_data <- blocked_agg %>%
  filter(gene %in% top_blocked) %>%
  pivot_wider(id_cols = gene, names_from = cond_label, values_from = z_expr) %>%
  column_to_rownames("gene")

# Drop rows/cols with all NAs, replace remaining NAs with 0
blocked_heatmap_data <- blocked_heatmap_data[rowSums(!is.na(blocked_heatmap_data)) > 0, , drop = FALSE]
blocked_heatmap_data[is.na(blocked_heatmap_data)] <- 0

if (nrow(blocked_heatmap_data) > 3) {
  library(pheatmap)

  cond_order <- c("Baseline", "Activin", "SB50 60min", "SB50 120min", "SB50 180min")
  # Only keep columns that exist in the data
  cond_order <- intersect(cond_order, colnames(blocked_heatmap_data))
  blocked_heatmap_data <- blocked_heatmap_data[, cond_order, drop = FALSE]

  heatmap_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

  pdf(results_path("q3_blocked_genes_expression.pdf"), width = 6, height = 10,
      family = "Helvetica")
  tryCatch(
    pheatmap(
      as.matrix(blocked_heatmap_data),
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      color = heatmap_colors,
      fontsize = 9,
      fontsize_row = 8,
      main = "Blocked genes: expression across conditions (z-score)",
      border_color = NA,
      cellwidth = 30,
      cellheight = 10
    ),
    error = function(e) message("Heatmap error: ", e$message),
    finally = dev.off()
  )
  cat("Saved:", results_path("q3_blocked_genes_expression.pdf"), "\n")
} else {
  cat("Skipped heatmap: too few blocked genes (", nrow(blocked_heatmap_data), ")\n")
}

# ============================================================================
# 6. BLOCKING WITH SHARED DEGs ONLY (Exp1 ∩ Exp2)
# ============================================================================

cat("\n========== BLOCKING WITH SHARED DEGs ONLY ==========\n")

# Load shared DEGs from q1
shared_degs <- read_csv(results_path("q1_shared_de_genes.csv"), show_col_types = FALSE)
shared_genes <- shared_degs$gene
cat(sprintf("Shared DEGs (Exp1 ∩ Exp2): %d\n", length(shared_genes)))

get_blocking_shared <- function(res_activin, res_sb50, tp, shared) {
  inner_join(
    res_activin %>% dplyr::select(gene, lfc_activin = log2FoldChange, padj_activin = padj),
    res_sb50 %>% dplyr::select(gene, lfc_sb50 = log2FoldChange, padj_sb50 = padj),
    by = "gene"
  ) %>%
    filter(gene %in% shared) %>%
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
      timepoint = tp
    )
}

data_60_shared <- get_blocking_shared(res_activin_vs_baseline, res_sb50_60_vs_baseline, 60, shared_genes)
data_120_shared <- get_blocking_shared(res_activin_vs_baseline, res_sb50_120_vs_baseline, 120, shared_genes)
data_180_shared <- get_blocking_shared(res_activin_vs_baseline, res_sb50_180_vs_baseline, 180, shared_genes)

# Scatter plots
make_scatter_shared <- function(data, timepoint) {
  sb50_duration <- 240 - timepoint
  ggplot(data, aes(x = lfc_activin, y = lfc_sb50, color = category)) +
    geom_point(data = filter(data, category == "NS"), alpha = 0.15, size = 0.8) +
    geom_point(data = filter(data, category != "NS"), alpha = 0.7, size = 1.8) +
    geom_hline(yintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted", linewidth = 0.6) +
    scale_color_manual(values = cat_colors,
                       breaks = c("Blocked", "Not blocked", "Reversed", "SB50-specific")) +
    labs(
      title = sprintf("SB50 at %dmin (shared DEGs)", timepoint),
      subtitle = sprintf("(%dmin SB50 exposure)", sb50_duration),
      x = "Activin effect (log2FC)", y = "SB50+Activin effect (log2FC)",
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

p_sh60 <- make_scatter_shared(data_60_shared, 60)
p_sh120 <- make_scatter_shared(data_120_shared, 120)
p_sh180 <- make_scatter_shared(data_180_shared, 180)

# Bar
shared_all <- bind_rows(data_60_shared, data_120_shared, data_180_shared) %>%
  filter(category != "NS") %>%
  mutate(timepoint = factor(timepoint))

shared_prop <- shared_all %>%
  count(timepoint, category) %>%
  group_by(timepoint) %>%
  mutate(total = sum(n), pct = 100 * n / total) %>%
  ungroup() %>%
  mutate(category = factor(category, levels = c("SB50-specific", "Reversed", "Not blocked", "Blocked")))

p_bar_shared <- ggplot(shared_prop, aes(x = timepoint, y = pct, fill = category)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = ifelse(pct > 5, sprintf("%d\n(%.0f%%)", n, pct), "")),
            position = position_stack(vjust = 0.5),
            size = 2.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = cat_colors,
                    breaks = c("Blocked", "Not blocked", "Reversed", "SB50-specific")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  labs(title = "Shared DEGs only", x = "SB50 added at", y = "% DE genes", fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    aspect.ratio = 1
  )

combined_shared <- (p_sh60 + p_sh120) / (p_sh180 + p_bar_shared) +
  plot_annotation(
    title = "SB50 blocking — shared DEGs only (Exp1 ∩ Exp2)",
    subtitle = "Restricted to genes significantly DE in both experiments",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
    )
  )

ggsave(results_path("q3_blocking_shared_degs.pdf"), combined_shared, width = 8, height = 8)
cat("Saved:", results_path("q3_blocking_shared_degs.pdf"), "\n")

# Also save shared gene lists
for (tp in c(60, 120, 180)) {
  d <- get(paste0("data_", tp, "_shared"))
  for (cat_name in c("Blocked", "Not blocked")) {
    genes_df <- d %>% filter(category == cat_name) %>%
      dplyr::select(gene, lfc_activin, padj_activin, lfc_sb50, padj_sb50) %>%
      arrange(padj_activin)
    write_csv(genes_df, results_path(sprintf("q3_shared_%s_%dmin.csv",
                                              gsub(" ", "_", tolower(cat_name)), tp)))
  }
}

# ============================================================================
# 7. GO ENRICHMENT ON SB50-SPECIFIC GENES
# ============================================================================
#
# The SB50-specific (green) genes are those NOT altered by Activin alone
# but significantly changed when SB50 is added. What biological functions
# do they represent?
#
# Run all three GO ontologies (BP, MF, CC) per timepoint to find the most
# informative categories. The per-timepoint approach avoids the dilution
# problem of pooling: a pathway enriched at 60 min can lose significance
# when 120+180 min genes (from different pathways) are added to the pool.
#

cat("\n========== GO ON SB50-SPECIFIC GENES ==========\n")

bg_map <- bitr(rownames(counts_filtered), fromType = "SYMBOL",
                toType = "ENTREZID", OrgDb = org.Dr.eg.db)

# Per-timepoint SB50-specific genes
sb50_spec_by_tp <- list(
  "60 min"  = data_60  %>% filter(category == "SB50-specific") %>% pull(gene),
  "120 min" = data_120 %>% filter(category == "SB50-specific") %>% pull(gene),
  "180 min" = data_180 %>% filter(category == "SB50-specific") %>% pull(gene)
)

# Pooled across all timepoints (for CSV reference only — not plotted)
sb50_spec_pooled <- bind_rows(data_60, data_120, data_180) %>%
  filter(category == "SB50-specific") %>%
  distinct(gene) %>%
  pull(gene)

cat(sprintf("SB50-specific genes: 60min=%d, 120min=%d, 180min=%d, pooled=%d\n",
            length(sb50_spec_by_tp[["60 min"]]),
            length(sb50_spec_by_tp[["120 min"]]),
            length(sb50_spec_by_tp[["180 min"]]),
            length(sb50_spec_pooled)))

# Run GO for all ontologies × all timepoints
sb50_go_all <- list()
ont_labels <- c("BP" = "Biological Process", "MF" = "Molecular Function", "CC" = "Cellular Component")

for (ont in c("BP", "MF", "CC")) {
  for (tp_name in c("60 min", "120 min", "180 min")) {
    genes <- sb50_spec_by_tp[[tp_name]]
    if (length(genes) < 5) {
      cat(sprintf("  %s %s: too few genes (%d), skipping\n", ont, tp_name, length(genes)))
      next
    }

    gene_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = org.Dr.eg.db)
    if (nrow(gene_ids) < 3) next

    ego <- enrichGO(
      gene = gene_ids$ENTREZID,
      universe = bg_map$ENTREZID,
      OrgDb = org.Dr.eg.db,
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      ego_df <- as.data.frame(ego)
      key <- paste(ont, tp_name)
      sb50_go_all[[key]] <- ego_df %>%
        head(8) %>%
        mutate(timepoint = tp_name, ontology = ont)
      cat(sprintf("  %s %s: %d significant terms (showing top 8)\n",
                  ont, tp_name, sum(ego_df$p.adjust < 0.05)))

      # Save per-timepoint per-ontology CSV
      fname <- gsub(" ", "_", tolower(tp_name))
      write_csv(ego_df, results_path(sprintf("q3_sb50_specific_go_%s_%s.csv",
                                              tolower(ont), fname)))
    } else {
      cat(sprintf("  %s %s: no significant GO terms\n", ont, tp_name))
    }
  }
}

# Also save pooled results as CSVs (all ontologies)
if (length(sb50_spec_pooled) >= 5) {
  gene_ids_pooled <- bitr(sb50_spec_pooled, fromType = "SYMBOL",
                           toType = "ENTREZID", OrgDb = org.Dr.eg.db)
  for (ont in c("BP", "MF", "CC")) {
    ego <- enrichGO(
      gene = gene_ids_pooled$ENTREZID,
      universe = bg_map$ENTREZID,
      OrgDb = org.Dr.eg.db,
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    ego_df <- as.data.frame(ego)
    write_csv(ego_df, results_path(sprintf("q3_sb50_specific_go_%s_pooled.csv", tolower(ont))))
    cat(sprintf("  GO %s (pooled): %d significant terms\n", ont, sum(ego_df$p.adjust < 0.05)))
  }
}

# Produce per-timepoint GO dot plot (faceted by timepoint, each with own terms)
# Since terms rarely overlap between timepoints, this avoids sparse dot matrices
if (length(sb50_go_all) > 0) {
  sb50_go_combined <- bind_rows(sb50_go_all) %>%
    mutate(
      GeneRatioNum = sapply(GeneRatio, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]]); parts[1] / parts[2]
      }),
      neg_log_p = -log10(p.adjust),
      Description_short = str_wrap(Description, width = 40),
      timepoint = factor(timepoint, levels = c("60 min", "120 min", "180 min")),
      ontology = factor(ontology, levels = c("BP", "MF", "CC"),
                        labels = c("Biological Process", "Molecular Function",
                                   "Cellular Component"))
    )

  # Keep top 8 per timepoint×ontology, order by p-value within each facet
  sb50_go_plot <- sb50_go_combined %>%
    group_by(timepoint, ontology) %>%
    arrange(p.adjust) %>%
    slice_head(n = 8) %>%
    mutate(Description_short = factor(Description_short,
                                      levels = Description_short[order(GeneRatioNum)])) %>%
    ungroup()

  # Determine best faceting: timepoint × ontology grid
  # Only include ontologies with results
  ontologies_present <- sb50_go_plot %>% distinct(ontology) %>% pull(ontology)
  n_ont <- length(ontologies_present)
  n_tp <- length(unique(sb50_go_plot$timepoint))

  p_sb50_go <- ggplot(sb50_go_plot,
                       aes(x = GeneRatioNum, y = Description_short,
                           size = Count, fill = neg_log_p)) +
    geom_point(shape = 21, color = "white", stroke = 0.6) +
    scale_fill_gradient(low = "#C7E9C0", high = "#00441B",
                        name = expression(-log[10](p[adj]))) +
    scale_size_continuous(range = c(3, 10), name = "Genes",
                          breaks = scales::breaks_pretty(n = 4)) +
    guides(size = guide_legend(override.aes = list(fill = "grey40"))) +
    facet_grid(timepoint ~ ontology, scales = "free_y", space = "free_y") +
    labs(x = "Gene Ratio", y = NULL,
         title = "GO enrichment: SB50-specific genes by activin exposure time",
         subtitle = "Activin-independent SB50 effects; each panel shows top enriched terms for that timepoint") +
    theme_minimal(base_size = 11, base_family = "Helvetica") +
    theme(
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.6),
      axis.text.y = element_text(size = 8, color = "grey15"),
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y = element_text(size = 10, face = "bold", angle = 0),
      legend.position = "right",
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey50")
    )

  # Height: ~3 per timepoint-row + 2 for legend/titles
  plot_height <- 2 + n_tp * 3.5

  ggsave(results_path("q3_sb50_specific_go.pdf"), p_sb50_go,
         width = 12, height = plot_height)
  cat("Saved:", results_path("q3_sb50_specific_go.pdf"), "\n")

  write_csv(sb50_go_combined, results_path("q3_sb50_specific_go_combined.csv"))
} else {
  cat("No GO results for SB50-specific genes — skipping plot\n")
}

cat("\n========== Q3 EXTENDED ANALYSIS COMPLETE ==========\n")
