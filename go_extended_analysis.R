# ============================================================================
# GO EXTENDED ANALYSIS
# ============================================================================
#
# Additional GO analyses requested:
#   1. Export ALL GO results as CSVs (from existing xlsx)
#   2. Comparative GO figure: Exp2 showing same terms as Exp1
#   3. GO for motility and adhesion terms (BP ontology)
#   4. GO for SB50-specific genes (green category from q3 blocking)
#   5. GO for all SB50 conditions (BP + CC) — complete listing
#

source("preprocess.R")
library(clusterProfiler)
library(org.Dr.eg.db)
library(openxlsx)

counts_filt <- counts_filtered

run_deseq_go <- function(test_meta, ref_meta) {
  run_deseq(test_meta, ref_meta, counts_mat = counts_filt)
}

run_go <- function(genes, ont = "MF", universe = NULL) {
  if (length(genes) < 3) return(NULL)
  entrez <- bitr(genes, "SYMBOL", "ENTREZID", org.Dr.eg.db, drop = TRUE)
  if (nrow(entrez) < 3) return(NULL)

  uni_entrez <- NULL
  if (!is.null(universe)) {
    uni <- bitr(universe, "SYMBOL", "ENTREZID", org.Dr.eg.db, drop = TRUE)
    uni_entrez <- uni$ENTREZID
  }

  ego <- enrichGO(entrez$ENTREZID, OrgDb = org.Dr.eg.db, ont = ont,
                  pAdjustMethod = "BH", pvalueCutoff = 0.05,
                  qvalueCutoff = 0.1, readable = TRUE,
                  universe = uni_entrez)
  if (is.null(ego) || nrow(ego) == 0) return(NULL)
  ego@result
}

# ============================================================================
# 1. EXPORT ALL GO RESULTS AS CSVs
# ============================================================================

cat("\n========== EXPORT ALL GO RESULTS AS CSVs ==========\n")

if (file.exists(results_path("exp1_go_mf_results.xlsx"))) {
  exp1_go <- read.xlsx(results_path("exp1_go_mf_results.xlsx"), sheet = 1)
  write_csv(exp1_go, results_path("go_exp1_mf_all.csv"))
  cat("Saved:", results_path("go_exp1_mf_all.csv"),
      sprintf("(%d terms)\n", nrow(exp1_go)))
} else {
  cat("  exp1_go_mf_results.xlsx not found — run go_analysis_simple.R first\n")
  exp1_go <- NULL
}

if (file.exists(results_path("exp2_go_mf_results.xlsx"))) {
  exp2_go <- read.xlsx(results_path("exp2_go_mf_results.xlsx"), sheet = 1)
  write_csv(exp2_go, results_path("go_exp2_mf_all.csv"))
  cat("Saved:", results_path("go_exp2_mf_all.csv"),
      sprintf("(%d terms)\n", nrow(exp2_go)))
} else {
  cat("  exp2_go_mf_results.xlsx not found — run go_analysis_simple.R first\n")
  exp2_go <- NULL
}

# ============================================================================
# 2. COMPARATIVE GO FIGURE: Exp2 with same terms as Exp1
# ============================================================================

cat("\n========== COMPARATIVE GO: Exp2 with Exp1 terms ==========\n")

if (!is.null(exp1_go) && !is.null(exp2_go)) {
  # Get top Exp1 terms
  exp1_top <- exp1_go %>%
    group_by(Description) %>%
    summarise(min_p = min(p.adjust), n_cond = n(), .groups = "drop") %>%
    arrange(min_p) %>%
    slice_head(n = 15) %>%
    pull(Description)

  cat("  Top 15 Exp1 GO MF terms:", paste(head(exp1_top, 5), collapse = "; "), "...\n")

  # Now run Exp2 GO for ALL SB50 conditions, keeping ONLY the Exp1 terms
  exp2 <- metadata[metadata$experiment == "Exp2", ]
  baseline <- exp2[exp2$concentration == "0ngml_DMSO", ]
  activin <- exp2[exp2$concentration == "15ngml_DMSO", ]

  comparisons <- list()
  for (time in c(60, 120, 180)) {
    sb50 <- exp2[exp2$concentration == "SB50" & exp2$time_min == time, ]
    if (nrow(sb50) < 2) next

    # vs Baseline
    name_bl <- paste0("SB50_", time, "min_vs_Baseline")
    res_bl <- run_deseq_go(sb50, baseline)
    de_bl <- res_bl$gene[res_bl$padj < 0.05 & abs(res_bl$log2FoldChange) >= 1]
    go_bl <- run_go(de_bl, "MF")
    if (!is.null(go_bl)) {
      comparisons[[name_bl]] <- go_bl %>% mutate(comparison = name_bl, time_min = time)
    }

    # vs Activin
    name_act <- paste0("SB50_", time, "min_vs_Activin")
    res_act <- run_deseq_go(sb50, activin)
    de_act <- res_act$gene[res_act$padj < 0.05 & abs(res_act$log2FoldChange) >= 1]
    go_act <- run_go(de_act, "MF")
    if (!is.null(go_act)) {
      comparisons[[name_act]] <- go_act %>% mutate(comparison = name_act, time_min = time)
    }
  }

  # Activin vs Baseline
  res_ab <- run_deseq_go(activin, baseline)
  de_ab <- res_ab$gene[res_ab$padj < 0.05 & abs(res_ab$log2FoldChange) >= 1]
  go_ab <- run_go(de_ab, "MF")
  if (!is.null(go_ab)) {
    comparisons[["Activin_vs_Baseline"]] <- go_ab %>%
      mutate(comparison = "Activin_vs_Baseline", time_min = NA)
  }

  exp2_go_full <- bind_rows(comparisons) %>% filter(p.adjust < 0.05)

  # Filter to Exp1 top terms ONLY
  exp2_filtered <- exp2_go_full %>%
    filter(Description %in% exp1_top) %>%
    mutate(
      neg_log_p = -log10(p.adjust),
      comparison = factor(comparison, levels = c(
        "Activin_vs_Baseline",
        "SB50_60min_vs_Baseline", "SB50_120min_vs_Baseline", "SB50_180min_vs_Baseline",
        "SB50_60min_vs_Activin", "SB50_120min_vs_Activin", "SB50_180min_vs_Activin"
      )),
      Description_short = str_wrap(Description, width = 40)
    )

  # Order terms same as in Exp1
  exp2_filtered$Description_short <- factor(
    exp2_filtered$Description_short,
    levels = unique(str_wrap(exp1_top, width = 40))
  )

  x_labels <- c("Activin", "60'", "120'", "180'", "60'", "120'", "180'")

  p_range <- range(exp2_filtered$neg_log_p, na.rm = TRUE)
  p_mid <- mean(p_range)

  p_comp <- ggplot(exp2_filtered, aes(x = comparison, y = Description_short)) +
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
             fill = "#B2182B", alpha = 0.12) +
    annotate("rect", xmin = 1.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
             fill = "#2166AC", alpha = 0.12) +
    annotate("rect", xmin = 4.5, xmax = 7.5, ymin = -Inf, ymax = Inf,
             fill = "#1B7837", alpha = 0.12) +
    annotate("label", x = 1, y = Inf, label = "Activin\neffect", vjust = 0,
             fontface = "bold", size = 3, fill = "#B2182B", color = "white",
             label.size = 0, label.padding = unit(0.2, "lines")) +
    annotate("label", x = 3, y = Inf, label = "SB50 vs Baseline", vjust = 0,
             fontface = "bold", size = 3, fill = "#2166AC", color = "white",
             label.size = 0, label.padding = unit(0.2, "lines")) +
    annotate("label", x = 6, y = Inf, label = "SB50 vs Activin", vjust = 0,
             fontface = "bold", size = 3, fill = "#1B7837", color = "white",
             label.size = 0, label.padding = unit(0.2, "lines")) +
    geom_point(aes(size = Count, fill = neg_log_p),
               shape = 21, color = "white", stroke = 0.6) +
    scale_fill_gradient2(
      low = "#4393C3", mid = "#F7F7F7", high = "#B2182B",
      midpoint = p_mid,
      name = expression(-log[10](p[adj]))
    ) +
    scale_size_continuous(range = c(3, 10), name = "Genes", breaks = c(5, 10, 20)) +
    guides(size = guide_legend(override.aes = list(fill = "grey40"))) +
    scale_x_discrete(labels = x_labels) +
    coord_cartesian(clip = "off") +
    labs(x = "", y = "",
         title = "Exp 2: GO Molecular Function (same terms as Exp 1)",
         subtitle = "Terms selected from top-15 Exp1 enrichments") +
    theme_minimal(base_size = 11, base_family = "Helvetica") +
    theme(
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.6),
      axis.text.y = element_text(size = 8.5, color = "grey15"),
      axis.text.x = element_text(size = 9, color = "grey30", face = "bold"),
      legend.position = "right",
      plot.margin = margin(t = 35, r = 10, b = 10, l = 5),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey50")
    )

  ggsave(results_path("go_exp2_comparative_exp1_terms.pdf"), p_comp,
         width = 11, height = 6)
  cat("Saved:", results_path("go_exp2_comparative_exp1_terms.pdf"), "\n")

  write_csv(exp2_filtered, results_path("go_exp2_comparative_terms.csv"))
}

# ============================================================================
# 3. GO FOR MOTILITY AND ADHESION (BP)
# ============================================================================

cat("\n========== GO: MOTILITY & ADHESION (BP) ==========\n")

exp2 <- metadata[metadata$experiment == "Exp2", ]
baseline <- exp2[exp2$concentration == "0ngml_DMSO", ]
activin_exp2 <- exp2[exp2$concentration == "15ngml_DMSO", ]

motility_adhesion_go <- list()
de_results_store <- list()   # store DESeq2 results for direction analysis

cat("Running BP enrichment for motility/adhesion terms...\n")

for (time in c(60, 120, 180)) {
  sb50 <- exp2[exp2$concentration == "SB50" & exp2$time_min == time, ]
  if (nrow(sb50) < 2) next

  # vs Baseline: presence of terms = motility programme still active
  name_bl <- paste0("SB50_", time, "min_vs_Baseline")
  res_bl <- run_deseq_go(sb50, baseline)
  de_results_store[[name_bl]] <- res_bl
  de_bl <- res_bl$gene[res_bl$padj < 0.05 & abs(res_bl$log2FoldChange) >= 1]
  go_bp_bl <- run_go(de_bl, "BP")
  if (!is.null(go_bp_bl)) {
    motility_adhesion_go[[name_bl]] <- go_bp_bl %>%
      mutate(comparison = name_bl, time_min = time)
  }

  # vs Activin: presence of terms = SB50 successfully blocked those pathways
  name_act <- paste0("SB50_", time, "min_vs_Activin")
  res_act <- run_deseq_go(sb50, activin_exp2)
  de_results_store[[name_act]] <- res_act
  de_act <- res_act$gene[res_act$padj < 0.05 & abs(res_act$log2FoldChange) >= 1]
  go_bp_act <- run_go(de_act, "BP")
  if (!is.null(go_bp_act)) {
    motility_adhesion_go[[name_act]] <- go_bp_act %>%
      mutate(comparison = name_act, time_min = time)
  }
}

# Activin vs Baseline
res_ab <- run_deseq_go(activin_exp2, baseline)
de_ab <- res_ab$gene[res_ab$padj < 0.05 & abs(res_ab$log2FoldChange) >= 1]
go_bp_ab <- run_go(de_ab, "BP")
if (!is.null(go_bp_ab)) {
  motility_adhesion_go[["Activin_vs_Baseline"]] <- go_bp_ab %>%
    mutate(comparison = "Activin_vs_Baseline", time_min = NA)
}

# Also run for Exp1 15ngml 240min for comparison
exp1 <- metadata[metadata$experiment == "Exp1", ]
exp1_15_240 <- exp1[exp1$concentration == "15ngml" & exp1$time_min == 240, ]
exp1_0_240 <- exp1[exp1$concentration == "0ngml" & exp1$time_min == 240, ]
if (nrow(exp1_15_240) >= 2 && nrow(exp1_0_240) >= 2) {
  res_exp1 <- run_deseq_go(exp1_15_240, exp1_0_240)
  de_exp1 <- res_exp1$gene[res_exp1$padj < 0.05 & abs(res_exp1$log2FoldChange) >= 1]
  go_bp_exp1 <- run_go(de_exp1, "BP")
  if (!is.null(go_bp_exp1)) {
    motility_adhesion_go[["Exp1_15ngml_240min"]] <- go_bp_exp1 %>%
      mutate(comparison = "Exp1_15ngml_240min", time_min = 240)
  }
}

all_bp_go <- bind_rows(motility_adhesion_go) %>% filter(p.adjust < 0.05)

# Filter for motility/adhesion/migration terms
motility_keywords <- c("motil", "migrat", "adhes", "locomot", "chemotax",
                        "cell moving", "movement")
pattern <- paste(motility_keywords, collapse = "|")

motility_terms <- all_bp_go %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))

if (nrow(motility_terms) > 0) {
  cat(sprintf("  Found %d motility/adhesion GO BP terms across conditions\n",
              n_distinct(motility_terms$Description)))

  write_csv(motility_terms, results_path("go_motility_adhesion_bp.csv"))
  cat("Saved:", results_path("go_motility_adhesion_bp.csv"), "\n")

  # Single combined plot: all 8 comparisons in one figure
  # Left group (vs Baseline): dots = motility programme still active
  # Right group (vs Activin): dots = SB50 successfully blocked motility
  motility_terms <- motility_terms %>%
    mutate(
      neg_log_p = -log10(p.adjust),
      Description_short = str_wrap(Description, width = 45),
      comparison = factor(comparison, levels = c(
        "Exp1_15ngml_240min", "Activin_vs_Baseline",
        "SB50_60min_vs_Baseline", "SB50_120min_vs_Baseline", "SB50_180min_vs_Baseline",
        "SB50_60min_vs_Activin", "SB50_120min_vs_Activin", "SB50_180min_vs_Activin"
      ))
    )

  term_order <- motility_terms %>%
    group_by(Description_short) %>%
    summarise(mean_p = mean(neg_log_p), .groups = "drop") %>%
    arrange(mean_p) %>%
    pull(Description_short)
  motility_terms$Description_short <- factor(motility_terms$Description_short, levels = term_order)

  # Background rectangles: highlight the two SB50 comparison groups
  # vs Baseline columns: positions 3-5, vs Activin columns: positions 6-8
  p_motility <- ggplot(motility_terms, aes(x = comparison, y = Description_short)) +
    # Shaded backgrounds to group comparison types
    annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
             fill = "grey90", alpha = 0.4) +
    annotate("rect", xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = Inf,
             fill = "#E8F5E9", alpha = 0.4) +
    annotate("rect", xmin = 5.5, xmax = 8.5, ymin = -Inf, ymax = Inf,
             fill = "#FFF3E0", alpha = 0.4) +
    # Group labels at top
    annotate("text", x = 1.5, y = Inf, label = "Activin effect",
             vjust = -1.5, fontface = "bold", size = 3.5, color = "grey30") +
    annotate("text", x = 4, y = Inf,
             label = "SB50 vs Baseline\n(dots = still active)",
             vjust = -0.8, fontface = "bold", size = 3, color = "#2E7D32") +
    annotate("text", x = 7, y = Inf,
             label = "SB50 vs Activin\n(dots = blocked by SB50)",
             vjust = -0.8, fontface = "bold", size = 3, color = "#E65100") +
    geom_point(aes(size = Count, fill = neg_log_p),
               shape = 21, color = "white", stroke = 0.6) +
    scale_fill_gradient(low = "#80CDC1", high = "#01665E",
                         name = expression(-log[10](p[adj]))) +
    scale_size_continuous(range = c(3, 10), name = "Genes",
                          breaks = scales::breaks_pretty(n = 4)) +
    guides(size = guide_legend(override.aes = list(fill = "grey40"))) +
    scale_x_discrete(drop = FALSE, labels = c(
      "Exp1_15ngml_240min"      = "Exp1\nActivin 240'\nvs ctrl",
      "Activin_vs_Baseline"     = "Exp2\nActivin 240'\nvs ctrl",
      "SB50_60min_vs_Baseline"  = "SB50 @ 60'\nvs ctrl",
      "SB50_120min_vs_Baseline" = "SB50 @ 120'\nvs ctrl",
      "SB50_180min_vs_Baseline" = "SB50 @ 180'\nvs ctrl",
      "SB50_60min_vs_Activin"   = "SB50 @ 60'\nvs Activin",
      "SB50_120min_vs_Activin"  = "SB50 @ 120'\nvs Activin",
      "SB50_180min_vs_Activin"  = "SB50 @ 180'\nvs Activin"
    )) +
    labs(x = NULL, y = "",
         title = "GO Biological Process: Motility & Cell Adhesion",
         subtitle = "Left: motility programme activation (vs ctrl)  |  Right: motility blocked by SB50 (vs Activin)") +
    theme_minimal(base_size = 11, base_family = "Helvetica") +
    theme(
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.6),
      axis.text.y = element_text(size = 8.5),
      axis.text.x = element_text(size = 8, face = "bold"),
      legend.position = "right",
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
      plot.margin = margin(t = 30, r = 10, b = 10, l = 10)
    )

  ggsave(results_path("go_motility_adhesion_bp.pdf"), p_motility,
         width = 13, height = 8)
  cat("Saved:", results_path("go_motility_adhesion_bp.pdf"), "\n")

  # ===========================================================================
  # 3b. DIRECTION ANALYSIS: are motility GO genes up or downregulated?
  # ===========================================================================
  # Hypothesis: genes enriched in SB50@60 vs Activin should be DOWNREGULATED
  # (lower in SB50 than in Activin → SB50 successfully blocked them).

  cat("\n========== MOTILITY GENE DIRECTION ANALYSIS ==========\n")

  # Extract unique genes from motility GO terms per comparison
  motility_gene_direction <- motility_terms %>%
    filter(grepl("vs_Activin|vs_Baseline|Activin_vs_Baseline", comparison)) %>%
    mutate(comparison = as.character(comparison)) %>%
    dplyr::select(comparison, geneID) %>%
    mutate(genes = strsplit(geneID, "/")) %>%
    tidyr::unnest(genes) %>%
    dplyr::select(comparison, gene = genes) %>%
    distinct()

  # Join with DESeq2 results to get log2FC and direction
  # Also store Activin vs Baseline results
  de_results_store[["Activin_vs_Baseline"]] <- res_ab

  motility_gene_direction <- motility_gene_direction %>%
    rowwise() %>%
    mutate(
      log2FC = {
        de <- de_results_store[[comparison]]
        if (!is.null(de) && gene %in% de$gene) {
          de$log2FoldChange[de$gene == gene]
        } else NA_real_
      },
      padj = {
        de <- de_results_store[[comparison]]
        if (!is.null(de) && gene %in% de$gene) {
          de$padj[de$gene == gene]
        } else NA_real_
      }
    ) %>%
    ungroup() %>%
    filter(!is.na(log2FC)) %>%
    mutate(
      direction = ifelse(log2FC > 0, "Upregulated", "Downregulated"),
      significant = padj < 0.05 & abs(log2FC) >= 1
    )

  # Summarise direction per comparison
  dir_summary <- motility_gene_direction %>%
    filter(significant) %>%
    group_by(comparison, direction) %>%
    summarise(n = n(), .groups = "drop") %>%
    tidyr::complete(comparison, direction, fill = list(n = 0)) %>%
    group_by(comparison) %>%
    mutate(total = sum(n), pct = n / total * 100) %>%
    ungroup()

  cat("\nDirection of significant motility GO genes per comparison:\n")
  print(dir_summary)

  # Also compute median log2FC per comparison
  median_lfc <- motility_gene_direction %>%
    filter(significant) %>%
    group_by(comparison) %>%
    summarise(median_log2FC = median(log2FC), mean_log2FC = mean(log2FC),
              n_genes = n(), .groups = "drop")
  cat("\nMedian log2FC of motility genes per comparison:\n")
  print(median_lfc)

  # ---- Plot: direction bar chart (mirrored style from go_term_direction_analysis.R) ----

  comparison_order <- c(
    "Activin_vs_Baseline",
    "SB50_60min_vs_Baseline", "SB50_120min_vs_Baseline", "SB50_180min_vs_Baseline",
    "SB50_60min_vs_Activin", "SB50_120min_vs_Activin", "SB50_180min_vs_Activin"
  )
  comparison_labels <- c(
    "Activin_vs_Baseline"      = "Activin\nvs ctrl",
    "SB50_60min_vs_Baseline"   = "SB50 @ 60'\nvs ctrl",
    "SB50_120min_vs_Baseline"  = "SB50 @ 120'\nvs ctrl",
    "SB50_180min_vs_Baseline"  = "SB50 @ 180'\nvs ctrl",
    "SB50_60min_vs_Activin"    = "SB50 @ 60'\nvs Activin",
    "SB50_120min_vs_Activin"   = "SB50 @ 120'\nvs Activin",
    "SB50_180min_vs_Activin"   = "SB50 @ 180'\nvs Activin"
  )

  dir_plot_data <- dir_summary %>%
    filter(comparison %in% comparison_order) %>%
    mutate(
      comparison = factor(comparison, levels = comparison_order),
      n_plot = ifelse(direction == "Downregulated", -n, n)
    )

  # Total gene counts per comparison for annotation
  total_counts <- dir_plot_data %>%
    group_by(comparison) %>%
    summarise(total = total[1], .groups = "drop")

  # Dynamic y-axis limits based on max gene count
  max_n <- max(abs(dir_plot_data$n_plot), na.rm = TRUE)
  y_lim <- ceiling(max_n * 1.3)

  p_direction <- ggplot(dir_plot_data, aes(x = comparison, y = n_plot, fill = direction)) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.5) +
    # Background shading
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
             fill = "#B2182B", alpha = 0.06) +
    annotate("rect", xmin = 1.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
             fill = "#2166AC", alpha = 0.06) +
    annotate("rect", xmin = 4.5, xmax = 7.5, ymin = -Inf, ymax = Inf,
             fill = "#1B7837", alpha = 0.06) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    geom_text(aes(label = abs(n)),
              size = 3.5, fontface = "bold", color = "grey20",
              vjust = ifelse(dir_plot_data$direction == "Upregulated", -0.5, 1.5)) +
    # Total gene count at top
    geom_text(data = total_counts,
              aes(x = comparison, y = y_lim * 0.95, label = paste0("n=", total)),
              inherit.aes = FALSE, size = 3, fontface = "italic", color = "grey50") +
    scale_fill_manual(
      values = c("Upregulated" = "#B2182B", "Downregulated" = "#2166AC"),
      name = ""
    ) +
    scale_x_discrete(labels = comparison_labels) +
    scale_y_continuous(limits = c(-y_lim, y_lim),
                       labels = function(x) abs(x)) +
    labs(
      x = NULL,
      y = "Number of significant motility genes",
      title = "Direction of motility GO genes across comparisons",
      subtitle = "Genes from enriched motility/migration/adhesion GO BP terms\nDownregulated in SB50 vs Activin = successfully blocked by SB50"
    ) +
    theme_minimal(base_size = 11, base_family = "Helvetica") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.6),
      axis.text.x = element_text(size = 9, face = "bold"),
      legend.position = "top",
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )

  ggsave(results_path("go_motility_direction.pdf"), p_direction,
         width = 10, height = 7)
  cat("Saved:", results_path("go_motility_direction.pdf"), "\n")

  write_csv(motility_gene_direction, results_path("go_motility_gene_direction.csv"))
  cat("Saved:", results_path("go_motility_gene_direction.csv"), "\n")

} else {
  cat("  No motility/adhesion GO BP terms found with current DEG thresholds.\n")
}

# Save ALL BP GO results as well
write_csv(all_bp_go, results_path("go_all_bp_results.csv"))
cat("Saved:", results_path("go_all_bp_results.csv"),
    sprintf("(%d terms)\n", nrow(all_bp_go)))

# ============================================================================
# 4. FULL GO LISTING: BP, CC, MF for all conditions
# ============================================================================

cat("\n========== COMPLETE GO LISTING (ALL ONTOLOGIES) ==========\n")

# Run all three ontologies for each Exp2 comparison
ontologies <- c("BP", "CC", "MF")
all_go_results <- list()

for (time in c(60, 120, 180)) {
  sb50 <- exp2[exp2$concentration == "SB50" & exp2$time_min == time, ]
  if (nrow(sb50) < 2) next

  for (ref_name in c("Baseline", "Activin")) {
    ref_meta <- if (ref_name == "Baseline") baseline else activin_exp2
    name <- paste0("SB50_", time, "min_vs_", ref_name)

    res <- run_deseq_go(sb50, ref_meta)
    de_genes <- res$gene[res$padj < 0.05 & abs(res$log2FoldChange) >= 1]

    # Split up/down
    de_up <- res$gene[res$padj < 0.05 & res$log2FoldChange >= 1]
    de_dn <- res$gene[res$padj < 0.05 & res$log2FoldChange <= -1]

    for (ont in ontologies) {
      # All DEGs
      go_all <- run_go(de_genes, ont)
      if (!is.null(go_all)) {
        all_go_results[[paste0(name, "_", ont, "_all")]] <- go_all %>%
          mutate(comparison = name, ontology = ont, direction = "all",
                 time_min = time)
      }
      # Up only
      go_up <- run_go(de_up, ont)
      if (!is.null(go_up)) {
        all_go_results[[paste0(name, "_", ont, "_up")]] <- go_up %>%
          mutate(comparison = name, ontology = ont, direction = "up",
                 time_min = time)
      }
      # Down only
      go_dn <- run_go(de_dn, ont)
      if (!is.null(go_dn)) {
        all_go_results[[paste0(name, "_", ont, "_dn")]] <- go_dn %>%
          mutate(comparison = name, ontology = ont, direction = "down",
                 time_min = time)
      }
    }
    cat(sprintf("  %s: done\n", name))
  }
}

# Activin vs Baseline
de_ab_up <- res_ab$gene[res_ab$padj < 0.05 & res_ab$log2FoldChange >= 1]
de_ab_dn <- res_ab$gene[res_ab$padj < 0.05 & res_ab$log2FoldChange <= -1]

for (ont in ontologies) {
  go_all <- run_go(de_ab, ont)
  if (!is.null(go_all)) {
    all_go_results[[paste0("Activin_vs_Baseline_", ont, "_all")]] <- go_all %>%
      mutate(comparison = "Activin_vs_Baseline", ontology = ont, direction = "all",
             time_min = NA)
  }
  go_up <- run_go(de_ab_up, ont)
  if (!is.null(go_up)) {
    all_go_results[[paste0("Activin_vs_Baseline_", ont, "_up")]] <- go_up %>%
      mutate(comparison = "Activin_vs_Baseline", ontology = ont, direction = "up",
             time_min = NA)
  }
  go_dn <- run_go(de_ab_dn, ont)
  if (!is.null(go_dn)) {
    all_go_results[[paste0("Activin_vs_Baseline_", ont, "_dn")]] <- go_dn %>%
      mutate(comparison = "Activin_vs_Baseline", ontology = ont, direction = "down",
             time_min = NA)
  }
}
cat("  Activin_vs_Baseline: done\n")

complete_go <- bind_rows(all_go_results) %>% filter(p.adjust < 0.05)
write_csv(complete_go, results_path("go_complete_all_ontologies.csv"))
cat("Saved:", results_path("go_complete_all_ontologies.csv"),
    sprintf("(%d significant terms across all conditions x ontologies)\n", nrow(complete_go)))

# Save by ontology
for (ont in ontologies) {
  ont_df <- complete_go %>% filter(ontology == ont)
  write_csv(ont_df, results_path(paste0("go_complete_", ont, ".csv")))
  cat(sprintf("  %s: %d terms → %s\n", ont, nrow(ont_df),
              results_path(paste0("go_complete_", ont, ".csv"))))
}

# ============================================================================
# 5. GO ON SB50-SPECIFIC GENES (from q3 blocking analysis)
# ============================================================================

cat("\n========== GO ON SB50-SPECIFIC GENES ==========\n")

# Read the blocking gene lists if available
sb50_specific_files <- list.files(results_path(""), pattern = "q3_gene_lists_", full.names = TRUE)

if (length(sb50_specific_files) > 0) {
  sb50_specific_go_list <- list()

  for (f in sb50_specific_files) {
    gene_df <- read_csv(f, show_col_types = FALSE)
    if ("SB50-specific" %in% unique(gene_df$category)) {
      sb50_genes <- gene_df %>% filter(category == "SB50-specific") %>% pull(gene)
      tp_name <- tools::file_path_sans_ext(basename(f)) %>%
        gsub("q3_gene_lists_", "", .)

      cat(sprintf("  %s: %d SB50-specific genes\n", tp_name, length(sb50_genes)))

      for (ont in c("BP", "MF")) {
        go_res <- run_go(sb50_genes, ont)
        if (!is.null(go_res) && sum(go_res$p.adjust < 0.05) > 0) {
          sb50_specific_go_list[[paste0(tp_name, "_", ont)]] <- go_res %>%
            filter(p.adjust < 0.05) %>%
            mutate(timepoint = tp_name, ontology = ont)
        }
      }
    }
  }

  if (length(sb50_specific_go_list) > 0) {
    sb50_specific_go <- bind_rows(sb50_specific_go_list)
    write_csv(sb50_specific_go, results_path("go_sb50_specific_genes.csv"))
    cat("Saved:", results_path("go_sb50_specific_genes.csv"),
        sprintf("(%d terms)\n", nrow(sb50_specific_go)))

    # Plot top terms
    top_sb50_terms <- sb50_specific_go %>%
      filter(ontology == "BP") %>%
      group_by(Description) %>%
      summarise(min_p = min(p.adjust), .groups = "drop") %>%
      slice_min(min_p, n = 15) %>%
      pull(Description)

    plot_sb50 <- sb50_specific_go %>%
      filter(ontology == "BP", Description %in% top_sb50_terms) %>%
      mutate(
        neg_log_p = -log10(p.adjust),
        Description_short = str_wrap(Description, width = 40)
      )

    if (nrow(plot_sb50) > 0) {
      p_sb50_go <- ggplot(plot_sb50,
                           aes(x = timepoint,
                               y = reorder(Description_short, neg_log_p))) +
        geom_point(aes(size = Count, fill = neg_log_p),
                   shape = 21, color = "white", stroke = 0.6) +
        scale_fill_gradient(low = "#B8E186", high = "#276419",
                             name = expression(-log[10](p[adj]))) +
        scale_size_continuous(range = c(3, 10), name = "Genes",
                              breaks = scales::breaks_pretty(n = 4)) +
        guides(size = guide_legend(override.aes = list(fill = "grey40"))) +
        labs(x = "SB50 timepoint", y = "",
             title = "GO BP: SB50-specific genes (not present in Activin)",
             subtitle = "Green category from blocking analysis") +
        theme_minimal(base_size = 11, base_family = "Helvetica") +
        theme(
          panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.6),
          axis.text.y = element_text(size = 8.5),
          legend.position = "right",
          plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey50")
        )

      ggsave(results_path("go_sb50_specific_bp.pdf"), p_sb50_go,
             width = 10, height = 7)
      cat("Saved:", results_path("go_sb50_specific_bp.pdf"), "\n")
    }
  } else {
    cat("  No significant GO terms for SB50-specific genes.\n")
  }
} else {
  cat("  No q3_gene_lists files found — run q3_extended_analysis.R first.\n")
}

cat("\n========== GO EXTENDED ANALYSIS COMPLETE ==========\n")
cat("\nOutput files:\n")
cat("  - go_exp1_mf_all.csv\n")
cat("  - go_exp2_mf_all.csv\n")
cat("  - go_exp2_comparative_exp1_terms.pdf + csv\n")
cat("  - go_motility_adhesion_bp.pdf + csv\n")
cat("  - go_all_bp_results.csv\n")
cat("  - go_complete_all_ontologies.csv\n")
cat("  - go_complete_BP.csv, go_complete_CC.csv, go_complete_MF.csv\n")
cat("  - go_sb50_specific_genes.csv + go_sb50_specific_bp.pdf\n")
