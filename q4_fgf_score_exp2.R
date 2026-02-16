# Question 4: FGF Signaling Score - Cumulative Expression with Variance Ribbons
# Direct FGF targets, FGF ligands, DUSP family, and key mesoderm genes
# Publication-ready figure for LaTeX (Helvetica)
# Inspired by the Nodal Score (q2_nodal_score_exp2.R)

source("preprocess.R")
library(patchwork)

# ============================================================================
# 1. LOAD AND PREPARE DATA (via preprocess.R — outlier S343239 already excluded)
# ============================================================================

metadata_q4 <- metadata %>%
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

counts_exp2 <- counts_filtered[, rownames(metadata_q4)]
counts_q4 <- counts_exp2[rowSums(counts_exp2 >= 10) >= 3, ]

metadata_q4$condition <- factor(metadata_q4$condition)
dds <- DESeqDataSetFromMatrix(counts_q4, metadata_q4, ~ 1)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

all_genes <- tolower(rownames(norm_counts))

# Use metadata_q4 as 'metadata' locally for the scoring function
metadata <- metadata_q4

# ============================================================================
# 2. DEFINE GENE SETS
# ============================================================================

# --- FGF direct target genes (well-established in zebrafish literature) ---
# ETS transcription factors (Raible & Brand 2001; Roehl & Nüsslein-Volhard 2001)
# Sprouty negative feedback (Fürthauer et al. 2001; Hacohen et al. 1998)
# DUSP6/MKP3 (Echevarria et al. 2005; Tsang et al. 2004)
# Sef/IL17RD (Tsang et al. 2002; Furthauer et al. 2002)
# Spred (Sivak et al. 2005)
fgf_direct_targets <- c(
  # ETS family - canonical FGF readout
  "etv4",         # pea3 - direct FGF target

  "etv5a",        # erm - direct FGF target
  "etv5b",        # erm-like
  # Sprouty family - negative feedback regulators
  "spry1",
  "spry2",
  "spry4",        # most responsive to FGF
  # DUSP6 (MKP3) - direct transcriptional target of FGF/ERK
  "dusp6",
  # Sef - specific FGF inhibitor
  "il17rd",       # sef in zebrafish
  # Spred family
  "spred1",
  "spred2a",
  "spred2b",
  "spred3",
  # FLRT - fibronectin leucine rich repeat
  "flrt3"
)

# --- FGF ligands family ---
fgf_ligands <- c(
  "fgf3",  "fgf4",  "fgf5",  "fgf6a", "fgf6b",
  "fgf7",  "fgf8a", "fgf8b", "fgf10a", "fgf10b",
  "fgf16", "fgf17", "fgf18a", "fgf18b", "fgf19",
  "fgf1a", "fgf1b", "fgf2",
  "fgf20a", "fgf20b", "fgf21", "fgf22", "fgf23", "fgf24"
)

# --- DUSP (dual-specificity phosphatase) family ---
dusp_family <- c(
  "dusp1", "dusp2", "dusp3a", "dusp3b", "dusp4", "dusp5", "dusp6",
  "dusp7", "dusp8a", "dusp8b", "dusp10", "dusp11", "dusp12",
  "dusp13a", "dusp14", "dusp16", "dusp19a", "dusp19b",
  "dusp22a", "dusp22b", "dusp23a", "dusp23b", "dusp26", "dusp27",
  "dusp28", "dusp29"
)

# --- Key mesoderm / FGF-responsive genes ---
mesoderm_fgf_genes <- c(
  "tbx16",        # spadetail - FGF responsive mesoderm TF
  "tbx6",         # related T-box factor
  "tbxta",        # no tail / brachyury
  "noto",         # notochord homeobox
  "msgn1",        # mesogenin - paraxial mesoderm
  "her1",         # hairy/enhancer of split - somitogenesis clock
  "her7",         # somitogenesis clock gene
  "hes6",         # hairy/enhancer of split family
  "mespaa",       # mesoderm posterior
  "mespab",
  "mespba",
  "mespbb",
  "ripply1",      # somite boundary
  "ripply2"
)

# ============================================================================
# 3. SETTINGS
# ============================================================================

conditions <- c("0ngml_DMSO", "15ngml_DMSO", "SB50_60min", "SB50_120min", "SB50_180min")
cond_labels <- c("0 ng/ml (Control)", "15 ng/ml Activin", "SB50 60 min", "SB50 120 min", "SB50 180 min")
names(cond_labels) <- conditions
cond_colors <- c("0ngml_DMSO" = "#BDBDBD", "15ngml_DMSO" = "#5FB358",
                 "SB50_60min" = "#FDBF6F", "SB50_120min" = "#FF7F00", "SB50_180min" = "#E31A1C")

# ============================================================================
# 4. HELPER FUNCTION: Compute cumulative score + stats + permutation tests
# ============================================================================

compute_cumulative_score <- function(gene_list, norm_counts, metadata, conditions,
                                     score_name = "Score", cond_labels, cond_colors) {

  # Match genes (case-insensitive)
  available <- gene_list[tolower(gene_list) %in% tolower(rownames(norm_counts))]
  gene_counts <- norm_counts[match(tolower(available), tolower(rownames(norm_counts))), , drop = FALSE]
  rownames(gene_counts) <- available
  n_genes <- nrow(gene_counts)

  cat(sprintf("\n=== %s ===\n", score_name))
  cat(sprintf("  Requested: %d | Found: %d | Missing: %s\n",
              length(gene_list), n_genes,
              paste(setdiff(gene_list, available), collapse = ", ")))
  cat(sprintf("  Genes used: %s\n", paste(available, collapse = ", ")))

  if (n_genes < 2) {
    cat("  Skipping: fewer than 2 genes found.\n")
    return(NULL)
  }

  # Order genes by control expression
  ctrl_cols <- metadata$condition == "0ngml_DMSO"
  gene_order <- order(rowMeans(gene_counts[, ctrl_cols, drop = FALSE]))
  ordered_genes <- rownames(gene_counts)[gene_order]

  # Cumulative expression per sample
  cumsum_mat <- sapply(colnames(gene_counts), function(s) {
    cumsum(log2(gene_counts[ordered_genes, s] + 1))
  })

  # Stats per condition at each rank
  cumsum_stats <- expand.grid(gene_rank = 1:n_genes, condition = conditions) %>%
    rowwise() %>%
    mutate(
      vals = list(cumsum_mat[gene_rank, rownames(metadata)[metadata$condition == condition]]),
      mean = mean(unlist(vals)),
      sd = sd(unlist(vals))
    ) %>%
    dplyr::select(-vals) %>%
    ungroup() %>%
    mutate(condition = factor(condition, levels = conditions))

  # Per-sample final cumulative score
  final_cumsum <- data.frame(
    sample = colnames(cumsum_mat),
    condition = factor(metadata[colnames(cumsum_mat), "condition"], levels = conditions),
    cumsum_score = cumsum_mat[n_genes, ]
  )

  # Permutation tests
  set.seed(42)
  n_perm <- 10000

  calc_mean_curve <- function(idx, mat) rowMeans(mat[, idx, drop = FALSE])

  ctrl_idx <- which(metadata$condition == "0ngml_DMSO")
  ctrl_curve <- calc_mean_curve(ctrl_idx, cumsum_mat)

  pairwise <- map_dfr(conditions[-1], function(cond) {
    test_idx <- which(metadata$condition == cond)
    test_curve <- calc_mean_curve(test_idx, cumsum_mat)
    obs_diff <- mean(test_curve - ctrl_curve)
    obs_max_dev <- max(abs(test_curve - ctrl_curve))

    all_idx <- c(ctrl_idx, test_idx)
    n_ctrl <- length(ctrl_idx)
    n_test <- length(test_idx)

    perm_diffs <- replicate(n_perm, {
      perm <- sample(all_idx)
      mean(calc_mean_curve(perm[(n_ctrl+1):(n_ctrl+n_test)], cumsum_mat) -
             calc_mean_curve(perm[1:n_ctrl], cumsum_mat))
    })

    p_perm <- (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1)
    sig <- case_when(p_perm < 0.001 ~ "***", p_perm < 0.01 ~ "**", p_perm < 0.05 ~ "*", TRUE ~ "ns")

    cat(sprintf("  %s: mean diff = %+.2f, p = %.4f %s\n", cond_labels[cond], obs_diff, p_perm, sig))

    tibble(condition = cond, mean_diff = obs_diff, max_dev = obs_max_dev, p_perm = p_perm, sig = sig)
  })

  # Build final stats for annotation
  final_stats <- cumsum_stats %>% filter(gene_rank == n_genes)
  ctrl_final <- final_stats$mean[final_stats$condition == "0ngml_DMSO"]

  annot <- pairwise %>%
    filter(sig != "ns") %>%
    arrange(desc(abs(mean_diff))) %>%
    mutate(
      y_test = final_stats$mean[match(condition, final_stats$condition)],
      bar_x = n_genes * 1.02
    )

  # Plot
  p <- ggplot(cumsum_stats, aes(x = gene_rank, y = mean, color = condition, fill = condition)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = cond_colors, labels = cond_labels, name = "") +
    scale_fill_manual(values = cond_colors, labels = cond_labels, name = "") +
    labs(
      x = bquote(.(score_name) * " genes (ranked, " * italic(n) * " = " * .(n_genes) * ")"),
      y = expression(Sigma ~ log[2](counts + 1)),
      caption = "Shading: \u00b11 SD | Permutation test (10,000 iter) vs control"
    ) +
    coord_cartesian(xlim = c(1, n_genes * 1.18)) +
    theme_bw(base_size = 11, base_family = "Helvetica") +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      plot.caption = element_text(size = 7, hjust = 0, color = "#666666"),
      plot.margin = margin(10, 40, 10, 10)
    )

  # Add significance brackets
  bar_width <- 0.4
  for (i in seq_len(nrow(annot))) {
    x_offset <- (i - 1) * 1.2
    bar_x <- n_genes + 0.5 + x_offset
    bar_color <- cond_colors[annot$condition[i]]

    p <- p +
      annotate("segment", x = bar_x, xend = bar_x,
               y = ctrl_final, yend = annot$y_test[i],
               color = bar_color, linewidth = 0.8) +
      annotate("segment", x = bar_x - bar_width, xend = bar_x,
               y = ctrl_final, yend = ctrl_final,
               color = bar_color, linewidth = 0.8) +
      annotate("segment", x = bar_x - bar_width, xend = bar_x,
               y = annot$y_test[i], yend = annot$y_test[i],
               color = bar_color, linewidth = 0.8) +
      annotate("text", x = bar_x + 0.3, y = (ctrl_final + annot$y_test[i]) / 2,
               label = annot$sig[i], size = 4, hjust = 0, fontface = "bold", color = "black")
  }

  list(plot = p, stats = pairwise, final_scores = final_cumsum,
       genes_used = available, n_genes = n_genes, score_name = score_name)
}

# ============================================================================
# 5. RUN ALL GENE SET SCORES
# ============================================================================

results <- list()

# A) FGF Direct Targets
results$fgf_targets <- compute_cumulative_score(
  fgf_direct_targets, norm_counts, metadata, conditions,
  "FGF Direct Targets", cond_labels, cond_colors
)

# B) FGF Ligands
results$fgf_ligands <- compute_cumulative_score(
  fgf_ligands, norm_counts, metadata, conditions,
  "FGF Ligands", cond_labels, cond_colors
)

# C) DUSP Family
results$dusp_family <- compute_cumulative_score(
  dusp_family, norm_counts, metadata, conditions,
  "DUSP Family", cond_labels, cond_colors
)

# D) Mesoderm / FGF-responsive genes (tbx16, noto, msgn1, etc.)
results$mesoderm <- compute_cumulative_score(
  mesoderm_fgf_genes, norm_counts, metadata, conditions,
  "Mesoderm / FGF-responsive", cond_labels, cond_colors
)

# ============================================================================
# 6. INDIVIDUAL PLOTS
# ============================================================================

for (name in names(results)) {
  res <- results[[name]]
  if (is.null(res)) next

  fname_base <- paste0("q4_", name, "_score")
  p_titled <- res$plot +
    plot_annotation(
      title = paste0(res$score_name, " Score after SB50 inhibition"),
      subtitle = "Cumulative expression of gene set (ranked by control)",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "Helvetica"),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#666666", family = "Helvetica")
      )
    )

  pdf(results_path(paste0(fname_base, ".pdf")), width = 10, height = 7, family = "Helvetica")
  tryCatch(print(p_titled),
           error = function(e) message("Plot error: ", e$message),
           finally = dev.off())

  png(results_path(paste0(fname_base, ".png")), width = 10, height = 5, units = "in", res = 300)
  tryCatch(print(p_titled),
           error = function(e) message("Plot error: ", e$message),
           finally = dev.off())

  write_csv(res$stats, results_path(paste0(fname_base, "_stats.csv")))
  write_csv(res$final_scores, results_path(paste0(fname_base, "_per_sample.csv")))

  cat(sprintf("Saved: %s.pdf/png\n", results_path(fname_base)))
}

# ============================================================================
# 7. COMBINED PANEL FIGURE
# ============================================================================

valid_results <- Filter(Negate(is.null), results)

if (length(valid_results) >= 2) {
  panel_plots <- imap(valid_results, function(res, name) {
    res$plot + ggtitle(res$score_name) +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "none"
      )
  })

  # Combine into panel
  combined <- wrap_plots(panel_plots, ncol = 2) +
    plot_annotation(
      title = "Gene Set Scores after SB50 inhibition",
      subtitle = "Cumulative log2(counts+1) ranked by control expression",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Helvetica"),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "#666666", family = "Helvetica")
      )
    )

  # Extract one legend
  legend_plot <- valid_results[[1]]$plot +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 10))
  legend_grob <- cowplot::get_legend(legend_plot)

  final_panel <- combined / legend_grob + plot_layout(heights = c(20, 1))

  pdf(results_path("q4_fgf_scores_combined.pdf"), width = 14, height = 12, family = "Helvetica")
  tryCatch(print(final_panel),
           error = function(e) message("Plot error: ", e$message),
           finally = dev.off())

  png(results_path("q4_fgf_scores_combined.png"), width = 14, height = 12, units = "in", res = 300)
  tryCatch(print(final_panel),
           error = function(e) message("Plot error: ", e$message),
           finally = dev.off())

  cat("\nSaved:", results_path("q4_fgf_scores_combined.pdf/png"), "\n")
}

# ============================================================================
# 8. SUMMARY TABLE
# ============================================================================

summary_table <- map_dfr(valid_results, function(res) {
  res$stats %>% mutate(gene_set = res$score_name)
})
write_csv(summary_table, results_path("q4_fgf_all_stats.csv"))

cat("\n============================================================\n")
cat("FGF SCORE ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("\nGene sets analyzed:\n")
for (res in valid_results) {
  cat(sprintf("  %s: %d genes\n", res$score_name, res$n_genes))
}
cat("\nOutput files (all in results/):\n")
cat("  - q4_fgf_targets_score.pdf/png + stats\n")
cat("  - q4_fgf_ligands_score.pdf/png + stats\n")
cat("  - q4_dusp_family_score.pdf/png + stats\n")
cat("  - q4_mesoderm_score.pdf/png + stats\n")
cat("  - q4_fgf_scores_combined.pdf/png (panel figure)\n")
cat("  - q4_fgf_all_stats.csv\n")
