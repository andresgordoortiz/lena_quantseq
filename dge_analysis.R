# QuantSeq DGE Analysis
# For 3' QuantSeq data, use raw counts (salmon.merged.gene_counts.tsv) for DESeq2
# as DESeq2 handles normalization internally

source("preprocess.R")
library(pheatmap)
library(RColorBrewer)

# Create condition labels
metadata$condition <- paste0(metadata$concentration, "_", metadata$time_min, "min")

cat("Samples:", ncol(counts_filtered), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")

# ============================================================================
# 2. LOAD NODAL SCORE GENES
# ============================================================================

nodal_score_data <- read_excel("docs/nodal-score-genes_complete.xlsx", skip = 1)
nodal_genes <- tolower(na.omit(nodal_score_data$`Nodal score`))


# ============================================================================
# 3. ANALYSIS
# ============================================================================

cat("\n============================================================\n")
cat("ANALYSIS 1: 0 ng/ml vs 15 ng/ml Activin (Exp1 vs Exp2)\n")
cat("============================================================\n")

# --- Experiment 1: Activin time course ---
exp1_samples <- metadata[metadata$experiment == "Exp1" &
                           metadata$concentration %in% c("0ngml", "15ngml"), ]

# For a cleaner comparison, use 60min timepoint (common between conditions)
exp1_60min <- metadata[metadata$experiment == "Exp1" &
                         metadata$concentration %in% c("0ngml", "15ngml") &
                         metadata$time_min == 60, ]

if(nrow(exp1_60min) >= 2) {
  dds_exp1 <- DESeqDataSetFromMatrix(
    countData = counts_filtered[, rownames(exp1_60min)],
    colData = exp1_60min,
    design = ~ concentration
  )
  dds_exp1$concentration <- relevel(factor(dds_exp1$concentration), ref = "0ngml")
  dds_exp1 <- DESeq(dds_exp1)
  res_exp1 <- lfcShrink(dds_exp1, coef = "concentration_15ngml_vs_0ngml", type = "normal", quiet = TRUE)
  res_exp1_df <- as.data.frame(res_exp1) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    arrange(padj)

  cat("\nExp1 - 15ngml vs 0ngml Activin at 60min:\n")
  cat("  Total DE genes (padj < 0.05):", sum(res_exp1_df$padj < 0.05, na.rm = TRUE), "\n")
  cat("  Upregulated:", sum(res_exp1_df$padj < 0.05 & res_exp1_df$log2FoldChange > 0, na.rm = TRUE), "\n")
  cat("  Downregulated:", sum(res_exp1_df$padj < 0.05 & res_exp1_df$log2FoldChange < 0, na.rm = TRUE), "\n")
  cat("\nTop 10 DE genes (Exp1):\n")
  print(head(res_exp1_df[, c("gene", "log2FoldChange", "padj")], 10))
}

# --- Experiment 2: DMSO controls ---
exp2_dmso <- metadata[metadata$experiment == "Exp2" &
                        metadata$concentration %in% c("0ngml_DMSO", "15ngml_DMSO"), ]

if(nrow(exp2_dmso) >= 2) {
  dds_exp2 <- DESeqDataSetFromMatrix(
    countData = counts_filtered[, rownames(exp2_dmso)],
    colData = exp2_dmso,
    design = ~ concentration
  )
  dds_exp2$concentration <- relevel(factor(dds_exp2$concentration), ref = "0ngml_DMSO")
  dds_exp2 <- DESeq(dds_exp2)
  res_exp2 <- lfcShrink(dds_exp2, coef = "concentration_15ngml_DMSO_vs_0ngml_DMSO", type = "normal", quiet = TRUE)
  res_exp2_df <- as.data.frame(res_exp2) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    arrange(padj)

  cat("\nExp2 - 15ngml DMSO vs 0ngml DMSO:\n")
  cat("  Total DE genes (padj < 0.05):", sum(res_exp2_df$padj < 0.05, na.rm = TRUE), "\n")
  cat("  Upregulated:", sum(res_exp2_df$padj < 0.05 & res_exp2_df$log2FoldChange > 0, na.rm = TRUE), "\n")
  cat("  Downregulated:", sum(res_exp2_df$padj < 0.05 & res_exp2_df$log2FoldChange < 0, na.rm = TRUE), "\n")
  cat("\nTop 10 DE genes (Exp2):\n")
  print(head(res_exp2_df[, c("gene", "log2FoldChange", "padj")], 10))
}

# Compare overlap between experiments
if(exists("res_exp1_df") & exists("res_exp2_df")) {
  de_exp1 <- res_exp1_df$gene[res_exp1_df$padj < 0.05]
  de_exp2 <- res_exp2_df$gene[res_exp2_df$padj < 0.05]
  overlap <- intersect(de_exp1, de_exp2)
  cat("\nComparison Exp1 vs Exp2:\n")
  cat("  DE genes unique to Exp1:", length(setdiff(de_exp1, de_exp2)), "\n")
  cat("  DE genes unique to Exp2:", length(setdiff(de_exp2, de_exp1)), "\n")
  cat("  Shared DE genes:", length(overlap), "\n")

  # Save results
  write_csv(res_exp1_df, results_path("results_exp1_15vs0_activin.csv"))
  write_csv(res_exp2_df, results_path("results_exp2_15vs0_dmso.csv"))
}

# ============================================================================
# 4. ANALYSIS 2: Nodal Score between conditions in Experiment 2
# ============================================================================

cat("\n============================================================\n")
cat("ANALYSIS 2: Nodal Score - Experiment 2 conditions\n")
cat("============================================================\n")

# Get normalized counts for nodal genes
exp2_all <- metadata[metadata$experiment == "Exp2", ]
if(nrow(exp2_all) >= 2) {
  dds_exp2_all <- DESeqDataSetFromMatrix(
    countData = counts_filtered[, rownames(exp2_all)],
    colData = exp2_all,
    design = ~ concentration
  )
  dds_exp2_all <- estimateSizeFactors(dds_exp2_all)
  norm_counts <- counts(dds_exp2_all, normalized = TRUE)

  # Find nodal genes in our dataset
  available_nodal <- intersect(tolower(rownames(norm_counts)), nodal_genes)
  cat("Nodal score genes found in dataset:", length(available_nodal), "of", length(nodal_genes), "\n")

  # Calculate nodal score per sample (mean of z-scored nodal gene expression)
  nodal_counts <- norm_counts[rownames(norm_counts) %in% available_nodal, ]

  # Z-score normalization across samples
  nodal_zscore <- t(scale(t(nodal_counts)))
  nodal_score <- colMeans(nodal_zscore, na.rm = TRUE)

  exp2_all$nodal_score <- nodal_score[rownames(exp2_all)]

  cat("\nNodal Score by condition (Exp2):\n")
  nodal_summary <- exp2_all %>%
    group_by(concentration) %>%
    summarise(
      n = n(),
      mean_score = mean(nodal_score, na.rm = TRUE),
      sd_score = sd(nodal_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_score))
  print(nodal_summary)

  # Statistical test between conditions
  if(length(unique(exp2_all$concentration)) > 1) {
    cat("\nKruskal-Wallis test for Nodal score differences:\n")
    kw_test <- kruskal.test(nodal_score ~ concentration, data = exp2_all)
    print(kw_test)

    # Pairwise comparisons
    if(kw_test$p.value < 0.1) {
      cat("\nPairwise Wilcoxon tests (15ngml_DMSO as reference):\n")
      dmso_15_score <- exp2_all$nodal_score[exp2_all$concentration == "15ngml_DMSO"]
      for(cond in unique(exp2_all$concentration)) {
        if(cond != "15ngml_DMSO") {
          other_score <- exp2_all$nodal_score[exp2_all$concentration == cond]
          if(length(other_score) >= 2 & length(dmso_15_score) >= 2) {
            wt <- wilcox.test(other_score, dmso_15_score)
            cat(sprintf("  %s vs 15ngml_DMSO: p = %.4f\n", cond, wt$p.value))
          }
        }
      }
    }
  }

  # Heatmap of nodal genes
  if(nrow(nodal_counts) > 2) {
    annotation_col <- data.frame(
      Condition = exp2_all$concentration,
      row.names = rownames(exp2_all)
    )

    pdf(results_path("nodal_score_heatmap_exp2.pdf"), width = 12, height = 10)
    tryCatch(
      pheatmap(nodal_zscore,
               annotation_col = annotation_col,
               cluster_cols = TRUE,
               cluster_rows = TRUE,
               show_rownames = TRUE,
               show_colnames = FALSE,
               main = "Nodal Score Genes - Experiment 2",
               fontsize_row = 6),
      error = function(e) message("Heatmap error: ", e$message),
      finally = dev.off()
    )
    cat("\nHeatmap saved:", results_path("nodal_score_heatmap_exp2.pdf"), "\n")
  }

  write_csv(exp2_all, results_path("nodal_scores_exp2.csv"))
}

# ============================================================================
# 5. ANALYSIS 3: SB50 conditions vs controls
# ============================================================================

cat("\n============================================================\n")
cat("ANALYSIS 3: SB50 treatments vs controls\n")
cat("============================================================\n")

# Get all SB50 samples and their timepoints
sb50_samples <- metadata[metadata$concentration == "SB50", ]
sb50_times <- unique(sb50_samples$time_min)

cat("SB50 timepoints available:", paste(sb50_times, "min", collapse = ", "), "\n\n")

# Reference: 15ngml DMSO 4h (240min) from Exp2
dmso_15_4h <- metadata[metadata$concentration == "15ngml_DMSO" &
                         metadata$time_min == 240, ]

# If 240min not available, use available 15ngml_DMSO
if(nrow(dmso_15_4h) == 0) {
  dmso_15_ref <- metadata[metadata$concentration == "15ngml_DMSO", ]
  cat("Note: Using available 15ngml_DMSO samples as reference (240min not found)\n")
  dmso_15_4h <- dmso_15_ref
}

# Create comparison function
run_deseq_comparison <- function(test_samples, ref_samples, counts_mat, comparison_name) {
  all_samples <- rbind(
    data.frame(test_samples, group = "test"),
    data.frame(ref_samples, group = "ref")
  )
  all_samples$group <- factor(all_samples$group, levels = c("ref", "test"))

  if(nrow(all_samples) < 4) {
    cat(sprintf("  %s: Insufficient samples\n", comparison_name))
    return(NULL)
  }

  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat[, rownames(all_samples)],
    colData = all_samples,
    design = ~ group
  )
  dds <- DESeq(dds, quiet = TRUE)
  res <- lfcShrink(dds, coef = "group_test_vs_ref", type = "normal", quiet = TRUE)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    arrange(padj)

  n_de <- sum(res_df$padj < 0.05, na.rm = TRUE)
  n_up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)

  cat(sprintf("  %s: %d DE genes (↑%d, ↓%d)\n", comparison_name, n_de, n_up, n_down))

  return(res_df)
}

# Run comparisons for each SB50 timepoint
sb50_results <- list()

for(time_pt in sb50_times) {
  cat(sprintf("\nSB50 %dmin comparisons:\n", time_pt))

  sb50_time <- sb50_samples[sb50_samples$time_min == time_pt, ]

  # Comparison 1: SB50 vs 15ngml DMSO (Exp2)
  comp_name <- sprintf("SB50_%dmin_vs_15ngml_DMSO", time_pt)
  res1 <- run_deseq_comparison(sb50_time, dmso_15_4h, counts_filtered, comp_name)
  if(!is.null(res1)) {
    sb50_results[[comp_name]] <- res1
    write_csv(res1, results_path(sprintf("results_%s.csv", comp_name)))
  }

  # Comparison 2: SB50 vs 15ngml Activin at same timepoint (Exp1)
  exp1_15ngml_time <- metadata[metadata$experiment == "Exp1" &
                                 metadata$concentration == "15ngml" &
                                 metadata$time_min == time_pt, ]

  if(nrow(exp1_15ngml_time) > 0) {
    comp_name2 <- sprintf("SB50_%dmin_vs_15ngml_%dmin_Exp1", time_pt, time_pt)
    res2 <- run_deseq_comparison(sb50_time, exp1_15ngml_time, counts_filtered, comp_name2)
    if(!is.null(res2)) {
      sb50_results[[comp_name2]] <- res2
      write_csv(res2, results_path(sprintf("results_%s.csv", comp_name2)))
    }
  } else {
    cat(sprintf("  No 15ngml Activin at %dmin in Exp1\n", time_pt))
  }
}

# ============================================================================
# 6. SUMMARY VISUALIZATIONS
# ============================================================================

cat("\n============================================================\n")
cat("GENERATING SUMMARY PLOTS\n")
cat("============================================================\n")

# PCA of all samples
dds_all <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ experiment
)
vsd <- vst(dds_all, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = c("experiment", "concentration"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pdf(results_path("pca_all_samples.pdf"), width = 10, height = 8)
tryCatch(
  print(
    ggplot(pca_data, aes(x = PC1, y = PC2, color = concentration, shape = experiment)) +
      geom_point(size = 3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      theme_bw() +
      ggtitle("PCA - All Samples") +
      theme(legend.position = "right")
  ),
  error = function(e) message("Plot error: ", e$message),
  finally = dev.off()
)
cat("Saved:", results_path("pca_all_samples.pdf"), "\n")

# Volcano plots for key comparisons
create_volcano <- function(res_df, title, filename) {
  if(is.null(res_df)) return(NULL)

  res_df$significance <- case_when(
    res_df$padj < 0.05 & res_df$log2FoldChange > 1 ~ "Up",
    res_df$padj < 0.05 & res_df$log2FoldChange < -1 ~ "Down",
    TRUE ~ "NS"
  )

  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    theme_bw() +
    ggtitle(title) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40")

  ggsave(filename, p, width = 8, height = 6)
  cat(sprintf("Saved: %s\n", filename))
}

if(exists("res_exp1_df")) {
  create_volcano(res_exp1_df, "Exp1: 15ngml vs 0ngml Activin (60min)", results_path("volcano_exp1.pdf"))
}
if(exists("res_exp2_df")) {
  create_volcano(res_exp2_df, "Exp2: 15ngml DMSO vs 0ngml DMSO", results_path("volcano_exp2.pdf"))
}

# Nodal score barplot
if(exists("exp2_all") && "nodal_score" %in% colnames(exp2_all)) {
  pdf(results_path("nodal_score_barplot.pdf"), width = 8, height = 6)
  tryCatch({
    nodal_plot_data <- exp2_all %>%
      group_by(concentration) %>%
      summarise(
        mean = mean(nodal_score, na.rm = TRUE),
        se = sd(nodal_score, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )

    p <- ggplot(nodal_plot_data, aes(x = reorder(concentration, -mean), y = mean)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
      theme_bw() +
      labs(x = "Condition", y = "Nodal Score (mean z-score)",
           title = "Nodal Score by Condition (Exp2)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
  },
  error = function(e) message("Plot error: ", e$message),
  finally = dev.off()
  )
  cat("Saved:", results_path("nodal_score_barplot.pdf"), "\n")
}

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("\nOutput files generated (all in results/):\n")
cat("  - results_exp1_15vs0_activin.csv\n")
cat("  - results_exp2_15vs0_dmso.csv\n")
cat("  - nodal_scores_exp2.csv\n")
cat("  - results_SB50_*_vs_*.csv (multiple comparisons)\n")
cat("  - pca_all_samples.pdf\n")
cat("  - volcano_exp1.pdf, volcano_exp2.pdf\n")
cat("  - nodal_score_heatmap_exp2.pdf\n")
cat("  - nodal_score_barplot.pdf\n")
