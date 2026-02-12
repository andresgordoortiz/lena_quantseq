# ============================================================================
# SHARED PREPROCESSING MODULE
# ============================================================================
#
# Source this file at the top of every analysis script:
#   source("preprocess.R")
#
# PROVIDES:
#   - counts_filtered  : integer count matrix, genes with ≥10 counts in ≥3 samples
#   - metadata         : data.frame with sample, treatment, experiment, concentration, time_min
#   - RESULTS_DIR      : "results" (all outputs go here)
#   - results_path()   : helper to build output paths inside results/
#   - run_deseq()      : DESeq2 + lfcShrink (normal prior) wrapper
#   - print_summary()  : quick summary printer
#
# CHANGES FROM PREVIOUS SCRIPTS:
#   1. lfcShrink with normal prior shrinks extreme log2FC values (fixes log2FC > 10)
#   2. All outputs go to results/ via results_path()
#   3. Consistent metadata parsing (q3-style: SB50 first, DMSO, then Activin)
#

library(DESeq2)
library(readr)
library(tidyverse)
library(readxl)
library(ggplot2)
library(patchwork)

# ============================================================================
# OUTPUT DIRECTORY
# ============================================================================

RESULTS_DIR <- "results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

results_path <- function(...) file.path(RESULTS_DIR, ...)

# ============================================================================
# LOAD COUNTS
# ============================================================================

counts_raw <- read.table("salmon.merged.gene_counts.tsv", header = TRUE, row.names = 1)
counts_raw <- counts_raw[, -1]
counts_int <- round(counts_raw)

# ============================================================================
# BUILD METADATA (consistent across all scripts)
# ============================================================================

samples <- read_csv("samples.csv", show_col_types = FALSE)

metadata <- data.frame(
  sample = paste0("S", samples$requests_sample_sample_id),
  treatment = sub("^\\d{8}_R\\d+_", "", samples$sample_description),
  row.names = paste0("S", samples$requests_sample_sample_id)
)

metadata$experiment <- ifelse(grepl("DMSO|SB50", metadata$treatment), "Exp2", "Exp1")

# Classification: most specific patterns first (SB50 → DMSO → plain Activin)
metadata$concentration <- case_when(
  grepl("^50uMSB50", metadata$treatment)   ~ "SB50",
  grepl("DMSO", metadata$treatment) & grepl("^15ngml", metadata$treatment) ~ "15ngml_DMSO",
  grepl("DMSO", metadata$treatment) & grepl("^0ngml", metadata$treatment)  ~ "0ngml_DMSO",
  grepl("^15ngmlActivin", metadata$treatment) & !grepl("DMSO", metadata$treatment) ~ "15ngml",
  grepl("^10ngmlActivin", metadata$treatment) ~ "10ngml",
  grepl("^5ngmlActivin", metadata$treatment)  ~ "5ngml",
  grepl("^0ngmlActivin", metadata$treatment) & !grepl("DMSO", metadata$treatment) ~ "0ngml",
  TRUE ~ "other"
)

metadata$time_min <- as.numeric(str_extract(metadata$treatment, "\\d+(?=min$)"))

# ============================================================================
# FILTER COUNTS
# ============================================================================

common_samples <- intersect(colnames(counts_int), rownames(metadata))
counts_int <- counts_int[, common_samples]
metadata <- metadata[common_samples, ]

keep <- rowSums(counts_int >= 10) >= 3
counts_filtered <- counts_int[keep, ]

cat("Preprocessing loaded:\n")
cat("  Samples:", ncol(counts_filtered), "\n")
cat("  Genes after filtering:", nrow(counts_filtered), "\n")
cat("  Output directory:", RESULTS_DIR, "\n\n")

# ============================================================================
# DESeq2 HELPER WITH LFC SHRINKAGE
# ============================================================================
#
# lfcShrink with type = "normal" applies a zero-centred Normal prior to
# log2 fold-change estimates. This:
#   - Shrinks extreme LFC values toward zero (fixes log2FC > 10 artefacts)
#   - Stabilises variance for low-count genes
#   - Produces more reliable effect-size rankings
#   - Does NOT affect p-values or padj (those are from the original DESeq2 test)
#
# Note: "apeglm" was tested but its adaptive Cauchy prior actually *increases*
# estimates for strongly significant genes (pitx2: 11.8 → 13.8). The normal
# prior instead shrinks them effectively (pitx2: 11.8 → 7.5).
#

run_deseq <- function(test_samples, ref_samples, counts_mat = counts_filtered,
                      shrink = TRUE) {
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

  if (shrink) {
    # Normal prior shrinkage — always pulls extreme LFC toward zero
    res <- lfcShrink(dds, coef = "group_test_vs_ref", type = "normal", quiet = TRUE)
  } else {
    res <- results(dds, contrast = c("group", "test", "ref"))
  }

  as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    arrange(padj)
}

print_summary <- function(res_df, name, lfc_thresh = 1.0) {
  sig <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= lfc_thresh
  n_de <- sum(sig, na.rm = TRUE)
  n_up <- sum(sig & res_df$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(sig & res_df$log2FoldChange < 0, na.rm = TRUE)
  cat(sprintf("%s: %d DE (↑%d ↓%d) [padj<0.05, |log2FC|≥%.1f]\n",
              name, n_de, n_up, n_down, lfc_thresh))
  data.frame(comparison = name, de_genes = n_de, up = n_up, down = n_down)
}
