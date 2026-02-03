# GO Molecular Function Analysis - Simplified
# Danio rerio reference genome

library(DESeq2)
library(readr)
library(tidyverse)
library(clusterProfiler)
library(org.Dr.eg.db)
library(openxlsx)
library(ggplot2)

# ============================================================================
# DATA PREP
# ============================================================================

counts_raw <- read.table("salmon.merged.gene_counts.tsv", header = TRUE, row.names = 1)[, -1]
counts_int <- round(counts_raw)
samples <- read_csv("samples.csv", show_col_types = FALSE)

metadata <- data.frame(
  sample = paste0("S", samples$requests_sample_sample_id),
  treatment = sub("^\\d{8}_R\\d+_", "", samples$sample_description),
  row.names = paste0("S", samples$requests_sample_sample_id)
) %>%
  mutate(
    experiment = ifelse(grepl("DMSO|SB50", treatment), "Exp2", "Exp1"),
    concentration = case_when(
      grepl("^0ngmlActivin", treatment) & !grepl("DMSO", treatment) ~ "0ngml",
      grepl("^5ngmlActivin", treatment) ~ "5ngml",
      grepl("^10ngmlActivin", treatment) ~ "10ngml",
      grepl("^15ngmlActivin", treatment) & !grepl("DMSO", treatment) ~ "15ngml",
      grepl("^0ngml.*DMSO", treatment) ~ "0ngml_DMSO",
      grepl("^15ngml.*DMSO", treatment) ~ "15ngml_DMSO",
      grepl("^50uMSB50", treatment) ~ "SB50"
    ),
    time_min = as.numeric(str_extract(treatment, "\\d+(?=min$)"))
  )

common <- intersect(colnames(counts_int), rownames(metadata))
counts_int <- counts_int[, common]
metadata <- metadata[common, ]
counts_filt <- counts_int[rowSums(counts_int >= 10) >= 3, ]

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

run_deseq <- function(test_meta, ref_meta) {
  combined <- rbind(
    data.frame(test_meta, group = "test"),
    data.frame(ref_meta, group = "ref")
  )
  combined$group <- factor(combined$group, levels = c("ref", "test"))
  dds <- DESeqDataSetFromMatrix(counts_filt[, rownames(combined)], combined, ~ group)
  dds <- DESeq(dds, quiet = TRUE)
  results(dds, contrast = c("group", "test", "ref")) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    arrange(padj)
}

run_go_mf <- function(genes, universe = NULL) {
  if (length(genes) < 3) return(NULL)
  entrez <- bitr(genes, "SYMBOL", "ENTREZID", org.Dr.eg.db, drop = TRUE)
  if (nrow(entrez) < 3) return(NULL)

  ego <- enrichGO(entrez$ENTREZID, OrgDb = org.Dr.eg.db, ont = "MF",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1,
                  readable = TRUE)
  if (is.null(ego) || nrow(ego) == 0) return(NULL)
  ego@result
}

# ============================================================================
# EXPERIMENT 1: Activin dose-response
# ============================================================================

cat("=== EXPERIMENT 1: GO MF Analysis ===\n")
exp1 <- metadata[metadata$experiment == "Exp1", ]

exp1_results <- list()
exp1_go <- list()

for (conc in c("5ngml", "10ngml", "15ngml")) {
  for (time in c(60, 120, 180, 240)) {
    test <- exp1[exp1$concentration == conc & exp1$time_min == time, ]
    ref <- exp1[exp1$concentration == "0ngml" & exp1$time_min == time, ]

    if (nrow(test) < 2 || nrow(ref) < 2) next

    name <- paste0(conc, "_", time, "min")
    cat(sprintf("  %s: ", name))

    res <- run_deseq(test, ref)
    de_genes <- res$gene[res$padj < 0.05 & abs(res$log2FoldChange) >= 1]
    cat(sprintf("%d DE genes -> ", length(de_genes)))

    go_res <- run_go_mf(de_genes)
    if (!is.null(go_res) && nrow(go_res) > 0) {
      go_res$comparison <- name
      go_res$concentration <- conc
      go_res$time_min <- time
      exp1_go[[name]] <- go_res
      cat(sprintf("%d GO terms\n", sum(go_res$p.adjust < 0.05)))
    } else {
      cat("0 GO terms\n")
    }
    exp1_results[[name]] <- res
  }
}

# Combine and save Exp1 GO results
exp1_go_df <- bind_rows(exp1_go) %>%
  filter(p.adjust < 0.05) %>%
  select(comparison, concentration, time_min, ID, Description,
         GeneRatio, pvalue, p.adjust, Count, geneID)

# ============================================================================
# EXPERIMENT 2: SB50 inhibition
# ============================================================================

cat("\n=== EXPERIMENT 2: GO MF Analysis ===\n")
exp2 <- metadata[metadata$experiment == "Exp2", ]
baseline <- exp2[exp2$concentration == "0ngml_DMSO", ]
activin <- exp2[exp2$concentration == "15ngml_DMSO", ]

exp2_go <- list()

for (time in c(60, 120, 180)) {
  sb50 <- exp2[exp2$concentration == "SB50" & exp2$time_min == time, ]
  if (nrow(sb50) < 2) next

  # vs Baseline (0ngml_DMSO)
  name_bl <- paste0("SB50_", time, "min_vs_Baseline")
  cat(sprintf("  %s: ", name_bl))
  res_bl <- run_deseq(sb50, baseline)
  de_bl <- res_bl$gene[res_bl$padj < 0.05 & abs(res_bl$log2FoldChange) >= 1]
  cat(sprintf("%d DE -> ", length(de_bl)))

  go_bl <- run_go_mf(de_bl)
  if (!is.null(go_bl) && nrow(go_bl) > 0) {
    go_bl$comparison <- name_bl
    go_bl$reference <- "0ngml_DMSO"
    go_bl$time_min <- time
    exp2_go[[name_bl]] <- go_bl
    cat(sprintf("%d GO\n", sum(go_bl$p.adjust < 0.05)))
  } else cat("0 GO\n")

  # vs Activin (15ngml_DMSO)
  name_act <- paste0("SB50_", time, "min_vs_Activin")
  cat(sprintf("  %s: ", name_act))
  res_act <- run_deseq(sb50, activin)
  de_act <- res_act$gene[res_act$padj < 0.05 & abs(res_act$log2FoldChange) >= 1]
  cat(sprintf("%d DE -> ", length(de_act)))

  go_act <- run_go_mf(de_act)
  if (!is.null(go_act) && nrow(go_act) > 0) {
    go_act$comparison <- name_act
    go_act$reference <- "15ngml_DMSO"
    go_act$time_min <- time
    exp2_go[[name_act]] <- go_act
    cat(sprintf("%d GO\n", sum(go_act$p.adjust < 0.05)))
  } else cat("0 GO\n")
}

# Also: Activin vs Baseline
cat("  Activin_vs_Baseline: ")
res_act_bl <- run_deseq(activin, baseline)
de_act_bl <- res_act_bl$gene[res_act_bl$padj < 0.05 & abs(res_act_bl$log2FoldChange) >= 1]
cat(sprintf("%d DE -> ", length(de_act_bl)))
go_act_bl <- run_go_mf(de_act_bl)
if (!is.null(go_act_bl) && nrow(go_act_bl) > 0) {
  go_act_bl$comparison <- "Activin_vs_Baseline"
  go_act_bl$reference <- "0ngml_DMSO"
  go_act_bl$time_min <- NA
  exp2_go[["Activin_vs_Baseline"]] <- go_act_bl
  cat(sprintf("%d GO\n", sum(go_act_bl$p.adjust < 0.05)))
} else cat("0 GO\n")

exp2_go_df <- bind_rows(exp2_go) %>%
  filter(p.adjust < 0.05) %>%
  select(comparison, reference, time_min, ID, Description,
         GeneRatio, pvalue, p.adjust, Count, geneID)

# ============================================================================
# SAVE TO EXCEL
# ============================================================================

cat("\n=== Saving Excel files ===\n")

# Exp1 Excel
wb1 <- createWorkbook()
addWorksheet(wb1, "GO_MF_Results")
writeData(wb1, 1, exp1_go_df)
addWorksheet(wb1, "Summary")
exp1_summary <- exp1_go_df %>%
  group_by(concentration, time_min) %>%
  summarize(n_terms = n(), top_terms = paste(head(Description, 3), collapse = "; "), .groups = "drop")
writeData(wb1, 2, exp1_summary)
saveWorkbook(wb1, "exp1_go_mf_results.xlsx", overwrite = TRUE)
cat("Saved: exp1_go_mf_results.xlsx\n")

# Exp2 Excel
wb2 <- createWorkbook()
addWorksheet(wb2, "GO_MF_Results")
writeData(wb2, 1, exp2_go_df)
addWorksheet(wb2, "Summary")
exp2_summary <- exp2_go_df %>%
  group_by(comparison, reference) %>%
  summarize(n_terms = n(), top_terms = paste(head(Description, 3), collapse = "; "), .groups = "drop")
writeData(wb2, 2, exp2_summary)
saveWorkbook(wb2, "exp2_go_mf_results.xlsx", overwrite = TRUE)
cat("Saved: exp2_go_mf_results.xlsx\n")

# ============================================================================
# PLOTS - EXPERIMENT 1
# ============================================================================

cat("\n=== Generating plots ===\n")

if (nrow(exp1_go_df) > 0) {
  # Get top GO terms across all conditions
  top_terms <- exp1_go_df %>%
    group_by(Description) %>%
    summarize(min_p = min(p.adjust), n = n(), .groups = "drop") %>%
    slice_min(min_p, n = 20) %>%
    pull(Description)

  plot_df <- exp1_go_df %>%
    filter(Description %in% top_terms) %>%
    mutate(
      neg_log_p = -log10(p.adjust),
      concentration = factor(concentration, levels = c("5ngml", "10ngml", "15ngml")),
      time_min = factor(time_min)
    )

  # Heatmap-style dotplot
  p1 <- ggplot(plot_df, aes(x = interaction(concentration, time_min, sep = "\n"),
                             y = reorder(Description, neg_log_p))) +
    geom_point(aes(size = Count, color = neg_log_p)) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(padj)") +
    scale_size_continuous(range = c(2, 8), name = "Gene Count") +
    labs(x = "Concentration / Time (min)", y = "",
         title = "Exp1: GO Molecular Function across Activin conditions") +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          legend.position = "right")

  ggsave("exp1_go_mf_dotplot.pdf", p1, width = 12, height = 8)
  cat("Saved: exp1_go_mf_dotplot.pdf\n")

  # Faceted by concentration
  p1b <- ggplot(plot_df, aes(x = time_min, y = reorder(Description, neg_log_p))) +
    geom_point(aes(size = Count, color = neg_log_p)) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(padj)") +
    scale_size_continuous(range = c(2, 6), name = "Genes") +
    facet_wrap(~concentration, nrow = 1) +
    labs(x = "Time (min)", y = "", title = "Exp1: GO MF by Activin concentration") +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 7),
          strip.background = element_rect(fill = "gray90"))

  ggsave("exp1_go_mf_faceted.pdf", p1b, width = 14, height = 8)
  cat("Saved: exp1_go_mf_faceted.pdf\n")
}

# ============================================================================
# PLOTS - EXPERIMENT 2
# ============================================================================

if (nrow(exp2_go_df) > 0) {
  top_terms2 <- exp2_go_df %>%
    group_by(Description) %>%
    summarize(min_p = min(p.adjust), .groups = "drop") %>%
    slice_min(min_p, n = 20) %>%
    pull(Description)

  plot_df2 <- exp2_go_df %>%
    filter(Description %in% top_terms2) %>%
    mutate(
      neg_log_p = -log10(p.adjust),
      comparison = factor(comparison, levels = c("Activin_vs_Baseline",
                                                   "SB50_60min_vs_Baseline", "SB50_120min_vs_Baseline", "SB50_180min_vs_Baseline",
                                                   "SB50_60min_vs_Activin", "SB50_120min_vs_Activin", "SB50_180min_vs_Activin")),
      time_min = factor(time_min, levels = c(60, 120, 180))
    )

  # Dotplot comparing SB50 timepoints
  p2 <- ggplot(plot_df2, aes(x = comparison, y = reorder(Description, neg_log_p))) +
    geom_point(aes(size = Count, color = neg_log_p)) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(padj)") +
    scale_size_continuous(range = c(2, 8), name = "Gene Count") +
    labs(x = "", y = "", title = "Exp2: GO Molecular Function - SB50 treatment effects") +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 9, angle = 45, hjust = 1))

  ggsave("exp2_go_mf_dotplot.pdf", p2, width = 11, height = 8)
  cat("Saved: exp2_go_mf_dotplot.pdf\n")

  # Split by reference
  p2b <- ggplot(plot_df2 %>% filter(!is.na(time_min)),
                aes(x = time_min, y = reorder(Description, neg_log_p))) +
    geom_point(aes(size = Count, color = neg_log_p)) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(padj)") +
    scale_size_continuous(range = c(2, 6), name = "Genes") +
    facet_wrap(~reference, labeller = labeller(reference = c("0ngml_DMSO" = "vs Baseline",
                                                              "15ngml_DMSO" = "vs Activin"))) +
    labs(x = "SB50 addition time (min)", y = "",
         title = "Exp2: GO MF - SB50 effects over time") +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 7),
          strip.background = element_rect(fill = "gray90"))

  ggsave("exp2_go_mf_faceted.pdf", p2b, width = 12, height = 8)
  cat("Saved: exp2_go_mf_faceted.pdf\n")
}

cat("\n=== Done! ===\n")
cat("Excel files: exp1_go_mf_results.xlsx, exp2_go_mf_results.xlsx\n")
cat("Plots: exp1_go_mf_dotplot.pdf, exp1_go_mf_faceted.pdf\n")
cat("       exp2_go_mf_dotplot.pdf, exp2_go_mf_faceted.pdf\n")
