# GO Term Gene Direction Analysis - Simple Boxplot
# Shows log2FC direction for genes in signaling receptor GO terms

source("preprocess.R")
library(org.Dr.eg.db)
library(openxlsx)
library(AnnotationDbi)

# ============================================================================
# DATA PREP (uses shared counts_filtered, metadata from preprocess.R)
# ============================================================================

counts_filt <- counts_filtered

# Use shared run_deseq from preprocess.R (with lfcShrink)
run_deseq_go <- function(test_meta, ref_meta) {
  run_deseq(test_meta, ref_meta, counts_mat = counts_filt)
}

# ============================================================================
# GET ALL GENES FROM TARGET GO TERMS (from org.Dr.eg.db)
# ============================================================================

library(AnnotationDbi)

# Target GO term IDs and their names
target_go_info <- data.frame(
  GO_ID = c("GO:0030545", "GO:0030546", "GO:0048018"),
  GO_term = c("signaling receptor regulator activity",
              "signaling receptor activator activity",
              "receptor ligand activity")
)

cat("=== Extracting genes from org.Dr.eg.db ===\n")

# Get all genes annotated to these GO terms
go_genes_list <- lapply(1:nrow(target_go_info), function(i) {
  go_id <- target_go_info$GO_ID[i]
  go_name <- target_go_info$GO_term[i]

  # Get Entrez IDs for this GO term
  entrez_ids <- tryCatch({
    AnnotationDbi::select(org.Dr.eg.db,
                          keys = go_id,
                          columns = "ENTREZID",
                          keytype = "GOALL")$ENTREZID
  }, error = function(e) character(0))

  if (length(entrez_ids) == 0) return(data.frame())

  # Convert to gene symbols
  symbols <- AnnotationDbi::select(org.Dr.eg.db,
                                   keys = unique(entrez_ids),
                                   columns = "SYMBOL",
                                   keytype = "ENTREZID")

  data.frame(
    GO_ID = go_id,
    GO_term = go_name,
    gene = symbols$SYMBOL
  )
})

all_go_genes_df <- bind_rows(go_genes_list) %>%
  filter(!is.na(gene)) %>%
  distinct()

# Get unique genes (some may be in multiple terms)
all_go_genes <- unique(all_go_genes_df$gene)

cat("\nGO terms and gene counts:\n")
all_go_genes_df %>%
  group_by(GO_ID, GO_term) %>%
  summarize(n_genes = n_distinct(gene), .groups = "drop") %>%
  print()

cat("\nTotal unique genes across all 3 terms:", length(all_go_genes), "\n")

# Filter to genes that are in our dataset
genes_in_data <- intersect(all_go_genes, rownames(counts_filt))
cat("Genes present in our dataset:", length(genes_in_data), "\n")

# ============================================================================
# RUN DE ANALYSIS FOR ALL EXP2 COMPARISONS
# ============================================================================

cat("\n=== Running DE analysis ===\n")
exp2 <- metadata[metadata$experiment == "Exp2", ]
baseline <- exp2[exp2$concentration == "0ngml_DMSO", ]
activin <- exp2[exp2$concentration == "15ngml_DMSO", ]

de_list <- list()

# Activin vs Baseline
cat("Activin vs Baseline\n")
de_list[["Activin"]] <- run_deseq_go(activin, baseline) %>%
  mutate(condition = "Activin")

# SB50 vs Baseline at each timepoint
for (time in c(60, 120, 180)) {
  sb50 <- exp2[exp2$concentration == "SB50" & exp2$time_min == time, ]

  name <- paste0("SB50_", time, "'\nvs Baseline")
  cat(name, "\n")
  de_list[[name]] <- run_deseq_go(sb50, baseline) %>%
    mutate(condition = name)
}

# SB50 vs Activin at each timepoint
for (time in c(60, 120, 180)) {
  sb50 <- exp2[exp2$concentration == "SB50" & exp2$time_min == time, ]

  name <- paste0("SB50_", time, "'\nvs Activin")
  cat(name, "\n")
  de_list[[name]] <- run_deseq_go(sb50, activin) %>%
    mutate(condition = name)
}

all_de <- bind_rows(de_list)

# ============================================================================
# FILTER TO GO TERM GENES (using all genes from database)
# ============================================================================

plot_data <- all_de %>%
  filter(gene %in% genes_in_data) %>%
  mutate(
    condition = factor(condition, levels = c(
      "Activin",
      "SB50_60'\nvs Baseline", "SB50_120'\nvs Baseline", "SB50_180'\nvs Baseline",
      "SB50_60'\nvs Activin", "SB50_120'\nvs Activin", "SB50_180'\nvs Activin"
    )),
    significant = padj < 0.05 & abs(log2FoldChange) >= 1
  )

cat("\n=== Summary ===\n")
plot_data %>%
  group_by(condition) %>%
  summarize(
    n_genes = n(),
    mean_log2FC = round(mean(log2FoldChange), 2),
    n_up = sum(log2FoldChange > 0 & significant),
    n_down = sum(log2FoldChange < 0 & significant),
    .groups = "drop"
  ) %>%
  print()

# ============================================================================
# CALCULATE ENRICHMENT RATIO (GO term genes / total DE genes)
# ============================================================================

# Define condition order explicitly
condition_order <- c(
  "Activin",
  "SB50_60'\nvs Baseline", "SB50_120'\nvs Baseline", "SB50_180'\nvs Baseline",
  "SB50_60'\nvs Activin", "SB50_120'\nvs Activin", "SB50_180'\nvs Activin"
)

# First, get TOTAL significant DE genes per condition (not just GO term genes)
total_de_per_condition <- all_de %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1) %>%
  group_by(condition) %>%
  summarize(total_de = n(), .groups = "drop") %>%
  mutate(condition = factor(condition, levels = condition_order))

# Count significant GO term genes per condition and direction
sig_counts <- plot_data %>%
  filter(significant) %>%
  mutate(
    direction = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated")
  ) %>%
  group_by(condition, direction) %>%
  summarize(n_go_genes = n(), .groups = "drop") %>%
  complete(condition, direction, fill = list(n_go_genes = 0)) %>%
  left_join(total_de_per_condition, by = "condition") %>%
  mutate(
    condition = factor(condition, levels = condition_order),
    # Proportion of DE genes that are from GO term
    prop = n_go_genes / total_de * 100,
    direction = factor(direction, levels = c("Upregulated", "Downregulated")),
    # For plotting downregulated as negative
    prop_plot = ifelse(direction == "Downregulated", -prop, prop)
  )

# Calculate cumulative (total) proportion per condition for the line
cumulative_data <- sig_counts %>%
  group_by(condition) %>%
  summarize(
    total_prop = sum(prop),
    .groups = "drop"
  ) %>%
  mutate(condition = factor(condition, levels = condition_order))

cat("\n=== Enrichment analysis ===\n")
cat("Total genes in GO terms (in dataset):", length(genes_in_data), "\n")
cat("Total genes tested:", nrow(counts_filt), "\n")
cat("Expected proportion if random:", round(length(genes_in_data)/nrow(counts_filt)*100, 2), "%\n\n")

summary_table <- sig_counts %>%
  group_by(condition) %>%
  summarize(
    total_de = first(total_de),
    go_term_de = sum(n_go_genes),
    pct_go_term = round(sum(n_go_genes) / first(total_de) * 100, 1),
    n_up = sum(n_go_genes[direction == "Upregulated"]),
    n_down = sum(n_go_genes[direction == "Downregulated"]),
    .groups = "drop"
  ) %>%
  mutate(condition = factor(condition, levels = condition_order)) %>%
  arrange(condition)
print(summary_table)

# ============================================================================
# PUBLICATION PLOT: Proportion of DE genes from GO terms (shows enrichment)
# ============================================================================

# Calculate scale factor for secondary axis (to fit cumulative line nicely)
max_bar <- max(abs(sig_counts$prop_plot), na.rm = TRUE)
max_cumulative <- max(cumulative_data$total_prop, na.rm = TRUE)
scale_factor <- max_bar / max_cumulative * 1.5

p <- ggplot(sig_counts, aes(x = condition, y = prop_plot)) +
  # Reference line at 0
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.5) +
  # Background rectangles for grouping
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = "#B2182B", alpha = 0.08) +
  annotate("rect", xmin = 1.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = "#2166AC", alpha = 0.08) +
  annotate("rect", xmin = 4.5, xmax = 7.5, ymin = -Inf, ymax = Inf,
           fill = "#1B7837", alpha = 0.08) +
  # Bars
  geom_col(aes(fill = direction), width = 0.7, color = "white", linewidth = 0.3) +
  # Add cumulative line (scaled to fit on same axis)
  geom_line(data = cumulative_data,
            aes(x = condition, y = total_prop * scale_factor, group = 1),
            color = "#7B3294", linewidth = 1.2) +
  geom_point(data = cumulative_data,
             aes(x = condition, y = total_prop * scale_factor),
             color = "#7B3294", size = 3) +
  # Add cumulative percentage labels
  geom_text(data = cumulative_data,
            aes(x = condition, y = total_prop * scale_factor + 0.4,
                label = paste0(round(total_prop, 1), "%")),
            color = "#7B3294", size = 3, fontface = "bold") +
  # Add percentage labels on bars
  geom_text(aes(label = paste0(round(abs(prop), 1), "%\n(", n_go_genes, ")"),
                y = ifelse(direction == "Upregulated", prop_plot + 0.3, prop_plot - 0.3)),
            size = 2.8, fontface = "bold", color = "grey20", lineheight = 0.8) +
  # Colors
  scale_fill_manual(
    values = c("Upregulated" = "#B2182B", "Downregulated" = "#2166AC"),
    name = ""
  ) +
  # Secondary axis for cumulative
  scale_y_continuous(
    sec.axis = sec_axis(~ . / scale_factor, name = "Total % (cumulative)")
  ) +
  # Labels
  labs(
    x = "",
    y = "% of DE genes from GO terms",
    title = "Signaling receptor GO terms: proportion of DE genes & direction",
    subtitle = "Purple line = total enrichment (up + down). Higher % = stronger enrichment."
  ) +
  # Theme
  theme_minimal(base_size = 11, base_family = "Helvetica") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.6),

    axis.text.x = element_text(size = 9, color = "grey20", face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 11, face = "bold", color = "#7B3294", margin = margin(l = 10)),
    axis.text.y.right = element_text(color = "#7B3294"),

    legend.position = "top",

    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
  )

ggsave(results_path("exp2_go_genes_direction.pdf"), p, width = 8, height = 8, device = pdf)
cat("\nSaved:", results_path("exp2_go_genes_direction.pdf"), "\n")

# Also save a version showing the actual significant genes
sig_genes <- plot_data %>%
  filter(significant) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
  dplyr::select(gene, condition, log2FoldChange, padj, direction) %>%
  arrange(condition, direction, desc(abs(log2FoldChange)))

write.csv(sig_genes, results_path("exp2_go_genes_significant.csv"), row.names = FALSE)
cat("Saved:", results_path("exp2_go_genes_significant.csv"), "\n")

# Save the summary table
write.csv(summary_table, results_path("exp2_go_enrichment_summary.csv"), row.names = FALSE)
cat("Saved:", results_path("exp2_go_enrichment_summary.csv"), "\n")
