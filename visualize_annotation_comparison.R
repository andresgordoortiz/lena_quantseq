# ============================================================================
# Visualize Genome Annotation Comparison: GRCz11 vs GRCz12tu
# ============================================================================
#
# This script visualizes the output from comapre_annotations.sh
# Run the bash script first to generate the required input files
#
# Usage:
#   1. Set ANNOTATION_DIR to the output directory from the bash script
#   2. (Optional) Add your genes of interest to MY_GENES vector
#   3. Run the script AFTER RUNNING Q3 SCRIPT
#

source("preprocess.R")

# ============================================================================
# CONFIGURATION
# ============================================================================

all_activin<-rbind(res_sb50_60_vs_activin_only,res_sb50_120_vs_activin_only,
                   res_sb50_180_vs_activin_only,res_activin_vs_baseline,
                   res_sb50_60_vs_baseline,res_sb50_120_vs_baseline,
                   res_sb50_180_vs_baseline)

activin_affected_genes<-all_activin %>%
  filter(abs(log2FoldChange)>=1.5 & padj <=0.05) %>%
  dplyr::select(gene) %>%
  as.list() %>%
  unique()
# Output directory for annotation figures
ANNOTATION_OUT_DIR <- results_path("compare_annotations")
if (!dir.exists(ANNOTATION_OUT_DIR)) dir.create(ANNOTATION_OUT_DIR, recursive = TRUE)

# Your genes of interest (add gene names or IDs to check against new annotations)
# Example: MY_GENES <- c("nodal", "lefty1", "lefty2", "pitx2", "gsc", "mixl1")
MY_GENES <- activin_affected_genes

# Set input directory to the annotation comparison output
ANNOTATION_DIR <- "compare_annotations"

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading annotation comparison data from:", ANNOTATION_DIR, "\n\n")

# Check if directory exists
if (!dir.exists(ANNOTATION_DIR)) {
  stop("Directory not found: ", ANNOTATION_DIR,
       "\nPlease run comapre_annotations.sh first and update ANNOTATION_DIR")
}

# Load gene lists
grcz11_genes <- read_tsv(
  file.path(ANNOTATION_DIR, "GRCz11_genes.txt"),
  col_names = c("gene_id", "gene_name", "biotype", "chr", "start", "end", "strand"),
  show_col_types = FALSE
)

grcz12tu_genes <- read_tsv(
  file.path(ANNOTATION_DIR, "GRCz12tu_genes.txt"),
  col_names = c("gene_id", "gene_name", "biotype", "chr", "start", "end", "strand"),
  show_col_types = FALSE
)

genes_removed <- read_tsv(
  file.path(ANNOTATION_DIR, "genes_removed_in_GRCz12tu.txt"),
  col_names = c("gene_id", "gene_name", "biotype", "chr", "start", "end", "strand"),
  show_col_types = FALSE
)

genes_added <- read_tsv(
  file.path(ANNOTATION_DIR, "genes_added_in_GRCz12tu.txt"),
  col_names = c("gene_id", "gene_name", "biotype", "chr", "start", "end", "strand"),
  show_col_types = FALSE
)

genes_new_names <- read_tsv(
  file.path(ANNOTATION_DIR, "genes_with_new_names.txt"),
  col_names = c("gene_id", "old_name", "arrow", "new_name", "location"),
  show_col_types = FALSE
) %>% dplyr::select(-arrow)

genes_biotype_changes <- read_tsv(
  file.path(ANNOTATION_DIR, "genes_with_biotype_changes.txt"),
  col_names = c("gene_id", "gene_name", "old_biotype", "arrow", "new_biotype"),
  show_col_types = FALSE
) %>% dplyr::select(-arrow)

newly_annotated_named <- read_tsv(
  file.path(ANNOTATION_DIR, "newly_annotated_with_names.txt"),
  col_names = c("gene_id", "gene_name", "biotype", "location"),
  show_col_types = FALSE
)

newly_annotated_computational <- read_tsv(
  file.path(ANNOTATION_DIR, "newly_annotated_computational.txt"),
  col_names = c("gene_id", "gene_name", "biotype", "location"),
  show_col_types = FALSE
)

# ============================================================================
# COLOR PALETTE (consistent with q3 style)
# ============================================================================

colors <- list(
  primary = "#2166AC",
  secondary = "#B2182B",
  accent1 = "#1B9E77",
  accent2 = "#D95F02",
  accent3 = "#7570B3",
  neutral = "grey70",
  light = "grey90"
)

biotype_colors <- c(

  "protein_coding" = "#2166AC",
  "lncRNA" = "#B2182B",
  "pseudogene" = "#7570B3",
  "miRNA" = "#1B9E77",
  "snRNA" = "#D95F02",
  "snoRNA" = "#E6AB02",
  "rRNA" = "#A6761D",
  "tRNA" = "#666666",
  "other" = "grey70"
)

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("=== ANNOTATION COMPARISON SUMMARY ===\n\n")
cat(paste0("GRCz11 total genes:     ", format(nrow(grcz11_genes), big.mark = ","), "\n"))
cat(paste0("GRCz12tu total genes:   ", format(nrow(grcz12tu_genes), big.mark = ","), "\n"))
net_change <- nrow(grcz12tu_genes) - nrow(grcz11_genes)
cat(paste0("Net change:             ", ifelse(net_change >= 0, "+", ""), format(net_change, big.mark = ","), "\n\n"))

cat(paste0("Genes removed:          ", format(nrow(genes_removed), big.mark = ","), "\n"))
cat(paste0("Genes added:            ", format(nrow(genes_added), big.mark = ","), "\n"))
cat(paste0("Genes with new names:   ", format(nrow(genes_new_names), big.mark = ","), "\n"))
cat(paste0("Biotype changes:        ", format(nrow(genes_biotype_changes), big.mark = ","), "\n\n"))

# ============================================================================
# PLOT 1: Overall Gene Count Comparison
# ============================================================================

summary_data <- tibble(
  category = c("Total genes", "Removed", "Added", "New names", "Biotype changed"),
  GRCz11 = c(nrow(grcz11_genes), nrow(genes_removed), 0, 0, nrow(genes_biotype_changes)),
  GRCz12tu = c(nrow(grcz12tu_genes), 0, nrow(genes_added), nrow(genes_new_names), nrow(genes_biotype_changes))
)

comparison_data <- tibble(
  category = factor(
    c("Removed in\nGRCz12tu", "Added in\nGRCz12tu", "Got official\nname", "Biotype\nchanged"),
    levels = c("Removed in\nGRCz12tu", "Added in\nGRCz12tu", "Got official\nname", "Biotype\nchanged")
  ),
  count = c(nrow(genes_removed), nrow(genes_added), nrow(genes_new_names), nrow(genes_biotype_changes)),
  type = c("Removed", "Added", "Improved", "Changed")
)

p1 <- ggplot(comparison_data, aes(x = category, y = count, fill = type)) +

geom_col(width = 0.7) +
  geom_text(aes(label = format(count, big.mark = ",")),
            vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c(
    "Removed" = colors$secondary,
    "Added" = colors$accent1,
    "Improved" = colors$primary,
    "Changed" = colors$accent3
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Annotation Changes: GRCz11 → GRCz12tu",
    subtitle = paste0("Net change: ", ifelse(net_change >= 0, "+", ""), format(net_change, big.mark = ","), " genes"),
    x = NULL,
    y = "Number of genes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 9)
  )

# ============================================================================
# PLOT 2: Biotype Distribution Comparison
# ============================================================================

# Prepare biotype data
simplify_biotype <- function(biotype) {
  case_when(
    biotype == "protein_coding" ~ "protein_coding",
    biotype %in% c("lncRNA", "lnc_RNA") ~ "lncRNA",
    grepl("pseudogene", biotype) ~ "pseudogene",
    biotype == "miRNA" ~ "miRNA",
    biotype == "snRNA" ~ "snRNA",
    biotype == "snoRNA" ~ "snoRNA",
    biotype == "rRNA" ~ "rRNA",
    biotype == "tRNA" ~ "tRNA",
    TRUE ~ "other"
  )
}

biotype_grcz11 <- grcz11_genes %>%
  mutate(biotype_simple = simplify_biotype(biotype)) %>%
  count(biotype_simple) %>%
  mutate(version = "GRCz11", pct = 100 * n / sum(n))

biotype_grcz12tu <- grcz12tu_genes %>%
  mutate(biotype_simple = simplify_biotype(biotype)) %>%
  count(biotype_simple) %>%
  mutate(version = "GRCz12tu", pct = 100 * n / sum(n))

biotype_combined <- bind_rows(biotype_grcz11, biotype_grcz12tu) %>%
  mutate(
    biotype_simple = factor(biotype_simple,
      levels = c("protein_coding", "lncRNA", "pseudogene", "miRNA",
                 "snRNA", "snoRNA", "rRNA", "tRNA", "other"))
  )

p2 <- ggplot(biotype_combined, aes(x = version, y = pct, fill = biotype_simple)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(
    aes(label = ifelse(pct > 3, sprintf("%.1f%%", pct), "")),
    position = position_stack(vjust = 0.5),
    size = 2.5, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(values = biotype_colors, name = "Biotype") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "Biotype Distribution",
    x = NULL,
    y = "% of genes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major.x = element_blank()
  )

# ============================================================================
# PLOT 3: New Genes by Biotype
# ============================================================================

new_genes_biotype <- genes_added %>%
  mutate(
    biotype_simple = simplify_biotype(biotype),
    annotation_type = ifelse(
      grepl("^LOC|^ENSDARG", gene_name) | gene_name == "NA",
      "Computational ID",
      "Official name"
    )
  ) %>%
  count(biotype_simple, annotation_type) %>%
  mutate(biotype_simple = factor(biotype_simple,
    levels = c("protein_coding", "lncRNA", "pseudogene", "miRNA",
               "snRNA", "snoRNA", "rRNA", "tRNA", "other")))

p3 <- ggplot(new_genes_biotype, aes(x = biotype_simple, y = n, fill = annotation_type)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(
    "Official name" = colors$primary,
    "Computational ID" = colors$neutral
  ), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Newly Added Genes in GRCz12tu",
    subtitle = paste0(format(nrow(genes_added), big.mark = ","), " new genes total"),
    x = "Biotype",
    y = "Number of genes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ============================================================================
# PLOT 4: Biotype Transitions (Sankey-style simplified)
# ============================================================================

if (nrow(genes_biotype_changes) > 0) {
  biotype_transitions <- genes_biotype_changes %>%
    mutate(
      old_simple = simplify_biotype(old_biotype),
      new_simple = simplify_biotype(new_biotype)
    ) %>%
    count(old_simple, new_simple, name = "count") %>%
    arrange(desc(count)) %>%
    head(15)  # Top 15 transitions

  biotype_transitions <- biotype_transitions %>%
    mutate(transition = paste(old_simple, "→", new_simple)) %>%
    mutate(transition = fct_reorder(transition, count))

  p4 <- ggplot(biotype_transitions, aes(x = count, y = transition, fill = old_simple)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = count), hjust = -0.2, size = 3) +
    scale_fill_manual(values = biotype_colors, name = "Original biotype") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Top Biotype Transitions",
      subtitle = paste0(format(nrow(genes_biotype_changes), big.mark = ","), " genes changed biotype"),
      x = "Number of genes",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
      legend.position = "none",
      panel.grid.major.y = element_blank()
    )
} else {
  p4 <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No biotype changes", size = 5) +
    theme_void()
}

# ============================================================================
# COMBINED FIGURE
# ============================================================================

combined_fig <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Zebrafish Genome Annotation Comparison",
    subtitle = "GRCz11 → GRCz12tu",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40")
    )
  )

ggsave(
  file.path(ANNOTATION_OUT_DIR, "annotation_comparison_figure.pdf"),
  combined_fig, width = 12, height = 10
)
cat("\nSaved:", file.path(ANNOTATION_OUT_DIR, "annotation_comparison_figure.pdf"), "\n")

# Also save as PNG for quick viewing
ggsave(
  file.path(ANNOTATION_OUT_DIR, "annotation_comparison_figure.png"),
  combined_fig, width = 12, height = 10, dpi = 150
)
cat("Saved:", file.path(ANNOTATION_OUT_DIR, "annotation_comparison_figure.png"), "\n")

# ============================================================================
# ============================================================================
# CHECK YOUR GENES OF INTEREST
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("CHECKING YOUR GENES OF INTEREST\n")
cat(strrep("=", 70), "\n\n")

if (length(MY_GENES) == 0) {
  cat("No genes specified in MY_GENES vector.\n")
  cat("Add gene names/IDs to MY_GENES at the top of the script to check them.\n\n")
} else {
  cat(sprintf("Checking %d genes: %s\n\n", length(MY_GENES), paste(MY_GENES, collapse = ", ")))

  # Combine all new genes for searching
  all_new_genes <- bind_rows(
    newly_annotated_named %>% mutate(annotation_status = "Official name"),
    newly_annotated_computational %>% mutate(annotation_status = "Computational ID")
  )

  # Search by gene name (case-insensitive)
  found_in_new <- all_new_genes %>%
    filter(
      tolower(gene_name) %in% tolower(MY_GENES) |
      tolower(gene_id) %in% tolower(MY_GENES)
    )

  # Search in genes that got new official names
  found_new_names <- genes_new_names %>%
    filter(
      tolower(old_name) %in% tolower(MY_GENES) |
      tolower(new_name) %in% tolower(MY_GENES) |
      tolower(gene_id) %in% tolower(MY_GENES)
    )

  # Search in removed genes
  found_removed <- genes_removed %>%
    filter(
      tolower(gene_name) %in% tolower(MY_GENES) |
      tolower(gene_id) %in% tolower(MY_GENES)
    )

  # Report findings
  cat("=== IN NEWLY ADDED GENES ===\n")
  if (nrow(found_in_new) > 0) {
    cat(sprintf("Found %d of your genes in newly added annotations:\n", nrow(found_in_new)))
    print(found_in_new %>% dplyr::select(gene_id, gene_name, biotype, annotation_status))
  } else {
    cat("None of your genes are in the newly added genes.\n")
  }

  cat("\n=== IN GENES WITH NEW OFFICIAL NAMES ===\n")
  if (nrow(found_new_names) > 0) {
    cat(sprintf("Found %d of your genes that received new official names:\n", nrow(found_new_names)))
    print(found_new_names %>% dplyr::select(gene_id, old_name, new_name, location))
  } else {
    cat("None of your genes received new official names.\n")
  }

  cat("\n=== IN REMOVED GENES ===\n")
  if (nrow(found_removed) > 0) {
    cat(sprintf("WARNING: %d of your genes were REMOVED in GRCz12tu:\n", nrow(found_removed)))
    print(found_removed %>% dplyr::select(gene_id, gene_name, biotype))
  } else {
    cat("None of your genes were removed (good!).\n")
  }

  # Summary
  genes_found <- unique(c(
    found_in_new$gene_name, found_in_new$gene_id,
    found_new_names$old_name, found_new_names$new_name, found_new_names$gene_id,
    found_removed$gene_name, found_removed$gene_id
  ))
  genes_found <- genes_found[tolower(genes_found) %in% tolower(MY_GENES)]

  genes_not_found <- MY_GENES[!tolower(MY_GENES) %in% tolower(genes_found)]

  cat("\n=== SUMMARY ===\n")
  cat(sprintf("Genes found in annotation changes: %d / %d\n",
              length(unique(genes_found)), length(MY_GENES)))

  if (length(genes_not_found) > 0) {
    cat(sprintf("Genes not affected by annotation changes: %s\n",
                paste(genes_not_found, collapse = ", ")))
    cat("(These genes likely exist in both versions with no changes)\n")
  }

  # Save results
  gene_check_results <- tibble(
    query_gene = MY_GENES,
    found_in_new_genes = tolower(MY_GENES) %in% tolower(c(found_in_new$gene_name, found_in_new$gene_id)),
    got_new_name = tolower(MY_GENES) %in% tolower(c(found_new_names$old_name, found_new_names$new_name, found_new_names$gene_id)),
    was_removed = tolower(MY_GENES) %in% tolower(c(found_removed$gene_name, found_removed$gene_id))
  )

  write_csv(gene_check_results, file.path(ANNOTATION_OUT_DIR, "my_genes_annotation_check.csv"))
  cat("\nSaved:", file.path(ANNOTATION_OUT_DIR, "my_genes_annotation_check.csv"), "\n")
}

# ============================================================================
# EXPORT SUMMARY TABLES
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("EXPORTING SUMMARY TABLES\n")
cat(strrep("=", 70), "\n\n")

# Top newly annotated genes with official names
if (nrow(newly_annotated_named) > 0) {
  write_csv(
    newly_annotated_named %>% arrange(gene_name),
    file.path(ANNOTATION_OUT_DIR, "newly_annotated_genes_with_names.csv")
  )
  cat(paste0("Exported ", format(nrow(newly_annotated_named), big.mark = ","), " newly annotated genes with official names\n"))
}

# Genes that got new official names
if (nrow(genes_new_names) > 0) {
  write_csv(
    genes_new_names %>% arrange(new_name),
    file.path(ANNOTATION_OUT_DIR, "genes_renamed_to_official.csv")
  )
  cat(paste0("Exported ", format(nrow(genes_new_names), big.mark = ","), " genes that received official names\n"))
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("DONE!\n")
cat(strrep("=", 70), "\n")
