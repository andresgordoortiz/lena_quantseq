# Heatmap of Nodal Score Genes - Experiment 1
# Z-score expression grouped by time and concentration

library(DESeq2)
library(readr)
library(tidyverse)
library(readxl)
library(pheatmap)
library(RColorBrewer)

# Load counts
counts_raw <- read.table("salmon.merged.gene_counts.tsv", header = TRUE, row.names = 1)
counts_raw <- counts_raw[, -1]
counts_int <- round(counts_raw)

# Load sample metadata
samples <- read_csv("samples.csv")
parse_sample_names <- function(sample_name) {
  sub("^\\d{8}_R\\d+_", "", sample_name)
}
clean_treatments <- sapply(samples$sample_description, parse_sample_names, USE.NAMES = FALSE)

metadata <- data.frame(
  sample = paste0("S", samples$requests_sample_sample_id),
  treatment = clean_treatments,
  row.names = paste0("S", samples$requests_sample_sample_id)
)

metadata$experiment <- ifelse(grepl("DMSO|SB50", metadata$treatment), "Exp2", "Exp1")
metadata$concentration <- case_when(
  grepl("^0ngmlActivin", metadata$treatment) ~ "0ngml",
  grepl("^5ngmlActivin", metadata$treatment) ~ "5ngml",
  grepl("^10ngmlActivin", metadata$treatment) ~ "10ngml",
  grepl("^15ngmlActivin", metadata$treatment) ~ "15ngml",
  TRUE ~ NA_character_
)
metadata$time_min <- as.numeric(str_extract(metadata$treatment, "\\d+(?=min$)"))

# Filter to Experiment 1 only
exp1_metadata <- metadata %>%
  filter(experiment == "Exp1", !is.na(concentration), !is.na(time_min))

cat("Experiment 1 samples:", nrow(exp1_metadata), "\n")

# Get counts for Exp1 samples
common_samples <- intersect(colnames(counts_int), rownames(exp1_metadata))
counts_exp1 <- counts_int[, common_samples]
exp1_metadata <- exp1_metadata[common_samples, ]

# Filter low-expressed genes
keep <- rowSums(counts_exp1 >= 10) >= 3
counts_filtered <- counts_exp1[keep, ]

# Normalize counts (using simple design for normalization only)
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = exp1_metadata,
  design = ~ 1
)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

# Load nodal score genes
nodal_score_data <- read_excel("docs/nodal-score-genes_complete.xlsx", skip = 1)
nodal_genes <- unique(tolower(na.omit(nodal_score_data$`Nodal score`)))

# Find nodal genes in dataset (case-insensitive matching)
available_nodal <- intersect(tolower(rownames(norm_counts)), nodal_genes)
cat("Nodal score genes found:", length(available_nodal), "of", length(nodal_genes), "\n")

# Filter to nodal genes only (match case-insensitively but keep original names)
nodal_counts <- norm_counts[tolower(rownames(norm_counts)) %in% available_nodal, , drop = FALSE]

# Calculate z-scores across all samples
nodal_zscore <- t(scale(t(nodal_counts)))

# Remove any genes with all NA values after scaling
nodal_zscore <- nodal_zscore[!apply(nodal_zscore, 1, function(x) all(is.na(x))), , drop = FALSE]

# Create grouping variable
exp1_metadata$group <- paste0(exp1_metadata$concentration, "_", exp1_metadata$time_min, "min")

# Average z-scores by group
group_means <- data.frame(nodal_zscore) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "zscore") %>%
  left_join(exp1_metadata %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  group_by(gene, group) %>%
  summarise(mean_zscore = mean(zscore, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_zscore) %>%
  column_to_rownames("gene")

# Order columns: first by time, then by concentration within time
unique_groups <- exp1_metadata %>%
  select(group, time_min, concentration) %>%
  filter(time_min %in% c(60,120,180,240)) %>%
  distinct() %>%
  mutate(
    conc_order = match(concentration, c("0ngml", "5ngml", "10ngml", "15ngml"))
  ) %>%
  arrange(time_min, conc_order) %>%
  select(-conc_order)

ordered_cols <- intersect(unique_groups$group, colnames(group_means))
group_means_ordered <- group_means[, ordered_cols, drop = FALSE]

# Create annotation for columns
annotation_col <- unique_groups %>%
  filter(group %in% ordered_cols) %>%
  ungroup() %>%
  select(group, time_min, concentration)

# Convert to data frame with row names
annotation_col <- data.frame(
  time_min = annotation_col$time_min,
  concentration = annotation_col$concentration,
  row.names = annotation_col$group
)

annotation_col$time_min <- factor(annotation_col$time_min)
annotation_col$concentration <- factor(annotation_col$concentration,
                                       levels = c("0ngml", "5ngml", "10ngml", "15ngml"))

# Create clean, academic color palettes
# Heatmap colors: refined blue-white-red
heatmap_colors <- colorRampPalette(c(
  "#2166AC",  # Deep blue
  "#4393C3",  # Medium blue
  "#92C5DE",  # Light blue
  "#F7F7F7",  # White
  "#F4A582",  # Light red
  "#D6604D",  # Medium red
  "#B2182B"   # Deep red
))(100)

# Concentration colors: grey for 0, shades of green for 5/10/15
green_base <- c(95, 179, 88) / 255  # Normalize RGB to 0-1
concentration_colors <- c(
  "0ngml" = "#BDBDBD",  # Grey for control
  "5ngml" = rgb(green_base[1] * 0.6, green_base[2] * 0.6, green_base[3] * 0.6),  # Light green
  "10ngml" = rgb(green_base[1] * 0.8, green_base[2] * 0.8, green_base[3] * 0.8),  # Medium green
  "15ngml" = rgb(green_base[1], green_base[2], green_base[3])  # Full green
)

# Time colors: subtle gradient
time_colors <- colorRampPalette(c("#F0F0F0", "#636363"))(length(unique(annotation_col$time_min)))
names(time_colors) <- levels(annotation_col$time_min)

ann_colors <- list(
  time_min = time_colors,
  concentration = concentration_colors
)



# Remove any rows with missing values
group_means_clean <- group_means_ordered[complete.cases(group_means_ordered), , drop = FALSE]
cat("Genes after removing NA rows:", nrow(group_means_clean), "\n")

pdf("nodal_heatmap_exp1_averaged.pdf", width = 10, height = 8, family = "Helvetica")
pheatmap(
  group_means_clean,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = FALSE,  # Cleaner without individual column labels
  annotation_names_col = FALSE,  # Cleaner annotation
  color = heatmap_colors,
  breaks = seq(-2.5, 2.5, length.out = 101),
  border_color = NA,  # No borders for cleaner look
  fontsize = 10,
  fontsize_row = 9,
  main = "",  # No title for cleaner look
  annotation_legend = TRUE,
  legend = TRUE,
  cellwidth = 12,
  cellheight = 10,
  treeheight_row = 30
)
dev.off()

# Save the averaged z-score matrix
#write_csv(group_means_clean %>% rownames_to_column("gene"))

