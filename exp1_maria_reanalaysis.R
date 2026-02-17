# Heatmap of Nodal Score Genes - Experiment 1
# Z-score expression grouped by time and concentration

source("preprocess.R")
library(pheatmap)
library(RColorBrewer)
library(gridExtra)

# Filter to Experiment 1 only
exp1_metadata <- metadata %>%
  filter(experiment == "Exp1", !is.na(concentration), !is.na(time_min))

cat("Experiment 1 samples:", nrow(exp1_metadata), "\n")

# Get counts for Exp1 samples
common_samples <- intersect(colnames(counts_filtered), rownames(exp1_metadata))
counts_exp1 <- counts_filtered[, common_samples]
exp1_metadata <- exp1_metadata[common_samples, ]

# Normalize counts (using simple design for normalization only)
dds <- DESeqDataSetFromMatrix(
  countData = counts_exp1,
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

ann_colors <- list(
  concentration = concentration_colors
)

# ── Prepare per-sample heatmap ──────────────────────────────────────────────
# Order samples by time then concentration
sample_order_df <- exp1_metadata %>%
  filter(time_min %in% c(60, 120, 180, 240)) %>%
  mutate(conc_order = match(concentration, c("0ngml", "5ngml", "10ngml", "15ngml"))) %>%
  arrange(time_min, conc_order)

ordered_samples <- intersect(rownames(sample_order_df), colnames(nodal_zscore))
sample_order_df <- sample_order_df[ordered_samples, ]

# Subset z-scores and remove genes with any NA values
nodal_zscore_clean <- nodal_zscore[, ordered_samples, drop = FALSE]
nodal_zscore_clean <- nodal_zscore_clean[complete.cases(nodal_zscore_clean), , drop = FALSE]
cat("Genes after removing NA rows:", nrow(nodal_zscore_clean), "\n")

# Per-sample column annotation
annotation_col_samples <- data.frame(
  time_min = factor(sample_order_df$time_min),
  concentration = factor(sample_order_df$concentration,
                         levels = c("0ngml", "5ngml", "10ngml", "15ngml")),
  row.names = rownames(sample_order_df)
)

# Time annotation colors (match present levels)
time_colors <- colorRampPalette(c("#F0F0F0", "#636363"))(length(levels(annotation_col_samples$time_min)))
names(time_colors) <- levels(annotation_col_samples$time_min)
ann_colors$time_min <- time_colors

# ── Calculate PDF dimensions ────────────────────────────────────────────────
n_genes <- nrow(nodal_zscore_clean)
n_samples_plot <- ncol(nodal_zscore_clean)

pdf_w <- n_samples_plot * 8 / 72 + 3   # cells + rownames + legend + margins
pdf_h <- n_genes * 10 / 72 + 1.5       # cells + annotation + margins

# ── Plot: Per-sample heatmap ────────────────────────────────────────────────
pdf(results_path("nodal_heatmap_exp1.pdf"), width = pdf_w, height = pdf_h, family = "Helvetica")
tryCatch(
  pheatmap(
    nodal_zscore_clean,
    annotation_col = annotation_col_samples,
    annotation_colors = ann_colors,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_names_col = FALSE,
    color = heatmap_colors,
    breaks = seq(-2.5, 2.5, length.out = 101),
    border_color = NA,
    fontsize = 10,
    fontsize_row = 9,
    main = "",
    annotation_legend = TRUE,
    legend = TRUE,
    cellwidth = 8,
    cellheight = 10,
    treeheight_row = 20
  ),
  error = function(e) message("Heatmap error: ", e$message),
  finally = dev.off()
)

