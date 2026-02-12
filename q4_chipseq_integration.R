# Question 4: ChIP-seq Integration Analysis
#
# Integrate Smad2 and EomesA ChIP-seq data (from Nelson et al. 2014,
# BMC Genomics) with the Quant-Seq RNA-seq results.
#
# CHIPSEQ DATA:
#   - 12915_2014_81_MOESM4_ESM.xlsx  → Smad2 ChIP-seq peaks (Table S2)
#   - 12915_2014_81_MOESM10_ESM.xlsx → EomesA ChIP-seq peaks (Table S5)
#   Both contain: Chromosome, Peak Start/Stop, Peak ID, Proximal gene,
#                 Foxh1 binding site (motif presence)
#   Genes marked with ^ in Smad2 data are upregulated by ndr1 overexpression
#
# QUESTIONS:
#   1. For each candidate gene, report Smad2, EomesA, and Foxh1 ChIP-seq binding
#   2. Compare temporal dynamics of candidate genes vs. Nodal score (direct targets)
#   3. Compute reversibility in Exp2 (SB50 conditions) for candidates vs. Nodal score
#

library(DESeq2)
library(readr)
library(tidyverse)
library(readxl)
library(ggplot2)
library(patchwork)
library(ggrepel)

# ============================================================================
# LOAD AND PREPARE DATA (consistent with previous scripts)
# ============================================================================

counts_raw <- read.table("salmon.merged.gene_counts.tsv", header = TRUE, row.names = 1)
counts_raw <- counts_raw[, -1]
counts_int <- round(counts_raw)

samples <- read_csv("samples.csv", show_col_types = FALSE)
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
  grepl("DMSO", metadata$treatment) & grepl("^0ngml", metadata$treatment) ~ "0ngml_DMSO",
  grepl("DMSO", metadata$treatment) & grepl("^15ngml", metadata$treatment) ~ "15ngml_DMSO",
  grepl("SB50", metadata$treatment) ~ "SB50",
  TRUE ~ NA_character_
)
metadata$time_min <- as.numeric(str_extract(metadata$treatment, "\\d+(?=min$)"))

common_samples <- intersect(colnames(counts_int), rownames(metadata))
counts_int <- counts_int[, common_samples]
metadata <- metadata[common_samples, ]

keep <- rowSums(counts_int >= 10) >= 3
counts_filtered <- counts_int[keep, ]

cat("Total samples:", ncol(counts_filtered), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")

# ============================================================================
# DEFINE GENE LISTS
# ============================================================================
#
# TWO GENE SETS:
#   1. NODAL SCORE (27 genes): Curated from the literature as established
#      *direct* targets of Nodal/Activin signalling (e.g. gsc, lft1, ndr1,
#      sox32). These are the canonical downstream effectors, many of which
#      are known to be bound by Smad2 and its co-factors (Foxh1, EomesA).
#
#   2. CANDIDATE LIST (18 genes): Genes of interest identified in Lena's
#      experiments as Activin-responsive, but whose status as direct vs
#      indirect targets is uncertain. By overlaying ChIP-seq binding data
#      we can assess which candidates have evidence of direct TF binding
#      (Smad2/EomesA peaks ± Foxh1 motifs), and by comparing their temporal
#      dynamics and reversibility to the Nodal score we can gauge whether
#      they behave like primary or secondary response genes.
#
#   OVERLAP: Some genes appear in both lists (dkk1b, kirrel3l, flrt3,
#   efnb2a). These are handled by assigning each gene to ONE group only:
#   if it is in the Nodal score it is labelled "Nodal score", otherwise
#   "Candidate". This avoids duplication artifacts in downstream analyses.
#

candidate_genes <- c("rhov", "net1", "flrt3", "rnd1b", "dkk1b", "plekha5b",
                     "rasgef1ba", "efna1a", "osr1", "frmd4ba", "abi1b",
                     "snai1b", "snai1a", "kirrel3l", "jcada", "pfkfb3",
                     "prickle1b", "efnb2a")

nodal_score_data <- read_excel("docs/nodal-score-genes_complete.xlsx", skip = 1)
nodal_genes <- unique(tolower(na.omit(nodal_score_data$`Nodal score`)))

# Genes shared between both lists
overlapping_genes <- intersect(candidate_genes, nodal_genes)
cat("Overlapping genes:", paste(overlapping_genes, collapse = ", "), "\n")

# For expression analyses, assign each gene to ONE group to avoid duplication:
# overlapping genes go to "Nodal score" (the established set)
candidate_genes_unique <- setdiff(candidate_genes, nodal_genes)

cat("Candidate genes:", length(candidate_genes), "(", length(candidate_genes_unique), "unique)\n")
cat("Nodal score genes:", length(nodal_genes), "\n")

# ============================================================================
# PART 1: ChIP-seq binding summary for candidate genes
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("PART 1: ChIP-seq Binding Summary for Candidate Genes\n")
cat(strrep("=", 70), "\n\n")

# Read ChIP-seq data
smad2_chipseq <- read_excel("12915_2014_81_MOESM4_ESM.xlsx", skip = 3)
eomesa_chipseq <- read_excel("12915_2014_81_MOESM10_ESM.xlsx", skip = 2)

# Function to search for a gene in ChIP-seq proximal gene lists
search_chipseq <- function(gene_name, chipseq_df, chip_factor) {
  # Match gene name as whole word in semicolon-separated list
  # Gene names may have ^ suffix (ndr1-upregulated marker in Smad2 data)
  pattern <- paste0("(?i)(^|;)", gene_name, "(\\^)?[;]")
  matches <- chipseq_df %>%
    filter(grepl(pattern, `Proximal gene`, perl = TRUE))

  if (nrow(matches) == 0) {
    return(tibble(
      gene = gene_name,
      chip_factor = chip_factor,
      n_peaks = 0,
      peak_ids = NA_character_,
      has_foxh1_motif = FALSE,
      foxh1_motifs = NA_character_,
      ndr1_upregulated = NA
    ))
  }

  has_foxh1 <- any(!is.na(matches$`Foxh1 binding site`))
  foxh1_seqs <- paste(na.omit(matches$`Foxh1 binding site`), collapse = "; ")
  if (foxh1_seqs == "") foxh1_seqs <- NA_character_

  # Check for ^ marker (only in Smad2 data)
  ndr1_up <- if (chip_factor == "Smad2") {
    any(grepl(paste0("(?i)", gene_name, "\\^"), matches$`Proximal gene`, perl = TRUE))
  } else {
    NA
  }

  tibble(
    gene = gene_name,
    chip_factor = chip_factor,
    n_peaks = nrow(matches),
    peak_ids = paste(matches$`Peak ID`, collapse = "; "),
    has_foxh1_motif = has_foxh1,
    foxh1_motifs = foxh1_seqs,
    ndr1_upregulated = ndr1_up
  )
}

# Search all candidate genes in both ChIP-seq datasets
all_genes_to_check <- unique(c(candidate_genes, nodal_genes))

chipseq_results <- bind_rows(
  map_dfr(all_genes_to_check, ~ search_chipseq(.x, smad2_chipseq, "Smad2")),
  map_dfr(all_genes_to_check, ~ search_chipseq(.x, eomesa_chipseq, "EomesA"))
)

# Create summary table for candidate genes
candidate_chip_summary <- chipseq_results %>%
  filter(gene %in% candidate_genes) %>%
  pivot_wider(
    id_cols = gene,
    names_from = chip_factor,
    values_from = c(n_peaks, has_foxh1_motif, ndr1_upregulated),
    names_sep = "_"
  ) %>%
  mutate(
    Smad2_bound = n_peaks_Smad2 > 0,
    EomesA_bound = n_peaks_EomesA > 0,
    Foxh1_at_Smad2 = has_foxh1_motif_Smad2,
    Foxh1_at_EomesA = has_foxh1_motif_EomesA,
    ndr1_up = ndr1_upregulated_Smad2
  ) %>%
  dplyr::select(gene, Smad2_bound, n_peaks_Smad2, Foxh1_at_Smad2, ndr1_up,
                EomesA_bound, n_peaks_EomesA, Foxh1_at_EomesA) %>%
  arrange(gene)

cat("=== ChIP-seq Binding Summary: Candidate Genes ===\n")
print(as.data.frame(candidate_chip_summary), row.names = FALSE)

# Same for Nodal score genes (for reference)
nodal_chip_summary <- chipseq_results %>%
  filter(gene %in% nodal_genes) %>%
  pivot_wider(
    id_cols = gene,
    names_from = chip_factor,
    values_from = c(n_peaks, has_foxh1_motif, ndr1_upregulated),
    names_sep = "_"
  ) %>%
  mutate(
    Smad2_bound = n_peaks_Smad2 > 0,
    EomesA_bound = n_peaks_EomesA > 0,
    Foxh1_at_Smad2 = has_foxh1_motif_Smad2,
    Foxh1_at_EomesA = has_foxh1_motif_EomesA,
    ndr1_up = ndr1_upregulated_Smad2
  ) %>%
  dplyr::select(gene, Smad2_bound, n_peaks_Smad2, Foxh1_at_Smad2, ndr1_up,
                EomesA_bound, n_peaks_EomesA, Foxh1_at_EomesA) %>%
  arrange(gene)

cat("\n=== ChIP-seq Binding Summary: Nodal Score Genes ===\n")
print(as.data.frame(nodal_chip_summary), row.names = FALSE)

# Save combined summary
write_csv(candidate_chip_summary, "q4_chipseq_candidates.csv")
write_csv(nodal_chip_summary, "q4_chipseq_nodal.csv")

# ============================================================================
# PART 1 FIGURE: ChIP-seq binding heatmap for candidate genes
# ============================================================================

# Prepare data for heatmap-style dot plot
chip_plot_data <- bind_rows(
  candidate_chip_summary %>% mutate(gene_group = "Candidate"),
  nodal_chip_summary %>% mutate(gene_group = "Nodal score")
) %>%
  pivot_longer(
    cols = c(Smad2_bound, EomesA_bound, Foxh1_at_Smad2, Foxh1_at_EomesA),
    names_to = "feature",
    values_to = "present"
  ) %>%
  mutate(
    feature = case_when(
      feature == "Smad2_bound" ~ "Smad2",
      feature == "EomesA_bound" ~ "EomesA",
      feature == "Foxh1_at_Smad2" ~ "Foxh1 motif\n(Smad2 peak)",
      feature == "Foxh1_at_EomesA" ~ "Foxh1 motif\n(EomesA peak)"
    ),
    feature = factor(feature, levels = c("Smad2", "Foxh1 motif\n(Smad2 peak)",
                                         "EomesA", "Foxh1 motif\n(EomesA peak)")),
    gene = factor(gene, levels = rev(unique(c(
      sort(candidate_genes[candidate_genes %in% gene]),
      sort(nodal_genes[nodal_genes %in% gene])
    ))))
  )

p_chip <- ggplot(chip_plot_data, aes(x = feature, y = gene)) +
  geom_tile(aes(fill = present), color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c("TRUE" = "#2166AC", "FALSE" = "#F0F0F0"),
    labels = c("TRUE" = "Bound", "FALSE" = "Not bound"),
    name = ""
  ) +
  facet_grid(gene_group ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = "", y = "",
    title = "ChIP-seq binding evidence",
    subtitle = "Smad2 & EomesA peaks ± Foxh1 motif (Nelson et al. 2014)"
  ) +
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8, face = "italic"),
    strip.text.y = element_text(size = 10, face = "bold", angle = 0),
    legend.position = "bottom",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "#666666")
  )

pdf("q4_chipseq_binding_heatmap.pdf", width = 6, height = 12, family = "Helvetica")
print(p_chip)
dev.off()
cat("\nSaved: q4_chipseq_binding_heatmap.pdf\n")

# ============================================================================
# PART 2: Temporal Dynamics – Candidate Genes vs. Nodal Score
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("PART 2: Temporal Dynamics in Experiment 1\n")
cat(strrep("=", 70), "\n\n")

# Use Experiment 1: 15 ng/ml Activin at different timepoints
# Normalize all Exp1 samples together
exp1_meta <- metadata %>%
  filter(experiment == "Exp1")

exp1_counts <- counts_filtered[, rownames(exp1_meta)]

dds_exp1 <- DESeqDataSetFromMatrix(
  countData = exp1_counts,
  colData = exp1_meta,
  design = ~ 1
)
dds_exp1 <- estimateSizeFactors(dds_exp1)
norm_counts_exp1 <- counts(dds_exp1, normalized = TRUE)

# Focus on 15 ng/ml Activin and 0 ng/ml control across timepoints
exp1_15ngml <- exp1_meta %>% filter(concentration == "15ngml")
exp1_0ngml  <- exp1_meta %>% filter(concentration == "0ngml")

# For each gene: compute log2FC at each timepoint (15 ng/ml vs mean of 0 ng/ml at same time)
timepoints <- sort(unique(exp1_15ngml$time_min))

# Genes available in the count matrix (use unique sets to avoid duplication)
available_candidates <- candidate_genes_unique[candidate_genes_unique %in% tolower(rownames(norm_counts_exp1))]
available_nodal <- nodal_genes[nodal_genes %in% tolower(rownames(norm_counts_exp1))]

cat("Available candidate genes:", length(available_candidates), "of", length(candidate_genes), "\n")
cat("  Missing:", paste(setdiff(candidate_genes, available_candidates), collapse = ", "), "\n")
cat("Available Nodal score genes:", length(available_nodal), "of", length(nodal_genes), "\n")
cat("  Missing:", paste(setdiff(nodal_genes, available_nodal), collapse = ", "), "\n\n")

# Compute temporal expression profiles (normalized counts, log2-transformed)
compute_temporal_profile <- function(genes, norm_mat, meta_15, meta_0, timepoints) {
  profiles <- map_dfr(genes, function(gene) {
    # Case-insensitive row matching
    row_idx <- which(tolower(rownames(norm_mat)) == gene)
    if (length(row_idx) == 0) return(NULL)

    map_dfr(timepoints, function(tp) {
      treated_samples <- rownames(meta_15)[meta_15$time_min == tp]
      control_samples <- rownames(meta_0)[meta_0$time_min == tp]

      if (length(treated_samples) == 0 || length(control_samples) == 0) return(NULL)

      treated_vals <- norm_mat[row_idx, treated_samples]
      control_vals <- norm_mat[row_idx, control_samples]

      # Log2 fold change (pseudocount of 1)
      log2fc <- log2(mean(treated_vals) + 1) - log2(mean(control_vals) + 1)

      tibble(
        gene = gene,
        time_min = tp,
        log2fc = log2fc,
        mean_treated = mean(treated_vals),
        mean_control = mean(control_vals)
      )
    })
  })
  return(profiles)
}

profiles_candidates <- compute_temporal_profile(available_candidates, norm_counts_exp1,
                                                 exp1_15ngml, exp1_0ngml, timepoints)
profiles_nodal <- compute_temporal_profile(available_nodal, norm_counts_exp1,
                                            exp1_15ngml, exp1_0ngml, timepoints)

profiles_candidates$gene_group <- "Candidate"
profiles_nodal$gene_group <- "Nodal score"

all_profiles <- bind_rows(profiles_candidates, profiles_nodal)

# Compute peak expression time for each gene
peak_times <- all_profiles %>%
  group_by(gene, gene_group) %>%
  summarise(
    peak_time = time_min[which.max(abs(log2fc))],
    peak_log2fc = log2fc[which.max(abs(log2fc))],
    direction = ifelse(log2fc[which.max(abs(log2fc))] > 0, "Up", "Down"),
    .groups = "drop"
  )

cat("=== Peak Expression Times ===\n")
cat("\nCandidate genes:\n")
peak_cand <- peak_times %>% filter(gene_group == "Candidate") %>% arrange(peak_time)
print(as.data.frame(peak_cand), row.names = FALSE)

cat("\nNodal score genes:\n")
peak_nodal <- peak_times %>% filter(gene_group == "Nodal score") %>% arrange(peak_time)
print(as.data.frame(peak_nodal), row.names = FALSE)

write_csv(peak_times, "q4_peak_expression_times.csv")

# ============================================================================
# PART 2 FIGURE: Temporal dynamics line plots
# ============================================================================

# Panel A: Faceted small multiples – one mini-panel per gene
# (avoids colour-indistinguishability with 17+ genes)
p_temporal_candidates <- ggplot(profiles_candidates,
       aes(x = time_min, y = log2fc)) +
  geom_area(fill = "#E66101", alpha = 0.15) +
  geom_line(color = "#E66101", linewidth = 0.7) +
  geom_point(color = "#E66101", size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#999999", linewidth = 0.3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  scale_x_continuous(breaks = c(min(timepoints), max(timepoints))) +
  labs(
    x = "Time (min)", y = expression(log[2]~"FC"),
    title = "Candidate genes"
  ) +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7, face = "italic"),
    strip.background = element_rect(fill = "#FFF5EB"),
    plot.title = element_text(size = 11, face = "bold")
  )

p_temporal_nodal <- ggplot(profiles_nodal,
       aes(x = time_min, y = log2fc)) +
  geom_area(fill = "#5E3C99", alpha = 0.15) +
  geom_line(color = "#5E3C99", linewidth = 0.7) +
  geom_point(color = "#5E3C99", size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#999999", linewidth = 0.3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  scale_x_continuous(breaks = c(min(timepoints), max(timepoints))) +
  labs(
    x = "Time (min)", y = expression(log[2]~"FC"),
    title = "Nodal score (direct targets)"
  ) +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7, face = "italic"),
    strip.background = element_rect(fill = "#F0EBF8"),
    plot.title = element_text(size = 11, face = "bold")
  )

# Panel B: Summary – mean |log2FC| trajectory per group with ribbon
group_summary <- all_profiles %>%
  group_by(gene_group, time_min) %>%
  summarise(
    mean_abs_log2fc = mean(abs(log2fc)),
    sd_abs_log2fc = sd(abs(log2fc)),
    mean_log2fc = mean(log2fc),
    sd_log2fc = sd(log2fc),
    .groups = "drop"
  )

group_colors <- c("Candidate" = "#E66101", "Nodal score" = "#5E3C99")

p_temporal_summary <- ggplot(group_summary,
       aes(x = time_min, y = mean_abs_log2fc, color = gene_group, fill = gene_group)) +
  geom_ribbon(aes(ymin = pmax(0, mean_abs_log2fc - sd_abs_log2fc),
                  ymax = mean_abs_log2fc + sd_abs_log2fc),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = group_colors, name = "") +
  scale_fill_manual(values = group_colors, name = "") +
  scale_x_continuous(breaks = timepoints) +
  labs(
    x = "Time (min)", y = expression("Mean |" * log[2] * " FC|"),
    title = "Response magnitude over time"
  ) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 11, face = "bold")
  )

# Panel C: Peak time distribution comparison
p_peak_dist <- ggplot(peak_times, aes(x = factor(peak_time), fill = gene_group)) +
  geom_bar(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = group_colors, name = "") +
  labs(
    x = "Time of peak |log₂FC| (min)", y = "Number of genes",
    title = "Peak response timing"
  ) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 11, face = "bold")
  )

# Save each temporal figure separately (small multiples need more space)
pdf("q4_temporal_candidates.pdf", width = 10, height = 6, family = "Helvetica")
print(p_temporal_candidates)
dev.off()

pdf("q4_temporal_nodal.pdf", width = 10, height = 8, family = "Helvetica")
print(p_temporal_nodal)
dev.off()

# Summary figure: mean response + peak timing
fig_temporal_summary <- (p_temporal_summary + p_peak_dist) +
  plot_annotation(
    title = "Temporal dynamics: Candidate genes vs. Nodal score targets",
    subtitle = "Experiment 1 – 15 ng/ml Activin stimulation across time",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                                family = "Helvetica"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#666666",
                                   family = "Helvetica")
    )
  )

pdf("q4_temporal_dynamics.pdf", width = 12, height = 5, family = "Helvetica")
print(fig_temporal_summary)
dev.off()
cat("\nSaved: q4_temporal_candidates.pdf, q4_temporal_nodal.pdf, q4_temporal_dynamics.pdf\n")

# ============================================================================
# PART 3: Reversibility Analysis (Exp2)
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("PART 3: Reversibility in Experiment 2\n")
cat(strrep("=", 70), "\n\n")

# DESIGN:
#   Reference: 15 ng/ml Activin + DMSO (maximum Activin stimulation)
#   Test conditions: 0 ng/ml DMSO (baseline), SB50 at 60/120/180 min
#   Reversibility = how much the gene returns toward baseline when SB50 is added
#
#   Reversibility score = (Activin - SB50) / (Activin - Baseline)
#     = 1 if fully reversed (gene back to baseline)
#     = 0 if not reversed at all (gene stays at Activin level)
#     < 0 if SB50 amplifies the Activin effect
#     > 1 if SB50 overcorrects past baseline

# Normalize Exp2 samples
exp2_meta <- metadata %>%
  filter(experiment == "Exp2") %>%
  filter(!is.na(concentration))

exp2_counts <- counts_filtered[, rownames(exp2_meta)]

dds_exp2 <- DESeqDataSetFromMatrix(
  countData = exp2_counts,
  colData = exp2_meta,
  design = ~ 1
)
dds_exp2 <- estimateSizeFactors(dds_exp2)
norm_counts_exp2 <- counts(dds_exp2, normalized = TRUE)

# Define conditions
conditions <- list(
  baseline = exp2_meta %>% filter(concentration == "0ngml_DMSO"),
  activin  = exp2_meta %>% filter(concentration == "15ngml_DMSO"),
  sb50_60  = exp2_meta %>% filter(concentration == "SB50", time_min == 60),
  sb50_120 = exp2_meta %>% filter(concentration == "SB50", time_min == 120),
  sb50_180 = exp2_meta %>% filter(concentration == "SB50", time_min == 180)
)

# Compute mean normalized expression per condition per gene
compute_mean_expr <- function(genes, norm_mat, conditions) {
  map_dfr(genes, function(gene) {
    row_idx <- which(tolower(rownames(norm_mat)) == gene)
    if (length(row_idx) == 0) return(NULL)

    map_dfr(names(conditions), function(cond_name) {
      samps <- rownames(conditions[[cond_name]])
      vals <- norm_mat[row_idx, samps]
      tibble(
        gene = gene,
        condition = cond_name,
        mean_expr = mean(vals),
        sd_expr = sd(vals),
        log2_mean = log2(mean(vals) + 1)
      )
    })
  })
}

expr_candidates <- compute_mean_expr(available_candidates, norm_counts_exp2, conditions)
expr_nodal <- compute_mean_expr(available_nodal, norm_counts_exp2, conditions)

# Compute reversibility score for each gene at each SB50 timepoint
compute_reversibility <- function(expr_df) {
  # Pivot to wide format
  expr_wide <- expr_df %>%
    dplyr::select(gene, condition, log2_mean) %>%
    pivot_wider(names_from = condition, values_from = log2_mean)

  # Compute reversibility at each SB50 timepoint
  expr_wide %>%
    mutate(
      # Delta = Activin effect = log2(Activin) - log2(Baseline)
      activin_effect = activin - baseline,
      # Reversibility = fraction of Activin effect reversed by SB50
      rev_60  = ifelse(abs(activin_effect) > 0.1,
                       (activin - sb50_60) / activin_effect, NA_real_),
      rev_120 = ifelse(abs(activin_effect) > 0.1,
                       (activin - sb50_120) / activin_effect, NA_real_),
      rev_180 = ifelse(abs(activin_effect) > 0.1,
                       (activin - sb50_180) / activin_effect, NA_real_)
    ) %>%
    dplyr::select(gene, activin_effect, rev_60, rev_120, rev_180)
}

rev_candidates <- compute_reversibility(expr_candidates) %>% mutate(gene_group = "Candidate")
rev_nodal <- compute_reversibility(expr_nodal) %>% mutate(gene_group = "Nodal score")
rev_all <- bind_rows(rev_candidates, rev_nodal)

cat("=== Reversibility Scores (1 = fully reversed, 0 = not reversed) ===\n")
cat("\nCandidate genes:\n")
rev_cand_print <- rev_candidates %>%
  dplyr::select(gene, activin_effect, rev_60, rev_120, rev_180) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  arrange(desc(abs(activin_effect)))
print(as.data.frame(rev_cand_print), row.names = FALSE)

cat("\nNodal score genes:\n")
rev_nod_print <- rev_nodal %>%
  dplyr::select(gene, activin_effect, rev_60, rev_120, rev_180) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  arrange(desc(abs(activin_effect)))
print(as.data.frame(rev_nod_print), row.names = FALSE)

write_csv(rev_all, "q4_reversibility_scores.csv")

# Statistical comparison of reversibility between groups
cat("\n=== Statistical Comparison: Reversibility (Candidates vs Nodal Score) ===\n")
for (tp in c("rev_60", "rev_120", "rev_180")) {
  tp_label <- gsub("rev_", "SB50 @ ", tp)
  vals_cand <- rev_candidates[[tp]]
  vals_nod <- rev_nodal[[tp]]
  vals_cand <- vals_cand[!is.na(vals_cand)]
  vals_nod <- vals_nod[!is.na(vals_nod)]

  if (length(vals_cand) >= 3 && length(vals_nod) >= 3) {
    wt <- wilcox.test(vals_cand, vals_nod)
    cat(sprintf("  %smin: Candidates median=%.2f, Nodal median=%.2f, Wilcoxon p=%.4f\n",
                tp_label, median(vals_cand), median(vals_nod), wt$p.value))
  }
}

# ============================================================================
# PART 3 FIGURE: Reversibility plots
# ============================================================================

# Panel A: Reversibility per gene (heatmap-style)
rev_long <- rev_all %>%
  pivot_longer(cols = c(rev_60, rev_120, rev_180),
               names_to = "sb50_time", values_to = "reversibility") %>%
  mutate(
    sb50_time = case_when(
      sb50_time == "rev_60" ~ "SB50 @ 60 min",
      sb50_time == "rev_120" ~ "SB50 @ 120 min",
      sb50_time == "rev_180" ~ "SB50 @ 180 min"
    ),
    sb50_time = factor(sb50_time, levels = c("SB50 @ 60 min", "SB50 @ 120 min", "SB50 @ 180 min")),
    gene = factor(gene, levels = rev(
      rev_all %>% arrange(gene_group, desc(abs(activin_effect))) %>% pull(gene) %>% unique()
    ))
  )

p_rev_heatmap <- ggplot(rev_long %>% filter(!is.na(reversibility)),
       aes(x = sb50_time, y = gene, fill = reversibility)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", reversibility)),
            size = 2.2, color = "black") +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0.5, limits = c(-0.5, 1.5),
    oob = scales::squish,
    name = "Reversibility\nscore"
  ) +
  facet_grid(gene_group ~ ., scales = "free_y", space = "free_y") +
  labs(x = "", y = "", title = "Reversibility per gene") +
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 7, face = "italic"),
    strip.text.y = element_text(size = 10, face = "bold", angle = 0),
    plot.title = element_text(size = 11, face = "bold")
  )

# Panel B: Box plot comparison of reversibility between groups
p_rev_box <- ggplot(rev_long %>% filter(!is.na(reversibility)),
       aes(x = sb50_time, y = reversibility, fill = gene_group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1, width = 0.6,
               position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = c(0, 1), linetype = c("dashed", "solid"),
             color = c("#B2182B", "#2166AC"), linewidth = 0.4) +
  scale_fill_manual(values = group_colors, name = "") +
  labs(
    x = "", y = "Reversibility score",
    title = "Reversibility comparison",
    caption = "Score = 1: fully reversed | Score = 0: not reversed"
  ) +
  coord_cartesian(ylim = c(-1, 2)) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 11, face = "bold"),
    plot.caption = element_text(size = 8, hjust = 0, color = "#666666")
  )

# Panel C: Expression profiles across Exp2 conditions (selected genes)
expr_all <- bind_rows(
  expr_candidates %>% mutate(gene_group = "Candidate"),
  expr_nodal %>% mutate(gene_group = "Nodal score")
) %>%
  mutate(
    condition_label = case_when(
      condition == "baseline" ~ "0 ng/ml\n(Baseline)",
      condition == "activin" ~ "15 ng/ml\n(Activin)",
      condition == "sb50_60" ~ "SB50\n60 min",
      condition == "sb50_120" ~ "SB50\n120 min",
      condition == "sb50_180" ~ "SB50\n180 min"
    ),
    condition_label = factor(condition_label,
                             levels = c("0 ng/ml\n(Baseline)", "15 ng/ml\n(Activin)",
                                        "SB50\n60 min", "SB50\n120 min", "SB50\n180 min"))
  )

# Normalize expression relative to baseline for comparability
# Group by gene AND gene_group to handle genes present in both lists
expr_norm <- expr_all %>%
  group_by(gene, gene_group) %>%
  mutate(
    baseline_expr = log2_mean[condition == "baseline"][1],
    delta_log2 = log2_mean - baseline_expr
  ) %>%
  ungroup()

# Mean delta per group
group_expr_summary <- expr_norm %>%
  group_by(gene_group, condition_label) %>%
  summarise(
    mean_delta = mean(delta_log2),
    sd_delta = sd(delta_log2),
    se_delta = sd(delta_log2) / sqrt(n()),
    .groups = "drop"
  )

p_expr_profile <- ggplot(group_expr_summary,
       aes(x = condition_label, y = mean_delta, color = gene_group, group = gene_group)) +
  geom_ribbon(aes(ymin = mean_delta - se_delta, ymax = mean_delta + se_delta,
                  fill = gene_group),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#999999", linewidth = 0.3) +
  scale_color_manual(values = group_colors, name = "") +
  scale_fill_manual(values = group_colors, name = "") +
  labs(
    x = "", y = expression(Delta ~ log[2] ~ "(norm. counts) vs baseline"),
    title = "Mean expression shift across Exp2 conditions",
    caption = "Shading: ± SE"
  ) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 11, face = "bold"),
    plot.caption = element_text(size = 8, hjust = 0, color = "#666666")
  )

# Combine reversibility figure
fig_reversibility <- (p_rev_heatmap | (p_rev_box / p_expr_profile)) +
  plot_annotation(
    title = "Reversibility of Activin-induced expression by SB50",
    subtitle = "Experiment 2 – Candidate genes vs. direct Nodal targets",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                                family = "Helvetica"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#666666",
                                   family = "Helvetica")
    )
  ) +
  plot_layout(widths = c(1.2, 1))

pdf("q4_reversibility.pdf", width = 14, height = 12, family = "Helvetica")
print(fig_reversibility)
dev.off()
cat("\nSaved: q4_reversibility.pdf\n")

# ============================================================================
# COMBINED OVERVIEW FIGURE: ChIP-seq + Dynamics + Reversibility
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("Creating combined overview figure\n")
cat(strrep("=", 70), "\n\n")

# Summary dot plot: combine ChIP-seq binding, peak time, and reversibility
combined_summary <- peak_times %>%
  left_join(
    rev_all %>% dplyr::select(gene, gene_group, activin_effect, rev_60),
    by = c("gene", "gene_group")
  ) %>%
  left_join(
    chipseq_results %>%
      dplyr::select(gene, chip_factor, n_peaks, has_foxh1_motif) %>%
      pivot_wider(
        id_cols = gene,
        names_from = chip_factor,
        values_from = c(n_peaks, has_foxh1_motif),
        names_sep = "_"
      ),
    by = "gene"
  ) %>%
  mutate(
    any_chipseq = (n_peaks_Smad2 > 0) | (n_peaks_EomesA > 0),
    both_chipseq = (n_peaks_Smad2 > 0) & (n_peaks_EomesA > 0),
    chip_label = case_when(
      n_peaks_Smad2 > 0 & n_peaks_EomesA > 0 ~ "Smad2 + EomesA",
      n_peaks_Smad2 > 0 ~ "Smad2 only",
      n_peaks_EomesA > 0 ~ "EomesA only",
      TRUE ~ "No binding"
    )
  )

chip_fill_colors <- c("Smad2 + EomesA" = "#7570B3", "Smad2 only" = "#D95F02",
                       "EomesA only" = "#1B9E77", "No binding" = "#E0E0E0")

p_overview <- ggplot(combined_summary,
       aes(x = peak_time, y = rev_60, size = abs(peak_log2fc))) +
  geom_point(aes(fill = chip_label, shape = gene_group),
             alpha = 0.8, stroke = 0.5, color = "black") +
  geom_text_repel(
    aes(label = gene),
    size = 2.8, fontface = "italic", color = "#333333",
    max.overlaps = 25, seed = 42, box.padding = 0.4,
    min.segment.length = 0.2
  ) +
  geom_hline(yintercept = c(0, 1), linetype = c("dashed", "solid"),
             color = c("#B2182B", "#2166AC"), linewidth = 0.3) +
  scale_shape_manual(values = c("Candidate" = 21, "Nodal score" = 24), name = "Gene group") +
  scale_fill_manual(values = chip_fill_colors, name = "ChIP-seq binding") +
  scale_size_continuous(range = c(2, 8), name = "|Peak log₂FC|") +
  scale_x_continuous(breaks = timepoints) +
  coord_cartesian(ylim = c(-0.5, 1.8)) +
  labs(
    x = "Time of peak response in Exp1 (min)",
    y = "Reversibility at SB50 60 min",
    title = "Integrated view: ChIP-seq binding, response timing, and reversibility",
    subtitle = "Direct Nodal targets tend to respond earlier and are more reversible",
    caption = paste0(
      "Reversibility = 1: fully blocked by SB50 | = 0: not blocked\n",
      "ChIP-seq data from Nelson et al. (2014) BMC Genomics"
    )
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 4)),
    shape = guide_legend(override.aes = list(size = 4))
  ) +
  theme_bw(base_size = 11, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "#666666"),
    plot.caption = element_text(size = 8, hjust = 0, color = "#666666")
  )

pdf("q4_integrated_overview.pdf", width = 12, height = 8, family = "Helvetica")
print(p_overview)
dev.off()
cat("\nSaved: q4_integrated_overview.pdf\n")

write_csv(combined_summary, "q4_combined_summary.csv")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n\n")

cat("ChIP-seq binding among candidate genes:\n")
cat(sprintf("  Smad2 peaks: %d / %d genes have nearby peaks\n",
            sum(candidate_chip_summary$Smad2_bound), nrow(candidate_chip_summary)))
cat(sprintf("  EomesA peaks: %d / %d genes have nearby peaks\n",
            sum(candidate_chip_summary$EomesA_bound), nrow(candidate_chip_summary)))
cat(sprintf("  Foxh1 motif at Smad2 peaks: %d\n",
            sum(candidate_chip_summary$Foxh1_at_Smad2, na.rm = TRUE)))
cat(sprintf("  Foxh1 motif at EomesA peaks: %d\n",
            sum(candidate_chip_summary$Foxh1_at_EomesA, na.rm = TRUE)))

cat("\nChIP-seq binding among Nodal score genes:\n")
cat(sprintf("  Smad2 peaks: %d / %d genes have nearby peaks\n",
            sum(nodal_chip_summary$Smad2_bound), nrow(nodal_chip_summary)))
cat(sprintf("  EomesA peaks: %d / %d genes have nearby peaks\n",
            sum(nodal_chip_summary$EomesA_bound), nrow(nodal_chip_summary)))

cat("\nTemporal dynamics:\n")
cat(sprintf("  Candidates: median peak at %d min\n",
            median(peak_cand$peak_time)))
cat(sprintf("  Nodal score: median peak at %d min\n",
            median(peak_nodal$peak_time)))

rev_60_cand <- rev_candidates$rev_60[!is.na(rev_candidates$rev_60)]
rev_60_nod <- rev_nodal$rev_60[!is.na(rev_nodal$rev_60)]
cat("\nReversibility (SB50 @ 60 min):\n")
cat(sprintf("  Candidates: median = %.2f (range: %.2f – %.2f)\n",
            median(rev_60_cand), min(rev_60_cand), max(rev_60_cand)))
cat(sprintf("  Nodal score: median = %.2f (range: %.2f – %.2f)\n",
            median(rev_60_nod), min(rev_60_nod), max(rev_60_nod)))

cat("\n=== Output Files ===\n")
cat("  q4_chipseq_candidates.csv      – ChIP-seq binding for candidate genes\n")
cat("  q4_chipseq_nodal.csv           – ChIP-seq binding for Nodal score genes\n")
cat("  q4_peak_expression_times.csv   – Peak expression timing per gene\n")
cat("  q4_reversibility_scores.csv    – Reversibility scores per gene\n")
cat("  q4_combined_summary.csv        – Integrated summary table\n")
cat("  q4_chipseq_binding_heatmap.pdf – ChIP-seq binding heatmap\n")
cat("  q4_temporal_candidates.pdf     – Candidate gene temporal facets\n")
cat("  q4_temporal_nodal.pdf          – Nodal score gene temporal facets\n")
cat("  q4_temporal_dynamics.pdf       – Summary temporal dynamics figure\n")
cat("  q4_reversibility.pdf           – Reversibility figure\n")
cat("  q4_integrated_overview.pdf     – Combined overview figure\n")
