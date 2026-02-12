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
# OUTPUT FIGURES:
#   1. ChIP-seq binding heatmap (candidates + Nodal score)
#   2. Temporal behaviour clusters (2×2 facet: 4 response patterns)
#   3. Reversibility expression profile (Exp2 mean expression per condition)
#   4. Integrated overview scatter (split: candidates / Nodal score)
#   5. Candidate validation volcano
#

source("preprocess.R")
library(ggrepel)

# ============================================================================
# DATA LOADED FROM preprocess.R: counts_filtered, metadata, results_path()
# ============================================================================

cat("Total samples:", ncol(counts_filtered), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")
cat("\nConcentration breakdown:\n")
print(table(metadata$concentration, metadata$experiment))

# ============================================================================
# DEFINE GENE LISTS
# ============================================================================
#
#   NODAL SCORE (27 genes): Established direct Nodal/Activin targets from
#     the literature (e.g. gsc, lft1, ndr1, sox32).
#
#   CANDIDATE LIST (18 genes): Genes identified as Activin-responsive in
#     Lena's experiments, but whose direct/indirect status is unknown.
#
#   OVERLAP: dkk1b, kirrel3l, flrt3, efnb2a appear in both. They are
#     assigned to the Nodal score group to avoid duplication.

candidate_genes <- c("rhov", "net1", "flrt3", "rnd1b", "dkk1b", "plekha5b",
                     "rasgef1ba", "efna1a", "osr1", "frmd4ba", "abi1b",
                     "snai1b", "snai1a", "kirrel3l", "jcada", "pfkfb3",
                     "prickle1b", "efnb2a")

nodal_score_data <- read_excel("docs/nodal-score-genes_complete.xlsx", skip = 1)
nodal_genes <- unique(tolower(na.omit(nodal_score_data$`Nodal score`)))

overlapping_genes <- intersect(candidate_genes, nodal_genes)
candidate_genes_unique <- setdiff(candidate_genes, nodal_genes)

cat("Candidate genes:", length(candidate_genes), "(", length(candidate_genes_unique), "unique)\n")
cat("Nodal score genes:", length(nodal_genes), "\n")
cat("Overlapping:", paste(overlapping_genes, collapse = ", "), "\n")

# Colour palette (consistent throughout)
group_colors <- c("Candidate" = "#E66101", "Nodal score" = "#5E3C99")

# ============================================================================
# PART 1: ChIP-seq Binding Summary + Heatmap
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("PART 1: ChIP-seq Binding Summary\n")
cat(strrep("=", 70), "\n\n")

smad2_chipseq <- read_excel("12915_2014_81_MOESM4_ESM.xlsx", skip = 3)
eomesa_chipseq <- read_excel("12915_2014_81_MOESM10_ESM.xlsx", skip = 2)

search_chipseq <- function(gene_name, chipseq_df, chip_factor) {
  pattern <- paste0("(?i)(^|;)", gene_name, "(\\^)?[;]")
  matches <- chipseq_df %>%
    filter(grepl(pattern, `Proximal gene`, perl = TRUE))

  if (nrow(matches) == 0) {
    return(tibble(gene = gene_name, chip_factor = chip_factor,
                  n_peaks = 0, has_foxh1_motif = FALSE, ndr1_upregulated = NA))
  }

  tibble(
    gene = gene_name,
    chip_factor = chip_factor,
    n_peaks = nrow(matches),
    has_foxh1_motif = any(!is.na(matches$`Foxh1 binding site`)),
    ndr1_upregulated = if (chip_factor == "Smad2") {
      any(grepl(paste0("(?i)", gene_name, "\\^"), matches$`Proximal gene`, perl = TRUE))
    } else NA
  )
}

all_genes_to_check <- unique(c(candidate_genes, nodal_genes))

chipseq_results <- bind_rows(
  map_dfr(all_genes_to_check, ~ search_chipseq(.x, smad2_chipseq, "Smad2")),
  map_dfr(all_genes_to_check, ~ search_chipseq(.x, eomesa_chipseq, "EomesA"))
)

# Build ChIP-seq summary tables
make_chip_summary <- function(genes) {
  chipseq_results %>%
    filter(gene %in% genes) %>%
    pivot_wider(id_cols = gene, names_from = chip_factor,
                values_from = c(n_peaks, has_foxh1_motif, ndr1_upregulated),
                names_sep = "_") %>%
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
}

candidate_chip_summary <- make_chip_summary(candidate_genes)
nodal_chip_summary <- make_chip_summary(nodal_genes)

cat("=== Candidate Genes ===\n")
print(as.data.frame(candidate_chip_summary), row.names = FALSE)
cat("\n=== Nodal Score Genes ===\n")
print(as.data.frame(nodal_chip_summary), row.names = FALSE)

write_csv(candidate_chip_summary, results_path("q4_chipseq_candidates.csv"))
write_csv(nodal_chip_summary, results_path("q4_chipseq_nodal.csv"))

# --- FIGURE 1: ChIP-seq binding heatmap ---

chip_plot_data <- bind_rows(
  candidate_chip_summary %>% mutate(gene_group = "Candidate"),
  nodal_chip_summary %>% mutate(gene_group = "Nodal score")
) %>%
  pivot_longer(
    cols = c(Smad2_bound, EomesA_bound, Foxh1_at_Smad2, Foxh1_at_EomesA),
    names_to = "feature", values_to = "present"
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
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8, face = "italic"),
    strip.text.y = element_text(size = 10, face = "bold", angle = 0),
    legend.position = "bottom",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey50")
  )

ggsave(results_path("q4_chipseq_binding_heatmap.pdf"), p_chip, width = 6, height = 12)
cat("\nSaved:", results_path("q4_chipseq_binding_heatmap.pdf"), "\n")

# ============================================================================
# PART 2: Temporal Dynamics & Behaviour Clustering
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("PART 2: Temporal Dynamics in Experiment 1\n")
cat(strrep("=", 70), "\n\n")

# Normalize Exp1 samples
exp1_meta <- metadata %>% filter(experiment == "Exp1")
exp1_counts <- counts_filtered[, rownames(exp1_meta)]

dds_exp1 <- DESeqDataSetFromMatrix(
  countData = exp1_counts, colData = exp1_meta, design = ~ 1
)
dds_exp1 <- estimateSizeFactors(dds_exp1)
norm_counts_exp1 <- counts(dds_exp1, normalized = TRUE)

exp1_15ngml <- exp1_meta %>% filter(concentration == "15ngml")
exp1_0ngml  <- exp1_meta %>% filter(concentration == "0ngml")
timepoints <- sort(unique(exp1_15ngml$time_min))

available_candidates <- candidate_genes_unique[candidate_genes_unique %in% tolower(rownames(norm_counts_exp1))]
available_nodal <- nodal_genes[nodal_genes %in% tolower(rownames(norm_counts_exp1))]

cat("Available candidates:", length(available_candidates), "of", length(candidate_genes_unique), "\n")
cat("Available Nodal genes:", length(available_nodal), "of", length(nodal_genes), "\n\n")

# Compute log2FC at each timepoint
compute_temporal_profile <- function(genes, norm_mat, meta_15, meta_0, tp) {
  map_dfr(genes, function(gene) {
    row_idx <- which(tolower(rownames(norm_mat)) == gene)
    if (length(row_idx) == 0) return(NULL)
    map_dfr(tp, function(t) {
      treated <- norm_mat[row_idx, rownames(meta_15)[meta_15$time_min == t]]
      control <- norm_mat[row_idx, rownames(meta_0)[meta_0$time_min == t]]
      if (length(treated) == 0 || length(control) == 0) return(NULL)
      tibble(gene = gene, time_min = t,
             log2fc = log2(mean(treated) + 1) - log2(mean(control) + 1))
    })
  })
}

profiles_candidates <- compute_temporal_profile(available_candidates, norm_counts_exp1,
                                                 exp1_15ngml, exp1_0ngml, timepoints) %>%
  mutate(gene_group = "Candidate")
profiles_nodal <- compute_temporal_profile(available_nodal, norm_counts_exp1,
                                            exp1_15ngml, exp1_0ngml, timepoints) %>%
  mutate(gene_group = "Nodal score")
all_profiles <- bind_rows(profiles_candidates, profiles_nodal)

# Peak times
peak_times <- all_profiles %>%
  group_by(gene, gene_group) %>%
  summarise(
    peak_time = time_min[which.max(abs(log2fc))],
    peak_log2fc = log2fc[which.max(abs(log2fc))],
    direction = ifelse(log2fc[which.max(abs(log2fc))] > 0, "Up", "Down"),
    .groups = "drop"
  )
write_csv(peak_times, results_path("q4_peak_expression_times.csv"))

# --- Cluster genes into 4 categories by temporal profile ---
#
# Four behaviour classes:
#   1. Down-regulated: Mean log2FC is negative across most timepoints
#   2. Time-exposure sensitive: Peak response at early/mid timepoints, declines later
#   3. Sustained activation: Steady increase or plateau at high log2FC
#   4. Other: Mixed or weak patterns
#

classify_temporal <- function(gene_name, gene_df) {
  fc <- gene_df$log2fc
  tp <- gene_df$time_min
  n_tp <- length(fc)
  mean_fc <- mean(fc)
  max_fc <- max(fc)
  min_fc <- min(fc)
  last_fc <- fc[n_tp]
  peak_idx <- which.max(abs(fc))
  peak_fc <- fc[peak_idx]

  # Monotonic increase check
  diffs <- diff(fc)
  mostly_up <- sum(diffs > -0.1) >= (n_tp - 2)

  if (mean_fc < -0.2 & sum(fc < 0) >= ceiling(n_tp / 2)) {
    return("Down-regulated")
  } else if (peak_fc > 0.3 & peak_idx < n_tp & last_fc < 0.7 * max_fc) {
    return("Time-exposure sensitive")
  } else if (mean_fc > 0.3 & mostly_up) {
    return("Sustained activation")
  } else {
    return("Other")
  }
}

gene_clusters <- all_profiles %>%
  group_by(gene, gene_group) %>%
  group_modify(~ tibble(cluster = classify_temporal(.y$gene, .x))) %>%
  ungroup()

gene_clusters$cluster <- factor(gene_clusters$cluster,
  levels = c("Down-regulated", "Time-exposure sensitive", "Sustained activation", "Other"))

cat("=== Temporal Behaviour Clusters ===\n")
for (cl in levels(gene_clusters$cluster)) {
  cl_genes <- gene_clusters %>% filter(cluster == cl)
  cat(sprintf("  %s (%d): %s\n", cl, nrow(cl_genes),
              paste(cl_genes$gene, collapse = ", ")))
}

write_csv(gene_clusters, results_path("q4_temporal_clusters.csv"))

# --- FIGURE 2: 2×2 cluster facet plot with gene labels ---

cluster_plot_data <- all_profiles %>%
  inner_join(gene_clusters, by = c("gene", "gene_group"))

# Build facet labels with gene count
cluster_gene_labels <- gene_clusters %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    gene_list = paste(gene, collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(facet_label = paste0(cluster, " (n = ", n, ")"))

cluster_plot_data <- cluster_plot_data %>%
  left_join(cluster_gene_labels %>% dplyr::select(cluster, facet_label), by = "cluster")

# Cluster mean per group (thick line)
cluster_ribbon <- cluster_plot_data %>%
  group_by(facet_label, gene_group, time_min) %>%
  summarise(mean_fc = mean(log2fc), sd_fc = sd(log2fc), .groups = "drop")

# Prepare gene-name annotations to display inside each panel
gene_annotations <- cluster_plot_data %>%
  distinct(gene, gene_group, facet_label) %>%
  group_by(facet_label, gene_group) %>%
  summarise(genes_text = paste(sort(gene), collapse = ", "), .groups = "drop")

p_clusters <- ggplot(cluster_plot_data, aes(x = time_min, y = log2fc)) +
  # Individual gene traces (thin, transparent)
  geom_line(aes(group = gene, color = gene_group), alpha = 0.35, linewidth = 0.4) +
  # Cluster mean per group
  geom_line(data = cluster_ribbon,
            aes(y = mean_fc, color = gene_group, group = gene_group),
            linewidth = 1.2) +
  geom_point(data = cluster_ribbon,
             aes(y = mean_fc, color = gene_group), size = 2) +
  # Gene names as text annotation at bottom of each panel
  geom_text(data = gene_annotations,
            aes(x = mean(timepoints), y = -Inf, label = genes_text, color = gene_group),
            vjust = -0.3, size = 2.2, fontface = "italic",
            position = position_dodge(width = 40), show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 2) +
  scale_color_manual(values = group_colors, name = NULL) +
  scale_x_continuous(breaks = timepoints) +
  labs(
    x = "Time (min)", y = "log2FC (15 vs 0 ng/ml Activin)",
    title = "Temporal behaviour clusters",
    subtitle = "Experiment 1 – gene expression response over time at 15 ng/ml Activin"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50")
  )

ggsave(results_path("q4_temporal_clusters.pdf"), p_clusters, width = 10, height = 8)
cat("\nSaved:", results_path("q4_temporal_clusters.pdf"), "\n")

# ============================================================================
# PART 3: Reversibility (Exp2) – Expression Profile
# ============================================================================
#
# REVERSIBILITY SCORE CALCULATION
# ================================
# All Exp2 samples are collected at 240 min.
#
#   Baseline = 0 ng/ml Activin + DMSO    (no stimulation)
#   Activin  = 15 ng/ml Activin + DMSO   (max stimulation, no inhibitor)
#   SB50     = 15 ng/ml Activin + 50 µM SB-431542 added at 60 / 120 / 180 min
#
# For each gene:
#
#   activin_effect = log2(Activin + 1) − log2(Baseline + 1)
#     → how much Activin changed the gene from baseline
#
#   reversal = log2(Activin + 1) − log2(SB50 + 1)
#     → how much SB50 brought expression back toward baseline
#
#   reversibility = reversal / activin_effect
#     → fraction of the Activin effect that is reversed
#
#     = 1.0 → SB50 *fully* reversed the effect (expression returns to Baseline)
#     = 0.0 → SB50 had *no* effect (expression stays at Activin level)
#     = between 0 and 1 → partial reversal
#
# Genes with |activin_effect| < 0.1 log2 units are excluded — the Activin
# effect is too small to meaningfully assess reversibility.
#

cat("\n")
cat(strrep("=", 70), "\n")
cat("PART 3: Reversibility in Experiment 2\n")
cat(strrep("=", 70), "\n\n")

exp2_meta <- metadata %>%
  filter(experiment == "Exp2") %>%
  filter(concentration != "other")

cat("Exp2 samples:", nrow(exp2_meta), "\n")
cat("  Baseline (0ngml_DMSO):", sum(exp2_meta$concentration == "0ngml_DMSO"), "\n")
cat("  Activin (15ngml_DMSO):", sum(exp2_meta$concentration == "15ngml_DMSO"), "\n")
cat("  SB50:", sum(exp2_meta$concentration == "SB50"), "\n\n")

exp2_counts <- counts_filtered[, rownames(exp2_meta)]

dds_exp2 <- DESeqDataSetFromMatrix(
  countData = exp2_counts, colData = exp2_meta, design = ~ 1
)
dds_exp2 <- estimateSizeFactors(dds_exp2)
norm_counts_exp2 <- counts(dds_exp2, normalized = TRUE)

conditions <- list(
  baseline = exp2_meta %>% filter(concentration == "0ngml_DMSO"),
  activin  = exp2_meta %>% filter(concentration == "15ngml_DMSO"),
  sb50_60  = exp2_meta %>% filter(concentration == "SB50", time_min == 60),
  sb50_120 = exp2_meta %>% filter(concentration == "SB50", time_min == 120),
  sb50_180 = exp2_meta %>% filter(concentration == "SB50", time_min == 180)
)

compute_mean_expr <- function(genes, norm_mat, conds) {
  map_dfr(genes, function(gene) {
    row_idx <- which(tolower(rownames(norm_mat)) == gene)
    if (length(row_idx) == 0) return(NULL)
    map_dfr(names(conds), function(cn) {
      vals <- norm_mat[row_idx, rownames(conds[[cn]])]
      tibble(gene = gene, condition = cn,
             mean_expr = mean(vals), sd_expr = sd(vals),
             log2_mean = log2(mean(vals) + 1))
    })
  })
}

expr_candidates <- compute_mean_expr(available_candidates, norm_counts_exp2, conditions)
expr_nodal <- compute_mean_expr(available_nodal, norm_counts_exp2, conditions)

# Reversibility scores
compute_reversibility <- function(expr_df) {
  expr_df %>%
    dplyr::select(gene, condition, log2_mean) %>%
    pivot_wider(names_from = condition, values_from = log2_mean) %>%
    mutate(
      activin_effect = activin - baseline,
      rev_60  = ifelse(abs(activin_effect) > 0.1, (activin - sb50_60) / activin_effect, NA_real_),
      rev_120 = ifelse(abs(activin_effect) > 0.1, (activin - sb50_120) / activin_effect, NA_real_),
      rev_180 = ifelse(abs(activin_effect) > 0.1, (activin - sb50_180) / activin_effect, NA_real_)
    ) %>%
    dplyr::select(gene, activin_effect, rev_60, rev_120, rev_180)
}

rev_candidates <- compute_reversibility(expr_candidates) %>% mutate(gene_group = "Candidate")
rev_nodal <- compute_reversibility(expr_nodal) %>% mutate(gene_group = "Nodal score")
rev_all <- bind_rows(rev_candidates, rev_nodal)

cat("=== Reversibility Scores ===\n\nCandidate genes:\n")
print(as.data.frame(rev_candidates %>%
  dplyr::select(gene, activin_effect, rev_60, rev_120, rev_180) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  arrange(desc(abs(activin_effect)))), row.names = FALSE)

cat("\nNodal score genes:\n")
print(as.data.frame(rev_nodal %>%
  dplyr::select(gene, activin_effect, rev_60, rev_120, rev_180) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  arrange(desc(abs(activin_effect)))), row.names = FALSE)

write_csv(rev_all, results_path("q4_reversibility_scores.csv"))

# --- FIGURE 3: Mean expression profile across Exp2 conditions ---

expr_all <- bind_rows(
  expr_candidates %>% mutate(gene_group = "Candidate"),
  expr_nodal %>% mutate(gene_group = "Nodal score")
) %>%
  mutate(
    condition_label = case_when(
      condition == "baseline" ~ "Baseline\n(0 ng/ml)",
      condition == "activin"  ~ "Activin\n(15 ng/ml)",
      condition == "sb50_60"  ~ "SB50\n@ 60 min",
      condition == "sb50_120" ~ "SB50\n@ 120 min",
      condition == "sb50_180" ~ "SB50\n@ 180 min"
    ),
    condition_label = factor(condition_label,
      levels = c("Baseline\n(0 ng/ml)", "Activin\n(15 ng/ml)",
                 "SB50\n@ 60 min", "SB50\n@ 120 min", "SB50\n@ 180 min"))
  )

expr_norm <- expr_all %>%
  group_by(gene, gene_group) %>%
  mutate(
    baseline_expr = log2_mean[condition == "baseline"][1],
    delta_log2 = log2_mean - baseline_expr
  ) %>%
  ungroup()

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
                  fill = gene_group), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = group_colors, name = NULL) +
  scale_fill_manual(values = group_colors, name = NULL) +
  labs(
    x = "", y = expression(Delta ~ "log2(norm. counts) vs baseline"),
    title = "Reversibility: mean expression shift across Exp2 conditions",
    subtitle = "If SB50 reverses Activin effect, expression should return toward baseline (zero line)",
    caption = paste0(
      "Reversibility score = (log2 Activin − log2 SB50) / (log2 Activin − log2 Baseline)\n",
      "Score 1 = fully reversed | Score 0 = not reversed | Shading = ± SE | All samples at 240 min"
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50"),
    plot.caption = element_text(size = 8, hjust = 0, color = "grey50")
  )

ggsave(results_path("q4_reversibility_profile.pdf"), p_expr_profile, width = 8, height = 6)
cat("\nSaved:", results_path("q4_reversibility_profile.pdf"), "\n")

# ============================================================================
# PART 4: Integrated Overview (split by gene group)
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("PART 4: Integrated Overview\n")
cat(strrep("=", 70), "\n\n")

# Combine peak timing, reversibility at 60 min, and ChIP-seq info
combined_summary <- peak_times %>%
  left_join(
    rev_all %>% dplyr::select(gene, gene_group, activin_effect, rev_60),
    by = c("gene", "gene_group")
  ) %>%
  left_join(
    chipseq_results %>%
      dplyr::select(gene, chip_factor, n_peaks, has_foxh1_motif) %>%
      pivot_wider(id_cols = gene, names_from = chip_factor,
                  values_from = c(n_peaks, has_foxh1_motif), names_sep = "_"),
    by = "gene"
  ) %>%
  mutate(
    chip_label = case_when(
      n_peaks_Smad2 > 0 & n_peaks_EomesA > 0 ~ "Smad2 + EomesA",
      n_peaks_Smad2 > 0 ~ "Smad2 only",
      n_peaks_EomesA > 0 ~ "EomesA only",
      TRUE ~ "No binding"
    ),
    # Clamp reversibility to [0, 1] for display
    rev_60_clamped = pmin(pmax(rev_60, 0), 1)
  )

write_csv(combined_summary, results_path("q4_combined_summary.csv"))

chip_fill_colors <- c("Smad2 + EomesA" = "#7570B3", "Smad2 only" = "#D95F02",
                       "EomesA only" = "#1B9E77", "No binding" = "#E0E0E0")

# --- FIGURE 4: Two separate overview panels ---
# Following q3 styling conventions:
#   theme_minimal(base_size = 10), ggsave(), consistent text sizing

make_overview_panel <- function(data, title_text, title_color) {
  ggplot(data %>% filter(!is.na(rev_60_clamped)),
         aes(x = peak_time, y = rev_60_clamped, size = abs(peak_log2fc))) +
    geom_point(aes(fill = chip_label), shape = 21,
               alpha = 0.85, stroke = 0.5, color = "black") +
    geom_text_repel(
      aes(label = gene),
      size = 3.2, fontface = "italic", color = "#333333",
      max.overlaps = 25, seed = 42, box.padding = 0.5,
      min.segment.length = 0.2
    ) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    scale_fill_manual(values = chip_fill_colors, name = "ChIP-seq binding",
                      drop = FALSE) +
    scale_size_continuous(range = c(2, 8), name = "|Peak log\u2082FC|") +
    scale_x_continuous(breaks = timepoints) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = c("0\n(not reversed)", "0.25", "0.5", "0.75", "1\n(fully reversed)")) +
    labs(x = "Time of peak response (min)", y = "Reversibility (SB50 @ 60 min)",
         title = title_text) +
    guides(fill = guide_legend(override.aes = list(size = 4))) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 9),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
                                color = title_color)
    )
}

p_overview_cand <- make_overview_panel(
  combined_summary %>% filter(gene_group == "Candidate"),
  "Candidate genes", "#E66101"
)
p_overview_nodal <- make_overview_panel(
  combined_summary %>% filter(gene_group == "Nodal score"),
  "Nodal score genes (direct targets)", "#5E3C99"
)

fig_overview <- p_overview_cand / p_overview_nodal +
  plot_annotation(
    title = "Integrated view: response timing vs reversibility",
    subtitle = "Bubble size = peak fold-change magnitude | Fill = ChIP-seq binding evidence",
    caption = "Y = 1: fully reversed by SB50 | Y = 0: not reversed | Data: Nelson et al. (2014) BMC Genomics",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50"),
      plot.caption = element_text(size = 8, hjust = 0, color = "grey50")
    )
  )

ggsave(results_path("q4_integrated_overview.pdf"), fig_overview, width = 10, height = 14)
cat("Saved:", results_path("q4_integrated_overview.pdf"), "\n")

# ============================================================================
# PART 5: Candidate Gene Validation & Alternatives
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("PART 5: Candidate Gene Validation\n")
cat(strrep("=", 70), "\n\n")

de_exp1 <- read_csv(results_path("q1_exp1_0vs15_activin_240min.csv"), show_col_types = FALSE)
de_exp2 <- read_csv(results_path("q3_Activin_vs_Baseline.csv"), show_col_types = FALSE)

validate_candidate <- function(gene_name, de1, de2, chip_res) {
  r1 <- de1 %>% filter(gene == gene_name)
  r2 <- de2 %>% filter(gene == gene_name)

  chip <- chip_res %>% filter(gene == gene_name) %>%
    group_by(gene) %>%
    summarise(sm = sum(n_peaks[chip_factor == "Smad2"]),
              eo = sum(n_peaks[chip_factor == "EomesA"]),
              fx = any(has_foxh1_motif, na.rm = TRUE), .groups = "drop")

  tibble(
    gene = gene_name,
    in_count_matrix = gene_name %in% tolower(rownames(counts_filtered)),
    exp1_log2fc = if (nrow(r1) > 0) r1$log2FoldChange[1] else NA_real_,
    exp1_padj   = if (nrow(r1) > 0) r1$padj[1] else NA_real_,
    exp1_sig = if (nrow(r1) > 0) r1$padj[1] < 0.05 & abs(r1$log2FoldChange[1]) >= 1.0 else FALSE,
    exp2_log2fc = if (nrow(r2) > 0) r2$log2FoldChange[1] else NA_real_,
    exp2_padj   = if (nrow(r2) > 0) r2$padj[1] else NA_real_,
    exp2_sig = if (nrow(r2) > 0) r2$padj[1] < 0.05 & abs(r2$log2FoldChange[1]) >= 1.0 else FALSE,
    smad2_bound  = if (nrow(chip) > 0) chip$sm[1] > 0 else FALSE,
    eomesa_bound = if (nrow(chip) > 0) chip$eo[1] > 0 else FALSE,
    foxh1_motif  = if (nrow(chip) > 0) chip$fx[1] else FALSE,
    de_consistent = ifelse(!is.na(exp1_sig) & !is.na(exp2_sig), exp1_sig & exp2_sig, FALSE),
    has_chipseq = smad2_bound | eomesa_bound,
    quality = case_when(
      !in_count_matrix ~ "NOT IN DATA",
      de_consistent & has_chipseq ~ "STRONG",
      de_consistent & !has_chipseq ~ "DE only",
      !de_consistent & has_chipseq ~ "ChIP only",
      TRUE ~ "WEAK"
    )
  )
}

validation <- map_dfr(candidate_genes, ~ validate_candidate(.x, de_exp1, de_exp2, chipseq_results))

cat("Quality tiers:\n")
cat("  STRONG    = Significantly DE in both experiments + ChIP-seq binding\n")
cat("  DE only   = Significantly DE but no ChIP-seq peaks\n")
cat("  ChIP only = ChIP-seq peaks but not consistently DE\n")
cat("  WEAK      = Neither significant DE nor ChIP-seq binding\n")
cat("  NOT IN DATA = Absent from count matrix\n\n")

print(as.data.frame(validation %>%
  dplyr::select(gene, quality, exp1_log2fc, exp1_padj, smad2_bound, eomesa_bound) %>%
  mutate(exp1_log2fc = round(exp1_log2fc, 2),
         exp1_padj = formatC(exp1_padj, format = "e", digits = 1))),
  row.names = FALSE)

cat("\n")
for (q in c("STRONG", "DE only", "ChIP only", "WEAK", "NOT IN DATA")) {
  g <- validation %>% filter(quality == q) %>% pull(gene)
  if (length(g) > 0) cat(sprintf("  %s (%d): %s\n", q, length(g), paste(g, collapse = ", ")))
}

write_csv(validation, results_path("q4_candidate_validation.csv"))

# Propose alternatives: DE in both experiments AND ChIP-seq binding
all_known <- unique(c(candidate_genes, nodal_genes))

alt_candidates <- de_exp1 %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1.0, !gene %in% all_known,
         !grepl("^LOC|^si:", gene)) %>%
  arrange(padj) %>%
  head(100) %>%
  filter(gene %in% (de_exp2 %>% filter(padj < 0.05, abs(log2FoldChange) >= 1.0) %>% pull(gene)))

alt_chip <- map_dfr(alt_candidates$gene, function(g) {
  tibble(
    gene = g,
    smad2_peaks = sum(chipseq_results$n_peaks[chipseq_results$gene == g & chipseq_results$chip_factor == "Smad2"]),
    eomesa_peaks = sum(chipseq_results$n_peaks[chipseq_results$gene == g & chipseq_results$chip_factor == "EomesA"]),
    has_foxh1 = any(chipseq_results$has_foxh1_motif[chipseq_results$gene == g], na.rm = TRUE)
  )
})

alt_full <- alt_candidates %>%
  dplyr::select(gene, log2FoldChange, padj, baseMean) %>%
  left_join(alt_chip, by = "gene") %>%
  mutate(
    has_chipseq = (smad2_peaks > 0) | (eomesa_peaks > 0),
    score = -log10(padj) + abs(log2FoldChange) * 2 + has_chipseq * 10 + has_foxh1 * 5
  ) %>%
  arrange(desc(score))

cat("\n=== Top 20 Alternative Candidates ===\n")
cat("Score = -log10(padj) + 2*|log2FC| + 10*(ChIP) + 5*(Foxh1)\n\n")
print(as.data.frame(alt_full %>% head(20) %>%
  mutate(log2FC = round(log2FoldChange, 2),
         padj_str = formatC(padj, format = "e", digits = 1),
         baseMean = round(baseMean, 0), score = round(score, 1)) %>%
  dplyr::select(gene, log2FC, padj_str, baseMean, smad2_peaks, eomesa_peaks, has_foxh1, score)),
  row.names = FALSE)

write_csv(alt_full, results_path("q4_alternative_candidates.csv"))

# --- FIGURE 5: Candidate validation volcano ---

val_plot_data <- validation %>%
  filter(in_count_matrix) %>%
  mutate(
    gene = factor(gene, levels = gene[order(quality, exp1_padj)]),
    quality = factor(quality, levels = c("STRONG", "DE only", "ChIP only", "WEAK"))
  ) %>%
  filter(!is.na(quality))

quality_colors <- c("STRONG" = "#2166AC", "DE only" = "#E66101",
                     "ChIP only" = "#B2182B", "WEAK" = "grey75")

p_validation <- ggplot(val_plot_data,
       aes(x = exp1_log2fc, y = -log10(exp1_padj), fill = quality)) +
  geom_point(shape = 21, size = 4, stroke = 0.5, color = "black", alpha = 0.85) +
  geom_text_repel(aes(label = gene), size = 3, fontface = "italic", color = "#333333",
                  max.overlaps = 20, seed = 42, box.padding = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_fill_manual(values = quality_colors, name = "Validation tier") +
  labs(
    x = "log2FC (Exp1: 15 vs 0 ng/ml Activin @ 240 min)",
    y = "-log10(adj. p-value)",
    title = "Candidate gene validation",
    subtitle = "DE significance vs ChIP-seq binding evidence",
    caption = "Dashed lines: |log\u2082FC| = 1 and padj = 0.05"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50"),
    plot.caption = element_text(size = 8, hjust = 0, color = "grey50")
  )

ggsave(results_path("q4_candidate_validation.pdf"), p_validation, width = 9, height = 7)
cat("\nSaved:", results_path("q4_candidate_validation.pdf"), "\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n\n")

cat("ChIP-seq binding (candidates):\n")
cat(sprintf("  Smad2 peaks: %d / %d | EomesA peaks: %d / %d\n",
    sum(candidate_chip_summary$Smad2_bound), nrow(candidate_chip_summary),
    sum(candidate_chip_summary$EomesA_bound), nrow(candidate_chip_summary)))

cat("\nChIP-seq binding (Nodal score):\n")
cat(sprintf("  Smad2 peaks: %d / %d | EomesA peaks: %d / %d\n",
    sum(nodal_chip_summary$Smad2_bound), nrow(nodal_chip_summary),
    sum(nodal_chip_summary$EomesA_bound), nrow(nodal_chip_summary)))

strong <- validation %>% filter(quality == "STRONG") %>% pull(gene)
weak   <- validation %>% filter(quality %in% c("WEAK", "NOT IN DATA")) %>% pull(gene)
cat(sprintf("\nValid candidates (STRONG): %s\n", paste(strong, collapse = ", ")))
cat(sprintf("Weak/missing candidates: %s\n", paste(weak, collapse = ", ")))

cat("\n=== Output Files (all in results/) ===\n")
cat("  CSV:  q4_chipseq_candidates.csv, q4_chipseq_nodal.csv\n")
cat("        q4_peak_expression_times.csv, q4_temporal_clusters.csv\n")
cat("        q4_reversibility_scores.csv, q4_combined_summary.csv\n")
cat("        q4_candidate_validation.csv, q4_alternative_candidates.csv\n")
cat("  PDF:  q4_chipseq_binding_heatmap.pdf  (Fig 1)\n")
cat("        q4_temporal_clusters.pdf        (Fig 2)\n")
cat("        q4_reversibility_profile.pdf    (Fig 3)\n")
cat("        q4_integrated_overview.pdf      (Fig 4)\n")
cat("        q4_candidate_validation.pdf     (Fig 5)\n")
