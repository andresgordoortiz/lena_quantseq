# ============================================================================
# GENE REGULATORY NETWORK: FGF–NODAL SIGNALLING CROSSTALK
# ============================================================================
#
# Biological question:
#   Activin/Nodal activates FGF genes initially, but at high Nodal levels
#   (concentration or prolonged exposure), FGF signalling is blocked through
#   negative feedback — likely via DUSPs (especially DUSP4; van Boxtel et al.
#   2018) and Sprouty/Spred switches (Dorey & Amaya 2010).
#
#   This script builds a data-driven gene regulatory network to:
#     (a) Identify the temporal sequence of pathway activation
#     (b) Detect dose-dependent feedback (genes UP at low Activin, DOWN at high)
#     (c) Find correlation modules linking Nodal → FGF → feedback loops
#     (d) Score candidate negative feedback regulators
#     (e) Visualise the resulting GRN
#
# Data sources:
#   - Exp1: Activin dose-response time course
#       0/5/10/15 ng/ml Activin × 15/30/60/120/180/240 min
#   - Exp2: SB431542 inhibitor (via Q3 commitment dynamics)
#   - Two ChIP-seq datasets: Smad2 & EomesA binding (Nelson et al. 2014)
#   - Literature: Dorey & Amaya 2010, Annual Reviews Nodal 2024,
#     van Boxtel et al. 2018
#
# Outputs: results/grn_*
#

source("preprocess.R")
library(pheatmap)
library(ggrepel)
library(patchwork)
library(igraph)
library(scales)
library(readxl)

# ============================================================================
# 0. EXPANDED GENE UNIVERSE
# ============================================================================
# "Be massive, yet smart in your choices" — cast a wide net for candidates
# but use multiple orthogonal evidence streams to prioritise.

# ---- Nodal/TGFβ pathway core ----
nodal_core <- c(
  # Ligands
  "ndr1", "ndr2", "spaw",
  # Co-receptor
  "tdgf1",
  # Receptors / Type I/II
  "acvr1ba", "acvr2aa", "acvr2ba",
  # Smads
  "smad2", "smad3a", "smad3b", "smad4a",
  # TFs downstream of Nodal
  "foxh1", "eomesa", "mixl1", "gata5",
  # Feedback inhibitors
  "lefty1", "lefty2", "lft1", "lft2"
)

# ---- FGF ligands (comprehensive) ----
fgf_ligands <- c(
  "fgf3", "fgf4", "fgf8a", "fgf8b", "fgf10a", "fgf10b",
  "fgf16", "fgf17", "fgf18a", "fgf18b", "fgf19",
  "fgf20a", "fgf20b", "fgf24"
)

# ---- FGF receptors ----
fgf_receptors <- c("fgfr1a", "fgfr1b", "fgfr2", "fgfr3", "fgfr4")

# ---- MAPK cascade ----
mapk_cascade <- c(
  # RAS/RAF/MEK/ERK transcripts (not all will be DE, but survey them)
  "hras", "kras", "nras",
  "braf", "raf1a", "raf1b",
  "map2k1", "map2k2a",
  "mapk1", "mapk3"
)

# ---- FGF/MAPK feedback targets (known direct) ----
fgf_feedback <- c("etv4", "etv5a", "etv5b",
                   "spry1", "spry2", "spry4",
                   "dusp6", "il17rd")

# ---- DUSP family (broad) ----
dusp_family <- c("dusp1", "dusp2", "dusp4", "dusp5", "dusp6", "dusp7",
                  "dusp10", "dusp11", "dusp14", "dusp16")

# ---- Sprouty/Spred family ----
spry_spred <- c("spry1", "spry2", "spry4",
                 "spred1", "spred2a", "spred2b", "spred3")

# ---- Mesoderm / FGF-dependent ----
mesoderm_genes <- c("tbx16", "tbx6", "tbxta", "noto", "msgn1",
                     "her1", "her7", "hes6",
                     "mespaa", "mespab", "ripply1", "ripply2")

# ---- Endoderm markers (FGF opposes endoderm) ----
endoderm_genes <- c("sox32", "sox17", "gata5", "gata6", "foxa2", "foxa3",
                     "bon", "cas")

# ---- BMP crosstalk (FGF-ERK phosphorylates Smad1 linker → BMP inhibition) ----
bmp_crosstalk <- c("bmp2b", "bmp4", "bmp7a", "bmp7b",
                    "chrd", "nog1", "nog2",
                    "smad1", "smad5")

# ---- Wnt crosstalk ----
wnt_crosstalk <- c("wnt8a", "wnt3a", "wnt11f2",
                    "ctnnb1", "ctnnb2",
                    "axin1", "axin2", "dkk1b")

# ---- Candidate negative feedback regulators ----
# These are genes that could mediate Nodal → FGF inhibition at high doses.
# Selection criteria: known phosphatases, TF repressors, signalling
# modulators implicated in FGF or MAPK pathway regulation.
candidate_feedback <- c(
  # Known DUSPs not in the direct list
  "dusp4", "dusp5",
  # Phosphatase regulators
  "ptpn11", "pten",
  # Transcriptional repressors
  "snai1a", "snai1b",
  "foxd3", "sox9a",
  # Signalling modulators
  "socs3a", "socs3b",
  # Cell adhesion / migration modulators
  "flrt3", "rhov", "net1", "rnd1b",
  # Additional candidates from the experiment
  "efna1a", "efnb2a", "rasgef1ba",
  # PI3K/AKT axis (cross-talks with MAPK)
  "pik3ca", "akt1",
  # GAPs / GEFs for Ras
  "rasa1a", "rasa2", "rasal2",
  "sos1", "sos2",
  # MEK interacting inhibitors
  "rkip"
)

# ---- Nodal score genes (from literature) ----
nodal_score_data <- tryCatch(
  read_excel("docs/nodal-score-genes_complete.xlsx", skip = 1),
  error = function(e) NULL
)
nodal_score_genes <- if (!is.null(nodal_score_data)) {
  unique(tolower(na.omit(nodal_score_data$`Nodal score`)))
} else { character(0) }

# Combine all into a master gene universe (deduplicate)
grn_universe <- unique(tolower(c(
  nodal_core, fgf_ligands, fgf_receptors, mapk_cascade,
  fgf_feedback, dusp_family, spry_spred, mesoderm_genes,
  endoderm_genes, bmp_crosstalk, wnt_crosstalk,
  candidate_feedback, nodal_score_genes
)))

# Map to available genes
grn_available <- grn_universe[grn_universe %in% tolower(rownames(counts_filtered))]
cat(sprintf("\nGRN gene universe: %d defined, %d found in dataset\n",
            length(grn_universe), length(grn_available)))

# Gene group membership (for annotation)
gene_groups <- data.frame(gene = grn_available, stringsAsFactors = FALSE) %>%
  mutate(
    group = case_when(
      gene %in% tolower(nodal_core) ~ "Nodal pathway",
      gene %in% tolower(fgf_ligands) ~ "FGF ligands",
      gene %in% tolower(fgf_receptors) ~ "FGF receptors",
      gene %in% tolower(mapk_cascade) ~ "MAPK cascade",
      gene %in% tolower(fgf_feedback) ~ "FGF/MAPK feedback",
      gene %in% tolower(dusp_family) ~ "DUSP family",
      gene %in% tolower(spry_spred) ~ "Sprouty/Spred",
      gene %in% tolower(mesoderm_genes) ~ "Mesoderm (FGF-dependent)",
      gene %in% tolower(endoderm_genes) ~ "Endoderm markers",
      gene %in% tolower(bmp_crosstalk) ~ "BMP crosstalk",
      gene %in% tolower(wnt_crosstalk) ~ "Wnt crosstalk",
      gene %in% tolower(candidate_feedback) ~ "Candidate feedback",
      gene %in% nodal_score_genes ~ "Nodal target",
      TRUE ~ "Other"
    )
  )

write_csv(gene_groups, results_path("grn_gene_universe.csv"))
cat("Gene groups:\n")
print(table(gene_groups$group))

# ============================================================================
# 1. NORMALISE EXP1 DATA (dose–response time course)
# ============================================================================

cat("\n========== SECTION 1: NORMALISE EXP1 ==========\n")

exp1_meta <- metadata %>% filter(experiment == "Exp1")
exp1_counts <- counts_filtered[, rownames(exp1_meta)]

dds_exp1 <- DESeqDataSetFromMatrix(exp1_counts, exp1_meta, ~ 1)
dds_exp1 <- estimateSizeFactors(dds_exp1)
norm_exp1 <- counts(dds_exp1, normalized = TRUE)

# Extract concentrations and timepoints
concentrations <- c("0ngml", "5ngml", "10ngml", "15ngml")
timepoints <- sort(unique(exp1_meta$time_min))
cat("Concentrations:", paste(concentrations, collapse = ", "), "\n")
cat("Timepoints:", paste(timepoints, collapse = ", "), "min\n")

# ============================================================================
# 2. DOSE–RESPONSE TEMPORAL PROFILES
# ============================================================================
#
# For each gene × concentration × timepoint: compute log2FC vs 0 ng/ml.
# This lets us see how dose shapes the temporal response.

cat("\n========== SECTION 2: DOSE-RESPONSE PROFILES ==========\n")

exp1_0ngml <- exp1_meta %>% filter(concentration == "0ngml")

compute_dose_profile <- function(genes, norm_mat, meta_all, meta_0, concs, tp) {
  genes_found <- genes[tolower(genes) %in% tolower(rownames(norm_mat))]
  purrr::map_dfr(genes_found, function(gene) {
    row_idx <- which(tolower(rownames(norm_mat)) == gene)
    if (length(row_idx) == 0) return(NULL)
    purrr::map_dfr(concs[concs != "0ngml"], function(conc) {
      meta_conc <- meta_all %>% filter(concentration == conc)
      purrr::map_dfr(tp, function(t) {
        treated <- norm_mat[row_idx, rownames(meta_conc)[meta_conc$time_min == t]]
        control <- norm_mat[row_idx, rownames(meta_0)[meta_0$time_min == t]]
        if (length(treated) == 0 || length(control) == 0) return(NULL)
        tibble(gene = gene, concentration = conc, time_min = t,
               log2fc = log2(mean(treated) + 1) - log2(mean(control) + 1),
               mean_treated = mean(treated),
               mean_control = mean(control))
      })
    })
  })
}

dose_profiles <- compute_dose_profile(
  grn_available, norm_exp1, exp1_meta, exp1_0ngml, concentrations, timepoints
)

cat(sprintf("Dose-response profiles: %d genes × %d conditions\n",
            n_distinct(dose_profiles$gene), n_distinct(dose_profiles$concentration)))

write_csv(dose_profiles, results_path("grn_dose_response_profiles.csv"))

# ============================================================================
# 3. IDENTIFY DOSE-DEPENDENT PATTERNS
# ============================================================================
#
# Key pattern: "activated at low dose, repressed at high dose" → feedback.
# Also: monotonic increase, saturation, bell-shaped.

cat("\n========== SECTION 3: DOSE-DEPENDENT PATTERNS ==========\n")

# For each gene, compute mean log2FC across shared timepoints per concentration
shared_tp <- timepoints[timepoints >= 60]  # 60, 120, 180, 240 min — shared across all concs

dose_summary <- dose_profiles %>%
  filter(time_min %in% shared_tp) %>%
  group_by(gene, concentration) %>%
  summarise(mean_log2fc = mean(log2fc, na.rm = TRUE),
            max_log2fc = max(log2fc, na.rm = TRUE),
            .groups = "drop")

# Classify dose-response pattern
dose_pattern <- dose_summary %>%
  pivot_wider(id_cols = gene, names_from = concentration,
              values_from = mean_log2fc, names_prefix = "fc_") %>%
  mutate(
    # Monotonic increasing: 5 < 10 < 15
    monotonic_up = !is.na(fc_5ngml) & !is.na(fc_10ngml) & !is.na(fc_15ngml) &
                   fc_5ngml < fc_10ngml & fc_10ngml < fc_15ngml & fc_15ngml > 0.5,

    # Monotonic decreasing: 5 > 10 > 15, all negative
    monotonic_down = !is.na(fc_5ngml) & !is.na(fc_10ngml) & !is.na(fc_15ngml) &
                     fc_5ngml > fc_10ngml & fc_10ngml > fc_15ngml & fc_15ngml < -0.5,

    # Bell-shaped (feedback signature): up at 5, down at 15 relative to 5
    # OR: 10 > 5 and 10 > 15 (peak at middle dose)
    bell_shaped = !is.na(fc_5ngml) & !is.na(fc_10ngml) & !is.na(fc_15ngml) &
                  fc_10ngml > fc_5ngml & fc_10ngml > fc_15ngml &
                  fc_10ngml > 0.3,

    # Saturation: increases then plateaus (5 << 10 ≈ 15)
    saturating = !is.na(fc_5ngml) & !is.na(fc_10ngml) & !is.na(fc_15ngml) &
                 fc_10ngml > fc_5ngml + 0.3 &
                 abs(fc_15ngml - fc_10ngml) < 0.2 * abs(fc_10ngml),

    # Classify
    pattern = case_when(
      bell_shaped ~ "Bell-shaped (feedback)",
      monotonic_up ~ "Monotonic activation",
      monotonic_down ~ "Monotonic repression",
      saturating ~ "Saturating activation",
      TRUE ~ "Other/weak"
    )
  )

cat("\nDose-response patterns:\n")
print(table(dose_pattern$pattern))

# Annotate with gene groups
dose_pattern <- dose_pattern %>%
  left_join(gene_groups, by = "gene")

write_csv(dose_pattern, results_path("grn_dose_patterns.csv"))

# ============================================================================
# 4. TEMPORAL PROFILE AT 15 ng/ml (PRIMARY ACTIVATION DYNAMICS)
# ============================================================================

cat("\n========== SECTION 4: TEMPORAL DYNAMICS AT 15 ng/ml ==========\n")

exp1_15ngml <- exp1_meta %>% filter(concentration == "15ngml")

temporal_15 <- purrr::map_dfr(grn_available, function(gene) {
  row_idx <- which(tolower(rownames(norm_exp1)) == gene)
  if (length(row_idx) == 0) return(NULL)
  purrr::map_dfr(timepoints, function(t) {
    treated <- norm_exp1[row_idx, rownames(exp1_15ngml)[exp1_15ngml$time_min == t]]
    control <- norm_exp1[row_idx, rownames(exp1_0ngml)[exp1_0ngml$time_min == t]]
    if (length(treated) == 0 || length(control) == 0) return(NULL)
    tibble(gene = gene, time_min = t,
           log2fc = log2(mean(treated) + 1) - log2(mean(control) + 1))
  })
})

temporal_15 <- temporal_15 %>%
  left_join(gene_groups, by = "gene")

write_csv(temporal_15, results_path("grn_temporal_15ngml.csv"))

# Compute activation timing metrics
timing_metrics <- temporal_15 %>%
  group_by(gene) %>%
  summarise(
    peak_time = time_min[which.max(abs(log2fc))],
    peak_log2fc = log2fc[which.max(abs(log2fc))],
    # Time of first significant response (|log2fc| > 0.5)
    onset_time = {
      sig <- time_min[abs(log2fc) > 0.5]
      if (length(sig) > 0) min(sig) else NA_real_
    },
    # Direction
    direction = ifelse(sum(log2fc) >= 0, "up", "down"),
    # Area under curve (signed)
    auc = sum(log2fc),
    # Late vs early ratio (to detect delayed responders)
    late_vs_early = mean(log2fc[time_min >= 120]) - mean(log2fc[time_min <= 60]),
    .groups = "drop"
  ) %>%
  left_join(gene_groups, by = "gene")

write_csv(timing_metrics, results_path("grn_timing_metrics.csv"))

cat("Timing metrics computed for", nrow(timing_metrics), "genes\n")
cat("Direction: ", sum(timing_metrics$direction == "up"), "up,",
    sum(timing_metrics$direction == "down"), "down\n")

# ============================================================================
# 5. CORRELATION NETWORK
# ============================================================================
#
# Build a correlation matrix of temporal profiles (at 15 ng/ml) across the GRN
# universe. Genes with anti-correlated temporal profiles to FGF targets are
# candidate negative regulators.

cat("\n========== SECTION 5: CORRELATION NETWORK ==========\n")

# Build a gene × timepoint matrix for correlation
profile_mat <- temporal_15 %>%
  dplyr::select(gene, time_min, log2fc) %>%
  pivot_wider(names_from = time_min, values_from = log2fc) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Remove genes with zero variance (uninformative)
gene_var <- apply(profile_mat, 1, var, na.rm = TRUE)
profile_mat <- profile_mat[gene_var > 0.01, , drop = FALSE]
cat(sprintf("Correlation matrix: %d genes with sufficient variance\n", nrow(profile_mat)))

# Compute pairwise Pearson correlation
cor_mat <- cor(t(profile_mat), use = "pairwise.complete.obs", method = "pearson")

# ---- FGF score: mean correlation with FGF feedback genes ----
fgf_feedback_in <- intersect(fgf_feedback, rownames(cor_mat))
if (length(fgf_feedback_in) >= 3) {
  cor_with_fgf <- rowMeans(cor_mat[, fgf_feedback_in, drop = FALSE], na.rm = TRUE)
  fgf_cor_df <- data.frame(gene = names(cor_with_fgf), cor_with_fgf = cor_with_fgf,
                            stringsAsFactors = FALSE) %>%
    left_join(gene_groups, by = "gene") %>%
    arrange(cor_with_fgf)

  cat("\nMost ANTI-correlated with FGF feedback targets (candidate inhibitors):\n")
  print(head(fgf_cor_df %>% filter(!gene %in% fgf_feedback), 20))

  cat("\nMost POSITIVELY correlated with FGF feedback (co-activated):\n")
  print(tail(fgf_cor_df %>% filter(!gene %in% fgf_feedback), 20))

  write_csv(fgf_cor_df, results_path("grn_fgf_correlation_ranking.csv"))
}

# ---- Nodal–FGF correlation ----
nodal_in <- intersect(tolower(nodal_core), rownames(cor_mat))
if (length(nodal_in) >= 2 && length(fgf_feedback_in) >= 2) {
  nodal_fgf_cors <- cor_mat[nodal_in, fgf_feedback_in, drop = FALSE]
  cat("\nNodal ↔ FGF feedback correlations:\n")
  print(round(nodal_fgf_cors, 3))
}

# ============================================================================
# 6. CANDIDATE NEGATIVE FEEDBACK SCORING
# ============================================================================
#
# Score each gene as a candidate negative feedback regulator based on:
#   (a) Anti-correlation with FGF targets over time
#   (b) Delayed activation relative to FGF (onset_time > FGF onset)
#   (c) Dose-dependent repression at high activin (bell-shaped dose-response
#       or reduced FC at 15 vs 5/10 ng/ml)
#   (d) Smad2/EomesA ChIP-seq binding (direct Nodal target → can mediate feedback)
#   (e) Positive correlation with Nodal pathway genes (co-regulated with Nodal)
#   (f) Upregulated at 15 ng/ml + late timing = delayed negative regulator

cat("\n========== SECTION 6: CANDIDATE SCORING ==========\n")

# Load ChIP-seq data
smad2_chipseq <- tryCatch(
  read_excel("12915_2014_81_MOESM4_ESM.xlsx", skip = 3),
  error = function(e) { cat("Smad2 ChIP-seq file not found\n"); NULL }
)
eomesa_chipseq <- tryCatch(
  read_excel("12915_2014_81_MOESM10_ESM.xlsx", skip = 2),
  error = function(e) { cat("EomesA ChIP-seq file not found\n"); NULL }
)

# Check ChIP binding for each gene
check_chip <- function(gene_name, chipseq_df) {
  if (is.null(chipseq_df)) return(FALSE)
  pattern <- paste0("(?i)(^|;)", gene_name)
  any(grepl(pattern, chipseq_df$`Proximal gene`, perl = TRUE))
}

# Build scoring dataframe
candidate_scores <- timing_metrics %>%
  dplyr::select(gene, peak_time, peak_log2fc, onset_time, direction,
                auc, late_vs_early, group) %>%
  # (a) Anti-correlation with FGF targets
  left_join(fgf_cor_df %>% dplyr::select(gene, cor_with_fgf), by = "gene") %>%
  # (b) Delayed activation
  mutate(
    fgf_median_onset = median(timing_metrics$onset_time[
      timing_metrics$gene %in% fgf_feedback], na.rm = TRUE),
    delayed = !is.na(onset_time) & onset_time > fgf_median_onset
  ) %>%
  # (c) Dose-dependent feedback (from Section 3)
  left_join(dose_pattern %>% dplyr::select(gene, pattern, fc_5ngml, fc_10ngml, fc_15ngml),
            by = "gene") %>%
  # (d) ChIP-seq binding
  rowwise() %>%
  mutate(
    smad2_bound = check_chip(gene, smad2_chipseq),
    eomesa_bound = check_chip(gene, eomesa_chipseq),
    chip_bound = smad2_bound | eomesa_bound
  ) %>%
  ungroup()

# Correlation with Nodal core
if (length(nodal_in) >= 2) {
  cor_with_nodal <- rowMeans(cor_mat[, nodal_in, drop = FALSE], na.rm = TRUE)
  nodal_cor_df <- data.frame(gene = names(cor_with_nodal),
                              cor_with_nodal = cor_with_nodal,
                              stringsAsFactors = FALSE)
  candidate_scores <- candidate_scores %>%
    left_join(nodal_cor_df, by = "gene")
} else {
  candidate_scores$cor_with_nodal <- NA
}

# ---- Composite score ----
# Higher score = stronger candidate for negative feedback regulation.
# Components (each scaled 0–1):
#   1. Anti-FGF correlation: −cor_with_fgf → higher when anticorrelated
#   2. Positive Nodal correlation: cor_with_nodal → higher when co-expressed with Nodal
#   3. Delayed onset: 1 if delayed, 0 if not
#   4. Bell-shaped dose: 1 if bell-shaped, 0.3 if saturating, 0 otherwise
#   5. ChIP-seq bound: 1 if Smad2/EomesA, 0 otherwise
#   6. Upregulated + late: bonus for genes that are UP but peak late

scale01 <- function(x) {
  x[is.na(x)] <- 0
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / diff(rng)
}

candidate_scores <- candidate_scores %>%
  mutate(
    # Component 1: anti-FGF (invert)
    score_antifgf = scale01(-replace_na(cor_with_fgf, 0)),
    # Component 2: pro-Nodal
    score_pronodal = scale01(replace_na(cor_with_nodal, 0)),
    # Component 3: delayed onset
    score_delayed = ifelse(delayed, 1, 0),
    # Component 4: dose pattern
    score_dose = case_when(
      pattern == "Bell-shaped (feedback)" ~ 1.0,
      pattern == "Saturating activation" ~ 0.3,
      pattern == "Monotonic repression" ~ 0.5,
      TRUE ~ 0
    ),
    # Component 5: ChIP bound
    score_chip = ifelse(chip_bound, 1, 0),
    # Component 6: up + late
    score_uplate = ifelse(direction == "up" & late_vs_early > 0.3, 0.5, 0),

    # Composite (equal weights)
    composite_score = (score_antifgf + score_pronodal + score_delayed +
                        score_dose + score_chip + score_uplate) / 6
  ) %>%
  arrange(desc(composite_score))

write_csv(candidate_scores, results_path("grn_candidate_scores.csv"))

cat("\n=== TOP 30 CANDIDATE NEGATIVE FEEDBACK REGULATORS ===\n")
candidate_scores %>%
  filter(!gene %in% fgf_feedback) %>%
  head(30) %>%
  dplyr::select(gene, group, composite_score, cor_with_fgf, cor_with_nodal,
                delayed, pattern, chip_bound, peak_time, peak_log2fc) %>%
  print(n = 30, width = 200)

# ============================================================================
# 7. FIGURES
# ============================================================================

cat("\n========== SECTION 7: FIGURES ==========\n")

# ---- 7a. Dose-response heatmap for key gene families ----
key_families <- c("Nodal pathway", "FGF/MAPK feedback", "FGF ligands",
                   "DUSP family", "Mesoderm (FGF-dependent)", "Endoderm markers")

dose_heatmap_data <- dose_profiles %>%
  filter(time_min %in% shared_tp) %>%
  left_join(gene_groups, by = "gene") %>%
  filter(group %in% key_families) %>%
  group_by(gene, concentration, group) %>%
  summarise(mean_fc = mean(log2fc), .groups = "drop") %>%
  mutate(concentration = factor(concentration, levels = c("5ngml", "10ngml", "15ngml"),
                                labels = c("5", "10", "15")))

# Build matrix for pheatmap
dose_hm_wide <- dose_heatmap_data %>%
  unite("cond", concentration, sep = " ng/ml") %>%
  pivot_wider(names_from = cond, values_from = mean_fc) %>%
  arrange(group, gene)

dose_hm_mat <- dose_hm_wide %>%
  dplyr::select(-gene, -group) %>%
  as.matrix()
rownames(dose_hm_mat) <- dose_hm_wide$gene

# Annotation for gene group
row_ann <- data.frame(Group = dose_hm_wide$group, row.names = dose_hm_wide$gene)

group_colors <- c(
  "Nodal pathway" = "#7570B3",
  "FGF/MAPK feedback" = "#E66101",
  "FGF ligands" = "#1B9E77",
  "DUSP family" = "#D95F02",
  "Mesoderm (FGF-dependent)" = "#66A61E",
  "Endoderm markers" = "#E7298A"
)

ann_colors <- list(Group = group_colors[levels(factor(row_ann$Group))])

# Plot
pdf(results_path("grn_dose_response_heatmap.pdf"), width = 8, height = 14)
pheatmap(dose_hm_mat,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = row_ann,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = seq(-4, 4, length.out = 101),
         fontsize_row = 8,
         main = "Dose-response: mean log2FC vs 0 ng/ml (60-240 min)",
         gaps_row = cumsum(table(dose_hm_wide$group)))
dev.off()
cat("Saved:", results_path("grn_dose_response_heatmap.pdf"), "\n")

# ---- 7b. Temporal dynamics faceted by gene group ----
temporal_plot_data <- temporal_15 %>%
  filter(group %in% c("Nodal pathway", "FGF/MAPK feedback", "FGF ligands",
                       "DUSP family", "Mesoderm (FGF-dependent)", "Endoderm markers",
                       "Sprouty/Spred", "BMP crosstalk", "Wnt crosstalk"))

family_means <- temporal_plot_data %>%
  group_by(group, time_min) %>%
  summarise(mean_fc = mean(log2fc), .groups = "drop")

endpoint_labels <- temporal_plot_data %>%
  filter(time_min == max(time_min))

family_palette <- c(
  "Nodal pathway" = "#7570B3",
  "FGF/MAPK feedback" = "#E66101",
  "FGF ligands" = "#1B9E77",
  "DUSP family" = "#D95F02",
  "Sprouty/Spred" = "#E7298A",
  "Mesoderm (FGF-dependent)" = "#66A61E",
  "Endoderm markers" = "#A6761D",
  "BMP crosstalk" = "#666666",
  "Wnt crosstalk" = "#1F78B4"
)

p_temporal <- ggplot(temporal_plot_data, aes(x = time_min, y = log2fc)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_line(aes(group = gene, color = group), alpha = 0.3, linewidth = 0.4) +
  geom_line(data = family_means, aes(y = mean_fc, color = group),
            linewidth = 1.3) +
  geom_point(data = family_means, aes(y = mean_fc, color = group), size = 2) +
  geom_text_repel(data = endpoint_labels,
                  aes(label = gene, color = group),
                  size = 2, fontface = "italic", direction = "y",
                  hjust = 0, nudge_x = 8, segment.size = 0.2,
                  max.overlaps = 25, show.legend = FALSE) +
  facet_wrap(~ group, ncol = 3, scales = "free_y") +
  scale_color_manual(values = family_palette, guide = "none") +
  scale_x_continuous(breaks = timepoints, expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = "Time (min)", y = "log2FC (15 vs 0 ng/ml Activin)",
       title = "GRN pathway temporal dynamics",
       subtitle = "Thin = individual genes | Thick = group mean") +
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50")
  )

ggsave(results_path("grn_temporal_dynamics.pdf"), p_temporal, width = 16, height = 14)
cat("Saved:", results_path("grn_temporal_dynamics.pdf"), "\n")

# ---- 7c. Dose-response curves for selected genes ----
highlight_genes <- c("dusp4", "dusp6", "etv4", "etv5a", "spry2", "spry4",
                      "fgf3", "fgf17", "fgf24", "il17rd",
                      "ndr1", "ndr2", "lefty1", "lefty2",
                      "tbxta", "tbx16", "msgn1",
                      "sox32", "sox17",
                      "dusp1", "dusp5")
highlight_genes <- highlight_genes[highlight_genes %in% dose_profiles$gene]

dose_highlight <- dose_profiles %>%
  filter(gene %in% highlight_genes) %>%
  left_join(gene_groups, by = "gene") %>%
  mutate(concentration = factor(concentration,
                                levels = c("5ngml", "10ngml", "15ngml"),
                                labels = c("5 ng/ml", "10 ng/ml", "15 ng/ml")))

p_dose <- ggplot(dose_highlight, aes(x = time_min, y = log2fc, color = concentration)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  scale_color_manual(values = c("5 ng/ml" = "#FEE08B", "10 ng/ml" = "#FC8D59",
                                 "15 ng/ml" = "#D73027"),
                     name = "Activin dose") +
  scale_x_continuous(breaks = c(60, 120, 180, 240)) +
  labs(x = "Time (min)", y = "log2FC vs 0 ng/ml",
       title = "Dose-response curves: key GRN genes",
       subtitle = "Each panel shows one gene across three Activin concentrations") +
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 9, face = "bold.italic"),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50")
  )

ggsave(results_path("grn_dose_response_curves.pdf"), p_dose, width = 16, height = 12)
cat("Saved:", results_path("grn_dose_response_curves.pdf"), "\n")

# ---- 7d. Correlation heatmap: FGF targets vs candidate regulators ----
# Select genes for a focused correlation view
focus_genes <- unique(c(
  fgf_feedback, "dusp4", "dusp5", "dusp1",
  "fgf3", "fgf17", "fgf24", "fgf8a",
  "ndr1", "ndr2", "lefty1", "lefty2",
  "tbxta", "tbx16", "msgn1",
  "sox32", "sox17",
  "spred1", "spred2a",
  "flrt3", "rhov", "net1"
))
focus_genes <- focus_genes[focus_genes %in% rownames(cor_mat)]

if (length(focus_genes) >= 5) {
  cor_focus <- cor_mat[focus_genes, focus_genes]

  # Annotation for gene group
  focus_ann <- data.frame(
    Group = gene_groups$group[match(focus_genes, gene_groups$gene)],
    row.names = focus_genes
  )

  pdf(results_path("grn_correlation_heatmap.pdf"), width = 12, height = 11)
  pheatmap(cor_focus,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = seq(-1, 1, length.out = 101),
           annotation_row = focus_ann,
           annotation_col = focus_ann,
           fontsize = 9,
           main = "Temporal profile correlation (15 ng/ml Activin)")
  dev.off()
  cat("Saved:", results_path("grn_correlation_heatmap.pdf"), "\n")
}

# ---- 7e. Gene regulatory network graph ----
# Build a directed network based on:
#   - Literature edges (known regulatory relationships)
#   - Strong correlations / anti-correlations from data

cat("\n--- Building GRN graph ---\n")

# Literature-based edges (directed: source → target, with sign)
lit_edges <- tribble(
  ~from,     ~to,        ~type,       ~sign,
  # Nodal → FGF ligands (van Boxtel 2018, Annual Reviews 2024)
  "ndr1",    "fgf3",     "activates", +1,
  "ndr1",    "fgf17",    "activates", +1,
  "ndr1",    "fgf8a",    "activates", +1,
  "ndr2",    "fgf3",     "activates", +1,
  "ndr2",    "fgf17",    "activates", +1,
  "ndr2",    "fgf8a",    "activates", +1,
  # Nodal → Dusp4 (van Boxtel 2018: Nodal induces intracellular FGF inhibition via Dusp4)
  "ndr1",    "dusp4",    "activates", +1,
  "ndr2",    "dusp4",    "activates", +1,
  # Nodal → Lefty (negative feedback)
  "ndr1",    "lefty1",   "activates", +1,
  "ndr1",    "lefty2",   "activates", +1,
  "ndr2",    "lefty1",   "activates", +1,
  "ndr2",    "lefty2",   "activates", +1,
  # Lefty −| Nodal
  "lefty1",  "ndr1",     "inhibits",  -1,
  "lefty1",  "ndr2",     "inhibits",  -1,
  "lefty2",  "ndr1",     "inhibits",  -1,
  "lefty2",  "ndr2",     "inhibits",  -1,
  # FGF ligands → FGFR → MAPK
  "fgf3",    "fgfr1a",   "activates", +1,
  "fgf8a",   "fgfr1a",   "activates", +1,
  "fgf17",   "fgfr1a",   "activates", +1,
  "fgf24",   "fgfr1a",   "activates", +1,
  "fgfr1a",  "mapk3",    "activates", +1,
  "fgfr1a",  "mapk1",    "activates", +1,
  # MAPK → ETS targets
  "mapk3",   "etv4",     "activates", +1,
  "mapk3",   "etv5a",    "activates", +1,
  "mapk3",   "etv5b",    "activates", +1,
  "mapk1",   "etv4",     "activates", +1,
  # MAPK → Sprouty (negative feedback)
  "mapk3",   "spry2",    "activates", +1,
  "mapk3",   "spry4",    "activates", +1,
  "mapk3",   "dusp6",    "activates", +1,
  "mapk3",   "il17rd",   "activates", +1,
  # Sprouty −| MAPK (feedback)
  "spry2",   "mapk3",    "inhibits",  -1,
  "spry4",   "mapk3",    "inhibits",  -1,
  # DUSPs −| MAPK
  "dusp6",   "mapk3",    "inhibits",  -1,
  "dusp6",   "mapk1",    "inhibits",  -1,
  "dusp4",   "mapk3",    "inhibits",  -1,
  "dusp4",   "mapk1",    "inhibits",  -1,
  # il17rd (Sef) −| FGFR signalling
  "il17rd",  "fgfr1a",   "inhibits",  -1,
  # FGF → mesoderm (permissive)
  "mapk3",   "tbxta",    "activates", +1,
  "mapk3",   "tbx16",    "activates", +1,
  "mapk3",   "msgn1",    "activates", +1,
  # Nodal → endoderm (high dose)
  "ndr1",    "sox32",    "activates", +1,
  "ndr2",    "sox32",    "activates", +1,
  "sox32",   "sox17",    "activates", +1,
  # Nodal → Smad2 → Foxh1/Eomesa
  "ndr1",    "foxh1",    "activates", +1,
  "ndr1",    "eomesa",   "activates", +1,
  "foxh1",   "ndr1",     "activates", +1,  # Foxh1 amplifies Nodal
  "eomesa",  "ndr1",     "activates", +1
)

# Only keep edges where both nodes are in our dataset
lit_edges <- lit_edges %>%
  filter(from %in% grn_available, to %in% grn_available)

# Add data-driven edges: strong correlations (|r| > 0.85)
# between non-redundant gene groups
if (exists("cor_mat") && nrow(cor_mat) > 5) {
  data_edges <- which(abs(cor_mat) > 0.85 & upper.tri(cor_mat), arr.ind = TRUE)
  if (nrow(data_edges) > 0) {
    data_edge_df <- data.frame(
      from = rownames(cor_mat)[data_edges[, 1]],
      to = rownames(cor_mat)[data_edges[, 2]],
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        r = purrr::map2_dbl(from, to, ~ cor_mat[.x, .y]),
        type = ifelse(r > 0, "co-expressed", "anti-expressed"),
        sign = sign(r)
      ) %>%
      # Exclude edges within the same known group (redundant)
      left_join(gene_groups %>% rename(from_group = group), by = c("from" = "gene")) %>%
      left_join(gene_groups %>% rename(to_group = group), by = c("to" = "gene")) %>%
      filter(from_group != to_group) %>%
      dplyr::select(from, to, type, sign)

    cat(sprintf("Data-driven edges (|r| > 0.85, cross-group): %d\n", nrow(data_edge_df)))
  } else {
    data_edge_df <- tibble(from = character(), to = character(),
                           type = character(), sign = numeric())
  }
} else {
  data_edge_df <- tibble(from = character(), to = character(),
                         type = character(), sign = numeric())
}

all_edges <- bind_rows(
  lit_edges %>% mutate(source = "literature"),
  data_edge_df %>% mutate(source = "data")
)

write_csv(all_edges, results_path("grn_network_edges.csv"))

# Build igraph
g <- graph_from_data_frame(all_edges %>% dplyr::select(from, to),
                            directed = TRUE,
                            vertices = data.frame(name = unique(c(all_edges$from, all_edges$to))))

# Add node attributes
V(g)$group <- gene_groups$group[match(V(g)$name, gene_groups$gene)]
V(g)$group[is.na(V(g)$group)] <- "Other"

# Node colors by group
node_palette <- c(
  "Nodal pathway" = "#7570B3",
  "FGF/MAPK feedback" = "#E66101",
  "FGF ligands" = "#1B9E77",
  "DUSP family" = "#D95F02",
  "Sprouty/Spred" = "#E7298A",
  "Mesoderm (FGF-dependent)" = "#66A61E",
  "Endoderm markers" = "#A6761D",
  "FGF receptors" = "#984EA3",
  "MAPK cascade" = "#FF7F00",
  "BMP crosstalk" = "#999999",
  "Wnt crosstalk" = "#1F78B4",
  "Candidate feedback" = "#F781BF",
  "Nodal target" = "#B15928",
  "Other" = "#CCCCCC"
)
V(g)$color <- node_palette[V(g)$group]

# Edge colors by sign
E(g)$sign <- all_edges$sign[1:ecount(g)]
E(g)$type <- all_edges$type[1:ecount(g)]
E(g)$color <- ifelse(E(g)$sign > 0, "#2166AC80", "#B2182B80")
E(g)$lty <- ifelse(all_edges$source[1:ecount(g)] == "literature", 1, 3)
E(g)$arrow.size <- 0.4

# Layout
set.seed(42)
layout <- layout_with_fr(g, niter = 1000)

# Plot network
pdf(results_path("grn_network_graph.pdf"), width = 14, height = 14)
par(mar = c(1, 1, 3, 1))
plot(g, layout = layout,
     vertex.size = 8,
     vertex.label.cex = 0.6,
     vertex.label.font = 3,  # italic
     vertex.label.color = "black",
     vertex.frame.color = NA,
     edge.arrow.size = 0.4,
     edge.curved = 0.15,
     main = "FGF–Nodal Gene Regulatory Network")
# Legend
legend("bottomleft",
       legend = names(node_palette)[names(node_palette) %in% V(g)$group],
       fill = node_palette[names(node_palette) %in% V(g)$group],
       border = NA, cex = 0.7, bty = "n", title = "Gene group")
legend("bottomright",
       legend = c("Activation (literature)", "Inhibition (literature)",
                  "Co-expression (data)", "Anti-expression (data)"),
       col = c("#2166AC", "#B2182B", "#2166AC80", "#B2182B80"),
       lty = c(1, 1, 3, 3), lwd = 2, cex = 0.7, bty = "n", title = "Edge type")
dev.off()
cat("Saved:", results_path("grn_network_graph.pdf"), "\n")

# ---- 7f. Candidate scoring dot plot ----
top_candidates <- candidate_scores %>%
  filter(!gene %in% fgf_feedback, !gene %in% tolower(nodal_core)) %>%
  head(25)

if (nrow(top_candidates) > 0) {
  p_candidates <- ggplot(top_candidates,
                          aes(x = cor_with_fgf, y = reorder(gene, composite_score))) +
    geom_point(aes(size = composite_score,
                   fill = group),
               shape = 21, stroke = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_size_continuous(range = c(3, 10), name = "Composite\nscore") +
    scale_fill_manual(values = node_palette, name = "Gene group") +
    labs(x = "Correlation with FGF feedback targets",
         y = "",
         title = "Top candidate negative feedback regulators",
         subtitle = "Score integrates: anti-FGF correlation, delay, dose-response, ChIP-seq, Nodal co-expression") +
    theme_minimal(base_size = 11, base_family = "Helvetica") +
    theme(
      axis.text.y = element_text(size = 10, face = "italic"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey50")
    )

  ggsave(results_path("grn_candidate_dotplot.pdf"), p_candidates, width = 10, height = 9)
  cat("Saved:", results_path("grn_candidate_dotplot.pdf"), "\n")
}

# ---- 7g. Nodal → FGF activation timeline ----
# Show selected genes ordered by onset time
timeline_genes <- timing_metrics %>%
  filter(gene %in% c("ndr1", "ndr2", "lefty1", "lefty2",
                      "fgf3", "fgf17", "fgf24", "fgf8a",
                      "etv4", "etv5a", "spry2", "spry4",
                      "dusp4", "dusp6", "il17rd",
                      "tbxta", "tbx16", "msgn1", "sox32"),
         !is.na(onset_time)) %>%
  arrange(onset_time)

if (nrow(timeline_genes) > 0) {
  p_timeline <- ggplot(timeline_genes,
                        aes(x = onset_time, y = reorder(gene, -onset_time))) +
    geom_segment(aes(xend = peak_time, yend = gene, color = group),
                 linewidth = 2, alpha = 0.7) +
    geom_point(aes(x = onset_time, color = group), size = 3, shape = 16) +
    geom_point(aes(x = peak_time, color = group), size = 4, shape = 18) +
    scale_color_manual(values = node_palette, name = "Gene group") +
    scale_x_continuous(breaks = timepoints,
                       limits = c(min(timepoints), max(timepoints))) +
    labs(x = "Time (min)", y = "",
         title = "Activation timeline: onset → peak",
         subtitle = "Circle = onset (|log2FC| > 0.5) | Diamond = peak response time") +
    theme_minimal(base_size = 11, base_family = "Helvetica") +
    theme(
      axis.text.y = element_text(size = 10, face = "italic"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "#F0F0F0"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey50")
    )

  ggsave(results_path("grn_activation_timeline.pdf"), p_timeline, width = 10, height = 8)
  cat("Saved:", results_path("grn_activation_timeline.pdf"), "\n")
}

# ---- 7h. FGF feedback vs Nodal: dose-response comparison ----
# Does FGF feedback increase linearly with dose, or show saturation/inversion?
fgf_vs_nodal_genes <- c("etv4", "etv5a", "spry2", "dusp6",     # FGF feedback
                         "dusp4", "dusp5", "dusp1",              # candidate DUSPs
                         "ndr1", "ndr2", "lefty1", "lefty2")     # Nodal reference
fgf_vs_nodal_genes <- fgf_vs_nodal_genes[fgf_vs_nodal_genes %in% dose_profiles$gene]

if (length(fgf_vs_nodal_genes) >= 4) {
  dose_at_late <- dose_profiles %>%
    filter(gene %in% fgf_vs_nodal_genes, time_min >= 120) %>%
    group_by(gene, concentration) %>%
    summarise(mean_fc = mean(log2fc), .groups = "drop") %>%
    left_join(gene_groups, by = "gene") %>%
    mutate(conc_num = as.numeric(gsub("ngml", "", concentration)))

  p_dose_compare <- ggplot(dose_at_late, aes(x = conc_num, y = mean_fc, color = group)) +
    geom_line(aes(group = gene), linewidth = 0.8) +
    geom_point(size = 2) +
    geom_text_repel(data = dose_at_late %>% filter(conc_num == 15),
                    aes(label = gene), size = 3, fontface = "italic",
                    hjust = 0, nudge_x = 0.5, show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = node_palette, name = "Gene group") +
    scale_x_continuous(breaks = c(5, 10, 15)) +
    labs(x = "Activin concentration (ng/ml)",
         y = "Mean log2FC vs 0 ng/ml (120-240 min)",
         title = "Dose-response: FGF feedback vs Nodal pathway",
         subtitle = "Do FGF targets saturate or decrease at high Activin?") +
    theme_minimal(base_size = 11, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey50")
    )

  ggsave(results_path("grn_dose_fgf_vs_nodal.pdf"), p_dose_compare, width = 10, height = 7)
  cat("Saved:", results_path("grn_dose_fgf_vs_nodal.pdf"), "\n")
}

# ============================================================================
# 8. GENOME-WIDE SCREEN FOR NOVEL FEEDBACK CANDIDATES
# ============================================================================
#
# Strategy: scan ALL expressed genes for the "late-onset anti-FGF" pattern.
# This goes beyond the curated GRN universe to find novel candidates.

cat("\n========== SECTION 8: GENOME-WIDE FEEDBACK SCREEN ==========\n")

# Compute temporal profiles for ALL genes at 15 ng/ml
all_gene_names <- tolower(rownames(norm_exp1))

# Efficient batch computation
temporal_all <- purrr::map_dfr(all_gene_names, function(gene) {
  row_idx <- which(tolower(rownames(norm_exp1)) == gene)
  if (length(row_idx) == 0) return(NULL)
  purrr::map_dfr(timepoints, function(t) {
    treated <- norm_exp1[row_idx, rownames(exp1_15ngml)[exp1_15ngml$time_min == t]]
    control <- norm_exp1[row_idx, rownames(exp1_0ngml)[exp1_0ngml$time_min == t]]
    if (length(treated) == 0 || length(control) == 0) return(NULL)
    tibble(gene = gene, time_min = t,
           log2fc = log2(mean(treated) + 1) - log2(mean(control) + 1))
  })
})

cat(sprintf("Genome-wide temporal profiles: %d genes\n", n_distinct(temporal_all$gene)))

# Build genome-wide profile matrix
profile_all_mat <- temporal_all %>%
  pivot_wider(names_from = time_min, values_from = log2fc) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Filter: must have variance
gene_var_all <- apply(profile_all_mat, 1, var, na.rm = TRUE)
profile_all_mat <- profile_all_mat[gene_var_all > 0.05, , drop = FALSE]
cat(sprintf("After variance filter: %d genes\n", nrow(profile_all_mat)))

# Mean FGF feedback profile
fgf_in_all <- intersect(fgf_feedback, rownames(profile_all_mat))
if (length(fgf_in_all) >= 3) {
  fgf_mean_profile <- colMeans(profile_all_mat[fgf_in_all, , drop = FALSE], na.rm = TRUE)

  # Correlate every gene with the FGF mean profile
  cor_with_fgf_all <- apply(profile_all_mat, 1, function(x) {
    cor(x, fgf_mean_profile, use = "pairwise.complete.obs")
  })

  # Screen criteria for novel feedback candidates:
  # 1. Anti-correlated with FGF (r < -0.5)
  # 2. Upregulated overall (mean log2fc > 0.3) — it's an activated gene, not just noise
  # 3. Late onset (peak at >= 120 min)
  gene_means <- rowMeans(profile_all_mat, na.rm = TRUE)
  gene_peaks <- apply(profile_all_mat, 1, function(x) {
    tp_names <- as.numeric(colnames(profile_all_mat))
    tp_names[which.max(abs(x))]
  })

  novel_screen <- data.frame(
    gene = names(cor_with_fgf_all),
    cor_fgf = cor_with_fgf_all,
    mean_fc = gene_means[names(cor_with_fgf_all)],
    peak_time = gene_peaks[names(cor_with_fgf_all)],
    stringsAsFactors = FALSE
  ) %>%
    filter(cor_fgf < -0.5,           # anti-correlated with FGF
           abs(mean_fc) > 0.3,       # biologically meaningful
           peak_time >= 120) %>%     # delayed response
    # Exclude genes already in our curated universe
    mutate(in_grn = gene %in% grn_available) %>%
    arrange(cor_fgf)

  cat(sprintf("\nNovel anti-FGF candidates (r < -0.5, |mean FC| > 0.3, peak >= 120 min): %d\n",
              nrow(novel_screen)))
  cat(sprintf("  Of which NOT in curated GRN: %d\n", sum(!novel_screen$in_grn)))

  write_csv(novel_screen, results_path("grn_novel_feedback_candidates.csv"))

  # Show top hits
  cat("\nTop 30 novel candidates:\n")
  novel_screen %>%
    filter(!in_grn) %>%
    head(30) %>%
    print(n = 30)
}

# ============================================================================
# 9. SUMMARY STATISTICS AND MODEL
# ============================================================================

cat("\n========== SECTION 9: SUMMARY ==========\n")

# Print the biological model supported by the data
cat("\n===================================================================\n")
cat("FGF–NODAL GENE REGULATORY NETWORK: SUMMARY MODEL\n")
cat("===================================================================\n\n")

cat("LITERATURE-BASED MODEL (van Boxtel et al. 2018, Dorey & Amaya 2010):\n")
cat("  1. Nodal (ndr1/ndr2) induces FGF ligand transcription (fgf3, fgf8a, fgf17)\n")
cat("  2. FGF activates MAPK/ERK cascade → ETS targets (etv4, etv5a, etv5b)\n")
cat("  3. FGF/MAPK activates immediate negative feedback: spry2, spry4, dusp6, il17rd\n")
cat("  4. At HIGH Nodal levels, Dusp4 is induced → intracellular FGF/ERK inhibition\n")
cat("  5. Sprouty/Spred switch: different feedback mechanisms at different FGF levels\n")
cat("  6. FGF/ERK promotes mesoderm (tbxta, tbx16) at expense of endoderm (sox32)\n\n")

cat("DATA-SUPPORTED FINDINGS:\n")

# Timing
if (nrow(timing_metrics) > 0) {
  fgf_onset <- median(timing_metrics$onset_time[timing_metrics$gene %in% fgf_feedback], na.rm = TRUE)
  nodal_onset <- median(timing_metrics$onset_time[timing_metrics$gene %in% tolower(nodal_core)], na.rm = TRUE)
  dusp4_onset <- timing_metrics$onset_time[timing_metrics$gene == "dusp4"]
  cat(sprintf("  - Nodal pathway median onset: %.0f min\n", nodal_onset))
  cat(sprintf("  - FGF feedback median onset: %.0f min\n", fgf_onset))
  if (length(dusp4_onset) > 0 && !is.na(dusp4_onset))
    cat(sprintf("  - DUSP4 onset: %.0f min (delayed = %s)\n",
                dusp4_onset, ifelse(dusp4_onset > fgf_onset, "YES", "NO")))
}

# Dose-response
if (exists("dose_pattern") && nrow(dose_pattern) > 0) {
  cat("\nDose-response patterns:\n")
  dp_summary <- dose_pattern %>%
    group_by(group, pattern) %>%
    tally() %>%
    filter(n > 0, pattern != "Other/weak") %>%
    arrange(group, desc(n))
  print(dp_summary, n = 40)
}

# Top candidates
cat("\nTOP 10 NOVEL NEGATIVE FEEDBACK CANDIDATES:\n")
if (exists("candidate_scores")) {
  candidate_scores %>%
    filter(!gene %in% fgf_feedback, !gene %in% tolower(nodal_core)) %>%
    head(10) %>%
    dplyr::select(gene, group, composite_score, cor_with_fgf,
                  pattern, chip_bound, peak_time) %>%
    print(n = 10, width = 200)
}

cat("\n========== GRN ANALYSIS COMPLETE ==========\n")
cat("\nOutput files (all in results/):\n")
cat("  grn_gene_universe.csv           — all genes with group annotations\n")
cat("  grn_dose_response_profiles.csv  — dose × time × gene matrix\n")
cat("  grn_dose_patterns.csv           — classification of dose-response shapes\n")
cat("  grn_temporal_15ngml.csv         — temporal profiles at 15 ng/ml\n")
cat("  grn_timing_metrics.csv          — onset time, peak time, AUC per gene\n")
cat("  grn_fgf_correlation_ranking.csv — correlation with FGF targets\n")
cat("  grn_candidate_scores.csv        — composite scoring of feedback candidates\n")
cat("  grn_novel_feedback_candidates.csv — genome-wide screen results\n")
cat("  grn_network_edges.csv           — network edges (literature + data)\n")
cat("\nFigures:\n")
cat("  grn_dose_response_heatmap.pdf   — dose-response heatmap by gene family\n")
cat("  grn_temporal_dynamics.pdf       — temporal dynamics by pathway\n")
cat("  grn_dose_response_curves.pdf    — individual gene dose curves\n")
cat("  grn_correlation_heatmap.pdf     — temporal correlation matrix\n")
cat("  grn_network_graph.pdf           — gene regulatory network\n")
cat("  grn_candidate_dotplot.pdf       — top candidate scoring\n")
cat("  grn_activation_timeline.pdf     — onset-to-peak gantt chart\n")
cat("  grn_dose_fgf_vs_nodal.pdf       — FGF vs Nodal dose curves\n")
