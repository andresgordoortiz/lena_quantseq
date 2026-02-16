# ============================================================================
# GENE REGULATORY NETWORK: FGF–NODAL SIGNALLING CROSSTALK
# ============================================================================
#
# PURPOSE
# -------
# This script investigates how Activin/Nodal signalling regulates FGF
# pathway activity over TIME.  The central biological question:
#
#   "Nodal activates FGF ligand transcription, yet prolonged Nodal
#    exposure leads to FGF/ERK inhibition — what genes mediate this
#    delayed negative feedback, and what is the temporal order?"
#
# The key insight from this dataset is that EXPOSURE TIME — not Activin
# concentration — drives the transcriptomic changes.  All analyses here
# therefore focus on the 15 ng/ml time course (15–240 min) vs 0 ng/ml.
#
# BIOLOGICAL MODEL (literature)
# ----
#   1. Nodal (ndr1/ndr2) → FGF ligand transcription (fgf3, fgf8a, fgf17)
#   2. FGF/MAPK/ERK cascade → ETS targets (etv4, etv5a, etv5b)
#   3. FGF/MAPK → immediate negative feedback (spry2, spry4, dusp6, il17rd)
#   4. Prolonged Nodal → induces DUSP4 → intracellular FGF/ERK block
#      (van Boxtel et al. 2018)
#   5. FGF/ERK promotes mesoderm (tbxta, tbx16) vs endoderm (sox32)
#
# ANALYTICAL STRATEGY
# ----
#   A. Validate the model using POSITIVE CONTROLS: genes whose roles in
#      this network are established in the literature.
#   B. Apply the same temporal framework to CANDIDATE genes to identify
#      novel feedback regulators.
#   C. Build a gene regulatory network combining literature edges with
#      data-driven correlation/anti-correlation patterns.
#
# REFERENCES
#   • van Boxtel et al. 2018 (Dev Cell 2018): Nodal induces DUSP4
#   • Dorey & Amaya 2010 (Annu Rev Cell Dev Biol): FGF/Nodal crosstalk
#   • Furthauer et al. 2001, 2002, 2004: Sprouty/Sef FGF feedback
#   • Tsang et al. 2004: DUSP6 as FGF feedback target
#   • Nelson et al. 2014: Smad2/EomesA ChIP-seq in zebrafish
#
# OUTPUTS (all in results/)
#   CSVs:  grn_gene_universe.csv, grn_temporal_profiles.csv,
#          grn_timing_metrics.csv, grn_fgf_correlation_ranking.csv,
#          grn_candidate_scores.csv, grn_novel_feedback_candidates.csv
#   PDFs:  grn_positive_controls.pdf

source("preprocess.R")
library(ggrepel)
library(patchwork)
library(scales)
library(readxl)


# ============================================================================
# SECTION 0 — GENE UNIVERSE
# ============================================================================
#
# We define genes in two categories:
#
#   POSITIVE CONTROLS  — genes with well-established roles in the FGF–Nodal
#     network, drawn from primary literature.  These VALIDATE the analysis:
#     if our temporal profiling recovers their known behaviour, we can trust
#     the framework for novel candidates.
#
#   CANDIDATES — genes hypothesised to mediate negative feedback (Nodal →
#     FGF inhibition) or related pathway crosstalk.  These are the DISCOVERY
#     targets of this analysis.

cat("\n========== SECTION 0: GENE UNIVERSE ==========\n")

# ---- Nodal/TGFβ core pathway (POSITIVE CONTROL) ----
# Nodal ligands (ndr1/ndr2/spaw), co-receptor (tdgf1), receptors,
# effectors (Smad2/3/4, FoxH1, Eomesa), and auto-feedback (Lefty1/2).
# These should show early, strong activation.
nodal_core <- c(
  "ndr1", "ndr2", "spaw",           # ligands

"tdgf1",                            # co-receptor (one-eyed pinhead)
  "acvr1ba", "acvr2aa", "acvr2ba",  # type I/II receptors
  "smad2", "smad3a", "smad3b", "smad4a",  # signal transducers
  "foxh1", "eomesa", "mixl1", "gata5",    # downstream TFs
  "lefty1", "lefty2", "lft1", "lft2"      # feedback inhibitors
)

# ---- FGF/MAPK feedback targets (POSITIVE CONTROL) ----
# Genes DIRECTLY induced by FGF/MAPK signalling.  Their transcription is
# abolished by FGF inhibition (SU5402 or dn-FGFR).  This group serves as
# the "FGF readout" — if these are active, FGF/ERK is signalling.
#
# etv4/etv5a/etv5b: PEA3-family ETS TFs (Raible & Brand 2001; Roehl &
#   Nusslein-Volhard 2001)
# spry1/spry2/spry4: Sprouty family — RTK/MAPK feedback inhibitors
#   (Furthauer et al. 2001; Hacohen et al. 1998)
# dusp6 (mkp3): MAPK phosphatase, immediate FGF target
#   (Tsang et al. 2004; Kawakami et al. 2003)
# il17rd (sef): FGF-pathway-specific feedback inhibitor
#   (Furthauer et al. 2002; Tsang et al. 2002)
fgf_feedback <- c("etv4", "etv5a", "etv5b",
                   "spry1", "spry2", "spry4",
                   "dusp6", "il17rd")

# ---- FGF ligands (POSITIVE CONTROL) ----
# Which FGF ligands does Nodal/Activin induce?
# Literature: fgf3, fgf8a, fgf17 are Nodal-induced.  We survey all
# zebrafish fgf family members to see which respond.
fgf_ligands <- c("fgf3", "fgf4", "fgf8a", "fgf8b",
                  "fgf10a", "fgf10b", "fgf16", "fgf17",
                  "fgf18a", "fgf18b", "fgf19",
                  "fgf20a", "fgf20b", "fgf24")

# ---- FGF receptors (POSITIVE CONTROL) ----
# Zebrafish FGFR genes.  Expression should be relatively stable (receptors
# are often constitutive), but worth profiling to detect regulation.
fgf_receptors <- c("fgfr1a", "fgfr1b", "fgfr2", "fgfr3", "fgfr4")

# ---- Posterior mesoderm / FGF-dependent genes (POSITIVE CONTROL) ----
# These mesodermal TFs REQUIRE FGF for expression.  They are NOT direct
# FGF transcriptional targets but need the FGF signalling environment.
# If the Nodal–DUSP4 axis shuts down FGF/ERK late, these should decline
# or plateau after initial activation.
#
# tbxta/tbx16/tbx6: T-box TFs (Griffin & Bhatt 2004; Kimelman 2006)
# msgn1: presomitic mesoderm (Yoo et al. 2003)
# her1/her7/hes6: hairy/E(spl), somitogenesis clock (Holley et al. 2002)
# mespaa/mespab: somite boundary (Sawada et al. 2000)
# ripply1/ripply2: segmentation repressors (Kawamura et al. 2005)
# noto: notochord homeobox (Talbot et al. 1995)
mesoderm_genes <- c("tbx16", "tbx6", "tbxta", "noto", "msgn1",
                     "her1", "her7", "hes6",
                     "mespaa", "mespab", "ripply1", "ripply2")

# ---- Endoderm markers (POSITIVE CONTROL) ----
# FGF/ERK opposes endoderm fate.  If FGF is active, these should be low;
# if DUSP4 shuts off FGF, these may rise late.
endoderm_genes <- c("sox32", "sox17", "gata5", "gata6",
                     "foxa2", "foxa3", "bon", "cas")

# ---- MAPK cascade components (SURVEY) ----
# Core signalling cascade transcripts.  These encode the signalling
# proteins: Ras → Raf → MEK → ERK.  Transcriptional changes are
# unlikely (post-translational regulation dominates), but we survey.
mapk_cascade <- c("hras", "kras", "nras",
                   "braf", "raf1a", "raf1b",
                   "map2k1", "map2k2a",
                   "mapk1", "mapk3")

# ---- DUSP family (CANDIDATES) ----
# Dual-specificity phosphatases.  dusp6 is a known immediate FGF target
# (positive control); the rest are surveyed to find novel feedback
# regulators like dusp4 (van Boxtel et al. 2018).
dusp_family <- c("dusp1", "dusp2", "dusp4", "dusp5", "dusp6", "dusp7",
                  "dusp10", "dusp11", "dusp14", "dusp16")

# ---- BMP crosstalk (CANDIDATES) ----
# FGF/ERK phosphorylates Smad1 linker region → inhibits BMP signalling.
# Surveying whether BMP pathway genes change in our time course.
bmp_crosstalk <- c("bmp2b", "bmp4", "bmp7a", "bmp7b",
                    "chrd", "nog1", "nog2",
                    "smad1", "smad5")

# ---- Wnt crosstalk (CANDIDATES) ----
# Wnt and FGF cooperate in mesoderm.  Including to detect crosstalk.
wnt_crosstalk <- c("wnt8a", "wnt3a", "wnt11f2",
                    "ctnnb1", "ctnnb2",
                    "axin1", "axin2", "dkk1b")

# ---- Nodal score genes (from collaborator literature) ----
nodal_score_data <- tryCatch(
  read_excel("docs/nodal-score-genes_complete.xlsx", skip = 1),
  error = function(e) { cat("Note: nodal score spreadsheet not found\n"); NULL }
)
nodal_score_genes <- if (!is.null(nodal_score_data)) {
  unique(tolower(na.omit(nodal_score_data$`Nodal score`)))
} else { character(0) }

# ---- Combine into master universe (deduplicate) ----
grn_universe <- unique(tolower(c(
  nodal_core, fgf_feedback, fgf_ligands, fgf_receptors,
  mesoderm_genes, endoderm_genes, mapk_cascade,
  dusp_family, bmp_crosstalk, wnt_crosstalk,
  nodal_score_genes
)))

# Map to genes actually present in the count matrix
grn_available <- grn_universe[grn_universe %in% tolower(rownames(counts_filtered))]
cat(sprintf("GRN gene universe: %d defined → %d present in dataset\n",
            length(grn_universe), length(grn_available)))

# Assign each gene to its primary pathway group.  Order of case_when
# defines priority when a gene appears in multiple lists (e.g. gata5
# in both nodal_core and endoderm_genes → assigned to Nodal pathway).
gene_groups <- data.frame(gene = grn_available, stringsAsFactors = FALSE) %>%
  mutate(
    group = case_when(
      gene %in% tolower(nodal_core)      ~ "Nodal pathway",
      gene %in% tolower(fgf_feedback)    ~ "FGF/MAPK feedback",
      gene %in% tolower(fgf_ligands)     ~ "FGF ligands",
      gene %in% tolower(fgf_receptors)   ~ "FGF receptors",
      gene %in% tolower(mapk_cascade)    ~ "MAPK cascade",
      gene %in% tolower(mesoderm_genes)  ~ "Mesoderm (FGF-dependent)",
      gene %in% tolower(endoderm_genes)  ~ "Endoderm markers",
      gene %in% tolower(dusp_family)     ~ "DUSP family",
      gene %in% tolower(bmp_crosstalk)   ~ "BMP crosstalk",
      gene %in% tolower(wnt_crosstalk)   ~ "Wnt crosstalk",
      gene %in% nodal_score_genes        ~ "Nodal target",
      TRUE                               ~ "Other"
    ),
    # Classify as literature-validated ("control") vs discovery ("candidate")
    evidence_class = case_when(
      group %in% c("Nodal pathway", "FGF/MAPK feedback", "FGF ligands",
                    "FGF receptors", "Mesoderm (FGF-dependent)",
                    "Endoderm markers") ~ "Positive control",
      TRUE ~ "Candidate / survey"
    )
  )

write_csv(gene_groups, results_path("grn_gene_universe.csv"))

cat("\nGenes per group:\n")
print(table(gene_groups$group))
cat("\nEvidence class:\n")
print(table(gene_groups$evidence_class))


# ============================================================================
# SECTION 1 — NORMALISE EXPERIMENT 1
# ============================================================================
#
# All temporal analyses use Experiment 1 (Activin dose × time course).
# We normalise counts using DESeq2 size factors, then extract samples
# at 15 ng/ml (treated) and 0 ng/ml (control) for each timepoint.
#
# The key comparison is always:
#   log2FC = log2(mean_treated + 1) − log2(mean_control + 1)
# at each timepoint, giving a temporal response curve per gene.

cat("\n========== SECTION 1: NORMALISE EXP1 ==========\n")

exp1_meta <- metadata %>% filter(experiment == "Exp1")
exp1_counts <- counts_filtered[, rownames(exp1_meta)]

# Fit a minimal DESeq2 model (intercept-only) to get size factors.
# We do NOT run differential expression here — just normalisation.
dds_exp1 <- DESeqDataSetFromMatrix(exp1_counts, exp1_meta, ~ 1)
dds_exp1 <- estimateSizeFactors(dds_exp1)
norm_exp1 <- counts(dds_exp1, normalized = TRUE)

# Identify the sample groups we need
exp1_15ngml <- exp1_meta %>% filter(concentration == "15ngml")
exp1_0ngml  <- exp1_meta %>% filter(concentration == "0ngml")
timepoints  <- sort(unique(exp1_15ngml$time_min))

cat("Timepoints:", paste(timepoints, collapse = ", "), "min\n")
cat("15 ng/ml samples:", nrow(exp1_15ngml), "\n")
cat("0 ng/ml  samples:", nrow(exp1_0ngml), "\n")


# ============================================================================
# SECTION 2 — TEMPORAL PROFILES (15 ng/ml vs 0 ng/ml)
# ============================================================================
#
# For every gene in our GRN universe, compute the log2 fold-change at each
# timepoint.  This is the CORE data structure for all downstream analyses.
#
# Biological rationale: exposure time determines the transcriptomic state.
# By comparing 15 ng/ml (Activin-stimulated) vs 0 ng/ml (unstimulated
# control) at matched timepoints, we isolate the Activin-driven response
# at each point in the time course.

cat("\n========== SECTION 2: TEMPORAL PROFILES ==========\n")

compute_temporal_profile <- function(genes, norm_mat, meta_treated, meta_ctrl, tp) {
  # For each gene and timepoint, compute the mean-based log2FC.
  # Returns a long-format tibble with columns: gene, time_min, log2fc.
  purrr::map_dfr(genes, function(gene) {
    row_idx <- which(tolower(rownames(norm_mat)) == gene)
    if (length(row_idx) == 0) return(NULL)   # gene not in count matrix
    purrr::map_dfr(tp, function(t) {
      treated_vals <- norm_mat[row_idx, rownames(meta_treated)[meta_treated$time_min == t]]
      control_vals <- norm_mat[row_idx, rownames(meta_ctrl)[meta_ctrl$time_min == t]]
      if (length(treated_vals) == 0 || length(control_vals) == 0) return(NULL)
      tibble(gene  = gene,
             time_min = t,
             log2fc   = log2(mean(treated_vals) + 1) - log2(mean(control_vals) + 1))
    })
  })
}

# Compute for all GRN genes
temporal_profiles <- compute_temporal_profile(
  grn_available, norm_exp1, exp1_15ngml, exp1_0ngml, timepoints
)

# Annotate with pathway group and evidence class
temporal_profiles <- temporal_profiles %>%
  left_join(gene_groups, by = "gene")

cat(sprintf("Temporal profiles computed: %d genes × %d timepoints\n",
            n_distinct(temporal_profiles$gene), length(timepoints)))

write_csv(temporal_profiles, results_path("grn_temporal_profiles.csv"))


# ============================================================================
# SECTION 3 — ACTIVATION TIMING METRICS
# ============================================================================
#
# For each gene, compute summary statistics that characterise ITS temporal
# behaviour:
#
#   peak_time    — timepoint of strongest response (max |log2FC|)
#   peak_log2fc  — the log2FC at that peak
#   onset_time   — first timepoint where |log2FC| > 1.0 (conservative
#                  threshold; 2-fold change)
#   direction    — overall direction (up vs down, by sign of sum)
#   auc          — signed area under the curve (captures cumulative effect)
#   late_vs_early — difference between late (≥120 min) and early (≤60 min)
#                   mean log2FC.  Positive = delayed responder.
#
# These metrics are used downstream for:
#   - The activation timeline figure (Section 7c)
#   - Candidate scoring: "delayed" genes that rise AFTER FGF targets are
#     potential negative feedback regulators.

cat("\n========== SECTION 3: TIMING METRICS ==========\n")

timing_metrics <- temporal_profiles %>%
  group_by(gene) %>%
  summarise(
    peak_time = time_min[which.max(abs(log2fc))],
    peak_log2fc = log2fc[which.max(abs(log2fc))],
    onset_time = {
      sig <- time_min[abs(log2fc) > 1.0]
      if (length(sig) > 0) min(sig) else NA_real_
    },
    direction = ifelse(sum(log2fc) >= 0, "up", "down"),
    auc = sum(log2fc),
    late_vs_early = mean(log2fc[time_min >= 120]) - mean(log2fc[time_min <= 60]),
    .groups = "drop"
  ) %>%
  left_join(gene_groups, by = "gene")

write_csv(timing_metrics, results_path("grn_timing_metrics.csv"))

# Report key timing results
cat(sprintf("Timing metrics for %d genes\n", nrow(timing_metrics)))
cat(sprintf("  Direction: %d up, %d down\n",
            sum(timing_metrics$direction == "up"),
            sum(timing_metrics$direction == "down")))
cat(sprintf("  Genes with onset ≤ 30 min: %d (early responders)\n",
            sum(timing_metrics$onset_time <= 30, na.rm = TRUE)))
cat(sprintf("  Genes with onset ≥ 120 min: %d (delayed responders)\n",
            sum(timing_metrics$onset_time >= 120, na.rm = TRUE)))

# Flag how FGF feedback genes behave (validation)
fgf_timing <- timing_metrics %>% filter(gene %in% fgf_feedback)
cat("\nPositive control check — FGF feedback genes:\n")
cat(sprintf("  Median onset: %.0f min\n", median(fgf_timing$onset_time, na.rm = TRUE)))
cat(sprintf("  Median peak:  %.0f min\n", median(fgf_timing$peak_time)))
print(fgf_timing %>% dplyr::select(gene, onset_time, peak_time, peak_log2fc, direction))


# ============================================================================
# SECTION 4 — CORRELATION ANALYSIS WITH STATISTICAL TESTING
# ============================================================================
#
# Strategy: build a gene × gene Pearson correlation matrix using the
# temporal profiles as feature vectors.  Each gene has 6 values
# (log2FC at 15, 30, 60, 120, 180, 240 min).
#
# STATISTICAL RIGOUR: with only 6 timepoints per profile, individual
# pairwise correlations can be spurious.  We therefore:
#   (a) Compute p-values for each correlation using cor.test()
#   (b) Apply Benjamini-Hochberg FDR correction across ALL pairwise tests
#   (c) Only consider correlations with FDR-adjusted p < 0.05
#
# Biological rationale:
#   — Genes co-regulated by the same pathway will have POSITIVE correlation
#     (e.g. Nodal targets all rise together)
#   — A candidate negative regulator of FGF should be ANTI-CORRELATED
#     with FGF readouts (it rises when they fall, or vice versa)
#   — A gene activated by Nodal that inhibits FGF should correlate
#     positively with Nodal genes AND negatively with FGF genes

cat("\n========== SECTION 4: CORRELATION NETWORK ==========\n")

# Reshape to gene × timepoint matrix
profile_mat <- temporal_profiles %>%
  dplyr::select(gene, time_min, log2fc) %>%
  pivot_wider(names_from = time_min, values_from = log2fc) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Remove genes with negligible variance (no temporal structure to correlate)
gene_var <- apply(profile_mat, 1, var, na.rm = TRUE)
profile_mat <- profile_mat[gene_var > 0.01, , drop = FALSE]
cat(sprintf("Genes with sufficient variance for correlation: %d / %d\n",
            nrow(profile_mat), n_distinct(temporal_profiles$gene)))

# Pearson correlation of temporal profiles
cor_mat <- cor(t(profile_mat), use = "pairwise.complete.obs", method = "pearson")

# ---- Compute pairwise p-values and apply FDR correction ----
# For each pair of genes, test H0: r = 0 with cor.test().
# Then correct for multiple testing using Benjamini-Hochberg.
n_genes_cor <- nrow(profile_mat)
gene_names_cor <- rownames(profile_mat)

# Pre-allocate p-value matrix (only upper triangle computed)
pval_mat <- matrix(NA, n_genes_cor, n_genes_cor,
                   dimnames = list(gene_names_cor, gene_names_cor))

# Compute p-values for all unique pairs
pair_indices <- which(upper.tri(cor_mat), arr.ind = TRUE)
pvals_raw <- numeric(nrow(pair_indices))

for (k in seq_len(nrow(pair_indices))) {
  i <- pair_indices[k, 1]
  j <- pair_indices[k, 2]
  x <- profile_mat[i, ]
  y <- profile_mat[j, ]
  complete <- complete.cases(x, y)
  if (sum(complete) >= 4) {
    ct <- cor.test(x[complete], y[complete], method = "pearson")
    pvals_raw[k] <- ct$p.value
  } else {
    pvals_raw[k] <- NA
  }
}

# FDR correction (Benjamini-Hochberg) across ALL pairwise tests
pvals_adj <- p.adjust(pvals_raw, method = "BH")

# Fill the symmetric p-value matrix
for (k in seq_len(nrow(pair_indices))) {
  i <- pair_indices[k, 1]
  j <- pair_indices[k, 2]
  pval_mat[i, j] <- pvals_adj[k]
  pval_mat[j, i] <- pvals_adj[k]
}

# Report FDR results
n_tested <- sum(!is.na(pvals_raw))
n_sig <- sum(pvals_adj < 0.05, na.rm = TRUE)
n_strong_sig <- sum(pvals_adj < 0.05 & abs(cor_mat[upper.tri(cor_mat)]) > 0.85,
                     na.rm = TRUE)
cat(sprintf("Correlation tests: %d pairs tested\n", n_tested))
cat(sprintf("  FDR < 0.05: %d significant correlations\n", n_sig))
cat(sprintf("  FDR < 0.05 AND |r| > 0.85: %d strong significant correlations\n",
            n_strong_sig))

# ---- Mean correlation with FGF feedback targets (FDR-filtered) ----
# For each gene, how similar is its temporal trajectory to the average
# FGF feedback response?  We only use FDR-significant correlations.
fgf_feedback_in <- intersect(fgf_feedback, rownames(cor_mat))
cat(sprintf("FGF feedback genes in correlation matrix: %d / %d\n",
            length(fgf_feedback_in), length(fgf_feedback)))

if (length(fgf_feedback_in) >= 3) {
  # For each gene, compute mean correlation with FGF feedback targets,
  # but set non-significant correlations (FDR >= 0.05) to 0
  cor_fgf_fdr <- cor_mat[, fgf_feedback_in, drop = FALSE]
  pval_fgf_fdr <- pval_mat[, fgf_feedback_in, drop = FALSE]
  cor_fgf_fdr[is.na(pval_fgf_fdr) | pval_fgf_fdr >= 0.05] <- 0

  cor_with_fgf <- rowMeans(cor_fgf_fdr, na.rm = TRUE)

  # Also store the raw (unfiltered) mean correlation for comparison
  cor_with_fgf_raw <- rowMeans(cor_mat[, fgf_feedback_in, drop = FALSE], na.rm = TRUE)

  # Count how many FGF targets each gene significantly correlates with
  n_sig_fgf <- rowSums(!is.na(pval_fgf_fdr) & pval_fgf_fdr < 0.05)

  fgf_cor_df <- tibble(
    gene = names(cor_with_fgf),
    cor_with_fgf = cor_with_fgf,
    cor_with_fgf_raw = cor_with_fgf_raw,
    n_sig_fgf_targets = n_sig_fgf[names(cor_with_fgf)]
  ) %>%
    left_join(gene_groups, by = "gene") %>%
    arrange(cor_with_fgf)

  cat("\n— Most ANTI-correlated with FGF readouts (FDR < 0.05):\n")
  fgf_cor_df %>%
    filter(!gene %in% fgf_feedback, n_sig_fgf_targets > 0) %>%
    head(15) %>%
    dplyr::select(gene, group, cor_with_fgf, n_sig_fgf_targets) %>%
    print(n = 15)

  cat("\n— Most POSITIVELY correlated with FGF readouts (FDR < 0.05):\n")
  fgf_cor_df %>%
    filter(!gene %in% fgf_feedback, n_sig_fgf_targets > 0) %>%
    arrange(desc(cor_with_fgf)) %>%
    head(15) %>%
    dplyr::select(gene, group, cor_with_fgf, n_sig_fgf_targets) %>%
    print(n = 15)

  write_csv(fgf_cor_df, results_path("grn_fgf_correlation_ranking.csv"))
}

# ---- Nodal ↔ FGF cross-correlation ----
# Do FGF targets track or oppose Nodal targets?
nodal_in <- intersect(tolower(nodal_core), rownames(cor_mat))
if (length(nodal_in) >= 2 && length(fgf_feedback_in) >= 2) {
  cat("\nNodal ↔ FGF feedback temporal correlations:\n")
  print(round(cor_mat[nodal_in, fgf_feedback_in, drop = FALSE], 3))
}

# ---- Mean correlation with Nodal pathway (FDR-filtered) ----
if (length(nodal_in) >= 2) {
  cor_nodal_fdr <- cor_mat[, nodal_in, drop = FALSE]
  pval_nodal_fdr <- pval_mat[, nodal_in, drop = FALSE]
  cor_nodal_fdr[is.na(pval_nodal_fdr) | pval_nodal_fdr >= 0.05] <- 0

  cor_with_nodal <- rowMeans(cor_nodal_fdr, na.rm = TRUE)
  nodal_cor_df <- tibble(
    gene = names(cor_with_nodal),
    cor_with_nodal = cor_with_nodal
  )
} else {
  nodal_cor_df <- tibble(gene = character(), cor_with_nodal = numeric())
}


# ============================================================================
# SECTION 5 — CANDIDATE NEGATIVE FEEDBACK SCORING
# ============================================================================
#
# We score each gene as a potential mediator of the Nodal → FGF inhibition
# feedback loop.  The ideal candidate would:
#
#   1. Be ANTI-CORRELATED with FGF targets (it rises as FGF/ERK declines)
#   2. Be POSITIVELY correlated with Nodal pathway (co-regulated by Nodal)
#   3. Have DELAYED onset (activates AFTER direct FGF targets)
#   4. Be a direct Nodal target (Smad2/EomesA ChIP-seq binding)
#   5. Show late upregulation (peaks ≥120 min, after initial FGF wave)
#
# Each component is scaled 0–1 and equally weighted.  The composite score
# ranks candidates for follow-up validation.

cat("\n========== SECTION 5: CANDIDATE SCORING ==========\n")

# ---- Load ChIP-seq data (Nelson et al. 2014) ----
smad2_chipseq <- tryCatch(
  read_excel("12915_2014_81_MOESM4_ESM.xlsx", skip = 3),
  error = function(e) { cat("  Smad2 ChIP-seq file not found\n"); NULL }
)
eomesa_chipseq <- tryCatch(
  read_excel("12915_2014_81_MOESM10_ESM.xlsx", skip = 2),
  error = function(e) { cat("  EomesA ChIP-seq file not found\n"); NULL }
)

# Helper: check if a gene has a ChIP-seq peak nearby
check_chip <- function(gene_name, chipseq_df) {
  if (is.null(chipseq_df)) return(FALSE)
  pattern <- paste0("(?i)(^|;)", gene_name)
  any(grepl(pattern, chipseq_df$`Proximal gene`, perl = TRUE))
}

# Build the scoring table
# Start from timing metrics, add correlation and ChIP data
fgf_median_onset <- median(
  timing_metrics$onset_time[timing_metrics$gene %in% fgf_feedback],
  na.rm = TRUE
)

candidate_scores <- timing_metrics %>%
  dplyr::select(gene, peak_time, peak_log2fc, onset_time, direction,
                auc, late_vs_early, group, evidence_class) %>%
  # Component 1: anti-correlation with FGF targets
  left_join(fgf_cor_df %>% dplyr::select(gene, cor_with_fgf), by = "gene") %>%
  # Component 2: positive correlation with Nodal pathway
  left_join(nodal_cor_df %>% dplyr::select(gene, cor_with_nodal), by = "gene") %>%
  # Component 3: delayed onset (later than FGF feedback median)
  mutate(delayed = !is.na(onset_time) & onset_time > fgf_median_onset) %>%
  # Component 4: Smad2/EomesA ChIP-seq binding
  rowwise() %>%
  mutate(
    smad2_bound  = check_chip(gene, smad2_chipseq),
    eomesa_bound = check_chip(gene, eomesa_chipseq),
    chip_bound   = smad2_bound | eomesa_bound
  ) %>%
  ungroup()

# ---- Scale each component 0–1 and compute composite score ----
scale01 <- function(x) {
  x[is.na(x)] <- 0
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / diff(rng)
}

candidate_scores <- candidate_scores %>%
  mutate(
    # (1) Anti-FGF: invert correlation so anti-correlated genes score high
    score_antifgf  = scale01(-replace_na(cor_with_fgf, 0)),
    # (2) Pro-Nodal: positive correlation with Nodal pathway
    score_pronodal = scale01(replace_na(cor_with_nodal, 0)),
    # (3) Delayed onset
    score_delayed  = ifelse(delayed, 1, 0),
    # (4) ChIP-seq bound (direct Nodal target)
    score_chip     = ifelse(chip_bound, 1, 0),
    # (5) Late upregulation: upregulated AND response shifts late
    score_uplate   = ifelse(direction == "up" & late_vs_early > 0.3, 1, 0),

    # Composite: equally weighted mean of 5 components
    composite_score = (score_antifgf + score_pronodal + score_delayed +
                        score_chip + score_uplate) / 5
  ) %>%
  arrange(desc(composite_score))

write_csv(candidate_scores, results_path("grn_candidate_scores.csv"))

cat(sprintf("\nFGF feedback median onset: %.0f min\n", fgf_median_onset))
cat("\n=== TOP 25 CANDIDATE NEGATIVE FEEDBACK REGULATORS ===\n")
candidate_scores %>%
  filter(!gene %in% fgf_feedback, !gene %in% tolower(nodal_core)) %>%
  head(25) %>%
  dplyr::select(gene, group, composite_score, cor_with_fgf, cor_with_nodal,
                delayed, chip_bound, peak_time, peak_log2fc) %>%
  print(n = 25, width = 200)


# ============================================================================
# SECTION 6 — GENOME-WIDE SCREEN
# ============================================================================
#
# Extend beyond the curated GRN universe: compute temporal profiles for
# ALL expressed genes and screen for the "late-onset, anti-FGF" signature.
#
# Criteria:
#   1. Anti-correlated with FGF feedback mean profile (r < −0.5)
#   2. Biologically meaningful response (|mean log2FC| > 0.5)
#   3. Delayed responder (peak ≥ 120 min)
#
# This detects novel genes NOT in our curated lists that may participate
# in the Nodal → FGF inhibition axis.

cat("\n========== SECTION 6: GENOME-WIDE SCREEN ==========\n")

# Compute temporal profiles genome-wide (15 vs 0 ng/ml)
all_gene_names <- tolower(rownames(norm_exp1))

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

# Build genome-wide profile matrix and filter for variance
profile_all_mat <- temporal_all %>%
  pivot_wider(names_from = time_min, values_from = log2fc) %>%
  column_to_rownames("gene") %>%
  as.matrix()

gene_var_all <- apply(profile_all_mat, 1, var, na.rm = TRUE)
profile_all_mat <- profile_all_mat[gene_var_all > 0.05, , drop = FALSE]
cat(sprintf("After variance filter: %d genes\n", nrow(profile_all_mat)))

# Correlate each gene with the FGF feedback mean temporal profile
fgf_in_all <- intersect(fgf_feedback, rownames(profile_all_mat))

if (length(fgf_in_all) >= 3) {
  fgf_mean_profile <- colMeans(profile_all_mat[fgf_in_all, , drop = FALSE], na.rm = TRUE)

  cor_with_fgf_all <- apply(profile_all_mat, 1, function(x) {
    cor(x, fgf_mean_profile, use = "pairwise.complete.obs")
  })

  gene_means <- rowMeans(profile_all_mat, na.rm = TRUE)
  gene_peaks <- apply(profile_all_mat, 1, function(x) {
    tp_names <- as.numeric(colnames(profile_all_mat))
    tp_names[which.max(abs(x))]
  })

  novel_screen <- tibble(
    gene = names(cor_with_fgf_all),
    cor_fgf = cor_with_fgf_all,
    mean_fc = gene_means[names(cor_with_fgf_all)],
    peak_time = gene_peaks[names(cor_with_fgf_all)]
  ) %>%
    filter(cor_fgf < -0.5,
           abs(mean_fc) > 0.5,
           peak_time >= 120) %>%
    mutate(in_grn = gene %in% grn_available) %>%
    arrange(cor_fgf)

  cat(sprintf("\nNovel anti-FGF candidates: %d total, %d outside curated GRN\n",
              nrow(novel_screen), sum(!novel_screen$in_grn)))

  write_csv(novel_screen, results_path("grn_novel_feedback_candidates.csv"))

  cat("\nTop 20 novel candidates (not in curated GRN):\n")
  novel_screen %>%
    filter(!in_grn) %>%
    head(20) %>%
    print(n = 20)
}


# ============================================================================
# SECTION 7 — FIGURES
# ============================================================================
#
# Plotting conventions (matching q3/q4 scripts):
#   - theme_minimal, base_size 10, base_family "Helvetica"
#   - Titles: bold, size 14, hjust 0.5
#   - Subtitles: grey50, size 10, hjust 0.5
#   - Gene labels: italic, geom_text_repel
#   - panel.grid.minor = element_blank()
#   - Consistent colour palette per gene group

cat("\n========== SECTION 7: FIGURES ==========\n")

# ---- Shared colour palette for gene groups ----
# Used across all figures for consistent visual identity.
group_palette <- c(
  "Nodal pathway"           = "#7570B3",
  "FGF/MAPK feedback"       = "#E66101",
  "FGF ligands"             = "#1B9E77",
  "DUSP family"             = "#D95F02",
  "Mesoderm (FGF-dependent)" = "#66A61E",
  "Endoderm markers"        = "#A6761D",
  "FGF receptors"           = "#984EA3",
  "MAPK cascade"            = "#FF7F00",
  "BMP crosstalk"           = "#999999",
  "Wnt crosstalk"           = "#1F78B4",
  "Candidate feedback"      = "#F781BF",
  "Nodal target"            = "#B15928",
  "Other"                   = "#CCCCCC"
)

# ---- Shared ggplot2 theme to ensure visual consistency ----
theme_grn <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50"),
      strip.text    = element_text(size = 10, face = "bold"),
      panel.grid.minor  = element_blank(),
      axis.title    = element_text(size = 10),
      legend.title  = element_text(size = 9, face = "bold"),
      legend.text   = element_text(size = 8)
    )
}


# ============================================================================
# 7a. POSITIVE CONTROL TEMPORAL DYNAMICS
# ============================================================================
#
# Purpose: validate the analysis framework by showing how genes with
# KNOWN roles in the FGF–Nodal network respond over time.  If these
# match literature expectations, we can trust the approach for candidates.
#
# Expected patterns:
#   Nodal pathway → early, strong activation (ndr1/2 feedback up, lefty up)
#   FGF ligands   → activated BY Nodal (fgf3/8a/17 should rise)
#   FGF feedback  → activated BY FGF/ERK (etv4, spry2, dusp6 etc.)
#   Mesoderm      → requires FGF (tbxta, tbx16 — should follow FGF)
#   Endoderm      → opposes FGF (sox32 — may rise late if FGF declines)

cat("\n--- Figure 7a: Pathway temporal dynamics ---\n")

control_groups <- c("Nodal pathway", "FGF/MAPK feedback", "FGF ligands",
                     "Mesoderm (FGF-dependent)", "Endoderm markers",
                     "DUSP family", "MAPK cascade",
                     "BMP crosstalk", "Wnt crosstalk")

control_data <- temporal_profiles %>%
  filter(group %in% control_groups) %>%
  mutate(group = factor(group, levels = control_groups))

# Drop empty groups
populated <- control_data %>%
  group_by(group) %>%
  summarise(n = n_distinct(gene), .groups = "drop") %>%
  filter(n > 0) %>% pull(group)
control_data <- control_data %>% filter(group %in% populated)

# Compute group means for thick summary lines
control_means <- control_data %>%
  group_by(group, time_min) %>%
  summarise(mean_fc = mean(log2fc), sd_fc = sd(log2fc), .groups = "drop")

# Gene labels at the rightmost timepoint
control_labels <- control_data %>%
  filter(time_min == max(time_min))

p_controls <- ggplot(control_data, aes(x = time_min, y = log2fc)) +
  # Zero reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60",
             linewidth = 0.3) +
  # Individual gene trajectories (thin, semi-transparent)
  geom_line(aes(group = gene, color = group), alpha = 0.3, linewidth = 0.4) +
  # Group mean (thick)
  geom_line(data = control_means, aes(y = mean_fc, color = group),
            linewidth = 1.3) +
  geom_point(data = control_means, aes(y = mean_fc, color = group), size = 2) +
  # Gene name labels at last timepoint
  geom_text_repel(data = control_labels,
                  aes(label = gene, color = group),
                  size = 2.2, fontface = "italic", direction = "y",
                  hjust = 0, nudge_x = 8, segment.size = 0.2,
                  max.overlaps = 30, show.legend = FALSE) +
  facet_wrap(~ group, ncol = 3, scales = "free_y") +
  scale_color_manual(values = group_palette, guide = "none") +
  scale_x_continuous(breaks = timepoints,
                     expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = "Exposure time (min)",
       y = expression(log[2]*"FC (15 vs 0 ng/ml Activin)"),
       title = "Known pathway genes: temporal dynamics",
       subtitle = paste0("Thin lines = individual genes | Thick lines = group mean\n",
                         "All groups from the curated gene universe")) +
  theme_grn()

ggsave(results_path("grn_positive_controls.pdf"), p_controls,
       width = 16, height = 18)
cat("Saved:", results_path("grn_positive_controls.pdf"), "\n")


# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========== GRN ANALYSIS COMPLETE ==========\n")
cat("\nOUTPUT FILES (results/):\n")
cat("  Tables:\n")
cat("    grn_gene_universe.csv              gene list + group annotations\n")
cat("    grn_temporal_profiles.csv          temporal profiles (15 vs 0 ng/ml)\n")
cat("    grn_timing_metrics.csv             onset, peak, AUC per gene\n")
cat("    grn_fgf_correlation_ranking.csv    correlation with FGF readouts\n")
cat("    grn_candidate_scores.csv           composite feedback scoring\n")
cat("    grn_novel_feedback_candidates.csv  genome-wide screen hits\n")
cat("  Figures:\n")
cat("    grn_positive_controls.pdf          known pathway dynamics\n")
cat("\n")
