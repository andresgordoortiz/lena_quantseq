# ============================================================================
# Q4 EXTENDED ANALYSIS
# ============================================================================
#
# Additional analyses requested:
#   1. Add DUSP4, DUSP6 to ChIP-seq heatmap
#   2. Add FGF ligands, DUSP family, mesoderm/FGF-responsive genes to
#      temporal expression patterns (like q4_chipseq_integration Fig 2)
#   3. Gene lists for reversibility figures
#   4. Integrated overview at 120min (in addition to 60min)
#   5. Cell migration DEGs temporal patterns
#   6. FGF signaling receptor genes + ChIP-seq + TF integration
#   7. FGF score from literature
#

source("preprocess.R")
library(readxl)
library(pheatmap)
library(ggrepel)
library(patchwork)
library(clusterProfiler)
library(org.Dr.eg.db)

# ============================================================================
# SETUP: Load ChIP-seq data and gene lists
# ============================================================================

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
    gene = gene_name, chip_factor = chip_factor,
    n_peaks = nrow(matches),
    has_foxh1_motif = any(!is.na(matches$`Foxh1 binding site`)),
    ndr1_upregulated = if (chip_factor == "Smad2") {
      any(grepl(paste0("(?i)", gene_name, "\\^"), matches$`Proximal gene`, perl = TRUE))
    } else NA
  )
}

# Original gene lists
candidate_genes <- c("rhov", "net1", "flrt3", "rnd1b", "dkk1b", "plekha5b",
                     "rasgef1ba", "efna1a", "osr1", "frmd4ba", "abi1b",
                     "snai1b", "snai1a", "kirrel3l", "jcada", "pfkfb3",
                     "prickle1b", "efnb2a")
nodal_score_data <- read_excel("docs/nodal-score-genes_complete.xlsx", skip = 1)
nodal_genes <- unique(tolower(na.omit(nodal_score_data$`Nodal score`)))

# ============================================================================
# GENE FAMILY DEFINITIONS (used in Sections 1 and 2)
# ============================================================================

# ---- FGF/MAPK feedback genes ----
# These are genes whose transcription is DIRECTLY induced by FGF/MAPK
# signaling and are validated as FGF readouts in zebrafish.
# Evidence criteria: expression lost upon FGF inhibition (SU5402 or
# dominant-negative FGFR), and/or shown to be direct MAPK/ETS targets.
#
# References:
#   etv4 (pea3), etv5a (erm), etv5b: PEA3-family ETS transcription factors,
#     direct FGF/MAPK transcriptional targets.
#     — Raible & Brand 2001 (Mech Dev 107:105)
#     — Roehl & Nüsslein-Volhard 2001 (Curr Biol 11:503)
#     — Furthauer et al. 2004 (Development 131:2407)
#   spry1, spry2, spry4: Sprouty RTK/MAPK feedback inhibitors, expression
#     induced by FGF signaling.
#     — Furthauer et al. 2001 (Development 128:2175)
#     — Hacohen et al. 1998 (Cell 92:253)
#     — Minowada et al. 1999 (Development 126:4465)
#   dusp6 (mkp3): MAPK phosphatase, FGF-induced negative feedback.
#     — Tsang et al. 2004 (Nat Cell Biol 6:981)
#     — Kawakami et al. 2003 (J Biol Chem 278:14157)
#   il17rd (sef): FGF pathway-specific feedback inhibitor.
#     — Furthauer et al. 2002 (Nat Cell Biol 4:170)
#     — Tsang et al. 2002 (Nat Cell Biol 4:165)
#
# NOT included here (placed elsewhere or removed):
#   spred1/2a/2b/3 — Ras/MAPK signalling inhibitors but their transcription
#     is NOT well-documented as FGF-dependent in zebrafish. They function
#     constitutively as signalling modulators (Wakioka et al. 2001, Nature 412:647).
#   flrt3 — FGFR-interacting leucine-rich repeat protein. Functions primarily
#     in cell adhesion/repulsion and is already in the candidate gene list.
#     Its transcriptional regulation is Nodal/Activin-dependent, not FGF-dependent
#     (Böttcher et al. 2004, J Cell Biol 166:313).
#
fgf_feedback_genes <- c("etv4", "etv5a", "etv5b", "spry1", "spry2", "spry4",
                         "dusp6", "il17rd")

# ---- FGF ligands ----
# Zebrafish FGF gene family members. These are the signalling molecules,
# not downstream targets. Included to see which FGF ligands respond to Activin.
# Gene names follow ZFIN nomenclature.
fgf_ligands <- c("fgf3", "fgf4", "fgf8a", "fgf8b", "fgf10a", "fgf10b",
                  "fgf16", "fgf17", "fgf18a", "fgf18b", "fgf19",
                  "fgf20a", "fgf20b", "fgf24")

# ---- DUSP family ----
# Dual-specificity phosphatases that dephosphorylate MAPK family kinases.
# Only dusp6 is a validated direct FGF target (see above); the others are
# included to survey which DUSP family members respond to Activin.
dusp_family <- c("dusp1", "dusp2", "dusp4", "dusp5", "dusp6", "dusp7",
                  "dusp10", "dusp11", "dusp14", "dusp16")

# ---- Posterior mesoderm / FGF-dependent genes ----
# These are mesodermal TFs whose expression REQUIRES FGF signalling in the
# tailbud and presomitic mesoderm, but they are NOT direct FGF transcriptional
# targets. FGF maintains the permissive signalling environment (via Wnt and
# T-box loops) needed for their expression.
#
# References:
#   tbxta (ntl/brachyury), tbx16 (spadetail), tbx6:
#     — Griffin & Bhatt 2004 (Development 131:4981)
#     — Kimelman 2006 (Curr Opin Genet Dev 16:422)
#   msgn1 (mesogenin): presomitic mesoderm specification.
#     — Yoo et al. 2003 (Genes Dev 17:2385)
#   her1, her7, hes6: hairy/enhancer-of-split, somitogenesis clock.
#     — Holley et al. 2002 (Development 129:1175)
#   mespaa, mespab: Mesp family, somite boundary formation.
#     — Sawada et al. 2000 (Mech Dev 91:293)
#   ripply1, ripply2: segmentation repressors.
#     — Kawamura et al. 2005 (Development 132:4385)
#   noto: notochord homeobox.
#     — Talbot et al. 1995 (Nature 378:150)
#
mesoderm_fgf_dependent <- c("tbx16", "tbx6", "tbxta", "noto", "msgn1",
                             "her1", "her7", "hes6", "mespaa", "mespab",
                             "ripply1", "ripply2")

# ---- FGF receptors ----
# All five zebrafish FGFR genes (ZFIN nomenclature).
fgf_receptors <- c("fgfr1a", "fgfr1b", "fgfr2", "fgfr3", "fgfr4")

# ============================================================================
# 1. CHIPSEQ INTEGRATION HEATMAPS (2 focused heatmaps)
# ============================================================================

cat("\n========== ChIP-seq HEATMAPS ==========\n")

# Extended gene list
extra_genes <- c("dusp4", "dusp6")
all_genes_chip <- unique(c(candidate_genes, nodal_genes, extra_genes,
                           fgf_feedback_genes, fgf_ligands, dusp_family,
                           fgf_receptors))

chipseq_results <- bind_rows(
  map_dfr(all_genes_chip, ~ search_chipseq(.x, smad2_chipseq, "Smad2")),
  map_dfr(all_genes_chip, ~ search_chipseq(.x, eomesa_chipseq, "EomesA"))
)

make_chip_summary <- function(genes) {
  chipseq_results %>%
    filter(gene %in% genes) %>%
    group_by(gene, chip_factor) %>%
    summarise(n_peaks = sum(n_peaks), has_foxh1_motif = any(has_foxh1_motif),
              ndr1_upregulated = any(ndr1_upregulated, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(id_cols = gene, names_from = chip_factor,
                values_from = c(n_peaks, has_foxh1_motif, ndr1_upregulated),
                names_sep = "_") %>%
    mutate(
      Smad2_bound = n_peaks_Smad2 > 0,
      EomesA_bound = n_peaks_EomesA > 0,
      Foxh1_at_Smad2 = has_foxh1_motif_Smad2,
      Foxh1_at_EomesA = has_foxh1_motif_EomesA
    ) %>%
    dplyr::select(gene, Smad2_bound, n_peaks_Smad2, Foxh1_at_Smad2,
                  EomesA_bound, n_peaks_EomesA, Foxh1_at_EomesA) %>%
    arrange(gene)
}

# ---- Heatmap 1: Nodal score + Candidate genes ----
chip_hm1 <- bind_rows(
  make_chip_summary(nodal_genes) %>% mutate(gene_group = "Nodal score"),
  make_chip_summary(candidate_genes) %>% mutate(gene_group = "Candidate")
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
    gene_group = factor(gene_group, levels = c("Nodal score", "Candidate"))
  )

p_chip_hm1 <- ggplot(chip_hm1, aes(x = feature, y = gene)) +
  geom_tile(aes(fill = present), color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#F0F0F0"),
                    labels = c("TRUE" = "Bound", "FALSE" = "Not bound"), name = "") +
  facet_grid(gene_group ~ ., scales = "free_y", space = "free_y") +
  labs(x = "", y = "",
       title = "ChIP-seq binding: Nodal score & Candidate genes",
       subtitle = "Smad2 & EomesA peaks ± Foxh1 motif (Nelson et al. 2014)") +
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

ggsave(results_path("q4_chipseq_nodal_candidates.pdf"), p_chip_hm1,
       width = 6, height = 12)
cat("Saved:", results_path("q4_chipseq_nodal_candidates.pdf"), "\n")

# ---- Heatmap 2: DUSP family + FGF pathway ----
chip_hm2 <- bind_rows(
  make_chip_summary(dusp_family) %>% mutate(gene_group = "DUSP family"),
  make_chip_summary(fgf_feedback_genes) %>% mutate(gene_group = "FGF/MAPK feedback"),
  make_chip_summary(fgf_ligands) %>% mutate(gene_group = "FGF ligands"),
  make_chip_summary(fgf_receptors) %>% mutate(gene_group = "FGF receptors")
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
    gene_group = factor(gene_group, levels = c("DUSP family", "FGF/MAPK feedback",
                                                "FGF ligands", "FGF receptors"))
  )

p_chip_hm2 <- ggplot(chip_hm2, aes(x = feature, y = gene)) +
  geom_tile(aes(fill = present), color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("TRUE" = "#D95F02", "FALSE" = "#F0F0F0"),
                    labels = c("TRUE" = "Bound", "FALSE" = "Not bound"), name = "") +
  facet_grid(gene_group ~ ., scales = "free_y", space = "free_y") +
  labs(x = "", y = "",
       title = "ChIP-seq binding: DUSP & FGF pathway genes",
       subtitle = "Smad2 & EomesA peaks ± Foxh1 motif (Nelson et al. 2014)") +
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

ggsave(results_path("q4_chipseq_dusp_fgf.pdf"), p_chip_hm2,
       width = 6, height = 14)
cat("Saved:", results_path("q4_chipseq_dusp_fgf.pdf"), "\n")

# ============================================================================
# 2. FGF/DUSP/MESODERM TEMPORAL EXPRESSION PATTERNS (Exp1)
# ============================================================================

cat("\n========== FGF/DUSP/MESODERM TEMPORAL PATTERNS ==========\n")

# Gene families already defined above (before Section 1)

# Normalize Exp1
exp1_meta <- metadata %>% filter(experiment == "Exp1")
exp1_counts <- counts_filtered[, rownames(exp1_meta)]
dds_exp1 <- DESeqDataSetFromMatrix(exp1_counts, exp1_meta, ~ 1)
dds_exp1 <- estimateSizeFactors(dds_exp1)
norm_exp1 <- counts(dds_exp1, normalized = TRUE)

exp1_15ngml <- exp1_meta %>% filter(concentration == "15ngml")
exp1_0ngml  <- exp1_meta %>% filter(concentration == "0ngml")
timepoints <- sort(unique(exp1_15ngml$time_min))

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

# Compute profiles for each family
gene_families <- list(
  "FGF/MAPK feedback" = fgf_feedback_genes,
  "FGF Ligands" = fgf_ligands,
  "DUSP Family" = dusp_family,
  "Mesoderm (FGF-dependent)" = mesoderm_fgf_dependent,
  "FGF Receptors" = fgf_receptors
)

all_family_profiles <- map_dfr(names(gene_families), function(family) {
  compute_temporal_profile(gene_families[[family]], norm_exp1,
                           exp1_15ngml, exp1_0ngml, timepoints) %>%
    mutate(gene_family = family)
})

cat("Genes found per family:\n")
all_family_profiles %>%
  group_by(gene_family) %>%
  summarise(n_genes = n_distinct(gene), .groups = "drop") %>%
  print()

write_csv(all_family_profiles, results_path("q4_family_temporal_profiles.csv"))

# Combined plot — one facet per family
family_means <- all_family_profiles %>%
  group_by(gene_family, time_min) %>%
  summarise(mean_fc = mean(log2fc), sd_fc = sd(log2fc), .groups = "drop")

gene_endpoint_labels <- all_family_profiles %>%
  filter(time_min == max(time_min))

family_colors <- c(
  "FGF/MAPK feedback" = "#E66101",
  "FGF Ligands" = "#1B9E77",
  "DUSP Family" = "#D95F02",
  "Mesoderm (FGF-dependent)" = "#7570B3",
  "FGF Receptors" = "#E7298A"
)

p_families <- ggplot(all_family_profiles, aes(x = time_min, y = log2fc)) +
  geom_line(aes(group = gene, color = gene_family), alpha = 0.35, linewidth = 0.4) +
  geom_line(data = family_means, aes(y = mean_fc, color = gene_family),
            linewidth = 1.2) +
  geom_point(data = family_means, aes(y = mean_fc, color = gene_family), size = 2) +
  geom_text_repel(data = gene_endpoint_labels,
                  aes(label = gene, color = gene_family),
                  size = 2.2, fontface = "italic", direction = "y",
                  hjust = 0, nudge_x = 8, segment.size = 0.2,
                  max.overlaps = 30, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  facet_wrap(~ gene_family, ncol = 3, scales = "free_y") +
  scale_color_manual(values = family_colors, guide = "none") +
  scale_x_continuous(breaks = timepoints, expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = "Time (min)", y = "log2FC (15 vs 0 ng/ml Activin)",
       title = "Gene family temporal dynamics in Experiment 1",
       subtitle = "Thin lines: individual genes | Thick lines: family mean") +
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50")
  )

ggsave(results_path("q4_gene_families_temporal.pdf"), p_families, width = 16, height = 10)
cat("Saved:", results_path("q4_gene_families_temporal.pdf"), "\n")

# ============================================================================
# 3. GENE LISTS FOR REVERSIBILITY FIGURE
# ============================================================================

cat("\n========== REVERSIBILITY GENE LISTS ==========\n")

# Load reversibility scores if they exist; otherwise compute
if (file.exists(results_path("q4_reversibility_scores.csv"))) {
  rev_all <- read_csv(results_path("q4_reversibility_scores.csv"), show_col_types = FALSE)
} else {
  cat("Reversibility scores not found; run q4_chipseq_integration.R first.\n")
  rev_all <- NULL
}

if (!is.null(rev_all)) {
  # Classify by reversibility at 60min
  rev_classified <- rev_all %>%
    mutate(
      rev_class_60 = case_when(
        is.na(rev_60) ~ "Insufficient effect",
        rev_60 > 0.8 ~ "Fully reversed",
        rev_60 > 0.3 ~ "Partially reversed",
        rev_60 > -0.1 ~ "Not reversed",
        TRUE ~ "Over-reversed"
      ),
      rev_class_120 = case_when(
        is.na(rev_120) ~ "Insufficient effect",
        rev_120 > 0.8 ~ "Fully reversed",
        rev_120 > 0.3 ~ "Partially reversed",
        rev_120 > -0.1 ~ "Not reversed",
        TRUE ~ "Over-reversed"
      )
    )

  write_csv(rev_classified, results_path("q4_reversibility_classified.csv"))
  cat("Saved:", results_path("q4_reversibility_classified.csv"), "\n")

  cat("\nReversibility at SB50 60min:\n")
  print(table(rev_classified$rev_class_60))
  cat("\nReversibility at SB50 120min:\n")
  print(table(rev_classified$rev_class_120))
}

# ============================================================================
# 4. INTEGRATED OVERVIEW AT 120MIN (in addition to 60min)
# ============================================================================

cat("\n========== INTEGRATED OVERVIEW AT 120min ==========\n")

if (!is.null(rev_all) && "rev_120" %in% colnames(rev_all)) {
  # Load combined summary
  if (file.exists(results_path("q4_combined_summary.csv"))) {
    combined_summary <- read_csv(results_path("q4_combined_summary.csv"), show_col_types = FALSE)
  } else {
    cat("Combined summary not found; run q4_chipseq_integration.R first.\n")
    combined_summary <- NULL
  }

  if (!is.null(combined_summary)) {
    # Add rev_120 data
    combined_120 <- combined_summary %>%
      left_join(rev_all %>% dplyr::select(gene, rev_120), by = "gene") %>%
      mutate(rev_120_clamped = pmin(pmax(rev_120, 0), 1))

    chip_fill_colors <- c("Smad2 + EomesA" = "#7570B3", "Smad2 only" = "#D95F02",
                           "EomesA only" = "#1B9E77", "No binding" = "#E0E0E0")

    make_overview_120 <- function(data, title_text, title_color) {
      ggplot(data %>% filter(!is.na(rev_120_clamped)),
             aes(x = peak_time, y = rev_120_clamped, size = abs(peak_log2fc))) +
        geom_point(aes(fill = chip_label), shape = 21,
                   alpha = 0.85, stroke = 0.5, color = "black") +
        geom_text_repel(aes(label = gene), size = 3.2, fontface = "italic",
                        color = "#333333", max.overlaps = 25, seed = 42) +
        geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
        scale_fill_manual(values = chip_fill_colors, name = "ChIP-seq binding") +
        scale_size_continuous(range = c(2, 8), name = "|Peak log2FC|") +
        scale_x_continuous(breaks = c(60, 120, 180, 240)) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
        labs(x = "Time of peak response (min)", y = "Reversibility (SB50 @ 120 min)",
             title = title_text) +
        coord_cartesian(ylim = c(0, 1)) +
        theme_minimal(base_size = 10) +
        theme(
          panel.grid.minor = element_blank(),
          legend.position = "right",
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, color = title_color)
        )
    }

    group_colors_q4 <- c("Candidate" = "#E66101", "Nodal score" = "#5E3C99")

    p_120_cand <- make_overview_120(
      combined_120 %>% filter(gene_group == "Candidate"), "Candidate genes", "#E66101")
    p_120_nodal <- make_overview_120(
      combined_120 %>% filter(gene_group == "Nodal score"), "Nodal score genes", "#5E3C99")

    fig_120 <- p_120_cand / p_120_nodal +
      plot_annotation(
        title = "Integrated view: response timing vs reversibility (SB50 @ 120 min)",
        subtitle = "Bubble size = peak fold-change | Fill = ChIP-seq binding",
        theme = theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50")
        )
      )

    ggsave(results_path("q4_integrated_overview_120min.pdf"), fig_120, width = 10, height = 14)
    cat("Saved:", results_path("q4_integrated_overview_120min.pdf"), "\n")
  }
}

# ============================================================================
# 5. CELL MIGRATION & CELL ADHESION DEGs — TEMPORAL PATTERNS (2 SEPARATE PLOTS)
# ============================================================================
#
# Gene selection: intersect GO-annotated genes with shared DEGs from Q1.
# Two separate GO categories produce two separate figures:
#   - Cell migration:  GO:0016477 (cell migration), GO:0030334 (regulation
#                      of cell migration), GO:0040011 (locomotion)
#   - Cell adhesion:   GO:0007155 (cell adhesion)
#

cat("\n========== CELL MIGRATION & ADHESION DEGs ==========\n")

# Helper: query GO terms → gene symbols → intersect with shared DEGs
get_go_degs <- function(go_ids, label) {
  entrez <- tryCatch({
    AnnotationDbi::select(org.Dr.eg.db, keys = go_ids,
                           columns = "ENTREZID", keytype = "GOALL")$ENTREZID
  }, error = function(e) character(0))
  if (length(entrez) == 0) return(character(0))
  symbols <- AnnotationDbi::select(org.Dr.eg.db,
    keys = unique(entrez), columns = "SYMBOL", keytype = "ENTREZID")$SYMBOL
  symbols <- unique(na.omit(symbols))

  if (file.exists(results_path("q1_shared_de_genes.csv"))) {
    shared_degs <- read_csv(results_path("q1_shared_de_genes.csv"), show_col_types = FALSE)
    degs <- intersect(shared_degs$gene, symbols)
  } else {
    degs <- intersect(symbols, rownames(counts_filtered))
  }
  cat(sprintf("  %s: %d GO genes in dataset, %d shared DEGs\n",
              label, length(symbols), length(degs)))
  degs
}

migration_degs <- get_go_degs(
  c("GO:0016477", "GO:0030334", "GO:0040011"),  # migration + regulation + locomotion
  "Cell migration"
)
adhesion_degs <- get_go_degs(
  c("GO:0007155"),  # cell adhesion
  "Cell adhesion"
)

# Remove adhesion genes that overlap with migration (keep them clean)
adhesion_only <- setdiff(adhesion_degs, migration_degs)
cat(sprintf("  Adhesion-only (excluding migration overlap): %d\n", length(adhesion_only)))

# — Plotting helper (reusable for both categories) —
make_temporal_plot <- function(deg_list, label, color_up, color_down,
                               out_pdf, out_csv) {
  if (length(deg_list) < 3) {
    cat(sprintf("  Skipping %s plot: only %d genes\n", label, length(deg_list)))
    return(invisible(NULL))
  }

  profiles <- compute_temporal_profile(
    tolower(deg_list), norm_exp1, exp1_15ngml, exp1_0ngml, timepoints
  ) %>% mutate(gene_family = label)

  write_csv(profiles, out_csv)

  direction_df <- profiles %>%
    group_by(gene) %>%
    summarise(mean_fc = mean(log2fc), .groups = "drop") %>%
    mutate(direction = ifelse(mean_fc >= 0, "Upregulated", "Downregulated"))

  profiles <- profiles %>%
    left_join(direction_df %>% dplyr::select(gene, direction), by = "gene")

  dir_means <- profiles %>%
    group_by(direction, time_min) %>%
    summarise(mean_fc = mean(log2fc), .groups = "drop")

  endpoints <- profiles %>% filter(time_min == max(time_min))
  n_up   <- sum(direction_df$direction == "Upregulated")
  n_down <- sum(direction_df$direction == "Downregulated")

  dir_colors <- c("Upregulated" = color_up, "Downregulated" = color_down)

  p <- ggplot(profiles, aes(x = time_min, y = log2fc)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_line(aes(group = gene, color = direction), alpha = 0.25, linewidth = 0.3) +
    geom_line(data = dir_means, aes(y = mean_fc, color = direction),
              linewidth = 1.4) +
    geom_point(data = dir_means, aes(y = mean_fc, color = direction), size = 2.5) +
    geom_text_repel(data = endpoints,
                    aes(label = gene, color = direction), size = 2,
                    fontface = "italic", direction = "y", hjust = 0,
                    nudge_x = 8, segment.size = 0.2, max.overlaps = 30,
                    show.legend = FALSE) +
    facet_wrap(~ direction, ncol = 2, scales = "free_y") +
    scale_color_manual(values = dir_colors, guide = "none") +
    scale_x_continuous(breaks = timepoints, expand = expansion(mult = c(0.05, 0.25))) +
    labs(x = "Time (min)", y = "log2FC (15 vs 0 ng/ml Activin)",
         title = paste0(label, " DEGs — temporal response"),
         subtitle = sprintf("n = %d shared DEGs (%d up, %d down)",
                            length(deg_list), n_up, n_down)) +
    theme_minimal(base_size = 11, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.6),
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50")
    )

  ggsave(out_pdf, p, width = 12, height = 7)
  cat("Saved:", out_pdf, "\n")
  invisible(p)
}

# Plot 1: Cell migration
make_temporal_plot(
  migration_degs, "Cell migration",
  color_up = "#D73027", color_down = "#4575B4",
  out_pdf = results_path("q4_migration_degs_temporal.pdf"),
  out_csv = results_path("q4_migration_degs_temporal.csv")
)

# Plot 2: Cell adhesion
make_temporal_plot(
  adhesion_only, "Cell adhesion",
  color_up = "#E66101", color_down = "#5E3C99",
  out_pdf = results_path("q4_adhesion_degs_temporal.pdf"),
  out_csv = results_path("q4_adhesion_degs_temporal.csv")
)

# Gene lists for reference
write_csv(data.frame(gene = migration_degs), results_path("q4_migration_degs_list.csv"))
write_csv(data.frame(gene = adhesion_only), results_path("q4_adhesion_degs_list.csv"))

# NOTE: Section 6 (FGF signaling ChIP-seq heatmap) was removed.
# Its content is now fully covered by the combined DUSP+FGF heatmap
# produced in Section 1 (q4_chipseq_dusp_fgf.pdf).

# ============================================================================
# 6. FGF SCORE FROM LITERATURE
# ============================================================================
#
# Literature-based FGF score genes (Furthauer et al. 2004, Tsang et al. 2004,
# Sivak et al. 2005, Roehl & Nusslein-Volhard 2001):
#   Core readout: etv4 (pea3), etv5a (erm), etv5b, dusp6 (mkp3)
#   Negative feedback: spry2, spry4, il17rd (sef)
#   These are the most validated direct FGF/MAPK targets in zebrafish
#

cat("\n========== FGF SCORE (LITERATURE-BASED) ==========\n")

fgf_score_genes <- c("etv4", "etv5a", "etv5b", "dusp6", "spry2", "spry4", "il17rd")

cat("FGF score genes (literature):", paste(fgf_score_genes, collapse = ", "), "\n")

# Check availability
fgf_in_data <- fgf_score_genes[tolower(fgf_score_genes) %in% tolower(rownames(norm_exp1))]
cat(sprintf("Found %d / %d in dataset\n", length(fgf_in_data), length(fgf_score_genes)))

# Temporal profiles
fgf_score_profiles <- compute_temporal_profile(
  fgf_in_data, norm_exp1, exp1_15ngml, exp1_0ngml, timepoints
) %>% mutate(gene_family = "FGF Score (literature)")

write_csv(fgf_score_profiles, results_path("q4_fgf_score_literature_temporal.csv"))

cat("\n========== Q4 EXTENDED ANALYSIS COMPLETE ==========\n")
cat("\nOutput files (all in results/):\n")
cat("  - q4_chipseq_nodal_candidates.pdf (Nodal + Candidate genes ChIP-seq)\n")
cat("  - q4_chipseq_dusp_fgf.pdf (DUSP + FGF pathway ChIP-seq)\n")
cat("  - q4_gene_families_temporal.pdf (FGF feedback/DUSP/mesoderm/ligands/receptors)\n")
cat("  - q4_family_temporal_profiles.csv\n")
cat("  - q4_reversibility_classified.csv\n")
cat("  - q4_integrated_overview_120min.pdf\n")
cat("  - q4_migration_degs_temporal.pdf + csv (cell migration DEGs)\n")
cat("  - q4_adhesion_degs_temporal.pdf + csv (cell adhesion DEGs)\n")
cat("  - q4_fgf_score_literature_temporal.csv\n")
