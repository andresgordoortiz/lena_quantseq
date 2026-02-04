# GO Molecular Function Analysis - Publication-Ready Plots
# Uses styling from exp1_maria_reanalysis.R and q3_sb50_comparison.R

library(tidyverse)
library(ggplot2)
library(openxlsx)
library(patchwork)

# ============================================================================
# COLOR PALETTES
# ============================================================================

concentration_colors <- c(
  "5ngml" = "#7FBF7B",
  "10ngml" = "#4DAF4A",
  "15ngml" = "#1B7837"
)

# ============================================================================
# LOAD DATA
# ============================================================================

if (!file.exists("exp1_go_mf_results.xlsx") || !file.exists("exp2_go_mf_results.xlsx")) {
  stop("Please run go_analysis_simple.R first to generate the GO results")
}

exp1_go_df <- read.xlsx("exp1_go_mf_results.xlsx", sheet = 1)
exp2_go_df <- read.xlsx("exp2_go_mf_results.xlsx", sheet = 1)

cat("Loaded", nrow(exp1_go_df), "Exp1 GO terms and", nrow(exp2_go_df), "Exp2 GO terms\n")

# ============================================================================
# EXPERIMENT 1: PUBLICATION-READY PLOT
# ============================================================================

if (nrow(exp1_go_df) > 0) {

  # Filter to only 5, 10, 15 ng/ml
  exp1_filtered <- exp1_go_df %>%
    filter(concentration %in% c("5ngml", "10ngml", "15ngml"))

  top_terms <- exp1_filtered %>%
    group_by(Description) %>%
    summarize(min_p = min(p.adjust), n_conditions = n(), .groups = "drop") %>%
    arrange(min_p) %>%
    slice_head(n = 15) %>%
    pull(Description)

  plot_df <- exp1_filtered %>%
    filter(Description %in% top_terms) %>%
    mutate(
      neg_log_p = -log10(p.adjust),
      concentration = factor(concentration, levels = c("5ngml", "10ngml", "15ngml")),
      time_min = factor(time_min, levels = c(60, 120, 180, 240)),
      Description_short = str_wrap(Description, width = 40)
    )

  term_order <- plot_df %>%
    group_by(Description_short) %>%
    summarize(mean_neg_log_p = mean(neg_log_p), .groups = "drop") %>%
    arrange(mean_neg_log_p) %>%
    pull(Description_short)

  plot_df$Description_short <- factor(plot_df$Description_short, levels = term_order)

  # -------------------------------------------------------------------------
  # Combined single panel with colored backgrounds
  # -------------------------------------------------------------------------

  # Create combined x-axis variable
  plot_df_combined <- plot_df %>%
    mutate(
      x_var = interaction(time_min, concentration, sep = "_")
    )

  x_order <- expand.grid(
    time = c("60", "120", "180", "240"),
    conc = c("5ngml", "10ngml", "15ngml")
  ) %>%
    mutate(x_var = paste(time, conc, sep = "_")) %>%
    pull(x_var)

  plot_df_combined$x_var <- factor(plot_df_combined$x_var, levels = x_order)
  x_labels <- rep(c("60", "120", "180", "240"), 3)

  # Calculate p-value range for color scale
  p_range <- range(plot_df_combined$neg_log_p, na.rm = TRUE)
  p_mid <- mean(p_range)

  p_exp1_combined <- ggplot(plot_df_combined, aes(x = x_var, y = Description_short)) +
    # Colored background rectangles for concentration groups
    annotate("rect", xmin = 0.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
             fill = "#7FBF7B", alpha = 0.18) +
    annotate("rect", xmin = 4.5, xmax = 8.5, ymin = -Inf, ymax = Inf,
             fill = "#4DAF4A", alpha = 0.18) +
    annotate("rect", xmin = 8.5, xmax = 12.5, ymin = -Inf, ymax = Inf,
             fill = "#1B7837", alpha = 0.18) +
    # Concentration labels at top
    annotate("label", x = 2.5, y = Inf, label = "5 ng/ml", vjust = 0,
             fontface = "bold", size = 3.5, fill = "#7FBF7B", color = "white",
             label.size = 0, label.padding = unit(0.25, "lines")) +
    annotate("label", x = 6.5, y = Inf, label = "10 ng/ml", vjust = 0,
             fontface = "bold", size = 3.5, fill = "#4DAF4A", color = "white",
             label.size = 0, label.padding = unit(0.25, "lines")) +
    annotate("label", x = 10.5, y = Inf, label = "15 ng/ml", vjust = 0,
             fontface = "bold", size = 3.5, fill = "#1B7837", color = "white",
             label.size = 0, label.padding = unit(0.25, "lines")) +
    # Data points - colored by p-value
    geom_point(aes(size = Count, fill = neg_log_p),
               shape = 21, color = "white", stroke = 0.6) +
    scale_fill_gradient2(
      low = "#4393C3",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = p_mid,
      name = expression(-log[10](p[adj]))
    ) +
    scale_size_continuous(
      range = c(3, 10),
      name = "Genes",
      breaks = c(5, 10, 20, 30)
    ) +
    scale_x_discrete(labels = x_labels) +
    coord_cartesian(clip = "off") +
    labs(
      x = "Time (min)",
      y = ""
    ) +
    theme_minimal(base_size = 11, base_family = "Helvetica") +
    theme(
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.6),

      axis.text.y = element_text(size = 8.5, color = "grey15", lineheight = 0.85),
      axis.text.x = element_text(size = 9, color = "grey30", face = "bold"),
      axis.title.x = element_text(size = 10, margin = margin(t = 10), face = "bold"),
      axis.ticks = element_line(color = "grey50", linewidth = 0.3),
      axis.ticks.length = unit(0.12, "cm"),

      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),

      plot.margin = margin(t = 30, r = 10, b = 10, l = 5)
    ) +
    guides(
      fill = guide_colorbar(order = 1, barwidth = 0.8, barheight = 4),
      size = guide_legend(order = 2, override.aes = list(fill = "grey50"))
    )

  ggsave("exp1_go_mf_publication.pdf", p_exp1_combined, width = 13, height = 6, device = cairo_pdf)
  cat("Saved: exp1_go_mf_publication.pdf\n")
}

# ============================================================================
# EXPERIMENT 2: PUBLICATION-READY PLOT
# ============================================================================

if (nrow(exp2_go_df) > 0) {

  top_terms2 <- exp2_go_df %>%
    group_by(Description) %>%
    summarize(min_p = min(p.adjust), n_conditions = n(), .groups = "drop") %>%
    arrange(min_p) %>%
    slice_head(n = 15) %>%
    pull(Description)

  plot_df2 <- exp2_go_df %>%
    filter(Description %in% top_terms2) %>%
    mutate(
      neg_log_p = -log10(p.adjust),
      comparison = factor(comparison, levels = c(
        "Activin_vs_Baseline",
        "SB50_60min_vs_Baseline", "SB50_120min_vs_Baseline", "SB50_180min_vs_Baseline",
        "SB50_60min_vs_Activin", "SB50_120min_vs_Activin", "SB50_180min_vs_Activin"
      )),
      Description_short = str_wrap(Description, width = 40),
      comparison_label = case_when(
        comparison == "Activin_vs_Baseline" ~ "Activin",
        comparison == "SB50_60min_vs_Baseline" ~ "60'",
        comparison == "SB50_120min_vs_Baseline" ~ "120'",
        comparison == "SB50_180min_vs_Baseline" ~ "180'",
        comparison == "SB50_60min_vs_Activin" ~ "60'",
        comparison == "SB50_120min_vs_Activin" ~ "120'",
        comparison == "SB50_180min_vs_Activin" ~ "180'"
      ),
      reference_group = case_when(
        comparison == "Activin_vs_Baseline" ~ "Activin effect",
        grepl("vs_Baseline", comparison) ~ "SB50 vs Baseline",
        grepl("vs_Activin", comparison) ~ "SB50 vs Activin"
      )
    )

  term_order2 <- plot_df2 %>%
    group_by(Description_short) %>%
    summarize(mean_neg_log_p = mean(neg_log_p), .groups = "drop") %>%
    arrange(mean_neg_log_p) %>%
    pull(Description_short)

  plot_df2$Description_short <- factor(plot_df2$Description_short, levels = term_order2)

  # Order x-axis: Activin, then SB50 vs Baseline (60, 120, 180), then SB50 vs Activin (60, 120, 180)
  x_order2 <- c("Activin_vs_Baseline",
                "SB50_60min_vs_Baseline", "SB50_120min_vs_Baseline", "SB50_180min_vs_Baseline",
                "SB50_60min_vs_Activin", "SB50_120min_vs_Activin", "SB50_180min_vs_Activin")

  plot_df2$comparison <- factor(plot_df2$comparison, levels = x_order2)

  # Calculate p-value range
  p_range2 <- range(plot_df2$neg_log_p, na.rm = TRUE)
  p_mid2 <- mean(p_range2)

  # Create x labels
  x_labels2 <- c("Activin", "60'", "120'", "180'", "60'", "120'", "180'")

  p_exp2 <- ggplot(plot_df2, aes(x = comparison, y = Description_short)) +
    # Colored background rectangles
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
             fill = "#B2182B", alpha = 0.12) +
    annotate("rect", xmin = 1.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
             fill = "#2166AC", alpha = 0.12) +
    annotate("rect", xmin = 4.5, xmax = 7.5, ymin = -Inf, ymax = Inf,
             fill = "#1B7837", alpha = 0.12) +
    # Group labels at top
    annotate("label", x = 1, y = Inf, label = "Activin\neffect", vjust = 0,
             fontface = "bold", size = 3, fill = "#B2182B", color = "white",
             label.size = 0, label.padding = unit(0.2, "lines")) +
    annotate("label", x = 3, y = Inf, label = "SB50 vs Baseline", vjust = 0,
             fontface = "bold", size = 3, fill = "#2166AC", color = "white",
             label.size = 0, label.padding = unit(0.2, "lines")) +
    annotate("label", x = 6, y = Inf, label = "SB50 vs Activin", vjust = 0,
             fontface = "bold", size = 3, fill = "#1B7837", color = "white",
             label.size = 0, label.padding = unit(0.2, "lines")) +
    # Data points
    geom_point(aes(size = Count, fill = neg_log_p),
               shape = 21, color = "white", stroke = 0.6) +
    scale_fill_gradient2(
      low = "#4393C3", mid = "#F7F7F7", high = "#B2182B",
      midpoint = p_mid2,
      name = expression(-log[10](p[adj]))
    ) +
    scale_size_continuous(range = c(3, 10), name = "Genes", breaks = c(5, 10, 20)) +
    scale_x_discrete(labels = x_labels2) +
    coord_cartesian(clip = "off") +
    labs(x = "", y = "") +
    theme_minimal(base_size = 11, base_family = "Helvetica") +
    theme(
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.6),

      axis.text.y = element_text(size = 8.5, color = "grey15", lineheight = 0.85),
      axis.text.x = element_text(size = 9, color = "grey30", face = "bold"),
      axis.ticks = element_line(color = "grey50", linewidth = 0.3),

      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),

      plot.margin = margin(t = 35, r = 10, b = 10, l = 5)
    ) +
    guides(
      fill = guide_colorbar(order = 1, barwidth = 0.8, barheight = 4),
      size = guide_legend(order = 2, override.aes = list(fill = "grey50"))
    )

  ggsave("exp2_go_mf_publication.pdf", p_exp2, width = 11, height = 6, device = cairo_pdf)
  cat("Saved: exp2_go_mf_publication.pdf\n")
}

# ============================================================================
# COMBINED FIGURE: Both experiments
# ============================================================================

if (nrow(exp1_go_df) > 0 && nrow(exp2_go_df) > 0) {

  p_exp1_tagged <- p_exp1_combined +
    labs(tag = "Exp 1") +
    theme(plot.tag = element_text(size = 12, face = "bold"))

  p_exp2_tagged <- p_exp2 +
    labs(tag = "Exp 2") +
    theme(plot.tag = element_text(size = 12, face = "bold"))

  combined <- p_exp1_tagged / p_exp2_tagged

  ggsave("go_mf_combined_publication.pdf", combined, width = 13, height = 11, device = pdf)
  cat("Saved: go_mf_combined_publication.pdf\n")
}

cat("\n=== Publication plots generated ===\n")
