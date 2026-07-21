#!/usr/bin/env Rscript

# Submission-oriented supplementary figures for cohort flow, baseline
# robustness and longitudinal sensitivity. Rendering is R-only and every panel
# is backed by a panel-tagged source-data CSV.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
tables_dir <- file.path(project_root, "results", "tables", "common_state")
source_dir <- file.path(project_root, "results", "source_data")
figure_dir <- file.path(project_root, "results", "figures", "supplementary")
preview_dir <- file.path(project_root, "results", "figures", "previews")
logs_dir <- file.path(project_root, "logs")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(svglite)
  library(ragg)
})

COL <- c(
  "anti-TNF" = "#D34A5A", "anti-IL12/23" = "#7257C8", "anti-integrin" = "#159A8C",
  "Active" = "#D34A5A", "Placebo" = "#66707F", "Yes" = "#159A8C", "No" = "#C34F5B",
  "Inflammation_z" = "#D34A5A", "TNFaNFKB_z" = "#EE8B46", "IL6JAKSTAT3_z" = "#7257C8",
  "OxidativePhos_neg_z" = "#3E87B7", "Random_neg_z" = "#9AA2AC"
)

theme_paper <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.35, colour = "#26313B"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#26313B"),
      axis.title = element_text(colour = "#26313B", face = "bold"),
      axis.text = element_text(colour = "#37424D"),
      strip.text = element_text(face = "bold", colour = "#26313B"),
      strip.background = element_rect(fill = "#F1F3F5", colour = NA),
      plot.title = element_text(face = "bold", size = rel(1.10), colour = "#15202B",
                                margin = margin(b = 3)),
      plot.subtitle = element_text(colour = "#52606D", size = rel(0.88), margin = margin(b = 4)),
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(3.1, "mm"),
      legend.key.width = unit(4.0, "mm")
    )
}

save_figure <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  ggsave(file.path(figure_dir, paste0(stem, ".svg")), plot,
         width = width_in, height = height_in, device = svglite::svglite, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
         width = width_in, height = height_in, device = grDevices::pdf,
         useDingbats = FALSE, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot,
         width = width_in, height = height_in, device = ragg::agg_tiff,
         res = 600, compression = "lzw", bg = "white")
  ggsave(file.path(preview_dir, paste0(stem, ".png")), plot,
         width = width_in, height = height_in, device = ragg::agg_png,
         res = 180, bg = "white")
}

bind_source <- function(named_list) {
  all_names <- unique(unlist(lapply(named_list, names)))
  rows <- Map(function(x, panel) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x$panel <- panel
    x[, c("panel", setdiff(all_names, "panel")), drop = FALSE]
  }, named_list, names(named_list))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

short_cohort <- c(
  "GSE16879 IFX" = "IFX · 16879", "GSE73661 IFX" = "IFX · 73661",
  "GSE23597 IFX" = "IFX · 23597", "GSE92415 golimumab" = "GLM · 92415",
  "GSE206285 ustekinumab" = "UST · 206285", "GSE73661 VDZ trial" = "VDZ trial · 73661",
  "GSE73661 VDZ observational" = "VDZ obs · 73661"
)
module_names <- c(
  "Inflammation_z" = "Inflammatory response", "TNFaNFKB_z" = "TNF-NF-kB",
  "IL6JAKSTAT3_z" = "IL-6-JAK-STAT3", "OxidativePhos_neg_z" = "OxPhos control",
  "Random_neg_z" = "Locked random control"
)

primary <- read.csv(file.path(tables_dir, "table_primary_cross_mechanism_associations.csv"), check.names = FALSE)
long_eff <- read.csv(file.path(tables_dir, "table_longitudinal_primary_change_effects.csv"), check.names = FALSE)
rct_int <- read.csv(file.path(tables_dir, "table_rct_treatment_score_interactions.csv"), check.names = FALSE)

# ========================== Figure S1: cohort flow ===========================
exclusion <- read.csv(file.path(tables_dir, "table_sample_exclusion_audit.csv"), check.names = FALSE)

cohort_nodes <- primary[, c("cohort", "mechanism", "n", "events")]
cohort_nodes$accession <- sub(" .*", "", cohort_nodes$cohort)
cohort_nodes$id <- cohort_nodes$cohort
cohort_nodes$x <- 1
cohort_nodes$y <- rev(seq_len(nrow(cohort_nodes)))
cohort_nodes$label <- paste0(short_cohort[cohort_nodes$cohort], "\nn=", cohort_nodes$n)

accession_ids <- unique(cohort_nodes$accession)
accession_nodes <- data.frame(
  id = accession_ids, x = 0,
  y = seq(max(cohort_nodes$y) - 0.2, min(cohort_nodes$y) + 0.2, length.out = length(accession_ids)),
  label = accession_ids, stringsAsFactors = FALSE
)

branch_nodes <- data.frame(
  id = c("baseline", "rct", "longitudinal"), x = 2,
  y = c(6.1, 3.8, 1.5),
  label = c(
    sprintf("Baseline prognosis\n7 cohorts · n=%d", sum(primary$n)),
    sprintf("Randomized interaction\n3 comparisons · n=%d", sum(rct_int$n)),
    sprintf("Active longitudinal\n6 cohorts · n=%d", sum(long_eff$n))
  ), stringsAsFactors = FALSE
)

edges_accession <- merge(
  cohort_nodes[, c("accession", "id", "x", "y", "mechanism")],
  accession_nodes[, c("id", "x", "y")], by.x = "accession", by.y = "id",
  suffixes = c("_end", "_start"), sort = FALSE
)
edges_accession$edge_type <- "Accession to active cohort"

edge_rows <- list()
for (i in seq_len(nrow(cohort_nodes))) {
  z <- cohort_nodes[i, ]
  destinations <- "baseline"
  if (z$accession %in% rct_int$cohort) destinations <- c(destinations, "rct")
  if (z$cohort %in% long_eff$cohort) destinations <- c(destinations, "longitudinal")
  for (destination in destinations) {
    b <- branch_nodes[branch_nodes$id == destination, ]
    edge_rows[[length(edge_rows) + 1L]] <- data.frame(
      cohort = z$cohort, mechanism = z$mechanism,
      x_start = z$x, y_start = z$y, x_end = b$x, y_end = b$y,
      edge_type = "Active cohort to evidence branch", branch = destination,
      stringsAsFactors = FALSE
    )
  }
}
edges_branch <- do.call(rbind, edge_rows)

p1a <- ggplot() +
  geom_curve(data = edges_accession,
             aes(x = x_start, y = y_start, xend = x_end, yend = y_end, colour = mechanism),
             curvature = 0.10, linewidth = 0.45, alpha = 0.50,
             arrow = arrow(length = unit(1.1, "mm"), type = "closed")) +
  geom_curve(data = edges_branch,
             aes(x = x_start, y = y_start, xend = x_end, yend = y_end, colour = mechanism),
             curvature = 0.10, linewidth = 0.42, alpha = 0.40,
             arrow = arrow(length = unit(1.0, "mm"), type = "closed")) +
  geom_label(data = accession_nodes, aes(x, y, label = label),
             fill = "#F1F3F5", colour = "#26313B", linewidth = 0,
             size = 2.35, fontface = "bold", label.padding = unit(1.6, "mm")) +
  geom_label(data = cohort_nodes, aes(x, y, label = label, fill = mechanism),
             colour = "white", linewidth = 0, size = 2.05,
             lineheight = 0.92, label.padding = unit(1.2, "mm")) +
  geom_label(data = branch_nodes, aes(x, y, label = label),
             fill = "white", colour = "#26313B", linewidth = 0.35,
             size = 2.10, lineheight = 0.95, label.padding = unit(1.6, "mm")) +
  annotate("text", x = 0, y = 7.85, label = "GEO accession", fontface = "bold", size = 2.6) +
  annotate("text", x = 1, y = 7.85, label = "Eligible active cohort", fontface = "bold", size = 2.6) +
  annotate("text", x = 2, y = 7.85, label = "Evidence branch", fontface = "bold", size = 2.6) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_fill_manual(values = COL, guide = "none") +
  coord_cartesian(xlim = c(-0.25, 2.35), ylim = c(0.25, 8.0), clip = "off") +
  labs(title = "Five accessions support three non-equivalent evidence branches") +
  theme_void(base_family = "Helvetica", base_size = 7.5) +
  theme(plot.title = element_text(face = "bold", size = 8.4, colour = "#15202B"),
        plot.margin = margin(5, 18, 4, 18))

reason_order <- c(
  "Included independent participant", "Duplicate baseline biopsy collapsed within participant",
  "Missing or unusable outcome", "Other treatment arm or cohort branch",
  "Non-baseline/follow-up sample", "Non-UC/control sample"
)
reason_colours <- c(
  "Included independent participant" = "#26313B",
  "Duplicate baseline biopsy collapsed within participant" = "#E5A33C",
  "Missing or unusable outcome" = "#C34F5B",
  "Other treatment arm or cohort branch" = "#AEB6BF",
  "Non-baseline/follow-up sample" = "#D8DDE2",
  "Non-UC/control sample" = "#ECEFF2"
)
exclusion$reason <- factor(exclusion$reason, levels = reason_order)
exclusion$cohort_short <- factor(short_cohort[exclusion$cohort], levels = rev(short_cohort[primary$cohort]))
p1b <- ggplot(exclusion, aes(n_samples, cohort_short, fill = reason)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = reason_colours, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Sample disposition is fully reconciled",
       subtitle = "Mutually exclusive reasons; bars sum to accession expression samples",
       x = "Expression samples", y = NULL, fill = "Disposition") +
  theme_paper(6.8) +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.5),
        legend.key.width = unit(3.0, "mm")) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

branch_summary <- data.frame(
  branch = factor(c("Baseline prognosis", "Randomized interaction", "Active longitudinal"),
                  levels = c("Baseline prognosis", "Randomized interaction", "Active longitudinal")),
  cohorts = c(nrow(primary), nrow(rct_int), nrow(long_eff)),
  participants = c(sum(primary$n), sum(rct_int$n), sum(long_eff$n)),
  events = c(sum(primary$events), sum(rct_int$events), sum(long_eff$responders)),
  evidence_unit = c("treated participant", "randomized participant", "paired treated participant"),
  stringsAsFactors = FALSE
)
branch_summary$branch_display <- factor(
  c("Baseline\nprognosis", "Randomized\ninteraction", "Active\nlongitudinal"),
  levels = c("Baseline\nprognosis", "Randomized\ninteraction", "Active\nlongitudinal")
)
p1c <- ggplot(branch_summary, aes(branch_display, participants, fill = branch)) +
  geom_col(width = 0.62, colour = "white") +
  geom_text(aes(label = paste0("n=", participants)), vjust = -0.45, size = 2.35, fontface = "bold") +
  geom_text(aes(y = participants * 0.50, label = paste0(cohorts, " cohorts\n", events, " events")),
            colour = "white", size = 1.85, lineheight = 0.95) +
  scale_fill_manual(values = c("Baseline prognosis" = "#D34A5A", "Randomized interaction" = "#7257C8",
                               "Active longitudinal" = "#159A8C"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(title = "Analysis units remain branch-specific", x = NULL, y = "Participants") +
  theme_paper(6.8) + theme(axis.text.x = element_text(size = 6.2))

figS1 <- p1a / (p1b | p1c) +
  plot_layout(heights = c(1.52, 1), widths = c(1.32, 1)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(family = "Helvetica", face = "bold", size = 8.5)))
save_figure(figS1, "FigureS1_cohort_inclusion_flow", 183, 154)
write.csv(bind_source(list(
  a_accession_nodes = accession_nodes, a_cohort_nodes = cohort_nodes,
  a_branch_nodes = branch_nodes, a_accession_edges = edges_accession,
  a_branch_edges = edges_branch, b_sample_disposition = exclusion,
  c_evidence_branch_counts = branch_summary
)), file.path(source_dir, "FigureS1_source_data.csv"), row.names = FALSE)

# ================== Figure S2: baseline robustness audit =====================
nonlinear <- read.csv(file.path(tables_dir, "table_nonlinearity_tests.csv"), check.names = FALSE)
loo <- read.csv(file.path(tables_dir, "table_leave_one_out.csv"), check.names = FALSE)
curves <- read.csv(file.path(source_dir, "common_state_nonlinear_predictions.csv"), check.names = FALSE)
module_meta <- read.csv(file.path(tables_dir, "table_module_meta_analysis.csv"), check.names = FALSE)

primary$cohort_short <- factor(short_cohort[primary$cohort], levels = rev(short_cohort[primary$cohort]))
p2a <- ggplot(primary, aes(odds_ratio, cohort_short, colour = mechanism)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#7D8792", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0, linewidth = 0.75) +
  geom_point(size = 2.5) +
  geom_text(aes(label = paste0(events, "/", n)), x = 1.82, hjust = 1, size = 2.0, colour = "#52606D") +
  scale_colour_manual(values = COL, guide = "none") +
  scale_x_log10(limits = c(0.018, 2.0), breaks = c(0.03, 0.1, 0.3, 1, 2)) +
  labs(title = "All cohort estimates retain the adverse direction",
       x = "Response OR per 1-SD higher baseline score", y = NULL) +
  theme_paper()

curves$cohort_short <- factor(short_cohort[curves$cohort], levels = short_cohort[primary$cohort])
p2b <- ggplot(curves, aes(Inflammation_z, probability, colour = mechanism, fill = mechanism)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.75) +
  facet_wrap(~ cohort_short, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = COL, guide = "none") +
  scale_fill_manual(values = COL, guide = "none") +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(title = "Cohort-specific curvature does not define a common threshold",
       x = "Baseline inflammatory state (SD)", y = "Modelled outcome probability") +
  theme_paper(6.2) +
  theme(strip.text = element_text(size = 5.8), axis.text = element_text(size = 5.6))

nonlinear$delta_AIC <- nonlinear$spline_AIC - nonlinear$linear_AIC
nonlinear$minus_log10_p <- -log10(nonlinear$p_nonlinearity)
nonlinear$cohort_short <- short_cohort[nonlinear$cohort]
p2c <- ggplot(nonlinear, aes(delta_AIC, minus_log10_p, colour = mechanism, label = cohort_short)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#AEB6BF") +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.35, linetype = 2, colour = "#AEB6BF") +
  geom_point(size = 2.3) +
  ggrepel::geom_text_repel(size = 1.9, colour = "#37424D", max.overlaps = Inf,
                           box.padding = 0.22, point.padding = 0.18, min.segment.length = 0) +
  scale_colour_manual(values = COL, guide = "none") +
  labs(title = "Non-linearity signals are cohort-specific",
       subtitle = "Negative Delta AIC favours the spline; dashed line denotes P=0.05",
       x = "Spline AIC - linear AIC", y = "-log10(P for non-linearity)") +
  theme_paper(6.6)

loo$label <- factor(paste("omit", short_cohort[loo$omitted]),
                    levels = rev(paste("omit", short_cohort[loo$omitted])))
p2d <- ggplot(loo, aes(odds_ratio, label)) +
  geom_vline(xintercept = 0.4389942, linewidth = 2.2, colour = "#E1E5E9") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0,
                colour = "#66707F", linewidth = 0.65) +
  geom_point(size = 2.1, colour = "#15202B") +
  scale_x_log10(limits = c(0.28, 0.70), breaks = c(0.3, 0.4, 0.5, 0.6)) +
  labs(title = "Leave-one-out summaries remain below one",
       x = "Pooled response OR", y = NULL) +
  theme_paper(6.6)

module_meta$module <- factor(module_names[module_meta$score], levels = rev(module_names))
p2e <- ggplot(module_meta, aes(odds_ratio, module, colour = score)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#7D8792", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0, linewidth = 0.7) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_x_log10(limits = c(0.20, 2.2), breaks = c(0.25, 0.5, 1, 2)) +
  labs(title = "Prespecified module controls",
       x = "Response OR", y = NULL) +
  theme_paper(6.6)

s2_design <- "
AAAA
BBBB
CCDD
EEEE
"
figS2 <- p2a + p2b + p2c + p2d + p2e +
  plot_layout(design = s2_design, heights = c(0.72, 1.22, 0.92, 0.60)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(family = "Helvetica", face = "bold", size = 8.5)))
save_figure(figS2, "FigureS2_baseline_sensitivity", 183, 240)
write.csv(bind_source(list(
  a_complete_cohort_forest = primary, b_cohort_spline_predictions = curves,
  c_nonlinearity_tests = nonlinear, d_leave_one_out = loo,
  e_module_meta_analysis = module_meta
)), file.path(source_dir, "FigureS2_source_data.csv"), row.names = FALSE)

# =============== Figure S3: longitudinal sensitivity audit ==================
mixed <- read.csv(file.path(tables_dir, "table_longitudinal_mixed_sensitivity.csv"), check.names = FALSE)
module_change <- read.csv(file.path(tables_dir, "table_longitudinal_module_change_effects.csv"), check.names = FALSE)

analytic <- long_eff[, c("cohort", "mechanism", "n", "responders", "nonresponders",
                         "estimate", "ci_lower", "ci_upper", "p_value")]
analytic$method <- "Analytic linear model"
bootstrap <- long_eff[, c("cohort", "mechanism", "n", "responders", "nonresponders")]
bootstrap$estimate <- long_eff$bootstrap_estimate
bootstrap$ci_lower <- long_eff$bootstrap_ci_lower
bootstrap$ci_upper <- long_eff$bootstrap_ci_upper
bootstrap$p_value <- long_eff$bootstrap_p
bootstrap$method <- "4,000-resample bootstrap"
long_sensitivity <- rbind(analytic, bootstrap)
long_sensitivity$cohort_short <- factor(short_cohort[long_sensitivity$cohort],
                                        levels = rev(short_cohort[long_eff$cohort]))
p3a <- ggplot(long_sensitivity, aes(estimate, cohort_short, colour = method, shape = method)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#7D8792", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0,
                position = position_dodge(width = 0.42), linewidth = 0.65) +
  geom_point(position = position_dodge(width = 0.42), size = 2.1, fill = "white") +
  scale_colour_manual(values = c("Analytic linear model" = "#66707F",
                                 "4,000-resample bootstrap" = "#159A8C")) +
  scale_shape_manual(values = c("Analytic linear model" = 16, "4,000-resample bootstrap" = 21)) +
  labs(title = "Bootstrap intervals preserve the negative change contrast",
       x = "Responder - non-responder change (SD)", y = NULL, colour = NULL, shape = NULL) +
  theme_paper(6.8) + theme(legend.position = "bottom")

mixed$cohort_short <- factor(short_cohort[mixed$cohort], levels = rev(short_cohort[mixed$cohort]))
p3b <- ggplot(mixed, aes(estimate, cohort_short)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#7D8792", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0,
                colour = "#7257C8", linewidth = 0.7) +
  geom_point(size = 2.3, colour = "#7257C8") +
  geom_text(aes(label = paste0("n=", n)), x = 0.72, hjust = 1, size = 1.9, colour = "#52606D") +
  coord_cartesian(xlim = c(-1.9, 0.8), clip = "off") +
  labs(title = "Participant mixed models remain directionally concordant",
       subtitle = "Time x response coefficient; Wald 95% CI",
       x = "Interaction estimate (SD)", y = NULL) +
  theme_paper(6.8)

module_change$module <- factor(module_names[module_change$score], levels = module_names)
module_change$cohort_short <- factor(short_cohort[module_change$cohort],
                                     levels = rev(short_cohort[long_eff$cohort]))
heat_limit <- max(abs(module_change$estimate), na.rm = TRUE)
p3c <- ggplot(module_change, aes(module, cohort_short, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", estimate)), size = 1.85,
            colour = ifelse(abs(module_change$estimate) > 0.75 * heat_limit, "white", "#26313B")) +
  scale_fill_gradient2(low = "#3E87B7", mid = "white", high = "#D34A5A", midpoint = 0,
                       limits = c(-heat_limit, heat_limit), oob = squish) +
  labs(title = "Complete cohort-by-module change audit",
       subtitle = "Negative values indicate deeper reduction among responders",
       x = NULL, y = NULL, fill = "Change\ncontrast") +
  theme_paper(6.2) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 5.7),
        axis.text.y = element_text(size = 5.8), axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right")

p3d <- ggplot(module_change, aes(estimate, module, colour = mechanism)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#7D8792", linewidth = 0.4) +
  geom_boxplot(aes(group = module), width = 0.52, outlier.shape = NA,
               colour = "#AEB6BF", fill = "#F5F6F7", linewidth = 0.45) +
  geom_point(position = position_jitter(height = 0.10, width = 0), size = 1.75, alpha = 0.90) +
  scale_colour_manual(values = COL, name = "Mechanism") +
  labs(title = "Inflammatory signals separate from control modules",
       x = "Responder - non-responder change (SD)", y = NULL) +
  theme_paper(6.4) + theme(legend.position = "bottom", legend.text = element_text(size = 5.4))

figS3 <- (p3a | p3b) / (p3c | p3d) +
  plot_layout(heights = c(1, 1.1), widths = c(1.25, 1)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(family = "Helvetica", face = "bold", size = 8.5)))
save_figure(figS3, "FigureS3_longitudinal_sensitivity", 183, 184)
write.csv(bind_source(list(
  a_analytic_and_bootstrap = long_sensitivity, b_mixed_models = mixed,
  c_module_heatmap = module_change, d_module_distribution = module_change
)), file.path(source_dir, "FigureS3_source_data.csv"), row.names = FALSE)

# ============= Figure S4: high state and paired RCT audit ====================
transition <- read.csv(file.path(tables_dir, "table_high_state_transition_summary.csv"), check.names = FALSE)
transition_or <- read.csv(file.path(tables_dir, "table_high_state_reversal_response_or.csv"), check.names = FALSE)
rct_change <- read.csv(file.path(tables_dir, "table_rct_active_placebo_change.csv"), check.names = FALSE)

totals <- aggregate(cbind(high_baseline_n, reversal_n) ~ response, transition, sum)
high_total <- sum(totals$high_baseline_n)
response_total <- setNames(totals$high_baseline_n, totals$response)
reversal_total <- setNames(totals$reversal_n, totals$response)
flow_nodes <- data.frame(
  node = c("high", "yes", "no", "yes_rev", "yes_persist", "no_rev", "no_persist"),
  x = c(0, 1, 1, 2, 2, 2, 2),
  y = c(1.5, 2.1, 0.9, 2.55, 1.75, 1.15, 0.35),
  n = c(high_total, response_total["Yes"], response_total["No"],
        reversal_total["Yes"], response_total["Yes"] - reversal_total["Yes"],
        reversal_total["No"], response_total["No"] - reversal_total["No"]),
  label = c(
    sprintf("High baseline tertile\nn=%d", high_total),
    sprintf("Responders\nn=%d", response_total["Yes"]),
    sprintf("Non-responders\nn=%d", response_total["No"]),
    sprintf("Reversal\n%d/%d", reversal_total["Yes"], response_total["Yes"]),
    sprintf("Persistent high\n%d/%d", response_total["Yes"] - reversal_total["Yes"], response_total["Yes"]),
    sprintf("Reversal\n%d/%d", reversal_total["No"], response_total["No"]),
    sprintf("Persistent high\n%d/%d", response_total["No"] - reversal_total["No"], response_total["No"])
  ),
  class = c("Definition", "Yes", "No", "Reversal", "Persistent", "Reversal", "Persistent"),
  stringsAsFactors = FALSE
)
flow_edges <- data.frame(
  from = c("high", "high", "yes", "yes", "no", "no"),
  to = c("yes", "no", "yes_rev", "yes_persist", "no_rev", "no_persist"),
  n = c(response_total["Yes"], response_total["No"], reversal_total["Yes"],
        response_total["Yes"] - reversal_total["Yes"], reversal_total["No"],
        response_total["No"] - reversal_total["No"]),
  stringsAsFactors = FALSE
)
flow_edges <- merge(flow_edges, flow_nodes[, c("node", "x", "y")], by.x = "from", by.y = "node")
flow_edges <- merge(flow_edges, flow_nodes[, c("node", "x", "y")], by.x = "to", by.y = "node",
                    suffixes = c("_start", "_end"))
p4a <- ggplot() +
  geom_curve(data = flow_edges,
             aes(x = x_start, y = y_start, xend = x_end, yend = y_end, linewidth = sqrt(n)),
             curvature = 0.08, colour = "#B7BFC7", alpha = 0.75,
             arrow = arrow(length = unit(1.0, "mm"), type = "closed")) +
  geom_label(data = flow_nodes, aes(x, y, label = label, fill = class),
             colour = ifelse(flow_nodes$class %in% c("Yes", "No", "Definition"), "white", "#26313B"),
             linewidth = 0, size = 2.05, lineheight = 0.95,
             label.padding = unit(1.5, "mm")) +
  scale_fill_manual(values = c("Definition" = "#26313B", "Yes" = "#159A8C", "No" = "#C34F5B",
                               "Reversal" = "#DDEDE8", "Persistent" = "#F3E1E4"), guide = "none") +
  scale_linewidth(range = c(0.4, 2.4), guide = "none") +
  coord_cartesian(xlim = c(-0.25, 2.35), ylim = c(0.05, 2.85), clip = "off") +
  labs(title = "Threshold crossing occurs in both outcome groups",
       subtitle = "High = cohort upper tertile; reversal = follow-up below cohort baseline median") +
  theme_void(base_family = "Helvetica", base_size = 7.2) +
  theme(plot.title = element_text(face = "bold", size = 8.2, colour = "#15202B"),
        plot.subtitle = element_text(size = 6.2, colour = "#52606D"),
        plot.margin = margin(5, 10, 5, 10))

transition$cohort_short <- factor(short_cohort[transition$cohort],
                                  levels = rev(short_cohort[unique(transition$cohort)]))
transition$cohort_index <- as.numeric(transition$cohort_short)
transition$plot_y <- transition$cohort_index + ifelse(transition$response == "Yes", 0.10, -0.10)
p4b <- ggplot(transition, aes(reversal_rate, plot_y, colour = response, size = high_baseline_n)) +
  geom_segment(aes(x = 0, xend = reversal_rate, yend = plot_y), linewidth = 0.7, alpha = 0.25) +
  geom_point(alpha = 0.95) +
  geom_text(aes(label = paste0(reversal_n, "/", high_baseline_n)), hjust = -0.30,
            size = 1.85, colour = "#52606D") +
  scale_colour_manual(values = COL, name = "Outcome") +
  scale_size(range = c(1.8, 4.2), guide = "none") +
  scale_x_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1.18),
                     breaks = c(0, 0.25, 0.50, 0.75, 1)) +
  scale_y_continuous(breaks = seq_along(levels(transition$cohort_short)),
                     labels = levels(transition$cohort_short)) +
  labs(title = "Reversal rates remain variable across cohorts",
       x = "Molecular reversal rate", y = NULL) +
  theme_paper(6.4) + theme(legend.position = "bottom")

transition_or$cohort_short <- factor(short_cohort[transition_or$cohort],
                                     levels = rev(short_cohort[transition_or$cohort]))
estimable_or <- transition_or[is.finite(transition_or$odds_ratio), ]
nonestimable_or <- transition_or[!is.finite(transition_or$odds_ratio), ]
p4c <- ggplot() +
  geom_vline(xintercept = 1, linetype = 2, colour = "#7D8792", linewidth = 0.4) +
  geom_errorbar(data = estimable_or, aes(xmin = ci_lower, xmax = ci_upper, y = cohort_short),
                orientation = "y", width = 0, colour = "#66707F", linewidth = 0.65) +
  geom_point(data = estimable_or, aes(odds_ratio, cohort_short), colour = "#15202B", size = 2.2) +
  geom_text(data = nonestimable_or, aes(x = 0.035, y = cohort_short, label = "not estimable"),
            hjust = 0, size = 1.8, colour = "#8A939D") +
  scale_x_log10(limits = c(0.03, 500), breaks = c(0.03, 0.1, 1, 10, 100, 500)) +
  labs(title = "High-state reversal is not response-specific",
       subtitle = "Exploratory response OR among high-baseline participants",
       x = "Response odds ratio", y = NULL) +
  theme_paper(6.4)

rct_change$cohort_label <- factor(rct_change$cohort,
                                  levels = rev(rct_change$cohort))
p4d <- ggplot(rct_change, aes(estimate, cohort_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#7D8792", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0,
                colour = "#7257C8", linewidth = 0.75) +
  geom_point(shape = 18, size = 3.0, colour = "#7257C8") +
  coord_cartesian(xlim = c(-1.6, 0.85), clip = "off") +
  labs(title = "Paired active-placebo changes remain imprecise",
       subtitle = "Arm-specific means are retained in source data and the legend",
       x = "Active - placebo change (SD)", y = NULL) +
  theme_paper(6.4)

figS4 <- (p4a | p4b) / (p4c | p4d) +
  plot_layout(heights = c(1.08, 1), widths = c(1.08, 1)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(family = "Helvetica", face = "bold", size = 8.5)))
save_figure(figS4, "FigureS4_high_state_rct_audit", 183, 154)
write.csv(bind_source(list(
  a_flow_nodes = flow_nodes, a_flow_edges = flow_edges,
  b_transition_rates = transition, c_reversal_response_or = transition_or,
  d_active_placebo_change = rct_change
)), file.path(source_dir, "FigureS4_source_data.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(logs_dir, "19_make_supplementary_figures_sessionInfo.txt"))
message("SUPPLEMENTARY FIGURES S1–S4 EXPORTED")
