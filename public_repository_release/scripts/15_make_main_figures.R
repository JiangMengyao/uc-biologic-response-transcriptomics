#!/usr/bin/env Rscript

# R-only, submission-oriented grouped figures for the common prognostic-state
# analysis. Every main figure is exported as SVG, PDF, 600-dpi TIFF and a PNG
# preview, with a panel-tagged source-data CSV.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
tables_dir <- file.path(project_root, "results", "tables", "common_state")
extended_dir <- file.path(project_root, "results", "tables", "extended")
source_dir <- file.path(project_root, "results", "source_data")
figure_dir <- file.path(project_root, "results", "figures", "main")
preview_dir <- file.path(project_root, "results", "figures", "previews")
logs_dir <- file.path(project_root, "logs")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggridges)
  library(ggrepel)
  library(logistf)
  library(splines)
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

theme_paper <- function(base_size = 8.2) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#E7E9EC", linewidth = 0.28),
      axis.title = element_text(colour = "#26313B", face = "bold"),
      axis.text = element_text(colour = "#37424D"),
      strip.text = element_text(face = "bold", colour = "#26313B"),
      strip.background = element_rect(fill = "#F1F3F5", colour = NA),
      plot.title = element_text(face = "bold", size = rel(1.16), colour = "#15202B", margin = margin(b = 4)),
      plot.subtitle = element_text(colour = "#52606D", margin = margin(b = 5)),
      plot.caption = element_text(colour = "#66707F", hjust = 0, size = rel(0.82)),
      legend.title = element_text(face = "bold"), legend.key.height = unit(3.4, "mm")
    )
}

save_figure <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  ggsave(file.path(figure_dir, paste0(stem, ".svg")), plot, width = width_in, height = height_in,
         device = svglite::svglite, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width_in, height = height_in,
         device = grDevices::pdf, bg = "white", useDingbats = FALSE)
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width_in, height = height_in,
         device = ragg::agg_tiff, res = 600, compression = "lzw", bg = "white")
  ggsave(file.path(preview_dir, paste0(stem, ".png")), plot, width = width_in, height = height_in,
         device = ragg::agg_png, res = 180, bg = "white")
}

bind_source <- function(named_list) {
  all_names <- unique(unlist(lapply(named_list, names)))
  rows <- Map(function(x, panel) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x$panel <- panel
    x[, c("panel", setdiff(all_names, "panel")), drop = FALSE]
  }, named_list, names(named_list))
  do.call(rbind, rows)
}

mechanism_label <- c("anti-TNF" = "anti-TNF", "anti-IL12/23" = "UST · anti-IL-12/23",
                     "anti-integrin" = "VDZ · anti-integrin")
short_cohort <- c(
  "GSE16879 IFX" = "IFX · 16879", "GSE73661 IFX" = "IFX · 73661",
  "GSE23597 IFX" = "IFX · 23597", "GSE92415 golimumab" = "GLM · 92415",
  "GSE206285 ustekinumab" = "UST · 206285", "GSE73661 VDZ trial" = "VDZ trial · 73661",
  "GSE73661 VDZ observational" = "VDZ obs · 73661"
)

# ========================= Figure 1: evidence architecture ====================
primary <- read.csv(file.path(tables_dir, "table_primary_cross_mechanism_associations.csv"), check.names = FALSE)
registry <- read.csv(file.path(extended_dir, "table_cohort_registry_summary.csv"), check.names = FALSE)
long_eff <- read.csv(file.path(tables_dir, "table_longitudinal_primary_change_effects.csv"), check.names = FALSE)

cohort_arch <- primary[, c("cohort", "mechanism", "n", "events", "endpoint")]
cohort_arch$cohort_short <- short_cohort[cohort_arch$cohort]
cohort_arch$non_events <- cohort_arch$n - cohort_arch$events
cohort_arch$rct <- cohort_arch$cohort %in% c("GSE23597 IFX", "GSE92415 golimumab", "GSE206285 ustekinumab")
cohort_arch$longitudinal <- cohort_arch$cohort %in% long_eff$cohort
cohort_arch$y <- rev(seq_len(nrow(cohort_arch)))

# A: cohort lanes from baseline biopsy to endpoint.
p1a <- ggplot(cohort_arch) +
  geom_segment(aes(x = 0, xend = 1, y = y, yend = y, colour = mechanism, linewidth = sqrt(n)),
               alpha = 0.82, lineend = "round") +
  geom_point(aes(x = 0, y = y), shape = 21, size = 3.0, stroke = 0.7, fill = "white", colour = "#26313B") +
  geom_point(aes(x = 1, y = y, fill = mechanism), shape = 21, size = 3.6, stroke = 0.7, colour = "white") +
  geom_text(aes(x = -0.04, y = y, label = cohort_short), hjust = 1, size = 2.65, colour = "#26313B") +
  geom_text(aes(x = 1.055, y = y, label = paste0(events, "/", n, " response")), hjust = 0,
            size = 2.45, colour = "#52606D") +
  annotate("text", x = 0, y = max(cohort_arch$y) + 0.8, label = "Baseline mucosal state",
           fontface = "bold", size = 2.9, colour = "#26313B") +
  annotate("text", x = 1, y = max(cohort_arch$y) + 0.8, label = "Clinical endpoint",
           fontface = "bold", size = 2.9, colour = "#26313B") +
  scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
  scale_linewidth(range = c(2.0, 7.2), guide = "none") +
  coord_cartesian(xlim = c(-0.38, 1.33), clip = "off") +
  labs(title = "Cohorts span three biologic mechanisms", colour = "Mechanism") +
  theme_void(base_family = "Helvetica", base_size = 8.2) +
  theme(plot.title = element_text(face = "bold", colour = "#15202B", size = 9.2),
        plot.subtitle = element_text(colour = "#52606D"), legend.position = "none",
        plot.margin = margin(8, 26, 5, 36))

# B: evidence tile matrix.
evidence <- expand.grid(cohort = cohort_arch$cohort, evidence = c("Baseline association", "RCT arm contrast",
                                                                  "Longitudinal tissue", "Outcome interaction"),
                        stringsAsFactors = FALSE)
evidence$available <- with(evidence,
  evidence == "Baseline association" |
  (evidence %in% c("RCT arm contrast", "Outcome interaction") & cohort %in%
     c("GSE23597 IFX", "GSE92415 golimumab", "GSE206285 ustekinumab")) |
  (evidence == "Longitudinal tissue" & cohort %in% long_eff$cohort))
evidence$cohort_short <- short_cohort[evidence$cohort]
evidence$cohort_short <- factor(evidence$cohort_short, levels = rev(cohort_arch$cohort_short))
evidence$evidence <- factor(evidence$evidence, levels = c("Baseline association", "RCT arm contrast",
                                                          "Longitudinal tissue", "Outcome interaction"))
p1b <- ggplot(evidence, aes(evidence, cohort_short, fill = available)) +
  geom_tile(colour = "white", linewidth = 1.1, width = 0.91, height = 0.86) +
  geom_point(data = subset(evidence, available), shape = 21, size = 2.2, fill = "white", colour = "#26313B", stroke = 0.6) +
  scale_fill_manual(values = c(`TRUE` = "#26313B", `FALSE` = "#E9ECEF"), guide = "none") +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(title = "Evidence is deliberately modular", x = NULL, y = NULL) +
  theme_paper() + theme(axis.text.x = element_text(size = 7), panel.grid = element_blank())

# C: response-rate lollipop makes endpoint imbalance explicit.
count_long <- cohort_arch
count_long$response_rate <- count_long$events / count_long$n
count_long$cohort_short <- factor(count_long$cohort_short, levels = rev(count_long$cohort_short))
p1c <- ggplot(count_long, aes(response_rate, cohort_short, colour = mechanism)) +
  geom_segment(aes(x = 0, xend = response_rate, yend = cohort_short), linewidth = 1.6, alpha = .28) +
  geom_point(aes(size = n), alpha = .95) +
  geom_text(aes(label = paste0(events, "/", n)), hjust = -0.25, size = 2.2, colour = "#52606D") +
  scale_colour_manual(values = COL, guide = "none") +
  scale_size(range = c(2.1, 6.3), guide = "none") +
  scale_x_continuous(labels = label_percent(accuracy = 1), limits = c(0, .9), breaks = c(0, .25, .5, .75)) +
  labs(title = "Response rates vary widely",
       x = "Observed response rate", y = NULL) +
  theme_paper(7.5) + theme(plot.title = element_text(size = 9.2), plot.subtitle = element_text(size = 7.3))

# D: fixed claim hierarchy.
claims <- data.frame(
  x = 1:5, y = c(1.0, 1.17, 1.34, 1.51, 1.68),
  label = c("Cross-mechanism", "RCT split", "Paired change", "High-state flow", "Cell localization"),
  role = c("Primary", "Disambiguation", "Primary", "Descriptive", "Orthogonal")
)
p1d <- ggplot(claims, aes(x, y)) +
  geom_segment(aes(x = x, xend = dplyr::lead(x), y = y, yend = dplyr::lead(y)),
               linewidth = 4.5, colour = "#DDE2E6", lineend = "round", na.rm = TRUE) +
  geom_point(aes(fill = role), shape = 21, size = 6.4, colour = "white", stroke = 1) +
  geom_text(aes(label = label), nudge_y = 0.24, size = 2.2, colour = "#26313B") +
  geom_text(aes(label = role), nudge_y = -0.24, size = 1.9, colour = "#66707F") +
  scale_fill_manual(values = c("Primary" = "#D34A5A", "Disambiguation" = "#7257C8",
                               "Descriptive" = "#E5A33C", "Orthogonal" = "#159A8C")) +
  coord_cartesian(xlim = c(0.45, 5.55), ylim = c(0.72, 2.02), clip = "off") +
  labs(title = "Evidence strength bounds the claim", fill = NULL) +
  theme_void(base_family = "Helvetica", base_size = 8.2) +
  theme(plot.title = element_text(face = "bold", colour = "#15202B", size = 9.2),
        legend.position = "none")

fig1_top <- p1a + p1c + plot_layout(widths = c(1.55, 1))
fig1_bottom <- p1b + p1d + plot_layout(widths = c(1.25, 1))
fig1 <- fig1_top / fig1_bottom +
  plot_layout(heights = c(1.05, 1)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(family = "Helvetica", face = "bold", size = 10)))
save_figure(fig1, "Figure1_evidence_architecture", 183, 144)
write.csv(bind_source(list(a_cohort_lanes = cohort_arch, b_outcome_counts = count_long,
                           c_evidence_matrix = evidence, d_claim_hierarchy = claims)),
          file.path(source_dir, "Figure1_source_data.csv"), row.names = FALSE)

# ================= Figure 2: cross-mechanism prognostic state =================
mechanism_meta <- read.csv(file.path(tables_dir, "table_mechanism_meta_analysis.csv"), check.names = FALSE)
module_meta <- read.csv(file.path(tables_dir, "table_module_meta_analysis.csv"), check.names = FALSE)
loo <- read.csv(file.path(tables_dir, "table_leave_one_out.csv"), check.names = FALSE)
patient <- read.csv(file.path(source_dir, "common_state_baseline_patient_data.csv"), check.names = FALSE)

forest_cohort <- primary
forest_cohort$type <- "Cohort"
forest_cohort$label <- short_cohort[forest_cohort$cohort]
forest_meta <- data.frame(
  cohort = mechanism_meta$group, mechanism = ifelse(mechanism_meta$group == "All biologic mechanisms", "Overall", mechanism_meta$group),
  n = NA, events = NA, odds_ratio = mechanism_meta$odds_ratio, ci_lower = mechanism_meta$ci_lower,
  ci_upper = mechanism_meta$ci_upper, p_value = mechanism_meta$p_value, type = "Summary",
  label = c("anti-TNF summary", "UST", "VDZ summary", "Overall summary")
)
forest <- rbind(
  forest_cohort[, intersect(names(forest_cohort), names(forest_meta))],
  forest_meta
)
forest$label <- factor(forest$label, levels = rev(c(
  "IFX · 16879", "IFX · 73661", "IFX · 23597", "GLM · 92415", "anti-TNF summary",
  "UST · 206285", "UST", "VDZ trial · 73661", "VDZ obs · 73661", "VDZ summary", "Overall summary"
)))
forest$plot_colour <- ifelse(forest$mechanism == "Overall", "#15202B", COL[forest$mechanism])
p2a <- ggplot(forest, aes(odds_ratio, label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#7D8792", linewidth = 0.45) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper, colour = plot_colour), orientation = "y", width = 0, linewidth = 0.75) +
  geom_point(aes(colour = plot_colour, size = type, shape = type), fill = "white", stroke = 1) +
  scale_x_log10(limits = c(0.015, 2.0), breaks = c(0.03, 0.1, 0.3, 1, 2)) +
  scale_colour_identity() + scale_size_manual(values = c(Cohort = 2.2, Summary = 3.2), guide = "none") +
  scale_shape_manual(values = c(Cohort = 16, Summary = 18), guide = "none") +
  labs(title = "Cross-mechanism associations align",
       x = "Response odds ratio (log scale)", y = NULL) + theme_paper()

patient$cohort_short <- short_cohort[patient$cohort]
patient$cohort_short <- factor(patient$cohort_short, levels = rev(short_cohort[unique(primary$cohort)]))
p2b <- ggplot(patient, aes(Inflammation_z, cohort_short, fill = response, colour = response)) +
  ggridges::geom_density_ridges(alpha = 0.45, scale = 0.93, rel_min_height = 0.01,
                               position = "identity", linewidth = 0.35) +
  geom_point(aes(y = as.numeric(cohort_short) - 0.36), alpha = 0.22, size = 0.45,
             position = position_jitter(height = 0.06), show.legend = FALSE) +
  scale_fill_manual(values = COL) + scale_colour_manual(values = COL) +
  coord_cartesian(xlim = c(-3.1, 3.2)) +
  labs(title = "Cohort distributions remain distinct",
       x = "Baseline inflammatory state (within-cohort SD)", y = NULL, fill = "Outcome", colour = "Outcome") +
  theme_paper() + theme(legend.position = "bottom")

# Mechanism-level spline curves adjusted for cohort indicator.
curve_rows <- list()
for (mech in unique(patient$mechanism)) {
  x <- patient[patient$mechanism == mech, ]
  x$cohort <- factor(x$cohort)
  nonlinear_model <- nrow(x) >= 80L
  fit <- if (nonlinear_model && nlevels(x$cohort) > 1L) {
    glm(response_binary ~ cohort + splines::ns(Inflammation_z, df = 3), data = x, family = binomial())
  } else if (nonlinear_model) {
    glm(response_binary ~ splines::ns(Inflammation_z, df = 3), data = x, family = binomial())
  } else if (nlevels(x$cohort) > 1L) {
    glm(response_binary ~ cohort + Inflammation_z, data = x, family = binomial())
  } else {
    glm(response_binary ~ Inflammation_z, data = x, family = binomial())
  }
  grid_x <- seq(max(-2.5, quantile(x$Inflammation_z, 0.02)), min(2.5, quantile(x$Inflammation_z, 0.98)), length.out = 120)
  ref_cohort <- levels(x$cohort)[which.max(table(x$cohort))]
  nd <- data.frame(cohort = factor(ref_cohort, levels = levels(x$cohort)), Inflammation_z = grid_x)
  pr <- predict(fit, nd, type = "link", se.fit = TRUE)
  curve_rows[[mech]] <- data.frame(mechanism = mech, Inflammation_z = grid_x,
                                    probability = plogis(pr$fit), ci_lower = plogis(pr$fit - 1.96 * pr$se.fit),
                                    ci_upper = plogis(pr$fit + 1.96 * pr$se.fit), reference_cohort = ref_cohort,
                                    model = ifelse(nonlinear_model, "3-df natural spline", "linear trend"))
}
curves <- do.call(rbind, curve_rows)
p2c <- ggplot(curves, aes(Inflammation_z, probability, colour = mechanism, fill = mechanism)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.13, colour = NA) +
  geom_line(linewidth = 1.05) +
  facet_wrap(~ mechanism, scales = "free_y", labeller = as_labeller(mechanism_label), nrow = 1) +
  scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(title = "The relationship is graded rather than threshold-defined",
       x = "Baseline inflammatory state (SD)", y = "Modelled response probability") +
  theme_paper() + theme(legend.position = "none")

module_names <- c("Inflammation_z" = "Inflammatory response", "TNFaNFKB_z" = "TNF-NF-kB",
                  "IL6JAKSTAT3_z" = "IL6-JAK-STAT3", "OxidativePhos_neg_z" = "OxPhos control",
                  "Random_neg_z" = "Locked random control")
module_meta$module <- factor(module_names[module_meta$score], levels = rev(module_names))
p2d <- ggplot(module_meta, aes(odds_ratio, module, colour = score)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#7D8792") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0, linewidth = 0.75) +
  geom_point(size = 2.8) + scale_x_log10(limits = c(0.18, 2.8), breaks = c(0.25, 0.5, 1, 2)) +
  scale_colour_manual(values = COL, guide = "none") +
  labs(title = "Biological specificity controls the interpretation",
       x = "Response odds ratio (log scale)", y = NULL) + theme_paper()

loo$label <- factor(paste("omit", short_cohort[loo$omitted]), levels = rev(paste("omit", short_cohort[loo$omitted])))
p2e <- ggplot(loo, aes(odds_ratio, label)) +
  geom_vline(xintercept = mechanism_meta$odds_ratio[mechanism_meta$group == "All biologic mechanisms"],
             linewidth = 2.6, colour = "#DDE2E6") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0, colour = "#66707F") +
  geom_point(size = 2.2, colour = "#15202B") +
  scale_x_log10(limits = c(0.22, 0.78), breaks = c(0.25, 0.4, 0.6)) +
  labs(title = "No cohort drives the summary",
       x = "Response odds ratio", y = NULL) + theme_paper()

fig2_top <- p2a + p2b + plot_layout(widths = c(1.08, 1))
fig2_bottom <- p2d + p2e + plot_layout(widths = c(1.05, 1))
fig2 <- fig2_top / p2c / fig2_bottom +
  plot_layout(heights = c(1.45, 1.05, 0.95)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(family = "Helvetica", face = "bold", size = 10)))
save_figure(fig2, "Figure2_cross_mechanism_prognostic_state", 183, 220)
write.csv(bind_source(list(a_forest = forest, b_patient_distributions = patient,
                           c_mechanism_splines = curves, d_module_meta = module_meta,
                           e_leave_one_out = loo)),
          file.path(source_dir, "Figure2_source_data.csv"), row.names = FALSE)

# ================= Figure 3: prognostic versus predictive =====================
rct_arms <- read.csv(file.path(tables_dir, "table_rct_arm_prognostic_associations.csv"), check.names = FALSE)
rct_int <- read.csv(file.path(tables_dir, "table_rct_treatment_score_interactions.csv"), check.names = FALSE)
rct_arms$label <- paste(rct_arms$cohort, rct_arms$arm, sep = " · ")
rct_arms$label <- factor(rct_arms$label, levels = rev(c("GSE92415 · Active", "GSE92415 · Placebo",
                                                        "GSE23597 · Active", "GSE23597 · Placebo",
                                                        "GSE206285 · Active", "GSE206285 · Placebo")))
p3a <- ggplot(rct_arms, aes(odds_ratio, label, colour = arm)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#7D8792") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0, linewidth = 0.8) +
  geom_point(size = 2.8) + scale_x_log10(limits = c(0.09, 3.1), breaks = c(0.1, 0.3, 1, 3)) +
  scale_colour_manual(values = COL) +
  labs(title = "Arm-specific prognostic associations",
       x = "Response odds ratio per 1 SD", y = NULL, colour = NULL) + theme_paper() + theme(legend.position = "bottom")

rct_int$cohort <- factor(rct_int$cohort, levels = rev(c("GSE92415", "GSE23597", "GSE206285")))
p3b <- ggplot(rct_int, aes(odds_ratio, cohort, colour = mechanism)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#7D8792") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0, linewidth = 0.9) +
  geom_point(shape = 18, size = 3.4) + scale_x_log10(limits = c(0.08, 2.3), breaks = c(0.1, 0.3, 1, 2)) +
  scale_colour_manual(values = COL, guide = "none") +
  labs(title = "Treatment interactions remain uncertain",
       x = "Treatment × state interaction OR", y = NULL) + theme_paper()

# Firth absolute-risk curves from reconstructed RCT baseline data.
rct_patient <- patient[patient$cohort %in% c("GSE92415 golimumab", "GSE23597 IFX", "GSE206285 ustekinumab"), ]
# Add placebo participants from cohort-level source objects by rebuilding minimal frames.
rebuild_placebo <- function(accession, visit, active_label, disease_filter = NULL) {
  d <- readRDS(file.path(project_root, "data", "derived", "extended", paste0(tolower(accession), "_all_scores.rds")))
  x <- d[d$visit == visit & d$treatment == "placebo" & d$response %in% c("Yes", "No"), ]
  if (!is.null(disease_filter)) x <- x[grepl(disease_filter, x$disease), ]
  split_x <- split(x, x$subject)
  do.call(rbind, lapply(split_x, function(z) data.frame(subject = z$subject[1], response = z$response[1],
    response_binary = as.integer(z$response[1] == "Yes"), Inflammation_z = mean(z$Inflammation_z), arm = "Placebo")))
}
active_frames <- lapply(c("GSE92415 golimumab", "GSE23597 IFX", "GSE206285 ustekinumab"), function(nm) {
  x <- patient[patient$cohort == nm, c("subject", "response", "response_binary", "Inflammation_z")]
  x$arm <- "Active"; x$cohort <- sub(" .*", "", nm); x
})
placebo_frames <- list(
  transform(rebuild_placebo("GSE92415", "Week 0", "golimumab", "Ulcerative"), cohort = "GSE92415"),
  transform(rebuild_placebo("GSE23597", "W0", "infliximab"), cohort = "GSE23597"),
  transform(rebuild_placebo("GSE206285", "WEEK I-0", "ustekinumab"), cohort = "GSE206285")
)
rct_raw <- rbind(do.call(rbind, active_frames), do.call(rbind, placebo_frames))
rct_raw$arm <- factor(rct_raw$arm, levels = c("Placebo", "Active"))

firth_curve <- function(x) {
  fit <- logistf::logistf(response_binary ~ arm * Inflammation_z, data = x)
  grid <- expand.grid(Inflammation_z = seq(max(-2.4, quantile(x$Inflammation_z, .02)),
                                           min(2.4, quantile(x$Inflammation_z, .98)), length.out = 100),
                      arm = levels(x$arm))
  mm <- model.matrix(~ arm * Inflammation_z, grid)
  beta <- coef(fit)
  vc <- vcov(fit)
  eta <- as.vector(mm %*% beta[colnames(mm)])
  se <- sqrt(pmax(0, diag(mm %*% vc[colnames(mm), colnames(mm)] %*% t(mm))))
  grid$probability <- plogis(eta); grid$ci_lower <- plogis(eta - 1.96 * se); grid$ci_upper <- plogis(eta + 1.96 * se)
  grid$cohort <- x$cohort[1]
  grid
}
rct_curves <- do.call(rbind, lapply(split(rct_raw, rct_raw$cohort), firth_curve))
p3c <- ggplot(rct_curves, aes(Inflammation_z, probability, colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.11, colour = NA) +
  geom_line(linewidth = 0.95) + facet_wrap(~ cohort, nrow = 1) +
  scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(title = "Absolute-risk curves show arm-specific prognosis",
       x = "Baseline inflammatory state (SD)", y = "Modelled response probability", colour = NULL, fill = NULL) +
  theme_paper() + theme(legend.position = "bottom")

# Quadrant map of active versus placebo log-odds slopes.
active <- rct_arms[rct_arms$arm == "Active", c("cohort", "log_odds")]
placebo <- rct_arms[rct_arms$arm == "Placebo", c("cohort", "log_odds")]
quad <- merge(active, placebo, by = "cohort", suffixes = c("_active", "_placebo"))
quad$mechanism <- ifelse(quad$cohort == "GSE206285", "anti-IL12/23", "anti-TNF")
p3d <- ggplot(quad, aes(log_odds_placebo, log_odds_active, label = cohort, colour = mechanism)) +
  geom_hline(yintercept = 0, colour = "#AAB1B8", linetype = 2) + geom_vline(xintercept = 0, colour = "#AAB1B8", linetype = 2) +
  geom_abline(slope = 1, intercept = 0, linewidth = 3, colour = "#EEF0F2") +
  geom_point(size = 3.2) + ggrepel::geom_text_repel(size = 2.5, seed = 20, show.legend = FALSE) +
  annotate("text", x = -0.92, y = -0.06, label = "both arms\nnegative", hjust = 0, vjust = 1,
           size = 2.5, colour = "#66707F") +
  scale_colour_manual(values = COL, guide = "none") +
  labs(title = "Negative active-arm slopes are not automatically predictive",
       x = "Placebo-arm log odds / SD", y = "Active-arm log odds / SD") +
  coord_cartesian(xlim = c(-1, 0.2), ylim = c(-1, 0.2)) + theme_paper()

fig3_top <- p3a + p3b + plot_layout(widths = c(1.05, 1))
fig3 <- fig3_top / p3c / p3d +
  plot_layout(heights = c(1, 1.08, 0.75)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(family = "Helvetica", face = "bold", size = 10)))
save_figure(fig3, "Figure3_prognostic_vs_predictive_RCT", 183, 190)
write.csv(bind_source(list(a_arm_effects = rct_arms, b_interactions = rct_int,
                           c_absolute_risk_curves = rct_curves, d_slope_quadrant = quad)),
          file.path(source_dir, "Figure3_source_data.csv"), row.names = FALSE)

# ================= Figure 4: longitudinal overcoming ==========================
paired <- read.csv(file.path(source_dir, "longitudinal_paired_patient_data.csv"), check.names = FALSE)
change <- read.csv(file.path(tables_dir, "table_longitudinal_primary_change_effects.csv"), check.names = FALSE)
change_meta <- read.csv(file.path(tables_dir, "table_longitudinal_primary_meta.csv"), check.names = FALSE)
module_change <- read.csv(file.path(tables_dir, "table_longitudinal_module_change_effects.csv"), check.names = FALSE)
rct_change <- read.csv(file.path(tables_dir, "table_rct_active_placebo_change.csv"), check.names = FALSE)

paired$cohort_short <- short_cohort[paired$cohort]
paired_long <- rbind(
  data.frame(subject = paired$subject, cohort = paired$cohort, cohort_short = paired$cohort_short,
             mechanism = paired$mechanism, response = paired$response, time = "Baseline", score = paired$Inflammation_z_baseline),
  data.frame(subject = paired$subject, cohort = paired$cohort, cohort_short = paired$cohort_short,
             mechanism = paired$mechanism, response = paired$response, time = "Follow-up", score = paired$Inflammation_z_followup)
)
paired_long$time <- factor(paired_long$time, levels = c("Baseline", "Follow-up"))
p4a <- ggplot(paired_long, aes(time, score, group = interaction(cohort, subject), colour = response)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "#B8BEC5") +
  geom_line(alpha = 0.25, linewidth = 0.42) + geom_point(alpha = 0.46, size = 0.72) +
  stat_summary(aes(group = response), fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(aes(group = response), fun = mean, geom = "point", size = 2.3, shape = 21, fill = "white", stroke = 0.9) +
  facet_wrap(~ cohort_short, nrow = 2) + scale_colour_manual(values = COL) +
  labs(title = "Paired trajectories across six cohorts",
       x = NULL, y = "Inflammatory state (SD)", colour = "Clinical outcome") +
  theme_paper() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 15, hjust = 1))

change_plot <- change
change_plot$label <- factor(short_cohort[change_plot$cohort], levels = rev(short_cohort[change_plot$cohort]))
meta_row <- data.frame(label = factor("Overall", levels = c(levels(change_plot$label), "Overall")),
                       estimate = change_meta$estimate, ci_lower = change_meta$ci_lower, ci_upper = change_meta$ci_upper,
                       mechanism = "Overall")
change_forest <- rbind(change_plot[, c("label", "estimate", "ci_lower", "ci_upper", "mechanism")], meta_row)
change_forest$label <- factor(as.character(change_forest$label), levels = rev(c(short_cohort[change$cohort], "Overall")))
change_forest$plot_colour <- ifelse(change_forest$mechanism == "Overall", "#15202B", COL[change_forest$mechanism])
p4b <- ggplot(change_forest, aes(estimate, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#7D8792") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper, colour = plot_colour), orientation = "y", width = 0, linewidth = 0.8) +
  geom_point(aes(colour = plot_colour, shape = mechanism == "Overall"), size = 2.7) +
  scale_colour_identity() + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  labs(title = "Deeper decrease in responders",
       x = "Difference in change (SD)", y = NULL) + theme_paper()

module_change$module <- factor(module_names[module_change$score], levels = rev(module_names))
module_change$cohort_short <- factor(short_cohort[module_change$cohort], levels = rev(short_cohort[unique(module_change$cohort)]))
p4c <- ggplot(module_change, aes(module, cohort_short, fill = estimate, size = -log10(pmax(p_value, 1e-6)))) +
  geom_point(shape = 21, colour = "white", stroke = 0.45) +
  scale_fill_gradient2(low = "#159A8C", mid = "#F5F6F7", high = "#D34A5A", midpoint = 0,
                       limits = max(abs(module_change$estimate), na.rm = TRUE) * c(-1, 1), oob = squish) +
  scale_size(range = c(2.2, 7.0), name = expression(-log[10](P))) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(title = "Change patterns span multiple modules",
       x = NULL, y = NULL, fill = "Difference") + theme_paper() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 6.6), legend.position = "bottom")

# Individual-flow view for locked high baseline tertile.
high <- paired[paired$baseline_high, ]
high$reversal_label <- ifelse(high$molecular_reversal, "Reversed below\nbaseline median", "Persistently high")
high$response_label <- ifelse(high$response == "Yes", "Clinical response", "Low response")
stage_y <- function(stage, label) {
  if (stage == "Baseline") return(0)
  if (stage == "Molecular") return(ifelse(grepl("Reversed", label), 0.75, -0.75))
  ifelse(label == "Clinical response", 0.75, -0.75)
}
flow_rows <- lapply(seq_len(nrow(high)), function(i) {
  labels <- c("High baseline", high$reversal_label[i], high$response_label[i])
  stages <- c("Baseline", "Molecular", "Clinical")
  y <- c(0, ifelse(high$molecular_reversal[i], 0.75, -0.75), ifelse(high$response[i] == "Yes", 0.75, -0.75))
  # Smoothstep interpolation gives a ribbon-like alluvial trajectory without
  # adding a package dependency.
  pieces <- lapply(1:2, function(k) {
    tt <- seq(0, 1, length.out = 25); smooth <- tt * tt * (3 - 2 * tt)
    data.frame(subject = high$subject[i], cohort = high$cohort[i], response = high$response[i],
               x = (k - 1) + tt + 1, y = y[k] + (y[k + 1] - y[k]) * smooth)
  })
  do.call(rbind, pieces)
})
flow <- do.call(rbind, flow_rows)
p4d <- ggplot(flow, aes(x, y, group = interaction(cohort, subject), colour = response)) +
  geom_path(alpha = 0.34, linewidth = 1.0, lineend = "round") +
  annotate("label", x = 1, y = 0, label = paste0("High baseline\nn=", nrow(high)), size = 2.5,
           fill = "white", colour = "#26313B", linewidth = 0.25) +
  annotate("label", x = 2, y = 0.75, label = paste0("Reversed\nn=", sum(high$molecular_reversal)), size = 2.4,
           fill = "white", colour = "#26313B", linewidth = 0.25) +
  annotate("label", x = 2, y = -0.75, label = paste0("Persisted\nn=", sum(!high$molecular_reversal)), size = 2.4,
           fill = "white", colour = "#26313B", linewidth = 0.25) +
  annotate("label", x = 3, y = 0.75, label = paste0("Response\nn=", sum(high$response == "Yes")), size = 2.4,
           fill = "white", colour = "#26313B", linewidth = 0.25) +
  annotate("label", x = 3, y = -0.75, label = paste0("Low response\nn=", sum(high$response == "No")), size = 2.4,
           fill = "white", colour = "#26313B", linewidth = 0.25) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_x_continuous(breaks = 1:3, labels = c("Baseline", "Molecular", "Clinical")) +
  coord_cartesian(xlim = c(.82, 3.18), ylim = c(-1.12, 1.12), clip = "off") +
  labs(title = "High-state transition audit", x = NULL, y = NULL) +
  theme_void(base_family = "Helvetica", base_size = 8.2) +
  theme(plot.title = element_text(face = "bold", colour = "#15202B"), plot.subtitle = element_text(colour = "#52606D"),
        axis.text.x = element_text(colour = "#37424D", face = "bold"), plot.margin = margin(8, 16, 8, 8))

rct_change$cohort <- factor(rct_change$cohort, levels = rev(rct_change$cohort))
p4e <- ggplot(rct_change, aes(estimate, cohort)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#7D8792") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0, linewidth = 0.8, colour = "#7257C8") +
  geom_point(size = 2.8, colour = "#7257C8") +
  labs(title = "Active-placebo paired change remains uncertain",
       x = "Difference in change (SD)", y = NULL) + theme_paper()

fig4_top <- p4a + p4b + plot_layout(widths = c(1.15, 1))
fig4_middle <- p4c + p4d + plot_layout(widths = c(1.15, 1))
fig4 <- fig4_top / fig4_middle / p4e +
  plot_layout(heights = c(1.25, 1.0, 0.62)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(family = "Helvetica", face = "bold", size = 10)))
save_figure(fig4, "Figure4_longitudinal_overcoming", 183, 224)
write.csv(bind_source(list(a_paired_trajectories = paired_long, b_change_forest = change_forest,
                           c_module_change = module_change, d_high_state_flow = flow,
                           e_rct_change = rct_change)),
          file.path(source_dir, "Figure4_source_data.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(logs_dir, "15_make_main_figures_sessionInfo.txt"))
message("FIGURES 1–4 EXPORTED")
