#!/usr/bin/env Rscript

# R-only journal figure construction. Figure 5 is localization-only. Baseline
# remission models, aggregation sensitivities and Match=Yes longitudinal
# analyses are exported intact as Supplementary Figure S5.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
tables_dir <- file.path(project_root, "results", "tables", "common_state")
source_dir <- file.path(project_root, "results", "source_data")
main_dir <- file.path(project_root, "results", "figures", "main")
supp_dir <- file.path(project_root, "results", "figures", "supplementary")
preview_dir <- file.path(project_root, "results", "figures", "previews")
logs_dir <- file.path(project_root, "logs")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
})

compartment_levels <- c("Colonic epithelium", "Myeloid", "Fibroblast/pericyte")
compartment_labels <- c(
  "Colonic epithelium" = "Epithelium",
  "Myeloid" = "Myeloid",
  "Fibroblast/pericyte" = "Stroma"
)
compartment_colours <- c(
  "Colonic epithelium" = "#376B9F",
  "Myeloid" = "#A84A3F",
  "Fibroblast/pericyte" = "#3F7D64"
)
outcome_colours <- c("Non-remission" = "#777E86", "Remission" = "#2D6D9E")
ink <- "#20262E"
muted <- "#66717D"
grid <- "#E2E5E8"
low_colour <- "#315A78"
mid_colour <- "#F5F3EE"
high_colour <- "#9C453F"

theme_paper <- function(base_size = 7.4) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(colour = ink, linewidth = 0.32),
      axis.ticks = element_line(colour = ink, linewidth = 0.28),
      axis.title = element_text(colour = ink, size = rel(.96)),
      axis.text = element_text(colour = ink),
      plot.title = element_text(face = "bold", size = rel(1.04), colour = ink,
                                hjust = 0, margin = margin(b = 3)),
      strip.background = element_rect(fill = "#F1F2F3", colour = NA),
      strip.text = element_text(face = "bold", colour = ink, hjust = 0,
                                margin = margin(1.5, 2.5, 1.5, 2.5)),
      legend.title = element_text(face = "bold", size = rel(.92)),
      legend.text = element_text(size = rel(.9)),
      legend.key.height = unit(3.2, "mm"),
      plot.margin = margin(3.5, 4.5, 3.5, 4)
    )
}

save_bundle <- function(plot, stem, output_dir, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(file.path(output_dir, paste0(stem, ".svg")),
                   width = width_in, height = height_in, bg = "white")
  print(plot)
  grDevices::dev.off()
  # This macOS R build has no XQuartz Cairo libraries. The base PDF device is
  # still an R-only vector export and keeps text editable.
  grDevices::pdf(file.path(output_dir, paste0(stem, ".pdf")),
                 width = width_in, height = height_in,
                 family = "Helvetica", useDingbats = FALSE, bg = "white")
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(file.path(output_dir, paste0(stem, ".tiff")),
                 width = width_in, height = height_in, units = "in",
                 res = 600, compression = "lzw", background = "white")
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(file.path(preview_dir, paste0(stem, ".png")),
                width = width_in, height = height_in, units = "in",
                res = 180, background = "white")
  print(plot)
  grDevices::dev.off()
}

bind_source <- function(named_list) {
  all_names <- unique(unlist(lapply(named_list, names)))
  rows <- Map(function(x, panel) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    for (nm in setdiff(all_names, names(x))) x[[nm]] <- NA
    x$panel <- panel
    x[, c("panel", setdiff(all_names, "panel")), drop = FALSE]
  }, named_list, names(named_list))
  do.call(rbind, rows)
}

clean_state <- function(x) {
  x <- sub("^Non ileal ", "", x)
  x <- gsub("pos\\b", "+", x, perl = TRUE)
  x <- gsub("([[:alnum:]])hi\\b", "\\1-high", x, perl = TRUE)
  x <- gsub("([[:alnum:]])lo\\b", "\\1-low", x, perl = TRUE)
  x <- gsub("_", " ", x, fixed = TRUE)
  gsub("  +", " ", x)
}

localization <- read.csv(file.path(tables_dir, "table_singlecell_compartment_localization.csv"),
                         check.names = FALSE)
donor_map <- read.csv(file.path(source_dir, "singlecell_donor_compartment_baseline.csv"),
                      check.names = FALSE)
dot <- read.csv(file.path(source_dir, "singlecell_compartment_key_gene_dotplot.csv"),
                check.names = FALSE)
outcome <- read.csv(file.path(tables_dir, "table_singlecell_compartment_outcome_effects.csv"),
                    check.names = FALSE)
sensitivity <- read.csv(file.path(tables_dir, "table_singlecell_compartment_outcome_sensitivity.csv"),
                        check.names = FALSE)
longitudinal <- read.csv(file.path(tables_dir, "table_singlecell_compartment_longitudinal_effects.csv"),
                         check.names = FALSE)
paired <- read.csv(file.path(source_dir, "singlecell_compartment_paired_donors.csv"),
                   check.names = FALSE)

for (object_name in c("localization", "donor_map", "dot", "outcome",
                      "sensitivity", "longitudinal", "paired")) {
  x <- get(object_name)
  x$compartment <- factor(x$compartment, levels = compartment_levels)
  assign(object_name, x)
}

# Figure 5a: the six highest localized states per compartment. Selection is by
# localization score only and never by clinical direction.
eligible_loc <- localization[localization$donors >= 3L &
                               is.finite(localization$mean_module_score), ]
top_localization <- do.call(rbind, lapply(
  split(eligible_loc, eligible_loc$compartment), function(z) {
    head(z[order(z$mean_module_score, decreasing = TRUE), ], 6L)
  }
))
top_localization$state_label <- clean_state(top_localization$state)
top_localization$state_key <- paste(top_localization$compartment,
                                    top_localization$state_label, sep = "___")
state_levels <- unlist(lapply(split(top_localization, top_localization$compartment),
                              function(z) z$state_key[order(z$mean_module_score)]),
                       use.names = FALSE)
top_localization$state_key <- factor(top_localization$state_key,
                                     levels = unique(state_levels))

p5a <- ggplot(top_localization,
              aes(mean_module_score, state_key, colour = compartment)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .32, colour = muted) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y",
                width = 0, linewidth = .58, na.rm = TRUE) +
  geom_point(aes(size = donors), alpha = .98) +
  facet_wrap(~ compartment, ncol = 1, scales = "free_y",
             labeller = as_labeller(compartment_labels)) +
  scale_colour_manual(values = compartment_colours, guide = "none") +
  scale_size_continuous(range = c(1.4, 3.3), breaks = c(4, 8, 12, 16),
                        name = "Donors") +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  labs(title = "Cell-state localization",
       x = "Inflammatory-state score (mean, 95% CI)", y = NULL) +
  theme_paper(7.15) +
  theme(panel.grid.major.x = element_line(colour = grid, linewidth = .24),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal")

# Figure 5b: sample-balanced donor map. Clinical outcome is intentionally not
# encoded in the main figure.
donor_order <- aggregate(module_score ~ donor, donor_map, mean)
donor_order <- donor_order[order(donor_order$module_score), ]
donor_map$donor <- factor(donor_map$donor, levels = donor_order$donor)
donor_map$compartment_short <- factor(compartment_labels[as.character(donor_map$compartment)],
                                      levels = compartment_labels[compartment_levels])
p5b <- ggplot(donor_map, aes(compartment_short, donor, fill = module_score)) +
  geom_tile(colour = "white", linewidth = .55) +
  scale_fill_gradient2(low = low_colour, mid = mid_colour, high = high_colour,
                       midpoint = 0, name = "State score") +
  labs(title = "Donor map", x = NULL, y = "UC donor") +
  theme_paper(7.0) +
  theme(panel.grid = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

# Figure 5c: a fixed key-gene list displayed in the two most localized states
# per compartment. Both means and fractions are sample-balanced.
top_two <- do.call(rbind, lapply(split(top_localization, top_localization$compartment),
                                 function(z) head(z[order(z$mean_module_score,
                                                          decreasing = TRUE), ], 2L)))
dot_top <- merge(dot, top_two[, c("compartment", "state")],
                 by = c("compartment", "state"))
dot_top$compartment <- factor(dot_top$compartment, levels = compartment_levels)
dot_top$row_label <- paste(compartment_labels[as.character(dot_top$compartment)],
                           clean_state(dot_top$state), sep = " | ")
row_order <- unique(unlist(lapply(split(top_two, top_two$compartment), function(z) {
  paste(compartment_labels[as.character(z$compartment)], clean_state(z$state), sep = " | ")
}), use.names = FALSE))
dot_top$row_label <- factor(dot_top$row_label, levels = rev(row_order))
gene_order <- c("IL1B", "TNF", "IL6", "CXCL8", "OSM", "CCL2", "CXCL10",
                "PLAUR", "ICAM1", "NFKB1", "STAT3", "SOCS3", "FOS", "JUN", "JAK2")
dot_top$gene <- factor(dot_top$gene, levels = gene_order[gene_order %in% dot_top$gene])
p5c <- ggplot(dot_top,
              aes(gene, row_label, size = pct_expressing, fill = mean_expression_z)) +
  geom_point(shape = 21, colour = "white", stroke = .32) +
  scale_size_continuous(range = c(.8, 5.0), labels = label_percent(accuracy = 1),
                        name = "Samples: mean\npositive fraction") +
  scale_fill_gradient2(low = low_colour, mid = mid_colour, high = high_colour,
                       midpoint = 0, limits = c(-2.5, 2.5), oob = squish,
                       name = "Mean expression z") +
  labs(title = "Key-gene expression", x = NULL, y = NULL) +
  theme_paper(7.2) +
  theme(panel.grid = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "right")

figure5 <- (p5a | p5b) / p5c +
  plot_layout(widths = c(1.36, .64), heights = c(1.47, .78), guides = "keep") +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(family = "Helvetica", face = "bold",
                                           size = 9, colour = ink))
  )

save_bundle(figure5, "Figure5_cross_compartment_localization",
            main_dir, width_mm = 183, height_mm = 158)
write.csv(bind_source(list(
  a_state_localization_displayed = top_localization,
  a_state_localization_complete = localization,
  b_donor_compartment = donor_map,
  c_key_genes = dot_top
)), file.path(source_dir, "Figure5_source_data.csv"), row.names = FALSE)

# Supplementary Figure S5a: corrected sample-first baseline audit.
outcome$label <- factor(compartment_labels[as.character(outcome$compartment)],
                        levels = rev(compartment_labels[compartment_levels]))
pS5a <- ggplot(outcome, aes(odds_ratio, label, colour = compartment)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = .32, colour = muted) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y",
                width = 0, linewidth = .68, na.rm = TRUE) +
  geom_point(size = 2.35, na.rm = TRUE) +
  scale_x_log10(breaks = c(.25, .5, 1, 2, 4)) +
  scale_colour_manual(values = compartment_colours, guide = "none") +
  labs(title = "Baseline remission audit",
       x = "Remission OR per 1-SD score", y = NULL) +
  theme_paper(7.3) +
  theme(panel.grid.major.x = element_line(colour = grid, linewidth = .24),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

# Supplementary Figure S5b: every prespecified aggregation/anatomical
# sensitivity, including the unmodified legacy result.
method_levels <- c(
  "Sample-first | Inflamed only",
  "Sample-first | All baseline samples",
  "Sample-first | Non-inflamed only",
  "Sample-first | Inflamed rectum",
  "Sample-first | Inflamed distal colon",
  "Legacy state-cell-weighted"
)
method_labels <- c(
  "Sample-first | Inflamed only" = "Inflamed (primary)",
  "Sample-first | All baseline samples" = "All baseline",
  "Sample-first | Non-inflamed only" = "Non-inflamed",
  "Sample-first | Inflamed rectum" = "Inflamed rectum",
  "Sample-first | Inflamed distal colon" = "Inflamed distal colon",
  "Legacy state-cell-weighted" = "Legacy cell-weighted"
)
sensitivity$method <- factor(method_labels[sensitivity$analysis_label],
                             levels = rev(method_labels[method_levels]))
sensitivity$compartment_short <- factor(compartment_labels[as.character(sensitivity$compartment)],
                                        levels = compartment_labels[compartment_levels])
sensitivity$legacy <- sensitivity$analysis_label == "Legacy state-cell-weighted"
legacy_epi <- sensitivity[sensitivity$legacy &
                            sensitivity$compartment == "Colonic epithelium", ]
pS5b <- ggplot(sensitivity, aes(odds_ratio, method, colour = compartment)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = .3, colour = muted) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y",
                width = 0, linewidth = .55, na.rm = TRUE) +
  geom_point(aes(shape = legacy), size = 2.05, na.rm = TRUE) +
  geom_point(data = legacy_epi, shape = 21, size = 3.0, stroke = .8,
             fill = "#F5C6BF", colour = "#9C2F26") +
  geom_text(data = legacy_epi, aes(x = .019, y = method), inherit.aes = FALSE,
            label = "bulk-discordant\nP=0.041; FDR=0.123",
            hjust = 0, vjust = -.48, colour = "#8B2E27", size = 2.15,
            family = "Helvetica") +
  facet_wrap(~ compartment_short, nrow = 1) +
  scale_x_log10(limits = c(.015, 20),
                breaks = c(.02, .1, .5, 1, 5, 20),
                labels = label_number(accuracy = .01)) +
  scale_colour_manual(values = compartment_colours, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  labs(title = "Baseline audit across aggregation choices",
       x = "Remission OR per 1-SD score", y = NULL) +
  theme_paper(6.95) +
  theme(panel.grid.major.x = element_line(colour = grid, linewidth = .22),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text = element_text(hjust = .5))

# Supplementary Figure S5c: all donor-level Match=Yes changes retained.
paired$outcome <- factor(ifelse(paired$remission == "Remission",
                                "Remission", "Non-remission"),
                         levels = c("Non-remission", "Remission"))
paired$compartment_short <- factor(compartment_labels[as.character(paired$compartment)],
                                   levels = compartment_labels[compartment_levels])
pS5c <- ggplot(paired, aes(outcome, delta, colour = outcome)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = .3, colour = muted) +
  geom_boxplot(width = .48, outlier.shape = NA, colour = ink, fill = "white",
               linewidth = .42) +
  geom_point(position = position_jitter(width = .09, height = 0, seed = 20260720),
             size = 1.35, alpha = .9) +
  facet_wrap(~ compartment_short, nrow = 1) +
  scale_colour_manual(values = outcome_colours, guide = "none") +
  labs(title = "Match=Yes donor-level molecular change",
       x = NULL, y = "Post - Pre score") +
  theme_paper(7.0) +
  theme(panel.grid.major.y = element_line(colour = grid, linewidth = .22),
        axis.text.x = element_text(angle = 20, hjust = 1),
        strip.text = element_text(hjust = .5))

# Supplementary Figure S5d: donor-level group contrast for the exact pairs.
longitudinal$label <- factor(compartment_labels[as.character(longitudinal$compartment)],
                             levels = rev(compartment_labels[compartment_levels]))
pS5d <- ggplot(longitudinal, aes(estimate, label, colour = compartment)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .32, colour = muted) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y",
                width = 0, linewidth = .68, na.rm = TRUE) +
  geom_point(size = 2.35, na.rm = TRUE) +
  scale_colour_manual(values = compartment_colours, guide = "none") +
  labs(title = "Longitudinal effect estimate",
       x = "Remission - non-remission change", y = NULL) +
  theme_paper(7.3) +
  theme(panel.grid.major.x = element_line(colour = grid, linewidth = .24),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

figureS5 <- (pS5a | pS5d) / pS5b / pS5c +
  plot_layout(widths = c(1, 1), heights = c(.70, 1.15, .95), guides = "keep") +
  plot_annotation(
    tag_levels = "a",
    caption = paste0(
      "Clinical audit only. The legacy epithelial estimate is opposite to the bulk mainline ",
      "(nominal P=0.041; FDR=0.123) and is not reproduced by sample-first analyses."
    ),
    theme = theme(
      plot.tag = element_text(family = "Helvetica", face = "bold", size = 9, colour = ink),
      plot.caption = element_text(family = "Helvetica", size = 6.5, colour = muted,
                                  hjust = 0, margin = margin(t = 4))
    )
  )

save_bundle(figureS5, "FigureS5_singlecell_clinical_audit",
            supp_dir, width_mm = 183, height_mm = 181)
write.csv(bind_source(list(
  a_baseline_primary = outcome,
  b_match_yes_longitudinal_effects = longitudinal,
  c_baseline_sensitivity_complete = sensitivity,
  d_match_yes_donor_changes = paired
)), file.path(source_dir, "FigureS5_source_data.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(logs_dir, "16_make_singlecell_figure_sessionInfo.txt"))
message("FIGURE EXPORT COMPLETE: localization-only Figure 5 and clinical-audit Figure S5")
