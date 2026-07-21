#!/usr/bin/env Rscript

# Publication-oriented presentation and sensitivity analysis.  This script
# re-parameterizes pathway scores per baseline SD for interpretation only; it
# does not change the underlying pre-specified model or its p values.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
derived_dir <- file.path(project_root, "data", "derived")
results_dir <- file.path(project_root, "results")
tables_dir <- file.path(results_dir, "tables")
fig_dir <- file.path(results_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(logistf)
})

main <- readRDS(file.path(derived_dir, "gse92415_baseline_metadata_scores.rds"))
external <- readRDS(file.path(derived_dir, "gse16879_baseline_metadata_scores.rds"))

extract_glm_interaction <- function(data, score) {
  data$marker_z <- as.numeric(scale(data[[score]]))
  fit <- stats::glm(response_binary ~ treatment_binary * marker_z, data = data, family = stats::binomial())
  coefficient <- summary(fit)$coefficients["treatment_binary:marker_z", ]
  ci <- stats::confint.default(fit)["treatment_binary:marker_z", ]
  data.frame(
    model = "Standard logistic", term = score, estimate = coefficient["Estimate"],
    odds_ratio = exp(coefficient["Estimate"]), ci_lower = exp(ci[1]), ci_upper = exp(ci[2]),
    p_value = coefficient["Pr(>|z|)"]
  )
}

extract_firth_interaction <- function(data, adjusted = FALSE) {
  data$marker_z <- as.numeric(scale(data$Inflammation))
  formula <- if (adjusted) {
    response_binary ~ treatment_binary * marker_z + age + mayo_score
  } else {
    response_binary ~ treatment_binary * marker_z
  }
  fit <- logistf::logistf(formula, data = data)
  term <- "treatment_binary:marker_z"
  data.frame(
    model = if (adjusted) "Firth + age + baseline Mayo" else "Firth logistic",
    term = "Inflammation", estimate = unname(fit$coefficients[term]),
    odds_ratio = exp(unname(fit$coefficients[term])), ci_lower = exp(unname(fit$ci.lower[term])),
    ci_upper = exp(unname(fit$ci.upper[term])), p_value = unname(fit$prob[term])
  )
}

pathways <- c("Inflammation", "TNFaNFKB", "IL6JAKSTAT3", "OxidativePhos_neg", "Random_neg")
pathway_forest <- do.call(rbind, lapply(pathways, extract_glm_interaction, data = main))
pathway_forest$class <- ifelse(pathway_forest$term %in% c("Inflammation", "TNFaNFKB", "IL6JAKSTAT3"),
                               "Inflammatory program", "Negative control")
pathway_forest$label <- c("Inflammatory response", "TNFα–NFκB", "IL6–JAK–STAT3", "Oxidative phosphorylation", "Random 200-gene set")

sensitivity <- rbind(
  subset(pathway_forest, term == "Inflammation", select = c(model, term, estimate, odds_ratio, ci_lower, ci_upper, p_value)),
  extract_firth_interaction(main, adjusted = FALSE),
  extract_firth_interaction(main, adjusted = TRUE)
)

# Model-predicted probabilities in the RCT. Scores are displayed over the
# observed 5th–95th percentile to avoid extrapolating beyond the cohort.
main$score_z <- as.numeric(scale(main$Inflammation))
fit_raw <- stats::glm(response_binary ~ treatment_binary * Inflammation, data = main, family = stats::binomial())
score_grid <- seq(stats::quantile(main$Inflammation, 0.05), stats::quantile(main$Inflammation, 0.95), length.out = 120L)
prediction <- do.call(rbind, lapply(c(0L, 1L), function(treatment) {
  new_data <- data.frame(treatment_binary = treatment, Inflammation = score_grid)
  pred <- predict(fit_raw, newdata = new_data, type = "link", se.fit = TRUE)
  data.frame(
    treatment = ifelse(treatment == 1L, "Golimumab", "Placebo"), score = score_grid,
    probability = stats::plogis(pred$fit), lower = stats::plogis(pred$fit - 1.96 * pred$se.fit),
    upper = stats::plogis(pred$fit + 1.96 * pred$se.fit)
  )
}))
main$treatment_label <- ifelse(main$treatment_binary == 1L, "Golimumab", "Placebo")
main$response_label <- ifelse(main$response_binary == 1L, "Response", "No response")

# Bootstrap uncertainty for the single-arm external association. This is an
# association, not a diagnostic-performance or treatment-selection claim.
bootstrap_external <- function(data, B = 2000L, seed = 20260720L) {
  set.seed(seed)
  rho <- numeric(B); superiority <- numeric(B)
  for (i in seq_len(B)) {
    sampled <- data[sample.int(nrow(data), nrow(data), replace = TRUE), , drop = FALSE]
    rho[i] <- suppressWarnings(stats::cor(sampled$Inflammation, sampled$response_binary, method = "spearman"))
    nr <- sampled$Inflammation[sampled$response_binary == 0]
    r <- sampled$Inflammation[sampled$response_binary == 1]
    superiority[i] <- if (length(nr) && length(r)) mean(outer(nr, r, ">")) else NA_real_
  }
  list(rho_ci = stats::quantile(rho, c(0.025, 0.975), na.rm = TRUE),
       superiority_ci = stats::quantile(superiority, c(0.025, 0.975), na.rm = TRUE))
}
external_boot <- bootstrap_external(external)
nr <- external$Inflammation[external$response_binary == 0]
r <- external$Inflammation[external$response_binary == 1]
superiority <- mean(outer(nr, r, ">"))
rho <- stats::cor(external$Inflammation, external$response_binary, method = "spearman")
external_effect <- data.frame(
  cohort = "GSE16879, baseline UC before infliximab", n = nrow(external), responders = sum(external$response_binary),
  spearman_rho = rho, spearman_bootstrap_lower = external_boot$rho_ci[1], spearman_bootstrap_upper = external_boot$rho_ci[2],
  probability_nonresponder_score_higher = superiority,
  superiority_bootstrap_lower = external_boot$superiority_ci[1], superiority_bootstrap_upper = external_boot$superiority_ci[2]
)

# Standardized mean differences are preferred to baseline p values in the RCT.
smd <- function(x, group) {
  group1 <- x[group == 1L]; group0 <- x[group == 0L]
  (mean(group1, na.rm = TRUE) - mean(group0, na.rm = TRUE)) /
    sqrt((stats::var(group1, na.rm = TRUE) + stats::var(group0, na.rm = TRUE)) / 2)
}
baseline_balance <- data.frame(
  variable = c("Age (years)", "Baseline Mayo score"),
  standardized_mean_difference = c(smd(main$age, main$treatment_binary), smd(main$mayo_score, main$treatment_binary))
)

utils::write.csv(pathway_forest, file.path(tables_dir, "tableS3_standardized_pathway_interactions.csv"), row.names = FALSE)
utils::write.csv(sensitivity, file.path(tables_dir, "tableS4_model_sensitivity.csv"), row.names = FALSE)
utils::write.csv(external_effect, file.path(tables_dir, "tableS5_external_bootstrap_effect_size.csv"), row.names = FALSE)
utils::write.csv(baseline_balance, file.path(tables_dir, "tableS6_baseline_standardized_mean_differences.csv"), row.names = FALSE)

theme_publication <- function() {
  theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = "#444444"),
          legend.position = "top", axis.title = element_text(face = "bold"))
}

p2 <- ggplot() +
  geom_point(data = main, aes(x = Inflammation, y = response_binary, color = treatment_label), alpha = 0.6, size = 2) +
  geom_ribbon(data = prediction, aes(x = score, ymin = lower, ymax = upper, fill = treatment), alpha = 0.16, color = NA) +
  geom_line(data = prediction, aes(x = score, y = probability, color = treatment), linewidth = 1.1) +
  scale_color_manual(values = c("Golimumab" = "#0072B2", "Placebo" = "#D55E00")) +
  scale_fill_manual(values = c("Golimumab" = "#0072B2", "Placebo" = "#D55E00")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::label_percent()) +
  labs(title = "Week-6 clinical response across baseline inflammatory-program score",
       subtitle = "GSE92415 randomized cohort: n = 87; curves are logistic-model predictions with 95% Wald CI",
       x = "Baseline Hallmark inflammatory-response score", y = "Observed binary outcome / predicted response probability",
       color = NULL, fill = NULL,
       caption = "Exploratory interaction OR = 0.34 (95% CI 0.016–7.31; P = 0.493). Wide intervals preclude a predictive claim.") +
  theme_publication()
ggsave(file.path(fig_dir, "figure2_modelled_response_probability.png"), p2, width = 8.5, height = 5.7, dpi = 350)

p3 <- ggplot(pathway_forest, aes(x = odds_ratio, y = reorder(label, odds_ratio), color = class, shape = class)) +
  geom_vline(xintercept = 1, linetype = 2, color = "#444444") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0.17, linewidth = 0.65) +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_color_manual(values = c("Inflammatory program" = "#0072B2", "Negative control" = "#666666")) +
  scale_shape_manual(values = c("Inflammatory program" = 16, "Negative control" = 1)) +
  labs(title = "Treatment-by-program interaction across pre-specified gene sets",
       subtitle = "GSE92415; odds ratios are per 1 baseline-SD program score, with 95% Wald CI",
       x = "Interaction odds ratio (log scale)", y = NULL, color = NULL, shape = NULL,
       caption = "All confidence intervals include the null. The figure tests directional consistency, not multiplicity-adjusted discovery.") +
  theme_publication()
ggsave(file.path(fig_dir, "figure3_pathway_interaction_forest.png"), p3, width = 8.5, height = 5.2, dpi = 350)

pS1 <- ggplot(sensitivity, aes(x = odds_ratio, y = reorder(model, odds_ratio))) +
  geom_vline(xintercept = 1, linetype = 2, color = "#444444") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0.16, linewidth = 0.7, color = "#0072B2") +
  geom_point(size = 3, color = "#0072B2") + scale_x_log10() +
  labs(title = "Sensitivity of the inflammatory-program interaction estimate",
       subtitle = "Per 1 baseline-SD score; Firth models use profile penalized-likelihood 95% CI",
       x = "Interaction odds ratio (log scale)", y = NULL,
       caption = "Direction persists after small-sample penalization and covariate adjustment, but each interval remains compatible with no interaction.") +
  theme_publication()
ggsave(file.path(fig_dir, "figureS1_interaction_model_sensitivity.png"), pS1, width = 8.5, height = 4.4, dpi = 350)

external$response_label <- factor(ifelse(external$response_binary == 1L, "Mucosal healing", "No mucosal healing"),
                                  levels = c("Mucosal healing", "No mucosal healing"))
external_summary <- do.call(rbind, lapply(levels(external$response_label), function(group) {
  x <- external$Inflammation[external$response_label == group]
  data.frame(response_label = group, median = stats::median(x), q25 = stats::quantile(x, 0.25), q75 = stats::quantile(x, 0.75))
}))
external_annotation <- sprintf("Spearman ρ = %.2f (bootstrap 95%% CI %.2f to %.2f)\nP(NR score > R score) = %.2f (%.2f to %.2f)",
                               rho, external_boot$rho_ci[1], external_boot$rho_ci[2], superiority,
                               external_boot$superiority_ci[1], external_boot$superiority_ci[2])
p4 <- ggplot(external, aes(x = response_label, y = Inflammation, fill = response_label)) +
  geom_linerange(data = external_summary, aes(x = response_label, ymin = q25, ymax = q75), inherit.aes = FALSE, linewidth = 1.2, color = "#333333") +
  geom_point(data = external_summary, aes(x = response_label, y = median), inherit.aes = FALSE, shape = 95, size = 15, color = "#333333") +
  geom_jitter(aes(color = response_label, shape = response_label), width = 0.07, height = 0, size = 2.8, alpha = 0.9) +
  annotate("label", x = 1.5, y = max(external$Inflammation) + 0.18, label = external_annotation, size = 3.4, linewidth = 0.25) +
  scale_fill_manual(values = c("Mucosal healing" = "#009E73", "No mucosal healing" = "#CC79A7"), guide = "none") +
  scale_color_manual(values = c("Mucosal healing" = "#009E73", "No mucosal healing" = "#CC79A7"), guide = "none") +
  scale_shape_manual(values = c("Mucosal healing" = 16, "No mucosal healing" = 17), guide = "none") +
  labs(title = "External directional association of baseline inflammation with mucosal healing",
       subtitle = "GSE16879, baseline UC before infliximab: n = 24 (8 healing; 16 no healing); all individual patients shown",
       x = NULL, y = "Baseline Hallmark inflammatory-response score", fill = NULL,
       caption = "Single-arm cohort: this supports an association only and cannot replicate the randomized treatment interaction.") +
  coord_cartesian(ylim = c(min(external$Inflammation) - 0.08, max(external$Inflammation) + 0.33)) +
  theme_publication()
ggsave(file.path(fig_dir, "figure4_external_effect_size.png"), p4, width = 8.2, height = 5.7, dpi = 350)

writeLines(capture.output(sessionInfo()), file.path(project_root, "logs", "06_advanced_statistics_sessionInfo.txt"))
message("Advanced sensitivity tables and publication-oriented figures written to results/.")
