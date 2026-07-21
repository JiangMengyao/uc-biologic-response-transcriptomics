#!/usr/bin/env Rscript

# Main exploratory model, external directional validation, and pre-specified
# robustness analyses.  Claims are limited to what these designs support.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
derived_dir <- file.path(project_root, "data", "derived")
results_dir <- file.path(project_root, "results")
tables_dir <- file.path(results_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

meta_base <- readRDS(file.path(derived_dir, "gse92415_baseline_metadata_scores.rds"))
meta_uc <- readRDS(file.path(derived_dir, "gse16879_baseline_metadata_scores.rds"))

model_one_set <- function(marker) {
  dat <- meta_base
  dat$marker <- dat[[marker]]
  fit <- stats::glm(
    response_binary ~ treatment_binary + marker + treatment_binary:marker,
    data = dat, family = stats::binomial()
  )
  sm <- summary(fit)$coefficients
  ci <- suppressMessages(stats::confint.default(fit))
  model_table <- data.frame(
    term = rownames(sm), estimate = sm[, "Estimate"],
    standard_error = sm[, "Std. Error"], z_value = sm[, "z value"], p_value = sm[, "Pr(>|z|)"],
    odds_ratio = exp(sm[, "Estimate"]), ci_lower = exp(ci[, 1]), ci_upper = exp(ci[, 2]),
    row.names = NULL, check.names = FALSE
  )
  group_corr <- do.call(rbind, lapply(c("golimumab", "placebo"), function(group) {
    x <- dat[dat$treatment == group, , drop = FALSE]
    test <- stats::cor.test(x$marker, x$response_binary, method = "spearman", exact = FALSE)
    data.frame(set = marker, cohort = group, n = nrow(x), events = sum(x$response_binary),
               spearman_rho = unname(test$estimate), p_value = test$p.value)
  }))
  ext <- meta_uc
  ext$marker <- ext[[marker]]
  external_corr <- stats::cor.test(ext$marker, ext$response_binary, method = "spearman", exact = FALSE)
  # Default exact test reproduces the small-cohort result documented in Notion.
  external_wilcox <- stats::wilcox.test(marker ~ response_binary, data = ext)
  external <- data.frame(
    set = marker, cohort = "infliximab_external", n = nrow(ext), events = sum(ext$response_binary),
    spearman_rho = unname(external_corr$estimate), spearman_p_value = external_corr$p.value,
    wilcoxon_p_value = external_wilcox$p.value,
    responder_mean = mean(ext$marker[ext$response_binary == 1]),
    nonresponder_mean = mean(ext$marker[ext$response_binary == 0])
  )
  list(model = model_table, group_corr = group_corr, external = external)
}

all_sets <- c("Inflammation", "TNFaNFKB", "IL6JAKSTAT3", "OxidativePhos_neg", "Random_neg")
fits <- lapply(all_sets, model_one_set)
names(fits) <- all_sets

primary_model <- fits[["Inflammation"]]$model
primary_groups <- fits[["Inflammation"]]$group_corr
primary_external <- fits[["Inflammation"]]$external

robustness <- do.call(rbind, lapply(all_sets, function(set) {
  interaction <- fits[[set]]$model
  interaction <- interaction[interaction$term == "treatment_binary:marker", , drop = FALSE]
  glm_row <- fits[[set]]$group_corr[fits[[set]]$group_corr$cohort == "golimumab", , drop = FALSE]
  pbo_row <- fits[[set]]$group_corr[fits[[set]]$group_corr$cohort == "placebo", , drop = FALSE]
  ext_row <- fits[[set]]$external
  data.frame(
    set = set,
    golimumab_rho = glm_row$spearman_rho, golimumab_p_value = glm_row$p_value,
    placebo_rho = pbo_row$spearman_rho, placebo_p_value = pbo_row$p_value,
    interaction_or = interaction$odds_ratio, interaction_ci_lower = interaction$ci_lower,
    interaction_ci_upper = interaction$ci_upper, interaction_p_value = interaction$p_value,
    infliximab_rho = ext_row$spearman_rho, infliximab_p_value = ext_row$spearman_p_value,
    infliximab_wilcoxon_p_value = ext_row$wilcoxon_p_value
  )
}))

cross_tab <- as.data.frame.matrix(table(meta_base$treatment, meta_base$response))
cross_tab$treatment <- rownames(cross_tab)
rownames(cross_tab) <- NULL
cross_tab$response_rate <- cross_tab$Yes / (cross_tab$Yes + cross_tab$No)
raw_or <- with(cross_tab, (Yes[treatment == "golimumab"] * No[treatment == "placebo"]) /
  (No[treatment == "golimumab"] * Yes[treatment == "placebo"]))

table1_baseline <- do.call(rbind, lapply(c("golimumab", "placebo"), function(group) {
  dat <- meta_base[meta_base$treatment == group, , drop = FALSE]
  data.frame(
    treatment = group, n = nrow(dat), responders = sum(dat$response_binary),
    response_rate = mean(dat$response_binary),
    age_mean = mean(dat$age, na.rm = TRUE), age_sd = stats::sd(dat$age, na.rm = TRUE),
    mayo_mean = mean(dat$mayo_score, na.rm = TRUE), mayo_sd = stats::sd(dat$mayo_score, na.rm = TRUE)
  )
}))

utils::write.csv(table1_baseline, file.path(tables_dir, "table1_baseline_characteristics.csv"), row.names = FALSE)
utils::write.csv(cross_tab, file.path(tables_dir, "table2_response_cross_tab.csv"), row.names = FALSE)
utils::write.csv(primary_model, file.path(tables_dir, "table3_primary_logistic_model.csv"), row.names = FALSE)
utils::write.csv(primary_groups, file.path(tables_dir, "primary_group_directionality.csv"), row.names = FALSE)
utils::write.csv(primary_external, file.path(tables_dir, "external_directional_validation.csv"), row.names = FALSE)
utils::write.csv(robustness, file.path(tables_dir, "table4_robustness.csv"), row.names = FALSE)
saveRDS(list(primary_model = primary_model, primary_groups = primary_groups,
             primary_external = primary_external, robustness = robustness,
             cross_tab = cross_tab, table1_baseline = table1_baseline, raw_or = raw_or),
        file.path(results_dir, "analysis_results.rds"))

message(sprintf("Raw treatment OR: %.3f", raw_or))
print(primary_model)
print(primary_groups)
print(primary_external)
print(robustness)
writeLines(capture.output(sessionInfo()), file.path(project_root, "logs", "03_run_analysis_sessionInfo.txt"))
