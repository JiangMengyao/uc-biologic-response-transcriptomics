#!/usr/bin/env Rscript

# Build submission-facing, machine-readable Table S1-S7 CSV files from the
# frozen result objects. These tables reorganize existing estimates only; no
# model is refitted and no result is selected by significance.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
common_dir <- file.path(project_root, "results", "tables", "common_state")
output_dir <- file.path(project_root, "results", "supplementary_tables")
logs_dir <- file.path(project_root, "logs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_common <- function(name) {
  utils::read.csv(file.path(common_dir, name), stringsAsFactors = FALSE, check.names = FALSE)
}

bind_fill <- function(rows) {
  rows <- rows[vapply(rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))]
  all_names <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

accession_from <- function(x) sub(" .*", "", x)
treatment_from <- function(x) {
  ifelse(grepl("ustekinumab", x, ignore.case = TRUE), "ustekinumab",
         ifelse(grepl("golimumab", x, ignore.case = TRUE), "golimumab",
                ifelse(grepl("VDZ", x, ignore.case = TRUE), "vedolizumab",
                       ifelse(grepl("IFX", x, ignore.case = TRUE), "infliximab", NA_character_))))
}

primary <- read_common("table_primary_cross_mechanism_associations.csv")
rct_arms <- read_common("table_rct_arm_prognostic_associations.csv")
rct_interactions <- read_common("table_rct_treatment_score_interactions.csv")
longitudinal <- read_common("table_longitudinal_primary_change_effects.csv")
rct_change <- read_common("table_rct_active_placebo_change.csv")
exclusion <- read_common("table_sample_exclusion_audit.csv")

# Table S1 -------------------------------------------------------------------
active_s1 <- data.frame(
  accession = accession_from(primary$cohort),
  analysis_cohort = primary$cohort,
  treatment = treatment_from(primary$cohort),
  mechanism = primary$mechanism,
  arm = "Active",
  endpoint = primary$endpoint,
  baseline_n = primary$n,
  outcome_events = primary$events,
  outcome_non_events = primary$non_events,
  primary_baseline_cohort = TRUE,
  randomized_interaction_available = accession_from(primary$cohort) %in% rct_interactions$cohort,
  active_longitudinal_available = primary$cohort %in% longitudinal$cohort,
  active_longitudinal_n = longitudinal$n[match(primary$cohort, longitudinal$cohort)],
  paired_rct_change_n = rct_change$n[match(accession_from(primary$cohort), rct_change$cohort)],
  stringsAsFactors = FALSE
)

placebo <- rct_arms[rct_arms$arm == "Placebo", ]
placebo_s1 <- data.frame(
  accession = placebo$cohort,
  analysis_cohort = paste(placebo$cohort, "placebo"),
  treatment = "placebo",
  mechanism = placebo$mechanism,
  arm = "Placebo",
  endpoint = placebo$endpoint,
  baseline_n = placebo$n,
  outcome_events = placebo$events,
  outcome_non_events = placebo$non_events,
  primary_baseline_cohort = FALSE,
  randomized_interaction_available = TRUE,
  active_longitudinal_available = FALSE,
  active_longitudinal_n = NA_integer_,
  paired_rct_change_n = rct_change$n[match(placebo$cohort, rct_change$cohort)],
  stringsAsFactors = FALSE
)

table_s1 <- rbind(active_s1, placebo_s1)
site_map <- c(
  "GSE23597" = "Protocol-fixed distal/sigmoid, approximately 15-20 cm",
  "GSE92415" = "Protocol-fixed distal/sigmoid, approximately 15-20 cm",
  "GSE206285" = "Protocol-fixed distal/sigmoid, approximately 15-20 cm",
  "GSE16879" = "Colonic mucosa; exact site not released",
  "GSE73661" = "Colonic mucosa; exact site not released"
)
table_s1$biopsy_site_documentation <- unname(site_map[table_s1$accession])
raw_map <- tapply(exclusion$raw_accession_samples, exclusion$accession, unique)
table_s1$raw_accession_expression_samples <- as.integer(raw_map[table_s1$accession])
table_s1 <- table_s1[order(table_s1$accession, factor(table_s1$arm, levels = c("Active", "Placebo")),
                           table_s1$analysis_cohort), ]
rownames(table_s1) <- NULL

# Table S2 -------------------------------------------------------------------
mechanism_meta <- read_common("table_mechanism_meta_analysis.csv")
cohort_s2 <- data.frame(
  analysis_level = "Cohort",
  analysis_label = primary$cohort,
  mechanism = primary$mechanism,
  endpoint = primary$endpoint,
  model = "Firth logistic regression",
  k = 1L,
  n = primary$n,
  events = primary$events,
  non_events = primary$non_events,
  estimand = primary$estimand,
  log_odds = primary$log_odds,
  odds_ratio = primary$odds_ratio,
  ci_lower = primary$ci_lower,
  ci_upper = primary$ci_upper,
  p_value = primary$p_value,
  tau2 = NA_real_,
  I2_percent = NA_real_,
  stringsAsFactors = FALSE
)
meta_s2 <- do.call(rbind, lapply(seq_len(nrow(mechanism_meta)), function(i) {
  z <- mechanism_meta[i, ]
  idx <- if (z$group == "All biologic mechanisms") rep(TRUE, nrow(primary)) else primary$mechanism == z$group
  data.frame(
    analysis_level = ifelse(z$group == "All biologic mechanisms", "Overall", "Mechanism"),
    analysis_label = z$group,
    mechanism = ifelse(z$group == "All biologic mechanisms", "All three mechanisms", z$group),
    endpoint = "Study-specific published binary outcome",
    model = z$model,
    k = z$k,
    n = sum(primary$n[idx]),
    events = sum(primary$events[idx]),
    non_events = sum(primary$non_events[idx]),
    estimand = "Response OR per 1-SD baseline inflammatory score",
    log_odds = z$log_odds,
    odds_ratio = z$odds_ratio,
    ci_lower = z$ci_lower,
    ci_upper = z$ci_upper,
    p_value = z$p_value,
    tau2 = z$tau2,
    I2_percent = z$I2,
    stringsAsFactors = FALSE
  )
}))
table_s2 <- rbind(cohort_s2, meta_s2)

# Table S3 -------------------------------------------------------------------
arm_s3 <- data.frame(
  analysis_type = "Within-arm prognostic association",
  cohort = rct_arms$cohort,
  mechanism = rct_arms$mechanism,
  endpoint = rct_arms$endpoint,
  arm_or_contrast = rct_arms$arm,
  model = "Firth logistic regression",
  estimand = rct_arms$estimand,
  n = rct_arms$n,
  events = rct_arms$events,
  non_events = rct_arms$non_events,
  log_odds = rct_arms$log_odds,
  odds_ratio = rct_arms$odds_ratio,
  ci_lower = rct_arms$ci_lower,
  ci_upper = rct_arms$ci_upper,
  p_value = rct_arms$p_value,
  interpretation = "Arm-specific prognosis; not treatment-effect modification",
  stringsAsFactors = FALSE
)
interaction_s3 <- data.frame(
  analysis_type = "Randomized treatment-by-score interaction",
  cohort = rct_interactions$cohort,
  mechanism = rct_interactions$mechanism,
  endpoint = rct_interactions$endpoint,
  arm_or_contrast = "Active vs placebo interaction",
  model = "Randomized Firth logistic regression",
  estimand = rct_interactions$estimand,
  n = rct_interactions$n,
  events = rct_interactions$events,
  non_events = rct_interactions$n - rct_interactions$events,
  log_odds = rct_interactions$log_odds,
  odds_ratio = rct_interactions$odds_ratio,
  ci_lower = rct_interactions$ci_lower,
  ci_upper = rct_interactions$ci_upper,
  p_value = rct_interactions$p_value,
  interpretation = "Predictive effect modification not established when the 95% CI crosses 1",
  stringsAsFactors = FALSE
)
table_s3 <- rbind(arm_s3, interaction_s3)

# Table S4 -------------------------------------------------------------------
long_meta <- read_common("table_longitudinal_primary_meta.csv")
mixed <- read_common("table_longitudinal_mixed_sensitivity.csv")
module_change <- read_common("table_longitudinal_module_change_effects.csv")

primary_s4 <- data.frame(
  analysis_type = "Primary cohort change contrast",
  cohort_or_group = longitudinal$cohort,
  mechanism = longitudinal$mechanism,
  endpoint = longitudinal$endpoint,
  score = longitudinal$score,
  model = "Linear model: change ~ response",
  n = longitudinal$n,
  responders = longitudinal$responders,
  nonresponders = longitudinal$nonresponders,
  estimate = longitudinal$estimate,
  ci_lower = longitudinal$ci_lower,
  ci_upper = longitudinal$ci_upper,
  p_value = longitudinal$p_value,
  resamples = NA_integer_,
  tau2 = NA_real_, I2_percent = NA_real_,
  stringsAsFactors = FALSE
)
bootstrap_s4 <- data.frame(
  analysis_type = "Outcome-stratified bootstrap sensitivity",
  cohort_or_group = longitudinal$cohort,
  mechanism = longitudinal$mechanism,
  endpoint = longitudinal$endpoint,
  score = longitudinal$score,
  model = "Outcome-stratified nonparametric bootstrap",
  n = longitudinal$n,
  responders = longitudinal$responders,
  nonresponders = longitudinal$nonresponders,
  estimate = longitudinal$bootstrap_estimate,
  ci_lower = longitudinal$bootstrap_ci_lower,
  ci_upper = longitudinal$bootstrap_ci_upper,
  p_value = longitudinal$bootstrap_p,
  resamples = 4000L,
  tau2 = NA_real_, I2_percent = NA_real_,
  stringsAsFactors = FALSE
)
mixed_s4 <- data.frame(
  analysis_type = "Participant mixed-model sensitivity",
  cohort_or_group = mixed$cohort,
  mechanism = primary$mechanism[match(mixed$cohort, primary$cohort)],
  endpoint = primary$endpoint[match(mixed$cohort, primary$cohort)],
  score = "Inflammation_z",
  model = "Random-intercept mixed model; time x response",
  n = mixed$n,
  responders = longitudinal$responders[match(mixed$cohort, longitudinal$cohort)],
  nonresponders = longitudinal$nonresponders[match(mixed$cohort, longitudinal$cohort)],
  estimate = mixed$estimate,
  ci_lower = mixed$ci_lower,
  ci_upper = mixed$ci_upper,
  p_value = mixed$p_value_wald,
  resamples = NA_integer_,
  tau2 = NA_real_, I2_percent = NA_real_,
  singular_fit = mixed$singular,
  stringsAsFactors = FALSE
)
meta_s4 <- data.frame(
  analysis_type = "Primary overall meta-analysis",
  cohort_or_group = "Six active-treatment cohorts",
  mechanism = "anti-TNF and anti-integrin",
  endpoint = "Study-specific published binary outcome",
  score = "Inflammation_z",
  model = long_meta$model,
  n = sum(longitudinal$n),
  responders = sum(longitudinal$responders),
  nonresponders = sum(longitudinal$nonresponders),
  estimate = long_meta$estimate,
  ci_lower = long_meta$ci_lower,
  ci_upper = long_meta$ci_upper,
  p_value = long_meta$p_value,
  resamples = NA_integer_,
  tau2 = long_meta$tau2,
  I2_percent = long_meta$I2,
  stringsAsFactors = FALSE
)
module_s4 <- data.frame(
  analysis_type = "Complete prespecified module change contrast",
  cohort_or_group = module_change$cohort,
  mechanism = module_change$mechanism,
  endpoint = module_change$endpoint,
  score = module_change$score,
  model = "Linear model: change ~ response",
  n = module_change$n,
  responders = module_change$responders,
  nonresponders = module_change$nonresponders,
  estimate = module_change$estimate,
  ci_lower = module_change$ci_lower,
  ci_upper = module_change$ci_upper,
  p_value = module_change$p_value,
  resamples = NA_integer_,
  bootstrap_estimate = module_change$bootstrap_estimate,
  bootstrap_ci_lower = module_change$bootstrap_ci_lower,
  bootstrap_ci_upper = module_change$bootstrap_ci_upper,
  bootstrap_p = module_change$bootstrap_p,
  tau2 = NA_real_, I2_percent = NA_real_,
  stringsAsFactors = FALSE
)
rct_s4 <- data.frame(
  analysis_type = "Paired randomized active-placebo change contrast",
  cohort_or_group = rct_change$cohort,
  mechanism = ifelse(rct_change$cohort == "GSE92415", "anti-TNF", "anti-TNF"),
  endpoint = primary$endpoint[match(rct_change$cohort, accession_from(primary$cohort))],
  score = "Inflammation_z",
  model = "Linear model: paired change ~ randomized arm",
  n = rct_change$n,
  responders = NA_integer_, nonresponders = NA_integer_,
  estimate = rct_change$estimate,
  ci_lower = rct_change$ci_lower,
  ci_upper = rct_change$ci_upper,
  p_value = rct_change$p_value,
  resamples = NA_integer_, tau2 = NA_real_, I2_percent = NA_real_,
  active_n = rct_change$active_n,
  placebo_n = rct_change$placebo_n,
  mean_delta_active = rct_change$mean_delta_active,
  mean_delta_placebo = rct_change$mean_delta_placebo,
  stringsAsFactors = FALSE
)
table_s4 <- bind_fill(list(primary_s4, bootstrap_s4, mixed_s4, meta_s4, module_s4, rct_s4))

# Table S5 -------------------------------------------------------------------
table_s5 <- read_common("table_singlecell_compartment_outcome_sensitivity.csv")
table_s5 <- table_s5[, c(
  "analysis_label", "primary", "compartment", "aggregation", "stratum",
  "donors", "remission_donors", "nonremission_donors", "odds_ratio",
  "ci_lower", "ci_upper", "p_value", "fdr", "direction_vs_bulk", "interpretation"
)]

# Table S6 -------------------------------------------------------------------
table_s6 <- read_common("table_singlecell_state_outcome_effects.csv")
table_s6$estimability <- ifelse(is.finite(table_s6$odds_ratio), "Estimated",
                               "Not estimable under the prespecified donor-level model")
table_s6 <- table_s6[, c(
  "compartment", "state", "aggregation", "stratum", "donors",
  "remission_donors", "nonremission_donors", "odds_ratio", "ci_lower",
  "ci_upper", "p_value", "fdr", "estimability"
)]

# Table S7 -------------------------------------------------------------------
pairing <- read_common("table_singlecell_longitudinal_pairing_audit.csv")
compartment_long <- read_common("table_singlecell_compartment_longitudinal_effects.csv")
state_long <- read_common("table_singlecell_state_longitudinal_effects.csv")

pair_s7 <- data.frame(
  analysis_type = "Exact Match=Yes pairing audit",
  compartment = pairing$compartment,
  state = NA_character_,
  match_yes_pre_samples = pairing$match_yes_pre_samples,
  match_yes_post_samples = pairing$match_yes_post_samples,
  paired_donor_sites = pairing$paired_donor_sites,
  paired_donors = pairing$paired_donors,
  remission_donors = pairing$remission_donors,
  nonremission_donors = pairing$nonremission_donors,
  mean_delta_remission = NA_real_, mean_delta_nonremission = NA_real_,
  estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
  p_value = NA_real_, fdr = NA_real_,
  aggregation = "Match=Yes exact donor-site; duplicate site-time collapse; donor-equal",
  exact_pairing = pairing$exact_pairing,
  estimability = "Pairing audit; no effect estimate",
  stringsAsFactors = FALSE
)
compartment_s7 <- data.frame(
  analysis_type = "Compartment-level longitudinal effect",
  compartment = compartment_long$compartment,
  state = NA_character_,
  match_yes_pre_samples = NA_integer_, match_yes_post_samples = NA_integer_,
  paired_donor_sites = pairing$paired_donor_sites[match(compartment_long$compartment, pairing$compartment)],
  paired_donors = compartment_long$paired_donors,
  remission_donors = compartment_long$remission_donors,
  nonremission_donors = compartment_long$nonremission_donors,
  mean_delta_remission = compartment_long$mean_delta_remission,
  mean_delta_nonremission = compartment_long$mean_delta_nonremission,
  estimate = compartment_long$estimate,
  ci_lower = compartment_long$ci_lower,
  ci_upper = compartment_long$ci_upper,
  p_value = compartment_long$p_value,
  fdr = compartment_long$fdr,
  aggregation = compartment_long$aggregation,
  exact_pairing = TRUE,
  estimability = "Estimated",
  stringsAsFactors = FALSE
)
state_s7 <- data.frame(
  analysis_type = "State-level longitudinal effect",
  compartment = state_long$compartment,
  state = state_long$state,
  match_yes_pre_samples = NA_integer_, match_yes_post_samples = NA_integer_,
  paired_donor_sites = NA_integer_,
  paired_donors = state_long$paired_donors,
  remission_donors = state_long$remission_donors,
  nonremission_donors = state_long$nonremission_donors,
  mean_delta_remission = state_long$mean_delta_remission,
  mean_delta_nonremission = state_long$mean_delta_nonremission,
  estimate = state_long$estimate,
  ci_lower = state_long$ci_lower,
  ci_upper = state_long$ci_upper,
  p_value = state_long$p_value,
  fdr = state_long$fdr,
  aggregation = state_long$aggregation,
  exact_pairing = TRUE,
  estimability = ifelse(is.finite(state_long$estimate), "Estimated",
                        "Not estimable under the prespecified donor-level model"),
  stringsAsFactors = FALSE
)
table_s7 <- rbind(pair_s7, compartment_s7, state_s7)

tables <- list(
  TableS1_cohort_registry = table_s1,
  TableS2_baseline_associations = table_s2,
  TableS3_RCT_prognostic_predictive = table_s3,
  TableS4_longitudinal_sensitivities = table_s4,
  TableS5_singlecell_baseline_sensitivity = table_s5,
  TableS6_singlecell_state_baseline = table_s6,
  TableS7_singlecell_longitudinal = table_s7
)

for (nm in names(tables)) {
  utils::write.csv(tables[[nm]], file.path(output_dir, paste0(nm, ".csv")), row.names = FALSE, na = "")
}

manifest <- data.frame(
  item = paste0("Table S", seq_along(tables)),
  title = c(
    "Cohort registry and analysis inclusion",
    "Baseline cohort and mechanism associations",
    "Randomized prognosis and treatment-effect modification",
    "Longitudinal estimates and sensitivity analyses",
    "Single-cell compartment baseline sensitivity audit",
    "Single-cell state-level baseline audit",
    "Exact Match=Yes single-cell longitudinal audit"
  ),
  file = paste0(names(tables), ".csv"),
  rows = vapply(tables, nrow, integer(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(manifest, file.path(output_dir, "supplementary_tables_manifest.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(logs_dir, "20_build_supplementary_tables_sessionInfo.txt"))
print(manifest)
message("SUBMISSION TABLES S1-S7 CSV FILES EXPORTED")
