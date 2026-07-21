#!/usr/bin/env Rscript

# Paired longitudinal analysis of whether clinical responders show greater
# reversal of the baseline inflammatory state. The inferential unit is the
# participant, with one collapsed record per subject and visit.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
derived_dir <- file.path(project_root, "data", "derived", "extended")
tables_dir <- file.path(project_root, "results", "tables", "common_state")
source_dir <- file.path(project_root, "results", "source_data")
logs_dir <- file.path(project_root, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(lme4)
  library(logistf)
  library(metafor)
})

set.seed(20260720)
score_cols <- c("Inflammation_z", "TNFaNFKB_z", "IL6JAKSTAT3_z",
                "OxidativePhos_neg_z", "Random_neg_z")

read_scores <- function(accession) {
  readRDS(file.path(derived_dir, paste0(tolower(accession), "_all_scores.rds")))
}

collapse_subject_time <- function(data) {
  split_data <- split(data, interaction(data$subject, data$visit, drop = TRUE))
  rows <- lapply(split_data, function(x) {
    response <- unique(x$response[x$response %in% c("Yes", "No")])
    if (length(response) != 1L) stop("Inconsistent response for ", x$subject[1], " at ", x$visit[1])
    treatment <- unique(x$treatment)
    if (length(treatment) != 1L) stop("Inconsistent treatment for ", x$subject[1], " at ", x$visit[1])
    scores <- vapply(score_cols, function(v) mean(x[[v]], na.rm = TRUE), numeric(1))
    data.frame(subject = x$subject[1], visit = x$visit[1], response = response,
               treatment = treatment, as.list(scores), check.names = FALSE,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

pair_visits <- function(data, baseline_visit, followup_visit, cohort, mechanism,
                        endpoint, active_treatments, placebo_treatments = character()) {
  x <- collapse_subject_time(data[data$visit %in% c(baseline_visit, followup_visit) &
                                  data$response %in% c("Yes", "No"), ])
  baseline <- x[x$visit == baseline_visit, ]
  followup <- x[x$visit == followup_visit, ]
  keep_cols <- c("subject", "response", "treatment", score_cols)
  paired <- merge(baseline[, keep_cols], followup[, keep_cols], by = "subject",
                  suffixes = c("_baseline", "_followup"))
  if (any(paired$response_baseline != paired$response_followup)) {
    stop("Response mismatch across visits in ", cohort)
  }
  paired$response <- paired$response_baseline
  paired$response_binary <- as.integer(paired$response == "Yes")
  paired$arm <- ifelse(paired$treatment_baseline %in% active_treatments, "Active",
                       ifelse(paired$treatment_baseline %in% placebo_treatments, "Placebo", NA_character_))
  if (anyNA(paired$arm)) stop("Unclassified treatment arm in ", cohort)
  paired$cohort <- cohort
  paired$mechanism <- mechanism
  paired$endpoint <- endpoint
  for (s in score_cols) paired[[paste0("delta_", s)]] <- paired[[paste0(s, "_followup")]] - paired[[paste0(s, "_baseline")]]

  # Thresholds are calculated once from all paired baseline participants in the
  # cohort, before looking at response or randomized arm.
  median_base <- median(paired$Inflammation_z_baseline, na.rm = TRUE)
  high_cut <- as.numeric(quantile(paired$Inflammation_z_baseline, 2/3, na.rm = TRUE, type = 7))
  paired$baseline_tertile <- cut(
    paired$Inflammation_z_baseline,
    breaks = c(-Inf, as.numeric(quantile(paired$Inflammation_z_baseline, 1/3, na.rm = TRUE)), high_cut, Inf),
    labels = c("Low", "Middle", "High"), include.lowest = TRUE
  )
  paired$baseline_high <- paired$Inflammation_z_baseline >= high_cut
  paired$molecular_reversal <- paired$Inflammation_z_followup < median_base
  paired$baseline_median <- median_base
  paired$baseline_high_cut <- high_cut
  paired
}

bootstrap_mean_difference <- function(data, variable, n_boot = 4000L) {
  yes <- data[data$response_binary == 1L, variable]
  no <- data[data$response_binary == 0L, variable]
  observed <- mean(yes) - mean(no)
  boot <- replicate(n_boot, mean(sample(yes, length(yes), replace = TRUE)) -
                      mean(sample(no, length(no), replace = TRUE)))
  c(estimate = observed,
    ci_lower = unname(quantile(boot, 0.025, type = 6)),
    ci_upper = unname(quantile(boot, 0.975, type = 6)),
    bootstrap_p = 2 * min(mean(boot <= 0), mean(boot >= 0)))
}

change_effect <- function(data, score) {
  variable <- paste0("delta_", score)
  fit <- lm(stats::as.formula(paste(variable, "~ response_binary")), data = data)
  ci <- confint(fit, "response_binary", level = 0.95)
  boot <- bootstrap_mean_difference(data, variable)
  data.frame(
    cohort = data$cohort[1], mechanism = data$mechanism[1], endpoint = data$endpoint[1],
    score = score, n = nrow(data), responders = sum(data$response_binary == 1L),
    nonresponders = sum(data$response_binary == 0L),
    mean_delta_responders = mean(data[data$response_binary == 1L, variable]),
    mean_delta_nonresponders = mean(data[data$response_binary == 0L, variable]),
    estimate = unname(coef(fit)["response_binary"]),
    ci_lower = ci[1], ci_upper = ci[2], p_value = coef(summary(fit))["response_binary", "Pr(>|t|)"],
    bootstrap_estimate = boot["estimate"], bootstrap_ci_lower = boot["ci_lower"],
    bootstrap_ci_upper = boot["ci_upper"], bootstrap_p = boot["bootstrap_p"]
  )
}

mixed_effect <- function(data) {
  long <- rbind(
    data.frame(subject = data$subject, cohort = data$cohort, mechanism = data$mechanism,
               response_binary = data$response_binary, time = 0,
               Inflammation_z = data$Inflammation_z_baseline),
    data.frame(subject = data$subject, cohort = data$cohort, mechanism = data$mechanism,
               response_binary = data$response_binary, time = 1,
               Inflammation_z = data$Inflammation_z_followup)
  )
  fit <- lme4::lmer(Inflammation_z ~ time * response_binary + (1 | subject), data = long,
                    control = lmerControl(check.conv.singular = "ignore"))
  term <- "time:response_binary"
  se <- sqrt(diag(vcov(fit)))[term]
  estimate <- fixef(fit)[term]
  data.frame(cohort = data$cohort[1], n = nrow(data), term = term,
             estimate = estimate, se = se, ci_lower = estimate - 1.96 * se,
             ci_upper = estimate + 1.96 * se,
             p_value_wald = 2 * pnorm(-abs(estimate / se)), singular = isSingular(fit))
}

active_only <- function(x) x[x$arm == "Active", ]

# ---------- Build six treated longitudinal cohorts ----------
d16879 <- read_scores("GSE16879")
p16879 <- pair_visits(d16879[tolower(d16879$disease) == "uc", ], "Week 0", "Week 4-6",
                      "GSE16879 IFX", "anti-TNF", "Week 4–6 endoscopic response",
                      active_treatments = "infliximab")

d73661 <- read_scores("GSE73661")
p73661_ifx <- pair_visits(d73661[d73661$treatment == "IFX", ], "W0", "W4_W6",
                          "GSE73661 IFX", "anti-TNF", "Week 4–6 endoscopic response",
                          active_treatments = "IFX")
p73661_vdz_trial <- pair_visits(d73661[d73661$treatment %in% c("vdz_vdz4w", "vdz_vdz8w", "vdz_plac"), ],
                                "W0", "W6", "GSE73661 VDZ trial", "anti-integrin",
                                "Week 6 response",
                                active_treatments = c("vdz_vdz4w", "vdz_vdz8w", "vdz_plac"))
p73661_vdz_obs <- pair_visits(d73661[d73661$treatment == "vdz4w", ], "W0", "W12",
                              "GSE73661 VDZ observational", "anti-integrin",
                              "Week 12 response", active_treatments = "vdz4w")

d23597 <- read_scores("GSE23597")
p23597 <- pair_visits(d23597, "W0", "W8", "GSE23597", "anti-TNF",
                      "Week 8 clinical response", active_treatments = "infliximab",
                      placebo_treatments = "placebo")

d92415 <- read_scores("GSE92415")
p92415 <- pair_visits(d92415[grepl("Ulcerative", d92415$disease), ], "Week 0", "Week 6",
                      "GSE92415", "anti-TNF", "Week 6 clinical response",
                      active_treatments = "golimumab", placebo_treatments = "placebo")

treated_cohorts <- list(
  "GSE16879 IFX" = active_only(p16879),
  "GSE73661 IFX" = active_only(p73661_ifx),
  "GSE23597 IFX" = active_only(transform(p23597, cohort = ifelse(arm == "Active", "GSE23597 IFX", "GSE23597 placebo"))),
  "GSE92415 golimumab" = active_only(transform(p92415, cohort = ifelse(arm == "Active", "GSE92415 golimumab", "GSE92415 placebo"))),
  "GSE73661 VDZ trial" = active_only(p73661_vdz_trial),
  "GSE73661 VDZ observational" = active_only(p73661_vdz_obs)
)

paired_all <- do.call(rbind, treated_cohorts)
rownames(paired_all) <- NULL

# ---------- Change effects and mixed-model sensitivity ----------
effect_rows <- list()
mixed_rows <- list()
for (nm in names(treated_cohorts)) {
  x <- treated_cohorts[[nm]]
  for (score in score_cols) effect_rows[[paste(nm, score)]] <- change_effect(x, score)
  mixed_rows[[nm]] <- mixed_effect(x)
}
effects <- do.call(rbind, effect_rows)
rownames(effects) <- NULL
primary <- effects[effects$score == "Inflammation_z", ]
mixed <- do.call(rbind, mixed_rows)
rownames(mixed) <- NULL

primary$sei <- (primary$ci_upper - primary$ci_lower) / (2 * 1.96)
meta_fit <- metafor::rma.uni(yi = primary$estimate, sei = primary$sei,
                             method = "REML", test = "knha")
meta <- data.frame(
  model = "REML Hartung-Knapp", k = nrow(primary), estimate = unname(meta_fit$b[1]),
  ci_lower = meta_fit$ci.lb, ci_upper = meta_fit$ci.ub, p_value = meta_fit$pval,
  tau2 = meta_fit$tau2, I2 = meta_fit$I2
)

# Active-versus-placebo paired change in the two RCTs with post-treatment tissue.
rct_change <- lapply(list(GSE92415 = p92415, GSE23597 = p23597), function(x) {
  x$treated <- as.integer(x$arm == "Active")
  fit <- lm(delta_Inflammation_z ~ treated, data = x)
  ci <- confint(fit, "treated")
  data.frame(cohort = x$cohort[1], n = nrow(x), active_n = sum(x$treated), placebo_n = sum(!x$treated),
             estimate = coef(fit)["treated"], ci_lower = ci[1], ci_upper = ci[2],
             p_value = coef(summary(fit))["treated", "Pr(>|t|)"],
             mean_delta_active = mean(x$delta_Inflammation_z[x$treated == 1L]),
             mean_delta_placebo = mean(x$delta_Inflammation_z[x$treated == 0L]))
})
rct_change <- do.call(rbind, rct_change)
rownames(rct_change) <- NULL

# ---------- Locked high-state transition analysis ----------
high <- paired_all[paired_all$baseline_high, ]
transition_counts <- aggregate(
  rep(1L, nrow(high)),
  list(cohort = high$cohort, mechanism = high$mechanism, response = high$response,
       molecular_reversal = high$molecular_reversal), sum
)
names(transition_counts)[5] <- "n"

transition_summary <- do.call(rbind, lapply(split(high, high$cohort), function(x) {
  do.call(rbind, lapply(split(x, x$response), function(z) {
    data.frame(cohort = z$cohort[1], mechanism = z$mechanism[1], response = z$response[1],
               high_baseline_n = nrow(z), reversal_n = sum(z$molecular_reversal),
               reversal_rate = mean(z$molecular_reversal))
  }))
}))
rownames(transition_summary) <- NULL

transition_or <- lapply(split(high, high$cohort), function(x) {
  if (length(unique(x$response_binary)) < 2L || length(unique(x$molecular_reversal)) < 2L) {
    return(data.frame(cohort = x$cohort[1], n = nrow(x), odds_ratio = NA_real_,
                      ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_))
  }
  fit <- logistf::logistf(response_binary ~ molecular_reversal, data = x)
  term <- "molecular_reversalTRUE"
  data.frame(cohort = x$cohort[1], n = nrow(x), odds_ratio = exp(fit$coefficients[term]),
             ci_lower = exp(fit$ci.lower[term]), ci_upper = exp(fit$ci.upper[term]),
             p_value = fit$prob[term])
})
transition_or <- do.call(rbind, transition_or)
rownames(transition_or) <- NULL

# ---------- Locked longitudinal gate ----------
direction_pass_n <- sum(primary$estimate < 0)
class_direction <- aggregate(estimate ~ mechanism, primary, function(v) {
  if (length(v) >= 2L && all(v >= 0)) "coherent_opposite" else "not_coherent_opposite"
})
gate <- data.frame(
  gate = c("At least two treated cohorts have responder-minus-nonresponder delta < 0",
           "No mechanism class with >=2 cohorts is coherently opposite"),
  observed = c(paste0(direction_pass_n, " of ", nrow(primary), " cohorts"),
               paste(class_direction$mechanism[class_direction$estimate == "coherent_opposite"], collapse = "; ")),
  passed = c(direction_pass_n >= 2L, !any(class_direction$estimate == "coherent_opposite"))
)
if (!nzchar(gate$observed[2])) gate$observed[2] <- "none"

# ---------- Outputs ----------
utils::write.csv(paired_all, file.path(source_dir, "longitudinal_paired_patient_data.csv"), row.names = FALSE)
utils::write.csv(effects, file.path(tables_dir, "table_longitudinal_module_change_effects.csv"), row.names = FALSE)
utils::write.csv(primary, file.path(tables_dir, "table_longitudinal_primary_change_effects.csv"), row.names = FALSE)
utils::write.csv(meta, file.path(tables_dir, "table_longitudinal_primary_meta.csv"), row.names = FALSE)
utils::write.csv(mixed, file.path(tables_dir, "table_longitudinal_mixed_sensitivity.csv"), row.names = FALSE)
utils::write.csv(rct_change, file.path(tables_dir, "table_rct_active_placebo_change.csv"), row.names = FALSE)
utils::write.csv(transition_counts, file.path(tables_dir, "table_high_state_transition_counts.csv"), row.names = FALSE)
utils::write.csv(transition_summary, file.path(tables_dir, "table_high_state_transition_summary.csv"), row.names = FALSE)
utils::write.csv(transition_or, file.path(tables_dir, "table_high_state_reversal_response_or.csv"), row.names = FALSE)
utils::write.csv(gate, file.path(tables_dir, "table_longitudinal_gate.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "12_longitudinal_overcoming_sessionInfo.txt"))

print(primary)
print(meta)
print(rct_change)
print(transition_summary)
print(gate)
if (!all(gate$passed)) {
  stop("LONGITUDINAL DIRECTION GATE FAILED. Stop before single-cell and main-figure generation.")
}
message("LONGITUDINAL REVERSAL GATE PASSED")
