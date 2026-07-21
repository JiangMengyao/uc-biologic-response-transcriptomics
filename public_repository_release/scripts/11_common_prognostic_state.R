#!/usr/bin/env Rscript

# Cross-mechanism common prognostic-state analysis. This script treats the
# baseline mucosal inflammatory score as prognostic unless an RCT interaction
# demonstrates treatment-specific effect modification.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
derived_dir <- file.path(project_root, "data", "derived", "extended")
tables_dir <- file.path(project_root, "results", "tables", "common_state")
source_dir <- file.path(project_root, "results", "source_data")
logs_dir <- file.path(project_root, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(logistf)
  library(metafor)
  library(splines)
})

score_cols <- c("Inflammation_z", "TNFaNFKB_z", "IL6JAKSTAT3_z",
                "OxidativePhos_neg_z", "Random_neg_z")

read_scores <- function(accession) {
  readRDS(file.path(derived_dir, paste0(tolower(accession), "_all_scores.rds")))
}

collapse_subject <- function(data) {
  if (!nrow(data)) stop("Cannot collapse an empty cohort")
  split_data <- split(data, data$subject, drop = TRUE)
  rows <- lapply(split_data, function(x) {
    response <- unique(x$response[x$response %in% c("Yes", "No")])
    if (length(response) != 1L) stop("Inconsistent response within subject ", x$subject[1])
    scores <- vapply(score_cols, function(v) mean(x[[v]], na.rm = TRUE), numeric(1))
    data.frame(subject = x$subject[1], response = response,
               as.list(scores), check.names = FALSE, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$response_binary <- as.integer(out$response == "Yes")
  out
}

decorate <- function(data, accession, cohort, mechanism, endpoint, arm = "Active") {
  data$accession <- accession
  data$cohort <- cohort
  data$mechanism <- mechanism
  data$endpoint <- endpoint
  data$arm <- arm
  data
}

extract_firth <- function(data, score, cohort, mechanism, endpoint, arm,
                          estimand = "Response OR per 1-SD baseline score") {
  fit <- logistf::logistf(stats::as.formula(paste("response_binary ~", score)), data = data)
  data.frame(
    cohort = cohort, mechanism = mechanism, endpoint = endpoint, arm = arm,
    score = score, estimand = estimand, n = nrow(data),
    events = sum(data$response_binary == 1L), non_events = sum(data$response_binary == 0L),
    log_odds = unname(fit$coefficients[score]),
    odds_ratio = exp(unname(fit$coefficients[score])),
    ci_lower = exp(unname(fit$ci.lower[score])),
    ci_upper = exp(unname(fit$ci.upper[score])),
    p_value = unname(fit$prob[score]), stringsAsFactors = FALSE
  )
}

extract_interaction <- function(data, cohort, mechanism, endpoint) {
  fit <- logistf::logistf(response_binary ~ treated * Inflammation_z, data = data)
  term <- "treated:Inflammation_z"
  data.frame(
    cohort = cohort, mechanism = mechanism, endpoint = endpoint,
    score = "Inflammation_z", estimand = "Treatment-by-score interaction OR",
    n = nrow(data), events = sum(data$response_binary == 1L),
    log_odds = unname(fit$coefficients[term]),
    odds_ratio = exp(unname(fit$coefficients[term])),
    ci_lower = exp(unname(fit$ci.lower[term])),
    ci_upper = exp(unname(fit$ci.upper[term])),
    p_value = unname(fit$prob[term]), stringsAsFactors = FALSE
  )
}

add_se <- function(x) {
  x$sei <- (log(x$ci_upper) - log(x$ci_lower)) / (2 * 1.96)
  x
}

meta_row <- function(x, label, method = NULL) {
  x <- add_se(x)
  if (nrow(x) == 1L) {
    return(data.frame(group = label, model = "Single cohort", k = 1L,
                      log_odds = x$log_odds, odds_ratio = x$odds_ratio,
                      ci_lower = x$ci_lower, ci_upper = x$ci_upper,
                      p_value = x$p_value, tau2 = NA_real_, I2 = NA_real_))
  }
  if (is.null(method)) method <- if (nrow(x) >= 3L) "REML" else "FE"
  fit <- if (method == "REML") {
    metafor::rma.uni(yi = x$log_odds, sei = x$sei, method = "REML", test = "knha")
  } else {
    metafor::rma.uni(yi = x$log_odds, sei = x$sei, method = "FE")
  }
  data.frame(group = label,
             model = if (method == "REML") "REML Hartung-Knapp" else "Fixed-effect descriptive",
             k = nrow(x), log_odds = unname(fit$b[1]), odds_ratio = exp(unname(fit$b[1])),
             ci_lower = exp(fit$ci.lb), ci_upper = exp(fit$ci.ub), p_value = fit$pval,
             tau2 = ifelse(method == "REML", fit$tau2, NA_real_),
             I2 = ifelse(method == "REML", fit$I2, NA_real_))
}

# ---------- Cohort construction ----------
d16879 <- read_scores("GSE16879")
ifx16879 <- collapse_subject(d16879[d16879$visit == "Week 0" &
                                      tolower(d16879$disease) == "uc" &
                                      d16879$response %in% c("Yes", "No"), ])
ifx16879 <- decorate(ifx16879, "GSE16879", "GSE16879 IFX", "anti-TNF",
                     "Week 4–6 endoscopic response")

d73661 <- read_scores("GSE73661")
ifx73661 <- collapse_subject(d73661[d73661$visit == "W0" & d73661$treatment == "IFX" &
                                      d73661$response %in% c("Yes", "No"), ])
ifx73661 <- decorate(ifx73661, "GSE73661", "GSE73661 IFX", "anti-TNF",
                     "Week 4–6 endoscopic response")

vdz_trial <- collapse_subject(d73661[d73661$visit == "W0" &
  d73661$treatment %in% c("vdz_vdz4w", "vdz_vdz8w", "vdz_plac") &
  d73661$response %in% c("Yes", "No"), ])
vdz_trial <- decorate(vdz_trial, "GSE73661", "GSE73661 VDZ trial", "anti-integrin",
                      "Week 6 response")

vdz_obs <- collapse_subject(d73661[d73661$visit == "W0" & d73661$treatment == "vdz4w" &
                                    d73661$response %in% c("Yes", "No"), ])
vdz_obs <- decorate(vdz_obs, "GSE73661", "GSE73661 VDZ observational", "anti-integrin",
                    "Week 12 response")

d23597 <- read_scores("GSE23597")
ifx23597 <- collapse_subject(d23597[d23597$visit == "W0" & d23597$treatment == "infliximab" &
                                      d23597$response %in% c("Yes", "No"), ])
ifx23597 <- decorate(ifx23597, "GSE23597", "GSE23597 IFX", "anti-TNF",
                     "Week 8 clinical response")

d92415 <- read_scores("GSE92415")
glm92415 <- collapse_subject(d92415[d92415$visit == "Week 0" &
  grepl("Ulcerative", d92415$disease) & d92415$treatment == "golimumab" &
  d92415$response %in% c("Yes", "No"), ])
glm92415 <- decorate(glm92415, "GSE92415", "GSE92415 golimumab", "anti-TNF",
                     "Week 6 clinical response")

d206285 <- read_scores("GSE206285")
ust206285 <- collapse_subject(d206285[d206285$eligible &
  d206285$treatment == "ustekinumab" & d206285$response %in% c("Yes", "No"), ])
ust206285 <- decorate(ust206285, "GSE206285", "GSE206285 ustekinumab", "anti-IL12/23",
                      "Week 8 clinical remission")

treated <- do.call(rbind, list(ifx16879, ifx73661, ifx23597, glm92415,
                              ust206285, vdz_trial, vdz_obs))
rownames(treated) <- NULL
treated$cohort <- factor(treated$cohort, levels = unique(treated$cohort))

# ---------- Primary and supportive cohort associations ----------
effect_rows <- list()
for (cohort_name in levels(treated$cohort)) {
  x <- treated[treated$cohort == cohort_name, ]
  for (score in score_cols) {
    effect_rows[[paste(cohort_name, score)]] <- extract_firth(
      x, score, cohort_name, x$mechanism[1], x$endpoint[1], "Active"
    )
  }
}
effects <- do.call(rbind, effect_rows)
rownames(effects) <- NULL
primary <- effects[effects$score == "Inflammation_z", ]

mechanism_meta <- rbind(
  meta_row(primary[primary$mechanism == "anti-TNF", ], "anti-TNF"),
  meta_row(primary[primary$mechanism == "anti-IL12/23", ], "anti-IL12/23"),
  meta_row(primary[primary$mechanism == "anti-integrin", ], "anti-integrin", "FE"),
  meta_row(primary, "All biologic mechanisms")
)

module_meta <- do.call(rbind, lapply(score_cols, function(s) {
  out <- meta_row(effects[effects$score == s, ], s)
  out$score <- s
  out
}))

# Leave-one-cohort-out robustness for the overall primary synthesis.
loo_rows <- lapply(seq_len(nrow(primary)), function(i) {
  x <- add_se(primary[-i, ])
  fit <- metafor::rma.uni(yi = x$log_odds, sei = x$sei, method = "REML", test = "knha")
  data.frame(omitted = primary$cohort[i], k = nrow(x),
             log_odds = unname(fit$b[1]), odds_ratio = exp(unname(fit$b[1])),
             ci_lower = exp(fit$ci.lb), ci_upper = exp(fit$ci.ub),
             p_value = fit$pval, tau2 = fit$tau2, I2 = fit$I2)
})
leave_one_out <- do.call(rbind, loo_rows)

# ---------- Nonlinearity and absolute-risk curves (supportive) ----------
nonlinear_rows <- list()
prediction_rows <- list()
for (cohort_name in levels(treated$cohort)) {
  x <- treated[treated$cohort == cohort_name, ]
  linear <- glm(response_binary ~ Inflammation_z, data = x, family = binomial())
  spline_fit <- glm(response_binary ~ splines::ns(Inflammation_z, df = 3),
                    data = x, family = binomial())
  comparison <- suppressWarnings(anova(linear, spline_fit, test = "LRT"))
  p_nonlin <- comparison$`Pr(>Chi)`[2]
  nonlinear_rows[[cohort_name]] <- data.frame(
    cohort = cohort_name, mechanism = x$mechanism[1], n = nrow(x),
    events = sum(x$response_binary), p_nonlinearity = p_nonlin,
    linear_AIC = AIC(linear), spline_AIC = AIC(spline_fit)
  )
  range_x <- as.numeric(quantile(x$Inflammation_z, c(0.02, 0.98), na.rm = TRUE))
  grid <- data.frame(Inflammation_z = seq(range_x[1], range_x[2], length.out = 100L))
  pred <- predict(spline_fit, newdata = grid, type = "link", se.fit = TRUE)
  grid$probability <- plogis(pred$fit)
  grid$ci_lower <- plogis(pred$fit - 1.96 * pred$se.fit)
  grid$ci_upper <- plogis(pred$fit + 1.96 * pred$se.fit)
  grid$cohort <- cohort_name
  grid$mechanism <- x$mechanism[1]
  prediction_rows[[cohort_name]] <- grid
}
nonlinearity <- do.call(rbind, nonlinear_rows)
predictions <- do.call(rbind, prediction_rows)

# ---------- RCT prognostic-versus-predictive decomposition ----------
make_rct <- function(data, active_label, placebo_label, cohort, mechanism, endpoint) {
  active <- collapse_subject(data[data$treatment == active_label & data$response %in% c("Yes", "No"), ])
  placebo <- collapse_subject(data[data$treatment == placebo_label & data$response %in% c("Yes", "No"), ])
  active$treated <- 1L
  placebo$treated <- 0L
  all <- rbind(active, placebo)
  all$cohort <- cohort
  all$mechanism <- mechanism
  all$endpoint <- endpoint
  list(active = active, placebo = placebo, all = all)
}

rct92415 <- make_rct(d92415[d92415$visit == "Week 0" & grepl("Ulcerative", d92415$disease), ],
                     "golimumab", "placebo", "GSE92415", "anti-TNF", "Week 6 clinical response")
rct23597 <- make_rct(d23597[d23597$visit == "W0", ], "infliximab", "placebo",
                     "GSE23597", "anti-TNF", "Week 8 clinical response")
rct206285 <- make_rct(d206285[d206285$eligible, ], "ustekinumab", "placebo",
                      "GSE206285", "anti-IL12/23", "Week 8 clinical remission")
rct_list <- list(GSE92415 = rct92415, GSE23597 = rct23597, GSE206285 = rct206285)

rct_arm_rows <- list()
rct_interaction_rows <- list()
for (nm in names(rct_list)) {
  z <- rct_list[[nm]]
  mechanism <- z$all$mechanism[1]
  endpoint <- z$all$endpoint[1]
  rct_arm_rows[[paste0(nm, "_active")]] <- extract_firth(
    z$active, "Inflammation_z", nm, mechanism, endpoint, "Active")
  rct_arm_rows[[paste0(nm, "_placebo")]] <- extract_firth(
    z$placebo, "Inflammation_z", nm, mechanism, endpoint, "Placebo")
  rct_interaction_rows[[nm]] <- extract_interaction(z$all, nm, mechanism, endpoint)
}
rct_arms <- do.call(rbind, rct_arm_rows)
rct_interactions <- do.call(rbind, rct_interaction_rows)

# ---------- Locked baseline gate ----------
anti_meta <- mechanism_meta[mechanism_meta$group == "anti-TNF", ]
ust_meta <- mechanism_meta[mechanism_meta$group == "anti-IL12/23", ]
vdz_meta <- mechanism_meta[mechanism_meta$group == "anti-integrin", ]
gate <- data.frame(
  mechanism = c("anti-TNF", "anti-IL12/23", "anti-integrin"),
  odds_ratio = c(anti_meta$odds_ratio, ust_meta$odds_ratio, vdz_meta$odds_ratio),
  criterion = "Mechanism-level response OR per 1-SD higher inflammation < 1",
  passed = c(anti_meta$odds_ratio, ust_meta$odds_ratio, vdz_meta$odds_ratio) < 1
)

# ---------- Outputs ----------
utils::write.csv(treated, file.path(source_dir, "common_state_baseline_patient_data.csv"), row.names = FALSE)
utils::write.csv(effects, file.path(tables_dir, "table_cohort_module_associations.csv"), row.names = FALSE)
utils::write.csv(primary, file.path(tables_dir, "table_primary_cross_mechanism_associations.csv"), row.names = FALSE)
utils::write.csv(mechanism_meta, file.path(tables_dir, "table_mechanism_meta_analysis.csv"), row.names = FALSE)
utils::write.csv(module_meta, file.path(tables_dir, "table_module_meta_analysis.csv"), row.names = FALSE)
utils::write.csv(leave_one_out, file.path(tables_dir, "table_leave_one_out.csv"), row.names = FALSE)
utils::write.csv(nonlinearity, file.path(tables_dir, "table_nonlinearity_tests.csv"), row.names = FALSE)
utils::write.csv(predictions, file.path(source_dir, "common_state_nonlinear_predictions.csv"), row.names = FALSE)
utils::write.csv(rct_arms, file.path(tables_dir, "table_rct_arm_prognostic_associations.csv"), row.names = FALSE)
utils::write.csv(rct_interactions, file.path(tables_dir, "table_rct_treatment_score_interactions.csv"), row.names = FALSE)
utils::write.csv(gate, file.path(tables_dir, "table_cross_mechanism_gate.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "11_common_prognostic_state_sessionInfo.txt"))

print(primary)
print(mechanism_meta)
print(rct_arms)
print(rct_interactions)
print(gate)
if (!all(gate$passed)) {
  stop("CROSS-MECHANISM BASELINE GATE FAILED: ",
       paste(gate$mechanism[!gate$passed], collapse = ", "))
}
message("CROSS-MECHANISM COMMON-PROGNOSTIC-STATE GATE PASSED")
