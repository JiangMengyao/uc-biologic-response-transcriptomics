#!/usr/bin/env Rscript

# Prespecified confounding and heterogeneity audit for the shared baseline
# inflammatory-state analysis. Only unambiguously participant-linked public
# covariates are eligible for adjustment. The script writes all audit outputs
# before enforcing the stop gates.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
raw_dir <- file.path(project_root, "data", "raw")
derived_dir <- file.path(project_root, "data", "derived", "extended")
tables_dir <- file.path(project_root, "results", "tables", "common_state")
source_dir <- file.path(project_root, "results", "source_data")
logs_dir <- file.path(project_root, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(logistf)
  library(metafor)
})

baseline_path <- file.path(source_dir, "common_state_baseline_patient_data.csv")
primary_path <- file.path(tables_dir, "table_primary_cross_mechanism_associations.csv")
mechanism_path <- file.path(tables_dir, "table_mechanism_meta_analysis.csv")
loo_path <- file.path(tables_dir, "table_leave_one_out.csv")
metadata206_path <- file.path(raw_dir, "GSE206285_additional_sample_metadata.tsv.gz")
required <- c(baseline_path, primary_path, mechanism_path, loo_path, metadata206_path)
if (any(!file.exists(required))) stop("Missing required audit input: ", paste(basename(required[!file.exists(required)]), collapse = ", "))

baseline <- utils::read.csv(baseline_path, stringsAsFactors = FALSE, check.names = FALSE)
primary <- utils::read.csv(primary_path, stringsAsFactors = FALSE, check.names = FALSE)
mechanism_meta <- utils::read.csv(mechanism_path, stringsAsFactors = FALSE, check.names = FALSE)
leave_one_out <- utils::read.csv(loo_path, stringsAsFactors = FALSE, check.names = FALSE)
baseline$subject <- as.character(baseline$subject)

read_scores <- function(accession) {
  readRDS(file.path(derived_dir, paste0(tolower(accession), "_all_scores.rds")))
}

read_pdata <- function(accession) {
  object <- GEOquery::getGEO(
    filename = file.path(raw_dir, paste0(accession, "_series_matrix.txt.gz")),
    getGPL = FALSE
  )
  if (is.list(object)) object <- object[[1]]
  Biobase::pData(object)
}

as_number <- function(x) suppressWarnings(as.numeric(ifelse(x %in% c("", ".", "NA"), NA, x)))

fill_by_key <- function(target, rows, key, source, destination) {
  map <- rows[!duplicated(rows[[key]]), c(key, source), drop = FALSE]
  idx <- match(target[[key]], map[[key]])
  target[[destination]] <- map[[source]][idx]
  target
}

clinical <- baseline
clinical$age_years <- NA_real_
clinical$baseline_activity <- NA_real_
clinical$baseline_activity_scale <- NA_character_
clinical$disease_duration_years <- NA_real_
clinical$prior_anti_tnf_exposure <- NA_real_

# GSE92415: age and total Mayo are complete for the active baseline cohort.
p92415 <- read_pdata("GSE92415")
m92415 <- data.frame(
  subject = as.character(p92415[["subject:ch1"]]),
  visit = p92415[["visit:ch1"]],
  treatment = tolower(p92415[["treatment:ch1"]]),
  age_years = as_number(p92415[["age:ch1"]]),
  baseline_activity = as_number(p92415[["mayo score:ch1"]]),
  stringsAsFactors = FALSE
)
m92415 <- m92415[m92415$visit == "Week 0" & m92415$treatment == "golimumab", ]
idx <- clinical$cohort == "GSE92415 golimumab"
clinical$age_years[idx] <- m92415$age_years[match(clinical$subject[idx], m92415$subject)]
clinical$baseline_activity[idx] <- m92415$baseline_activity[match(clinical$subject[idx], m92415$subject)]
clinical$baseline_activity_scale[idx] <- "Total Mayo (0-12)"

# GSE73661: baseline endoscopic Mayo subscore is participant-linked.
p73661 <- read_pdata("GSE73661")
m73661 <- data.frame(
  subject = as.character(p73661[["study individual number:ch1"]]),
  visit = p73661[["week (w):ch1"]],
  treatment = p73661[["induction therapy_maintenance therapy:ch1"]],
  baseline_activity = as_number(p73661[["mayo endoscopic subscore:ch1"]]),
  stringsAsFactors = FALSE
)
m73661 <- m73661[m73661$visit == "W0", ]
m73661$cohort <- ifelse(
  m73661$treatment == "IFX", "GSE73661 IFX",
  ifelse(m73661$treatment == "vdz4w", "GSE73661 VDZ observational",
         ifelse(m73661$treatment %in% c("vdz_vdz4w", "vdz_vdz8w", "vdz_plac"),
                "GSE73661 VDZ trial", NA_character_))
)
for (cohort_name in unique(na.omit(m73661$cohort))) {
  idx <- clinical$cohort == cohort_name
  map <- m73661[m73661$cohort == cohort_name, ]
  clinical$baseline_activity[idx] <- map$baseline_activity[match(clinical$subject[idx], map$subject)]
  clinical$baseline_activity_scale[idx] <- "Mayo endoscopic subscore (0-3)"
}

# GSE206285: GEO-provided additional participant metadata.
p206285 <- read_pdata("GSE206285")
extra206 <- utils::read.delim(
  gzfile(metadata206_path), stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA", ".")
)
m206285 <- data.frame(
  gsm = as.character(extra206[["Accession"]]),
  baseline_activity = as_number(extra206[["Field[Mayo_score]"]]),
  disease_duration_years = as_number(extra206[["Field[Disease_duration]"]]),
  prior_anti_tnf_exposure = ifelse(extra206[["Field[Anti_TNF_history]"]] == "Y", 1,
                                   ifelse(extra206[["Field[Anti_TNF_history]"]] == "N", 0, NA)),
  stringsAsFactors = FALSE
)
m206285$subject <- as.character(p206285[["donor id:ch1"]][match(m206285$gsm, rownames(p206285))])
idx <- clinical$cohort == "GSE206285 ustekinumab"
clinical$baseline_activity[idx] <- m206285$baseline_activity[match(clinical$subject[idx], m206285$subject)]
clinical$disease_duration_years[idx] <- m206285$disease_duration_years[match(clinical$subject[idx], m206285$subject)]
clinical$prior_anti_tnf_exposure[idx] <- m206285$prior_anti_tnf_exposure[match(clinical$subject[idx], m206285$subject)]
clinical$baseline_activity_scale[idx] <- "Total Mayo (5-12 in released metadata)"

cohort_n <- setNames(primary$n, primary$cohort)
cohorts <- as.character(primary$cohort)
requested_covariates <- c(
  "Age", "Sex", "Disease duration", "Baseline disease activity",
  "Corticosteroid use", "Prior biologic failure", "Disease extent", "Biopsy site"
)

availability <- expand.grid(
  cohort = cohorts, covariate = requested_covariates,
  stringsAsFactors = FALSE
)
availability$n_total <- as.integer(cohort_n[availability$cohort])
availability$n_individual_nonmissing <- 0L
availability$missing_n <- availability$n_total
availability$availability <- "Not publicly released at participant level"
availability$source_field <- NA_character_
availability$within_cohort_variation <- "Unavailable"
availability$eligible_for_adjustment <- FALSE
availability$used_in_adjustment <- FALSE
availability$note <- "Not present in the GEO series matrix or accession-level supplementary metadata"

set_availability <- function(cohort_names, covariate, n_nonmissing, status, source_field,
                             variation, eligible, used, note) {
  for (cohort_name in cohort_names) {
    ii <- availability$cohort == cohort_name & availability$covariate == covariate
    nn <- if (length(n_nonmissing) == 1L && is.character(n_nonmissing)) {
      sum(!is.na(clinical[[n_nonmissing]][clinical$cohort == cohort_name]))
    } else if (length(n_nonmissing) == 1L) as.integer(n_nonmissing) else as.integer(n_nonmissing[[cohort_name]])
    availability$n_individual_nonmissing[ii] <<- nn
    availability$missing_n[ii] <<- availability$n_total[ii] - nn
    availability$availability[ii] <<- status
    availability$source_field[ii] <<- source_field
    availability$within_cohort_variation[ii] <<- variation
    availability$eligible_for_adjustment[ii] <<- eligible
    availability$used_in_adjustment[ii] <<- used
    availability$note[ii] <<- note
  }
}

set_availability("GSE92415 golimumab", "Age", "age_years", "Individual-level complete",
                 "GSE92415 age:ch1", "Yes", TRUE, TRUE, "Age was included as a precision/confounding adjustment")
set_availability("GSE92415 golimumab", "Baseline disease activity", "baseline_activity",
                 "Individual-level complete", "GSE92415 mayo score:ch1", "Yes", TRUE, TRUE,
                 "Total Mayo score")
set_availability(c("GSE73661 IFX", "GSE73661 VDZ trial", "GSE73661 VDZ observational"),
                 "Baseline disease activity", "baseline_activity", "Individual-level complete",
                 "GSE73661 mayo endoscopic subscore:ch1", "Yes", TRUE, TRUE,
                 "Endoscopic Mayo subscore; not interchangeable with total Mayo")
set_availability("GSE206285 ustekinumab", "Baseline disease activity", "baseline_activity",
                 "Individual-level complete", "GSE206285 additional metadata Field[Mayo_score]", "Yes", TRUE, TRUE,
                 "Total Mayo score")
set_availability("GSE206285 ustekinumab", "Disease duration", "disease_duration_years",
                 "Individual-level complete", "GSE206285 additional metadata Field[Disease_duration]", "Yes", TRUE, TRUE,
                 "Released in years")
set_availability("GSE206285 ustekinumab", "Prior biologic failure", "prior_anti_tnf_exposure",
                 "Proxy only", "GSE206285 additional metadata Field[Anti_TNF_history]", "Yes", TRUE, TRUE,
                 "Prior anti-TNF exposure, not documented treatment failure; model and report use exposure wording")
set_availability("GSE16879 IFX", "Baseline disease activity", 0L, "Study-level eligibility only",
                 "GSE16879 overall design", "No", FALSE, FALSE,
                 "Biopsies were from actively inflamed mucosa; no participant-level severity score")
set_availability("GSE23597 IFX", "Baseline disease activity", 0L, "Study-level eligibility only",
                 "GSE23597 overall design", "No", FALSE, FALSE,
                 "Moderate-to-severe active UC; no participant-level severity score")
set_availability("GSE16879 IFX", "Corticosteroid use", 0L, "Study-level eligibility only",
                 "GSE16879 overall design", "No", FALSE, FALSE,
                 "Refractory to corticosteroids and/or immunosuppression; current individual use unavailable")

distal_cohorts <- c("GSE23597 IFX", "GSE92415 golimumab", "GSE206285 ustekinumab")
unspecified_cohorts <- setdiff(cohorts, distal_cohorts)
set_availability(distal_cohorts, "Biopsy site", 0L, "Study-level constant",
                 "GEO series design/protocol", "No", FALSE, FALSE,
                 "Protocol-fixed distal/sigmoid colon at approximately 15-20 cm from the anal verge")
set_availability(unspecified_cohorts, "Biopsy site", 0L, "Study-level constant; exact site not released",
                 "GEO sample/series metadata", "No", FALSE, FALSE,
                 "Colonic mucosa; exact participant-level location unavailable")

# Model effects --------------------------------------------------------------
extract_effect <- function(fit, term = "Inflammation_z") {
  data.frame(
    log_odds = unname(fit$coefficients[term]),
    odds_ratio = exp(unname(fit$coefficients[term])),
    ci_lower = exp(unname(fit$ci.lower[term])),
    ci_upper = exp(unname(fit$ci.upper[term])),
    p_value = unname(fit$prob[term]), stringsAsFactors = FALSE
  )
}

z_if_variable <- function(x) {
  spread <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(spread) || spread <= 0) stop("Adjustment covariate has no variation")
  as.numeric(scale(x))
}

model_specs <- list(
  "GSE73661 IFX" = c("baseline_activity_z"),
  "GSE73661 VDZ trial" = c("baseline_activity_z"),
  "GSE73661 VDZ observational" = c("baseline_activity_z"),
  "GSE92415 golimumab" = c("age_z", "baseline_activity_z"),
  "GSE206285 ustekinumab" = c("baseline_activity_z", "disease_duration_z", "prior_anti_tnf_exposure")
)
model_labels <- list(
  "GSE73661 IFX" = "Baseline Mayo endoscopic subscore",
  "GSE73661 VDZ trial" = "Baseline Mayo endoscopic subscore",
  "GSE73661 VDZ observational" = "Baseline Mayo endoscopic subscore",
  "GSE92415 golimumab" = "Age + baseline total Mayo",
  "GSE206285 ustekinumab" = "Baseline total Mayo + disease duration + prior anti-TNF exposure"
)

adjusted_rows <- list()
for (cohort_name in cohorts) {
  x <- clinical[clinical$cohort == cohort_name, ]
  original <- primary[primary$cohort == cohort_name, ]
  base_row <- data.frame(
    cohort = cohort_name, mechanism = original$mechanism, endpoint = original$endpoint,
    n_primary = original$n, events_primary = original$events,
    unadjusted_or = original$odds_ratio, unadjusted_ci_lower = original$ci_lower,
    unadjusted_ci_upper = original$ci_upper, unadjusted_p = original$p_value,
    adjustment_status = ifelse(cohort_name %in% names(model_specs), "Estimated", "Not estimable from public participant-level covariates"),
    covariates = ifelse(cohort_name %in% names(model_specs), model_labels[[cohort_name]], NA_character_),
    n_complete = NA_integer_, events_complete = NA_integer_,
    complete_case_unadjusted_or = NA_real_, complete_case_unadjusted_ci_lower = NA_real_,
    complete_case_unadjusted_ci_upper = NA_real_, complete_case_unadjusted_p = NA_real_,
    adjusted_log_odds = NA_real_, adjusted_or = NA_real_, adjusted_ci_lower = NA_real_,
    adjusted_ci_upper = NA_real_, adjusted_p = NA_real_, direction_concordant = NA,
    relative_log_effect_change = NA_real_, stringsAsFactors = FALSE
  )
  if (cohort_name %in% names(model_specs)) {
    if (any(!is.na(x$age_years))) x$age_z <- z_if_variable(x$age_years)
    if (any(!is.na(x$baseline_activity))) x$baseline_activity_z <- z_if_variable(x$baseline_activity)
    if (any(!is.na(x$disease_duration_years))) x$disease_duration_z <- z_if_variable(x$disease_duration_years)
    vars <- c("response_binary", "Inflammation_z", model_specs[[cohort_name]])
    cc <- x[stats::complete.cases(x[, vars, drop = FALSE]), ]
    if (nrow(cc) < 10L || length(unique(cc$response_binary)) < 2L) stop("Adjusted model gate failed for ", cohort_name)
    unadjusted_fit <- logistf::logistf(response_binary ~ Inflammation_z, data = cc)
    formula <- stats::as.formula(paste("response_binary ~ Inflammation_z +", paste(model_specs[[cohort_name]], collapse = " + ")))
    adjusted_fit <- logistf::logistf(formula, data = cc)
    u <- extract_effect(unadjusted_fit)
    a <- extract_effect(adjusted_fit)
    base_row$n_complete <- nrow(cc)
    base_row$events_complete <- sum(cc$response_binary)
    base_row$complete_case_unadjusted_or <- u$odds_ratio
    base_row$complete_case_unadjusted_ci_lower <- u$ci_lower
    base_row$complete_case_unadjusted_ci_upper <- u$ci_upper
    base_row$complete_case_unadjusted_p <- u$p_value
    base_row$adjusted_log_odds <- a$log_odds
    base_row$adjusted_or <- a$odds_ratio
    base_row$adjusted_ci_lower <- a$ci_lower
    base_row$adjusted_ci_upper <- a$ci_upper
    base_row$adjusted_p <- a$p_value
    base_row$direction_concordant <- sign(a$log_odds) == sign(u$log_odds)
    base_row$relative_log_effect_change <- (a$log_odds - u$log_odds) / abs(u$log_odds)
  }
  adjusted_rows[[cohort_name]] <- base_row
}
adjusted <- do.call(rbind, adjusted_rows)
rownames(adjusted) <- NULL

add_se <- function(x, log_col = "log_odds", lower_col = "ci_lower", upper_col = "ci_upper") {
  x$sei <- (log(x[[upper_col]]) - log(x[[lower_col]])) / (2 * 1.96)
  x$yi <- x[[log_col]]
  x
}

meta_summary <- function(x, label, model_preference = NULL, log_col = "log_odds",
                         lower_col = "ci_lower", upper_col = "ci_upper") {
  x <- add_se(x, log_col, lower_col, upper_col)
  if (nrow(x) == 1L) {
    return(data.frame(group = label, model = "Single cohort", k = 1L,
                      n_total = sum(x$n), events = sum(x$events), log_odds = x$yi,
                      odds_ratio = exp(x$yi), ci_lower = x[[lower_col]], ci_upper = x[[upper_col]],
                      p_value = x$p_value, tau2 = NA_real_, I2 = NA_real_, stringsAsFactors = FALSE))
  }
  method <- if (is.null(model_preference)) ifelse(nrow(x) >= 3L, "REML", "FE") else model_preference
  fit <- if (method == "REML") {
    metafor::rma.uni(yi = x$yi, sei = x$sei, method = "REML", test = "knha")
  } else metafor::rma.uni(yi = x$yi, sei = x$sei, method = "FE")
  data.frame(group = label,
             model = ifelse(method == "REML", "REML Hartung-Knapp", "Fixed-effect descriptive"),
             k = nrow(x), n_total = sum(x$n), events = sum(x$events),
             log_odds = unname(fit$b[1]), odds_ratio = exp(unname(fit$b[1])),
             ci_lower = exp(fit$ci.lb), ci_upper = exp(fit$ci.ub), p_value = fit$pval,
             tau2 = ifelse(method == "REML", fit$tau2, NA_real_),
             I2 = ifelse(method == "REML", fit$I2, NA_real_), stringsAsFactors = FALSE)
}

available_adjusted <- adjusted[adjusted$adjustment_status == "Estimated", ]
adj_meta_input <- data.frame(
  cohort = available_adjusted$cohort, mechanism = available_adjusted$mechanism,
  n = available_adjusted$n_complete, events = available_adjusted$events_complete,
  log_odds = available_adjusted$adjusted_log_odds, ci_lower = available_adjusted$adjusted_ci_lower,
  ci_upper = available_adjusted$adjusted_ci_upper, p_value = available_adjusted$adjusted_p,
  stringsAsFactors = FALSE
)
unadj_available <- primary[match(available_adjusted$cohort, primary$cohort), ]
adjusted_meta <- rbind(
  cbind(analysis = "Unadjusted full cohort set", meta_summary(primary, "All mechanisms")),
  cbind(analysis = "Unadjusted covariate-available subset", meta_summary(unadj_available, "All mechanisms")),
  cbind(analysis = "Adjusted covariate-available subset", meta_summary(adj_meta_input, "All mechanisms"))
)
for (mech in unique(adj_meta_input$mechanism)) {
  z <- adj_meta_input[adj_meta_input$mechanism == mech, ]
  adjusted_meta <- rbind(adjusted_meta,
                         cbind(analysis = "Adjusted covariate-available subset", meta_summary(z, mech)))
}
rownames(adjusted_meta) <- NULL

# Explicit between-mechanism audit. The omnibus moderator test is paired with
# a magnitude gate because a non-significant P value alone is not evidence of
# equivalence when mechanism groups contain few cohorts.
mechanism_moderator_test <- function(x, analysis_label, mechanism_summaries,
                                     log_col = "log_odds", lower_col = "ci_lower",
                                     upper_col = "ci_upper") {
  x <- add_se(x, log_col, lower_col, upper_col)
  x$mechanism <- factor(x$mechanism)
  fit <- metafor::rma.uni(
    yi = yi, sei = sei, mods = ~ mechanism, data = x,
    method = "REML", test = "knha"
  )
  mechanism_or <- mechanism_summaries$odds_ratio
  data.frame(
    analysis = analysis_label,
    k = nrow(x),
    n_mechanisms = nlevels(x$mechanism),
    moderator_test = "REML meta-regression; Hartung-Knapp omnibus F test",
    moderator_statistic = unname(fit$QM),
    df1 = unname(fit$QMdf[1]),
    df2 = unname(fit$QMdf[2]),
    moderator_p = unname(fit$QMp),
    residual_I2 = unname(fit$I2),
    min_mechanism_or = min(mechanism_or),
    max_mechanism_or = max(mechanism_or),
    max_min_or_ratio = max(mechanism_or) / min(mechanism_or),
    stringsAsFactors = FALSE
  )
}

full_mechanism_summaries <- mechanism_meta[mechanism_meta$group != "All biologic mechanisms", ]
adjusted_mechanism_summaries <- adjusted_meta[
  adjusted_meta$analysis == "Adjusted covariate-available subset" &
    adjusted_meta$group != "All mechanisms",
]
mechanism_moderator <- rbind(
  mechanism_moderator_test(
    primary, "Unadjusted full cohort set", full_mechanism_summaries
  ),
  mechanism_moderator_test(
    adj_meta_input, "Adjusted covariate-available subset", adjusted_mechanism_summaries
  )
)

# Stratified sensitivity meta-analysis ---------------------------------------
primary$endpoint_type <- ifelse(
  primary$cohort %in% c("GSE16879 IFX", "GSE73661 IFX", "GSE73661 VDZ trial", "GSE73661 VDZ observational"),
  "Endoscopic/mucosal-healing response",
  ifelse(primary$cohort == "GSE206285 ustekinumab", "Clinical remission", "Clinical response")
)
primary$evaluation_time <- ifelse(
  primary$cohort %in% c("GSE16879 IFX", "GSE73661 IFX", "GSE92415 golimumab", "GSE73661 VDZ trial"),
  "Week 4-6", ifelse(primary$cohort %in% c("GSE23597 IFX", "GSE206285 ustekinumab"), "Week 8", "Week 12")
)
primary$biopsy_site_group <- ifelse(primary$cohort %in% distal_cohorts,
                                    "Protocol-fixed distal/sigmoid (15-20 cm)",
                                    "Colonic mucosa; exact site unspecified")

stratified_rows <- list()
for (dimension in c("endpoint_type", "evaluation_time", "mechanism", "biopsy_site_group")) {
  for (stratum in unique(primary[[dimension]])) {
    x <- primary[primary[[dimension]] == stratum, ]
    z <- meta_summary(x, stratum)
    z$dimension <- dimension
    z$stratum <- stratum
    stratified_rows[[paste(dimension, stratum)]] <- z
  }
}
stratified <- do.call(rbind, stratified_rows)
stratified <- stratified[, c("dimension", "stratum", "model", "k", "n_total", "events",
                             "log_odds", "odds_ratio", "ci_lower", "ci_upper", "p_value", "tau2", "I2")]
rownames(stratified) <- NULL

# Hierarchical, mutually exclusive sample-flow reasons -----------------------
cohort_rules <- list(
  "GSE16879 IFX" = list(accession = "GSE16879", uc = function(x) tolower(x$disease) == "uc",
                         baseline = function(x) x$visit == "Week 0", branch = function(x) rep(TRUE, nrow(x))),
  "GSE73661 IFX" = list(accession = "GSE73661", uc = function(x) x$disease == "ulcerative colitis",
                         baseline = function(x) x$visit == "W0", branch = function(x) x$treatment == "IFX"),
  "GSE73661 VDZ trial" = list(accession = "GSE73661", uc = function(x) x$disease == "ulcerative colitis",
                               baseline = function(x) x$visit == "W0",
                               branch = function(x) x$treatment %in% c("vdz_vdz4w", "vdz_vdz8w", "vdz_plac")),
  "GSE73661 VDZ observational" = list(accession = "GSE73661", uc = function(x) x$disease == "ulcerative colitis",
                                       baseline = function(x) x$visit == "W0", branch = function(x) x$treatment == "vdz4w"),
  "GSE23597 IFX" = list(accession = "GSE23597", uc = function(x) rep(TRUE, nrow(x)),
                         baseline = function(x) x$visit == "W0", branch = function(x) x$treatment == "infliximab"),
  "GSE92415 golimumab" = list(accession = "GSE92415", uc = function(x) grepl("Ulcerative", x$disease),
                               baseline = function(x) x$visit == "Week 0", branch = function(x) x$treatment == "golimumab"),
  "GSE206285 ustekinumab" = list(accession = "GSE206285", uc = function(x) x$disease == "ulcerative colitis",
                                  baseline = function(x) x$visit == "WEEK I-0", branch = function(x) x$treatment == "ustekinumab")
)

exclusion_rows <- list()
for (cohort_name in names(cohort_rules)) {
  rule <- cohort_rules[[cohort_name]]
  x <- read_scores(rule$accession)
  x$subject <- as.character(x$subject)
  is_uc <- rule$uc(x)
  is_base <- rule$baseline(x)
  is_branch <- rule$branch(x)
  known_outcome <- x$response %in% c("Yes", "No")
  reason <- ifelse(!is_uc, "Non-UC/control sample",
                   ifelse(!is_base, "Non-baseline/follow-up sample",
                          ifelse(!is_branch, "Other treatment arm or cohort branch",
                                 ifelse(!known_outcome, "Missing or unusable outcome", "Eligible before subject collapse"))))
  eligible_idx <- which(reason == "Eligible before subject collapse")
  if (length(eligible_idx)) {
    duplicate_idx <- eligible_idx[duplicated(x$subject[eligible_idx])]
    reason[eligible_idx] <- "Included independent participant"
    if (length(duplicate_idx)) reason[duplicate_idx] <- "Duplicate baseline biopsy collapsed within participant"
  }
  counts <- as.data.frame(table(reason), stringsAsFactors = FALSE)
  names(counts) <- c("reason", "n_samples")
  counts$cohort <- cohort_name
  counts$accession <- rule$accession
  counts$raw_accession_samples <- nrow(x)
  counts$final_participants <- sum(reason == "Included independent participant")
  counts$share_of_accession_samples <- counts$n_samples / nrow(x)
  exclusion_rows[[cohort_name]] <- counts
}
exclusion_audit <- do.call(rbind, exclusion_rows)
exclusion_audit <- exclusion_audit[, c("accession", "cohort", "reason", "n_samples",
                                       "share_of_accession_samples", "raw_accession_samples", "final_participants")]
rownames(exclusion_audit) <- NULL

# Stop gates -----------------------------------------------------------------
adjusted_overall <- adjusted_meta[adjusted_meta$analysis == "Adjusted covariate-available subset" &
                                    adjusted_meta$group == "All mechanisms", ]
full_overall <- adjusted_meta[adjusted_meta$analysis == "Unadjusted full cohort set", ]
reverse_supported <- any(available_adjusted$adjusted_or > 1 & available_adjusted$adjusted_ci_lower > 1, na.rm = TRUE)

adjusted_loo_rows <- lapply(seq_len(nrow(adj_meta_input)), function(i) {
  z <- meta_summary(adj_meta_input[-i, ], paste0("omit ", adj_meta_input$cohort[i]))
  z$omitted <- adj_meta_input$cohort[i]
  z
})
adjusted_loo <- do.call(rbind, adjusted_loo_rows)

mechanism_heterogeneity_pass <- all(mechanism_moderator$moderator_p >= 0.10) &&
  all(mechanism_moderator$max_min_or_ratio < 2)

gate <- data.frame(
  gate = c(
    "Adjusted pooled direction",
    "No supported reverse adjusted cohort",
    "Adjusted heterogeneity below substantive threshold",
    "Unadjusted heterogeneity below substantive threshold",
    "No substantial between-mechanism heterogeneity",
    "Mechanism-level directions remain adverse",
    "Full leave-one-out direction and significance",
    "Adjusted-subset leave-one-out direction"
  ),
  criterion = c(
    "Adjusted pooled OR < 1",
    "No adjusted cohort has OR > 1 with 95% CI wholly > 1",
    "Adjusted pooled I2 < 50%",
    "Unadjusted pooled I2 < 50%",
    "Full and adjusted moderator P >= 0.10 and max/min mechanism OR ratio < 2",
    "All mechanism-level ORs < 1",
    "Every full leave-one-out OR < 1 and upper 95% CI < 1",
    "Every adjusted-subset leave-one-out OR < 1"
  ),
  observed = c(
    sprintf("OR=%.4f; 95%% CI %.4f-%.4f", adjusted_overall$odds_ratio, adjusted_overall$ci_lower, adjusted_overall$ci_upper),
    paste(available_adjusted$cohort[available_adjusted$adjusted_or > 1 & available_adjusted$adjusted_ci_lower > 1], collapse = ";"),
    sprintf("I2=%.2f%%", adjusted_overall$I2),
    sprintf("I2=%.2f%%", full_overall$I2),
    paste(
      mechanism_moderator$analysis,
      sprintf("P=%.4f ratio=%.3f", mechanism_moderator$moderator_p,
              mechanism_moderator$max_min_or_ratio),
      collapse = "; "
    ),
    paste(mechanism_meta$group[mechanism_meta$group != "All biologic mechanisms"],
          sprintf("OR=%.3f", mechanism_meta$odds_ratio[mechanism_meta$group != "All biologic mechanisms"]), collapse = "; "),
    paste(leave_one_out$omitted, sprintf("OR=%.3f UCL=%.3f", leave_one_out$odds_ratio, leave_one_out$ci_upper), collapse = "; "),
    paste(adjusted_loo$omitted, sprintf("OR=%.3f", adjusted_loo$odds_ratio), collapse = "; ")
  ),
  passed = c(
    adjusted_overall$odds_ratio < 1,
    !reverse_supported,
    is.na(adjusted_overall$I2) || adjusted_overall$I2 < 50,
    is.na(full_overall$I2) || full_overall$I2 < 50,
    mechanism_heterogeneity_pass,
    all(mechanism_meta$odds_ratio[mechanism_meta$group != "All biologic mechanisms"] < 1),
    all(leave_one_out$odds_ratio < 1 & leave_one_out$ci_upper < 1),
    all(adjusted_loo$odds_ratio < 1)
  ),
  stop_if_failed = TRUE,
  stringsAsFactors = FALSE
)
gate$decision <- ifelse(gate$passed, "Continue", "STOP AND RE-PLAN")

utils::write.csv(clinical, file.path(source_dir, "covariate_augmented_baseline_patient_data.csv"), row.names = FALSE)
utils::write.csv(availability, file.path(tables_dir, "table_covariate_availability_matrix.csv"), row.names = FALSE)
utils::write.csv(adjusted, file.path(tables_dir, "table_adjusted_baseline_associations.csv"), row.names = FALSE)
utils::write.csv(adjusted_meta, file.path(tables_dir, "table_adjusted_meta_analysis.csv"), row.names = FALSE)
utils::write.csv(mechanism_moderator, file.path(tables_dir, "table_mechanism_moderator_test.csv"), row.names = FALSE)
utils::write.csv(stratified, file.path(tables_dir, "table_stratified_sensitivity_meta.csv"), row.names = FALSE)
utils::write.csv(exclusion_audit, file.path(tables_dir, "table_sample_exclusion_audit.csv"), row.names = FALSE)
utils::write.csv(adjusted_loo, file.path(tables_dir, "table_adjusted_subset_leave_one_out.csv"), row.names = FALSE)
utils::write.csv(gate, file.path(tables_dir, "table_confounding_heterogeneity_gate.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "18_confounding_heterogeneity_audit_sessionInfo.txt"))

print(adjusted[, c("cohort", "n_complete", "unadjusted_or", "adjusted_or", "adjusted_ci_lower", "adjusted_ci_upper", "adjusted_p")])
print(adjusted_meta)
print(stratified)
print(gate)
if (!all(gate$passed)) {
  stop("CONFOUNDING/HETEROGENEITY AUDIT STOP GATE: ", paste(gate$gate[!gate$passed], collapse = "; "))
}
message("CONFOUNDING/HETEROGENEITY AUDIT PASSED")
