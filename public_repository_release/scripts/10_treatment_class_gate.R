#!/usr/bin/env Rscript

# Treatment-class branch gate. Ustekinumab is evaluated as an RCT interaction;
# vedolizumab is evaluated in two separately defined treated cohorts. A non-TNF
# effect in the same direction and at least half the corresponding anti-TNF
# log-odds magnitude stops the workflow for scientific re-planning.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
raw_dir <- file.path(project_root, "data", "raw")
derived_dir <- file.path(project_root, "data", "derived", "extended")
tables_dir <- file.path(project_root, "results", "tables", "extended")
logs_dir <- file.path(project_root, "logs")

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(logistf)
  library(metafor)
})

read_scores <- function(accession) {
  readRDS(file.path(derived_dir, paste0(tolower(accession), "_all_scores.rds")))
}

extract_firth <- function(data, formula, term, cohort, estimand) {
  fit <- logistf::logistf(formula, data = data)
  data.frame(
    cohort = cohort, estimand = estimand, n = nrow(data), events = sum(data$response_binary == 1L),
    log_odds = unname(fit$coefficients[term]), odds_ratio = exp(unname(fit$coefficients[term])),
    ci_lower = exp(unname(fit$ci.lower[term])), ci_upper = exp(unname(fit$ci.upper[term])),
    p_value = unname(fit$prob[term]), stringsAsFactors = FALSE
  )
}

# Ustekinumab RCT: clinical remission is complete for all 550 UC participants.
ust <- read_scores("GSE206285")
ust <- ust[ust$eligible & ust$response %in% c("Yes", "No"), ]
ust$response_binary <- as.integer(ust$response == "Yes")
ust$treated <- as.integer(ust$treatment == "ustekinumab")
ust_interaction <- extract_firth(
  ust, response_binary ~ treated * Inflammation_z, "treated:Inflammation_z",
  "GSE206285 UST vs placebo", "Treatment-by-score interaction for Week 8 clinical remission"
)
ust_arms <- rbind(
  extract_firth(ust[ust$treated == 1L, ], response_binary ~ Inflammation_z, "Inflammation_z",
                "GSE206285 ustekinumab", "Remission OR per 1-SD baseline score"),
  extract_firth(ust[ust$treated == 0L, ], response_binary ~ Inflammation_z, "Inflammation_z",
                "GSE206285 placebo", "Remission OR per 1-SD baseline score")
)

# Secondary UST endpoint: mucosal healing, read directly from the public matrix.
g206 <- GEOquery::getGEO(filename = file.path(raw_dir, "GSE206285_series_matrix.txt.gz"), getGPL = FALSE)
if (is.list(g206)) g206 <- g206[[1]]
p206 <- Biobase::pData(g206)
mucosal <- p206[["mucosal healing at week 8:ch1"]]
names(mucosal) <- rownames(p206)
ust$mucosal_response <- mucosal[ust$gsm]
ust_mh <- ust[ust$mucosal_response %in% c("Y", "N"), ]
ust_mh$response_binary <- as.integer(ust_mh$mucosal_response == "Y")
ust_mh_interaction <- extract_firth(
  ust_mh, response_binary ~ treated * Inflammation_z, "treated:Inflammation_z",
  "GSE206285 UST vs placebo", "Treatment-by-score interaction for Week 8 mucosal healing"
)

# Vedolizumab cohorts remain separate because they have different protocols and
# outcome times. Only subjects with a title-derived response are included.
vdz <- read_scores("GSE73661")
vdz_trial_codes <- c("vdz_vdz4w", "vdz_vdz8w", "vdz_plac")
vdz_trial <- vdz[vdz$visit == "W0" & vdz$treatment %in% vdz_trial_codes &
                   vdz$response %in% c("Yes", "No"), ]
vdz_obs <- vdz[vdz$visit == "W0" & vdz$treatment == "vdz4w" &
                 vdz$response %in% c("Yes", "No"), ]
vdz_trial$response_binary <- as.integer(vdz_trial$response == "Yes")
vdz_obs$response_binary <- as.integer(vdz_obs$response == "Yes")
vdz_rows <- rbind(
  extract_firth(vdz_trial, response_binary ~ Inflammation_z, "Inflammation_z",
                "GSE73661 VDZ trial", "Week 6 response OR per 1-SD baseline score"),
  extract_firth(vdz_obs, response_binary ~ Inflammation_z, "Inflammation_z",
                "GSE73661 VDZ observational", "Week 12 response OR per 1-SD baseline score")
)
vdz_rows$sei <- (log(vdz_rows$ci_upper) - log(vdz_rows$ci_lower)) / (2 * 1.96)
vdz_fixed <- metafor::rma.uni(yi = vdz_rows$log_odds, sei = vdz_rows$sei, method = "FE")
vdz_summary <- data.frame(
  cohort = "GSE73661 VDZ fixed-effect descriptive synthesis",
  estimand = "Response OR per 1-SD baseline score", n = sum(vdz_rows$n), events = sum(vdz_rows$events),
  log_odds = unname(vdz_fixed$b[1]), odds_ratio = exp(unname(vdz_fixed$b[1])),
  ci_lower = exp(vdz_fixed$ci.lb), ci_upper = exp(vdz_fixed$ci.ub), p_value = vdz_fixed$pval
)

utils::write.csv(rbind(ust_interaction, ust_mh_interaction),
                 file.path(tables_dir, "table_ustekinumab_rct_interactions.csv"), row.names = FALSE)
utils::write.csv(ust_arms, file.path(tables_dir, "table_ustekinumab_arm_associations.csv"), row.names = FALSE)
utils::write.csv(vdz_rows, file.path(tables_dir, "table_vedolizumab_associations.csv"), row.names = FALSE)
utils::write.csv(vdz_summary, file.path(tables_dir, "table_vedolizumab_association_synthesis.csv"), row.names = FALSE)

anti_rct <- utils::read.csv(file.path(tables_dir, "table_antitnf_rct_interaction_synthesis.csv"))
anti_assoc <- utils::read.csv(file.path(tables_dir, "table_antitnf_association_synthesis.csv"))
ust_threshold <- 0.5 * abs(anti_rct$log_odds[1])
vdz_threshold <- 0.5 * abs(anti_assoc$log_odds[1])
ust_comparable <- ust_interaction$log_odds < 0 && abs(ust_interaction$log_odds) >= ust_threshold
vdz_comparable <- vdz_summary$log_odds < 0 && abs(vdz_summary$log_odds) >= vdz_threshold

decision <- data.frame(
  comparison = c("UST RCT interaction", "VDZ treated-arm association"),
  non_tnf_log_odds = c(ust_interaction$log_odds, vdz_summary$log_odds),
  anti_tnf_reference_log_odds = c(anti_rct$log_odds[1], anti_assoc$log_odds[1]),
  locked_half_magnitude_threshold = c(ust_threshold, vdz_threshold),
  same_direction_and_comparable = c(ust_comparable, vdz_comparable),
  stringsAsFactors = FALSE
)
utils::write.csv(decision, file.path(tables_dir, "table_treatment_class_branch_decision.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "10_treatment_class_gate_sessionInfo.txt"))

print(rbind(ust_interaction, ust_mh_interaction))
print(ust_arms)
print(vdz_rows)
print(vdz_summary)
print(decision)
if (any(decision$same_direction_and_comparable)) {
  triggered <- paste(decision$comparison[decision$same_direction_and_comparable], collapse = ", ")
  stop("TREATMENT-CLASS BRANCH GATE TRIGGERED: ", triggered,
       ". Stop before longitudinal and single-cell analysis; re-plan the paper as a potentially general biologic-resistance state.")
}
message("TREATMENT-CLASS GATE PASSED FOR THE EXPLORATORY ANTI-TNF-SPECIFIC BRANCH")
