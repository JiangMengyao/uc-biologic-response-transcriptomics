#!/usr/bin/env Rscript

# Locked anti-TNF direction analysis. Results are written before the direction
# gates are evaluated so a failed gate is auditable. No downstream class-
# specificity analysis should run if this script exits non-zero.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
derived_dir <- file.path(project_root, "data", "derived", "extended")
tables_dir <- file.path(project_root, "results", "tables", "extended")
logs_dir <- file.path(project_root, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(logistf)
  library(metafor)
})

read_scores <- function(accession) {
  readRDS(file.path(derived_dir, paste0(tolower(accession), "_all_scores.rds")))
}

collapse_subject_time <- function(data) {
  score_cols <- grep("(_z$|^Inflammation$|^TNFaNFKB$|^IL6JAKSTAT3$|^OxidativePhos_neg$|^Random_neg$)",
                     names(data), value = TRUE)
  keys <- c("accession", "subject", "disease", "treatment", "dose", "visit", "response", "endpoint")
  aggregate(data[, score_cols, drop = FALSE], data[keys], mean, na.rm = TRUE)
}

extract_firth <- function(data, formula, term, cohort, estimand) {
  fit <- logistf::logistf(formula, data = data)
  data.frame(
    cohort = cohort, estimand = estimand, n = nrow(data), events = sum(data$response_binary == 1L),
    term = term, log_odds = unname(fit$coefficients[term]),
    odds_ratio = exp(unname(fit$coefficients[term])),
    ci_lower = exp(unname(fit$ci.lower[term])), ci_upper = exp(unname(fit$ci.upper[term])),
    p_value = unname(fit$prob[term]), stringsAsFactors = FALSE
  )
}

prepare_binary <- function(data) {
  data <- data[data$response %in% c("Yes", "No"), ]
  data$response_binary <- as.integer(data$response == "Yes")
  data
}

# Baseline cohorts.
d16879 <- read_scores("GSE16879")
d16879 <- prepare_binary(d16879[d16879$visit == "Week 0" & tolower(d16879$disease) == "uc", ])

d73661 <- read_scores("GSE73661")
d73661_ifx <- prepare_binary(d73661[d73661$visit == "W0" & d73661$treatment == "IFX", ])

d23597 <- read_scores("GSE23597")
d23597 <- collapse_subject_time(d23597[d23597$visit == "W0" & d23597$response %in% c("Yes", "No"), ])
d23597 <- prepare_binary(d23597)
d23597$treated <- as.integer(d23597$treatment == "infliximab")

d92415 <- read_scores("GSE92415")
d92415 <- prepare_binary(d92415[d92415$visit == "Week 0" & grepl("Ulcerative", d92415$disease), ])
d92415$treated <- as.integer(d92415$treatment == "golimumab")

association_rows <- rbind(
  extract_firth(d16879, response_binary ~ Inflammation_z, "Inflammation_z", "GSE16879 IFX",
                "Response OR per 1-SD baseline score"),
  extract_firth(d73661_ifx, response_binary ~ Inflammation_z, "Inflammation_z", "GSE73661 IFX",
                "Response OR per 1-SD baseline score"),
  extract_firth(d23597[d23597$treated == 1L, ], response_binary ~ Inflammation_z, "Inflammation_z", "GSE23597 IFX",
                "Response OR per 1-SD baseline score"),
  extract_firth(d92415[d92415$treated == 1L, ], response_binary ~ Inflammation_z, "Inflammation_z", "GSE92415 golimumab",
                "Response OR per 1-SD baseline score")
)
utils::write.csv(association_rows, file.path(tables_dir, "table_antitnf_treated_arm_associations.csv"), row.names = FALSE)

rct_rows <- rbind(
  extract_firth(d92415, response_binary ~ treated * Inflammation_z, "treated:Inflammation_z", "GSE92415",
                "Treatment-by-score interaction OR"),
  extract_firth(d23597, response_binary ~ treated * Inflammation_z, "treated:Inflammation_z", "GSE23597",
                "Treatment-by-score interaction OR")
)
utils::write.csv(rct_rows, file.path(tables_dir, "table_antitnf_rct_interactions.csv"), row.names = FALSE)

# Descriptive inverse-variance synthesis. Profile intervals are converted to an
# approximate standard error; k=2 precludes a reliable random-effects claim.
rct_rows$sei <- (log(rct_rows$ci_upper) - log(rct_rows$ci_lower)) / (2 * 1.96)
fixed_rct <- metafor::rma.uni(yi = rct_rows$log_odds, sei = rct_rows$sei, method = "FE")
rct_summary <- data.frame(
  cohort = "Fixed-effect descriptive synthesis", estimand = "Treatment-by-score interaction OR",
  n = sum(rct_rows$n), events = sum(rct_rows$events), term = "treated:Inflammation_z",
  log_odds = unname(fixed_rct$b[1]), odds_ratio = exp(unname(fixed_rct$b[1])),
  ci_lower = exp(fixed_rct$ci.lb), ci_upper = exp(fixed_rct$ci.ub), p_value = fixed_rct$pval
)
utils::write.csv(rct_summary, file.path(tables_dir, "table_antitnf_rct_interaction_synthesis.csv"), row.names = FALSE)

# Treated-arm association synthesis: show cohort effects as primary; REML-HK is
# supportive because only four cohorts are available and endpoints differ.
association_rows$sei <- (log(association_rows$ci_upper) - log(association_rows$ci_lower)) / (2 * 1.96)
random_assoc <- metafor::rma.uni(yi = association_rows$log_odds, sei = association_rows$sei,
                                 method = "REML", test = "knha")
association_summary <- data.frame(
  model = "REML Hartung-Knapp", k = random_assoc$k, log_odds = unname(random_assoc$b[1]),
  odds_ratio = exp(unname(random_assoc$b[1])), ci_lower = exp(random_assoc$ci.lb),
  ci_upper = exp(random_assoc$ci.ub), p_value = random_assoc$pval,
  tau2 = random_assoc$tau2, I2 = random_assoc$I2
)
utils::write.csv(association_summary, file.path(tables_dir, "table_antitnf_association_synthesis.csv"), row.names = FALSE)

gate_table <- rbind(
  transform(association_rows[, c("cohort", "odds_ratio", "ci_lower", "ci_upper", "p_value")],
            gate = "treated-arm response OR < 1", passed = odds_ratio < 1),
  transform(rct_rows[, c("cohort", "odds_ratio", "ci_lower", "ci_upper", "p_value")],
            gate = "RCT interaction OR < 1", passed = odds_ratio < 1)
)
utils::write.csv(gate_table, file.path(tables_dir, "table_antitnf_direction_gates.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "09_antitnf_association_and_rct_sessionInfo.txt"))

print(association_rows)
print(rct_rows)
print(rct_summary)
if (!all(gate_table$passed)) {
  failed <- gate_table[!gate_table$passed, ]
  stop("DIRECTION GATE FAILED: ", paste(failed$cohort, collapse = ", "))
}
message("ALL LOCKED ANTI-TNF DIRECTION GATES PASSED")
