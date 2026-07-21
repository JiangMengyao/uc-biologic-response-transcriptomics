#!/usr/bin/env Rscript

# Build a subject-aware registry for all bulk cohorts and enforce the audited
# sample-count gates before any model is fit. A failed gate stops the workflow.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
raw_dir <- file.path(project_root, "data", "raw")
derived_dir <- file.path(project_root, "data", "derived", "extended")
tables_dir <- file.path(project_root, "results", "tables", "extended")
logs_dir <- file.path(project_root, "logs")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

read_eset <- function(accession) {
  path <- file.path(raw_dir, paste0(accession, "_series_matrix.txt.gz"))
  if (!file.exists(path)) stop("Missing matrix: ", path)
  object <- GEOquery::getGEO(filename = path, getGPL = FALSE)
  if (is.list(object)) object <- object[[1]]
  object
}

gate_equal <- function(observed, expected, label) {
  if (!identical(as.integer(observed), as.integer(expected))) {
    stop(sprintf("GATE FAILED [%s]: observed=%s expected=%s",
                 label, paste(observed, collapse = ","), paste(expected, collapse = ",")))
  }
  message("GATE PASS [", label, "]: ", paste(observed, collapse = ","))
}

registry_rows <- list()
append_registry <- function(data) {
  registry_rows[[length(registry_rows) + 1L]] <<- data
}

# GSE92415: randomized golimumab/placebo induction study.
g92415 <- read_eset("GSE92415")
p <- Biobase::pData(g92415)
d92415 <- data.frame(
  accession = "GSE92415", gsm = rownames(p), subject = p[["subject:ch1"]],
  disease = p[["disease:ch1"]], treatment = tolower(p[["treatment:ch1"]]),
  dose = NA_character_, visit = p[["visit:ch1"]],
  response = p[["wk6response:ch1"]], endpoint = "Week 6 clinical response",
  tissue = "colonic mucosal biopsy", platform = Biobase::annotation(g92415),
  stringsAsFactors = FALSE
)
d92415$eligible <- grepl("Ulcerative", d92415$disease) & d92415$response %in% c("Yes", "No")
b92415 <- d92415[d92415$eligible & d92415$visit == "Week 0", ]
w92415 <- d92415[d92415$eligible & d92415$visit == "Week 6", ]
gate_equal(nrow(b92415), 87L, "GSE92415 baseline n")
gate_equal(c(sum(b92415$treatment == "golimumab"), sum(b92415$treatment == "placebo")),
           c(59L, 28L), "GSE92415 randomized arms")
gate_equal(nrow(w92415), 75L, "GSE92415 Week 6 n")
gate_equal(length(intersect(b92415$subject, w92415$subject)), 65L, "GSE92415 paired Week 0-6")
append_registry(d92415)

# GSE16879: paired UC biopsies before/after first infliximab treatment.
g16879 <- read_eset("GSE16879")
p <- Biobase::pData(g16879)
d16879 <- data.frame(
  accession = "GSE16879", gsm = rownames(p),
  subject = sub("_(before|after)T$", "", p$title),
  disease = p[["disease:ch1"]], treatment = "infliximab", dose = NA_character_,
  visit = ifelse(grepl("Before", p[["before or after first infliximab treatment:ch1"]]),
                 "Week 0", "Week 4-6"),
  response = p[["response to infliximab:ch1"]], endpoint = "Week 4-6 mucosal healing",
  tissue = p[["tissue:ch1"]], platform = Biobase::annotation(g16879),
  stringsAsFactors = FALSE
)
d16879$eligible <- tolower(d16879$disease) == "uc" & d16879$response %in% c("Yes", "No")
b16879 <- d16879[d16879$eligible & d16879$visit == "Week 0", ]
a16879 <- d16879[d16879$eligible & d16879$visit == "Week 4-6", ]
gate_equal(c(nrow(b16879), sum(b16879$response == "Yes"), sum(b16879$response == "No")),
           c(24L, 8L, 16L), "GSE16879 baseline/R/NR")
gate_equal(length(intersect(b16879$subject, a16879$subject)), 24L, "GSE16879 paired")
append_registry(d16879)

# GSE23597: ACT1 randomized infliximab/placebo substudy.
g23597 <- read_eset("GSE23597")
p <- Biobase::pData(g23597)
d23597 <- data.frame(
  accession = "GSE23597", gsm = rownames(p), subject = p[["subject:ch1"]],
  disease = "ulcerative colitis", treatment = ifelse(p[["dose:ch1"]] == "placebo", "placebo", "infliximab"),
  dose = p[["dose:ch1"]], visit = p[["time:ch1"]], response = p[["wk8 response:ch1"]],
  endpoint = "Week 8 clinical response", tissue = "colonic mucosal biopsy",
  platform = Biobase::annotation(g23597), stringsAsFactors = FALSE
)
d23597$eligible <- d23597$response %in% c("Yes", "No")
b23597 <- d23597[d23597$eligible & d23597$visit == "W0", ]
subject_base_23597 <- unique(b23597[, c("subject", "treatment", "dose", "response")])
gate_equal(nrow(subject_base_23597), 44L, "GSE23597 independent baseline subjects")
gate_equal(c(sum(subject_base_23597$treatment == "infliximab"), sum(subject_base_23597$treatment == "placebo")),
           c(31L, 13L), "GSE23597 randomized arms")
gate_equal(c(sum(subject_base_23597$treatment == "infliximab" & subject_base_23597$response == "Yes"),
             sum(subject_base_23597$treatment == "infliximab" & subject_base_23597$response == "No"),
             sum(subject_base_23597$treatment == "placebo" & subject_base_23597$response == "Yes"),
             sum(subject_base_23597$treatment == "placebo" & subject_base_23597$response == "No")),
           c(24L, 7L, 5L, 8L), "GSE23597 Week 8 response by arm")
gate_equal(length(intersect(unique(b23597$subject), unique(d23597$subject[d23597$visit == "W8"]))),
           33L, "GSE23597 paired Week 0-8")
gate_equal(length(intersect(unique(b23597$subject), unique(d23597$subject[d23597$visit == "W30"]))),
           29L, "GSE23597 paired Week 0-30")
append_registry(d23597)

# GSE73661: independent IFX and VDZ treatment series. Response is encoded in
# the post-treatment sample title, so derive it at subject level and keep the
# W6 randomized VDZ series separate from the W12 observational VDZ series.
g73661 <- read_eset("GSE73661")
p <- Biobase::pData(g73661)
response_from_title <- ifelse(grepl("UC NR", p$title), "No",
                              ifelse(grepl("UC R", p$title), "Yes", NA_character_))
d73661 <- data.frame(
  accession = "GSE73661", gsm = rownames(p), subject = p[["study individual number:ch1"]],
  disease = ifelse(grepl("UC", p$title), "ulcerative colitis", "control"),
  treatment = p[["induction therapy_maintenance therapy:ch1"]], dose = NA_character_,
  visit = p[["week (w):ch1"]], response = response_from_title,
  endpoint = NA_character_, tissue = p[["tissue:ch1"]],
  platform = Biobase::annotation(g73661), stringsAsFactors = FALSE
)
ifx_post <- d73661[d73661$treatment == "IFX" & d73661$visit == "W4_W6", c("subject", "response")]
ifx_response <- setNames(ifx_post$response, ifx_post$subject)
ifx_base_idx <- d73661$treatment == "IFX" & d73661$visit == "W0"
d73661$response[ifx_base_idx] <- ifx_response[d73661$subject[ifx_base_idx]]
d73661$endpoint[ifx_base_idx | (d73661$treatment == "IFX" & d73661$visit == "W4_W6")] <- "Week 4-6 response"

vdz_trial_codes <- c("vdz_vdz4w", "vdz_vdz8w", "vdz_plac")
vdz_trial_post <- d73661[d73661$treatment %in% vdz_trial_codes & d73661$visit == "W6",
                         c("subject", "response")]
vdz_trial_response <- setNames(vdz_trial_post$response, vdz_trial_post$subject)
vdz_trial_base_idx <- d73661$treatment %in% vdz_trial_codes & d73661$visit == "W0"
d73661$response[vdz_trial_base_idx] <- vdz_trial_response[d73661$subject[vdz_trial_base_idx]]
d73661$endpoint[vdz_trial_base_idx | (d73661$treatment %in% vdz_trial_codes & d73661$visit == "W6")] <- "Week 6 response"

vdz_obs_post <- d73661[d73661$treatment == "vdz4w" & d73661$visit == "W12", c("subject", "response")]
vdz_obs_response <- setNames(vdz_obs_post$response, vdz_obs_post$subject)
vdz_obs_base_idx <- d73661$treatment == "vdz4w" & d73661$visit == "W0"
d73661$response[vdz_obs_base_idx] <- vdz_obs_response[d73661$subject[vdz_obs_base_idx]]
d73661$endpoint[vdz_obs_base_idx | (d73661$treatment == "vdz4w" & d73661$visit == "W12")] <- "Week 12 response"

b73661_ifx <- d73661[ifx_base_idx, ]
gate_equal(c(nrow(b73661_ifx), sum(b73661_ifx$response == "Yes"), sum(b73661_ifx$response == "No")),
           c(23L, 8L, 15L), "GSE73661 IFX baseline/R/NR")
b73661_vdz_trial <- d73661[vdz_trial_base_idx & !is.na(d73661$response), ]
gate_equal(c(nrow(b73661_vdz_trial), sum(b73661_vdz_trial$response == "Yes"), sum(b73661_vdz_trial$response == "No")),
           c(27L, 6L, 21L), "GSE73661 VDZ trial baseline with Week 6 outcome")
b73661_vdz_obs <- d73661[vdz_obs_base_idx, ]
gate_equal(c(nrow(b73661_vdz_obs), sum(b73661_vdz_obs$response == "Yes"), sum(b73661_vdz_obs$response == "No")),
           c(13L, 3L, 10L), "GSE73661 VDZ observational baseline/R/NR")
d73661$eligible <- !is.na(d73661$response)
append_registry(d73661)

# GSE206285: UNIFI ustekinumab/placebo randomized trial.
g206285 <- read_eset("GSE206285")
p <- Biobase::pData(g206285)
d206285 <- data.frame(
  accession = "GSE206285", gsm = rownames(p), subject = p[["donor id:ch1"]],
  disease = p[["diagnosis:ch1"]], treatment = ifelse(grepl("Ustekinumab", p[["treatment:ch1"]]),
                                                       "ustekinumab", "placebo"),
  dose = p[["treatment:ch1"]], visit = p[["visit:ch1"]],
  response = p[["clinical remission at week 8:ch1"]], endpoint = "Week 8 clinical remission",
  tissue = p[["tissue:ch1"]], platform = Biobase::annotation(g206285),
  stringsAsFactors = FALSE
)
d206285$response[d206285$response == "Y"] <- "Yes"
d206285$response[d206285$response == "N"] <- "No"
d206285$eligible <- d206285$disease == "ulcerative colitis" & d206285$response %in% c("Yes", "No")
b206285 <- d206285[d206285$eligible, ]
gate_equal(nrow(b206285), 550L, "GSE206285 UC baseline with remission outcome")
gate_equal(c(sum(b206285$treatment == "ustekinumab"), sum(b206285$treatment == "placebo")),
           c(364L, 186L), "GSE206285 randomized arms")
gate_equal(c(sum(b206285$response == "Yes"), sum(b206285$response == "No")),
           c(63L, 487L), "GSE206285 Week 8 remission")
append_registry(d206285)

registry <- do.call(rbind, registry_rows)
registry$randomized_arm <- registry$treatment %in% c("golimumab", "infliximab", "ustekinumab", "placebo") &
  registry$accession %in% c("GSE92415", "GSE23597", "GSE206285")
utils::write.csv(registry, file.path(derived_dir, "sample_registry.csv"), row.names = FALSE)

eligible_registry <- registry[!is.na(registry$eligible) & registry$eligible, ]
summary_rows <- do.call(rbind, lapply(split(eligible_registry, eligible_registry$accession), function(x) {
  data.frame(accession = x$accession[1], samples = nrow(x), subjects = length(unique(x$subject)),
             response_yes = sum(x$response == "Yes", na.rm = TRUE),
             response_no = sum(x$response == "No", na.rm = TRUE))
}))
utils::write.csv(summary_rows, file.path(tables_dir, "table_cohort_registry_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "07_build_cohort_registry_sessionInfo.txt"))
message("ALL REGISTRY GATES PASSED")
