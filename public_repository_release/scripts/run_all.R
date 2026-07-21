#!/usr/bin/env Rscript

# Reproducible, checkpointed entry point for the current cross-mechanism paper.
# Each stage runs in an isolated process, writes a dedicated console log and is
# skipped only when its outputs are newer than its scripts and dependencies.

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1])) else NA_character_
project_root <- if (nzchar(Sys.getenv("PROJECT_ROOT"))) {
  normalizePath(Sys.getenv("PROJECT_ROOT"))
} else if (!is.na(script_path)) {
  dirname(dirname(script_path))
} else {
  normalizePath(getwd())
}
Sys.setenv(PROJECT_ROOT = project_root)

args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(flag) flag %in% args
arg_value <- function(prefix, default = NULL) {
  hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", prefix, "="), "", hit[1]) else default
}

stage_names <- c(
  "preflight", "gene_sets", "bulk_download", "registry", "bulk_scores",
  "baseline_models", "confounding_audit", "longitudinal", "bulk_figures", "supplementary_figures",
  "singlecell_download", "singlecell_analysis", "singlecell_figure", "supplementary_tables", "validate"
)
if (has_flag("--list-stages")) {
  writeLines(stage_names)
  quit(status = 0L)
}

from_stage <- arg_value("--from", stage_names[1])
to_stage <- arg_value("--to", tail(stage_names, 1))
if (!from_stage %in% stage_names) stop("Unknown --from stage: ", from_stage)
if (!to_stage %in% stage_names) stop("Unknown --to stage: ", to_stage)
from_i <- match(from_stage, stage_names); to_i <- match(to_stage, stage_names)
if (from_i > to_i) stop("--from must precede --to")
force <- has_flag("--force")
offline <- has_flag("--offline")
skip_singlecell <- has_flag("--skip-single-cell")

p <- function(...) file.path(project_root, ...)
bulk_files <- p("data", "raw", paste0(c("GSE92415", "GSE16879", "GSE23597", "GSE73661", "GSE206285"),
                                        "_series_matrix.txt.gz"))
bulk_files <- c(bulk_files, p("data", "raw", "GSE206285_additional_sample_metadata.tsv.gz"))
singlecell_files <- p("data", "raw", "GSE282122",
                      c("paired_sample_list.csv", "UMAP_combined_objects.txt.gz", "myeloid_final.h5ad",
                        "epicolonic_final.h5ad", "fibperi_final.h5ad"))
figure_bundle <- function(stem, number) {
  c(p("results", "figures", "main", paste0(stem, c(".svg", ".pdf", ".tiff"))),
    p("results", "figures", "previews", paste0(stem, ".png")),
    p("results", "source_data", paste0("Figure", number, "_source_data.csv")))
}
supplement_bundle <- function(stem, source_name) {
  c(p("results", "figures", "supplementary", paste0(stem, c(".svg", ".pdf", ".tiff"))),
    p("results", "figures", "previews", paste0(stem, ".png")),
    p("results", "source_data", source_name))
}
bulk_figure_outputs <- unlist(Map(figure_bundle,
  c("Figure1_evidence_architecture", "Figure2_cross_mechanism_prognostic_state",
    "Figure3_prognostic_vs_predictive_RCT", "Figure4_longitudinal_overcoming"), 1:4), use.names = FALSE)
singlecell_analysis_outputs <- c(
  p("results", "source_data", "singlecell_sample_state_pseudobulk.csv"),
  p("results", "source_data", "singlecell_sample_compartment_pseudobulk.csv"),
  p("results", "source_data", "singlecell_compartment_pseudobulk.csv"),
  p("results", "source_data", "singlecell_donor_compartment_baseline.csv"),
  p("results", "source_data", "singlecell_compartment_paired_donors.csv"),
  p("results", "source_data", "singlecell_longitudinal_site_pairs.csv"),
  p("results", "source_data", "singlecell_compartment_key_gene_dotplot.csv"),
  p("results", "tables", "common_state", c(
    "table_singlecell_compartment_localization.csv",
    "table_singlecell_compartment_outcome_effects.csv",
    "table_singlecell_compartment_outcome_sensitivity.csv",
    "table_singlecell_compartment_longitudinal_effects.csv",
    "table_singlecell_longitudinal_pairing_audit.csv",
    "table_singlecell_state_outcome_effects.csv",
    "table_singlecell_state_longitudinal_effects.csv",
    "table_singlecell_dataset_level_audit.csv")),
  p("data", "derived", "singlecell", "cross_compartment_scaling.rds")
)
supplementary_figure_outputs <- unlist(Map(
  supplement_bundle,
  c("FigureS1_cohort_inclusion_flow", "FigureS2_baseline_sensitivity",
    "FigureS3_longitudinal_sensitivity", "FigureS4_high_state_rct_audit"),
  paste0("FigureS", 1:4, "_source_data.csv")
), use.names = FALSE)
supplementary_table_outputs <- c(
  p("results", "supplementary_tables", c(
    "TableS1_cohort_registry.csv",
    "TableS2_baseline_associations.csv",
    "TableS3_RCT_prognostic_predictive.csv",
    "TableS4_longitudinal_sensitivities.csv",
    "TableS5_singlecell_baseline_sensitivity.csv",
    "TableS6_singlecell_state_baseline.csv",
    "TableS7_singlecell_longitudinal.csv",
    "supplementary_tables_manifest.csv"
  ))
)

stage <- list(
  preflight = list(type = "R", scripts = "00_preflight.R", deps = character(),
                   outputs = p("environment", c("package_versions.csv", "system_info.csv", "sessionInfo.txt"))),
  gene_sets = list(type = "R", scripts = "00_define_gene_sets.R", deps = "preflight",
                   outputs = c(p("config", "predefined_gene_sets.csv"), p("data", "derived", "predefined_gene_sets.rds"))),
  bulk_download = list(type = "R", scripts = "01_download_data.R", deps = "preflight", outputs = bulk_files, download = TRUE),
  registry = list(type = "R", scripts = "07_build_cohort_registry.R", deps = "bulk_download",
                  outputs = p("data", "derived", "extended", "sample_registry.csv")),
  bulk_scores = list(type = "R", scripts = "08_prepare_extended_cohorts.R", deps = c("registry", "gene_sets"),
                     outputs = p("data", "derived", "extended",
                                 paste0(tolower(c("GSE92415", "GSE16879", "GSE23597", "GSE73661", "GSE206285")),
                                        "_all_scores.rds"))),
  baseline_models = list(type = "R", scripts = "11_common_prognostic_state.R", deps = "bulk_scores",
                         outputs = p("results", "tables", "common_state", c(
                           "table_mechanism_meta_analysis.csv", "table_cross_mechanism_gate.csv",
                           "table_primary_cross_mechanism_associations.csv",
                           "table_rct_treatment_score_interactions.csv"))),
  confounding_audit = list(type = "R", scripts = "18_confounding_heterogeneity_audit.R",
                           deps = c("baseline_models", "bulk_scores"),
                           outputs = p("results", "tables", "common_state", c(
                             "table_covariate_availability_matrix.csv",
                             "table_adjusted_baseline_associations.csv",
                             "table_adjusted_meta_analysis.csv",
                             "table_mechanism_moderator_test.csv",
                             "table_stratified_sensitivity_meta.csv",
                             "table_sample_exclusion_audit.csv",
                             "table_confounding_heterogeneity_gate.csv"))),
  longitudinal = list(type = "R", scripts = "12_longitudinal_overcoming.R", deps = "bulk_scores",
                      outputs = c(p("results", "tables", "common_state", c(
                        "table_longitudinal_primary_meta.csv", "table_longitudinal_gate.csv",
                        "table_longitudinal_primary_change_effects.csv")),
                        p("results", "source_data", "longitudinal_paired_patient_data.csv"))),
  singlecell_download = list(type = "shell", scripts = "13_download_singlecell_subset.sh", deps = "preflight",
                             outputs = singlecell_files, download = TRUE, singlecell = TRUE),
  singlecell_analysis = list(type = "R", scripts = "14_singlecell_compartment_localization.R",
                             deps = c("singlecell_download", "gene_sets"),
                             outputs = singlecell_analysis_outputs,
                             singlecell = TRUE),
  bulk_figures = list(type = "R", scripts = "15_make_main_figures.R",
                      deps = c("baseline_models", "longitudinal"),
                      outputs = bulk_figure_outputs),
  supplementary_figures = list(type = "R", scripts = "19_make_supplementary_figures.R",
                               deps = c("confounding_audit", "longitudinal"),
                               outputs = supplementary_figure_outputs),
  singlecell_figure = list(type = "R", scripts = "16_make_singlecell_figure.R", deps = "singlecell_analysis",
                           outputs = c(
                             figure_bundle("Figure5_cross_compartment_localization", 5),
                             supplement_bundle("FigureS5_singlecell_clinical_audit",
                                               "FigureS5_source_data.csv")
                           ),
                           singlecell = TRUE),
  supplementary_tables = list(type = "R", scripts = "20_build_supplementary_tables.R",
                              deps = c("confounding_audit", "longitudinal", "singlecell_analysis"),
                              outputs = supplementary_table_outputs,
                              singlecell = TRUE),
  validate = list(type = "R", scripts = "17_validate_rebuild.R",
                  deps = c("baseline_models", "confounding_audit", "longitudinal", "bulk_figures",
                           "supplementary_figures", "singlecell_analysis", "singlecell_figure",
                           "supplementary_tables"),
                  outputs = p("results", "reproducibility", c("validation_report.csv", "artifact_inventory.csv")))
)

logs_dir <- p("logs", "pipeline")
stamp_dir <- p(".pipeline", "stamps")
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stamp_dir, recursive = TRUE, showWarnings = FALSE)

all_outputs <- function(names) unique(unlist(lapply(names, function(x) stage[[x]]$outputs), use.names = FALSE))
is_fresh <- function(name) {
  info <- stage[[name]]
  if (!length(info$outputs) || any(!file.exists(info$outputs)) || any(file.info(info$outputs)$size <= 0)) return(FALSE)
  script_paths <- p("scripts", info$scripts)
  if (any(!file.exists(script_paths))) return(FALSE)
  dependency_outputs <- all_outputs(info$deps)
  inputs <- c(script_paths, dependency_outputs,
              if (name %in% c("gene_sets", "bulk_download", "singlecell_download")) p("config", "data_manifest.csv") else character())
  inputs <- inputs[file.exists(inputs)]
  if (!length(inputs)) return(TRUE)
  min(file.info(info$outputs)$mtime) >= max(file.info(inputs)$mtime)
}

selected <- stage_names[from_i:to_i]
if (skip_singlecell) selected <- selected[!vapply(selected, function(x) isTRUE(stage[[x]]$singlecell), logical(1))]
if (skip_singlecell && "validate" %in% selected) {
  message("--skip-single-cell also skips the full validation stage; use --to=supplementary_figures for an explicit bulk-only run.")
  selected <- setdiff(selected, "validate")
}

status_rows <- list()
for (name in selected) {
  info <- stage[[name]]
  start <- Sys.time()
  log_path <- file.path(logs_dir, paste0(name, ".log"))

  if (isTRUE(info$download) && offline) {
    missing <- info$outputs[!file.exists(info$outputs) | file.info(info$outputs)$size <= 0]
    if (length(missing)) stop("Offline mode: missing required inputs for ", name, ": ", paste(basename(missing), collapse = ", "))
    action <- "offline-existing"
    message(sprintf("[%s] OFFLINE PASS", name))
  } else if (!force && is_fresh(name)) {
    action <- "fresh-skip"
    message(sprintf("[%s] FRESH; skipped", name))
  } else {
    action <- "ran"
    script_file <- p("scripts", info$scripts)
    message(sprintf("[%s] RUN -> %s", name, log_path))
    executable <- if (info$type == "R") file.path(R.home("bin"), "Rscript") else Sys.which("bash")
    result <- system2(executable, script_file, stdout = log_path, stderr = log_path,
                      env = paste0("PROJECT_ROOT=", project_root))
    if (!identical(result, 0L)) {
      tail_lines <- if (file.exists(log_path)) tail(readLines(log_path, warn = FALSE), 30L) else "<no log>"
      writeLines(tail_lines, stderr())
      stop("Pipeline stage failed: ", name, ". See ", log_path)
    }
    missing <- info$outputs[!file.exists(info$outputs) | file.info(info$outputs)$size <= 0]
    if (length(missing)) stop("Stage ", name, " completed without required outputs: ", paste(missing, collapse = ", "))
    stamp <- data.frame(stage = name, completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                        scripts = paste(info$scripts, collapse = ";"),
                        script_md5 = paste(unname(tools::md5sum(script_file)), collapse = ";"),
                        stringsAsFactors = FALSE)
    utils::write.csv(stamp, file.path(stamp_dir, paste0(name, ".csv")), row.names = FALSE)
  }
  end <- Sys.time()
  status_rows[[length(status_rows) + 1L]] <- data.frame(
    stage = name, action = action, started = format(start, tz = "UTC", usetz = TRUE),
    finished = format(end, tz = "UTC", usetz = TRUE), duration_seconds = as.numeric(difftime(end, start, units = "secs")),
    log = ifelse(file.exists(log_path), log_path, NA_character_), stringsAsFactors = FALSE
  )
}

status <- do.call(rbind, status_rows)
dir.create(p("results", "reproducibility"), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(status, p("results", "reproducibility", "last_pipeline_run.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), p("logs", "run_all_sessionInfo.txt"))
message("PIPELINE COMPLETE: ", paste(selected, collapse = " -> "))
