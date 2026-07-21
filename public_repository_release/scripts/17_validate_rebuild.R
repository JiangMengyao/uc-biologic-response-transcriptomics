#!/usr/bin/env Rscript

# Machine-readable end-of-pipeline validation. This checks gates, artifact
# presence, source-data readability and the final figure export bundle.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
out_dir <- file.path(project_root, "results", "reproducibility")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

checks <- list()
add_check <- function(category, check, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(category = category, check = check,
                                               passed = isTRUE(passed), detail = as.character(detail),
                                               stringsAsFactors = FALSE)
}

baseline_gate <- utils::read.csv(file.path(project_root, "results", "tables", "common_state", "table_cross_mechanism_gate.csv"))
longitudinal_gate <- utils::read.csv(file.path(project_root, "results", "tables", "common_state", "table_longitudinal_gate.csv"))
confounding_gate <- utils::read.csv(file.path(project_root, "results", "tables", "common_state", "table_confounding_heterogeneity_gate.csv"))
add_check("scientific_gate", "cross_mechanism", all(baseline_gate$passed), paste(baseline_gate$mechanism, baseline_gate$passed, collapse = "; "))
add_check("scientific_gate", "longitudinal", all(longitudinal_gate$passed), paste(longitudinal_gate$observed, collapse = "; "))
add_check("scientific_gate", "confounding_heterogeneity", all(confounding_gate$passed),
          paste(confounding_gate$gate, confounding_gate$passed, collapse = "; "))

availability_path <- file.path(project_root, "results", "tables", "common_state",
                               "table_covariate_availability_matrix.csv")
adjusted_path <- file.path(project_root, "results", "tables", "common_state",
                           "table_adjusted_baseline_associations.csv")
stratified_path <- file.path(project_root, "results", "tables", "common_state",
                             "table_stratified_sensitivity_meta.csv")
moderator_path <- file.path(project_root, "results", "tables", "common_state",
                            "table_mechanism_moderator_test.csv")
exclusion_path <- file.path(project_root, "results", "tables", "common_state",
                            "table_sample_exclusion_audit.csv")
if (all(file.exists(c(availability_path, adjusted_path, stratified_path, moderator_path, exclusion_path)))) {
  availability <- utils::read.csv(availability_path, check.names = FALSE)
  adjusted_audit <- utils::read.csv(adjusted_path, check.names = FALSE)
  stratified_audit <- utils::read.csv(stratified_path, check.names = FALSE)
  moderator_audit <- utils::read.csv(moderator_path, check.names = FALSE)
  exclusion_audit <- utils::read.csv(exclusion_path, check.names = FALSE)
  adjusted_rows <- adjusted_audit[adjusted_audit$adjustment_status == "Estimated", ]
  add_check("confounding_audit", "covariate_matrix_complete",
            nrow(availability) == 56L && length(unique(availability$cohort)) == 7L &&
              length(unique(availability$covariate)) == 8L,
            paste(nrow(availability), "rows"))
  add_check("confounding_audit", "adjusted_models_direction",
            nrow(adjusted_rows) == 5L && all(adjusted_rows$n_complete == adjusted_rows$n_primary) &&
              all(adjusted_rows$adjusted_or < 1),
            paste(adjusted_rows$cohort, signif(adjusted_rows$adjusted_or, 4), collapse = "; "))
  add_check("confounding_audit", "requested_strata_complete",
            identical(sort(unique(stratified_audit$dimension)),
                      sort(c("endpoint_type", "evaluation_time", "mechanism", "biopsy_site_group"))),
            paste(unique(stratified_audit$dimension), collapse = "; "))
  add_check("confounding_audit", "between_mechanism_heterogeneity",
            nrow(moderator_audit) == 2L && all(moderator_audit$moderator_p >= 0.10) &&
              all(moderator_audit$max_min_or_ratio < 2),
            paste(moderator_audit$analysis,
                  sprintf("P=%.4f ratio=%.3f", moderator_audit$moderator_p,
                          moderator_audit$max_min_or_ratio), collapse = "; "))
  flow_totals <- stats::aggregate(n_samples ~ cohort + raw_accession_samples,
                                  data = exclusion_audit, FUN = sum)
  add_check("confounding_audit", "sample_flow_reconciled",
            nrow(flow_totals) == 7L && all(flow_totals$n_samples == flow_totals$raw_accession_samples),
            paste(flow_totals$cohort, flow_totals$n_samples, collapse = "; "))
} else {
  add_check("confounding_audit", "audit_outputs", FALSE,
            "covariate, adjusted, moderator, stratified or exclusion table missing")
}

main_dir <- file.path(project_root, "results", "figures", "main")
preview_dir <- file.path(project_root, "results", "figures", "previews")
figure_stems <- c("Figure1_evidence_architecture", "Figure2_cross_mechanism_prognostic_state",
                  "Figure3_prognostic_vs_predictive_RCT", "Figure4_longitudinal_overcoming",
                  "Figure5_cross_compartment_localization")
for (i in seq_along(figure_stems)) {
  stem <- figure_stems[i]
  files <- c(file.path(main_dir, paste0(stem, c(".svg", ".pdf", ".tiff"))),
             file.path(preview_dir, paste0(stem, ".png")))
  pass <- all(file.exists(files)) && all(file.info(files)$size > 0)
  add_check("figure_bundle", paste0("Figure", i), pass, paste(basename(files), collapse = "; "))
}

supp_stems <- c(
  "FigureS1_cohort_inclusion_flow",
  "FigureS2_baseline_sensitivity",
  "FigureS3_longitudinal_sensitivity",
  "FigureS4_high_state_rct_audit",
  "FigureS5_singlecell_clinical_audit"
)
for (i in seq_along(supp_stems)) {
  stem <- supp_stems[i]
  files <- c(file.path(project_root, "results", "figures", "supplementary",
                       paste0(stem, c(".svg", ".pdf", ".tiff"))),
             file.path(preview_dir, paste0(stem, ".png")))
  add_check("figure_bundle", paste0("FigureS", i),
            all(file.exists(files)) && all(file.info(files)$size > 0),
            paste(basename(files), collapse = "; "))
}

for (i in 1:5) {
  path <- file.path(project_root, "results", "source_data", paste0("Figure", i, "_source_data.csv"))
  ok <- file.exists(path) && file.info(path)$size > 0
  rows <- if (ok) nrow(utils::read.csv(path, check.names = FALSE)) else 0L
  add_check("source_data", paste0("Figure", i), ok && rows > 0, paste(rows, "rows"))
}
supp_sources <- file.path(project_root, "results", "source_data",
                          paste0("FigureS", 1:5, "_source_data.csv"))
for (i in seq_along(supp_sources)) {
  ok <- file.exists(supp_sources[i]) && file.info(supp_sources[i])$size > 0
  rows <- if (ok) nrow(utils::read.csv(supp_sources[i], check.names = FALSE)) else 0L
  add_check("source_data", paste0("FigureS", i), ok && rows > 0, paste(rows, "rows"))
}
supp_source <- supp_sources[5]

supplementary_table_names <- c(
  "TableS1_cohort_registry.csv",
  "TableS2_baseline_associations.csv",
  "TableS3_RCT_prognostic_predictive.csv",
  "TableS4_longitudinal_sensitivities.csv",
  "TableS5_singlecell_baseline_sensitivity.csv",
  "TableS6_singlecell_state_baseline.csv",
  "TableS7_singlecell_longitudinal.csv"
)
supplementary_table_expected_rows <- c(10L, 11L, 9L, 51L, 18L, 34L, 37L)
supplementary_table_dir <- file.path(project_root, "results", "supplementary_tables")
for (i in seq_along(supplementary_table_names)) {
  path <- file.path(supplementary_table_dir, supplementary_table_names[i])
  ok <- file.exists(path) && file.info(path)$size > 0
  rows <- if (ok) nrow(utils::read.csv(path, check.names = FALSE)) else 0L
  add_check("supplementary_tables", paste0("TableS", i, "_complete"),
            ok && rows == supplementary_table_expected_rows[i],
            paste(rows, "rows; expected", supplementary_table_expected_rows[i]))
}
supplementary_workbook <- file.path(supplementary_table_dir, "Supplementary_Tables_S1-S7.xlsx")
workbook_ok <- file.exists(supplementary_workbook) && file.info(supplementary_workbook)$size > 0
workbook_signature <- if (workbook_ok) {
  con <- file(supplementary_workbook, "rb")
  on.exit(close(con), add = TRUE)
  rawToChar(readBin(con, what = "raw", n = 2L))
} else ""
add_check("supplementary_tables", "formatted_workbook",
          workbook_ok && identical(workbook_signature, "PK"),
          if (workbook_ok) paste(file.info(supplementary_workbook)$size, "bytes") else "missing")

manifest <- utils::read.csv(file.path(project_root, "config", "data_manifest.csv"), stringsAsFactors = FALSE)
checksum_rows <- manifest[nzchar(manifest$expected_md5), ]
for (i in seq_len(nrow(checksum_rows))) {
  path <- if (checksum_rows$dataset[i] == "GSE282122") {
    file.path(project_root, "data", "raw", "GSE282122", checksum_rows$file[i])
  } else {
    file.path(project_root, "data", "raw", checksum_rows$file[i])
  }
  observed <- if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_
  add_check("checksum", checksum_rows$file[i], identical(observed, checksum_rows$expected_md5[i]),
            paste0("observed=", observed, "; expected=", checksum_rows$expected_md5[i]))
}

localization_path <- file.path(project_root, "results", "tables", "common_state",
                               "table_singlecell_compartment_localization.csv")
pseudobulk_path <- file.path(project_root, "results", "source_data", "singlecell_compartment_pseudobulk.csv")
sample_state_path <- file.path(project_root, "results", "source_data", "singlecell_sample_state_pseudobulk.csv")
sample_compartment_path <- file.path(project_root, "results", "source_data", "singlecell_sample_compartment_pseudobulk.csv")
pairing_path <- file.path(project_root, "results", "tables", "common_state",
                          "table_singlecell_longitudinal_pairing_audit.csv")
longitudinal_path <- file.path(project_root, "results", "tables", "common_state",
                               "table_singlecell_compartment_longitudinal_effects.csv")
sensitivity_path <- file.path(project_root, "results", "tables", "common_state",
                              "table_singlecell_compartment_outcome_sensitivity.csv")
scaling_path <- file.path(project_root, "data", "derived", "singlecell", "cross_compartment_scaling.rds")
if (all(file.exists(c(localization_path, pseudobulk_path, sample_state_path,
                      sample_compartment_path, pairing_path, longitudinal_path,
                      sensitivity_path, scaling_path)))) {
  localization <- utils::read.csv(localization_path, check.names = FALSE)
  pseudobulk <- utils::read.csv(pseudobulk_path, check.names = FALSE)
  sample_state <- utils::read.csv(sample_state_path, check.names = FALSE)
  sample_compartment <- utils::read.csv(sample_compartment_path, check.names = FALSE)
  pairing <- utils::read.csv(pairing_path, check.names = FALSE)
  longitudinal <- utils::read.csv(longitudinal_path, check.names = FALSE)
  sensitivity <- utils::read.csv(sensitivity_path, check.names = FALSE)
  scaling <- readRDS(scaling_path)
  expected_compartments <- c("Colonic epithelium", "Myeloid", "Fibroblast/pericyte")
  observed_compartments <- sort(unique(as.character(localization$compartment)))
  baseline_donors <- tapply(pseudobulk$donor[pseudobulk$treatment_time == "Pre"],
                            pseudobulk$compartment[pseudobulk$treatment_time == "Pre"],
                            function(z) length(unique(z)))
  add_check("singlecell_gate", "three_compartments",
            identical(sort(expected_compartments), observed_compartments),
            paste(observed_compartments, collapse = "; "))
  add_check("singlecell_gate", "baseline_donors", all(baseline_donors[expected_compartments] >= 8L),
            paste(names(baseline_donors), baseline_donors, collapse = "; "))
  add_check("singlecell_gate", "unique_pseudobulks", !anyDuplicated(pseudobulk$group_id),
            paste(nrow(pseudobulk), "groups"))
  add_check("singlecell_gate", "sample_first_units",
            !anyDuplicated(sample_state$group_id) && !anyDuplicated(sample_compartment$group_id) &&
              all(sample_state$n_samples == 1L) && all(sample_compartment$n_samples == 1L),
            paste(nrow(sample_state), "sample-state;", nrow(sample_compartment), "sample-compartment"))
  add_check("singlecell_gate", "mapped_signature", length(scaling$valid_genes) >= 100L,
            paste(length(scaling$valid_genes), "variable common genes"))
  add_check("singlecell_gate", "frozen_signature",
            identical(scaling$signature_md5, "e072b3b52623f871f902da42aaed04c9") &&
              identical(scaling$algorithm_version, "cross-compartment-sample-first-csr-v2"),
            paste(scaling$algorithm_version, scaling$signature_md5))
  add_check("singlecell_gate", "match_yes_exact_pairing",
            all(pairing$exact_pairing) && all(pairing$paired_donors >= 6L) &&
              all(pairing$remission_donors >= 2L) && all(pairing$nonremission_donors >= 2L),
            paste(pairing$compartment, pairing$paired_donors, collapse = "; "))
  add_check("singlecell_gate", "longitudinal_expected_direction",
            all(is.finite(longitudinal$estimate)) && all(longitudinal$estimate < 0),
            paste(longitudinal$compartment, signif(longitudinal$estimate, 4), collapse = "; "))
  legacy_epi <- sensitivity[sensitivity$analysis_label == "Legacy state-cell-weighted" &
                              sensitivity$compartment == "Colonic epithelium", ]
  add_check("singlecell_gate", "legacy_epithelial_preserved",
            nrow(legacy_epi) == 1L && abs(legacy_epi$odds_ratio - 2.89692486714913) < 1e-6 &&
              legacy_epi$direction_vs_bulk == "Discordant" && legacy_epi$fdr > .05,
            if (nrow(legacy_epi)) paste("OR", legacy_epi$odds_ratio, "FDR", legacy_epi$fdr) else "missing")
} else {
  add_check("singlecell_gate", "cross_compartment_outputs", FALSE,
            "sample-first localization, pairing, sensitivity or scaling output missing")
}

figure5_source <- file.path(project_root, "results", "source_data", "Figure5_source_data.csv")
figure1_source <- file.path(project_root, "results", "source_data", "Figure1_source_data.csv")
if (file.exists(figure1_source)) {
  fig1 <- utils::read.csv(figure1_source, check.names = FALSE)
  expected_fig1_panels <- c("a_cohort_lanes", "b_outcome_counts",
                            "c_evidence_matrix", "d_claim_hierarchy")
  observed_fig1_panels <- unique(as.character(fig1$panel))
  add_check("figure_contract", "Figure1_panel_source_mapping",
            identical(observed_fig1_panels, expected_fig1_panels),
            paste(observed_fig1_panels, collapse = "; "))
}
if (file.exists(figure5_source)) {
  fig5 <- utils::read.csv(figure5_source, check.names = FALSE)
  main_panels <- unique(as.character(fig5$panel))
  localization_only <- all(grepl("^(a_|b_|c_)", main_panels)) &&
    !any(grepl("outcome|longitudinal|remission|change", main_panels, ignore.case = TRUE))
  add_check("figure_contract", "Figure5_localization_only", localization_only,
            paste(main_panels, collapse = "; "))
} else {
  add_check("figure_contract", "Figure5_localization_only", FALSE, "Figure5 source data missing")
}

supplementary_panel_contracts <- list(
  FigureS1 = c("a_accession_nodes", "a_cohort_nodes", "a_branch_nodes",
               "a_accession_edges", "a_branch_edges", "b_sample_disposition",
               "c_evidence_branch_counts"),
  FigureS2 = c("a_complete_cohort_forest", "b_cohort_spline_predictions",
               "c_nonlinearity_tests", "d_leave_one_out", "e_module_meta_analysis"),
  FigureS3 = c("a_analytic_and_bootstrap", "b_mixed_models",
               "c_module_heatmap", "d_module_distribution"),
  FigureS4 = c("a_flow_nodes", "a_flow_edges", "b_transition_rates",
               "c_reversal_response_or", "d_active_placebo_change")
)
for (i in seq_along(supplementary_panel_contracts)) {
  path <- file.path(project_root, "results", "source_data",
                    paste0("FigureS", i, "_source_data.csv"))
  observed <- if (file.exists(path)) {
    unique(as.character(utils::read.csv(path, check.names = FALSE)$panel))
  } else character()
  expected <- supplementary_panel_contracts[[i]]
  add_check("figure_contract", paste0("FigureS", i, "_panel_source_mapping"),
            identical(observed, expected), paste(observed, collapse = "; "))
}

if (file.exists(supp_source)) {
  figs5 <- utils::read.csv(supp_source, check.names = FALSE)
  observed_supp_panels <- unique(as.character(figs5$panel))
  expected_supp_panels <- c(
    "a_baseline_primary",
    "b_match_yes_longitudinal_effects",
    "c_baseline_sensitivity_complete",
    "d_match_yes_donor_changes"
  )
  add_check("figure_contract", "FigureS5_panel_source_mapping",
            identical(observed_supp_panels, expected_supp_panels),
            paste(observed_supp_panels, collapse = "; "))
}

manuscript_dir <- file.path(project_root, "manuscript")
required_manuscript_files <- c(
  "README.md",
  "00_terminology_and_scope.md",
  "01_methods_and_results.md",
  "02_figure_legends.md",
  "03_claim_evidence_matrix.md",
  "04_supplementary_materials.md",
  "supplementary_manifest.csv",
  "05_presubmission_mock_review.md",
  "06_title_abstract_introduction_discussion.md"
)
manuscript_paths <- file.path(manuscript_dir, required_manuscript_files)
add_check(
  "manuscript_package", "required_files",
  all(file.exists(manuscript_paths)) && all(file.info(manuscript_paths)$size > 0),
  paste(required_manuscript_files, collapse = "; ")
)
if (file.exists(file.path(manuscript_dir, "supplementary_manifest.csv"))) {
  supplementary_manifest <- utils::read.csv(
    file.path(manuscript_dir, "supplementary_manifest.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  referenced_paths <- file.path(project_root, supplementary_manifest$primary_repository_path)
  add_check(
    "manuscript_package", "supplementary_manifest_paths",
    nrow(supplementary_manifest) >= 15L && all(file.exists(referenced_paths)),
    paste(nrow(supplementary_manifest), "items;",
          sum(file.exists(referenced_paths)), "repository objects found")
  )
}
methods_results_path <- file.path(manuscript_dir, "01_methods_and_results.md")
legends_path <- file.path(manuscript_dir, "02_figure_legends.md")
claims_path <- file.path(manuscript_dir, "03_claim_evidence_matrix.md")
if (all(file.exists(c(methods_results_path, legends_path, claims_path)))) {
  methods_results_text <- paste(readLines(methods_results_path, warn = FALSE), collapse = "\n")
  legends_text <- paste(readLines(legends_path, warn = FALSE), collapse = "\n")
  claims_text <- paste(readLines(claims_path, warn = FALSE), collapse = "\n")
  required_methods_terms <- c(
    "Sample-first single-cell localization and donor-level inference",
    "Match=Yes",
    "donor remained the inferential unit",
    "do not independently test clinical prediction",
    "OR of 2.90",
    "FDR=0.123"
  )
  add_check(
    "manuscript_package", "methods_results_claim_boundary",
    all(vapply(required_methods_terms, grepl, logical(1), x = methods_results_text,
               fixed = TRUE)),
    paste(required_methods_terms, collapse = "; ")
  )
  required_legend_terms <- c(
    "## Figure 1", "## Figure 2", "## Figure 3", "## Figure 4",
    "## Figure 5", "## Figure S1", "## Figure S2", "## Figure S3",
    "## Figure S4", "## Figure S5", "aggregation-sensitive",
    "Figure1_source_data.csv", "Figure2_source_data.csv",
    "Figure3_source_data.csv", "Figure4_source_data.csv",
    "Figure5_source_data.csv", "FigureS1_source_data.csv",
    "FigureS2_source_data.csv", "FigureS3_source_data.csv",
    "FigureS4_source_data.csv", "FigureS5_source_data.csv"
  )
  add_check(
    "manuscript_package", "figure_legends_complete",
    all(vapply(required_legend_terms, grepl, logical(1), x = legends_text,
               fixed = TRUE)),
    paste(required_legend_terms, collapse = "; ")
  )
  required_claim_terms <- c(
    "shared adverse prognostic state",
    "deeper molecular reversal",
    "localizes across epithelial, myeloid and fibroblast/pericyte compartments",
    "does not independently validate baseline clinical prediction"
  )
  add_check(
    "manuscript_package", "three_level_claim_contract",
    all(vapply(required_claim_terms, grepl, logical(1), x = claims_text,
               fixed = TRUE)),
    paste(required_claim_terms, collapse = "; ")
  )
}

remaining_text_path <- file.path(manuscript_dir, "06_title_abstract_introduction_discussion.md")
if (file.exists(remaining_text_path)) {
  remaining_text <- paste(readLines(remaining_text_path, warn = FALSE), collapse = "\n")
  required_remaining_terms <- c(
    "# Title",
    "shared adverse prognostic state",
    "# Abstract",
    "# Introduction — final study-aim paragraph",
    "# Discussion",
    "## Limitations",
    "# Data Availability",
    "# Code Availability",
    "do not establish a drug-selection biomarker",
    "aggregation-sensitive"
  )
  add_check(
    "manuscript_package", "remaining_sections_claim_boundary",
    all(vapply(required_remaining_terms, grepl, logical(1), x = remaining_text,
               fixed = TRUE)),
    paste(required_remaining_terms, collapse = "; ")
  )
}

report <- do.call(rbind, checks)
utils::write.csv(report, file.path(out_dir, "validation_report.csv"), row.names = FALSE)

# Submission/revision handoff inventory. Raw inputs remain covered by the
# locked data manifest; this inventory fingerprints code, documentation and
# generated analytical artifacts from the validated rebuild.
inventory_files <- unique(c(
  file.path(project_root, "README.md"),
  list.files(file.path(project_root, "scripts"), pattern = "\\.(R|sh)$", full.names = TRUE),
  list.files(file.path(project_root, "config"), full.names = TRUE),
  list.files(file.path(project_root, "docs"), pattern = "\\.(md|csv)$", full.names = TRUE),
  list.files(file.path(project_root, "manuscript"), pattern = "\\.(md|csv)$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(project_root, "environment"), full.names = TRUE),
  list.files(file.path(project_root, "results", "tables"), pattern = "\\.csv$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(project_root, "results", "supplementary_tables"),
             pattern = "\\.(csv|xlsx)$", full.names = TRUE),
  list.files(file.path(project_root, "results", "source_data"), pattern = "\\.csv$", full.names = TRUE),
  list.files(file.path(project_root, "results", "figures", "main"), pattern = "\\.(svg|pdf|tiff)$", full.names = TRUE),
  list.files(file.path(project_root, "results", "figures", "supplementary"), pattern = "\\.(svg|pdf|tiff)$", full.names = TRUE),
  list.files(file.path(project_root, "results", "figures", "previews"), pattern = "\\.png$", full.names = TRUE)
))
inventory_files <- inventory_files[file.exists(inventory_files) & !dir.exists(inventory_files)]
info <- file.info(inventory_files)
inventory <- data.frame(
  relative_path = substring(inventory_files, nchar(project_root) + 2L),
  bytes = info$size,
  modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
  md5 = unname(tools::md5sum(inventory_files)),
  stringsAsFactors = FALSE
)
inventory <- inventory[order(inventory$relative_path), ]
utils::write.csv(inventory, file.path(out_dir, "artifact_inventory.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "validation_sessionInfo.txt"))
if (!all(report$passed)) {
  failed <- report[!report$passed, ]
  stop("REBUILD VALIDATION FAILED: ", paste(failed$check, collapse = ", "))
}
message("REBUILD VALIDATION PASS: ", nrow(report), " checks")
