#!/usr/bin/env Rscript

# Harmonize probe-to-gene mapping and GSVA scoring across all bulk cohorts.
# Scores are standardized using each accession's baseline UC distribution and
# the same center/scale is then applied to longitudinal samples.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
raw_dir <- file.path(project_root, "data", "raw")
derived_dir <- file.path(project_root, "data", "derived", "extended")
tables_dir <- file.path(project_root, "results", "tables", "extended")
logs_dir <- file.path(project_root, "logs")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(AnnotationDbi)
  library(hthgu133pluspm.db)
  library(hgu133plus2.db)
  library(GSVA)
})

registry_path <- file.path(derived_dir, "sample_registry.csv")
if (!file.exists(registry_path)) stop("Run scripts/07_build_cohort_registry.R first")
registry <- utils::read.csv(registry_path, stringsAsFactors = FALSE, check.names = FALSE)
sets <- readRDS(file.path(project_root, "data", "derived", "predefined_gene_sets.rds"))

read_eset <- function(accession) {
  path <- file.path(raw_dir, paste0(accession, "_series_matrix.txt.gz"))
  object <- GEOquery::getGEO(filename = path, getGPL = FALSE)
  if (is.list(object)) object <- object[[1]]
  object
}

collapse_annotation_db <- function(expression, annotation_db) {
  valid <- intersect(rownames(expression), AnnotationDbi::keys(annotation_db, keytype = "PROBEID"))
  annotation <- AnnotationDbi::select(annotation_db, keys = valid, columns = "SYMBOL", keytype = "PROBEID")
  annotation$PROBEID <- as.character(annotation$PROBEID)
  annotation <- annotation[!is.na(annotation$SYMBOL) & nzchar(annotation$SYMBOL), ]
  annotation <- annotation[!duplicated(annotation$PROBEID), ]
  probe_iqr <- apply(expression[annotation$PROBEID, , drop = FALSE], 1L, stats::IQR, na.rm = TRUE)
  annotation$probe_iqr <- probe_iqr[annotation$PROBEID]
  annotation <- annotation[order(annotation$SYMBOL, -annotation$probe_iqr, annotation$PROBEID), ]
  annotation <- annotation[!duplicated(annotation$SYMBOL), ]
  result <- expression[annotation$PROBEID, , drop = FALSE]
  rownames(result) <- annotation$SYMBOL
  result
}

get_gpl6244_map <- function() {
  cache_path <- file.path(derived_dir, "GPL6244_symbol_map.csv")
  if (file.exists(cache_path)) {
    mapping <- utils::read.csv(cache_path, stringsAsFactors = FALSE, colClasses = "character")
    return(mapping)
  }
  platform <- GEOquery::getGEO("GPL6244")
  table <- GEOquery::Table(platform)
  first_assignment <- vapply(strsplit(table$gene_assignment, " /// ", fixed = TRUE),
                             function(x) x[1], character(1))
  symbol <- vapply(strsplit(first_assignment, " // ", fixed = TRUE),
                   function(x) if (length(x) >= 2L) x[2] else NA_character_, character(1))
  mapping <- data.frame(PROBEID = as.character(table$ID), SYMBOL = symbol, stringsAsFactors = FALSE)
  mapping <- mapping[!is.na(mapping$SYMBOL) & mapping$SYMBOL != "---" & nzchar(mapping$SYMBOL), ]
  utils::write.csv(mapping, cache_path, row.names = FALSE)
  mapping
}

collapse_gpl6244 <- function(expression) {
  annotation <- get_gpl6244_map()
  # GPL6244 probe identifiers contain digits only. read.csv() otherwise infers
  # them as numeric and matrix indexing silently switches from names to row
  # positions, which breaks a clean rebuild.
  annotation$PROBEID <- as.character(annotation$PROBEID)
  annotation <- annotation[annotation$PROBEID %in% rownames(expression), ]
  probe_iqr <- apply(expression[annotation$PROBEID, , drop = FALSE], 1L, stats::IQR, na.rm = TRUE)
  annotation$probe_iqr <- probe_iqr[annotation$PROBEID]
  annotation <- annotation[order(annotation$SYMBOL, -annotation$probe_iqr, annotation$PROBEID), ]
  annotation <- annotation[!duplicated(annotation$SYMBOL), ]
  result <- expression[annotation$PROBEID, , drop = FALSE]
  rownames(result) <- annotation$SYMBOL
  result
}

score_sets <- function(expression, gene_sets) {
  present_sets <- lapply(gene_sets, function(genes) intersect(unique(genes), rownames(expression)))
  mapped_n <- lengths(present_sets)
  if (mapped_n[["Inflammation"]] < 100L) {
    stop("GATE FAILED: fewer than 100 primary-signature genes mapped")
  }
  parameter <- GSVA::gsvaParam(exprData = expression, geneSets = present_sets, kcdf = "Gaussian")
  list(scores = t(GSVA::gsva(parameter, verbose = FALSE)), mapped_n = mapped_n)
}

baseline_index <- function(accession, metadata) {
  switch(accession,
         GSE92415 = metadata$visit == "Week 0" & grepl("Ulcerative", metadata$disease),
         GSE16879 = metadata$visit == "Week 0" & tolower(metadata$disease) == "uc",
         GSE23597 = metadata$visit == "W0",
         GSE73661 = metadata$visit == "W0" & metadata$disease == "ulcerative colitis",
         GSE206285 = metadata$disease == "ulcerative colitis",
         stop("Unknown accession: ", accession))
}

platform_method <- list(
  GSE92415 = function(x) collapse_annotation_db(x, hthgu133pluspm.db),
  GSE16879 = function(x) collapse_annotation_db(x, hgu133plus2.db),
  GSE23597 = function(x) collapse_annotation_db(x, hgu133plus2.db),
  GSE73661 = collapse_gpl6244,
  GSE206285 = function(x) collapse_annotation_db(x, hthgu133pluspm.db)
)

mapping_rows <- list()
for (accession in names(platform_method)) {
  message("PREPARING ", accession)
  eset <- read_eset(accession)
  expression <- Biobase::exprs(eset)
  expression_gene <- platform_method[[accession]](expression)
  scored <- score_sets(expression_gene, sets)

  metadata <- registry[registry$accession == accession, ]
  metadata <- metadata[match(colnames(expression_gene), metadata$gsm), ]
  if (anyNA(metadata$gsm)) stop("GATE FAILED [", accession, "]: expression/registry mismatch")
  score_frame <- data.frame(gsm = rownames(scored$scores), scored$scores,
                            check.names = FALSE, stringsAsFactors = FALSE)
  merged <- cbind(metadata, score_frame[match(metadata$gsm, score_frame$gsm), setdiff(names(score_frame), "gsm"), drop = FALSE])

  base <- baseline_index(accession, merged)
  for (set_name in names(sets)) {
    center <- mean(merged[[set_name]][base], na.rm = TRUE)
    spread <- stats::sd(merged[[set_name]][base], na.rm = TRUE)
    if (!is.finite(spread) || spread <= 0) stop("GATE FAILED [", accession, "/", set_name, "]: invalid SD")
    merged[[paste0(set_name, "_z")]] <- (merged[[set_name]] - center) / spread
    attr(merged[[paste0(set_name, "_z")]], "center") <- center
    attr(merged[[paste0(set_name, "_z")]], "scale") <- spread
  }

  saveRDS(merged, file.path(derived_dir, paste0(tolower(accession), "_all_scores.rds")))
  utils::write.csv(merged, file.path(derived_dir, paste0(tolower(accession), "_all_scores.csv")), row.names = FALSE)
  mapping_rows[[accession]] <- data.frame(
    accession = accession, expression_features = nrow(expression), mapped_symbols = nrow(expression_gene),
    gene_set = names(scored$mapped_n), mapped_gene_set_size = as.integer(scored$mapped_n),
    stringsAsFactors = FALSE
  )
  rm(eset, expression, expression_gene, scored, metadata, score_frame, merged)
  gc()
}

mapping_table <- do.call(rbind, mapping_rows)
utils::write.csv(mapping_table, file.path(tables_dir, "table_gene_mapping_audit.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "08_prepare_extended_cohorts_sessionInfo.txt"))
message("ALL PREPROCESSING GATES PASSED")
