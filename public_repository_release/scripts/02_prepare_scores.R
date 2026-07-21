#!/usr/bin/env Rscript

# Rebuild the two analysis datasets and pre-specified pathway scores directly
# from the public GEO matrices.  No temporary metadata files are required.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
raw_dir <- file.path(project_root, "data", "raw")
derived_dir <- file.path(project_root, "data", "derived")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(AnnotationDbi)
  library(hthgu133pluspm.db)
  library(hgu133plus2.db)
  library(GSVA)
  library(msigdbr)
})

read_eset <- function(accession) {
  file <- file.path(raw_dir, paste0(accession, "_series_matrix.txt.gz"))
  stopifnot(file.exists(file))
  eset <- GEOquery::getGEO(filename = file, getGPL = FALSE)
  if (is.list(eset)) eset <- eset[[1]]
  eset
}

map_to_symbols <- function(expression, annotation_db) {
  annotation <- AnnotationDbi::select(
    annotation_db, keys = rownames(expression), columns = "SYMBOL", keytype = "PROBEID"
  )
  annotation <- annotation[!duplicated(annotation$PROBEID), ]
  rownames(expression) <- annotation$SYMBOL[match(rownames(expression), annotation$PROBEID)]
  expression <- expression[!is.na(rownames(expression)) & !duplicated(rownames(expression)), , drop = FALSE]
  expression
}

score_gene_set <- function(expression, genes, label) {
  present <- intersect(genes, rownames(expression))
  if (length(present) < 10L) stop("Too few genes mapped for ", label, ": ", length(present))
  # GSVA 2.x API. The fallback is retained only for package/API failures.
  tryCatch({
    parameter <- GSVA::gsvaParam(
      exprData = expression, geneSets = stats::setNames(list(present), label), kcdf = "Gaussian"
    )
    as.numeric(GSVA::gsva(parameter, verbose = FALSE)[1, ])
  }, error = function(e) {
    message("GSVA failed for ", label, "; using documented fallback: ", conditionMessage(e))
    z <- t(scale(t(expression)))
    colMeans(z[present, , drop = FALSE], na.rm = TRUE)
  })
}

gse92415 <- read_eset("GSE92415")
p1 <- Biobase::pData(gse92415)
meta_base <- data.frame(
  gsm = rownames(p1),
  disease = p1[["disease:ch1"]],
  treatment = tolower(p1[["treatment:ch1"]]),
  visit = p1[["visit:ch1"]],
  response = p1[["wk6response:ch1"]],
  subject = p1[["subject:ch1"]],
  age = suppressWarnings(as.numeric(p1[["age:ch1"]])),
  mayo_score = suppressWarnings(as.numeric(p1[["mayo score:ch1"]])),
  stringsAsFactors = FALSE
)
meta_base <- meta_base[
  meta_base$visit == "Week 0" &
    grepl("Ulcerative", meta_base$disease) &
    meta_base$response %in% c("Yes", "No"), , drop = FALSE
]
meta_base$response_binary <- as.integer(meta_base$response == "Yes")
meta_base$treatment_binary <- as.integer(meta_base$treatment == "golimumab")
stopifnot(nrow(meta_base) == 87L, length(unique(meta_base$subject)) == 87L)
stopifnot(identical(as.integer(table(meta_base$treatment, meta_base$response)), c(27L, 17L, 32L, 11L)))

expr1 <- map_to_symbols(Biobase::exprs(gse92415), hthgu133pluspm.db)
expr1_base <- expr1[, meta_base$gsm, drop = FALSE]

gse16879 <- read_eset("GSE16879")
p2 <- Biobase::pData(gse16879)
meta_uc <- data.frame(
  gsm = rownames(p2),
  disease = p2[["disease:ch1"]],
  response = p2[["response to infliximab:ch1"]],
  visit = p2[["before or after first infliximab treatment:ch1"]],
  stringsAsFactors = FALSE
)
meta_uc <- meta_uc[
  tolower(meta_uc$disease) %in% c("uc", "ulcerative colitis") &
    meta_uc$visit == "Before first infliximab treatment" &
    meta_uc$response %in% c("Yes", "No"), , drop = FALSE
]
meta_uc$response_binary <- as.integer(meta_uc$response == "Yes")
stopifnot(nrow(meta_uc) == 24L, sum(meta_uc$response_binary) == 8L)

expr2 <- map_to_symbols(Biobase::exprs(gse16879), hgu133plus2.db)
expr2_base <- expr2[, meta_uc$gsm, drop = FALSE]

msig <- msigdbr::msigdbr(species = "Homo sapiens")
gene_set <- function(name) unique(msig$gene_symbol[msig$gs_name == name])
sets <- list(
  Inflammation = gene_set("HALLMARK_INFLAMMATORY_RESPONSE"),
  TNFaNFKB = gene_set("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  IL6JAKSTAT3 = gene_set("HALLMARK_IL6_JAK_STAT3_SIGNALING"),
  OxidativePhos_neg = gene_set("HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
  Random_neg = { set.seed(2026); sample(unique(msig$gene_symbol), 200L) }
)

score_both <- function(genes, label) {
  list(
    main = score_gene_set(expr1_base, genes, label),
    external = score_gene_set(expr2_base, genes, label)
  )
}
scores <- lapply(names(sets), function(name) score_both(sets[[name]], name))
names(scores) <- names(sets)
for (name in names(scores)) {
  meta_base[[name]] <- scores[[name]]$main
  meta_uc[[name]] <- scores[[name]]$external
}

utils::write.csv(meta_base, file.path(derived_dir, "gse92415_baseline_metadata_scores.csv"), row.names = FALSE)
utils::write.csv(meta_uc, file.path(derived_dir, "gse16879_baseline_metadata_scores.csv"), row.names = FALSE)
saveRDS(meta_base, file.path(derived_dir, "gse92415_baseline_metadata_scores.rds"))
saveRDS(meta_uc, file.path(derived_dir, "gse16879_baseline_metadata_scores.rds"))
saveRDS(sets, file.path(derived_dir, "predefined_gene_sets.rds"))

writeLines(capture.output(sessionInfo()), file.path(project_root, "logs", "02_prepare_scores_sessionInfo.txt"))
message("Prepared GSE92415 baseline n=", nrow(meta_base), " and GSE16879 UC baseline n=", nrow(meta_uc))
