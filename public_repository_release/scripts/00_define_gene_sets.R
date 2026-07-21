#!/usr/bin/env Rscript

# Freeze the exact gene-set membership used by every bulk and single-cell
# stage. The tracked CSV is the primary source; msigdbr is used only once to
# create it when rebuilding a fresh project copy.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
config_dir <- file.path(project_root, "config")
derived_dir <- file.path(project_root, "data", "derived")
logs_dir <- file.path(project_root, "logs")
dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

csv_path <- file.path(config_dir, "predefined_gene_sets.csv")
rds_path <- file.path(derived_dir, "predefined_gene_sets.rds")

if (file.exists(csv_path)) {
  long <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  required <- c("gene_set", "gene_symbol")
  if (!all(required %in% names(long))) stop("Malformed frozen gene-set CSV: ", csv_path)
  sets <- split(long$gene_symbol, long$gene_set)
} else {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Missing frozen gene-set CSV and package msigdbr is unavailable")
  }
  msig <- msigdbr::msigdbr(species = "Homo sapiens")
  gene_set <- function(name) unique(msig$gene_symbol[msig$gs_name == name])
  sets <- list(
    Inflammation = gene_set("HALLMARK_INFLAMMATORY_RESPONSE"),
    TNFaNFKB = gene_set("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
    IL6JAKSTAT3 = gene_set("HALLMARK_IL6_JAK_STAT3_SIGNALING"),
    OxidativePhos_neg = gene_set("HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
    Random_neg = { set.seed(2026); sample(unique(msig$gene_symbol), 200L) }
  )
  long <- do.call(rbind, lapply(names(sets), function(nm) {
    data.frame(gene_set = nm, gene_symbol = sets[[nm]], stringsAsFactors = FALSE)
  }))
  rownames(long) <- NULL
  utils::write.csv(long, csv_path, row.names = FALSE)
}

expected_sets <- c("Inflammation", "TNFaNFKB", "IL6JAKSTAT3", "OxidativePhos_neg", "Random_neg")
if (!all(expected_sets %in% names(sets))) stop("Frozen gene-set file is missing required sets")
if (length(unique(sets$Inflammation)) != 200L) stop("Inflammation gene set must contain 200 unique genes")
sets <- lapply(sets[expected_sets], unique)
saveRDS(sets, rds_path)
utils::write.csv(data.frame(gene_set = names(sets), genes = lengths(sets)),
                 file.path(config_dir, "predefined_gene_set_sizes.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "00_define_gene_sets_sessionInfo.txt"))
message("Frozen gene sets ready: ", paste(names(sets), lengths(sets), sep = "=", collapse = ", "))
