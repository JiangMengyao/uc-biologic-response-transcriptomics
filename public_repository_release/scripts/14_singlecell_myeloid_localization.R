#!/usr/bin/env Rscript

# Donor-aware localization of the locked inflammatory-response gene set in the
# GSE282122 myeloid compartment. The published h5ad stores a CSR log-normalized
# matrix. Only signature columns are materialized, avoiding a full dense object.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
h5_path <- file.path(project_root, "data", "raw", "GSE282122", "myeloid_final.h5ad")
derived_dir <- file.path(project_root, "data", "derived", "singlecell")
tables_dir <- file.path(project_root, "results", "tables", "common_state")
source_dir <- file.path(project_root, "results", "source_data")
logs_dir <- file.path(project_root, "logs")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(h5_path)) stop("Run scripts/13_download_singlecell_subset.sh first")

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(logistf)
})
set.seed(20260720)

decode_categorical <- function(name) {
  base <- paste0("/obs/", name)
  categories <- rhdf5::h5read(h5_path, paste0(base, "/categories"))
  codes <- rhdf5::h5read(h5_path, paste0(base, "/codes"))
  out <- rep(NA_character_, length(codes))
  ok <- codes >= 0L
  out[ok] <- categories[codes[ok] + 1L]
  out
}

obs <- data.frame(
  cell_id = rhdf5::h5read(h5_path, "/obs/_index"),
  donor = decode_categorical("Patient"),
  sample_id = decode_categorical("sample_id"),
  disease = decode_categorical("Disease"),
  treatment_time = decode_categorical("Treatment"),
  remission = decode_categorical("Remission_status"),
  inflammation = decode_categorical("Inflammation"),
  inflammation_score = rhdf5::h5read(h5_path, "/obs/Inflammation_score"),
  site = decode_categorical("Site"),
  state = decode_categorical("final_analysis"),
  stringsAsFactors = FALSE
)
obs$umap1 <- rhdf5::h5read(h5_path, "/obsm/X_umap/0")
obs$umap2 <- rhdf5::h5read(h5_path, "/obsm/X_umap/1")

gene_symbols <- rhdf5::h5read(h5_path, "/var/gene_symbol")
sets <- readRDS(file.path(project_root, "data", "derived", "predefined_gene_sets.rds"))
target_symbols <- intersect(unique(sets$Inflammation), unique(gene_symbols))
if (length(target_symbols) < 100L) stop("SINGLE-CELL MAPPING GATE FAILED: fewer than 100 genes mapped")
feature_pos <- match(target_symbols, gene_symbols)
target_symbols <- target_symbols[!is.na(feature_pos)]
feature_pos <- feature_pos[!is.na(feature_pos)]
feature_index0 <- feature_pos - 1L

message("Reading CSR index vectors and signature non-zero values")
indptr <- rhdf5::h5read(h5_path, "/X/indptr")
indices <- rhdf5::h5read(h5_path, "/X/indices")
keep_nz <- indices %in% feature_index0
nz_position0 <- which(keep_nz) - 1L
row_index <- findInterval(nz_position0, indptr)
column_index <- match(indices[keep_nz], feature_index0)
values <- as.numeric(rhdf5::h5read(h5_path, "/X/data")[keep_nz])
signature_matrix <- Matrix::sparseMatrix(
  i = row_index, j = column_index, x = values,
  dims = c(nrow(obs), length(target_symbols)), dimnames = list(obs$cell_id, target_symbols)
)
rm(indices, keep_nz, nz_position0, row_index, column_index, values, indptr)
gc()

# UC baseline and post-treatment cells with an interpretable clinical outcome.
eligible <- obs$disease == "UC" & obs$treatment_time %in% c("Pre", "Post") &
  obs$remission %in% c("Remission", "Non_Remission") & !is.na(obs$state)
cell_meta <- obs[eligible, ]
x_uc <- signature_matrix[eligible, , drop = FALSE]
rm(signature_matrix)
gc()

# One pseudobulk per donor x cell-state x time. Groups with fewer than 20 cells
# are excluded before scoring; no cell-level hypothesis test is performed.
group_key <- interaction(cell_meta$donor, cell_meta$state, cell_meta$treatment_time, drop = TRUE, lex.order = TRUE)
group_count <- table(group_key)
keep_group <- names(group_count)[group_count >= 20L]
keep_cell <- group_key %in% keep_group
cell_meta <- cell_meta[keep_cell, ]
x_uc <- x_uc[keep_cell, , drop = FALSE]
group_key <- droplevels(group_key[keep_cell])
group_index <- as.integer(group_key)
group_n <- tabulate(group_index, nbins = nlevels(group_key))
aggregator <- Matrix::sparseMatrix(i = group_index, j = seq_along(group_index),
                                   x = 1 / group_n[group_index],
                                   dims = c(nlevels(group_key), length(group_index)))
pseudobulk_expression <- as.matrix(aggregator %*% x_uc)
rownames(pseudobulk_expression) <- levels(group_key)
colnames(pseudobulk_expression) <- target_symbols

split_cells <- split(seq_len(nrow(cell_meta)), group_key)
pseudobulk_meta <- do.call(rbind, lapply(names(split_cells), function(k) {
  z <- cell_meta[split_cells[[k]], ]
  data.frame(group_id = k, donor = z$donor[1], state = z$state[1],
             treatment_time = z$treatment_time[1], remission = z$remission[1],
             n_cells = nrow(z), n_samples = length(unique(z$sample_id)),
             mean_inflammation_score = mean(z$inflammation_score, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
rownames(pseudobulk_meta) <- NULL
pseudobulk_expression <- pseudobulk_expression[match(pseudobulk_meta$group_id, rownames(pseudobulk_expression)), , drop = FALSE]

# Standardize every signature gene using baseline UC donor-state pseudobulks and
# apply those fixed parameters to post-treatment pseudobulks.
baseline_idx <- pseudobulk_meta$treatment_time == "Pre"
gene_center <- colMeans(pseudobulk_expression[baseline_idx, , drop = FALSE])
gene_scale <- apply(pseudobulk_expression[baseline_idx, , drop = FALSE], 2, sd)
valid_gene <- is.finite(gene_scale) & gene_scale > 0
z_expression <- sweep(pseudobulk_expression[, valid_gene, drop = FALSE], 2, gene_center[valid_gene], "-")
z_expression <- sweep(z_expression, 2, gene_scale[valid_gene], "/")
pseudobulk_meta$module_score <- rowMeans(z_expression)
pseudobulk_meta$module_score_z <- as.numeric(scale(pseudobulk_meta$module_score))

# Baseline state localization with donor-based standard errors.
baseline <- pseudobulk_meta[pseudobulk_meta$treatment_time == "Pre", ]
localization <- do.call(rbind, lapply(split(baseline, baseline$state), function(x) {
  n <- nrow(x); m <- mean(x$module_score); se <- sd(x$module_score) / sqrt(n)
  data.frame(state = x$state[1], donors = n, mean_module_score = m, se = se,
             ci_lower = m - qt(.975, df = max(1, n - 1)) * se,
             ci_upper = m + qt(.975, df = max(1, n - 1)) * se,
             median_module_score = median(x$module_score))
}))
rownames(localization) <- NULL

# Baseline outcome association per state: one score per donor, Firth logistic.
outcome_rows <- lapply(split(baseline, baseline$state), function(x) {
  x$remission_binary <- as.integer(x$remission == "Remission")
  x$score_within_state_z <- as.numeric(scale(x$module_score))
  if (nrow(x) < 8L || min(table(x$remission_binary)) < 2L || !all(is.finite(x$score_within_state_z))) {
    return(data.frame(state = x$state[1], donors = nrow(x), remission_donors = sum(x$remission_binary),
                      odds_ratio = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_))
  }
  fit <- logistf::logistf(remission_binary ~ score_within_state_z, data = x)
  term <- "score_within_state_z"
  data.frame(state = x$state[1], donors = nrow(x), remission_donors = sum(x$remission_binary),
             odds_ratio = exp(fit$coefficients[term]), ci_lower = exp(fit$ci.lower[term]),
             ci_upper = exp(fit$ci.upper[term]), p_value = fit$prob[term])
})
outcome_effect <- do.call(rbind, outcome_rows)
rownames(outcome_effect) <- NULL
outcome_effect$fdr <- p.adjust(outcome_effect$p_value, method = "BH")

# Donor-paired state change. The clinical contrast is remission minus
# nonremission change; negative values mean deeper resolution in remission.
pre <- pseudobulk_meta[pseudobulk_meta$treatment_time == "Pre", ]
post <- pseudobulk_meta[pseudobulk_meta$treatment_time == "Post", ]
paired <- merge(pre[, c("donor", "state", "remission", "module_score", "n_cells")],
                post[, c("donor", "state", "remission", "module_score", "n_cells")],
                by = c("donor", "state"), suffixes = c("_pre", "_post"))
paired <- paired[paired$remission_pre == paired$remission_post, ]
paired$remission <- paired$remission_pre
paired$remission_binary <- as.integer(paired$remission == "Remission")
paired$delta <- paired$module_score_post - paired$module_score_pre

longitudinal_rows <- lapply(split(paired, paired$state), function(x) {
  if (nrow(x) < 6L || length(unique(x$remission_binary)) < 2L || min(table(x$remission_binary)) < 2L) {
    return(data.frame(state = x$state[1], paired_donors = nrow(x), remission_donors = sum(x$remission_binary),
                      mean_delta_remission = ifelse(any(x$remission_binary == 1), mean(x$delta[x$remission_binary == 1]), NA),
                      mean_delta_nonremission = ifelse(any(x$remission_binary == 0), mean(x$delta[x$remission_binary == 0]), NA),
                      estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_))
  }
  fit <- lm(delta ~ remission_binary, data = x)
  ci <- confint(fit, "remission_binary")
  data.frame(state = x$state[1], paired_donors = nrow(x), remission_donors = sum(x$remission_binary),
             mean_delta_remission = mean(x$delta[x$remission_binary == 1]),
             mean_delta_nonremission = mean(x$delta[x$remission_binary == 0]),
             estimate = coef(fit)["remission_binary"], ci_lower = ci[1], ci_upper = ci[2],
             p_value = coef(summary(fit))["remission_binary", "Pr(>|t|)"])
})
longitudinal_effect <- do.call(rbind, longitudinal_rows)
rownames(longitudinal_effect) <- NULL
longitudinal_effect$fdr <- p.adjust(longitudinal_effect$p_value, method = "BH")

# Dataset-level donor checks: compare the paper's inflammation score and our
# myeloid module score with future remission. These are audit rows, not gates.
donor_baseline <- aggregate(cbind(module_score, mean_inflammation_score) ~ donor + remission,
                            baseline, mean)
donor_baseline$remission_binary <- as.integer(donor_baseline$remission == "Remission")
audit_effect <- function(variable, label) {
  donor_baseline$z <- as.numeric(scale(donor_baseline[[variable]]))
  fit <- logistf::logistf(remission_binary ~ z, data = donor_baseline)
  data.frame(score = label, donors = nrow(donor_baseline), remission_donors = sum(donor_baseline$remission_binary),
             odds_ratio = exp(fit$coefficients["z"]), ci_lower = exp(fit$ci.lower["z"]),
             ci_upper = exp(fit$ci.upper["z"]), p_value = fit$prob["z"])
}
dataset_audit <- rbind(audit_effect("mean_inflammation_score", "Published inflammation score"),
                       audit_effect("module_score", "Locked Hallmark module (myeloid average)"))

# Cell-level values below are descriptive only and used for UMAP/dot-plot
# rendering. They do not enter any p-value.
cell_meta$cell_module_raw <- Matrix::rowMeans(x_uc)
baseline_cells <- cell_meta$treatment_time == "Pre"
cell_meta$cell_module_z <- (cell_meta$cell_module_raw - mean(cell_meta$cell_module_raw[baseline_cells])) /
  sd(cell_meta$cell_module_raw[baseline_cells])

key_genes_requested <- c("IL1B", "TNF", "IL6", "CXCL10", "CCL2", "NFKB1",
                         "STAT3", "PLAUR", "ICAM1", "FOS", "JUN", "SOCS3")
key_genes <- intersect(key_genes_requested, colnames(x_uc))
dot_rows <- lapply(split(which(baseline_cells), cell_meta$state[baseline_cells]), function(idx) {
  m <- x_uc[idx, key_genes, drop = FALSE]
  data.frame(state = cell_meta$state[idx[1]], gene = key_genes,
             mean_expression = Matrix::colMeans(m), pct_expressing = Matrix::colMeans(m > 0),
             cells = length(idx))
})
dot_data <- do.call(rbind, dot_rows)
rownames(dot_data) <- NULL
dot_data$mean_expression_z <- ave(dot_data$mean_expression, dot_data$gene,
                                  FUN = function(z) as.numeric(scale(z)))

# Baseline donor composition is also donor-aware descriptive source data.
composition_counts <- as.data.frame(table(donor = cell_meta$donor[baseline_cells],
                                          state = cell_meta$state[baseline_cells]))
composition_counts <- composition_counts[composition_counts$Freq > 0, ]
composition_counts$proportion <- ave(composition_counts$Freq, composition_counts$donor,
                                     FUN = function(z) z / sum(z))
donor_outcome <- unique(cell_meta[baseline_cells, c("donor", "remission")])
composition <- merge(composition_counts, donor_outcome, by = "donor")

utils::write.csv(pseudobulk_meta, file.path(source_dir, "singlecell_myeloid_pseudobulk.csv"), row.names = FALSE)
utils::write.csv(cell_meta[baseline_cells, c("cell_id", "donor", "sample_id", "remission", "state",
                                             "umap1", "umap2", "cell_module_z")],
                 file.path(source_dir, "singlecell_myeloid_umap.csv"), row.names = FALSE)
utils::write.csv(localization, file.path(tables_dir, "table_singlecell_myeloid_localization.csv"), row.names = FALSE)
utils::write.csv(outcome_effect, file.path(tables_dir, "table_singlecell_myeloid_outcome_effects.csv"), row.names = FALSE)
utils::write.csv(paired, file.path(source_dir, "singlecell_myeloid_paired_donors.csv"), row.names = FALSE)
utils::write.csv(longitudinal_effect, file.path(tables_dir, "table_singlecell_myeloid_longitudinal_effects.csv"), row.names = FALSE)
utils::write.csv(dataset_audit, file.path(tables_dir, "table_singlecell_dataset_level_audit.csv"), row.names = FALSE)
utils::write.csv(dot_data, file.path(source_dir, "singlecell_myeloid_key_gene_dotplot.csv"), row.names = FALSE)
utils::write.csv(composition, file.path(source_dir, "singlecell_myeloid_donor_composition.csv"), row.names = FALSE)
saveRDS(list(gene_center = gene_center, gene_scale = gene_scale, mapped_genes = target_symbols,
             valid_genes = names(gene_scale)[valid_gene]),
        file.path(derived_dir, "myeloid_module_scaling.rds"))
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "14_singlecell_myeloid_localization_sessionInfo.txt"))

print(dataset_audit)
print(localization[order(localization$mean_module_score, decreasing = TRUE), ])
print(outcome_effect[order(outcome_effect$p_value), ])
print(longitudinal_effect[order(longitudinal_effect$p_value), ])
message("DONOR-AWARE MYELOID LOCALIZATION COMPLETE: ", length(unique(baseline$donor)),
        " baseline UC donors, ", sum(baseline$n_cells), " baseline myeloid cells, ",
        length(target_symbols), " mapped signature genes")
