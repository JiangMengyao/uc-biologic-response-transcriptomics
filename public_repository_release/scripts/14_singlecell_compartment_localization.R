#!/usr/bin/env Rscript

# Sample-first, donor-inference localization of the locked inflammatory-response
# program across GSE282122 colonic epithelium, myeloid and fibroblast/pericyte
# compartments. H5AD CSR matrices are streamed in bounded chunks. The clinical
# audit is deliberately separated from the localization result.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
raw_dir <- file.path(project_root, "data", "raw", "GSE282122")
derived_dir <- file.path(project_root, "data", "derived", "singlecell")
tables_dir <- file.path(project_root, "results", "tables", "common_state")
source_dir <- file.path(project_root, "results", "source_data")
logs_dir <- file.path(project_root, "logs")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(logistf)
})
set.seed(20260720)

algorithm_version <- "cross-compartment-sample-first-csr-v2"
minimum_cells <- 20L
minimum_baseline_donors <- 8L
minimum_paired_donors <- 6L
locked_signature_md5 <- "e072b3b52623f871f902da42aaed04c9"
chunk_nnz <- as.integer(Sys.getenv("H5AD_CHUNK_NNZ", "5000000"))
if (!is.finite(chunk_nnz) || chunk_nnz < 100000L) {
  stop("H5AD_CHUNK_NNZ must be at least 100000")
}

compartment_files <- c(
  "Colonic epithelium" = "epicolonic_final.h5ad",
  "Myeloid" = "myeloid_final.h5ad",
  "Fibroblast/pericyte" = "fibperi_final.h5ad"
)
h5_paths <- file.path(raw_dir, unname(compartment_files))
names(h5_paths) <- names(compartment_files)
missing_h5 <- h5_paths[!file.exists(h5_paths) | file.info(h5_paths)$size <= 0]
if (length(missing_h5)) {
  stop("SINGLE-CELL INPUT GATE FAILED: missing ", paste(basename(missing_h5), collapse = ", "),
       ". Run scripts/13_download_singlecell_subset.sh first.")
}

signature_file <- file.path(project_root, "config", "predefined_gene_sets.csv")
signature_md5 <- unname(tools::md5sum(signature_file))
if (!identical(signature_md5, locked_signature_md5)) {
  stop("FROZEN SIGNATURE GATE FAILED: observed MD5 ", signature_md5,
       "; locked amendment MD5 ", locked_signature_md5)
}

read_categorical <- function(path, name) {
  base <- paste0("/obs/", name)
  categories <- rhdf5::h5read(path, paste0(base, "/categories"))
  codes <- rhdf5::h5read(path, paste0(base, "/codes"))
  out <- rep(NA_character_, length(codes))
  ok <- !is.na(codes) & codes >= 0L
  out[ok] <- as.character(categories[codes[ok] + 1L])
  out
}

read_obs <- function(path) {
  out <- data.frame(
    cell_id = as.character(rhdf5::h5read(path, "/obs/_index")),
    donor = read_categorical(path, "Patient"),
    sample_id = read_categorical(path, "sample_id"),
    disease = read_categorical(path, "Disease"),
    treatment_time = read_categorical(path, "Treatment"),
    remission = read_categorical(path, "Remission_status"),
    inflammation = read_categorical(path, "Inflammation"),
    inflammation_score = as.numeric(rhdf5::h5read(path, "/obs/Inflammation_score")),
    site = read_categorical(path, "Site"),
    match = read_categorical(path, "Match"),
    state = read_categorical(path, "final_analysis"),
    stringsAsFactors = FALSE
  )
  required <- c("donor", "sample_id", "disease", "treatment_time", "remission",
                "inflammation", "site", "match", "state")
  if (any(vapply(out[required], length, integer(1)) != nrow(out))) {
    stop("H5AD METADATA GATE FAILED: obs vectors have inconsistent lengths in ", basename(path))
  }
  out
}

assert_sample_metadata <- function(obs, eligible, compartment) {
  fields <- c("sample_id", "donor", "disease", "treatment_time", "remission",
              "inflammation", "site", "match")
  sample_meta <- unique(obs[eligible, fields, drop = FALSE])
  conflict <- duplicated(sample_meta$sample_id) | duplicated(sample_meta$sample_id, fromLast = TRUE)
  if (any(conflict)) {
    bad <- unique(sample_meta$sample_id[conflict])
    stop("SAMPLE METADATA GATE FAILED in ", compartment, ": conflicting labels for ",
         paste(head(bad, 8L), collapse = ", "))
  }
}

sets <- readRDS(file.path(project_root, "data", "derived", "predefined_gene_sets.rds"))
signature <- unique(as.character(sets$Inflammation))
if (length(signature) != 200L) {
  stop("FROZEN SIGNATURE SIZE GATE FAILED: expected 200 genes, observed ", length(signature))
}
feature_symbols <- lapply(h5_paths, function(path) {
  as.character(rhdf5::h5read(path, "/var/gene_symbol"))
})
common_symbols <- Reduce(intersect, lapply(feature_symbols, unique))
target_symbols <- signature[signature %in% common_symbols]
if (length(target_symbols) < 100L) {
  stop("SINGLE-CELL MAPPING GATE FAILED: only ", length(target_symbols),
       " of 200 locked genes are common to all three compartments")
}

make_group_definition <- function(obs, eligible, keys, prefix, retain_fraction = FALSE) {
  key_ok <- Reduce(`&`, lapply(obs[, keys, drop = FALSE], function(z) {
    !is.na(z) & nzchar(as.character(z))
  }))
  keep <- eligible & key_ok
  group_key_all <- do.call(paste, c(obs[keys], sep = "|||"))
  counts <- table(group_key_all[keep])
  group_levels <- sort(names(counts)[counts >= minimum_cells])
  cell_to_group <- match(group_key_all, group_levels)
  cell_to_group[!keep] <- NA_integer_
  retained <- !is.na(cell_to_group)
  if (!length(group_levels) || !any(retained)) {
    stop("PSEUDOBULK GATE FAILED: no ", prefix, " groups with >= ", minimum_cells, " cells")
  }
  first_row <- match(group_levels, group_key_all)
  group_cells <- tabulate(cell_to_group[retained], nbins = length(group_levels))
  sample_split <- split(obs$sample_id[retained], cell_to_group[retained])
  mean_inflammation <- as.numeric(tapply(obs$inflammation_score[retained],
                                         cell_to_group[retained], mean, na.rm = TRUE))
  meta <- obs[first_row, c("donor", "sample_id", "treatment_time", "remission",
                           "inflammation", "site", "match", "state"), drop = FALSE]
  meta$group_id <- paste(prefix, group_levels, sep = "|||")
  meta$n_cells <- group_cells
  meta$n_samples <- vapply(seq_along(group_levels), function(i) {
    length(unique(sample_split[[as.character(i)]]))
  }, integer(1))
  meta$mean_inflammation_score <- mean_inflammation
  meta <- meta[, c("group_id", "donor", "sample_id", "state", "treatment_time",
                   "remission", "inflammation", "site", "match", "n_cells",
                   "n_samples", "mean_inflammation_score")]
  rownames(meta) <- meta$group_id
  list(meta = meta, map = cell_to_group, retain_fraction = retain_fraction)
}

process_compartment <- function(compartment, path, symbols) {
  slug <- gsub("[^a-z0-9]+", "_", tolower(compartment))
  cache_path <- file.path(derived_dir, paste0(slug, "_sample_first_signature_pseudobulk.rds"))
  input_md5 <- unname(tools::md5sum(path))
  if (file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    cache_ok <- identical(cached$algorithm_version, algorithm_version) &&
      identical(cached$input_md5, input_md5) &&
      identical(cached$signature_md5, signature_md5) &&
      identical(cached$target_symbols, target_symbols)
    if (cache_ok) {
      message("CACHE PASS: ", compartment, " sample-first pseudobulks")
      return(cached)
    }
  }

  message("READ OBS: ", compartment)
  obs <- read_obs(path)
  eligible <- obs$disease == "UC" & obs$treatment_time %in% c("Pre", "Post") &
    obs$remission %in% c("Remission", "Non_Remission") &
    !is.na(obs$donor) & nzchar(obs$donor) &
    !is.na(obs$sample_id) & nzchar(obs$sample_id) &
    !is.na(obs$site) & nzchar(obs$site) &
    !is.na(obs$match) & nzchar(obs$match)
  eligible[is.na(eligible)] <- FALSE
  if (!any(eligible)) stop("ELIGIBILITY GATE FAILED: no eligible UC cells in ", compartment)
  assert_sample_metadata(obs, eligible, compartment)

  state_eligible <- eligible & !is.na(obs$state) & nzchar(obs$state)
  definitions <- list(
    sample_state = make_group_definition(obs, state_eligible, c("sample_id", "state"),
                                         paste0(slug, "|||sample_state"), TRUE),
    sample_compartment = make_group_definition(obs, eligible, "sample_id",
                                               paste0(slug, "|||sample"), FALSE),
    legacy_state = make_group_definition(obs, state_eligible,
                                         c("donor", "state", "treatment_time"),
                                         paste0(slug, "|||legacy"), FALSE)
  )
  for (kind in names(definitions)) {
    definitions[[kind]]$meta$compartment <- compartment
    definitions[[kind]]$meta <- definitions[[kind]]$meta[, c(
      "group_id", "compartment", "donor", "sample_id", "state", "treatment_time",
      "remission", "inflammation", "site", "match", "n_cells", "n_samples",
      "mean_inflammation_score")]
    nr <- nrow(definitions[[kind]]$meta)
    definitions[[kind]]$sums <- matrix(0, nr, length(target_symbols),
                                       dimnames = list(definitions[[kind]]$meta$group_id, target_symbols))
    if (definitions[[kind]]$retain_fraction) {
      definitions[[kind]]$positive <- matrix(0, nr, length(target_symbols),
                                             dimnames = list(definitions[[kind]]$meta$group_id,
                                                             target_symbols))
    }
  }

  feature_pos <- match(target_symbols, symbols)
  if (anyNA(feature_pos)) stop("INTERNAL MAPPING ERROR in ", compartment)
  feature_lookup <- integer(length(symbols))
  feature_lookup[feature_pos] <- seq_along(feature_pos)
  indptr <- as.numeric(rhdf5::h5read(path, "/X/indptr"))
  if (length(indptr) != nrow(obs) + 1L) stop("CSR SHAPE GATE FAILED in ", compartment)
  nnz <- tail(indptr, 1)
  starts <- seq.int(1, nnz, by = chunk_nnz)
  message("STREAM CSR: ", compartment, "; ", format(nnz, big.mark = ","),
          " non-zero values in ", length(starts), " chunks")

  for (chunk_i in seq_along(starts)) {
    lo <- starts[chunk_i]
    hi <- min(nnz, lo + chunk_nnz - 1)
    count <- hi - lo + 1
    indices0 <- as.integer(rhdf5::h5read(path, "/X/indices", start = lo, count = count))
    gene_j <- feature_lookup[indices0 + 1L]
    target_nz <- which(gene_j > 0L)
    if (length(target_nz)) {
      values <- as.numeric(rhdf5::h5read(path, "/X/data", start = lo, count = count))
      positions0 <- (lo - 1) + (target_nz - 1)
      cell_i <- findInterval(positions0, indptr)
      gj_all <- gene_j[target_nz]
      gx_all <- values[target_nz]
      for (kind in names(definitions)) {
        group_i <- definitions[[kind]]$map[cell_i]
        keep <- !is.na(group_i)
        if (!any(keep)) next
        gi <- group_i[keep]
        gj <- gj_all[keep]
        gx <- gx_all[keep]
        dims <- c(nrow(definitions[[kind]]$meta), length(target_symbols))
        definitions[[kind]]$sums <- definitions[[kind]]$sums + as.matrix(
          Matrix::sparseMatrix(i = gi, j = gj, x = gx, dims = dims)
        )
        if (definitions[[kind]]$retain_fraction && any(gx > 0)) {
          positive <- gx > 0
          definitions[[kind]]$positive <- definitions[[kind]]$positive + as.matrix(
            Matrix::sparseMatrix(i = gi[positive], j = gj[positive],
                                 x = rep(1, sum(positive)), dims = dims)
          )
        }
      }
    }
    if (chunk_i %% 10L == 0L || chunk_i == length(starts)) {
      message("  ", compartment, ": chunk ", chunk_i, "/", length(starts))
    }
  }

  result_groups <- lapply(definitions, function(definition) {
    expression <- sweep(definition$sums, 1L, definition$meta$n_cells, "/")
    fraction <- if (definition$retain_fraction) {
      sweep(definition$positive, 1L, definition$meta$n_cells, "/")
    } else NULL
    list(meta = definition$meta, expression = expression,
         fraction_expressing = fraction)
  })
  baseline_donors <- length(unique(result_groups$sample_compartment$meta$donor[
    result_groups$sample_compartment$meta$treatment_time == "Pre" &
      result_groups$sample_compartment$meta$inflammation == "Inflamed"]))
  if (baseline_donors < minimum_baseline_donors) {
    stop("DONOR GATE FAILED: ", compartment, " has only ", baseline_donors,
         " retained baseline inflamed UC donors")
  }
  result <- c(list(
    algorithm_version = algorithm_version,
    input_md5 = input_md5,
    signature_md5 = signature_md5,
    target_symbols = target_symbols,
    audit = data.frame(compartment = compartment, input_cells = nrow(obs),
                       eligible_cells = sum(eligible),
                       eligible_samples = length(unique(obs$sample_id[eligible])),
                       sample_state_groups = nrow(result_groups$sample_state$meta),
                       sample_compartment_groups = nrow(result_groups$sample_compartment$meta),
                       legacy_state_groups = nrow(result_groups$legacy_state$meta),
                       stringsAsFactors = FALSE)
  ), result_groups)
  saveRDS(result, cache_path, compress = "gzip")
  message("CACHE WRITE: ", basename(cache_path), "; ",
          nrow(result$sample_state$meta), " sample-state and ",
          nrow(result$sample_compartment$meta), " sample-compartment pseudobulks")
  result
}

objects <- Map(process_compartment, names(h5_paths), h5_paths, feature_symbols)

combine_kind <- function(objects, kind, field) {
  values <- lapply(objects, function(x) x[[kind]][[field]])
  if (field == "meta") {
    out <- do.call(rbind, values)
    rownames(out) <- out$group_id
    return(out)
  }
  out <- do.call(rbind, values)
  out
}

sample_state_meta <- combine_kind(objects, "sample_state", "meta")
sample_state_expression <- combine_kind(objects, "sample_state", "expression")
sample_state_fraction <- combine_kind(objects, "sample_state", "fraction_expressing")
sample_comp_meta <- combine_kind(objects, "sample_compartment", "meta")
sample_comp_expression <- combine_kind(objects, "sample_compartment", "expression")
legacy_meta <- combine_kind(objects, "legacy_state", "meta")
legacy_expression <- combine_kind(objects, "legacy_state", "expression")
sample_state_expression <- sample_state_expression[sample_state_meta$group_id, , drop = FALSE]
sample_state_fraction <- sample_state_fraction[sample_state_meta$group_id, , drop = FALSE]
sample_comp_expression <- sample_comp_expression[sample_comp_meta$group_id, , drop = FALSE]
legacy_expression <- legacy_expression[legacy_meta$group_id, , drop = FALSE]

# Metadata must also agree when the same biological sample occurs in different
# compartment files.
cross_fields <- c("sample_id", "donor", "treatment_time", "remission",
                  "inflammation", "site", "match")
cross_sample <- unique(sample_comp_meta[, cross_fields])
cross_conflict <- duplicated(cross_sample$sample_id) |
  duplicated(cross_sample$sample_id, fromLast = TRUE)
if (any(cross_conflict)) {
  stop("CROSS-COMPARTMENT SAMPLE METADATA GATE FAILED for: ",
       paste(unique(cross_sample$sample_id[cross_conflict]), collapse = ", "))
}

collapse_expression <- function(meta, expression, keys, filter, prefix) {
  idx <- which(filter & stats::complete.cases(meta[, keys, drop = FALSE]))
  if (!length(idx)) stop("COLLAPSE GATE FAILED: no rows for ", prefix)
  group_key <- do.call(paste, c(meta[idx, keys, drop = FALSE], sep = "|||"))
  groups <- split(idx, factor(group_key, levels = unique(group_key)))
  rows <- vector("list", length(groups))
  matrices <- vector("list", length(groups))
  for (i in seq_along(groups)) {
    ii <- groups[[i]]
    z <- meta[ii, , drop = FALSE]
    row <- z[1, , drop = FALSE]
    row$group_id <- paste(prefix, names(groups)[i], sep = "|||")
    row$n_cells <- sum(z$n_cells)
    row$n_samples <- length(unique(z$sample_id))
    row$sample_id <- if (row$n_samples == 1L) unique(z$sample_id) else "multiple_equal_weight"
    row$site <- if (length(unique(z$site)) == 1L) unique(z$site) else "Multiple_sites"
    row$inflammation <- if (length(unique(z$inflammation)) == 1L) unique(z$inflammation) else "Mixed"
    row$match <- if (length(unique(z$match)) == 1L) unique(z$match) else "Mixed"
    row$mean_inflammation_score <- mean(z$mean_inflammation_score, na.rm = TRUE)
    rows[[i]] <- row
    matrices[[i]] <- colMeans(expression[ii, , drop = FALSE])
  }
  out_meta <- do.call(rbind, rows)
  rownames(out_meta) <- out_meta$group_id
  out_expression <- do.call(rbind, matrices)
  rownames(out_expression) <- out_meta$group_id
  list(meta = out_meta, expression = out_expression)
}

localization_filter <- sample_state_meta$treatment_time == "Pre" &
  sample_state_meta$inflammation == "Inflamed"
donor_state_primary <- collapse_expression(
  sample_state_meta, sample_state_expression,
  c("compartment", "donor", "state"), localization_filter, "sample_first_donor_state"
)

# Frozen cross-compartment scaling is learned only from the sample-balanced,
# baseline-inflamed donor-state reference and then applied unchanged elsewhere.
gene_center <- colMeans(donor_state_primary$expression)
gene_scale <- apply(donor_state_primary$expression, 2L, stats::sd)
valid_gene <- is.finite(gene_scale) & gene_scale > 0
if (sum(valid_gene) < 100L) {
  stop("SCALING GATE FAILED: fewer than 100 variable common signature genes")
}
score_matrix <- function(expression) {
  z <- sweep(expression[, valid_gene, drop = FALSE], 2L, gene_center[valid_gene], "-")
  z <- sweep(z, 2L, gene_scale[valid_gene], "/")
  rowMeans(z)
}
sample_state_meta$module_score <- score_matrix(sample_state_expression)
sample_comp_meta$module_score <- score_matrix(sample_comp_expression)
donor_state_primary$meta$module_score <- score_matrix(donor_state_primary$expression)

summarize_localization <- function(x) {
  n <- nrow(x)
  m <- mean(x$module_score)
  se <- if (n >= 2L) stats::sd(x$module_score) / sqrt(n) else NA_real_
  critical <- if (n >= 3L) stats::qt(.975, df = n - 1L) else NA_real_
  data.frame(
    compartment = x$compartment[1], state = x$state[1], donors = n,
    mean_module_score = m, median_module_score = stats::median(x$module_score),
    se = se, ci_lower = m - critical * se, ci_upper = m + critical * se,
    total_cells = sum(x$n_cells), total_samples = sum(x$n_samples),
    aggregation = "Sample-first; donor-equal", baseline_stratum = "Inflamed",
    stringsAsFactors = FALSE
  )
}
localization <- do.call(rbind, lapply(
  split(donor_state_primary$meta,
        interaction(donor_state_primary$meta$compartment,
                    donor_state_primary$meta$state, drop = TRUE)),
  summarize_localization
))
rownames(localization) <- NULL

collapse_scores <- function(x, keys, filter, prefix) {
  idx <- which(filter & stats::complete.cases(x[, keys, drop = FALSE]))
  if (!length(idx)) return(x[FALSE, , drop = FALSE])
  group_key <- do.call(paste, c(x[idx, keys, drop = FALSE], sep = "|||"))
  groups <- split(idx, factor(group_key, levels = unique(group_key)))
  out <- lapply(seq_along(groups), function(i) {
    ii <- groups[[i]]
    z <- x[ii, , drop = FALSE]
    row <- z[1, , drop = FALSE]
    row$group_id <- paste(prefix, names(groups)[i], sep = "|||")
    row$module_score <- mean(z$module_score)
    row$n_cells <- sum(z$n_cells)
    row$n_samples <- length(unique(z$sample_id))
    row$sample_id <- if (row$n_samples == 1L) unique(z$sample_id) else "multiple_equal_weight"
    row$site <- if (length(unique(z$site)) == 1L) unique(z$site) else "Multiple_sites"
    row$inflammation <- if (length(unique(z$inflammation)) == 1L) unique(z$inflammation) else "Mixed"
    row$match <- if (length(unique(z$match)) == 1L) unique(z$match) else "Mixed"
    row$mean_inflammation_score <- mean(z$mean_inflammation_score, na.rm = TRUE)
    row
  })
  result <- do.call(rbind, out)
  rownames(result) <- NULL
  result
}

baseline_inflamed_filter <- sample_comp_meta$treatment_time == "Pre" &
  sample_comp_meta$inflammation == "Inflamed"
donor_compartment_primary <- collapse_scores(
  sample_comp_meta, c("compartment", "donor"), baseline_inflamed_filter,
  "sample_first_baseline_inflamed"
)
baseline_counts <- tapply(donor_compartment_primary$donor,
                          donor_compartment_primary$compartment,
                          function(z) length(unique(z)))
if (any(baseline_counts[names(compartment_files)] < minimum_baseline_donors)) {
  stop("BASELINE DONOR GATE FAILED: ",
       paste(names(baseline_counts), baseline_counts, collapse = "; "))
}

firth_effect <- function(x, labels) {
  result <- as.data.frame(labels, stringsAsFactors = FALSE)
  x$remission_binary <- as.integer(x$remission == "Remission")
  x$score_z <- as.numeric(scale(x$module_score))
  result$donors <- nrow(x)
  result$remission_donors <- sum(x$remission_binary)
  result$nonremission_donors <- sum(x$remission_binary == 0L)
  result$odds_ratio <- result$ci_lower <- result$ci_upper <- result$p_value <- NA_real_
  if (nrow(x) < minimum_baseline_donors || length(unique(x$remission_binary)) < 2L ||
      min(table(x$remission_binary)) < 2L || !all(is.finite(x$score_z))) return(result)
  fit <- logistf::logistf(remission_binary ~ score_z, data = x)
  result$odds_ratio <- exp(fit$coefficients["score_z"])
  result$ci_lower <- exp(fit$ci.lower["score_z"])
  result$ci_upper <- exp(fit$ci.upper["score_z"])
  result$p_value <- fit$prob["score_z"]
  result
}

outcome_compartment <- do.call(rbind, lapply(
  split(donor_compartment_primary, donor_compartment_primary$compartment),
  function(x) firth_effect(x, list(compartment = x$compartment[1]))
))
rownames(outcome_compartment) <- NULL
outcome_compartment$aggregation <- "Sample-first; donor-equal"
outcome_compartment$stratum <- "Baseline inflamed"
outcome_compartment$fdr <- p.adjust(outcome_compartment$p_value, method = "BH")
outcome_compartment$direction_vs_bulk <- ifelse(
  outcome_compartment$odds_ratio < 1, "Concordant", "Discordant"
)

outcome_state <- do.call(rbind, lapply(
  split(donor_state_primary$meta,
        interaction(donor_state_primary$meta$compartment,
                    donor_state_primary$meta$state, drop = TRUE)),
  function(x) firth_effect(x, list(compartment = x$compartment[1], state = x$state[1]))
))
rownames(outcome_state) <- NULL
outcome_state$aggregation <- "Sample-first; donor-equal"
outcome_state$stratum <- "Baseline inflamed"
outcome_state$fdr <- ave(outcome_state$p_value, outcome_state$compartment,
                         FUN = function(z) p.adjust(z, method = "BH"))

# Prespecified sample-first sensitivity strata.
strata <- list(
  "Sample-first | Inflamed only" = sample_comp_meta$treatment_time == "Pre" &
    sample_comp_meta$inflammation == "Inflamed",
  "Sample-first | All baseline samples" = sample_comp_meta$treatment_time == "Pre",
  "Sample-first | Non-inflamed only" = sample_comp_meta$treatment_time == "Pre" &
    sample_comp_meta$inflammation == "Non_Inflamed",
  "Sample-first | Inflamed rectum" = sample_comp_meta$treatment_time == "Pre" &
    sample_comp_meta$inflammation == "Inflamed" & sample_comp_meta$site == "Rectum",
  "Sample-first | Inflamed distal colon" = sample_comp_meta$treatment_time == "Pre" &
    sample_comp_meta$inflammation == "Inflamed" &
    sample_comp_meta$site %in% c("Descending_Colon", "Sigmoid", "Rectum")
)
outcome_sensitivity <- do.call(rbind, lapply(names(strata), function(label) {
  collapsed <- collapse_scores(sample_comp_meta, c("compartment", "donor"),
                               strata[[label]], paste0("sensitivity_", gsub("[^a-z]+", "_", tolower(label))))
  result <- do.call(rbind, lapply(split(collapsed, collapsed$compartment), function(x) {
    firth_effect(x, list(compartment = x$compartment[1]))
  }))
  result$analysis_label <- label
  result$aggregation <- "Sample-first; donor-equal"
  result$stratum <- sub("^Sample-first \\| ", "", label)
  result$primary <- label == "Sample-first | Inflamed only"
  result
}))

# Exact reconstruction of the pre-amendment donor-state/cell-weighted clinical
# result. It remains supplementary and is never used to alter the signature.
legacy_baseline <- legacy_meta$treatment_time == "Pre"
legacy_center <- colMeans(legacy_expression[legacy_baseline, , drop = FALSE])
legacy_scale <- apply(legacy_expression[legacy_baseline, , drop = FALSE], 2L, stats::sd)
legacy_valid <- is.finite(legacy_scale) & legacy_scale > 0
legacy_z <- sweep(legacy_expression[, legacy_valid, drop = FALSE], 2L,
                  legacy_center[legacy_valid], "-")
legacy_z <- sweep(legacy_z, 2L, legacy_scale[legacy_valid], "/")
legacy_meta$module_score <- rowMeans(legacy_z)
legacy_groups <- split(seq_len(nrow(legacy_meta)),
                       interaction(legacy_meta$compartment, legacy_meta$donor,
                                   legacy_meta$treatment_time, drop = TRUE))
legacy_donor <- do.call(rbind, lapply(legacy_groups, function(ii) {
  z <- legacy_meta[ii, , drop = FALSE]
  data.frame(compartment = z$compartment[1], donor = z$donor[1],
             treatment_time = z$treatment_time[1], remission = z$remission[1],
             module_score = weighted.mean(z$module_score, z$n_cells),
             n_cells = sum(z$n_cells), n_states = nrow(z), stringsAsFactors = FALSE)
}))
rownames(legacy_donor) <- NULL
legacy_outcome <- do.call(rbind, lapply(
  split(legacy_donor[legacy_donor$treatment_time == "Pre", ],
        legacy_donor$compartment[legacy_donor$treatment_time == "Pre"]),
  function(x) firth_effect(x, list(compartment = x$compartment[1]))
))
legacy_outcome$analysis_label <- "Legacy state-cell-weighted"
legacy_outcome$aggregation <- "Donor-state first; cell-weighted donor summary"
legacy_outcome$stratum <- "All baseline cells; sites and inflammation mixed"
legacy_outcome$primary <- FALSE
legacy_epithelial <- legacy_outcome$odds_ratio[
  legacy_outcome$compartment == "Colonic epithelium"]
if (length(legacy_epithelial) != 1L || !is.finite(legacy_epithelial) ||
    abs(log(legacy_epithelial / 2.8969)) > 0.05) {
  stop("LEGACY RECONSTRUCTION GATE FAILED: epithelial OR=", legacy_epithelial,
       "; expected approximately 2.90")
}
outcome_sensitivity <- rbind(outcome_sensitivity, legacy_outcome)
rownames(outcome_sensitivity) <- NULL
outcome_sensitivity$fdr <- ave(outcome_sensitivity$p_value,
                               outcome_sensitivity$analysis_label,
                               FUN = function(z) p.adjust(z, method = "BH"))
outcome_sensitivity$direction_vs_bulk <- ifelse(
  outcome_sensitivity$odds_ratio < 1, "Concordant",
  ifelse(outcome_sensitivity$odds_ratio > 1, "Discordant", NA_character_)
)
outcome_sensitivity$interpretation <- ifelse(
  outcome_sensitivity$analysis_label == "Legacy state-cell-weighted" &
    outcome_sensitivity$compartment == "Colonic epithelium",
  "Direction inconsistent with bulk; nominal P<0.05 but FDR>0.05; aggregation-sensitive",
  "Exploratory single-cell clinical audit"
)

paired_effect <- function(x, labels) {
  result <- as.data.frame(labels, stringsAsFactors = FALSE)
  result$paired_donors <- nrow(x)
  result$remission_donors <- sum(x$remission_binary)
  result$nonremission_donors <- sum(x$remission_binary == 0L)
  result$mean_delta_remission <- if (any(x$remission_binary == 1L)) {
    mean(x$delta[x$remission_binary == 1L])
  } else NA_real_
  result$mean_delta_nonremission <- if (any(x$remission_binary == 0L)) {
    mean(x$delta[x$remission_binary == 0L])
  } else NA_real_
  result$estimate <- result$ci_lower <- result$ci_upper <- result$p_value <- NA_real_
  if (nrow(x) < minimum_paired_donors || length(unique(x$remission_binary)) < 2L ||
      min(table(x$remission_binary)) < 2L) return(result)
  fit <- stats::lm(delta ~ remission_binary, data = x)
  ci <- stats::confint(fit, "remission_binary")
  result$estimate <- unname(stats::coef(fit)["remission_binary"])
  result$ci_lower <- ci[1]
  result$ci_upper <- ci[2]
  result$p_value <- stats::coef(summary(fit))["remission_binary", "Pr(>|t|)"]
  result
}

make_exact_pairs <- function(x, state_level = FALSE) {
  keys <- c("compartment", "donor", "site")
  if (state_level) keys <- c(keys, "state")
  matched <- x$match == "Yes" & x$treatment_time %in% c("Pre", "Post")
  site_time <- collapse_scores(x, c(keys, "treatment_time"), matched,
                               if (state_level) "matched_site_time_state" else "matched_site_time")
  pre <- site_time[site_time$treatment_time == "Pre",
                   c(keys, "remission", "module_score", "n_cells", "n_samples")]
  post <- site_time[site_time$treatment_time == "Post",
                    c(keys, "remission", "module_score", "n_cells", "n_samples")]
  if (anyDuplicated(pre[, keys, drop = FALSE]) || anyDuplicated(post[, keys, drop = FALSE])) {
    stop("LONGITUDINAL PAIRING GATE FAILED: duplicate donor-site-time records")
  }
  paired <- merge(pre, post, by = keys, suffixes = c("_pre", "_post"))
  paired <- paired[paired$remission_pre == paired$remission_post, , drop = FALSE]
  paired$remission <- paired$remission_pre
  paired$remission_binary <- as.integer(paired$remission == "Remission")
  paired$delta <- paired$module_score_post - paired$module_score_pre
  paired$pairing_rule <- if (state_level) {
    "Match=Yes; exact donor-site-state"
  } else "Match=Yes; exact donor-site"
  paired
}

paired_site_compartment <- make_exact_pairs(sample_comp_meta, FALSE)
paired_donor_groups <- split(seq_len(nrow(paired_site_compartment)),
                             interaction(paired_site_compartment$compartment,
                                         paired_site_compartment$donor, drop = TRUE))
paired_compartment <- do.call(rbind, lapply(paired_donor_groups, function(ii) {
  z <- paired_site_compartment[ii, , drop = FALSE]
  data.frame(compartment = z$compartment[1], donor = z$donor[1],
             remission = z$remission[1],
             remission_binary = z$remission_binary[1], delta = mean(z$delta),
             mean_pre = mean(z$module_score_pre), mean_post = mean(z$module_score_post),
             n_sites = nrow(z), sites = paste(sort(unique(z$site)), collapse = ";"),
             pairing_rule = "Match=Yes; exact donor-site then donor-equal",
             stringsAsFactors = FALSE)
}))
rownames(paired_compartment) <- NULL
longitudinal_compartment <- do.call(rbind, lapply(
  split(paired_compartment, paired_compartment$compartment),
  function(x) paired_effect(x, list(compartment = x$compartment[1]))
))
rownames(longitudinal_compartment) <- NULL
longitudinal_compartment$aggregation <- "Sample-first; exact donor-site; donor-equal"
longitudinal_compartment$match_required <- "Yes"
longitudinal_compartment$fdr <- p.adjust(longitudinal_compartment$p_value, method = "BH")

pairing_audit <- do.call(rbind, lapply(names(compartment_files), function(compartment) {
  sc <- sample_comp_meta[sample_comp_meta$compartment == compartment &
                           sample_comp_meta$match == "Yes", ]
  ps <- paired_site_compartment[paired_site_compartment$compartment == compartment, ]
  pd <- paired_compartment[paired_compartment$compartment == compartment, ]
  data.frame(compartment = compartment,
             match_yes_pre_samples = length(unique(sc$sample_id[sc$treatment_time == "Pre"])),
             match_yes_post_samples = length(unique(sc$sample_id[sc$treatment_time == "Post"])),
             paired_donor_sites = nrow(ps), paired_donors = nrow(pd),
             remission_donors = sum(pd$remission == "Remission"),
             nonremission_donors = sum(pd$remission == "Non_Remission"),
             exact_pairing = TRUE, stringsAsFactors = FALSE)
}))
if (any(pairing_audit$paired_donors < minimum_paired_donors) ||
    any(pairing_audit$remission_donors < 2L) ||
    any(pairing_audit$nonremission_donors < 2L)) {
  stop("MATCH=YES LONGITUDINAL DONOR GATE FAILED: ",
       paste(pairing_audit$compartment, pairing_audit$paired_donors,
             pairing_audit$remission_donors, pairing_audit$nonremission_donors,
             sep = "/", collapse = "; "))
}
if (any(!is.finite(longitudinal_compartment$estimate)) ||
    any(longitudinal_compartment$estimate >= 0)) {
  stop("EXPECTED-DIRECTION GATE FAILED: Match=Yes remission-minus-nonremission ",
       "change was not negative in every compartment: ",
       paste(longitudinal_compartment$compartment,
             signif(longitudinal_compartment$estimate, 4), collapse = "; "))
}

paired_site_state <- make_exact_pairs(sample_state_meta, TRUE)
paired_state_groups <- split(seq_len(nrow(paired_site_state)),
                             interaction(paired_site_state$compartment,
                                         paired_site_state$donor,
                                         paired_site_state$state, drop = TRUE))
paired_state <- do.call(rbind, lapply(paired_state_groups, function(ii) {
  z <- paired_site_state[ii, , drop = FALSE]
  data.frame(compartment = z$compartment[1], donor = z$donor[1], state = z$state[1],
             remission = z$remission[1], remission_binary = z$remission_binary[1],
             delta = mean(z$delta), n_sites = nrow(z), stringsAsFactors = FALSE)
}))
longitudinal_state <- do.call(rbind, lapply(
  split(paired_state, interaction(paired_state$compartment, paired_state$state, drop = TRUE)),
  function(x) paired_effect(x, list(compartment = x$compartment[1], state = x$state[1]))
))
rownames(longitudinal_state) <- NULL
longitudinal_state$aggregation <- "Sample-first; Match=Yes exact donor-site-state"
longitudinal_state$fdr <- ave(longitudinal_state$p_value, longitudinal_state$compartment,
                              FUN = function(z) p.adjust(z, method = "BH"))

# Descriptive sample-balanced key-gene map. No cell-level p-value is calculated.
key_requested <- c("IL1B", "TNF", "IL6", "CXCL10", "CCL2", "NFKB1", "STAT3",
                   "PLAUR", "ICAM1", "FOS", "JUN", "SOCS3", "CXCL8", "OSM", "JAK2")
key_genes <- intersect(key_requested, colnames(sample_state_expression))
dot_idx <- which(localization_filter)
dot_groups <- split(dot_idx,
                    interaction(sample_state_meta$compartment[dot_idx],
                                sample_state_meta$state[dot_idx], drop = TRUE))
dot <- do.call(rbind, lapply(dot_groups, function(ii) {
  data.frame(
    compartment = sample_state_meta$compartment[ii[1]],
    state = sample_state_meta$state[ii[1]], gene = key_genes,
    donors = length(unique(sample_state_meta$donor[ii])), samples = length(ii),
    cells = sum(sample_state_meta$n_cells[ii]),
    mean_expression = colMeans(sample_state_expression[ii, key_genes, drop = FALSE]),
    pct_expressing = colMeans(sample_state_fraction[ii, key_genes, drop = FALSE]),
    aggregation = "Sample-first; samples equal-weighted",
    stringsAsFactors = FALSE
  )
}))
rownames(dot) <- NULL
dot$mean_expression_z <- ave(dot$mean_expression, dot$gene, FUN = function(z) {
  s <- stats::sd(z)
  if (!is.finite(s) || s == 0) rep(0, length(z)) else (z - mean(z)) / s
})

sample_state_export <- sample_state_meta[, c(
  "group_id", "compartment", "donor", "sample_id", "state", "treatment_time",
  "remission", "inflammation", "site", "match", "n_cells", "n_samples",
  "mean_inflammation_score", "module_score")]
sample_comp_export <- sample_comp_meta[, c(
  "group_id", "compartment", "donor", "sample_id", "treatment_time", "remission",
  "inflammation", "site", "match", "n_cells", "n_samples",
  "mean_inflammation_score", "module_score")]
donor_state_export <- donor_state_primary$meta[, c(
  "group_id", "compartment", "donor", "state", "treatment_time", "remission",
  "inflammation", "site", "match", "n_cells", "n_samples",
  "mean_inflammation_score", "module_score")]
dataset_audit <- do.call(rbind, lapply(objects, `[[`, "audit"))

utils::write.csv(sample_state_export,
                 file.path(source_dir, "singlecell_sample_state_pseudobulk.csv"), row.names = FALSE)
utils::write.csv(sample_comp_export,
                 file.path(source_dir, "singlecell_sample_compartment_pseudobulk.csv"), row.names = FALSE)
utils::write.csv(donor_state_export,
                 file.path(source_dir, "singlecell_compartment_pseudobulk.csv"), row.names = FALSE)
utils::write.csv(donor_compartment_primary,
                 file.path(source_dir, "singlecell_donor_compartment_baseline.csv"), row.names = FALSE)
utils::write.csv(localization,
                 file.path(tables_dir, "table_singlecell_compartment_localization.csv"), row.names = FALSE)
utils::write.csv(outcome_compartment,
                 file.path(tables_dir, "table_singlecell_compartment_outcome_effects.csv"), row.names = FALSE)
utils::write.csv(outcome_sensitivity,
                 file.path(tables_dir, "table_singlecell_compartment_outcome_sensitivity.csv"), row.names = FALSE)
utils::write.csv(longitudinal_compartment,
                 file.path(tables_dir, "table_singlecell_compartment_longitudinal_effects.csv"), row.names = FALSE)
utils::write.csv(pairing_audit,
                 file.path(tables_dir, "table_singlecell_longitudinal_pairing_audit.csv"), row.names = FALSE)
utils::write.csv(outcome_state,
                 file.path(tables_dir, "table_singlecell_state_outcome_effects.csv"), row.names = FALSE)
utils::write.csv(longitudinal_state,
                 file.path(tables_dir, "table_singlecell_state_longitudinal_effects.csv"), row.names = FALSE)
utils::write.csv(dataset_audit,
                 file.path(tables_dir, "table_singlecell_dataset_level_audit.csv"), row.names = FALSE)
utils::write.csv(paired_compartment,
                 file.path(source_dir, "singlecell_compartment_paired_donors.csv"), row.names = FALSE)
utils::write.csv(paired_site_compartment,
                 file.path(source_dir, "singlecell_longitudinal_site_pairs.csv"), row.names = FALSE)
utils::write.csv(dot,
                 file.path(source_dir, "singlecell_compartment_key_gene_dotplot.csv"), row.names = FALSE)
saveRDS(list(
  algorithm_version = algorithm_version,
  amendment = "docs/analysis_amendment_2026-07-20_singlecell_sample_first.md",
  signature_md5 = signature_md5,
  gene_center = gene_center, gene_scale = gene_scale,
  mapped_genes = target_symbols, valid_genes = names(gene_scale)[valid_gene],
  legacy_gene_center = legacy_center, legacy_gene_scale = legacy_scale,
  minimum_cells = minimum_cells, chunk_nnz = chunk_nnz,
  localization_stratum = "UC; Pre; Inflamed; sample-first; donor-equal",
  longitudinal_pairing = "Match=Yes; exact donor-site; donor-equal"
), file.path(derived_dir, "cross_compartment_scaling.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(logs_dir, "14_singlecell_compartment_localization_sessionInfo.txt"))

message("SAMPLE-FIRST CROSS-COMPARTMENT ANALYSIS COMPLETE: ",
        nrow(sample_state_export), " sample-state groups; ",
        nrow(sample_comp_export), " sample-compartment groups; ",
        length(target_symbols), " common mapped genes")
print(outcome_compartment)
print(longitudinal_compartment)
print(pairing_audit)
