#!/usr/bin/env Rscript

# Download and checksum-verify the public GEO bulk inputs used by the
# cross-mechanism workflow. A mismatched existing file is never silently used.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
project_root <- dirname(dirname(script_path))
raw_dir <- file.path(project_root, "data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- utils::read.csv(file.path(project_root, "config", "data_manifest.csv"),
                            stringsAsFactors = FALSE)

download_manifest_row <- function(row) {
  destination <- file.path(raw_dir, row$file)
  if (!nzchar(row$expected_md5) || !nzchar(row$source_url)) {
    stop("Manifest gate failed for ", row$dataset, "/", row$file)
  }
  if (file.exists(destination) && file.info(destination)$size > 0) {
    observed <- unname(tools::md5sum(destination))
    if (!identical(observed, row$expected_md5)) {
      stop("CHECKSUM FAILED for existing ", basename(destination),
           ": observed=", observed, "; expected=", row$expected_md5)
    }
    message("CHECKSUM PASS: ", basename(destination))
    return(invisible(destination))
  }
  partial <- paste0(destination, ".part")
  message("Downloading ", row$dataset, " from ", row$source_url)
  download.file(row$source_url, partial, mode = "wb", quiet = FALSE)
  observed <- unname(tools::md5sum(partial))
  if (!identical(observed, row$expected_md5)) {
    stop("CHECKSUM FAILED after download for ", basename(destination),
         ": observed=", observed, "; expected=", row$expected_md5)
  }
  if (!file.rename(partial, destination)) stop("Could not finalize download: ", destination)
  message("CHECKSUM PASS: ", basename(destination))
  invisible(destination)
}

bulk_manifest <- manifest[tolower(manifest$required) == "yes" & manifest$dataset != "GSE282122", ]
invisible(lapply(seq_len(nrow(bulk_manifest)), function(i) download_manifest_row(bulk_manifest[i, ])))
message("Public data are available in: ", raw_dir)
