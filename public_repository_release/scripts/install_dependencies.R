#!/usr/bin/env Rscript

# Install the packages required by the reproducible workflow. This script is
# intentionally separate from run_all.R: analysis runs never modify the R
# library as a side effect.

cran_packages <- c(
  "logistf", "metafor", "lme4", "ggplot2", "patchwork", "ggridges",
  "ggrepel", "scales", "svglite", "ragg", "dplyr", "msigdbr"
)
bioconductor_packages <- c(
  "GEOquery", "Biobase", "AnnotationDbi", "hthgu133pluspm.db",
  "hgu133plus2.db", "GSVA", "rhdf5"
)

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) utils::install.packages(missing_cran, repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
missing_bioc <- bioconductor_packages[!vapply(bioconductor_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
message("Dependency installation complete. Run Rscript scripts/00_preflight.R to verify.")

