#!/usr/bin/env Rscript

# Fail-fast dependency audit and machine-readable environment capture.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
environment_dir <- file.path(project_root, "environment")
logs_dir <- file.path(project_root, "logs")
dir.create(environment_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "GEOquery", "Biobase", "AnnotationDbi", "hthgu133pluspm.db", "hgu133plus2.db",
  "GSVA", "logistf", "metafor", "lme4", "ggplot2", "patchwork", "ggridges",
  "ggrepel", "scales", "svglite", "ragg", "rhdf5", "Matrix", "dplyr"
)
optional_bootstrap_package <- "msigdbr"
available <- vapply(c(required_packages, optional_bootstrap_package), requireNamespace,
                    logical(1), quietly = TRUE)
if (!all(available[required_packages])) {
  stop("Missing required R packages: ", paste(names(available)[!available & names(available) %in% required_packages],
                                               collapse = ", "),
       ". Run: Rscript scripts/install_dependencies.R")
}

required_tools <- c("bash", "curl")
tool_paths <- Sys.which(required_tools)
if (any(!nzchar(tool_paths))) stop("Missing command-line tools: ", paste(names(tool_paths)[!nzchar(tool_paths)], collapse = ", "))
md5_tool <- Sys.which(c("md5", "md5sum"))
if (!any(nzchar(md5_tool))) stop("Missing checksum tool: install md5 (macOS) or md5sum (Linux)")
tool_paths <- c(tool_paths, md5 = unname(md5_tool[nzchar(md5_tool)][1]))

installed <- utils::installed.packages()
package_rows <- data.frame(
  package = names(available),
  required = names(available) %in% required_packages,
  available = unname(available),
  version = vapply(names(available), function(x) {
    if (available[[x]]) as.character(utils::packageVersion(x)) else NA_character_
  }, character(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(package_rows, file.path(environment_dir, "package_versions.csv"), row.names = FALSE)

system_rows <- data.frame(
  field = c("R.version", "platform", "OS", "timezone", "locale", "project_root",
            paste0("tool.", names(tool_paths))),
  value = c(R.version.string, R.version$platform, paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "),
            Sys.timezone(), Sys.getlocale(), project_root, unname(tool_paths)),
  stringsAsFactors = FALSE
)
utils::write.csv(system_rows, file.path(environment_dir, "system_info.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(environment_dir, "sessionInfo.txt"))
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "00_preflight_sessionInfo.txt"))
message("PREFLIGHT PASS: R ", getRversion(), "; ", length(required_packages),
        " required packages; tools: ", paste(required_tools, collapse = ", "))
