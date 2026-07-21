#!/usr/bin/env Rscript

# Reconstructed simulation-based power analysis for the treatment-by-program
# interaction. It uses the observed 59:28 treatment allocation, the observed
# average treatment effect (OR 1.83), and a pure predictive interaction.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
tables_dir <- file.path(project_root, "results", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

simulate_interaction_power <- function(n, interaction_or, marker = c("continuous", "binary"),
                                       simulations = 4000L, seed = 2026L) {
  marker <- match.arg(marker)
  set.seed(seed)
  n_treated <- round(n * 59 / 87)
  treatment <- c(rep(1L, n_treated), rep(0L, n - n_treated))
  p_values <- replicate(simulations, {
    program <- if (marker == "continuous") stats::rnorm(n) else stats::rbinom(n, 1L, 0.5)
    probability <- stats::plogis(-0.44 + log(1.83) * treatment + log(interaction_or) * treatment * program)
    response <- stats::rbinom(n, 1L, probability)
    fit <- suppressWarnings(stats::glm(response ~ treatment * program, family = stats::binomial()))
    coefficients <- summary(fit)$coefficients
    if (!"treatment:program" %in% rownames(coefficients)) return(NA_real_)
    coefficients["treatment:program", "Pr(>|z|)"]
  })
  data.frame(
    n = n, marker = marker, interaction_or = interaction_or, simulations = simulations,
    usable_simulations = sum(!is.na(p_values)), power = mean(p_values < 0.05, na.rm = TRUE)
  )
}

power_grid <- do.call(rbind, c(
  lapply(c(2.0, 2.5, 3.0), function(or) simulate_interaction_power(87, or, "continuous")),
  lapply(c(2.0, 2.5, 3.0), function(or) simulate_interaction_power(87, or, "binary"))
))
sample_size_scan <- do.call(rbind, lapply(c(87L, 130L, 175L, 220L, 270L), function(n) {
  do.call(rbind, lapply(c(2.5, 3.0), function(or) simulate_interaction_power(n, or, "continuous")))
}))

utils::write.csv(power_grid, file.path(tables_dir, "tableS1_interaction_power.csv"), row.names = FALSE)
utils::write.csv(sample_size_scan, file.path(tables_dir, "tableS1_sample_size_scan.csv"), row.names = FALSE)
print(power_grid)
print(sample_size_scan)
writeLines(capture.output(sessionInfo()), file.path(project_root, "logs", "05_interaction_power_sessionInfo.txt"))
