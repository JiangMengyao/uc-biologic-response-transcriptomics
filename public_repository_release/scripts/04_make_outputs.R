#!/usr/bin/env Rscript

# Create transparent, manuscript-ready draft figures from the reproduced data.

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()))
derived_dir <- file.path(project_root, "data", "derived")
fig_dir <- file.path(project_root, "results", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

meta_base <- readRDS(file.path(derived_dir, "gse92415_baseline_metadata_scores.rds"))
meta_uc <- readRDS(file.path(derived_dir, "gse16879_baseline_metadata_scores.rds"))
result <- readRDS(file.path(project_root, "results", "analysis_results.rds"))

png(file.path(fig_dir, "figure1_sample_flow.png"), width = 1500, height = 900, res = 180)
plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, 10)); par(xpd = NA)
box <- function(x1, y1, x2, y2, label, col, cex = 1.05) {
  rect(x1, y1, x2, y2, col = col, border = "grey25", lwd = 2)
  text((x1 + x2) / 2, (y1 + y2) / 2, label, cex = cex)
}
arrow <- function(x1, y1, x2, y2) arrows(x1, y1, x2, y2, length = 0.08, lwd = 2)
box(0.7, 7.4, 4.3, 9.0, "GSE92415\n183 total samples", "#d9edf7")
box(5.7, 7.4, 9.3, 9.0, "GSE16879\n133 total samples", "#d9edf7")
box(0.7, 3.8, 4.3, 6.2, "Baseline UC biopsies\nn = 87; 87 unique subjects\ngolimumab 59 / placebo 28", "#dff0d8", cex = 0.82)
box(5.7, 3.8, 9.3, 6.2, "Baseline UC before infliximab\nn = 24\nresponders 8 / non-responders 16", "#dff0d8", cex = 0.82)
box(2.2, 0.8, 7.8, 2.5, "Primary: treatment × inflammatory-program interaction\nExternal: directional association only (single arm)", "#fcf8e3", cex = 0.82)
arrow(2.5, 7.3, 2.5, 6.3); arrow(7.5, 7.3, 7.5, 6.3); arrow(2.5, 3.7, 4.4, 2.6); arrow(7.5, 3.7, 5.6, 2.6)
title("Figure 1. Reproduced study sample flow")
dev.off()

png(file.path(fig_dir, "figure2_program_response_by_treatment.png"), width = 1800, height = 800, res = 180)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
for (trt in c("golimumab", "placebo")) {
  dat <- meta_base[meta_base$treatment == trt, , drop = FALSE]
  boxplot(Inflammation ~ response, data = dat, col = c("#d95f02", "#1b9e77"),
          ylab = "GSVA inflammatory-program score", xlab = "Week 6 clinical response",
          main = trt)
  stripchart(Inflammation ~ response, data = dat, method = "jitter", vertical = TRUE,
             pch = 16, add = TRUE, col = grDevices::adjustcolor("black", 0.55))
}
dev.off()

interaction <- result$primary_model[result$primary_model$term == "treatment_binary:marker", , drop = FALSE]
png(file.path(fig_dir, "figure3_interaction_effect.png"), width = 1300, height = 650, res = 180)
par(mar = c(4, 8, 3, 2))
plot(interaction$odds_ratio, 1, log = "x", xlim = range(c(interaction$ci_lower, interaction$ci_upper, 1)),
     ylim = c(0.5, 1.5), pch = 19, yaxt = "n", ylab = "", xlab = "Odds ratio (log scale)",
     main = "Treatment × inflammatory-program interaction")
segments(interaction$ci_lower, 1, interaction$ci_upper, 1, lwd = 3)
abline(v = 1, lty = 2, col = "grey40")
axis(2, at = 1, labels = "Treatment × program")
mtext(sprintf("OR %.2f (95%% CI %.3f–%.2f), P = %.3f", interaction$odds_ratio,
              interaction$ci_lower, interaction$ci_upper, interaction$p_value), side = 3, line = 0.2)
dev.off()

heat <- as.matrix(result$robustness[, c("golimumab_rho", "placebo_rho", "infliximab_rho")])
rownames(heat) <- result$robustness$set
png(file.path(fig_dir, "figure4_robustness_directionality.png"), width = 1400, height = 900, res = 180)
par(mar = c(6, 12, 3, 2))
image(x = seq_len(ncol(heat)), y = seq_len(nrow(heat)), z = t(heat[nrow(heat):1, , drop = FALSE]),
      col = colorRampPalette(c("#2166ac", "white", "#b2182b"))(101), zlim = c(-1, 1),
      xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "Directionality across pre-specified gene sets")
axis(1, at = seq_len(ncol(heat)), labels = colnames(heat), las = 2)
axis(2, at = seq_len(nrow(heat)), labels = rev(rownames(heat)), las = 2)
for (i in seq_len(nrow(heat))) for (j in seq_len(ncol(heat))) text(j, nrow(heat) - i + 1, sprintf("%.2f", heat[i, j]))
dev.off()

writeLines(capture.output(sessionInfo()), file.path(project_root, "logs", "04_make_outputs_sessionInfo.txt"))
message("Draft figures written to: ", fig_dir)
