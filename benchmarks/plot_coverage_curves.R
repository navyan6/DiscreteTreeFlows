#!/usr/bin/env Rscript
# Plot Table 3 coverage curves:
#   K vs coverage_obs@e=2, clade_recall, mean_min_edit (+ legacy site_recall).
#
# Usage:
#   Rscript benchmarks/plot_coverage_curves.R \
#       benchmarks/results/coverage_curves_h3n2_N16_eabs.csv \
#       benchmarks/results/plots/coverage_curves_h3n2_N16_eabs.pdf

args <- commandArgs(trailingOnly = TRUE)
in_csv <- if (length(args) >= 1) args[[1]] else "benchmarks/results/coverage_curves_h3n2.csv"
out_pdf <- if (length(args) >= 2) args[[2]] else {
  sub("\\.csv$", ".pdf", in_csv)
}

if (!file.exists(in_csv)) {
  stop("CSV not found: ", in_csv)
}

d <- read.csv(in_csv, stringsAsFactors = FALSE)
if (!nrow(d)) stop("empty CSV: ", in_csv)

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

methods <- unique(d$method)
cols <- c("#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9")
col_map <- setNames(cols[seq_along(methods)], methods)

# Prefer absolute-e coverage; fall back to legacy fractional coverage
y_cov <- if ("coverage_obs_e2" %in% names(d)) "coverage_obs_e2" else "coverage"
cov_lab <- if (y_cov == "coverage_obs_e2") {
  "Coverage obs within e=2"
} else {
  sprintf("Coverage (eps=%.3g)", d$eps_frac[1])
}

n_panels <- if ("clade_recall" %in% names(d)) 4 else 3
pdf(out_pdf, width = 3.1 * n_panels, height = 4.2)
op <- par(mfrow = c(1, n_panels), mar = c(4.2, 4.2, 2.5, 1.2), oma = c(0, 0, 2, 0))

plot_metric <- function(ycol, ylab, ylim = NULL) {
  if (!ycol %in% names(d)) {
    plot.new(); title(main = paste(ycol, "(missing)")); return(invisible())
  }
  ys <- d[[ycol]]
  ok <- is.finite(ys)
  if (!any(ok)) {
    plot.new(); title(main = paste(ycol, "(all NA)")); return(invisible())
  }
  if (is.null(ylim)) {
    ylim <- range(ys[ok], na.rm = TRUE)
    if (ycol %in% c("coverage", "coverage_obs_e2", "site_recall", "clade_recall",
                    "frac_gen_e2")) {
      ylim <- c(0, max(1, ylim[2]))
    }
  }
  plot(NA, xlim = range(d$K), ylim = ylim, xlab = "Number of trees (K)",
       ylab = ylab, main = ylab)
  for (m in methods) {
    sub <- d[d$method == m, ]
    sub <- sub[order(sub$K), ]
    lines(sub$K, sub[[ycol]], type = "b", pch = 16, lwd = 2, col = col_map[[m]])
  }
  legend("best", legend = methods, col = col_map[methods], lwd = 2, pch = 16,
         bty = "n", cex = 0.8)
}

plot_metric(y_cov, cov_lab)
plot_metric("mean_min_edit", "Mean min edit (Hamming)")
plot_metric("site_recall", "Site recall")
if ("clade_recall" %in% names(d)) {
  plot_metric("clade_recall", "Clade recall")
}

N_lab <- if ("N" %in% names(d)) paste0("N=", paste(unique(d$N), collapse = ",")) else ""
n_roots <- if ("n_roots" %in% names(d)) unique(d$n_roots)[1] else NA
root_lab <- if (!is.na(n_roots)) paste0(n_roots, " roots") else "held-out roots"
title(main = paste("Table 3 coverage curves — H3N2", root_lab, N_lab),
      outer = TRUE, cex.main = 1.05)
par(op)
dev.off()

# Optional second page: coverage vs e at K_max (if e columns present)
e_cols <- grep("^coverage_obs_e[0-9]+$", names(d), value = TRUE)
if (length(e_cols) >= 2) {
  e_vals <- as.integer(sub("^coverage_obs_e", "", e_cols))
  ord <- order(e_vals)
  e_cols <- e_cols[ord]
  e_vals <- e_vals[ord]
  k_max <- max(d$K)
  dK <- d[d$K == k_max, ]
  out_e <- sub("\\.pdf$", "_vs_e.pdf", out_pdf)
  pdf(out_e, width = 5.5, height = 4.2)
  ylim <- c(0, 1)
  plot(NA, xlim = range(e_vals), ylim = ylim, xlab = "Hamming radius e (AA)",
       ylab = "Coverage obs within e",
       main = sprintf("Coverage vs e at K=%s", k_max))
  for (m in methods) {
    sub <- dK[dK$method == m, ]
    if (!nrow(sub)) next
    ys <- as.numeric(sub[1, e_cols])
    lines(e_vals, ys, type = "b", pch = 16, lwd = 2, col = col_map[[m]])
  }
  legend("best", legend = methods, col = col_map[methods], lwd = 2, pch = 16,
         bty = "n", cex = 0.8)
  dev.off()
  message("Wrote ", out_e)
}

message("Wrote ", out_pdf)
