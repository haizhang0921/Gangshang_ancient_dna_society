#!/usr/bin/env Rscript

# ============================================================================
# Gangshang Southern Cemetery kinship sensitivity figure
# Formal R outputs are used as input; this script performs plotting only.
#
# Panels:
#   a, SM1 OxCal vs main 2nd-degree vs sibling posterior
#   b, SM2 OxCal vs main 2nd-degree vs sibling posterior
#   c, posterior distribution of |SM1-SM2| temporal separation
#   d, absolute changes in key posterior means under the sibling alternative
#
# Usage:
#   source("05_plot_south_kinship_sensitivity.R")
#
# Expected directories in working directory:
#   south_R_exact_main/
#   south_R_exact_alt_sibling/
#   south_two_scenarios_key_results.csv
# ============================================================================

main_dir <- "south_R_exact_main"
alt_dir  <- "south_R_exact_alt_sibling"
key_file <- "south_two_scenarios_key_results.csv"

required <- c(
  file.path(main_dir, "burial_posterior_density.csv"),
  file.path(alt_dir, "burial_posterior_density.csv"),
  file.path(main_dir, "relationship_gap_posterior_density.csv"),
  file.path(alt_dir, "relationship_gap_posterior_density.csv"),
  key_file
)

missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
  stop("Missing required file(s):\n", paste(missing, collapse = "\n"), call. = FALSE)
}

main <- read.csv(file.path(main_dir, "burial_posterior_density.csv"), check.names = FALSE)
alt  <- read.csv(file.path(alt_dir,  "burial_posterior_density.csv"), check.names = FALSE)
gap_main <- read.csv(file.path(main_dir, "relationship_gap_posterior_density.csv"), check.names = FALSE)
gap_alt  <- read.csv(file.path(alt_dir,  "relationship_gap_posterior_density.csv"), check.names = FALSE)
key <- read.csv(key_file, check.names = FALSE)

# restrained publication palette
col_oxcal <- "#4C78A8"
col_main  <- "#F58518"
col_alt   <- "#54A24B"
col_grid  <- "#E6E6E6"

plot_burial <- function(burial, panel_letter) {
  dm <- main[main$burial == burial, ]
  da <- alt[alt$burial == burial, ]
  if (nrow(dm) == 0L || nrow(da) == 0L) stop("Missing burial: ", burial)

  x <- -dm$calendar_year
  ymax <- max(c(dm$oxcal_density, dm$kinship_posterior_density, da$kinship_posterior_density), na.rm = TRUE) * 1.08

  plot(
    range(x), c(0, ymax), type = "n", xlim = rev(range(x)),
    xlab = "Calendar age (cal BCE)", ylab = "Posterior density",
    main = burial, las = 1
  )
  abline(v = axTicks(1), col = col_grid, lwd = 0.7)
  lines(x, dm$oxcal_density, col = col_oxcal, lwd = 1.5)
  lines(x, dm$kinship_posterior_density, col = col_main, lwd = 1.7)
  lines(-da$calendar_year, da$kinship_posterior_density, col = col_alt, lwd = 1.7)
  mtext(panel_letter, side = 3, adj = -0.12, line = 0.3, font = 2, cex = 1.2)
}

make_figure <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))

  par(mfrow = c(2, 2), mar = c(4.3, 4.5, 2.3, 1.1), oma = c(0.2, 0.2, 0.2, 0.2), family = "sans")

  plot_burial("SM1", "a")
  legend(
    "topleft",
    legend = c("OxCal archaeological posterior", "Main: 2nd-degree", "Alternative: sibling"),
    col = c(col_oxcal, col_main, col_alt), lwd = c(1.5, 1.7, 1.7),
    bty = "n", cex = 0.78
  )

  plot_burial("SM2", "b")

  ymax <- max(c(gap_main$probability_density, gap_alt$probability_density), na.rm = TRUE) * 1.08
  plot(
    gap_main$gap_years, gap_main$probability_density,
    type = "l", lwd = 1.7, col = col_main,
    xlab = "|SM1 - SM2| temporal separation (years)",
    ylab = "Posterior density", main = "SM1-SM2 absolute gap",
    ylim = c(0, ymax), las = 1
  )
  lines(gap_alt$gap_years, gap_alt$probability_density, lwd = 1.7, col = col_alt)
  abline(v = 35, lty = 2, lwd = 0.8, col = col_main)
  abline(v = 26, lty = 2, lwd = 0.8, col = col_alt)
  legend("topright", legend = c("Main: 2nd-degree", "Alternative: sibling"), col = c(col_main, col_alt), lwd = 1.7, bty = "n", cex = 0.78)
  mtext("c", side = 3, adj = -0.12, line = 0.3, font = 2, cex = 1.2)

  qnames <- c(
    "SM2 posterior mean",
    "SM1 posterior mean",
    "Cemetery Start",
    "Cemetery End",
    "Cemetery Duration"
  )
  vals <- abs(key$alt_minus_main[match(qnames, key$quantity)])
  ypos <- rev(seq_along(qnames))
  plot(
    c(0, max(vals, na.rm = TRUE) * 1.25), range(ypos), type = "n",
    xlab = "Absolute shift from main model (years)", ylab = "",
    yaxt = "n", main = "Sensitivity of posterior means", las = 1
  )
  axis(2, at = ypos, labels = qnames, las = 1, tick = FALSE, cex.axis = 0.78)
  abline(v = axTicks(1), col = col_grid, lwd = 0.7)
  segments(0, ypos, vals, ypos, lwd = 1.0, col = col_main)
  points(vals, ypos, pch = 16, cex = 0.9, col = col_main)
  text(vals + max(vals) * 0.025, ypos, labels = sprintf("%.2f", vals), pos = 4, cex = 0.75)
  mtext("d", side = 3, adj = -0.12, line = 0.3, font = 2, cex = 1.2)
}

pdf("Gangshang_South_kinship_sensitivity_figure_R.pdf", width = 10.2, height = 7.2, family = "sans")
make_figure()
dev.off()

png("Gangshang_South_kinship_sensitivity_figure_R.png", width = 3060, height = 2160, res = 300)
make_figure()
dev.off()

tiff("Gangshang_South_kinship_sensitivity_figure_R_600dpi.tiff", width = 10.2, height = 7.2, units = "in", res = 600, compression = "lzw")
make_figure()
dev.off()

cat("South sensitivity figure completed.\n")
