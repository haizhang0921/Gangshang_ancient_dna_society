#!/usr/bin/env Rscript

# ============================================================================
# Compare the two South Cemetery kinship scenarios
#   Main: SM1--SM2 = 2nd relationship, 35 +/- 26 years
#   Alt : SM1--SM2 = sibling,          26 +/- 22 years
#
# Expected result directories:
#   south_R_exact_main
#   south_R_boundaries_main
#   south_R_exact_alt_sibling
#   south_R_boundaries_alt_sibling
#
# Publication calendar ages are exported in cal BCE.
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1L) args[[1L]] else "."

main_exact <- file.path(base_dir, "south_R_exact_main")
alt_exact  <- file.path(base_dir, "south_R_exact_alt_sibling")
main_bound <- file.path(base_dir, "south_R_boundaries_main")
alt_bound  <- file.path(base_dir, "south_R_boundaries_alt_sibling")

needed <- c(
  file.path(main_exact, "burial_posterior_summary.csv"),
  file.path(alt_exact,  "burial_posterior_summary.csv"),
  file.path(main_exact, "relationship_posterior_diagnostics.csv"),
  file.path(alt_exact,  "relationship_posterior_diagnostics.csv"),
  file.path(main_bound, "boundary_duration_summary_calBCE.csv"),
  file.path(alt_bound,  "boundary_duration_summary_calBCE.csv")
)
miss <- needed[!file.exists(needed)]
if (length(miss) > 0L) {
  stop("Missing result file(s):\n", paste(miss, collapse = "\n"), call. = FALSE)
}

readu <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                                 fileEncoding = "UTF-8")

# ---------------------------------------------------------------------------
# 1. Burial-level comparison
# ---------------------------------------------------------------------------
m <- readu(file.path(main_exact, "burial_posterior_summary.csv"))
a <- readu(file.path(alt_exact,  "burial_posterior_summary.csv"))

keep_cols <- c("burial", "oxcal_mean", "posterior_mean", "posterior_median",
               "posterior_hpd95_width", "posterior_mean_change_years")
if (!all(keep_cols %in% names(m)) || !all(keep_cols %in% names(a))) {
  stop("Unexpected burial_posterior_summary.csv structure.", call. = FALSE)
}

mm <- m[, keep_cols]
aa <- a[, keep_cols]

names(mm)[-1] <- paste0("main_", names(mm)[-1])
names(aa)[-1] <- paste0("alt_sibling_", names(aa)[-1])

bur <- merge(mm, aa, by = "burial", all = TRUE, sort = FALSE)

# Internal model years are astronomical calendar years; publication values are BCE.
bur$oxcal_mean_cal_BCE <- -bur$main_oxcal_mean
bur$main_posterior_mean_cal_BCE <- -bur$main_posterior_mean
bur$alt_sibling_posterior_mean_cal_BCE <- -bur$alt_sibling_posterior_mean
bur$alt_minus_main_mean_years <- (
  bur$alt_sibling_posterior_mean_cal_BCE - bur$main_posterior_mean_cal_BCE
)
bur$absolute_mean_shift_years <- abs(bur$alt_minus_main_mean_years)

bur <- bur[order(-bur$absolute_mean_shift_years, bur$burial), ]

write.csv(
  bur,
  file.path(base_dir, "south_two_scenarios_burial_mean_comparison.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ---------------------------------------------------------------------------
# 2. Relationship diagnostics comparison
# ---------------------------------------------------------------------------
rm <- readu(file.path(main_exact, "relationship_posterior_diagnostics.csv"))
ra <- readu(file.path(alt_exact,  "relationship_posterior_diagnostics.csv"))
rm$scenario <- "Main: 2nd relationship"
ra$scenario <- "Alternative: sibling"
rel <- rbind(rm, ra)

write.csv(
  rel,
  file.path(base_dir, "south_two_scenarios_relationship_diagnostics.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ---------------------------------------------------------------------------
# 3. Boundary/duration comparison, publication scale
# ---------------------------------------------------------------------------
bm <- readu(file.path(main_bound, "boundary_duration_summary_calBCE.csv"))
ba <- readu(file.path(alt_bound,  "boundary_duration_summary_calBCE.csv"))

# Standardize parameter identifier
bm$parameter <- paste(bm$parameter_type, bm$name, sep = ": ")
ba$parameter <- paste(ba$parameter_type, ba$name, sep = ": ")

cols <- c("parameter", "parameter_type", "name", "scale", "mu", "median",
          "from_68_3", "to_68_3", "from_95_4", "to_95_4")
cols <- cols[cols %in% names(bm) & cols %in% names(ba)]

bm2 <- bm[, cols, drop = FALSE]
ba2 <- ba[, cols, drop = FALSE]

names(bm2)[names(bm2) %in% c("mu", "median", "from_68_3", "to_68_3", "from_95_4", "to_95_4")] <-
  paste0("main_", names(bm2)[names(bm2) %in% c("mu", "median", "from_68_3", "to_68_3", "from_95_4", "to_95_4")])
names(ba2)[names(ba2) %in% c("mu", "median", "from_68_3", "to_68_3", "from_95_4", "to_95_4")] <-
  paste0("alt_sibling_", names(ba2)[names(ba2) %in% c("mu", "median", "from_68_3", "to_68_3", "from_95_4", "to_95_4")])

# avoid duplicate metadata columns from merge
ba2 <- ba2[, setdiff(names(ba2), c("parameter_type", "name", "scale")), drop = FALSE]
bd <- merge(bm2, ba2, by = "parameter", all = TRUE, sort = FALSE)

if (all(c("main_mu", "alt_sibling_mu") %in% names(bd))) {
  bd$alt_minus_main_mean_years <- bd$alt_sibling_mu - bd$main_mu
  bd$absolute_mean_shift_years <- abs(bd$alt_minus_main_mean_years)
}

write.csv(
  bd,
  file.path(base_dir, "south_two_scenarios_boundary_duration_comparison.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ---------------------------------------------------------------------------
# 4. Compact key-results table
# ---------------------------------------------------------------------------
key <- data.frame(
  quantity = c(
    "SM1 posterior mean",
    "SM2 posterior mean",
    "SM1-SM2 posterior mean absolute gap",
    "Cemetery Start",
    "Cemetery End",
    "Cemetery Duration"
  ),
  unit = c("cal BCE", "cal BCE", "years", "cal BCE", "cal BCE", "years"),
  main = NA_real_,
  alt_sibling = NA_real_,
  stringsAsFactors = FALSE
)

get_bur <- function(tab, b, col) tab[tab$burial == b, col][1]
key$main[1] <- -get_bur(m, "SM1", "posterior_mean")
key$main[2] <- -get_bur(m, "SM2", "posterior_mean")
key$alt_sibling[1] <- -get_bur(a, "SM1", "posterior_mean")
key$alt_sibling[2] <- -get_bur(a, "SM2", "posterior_mean")
key$main[3] <- rm$posterior_gap_mean[1]
key$alt_sibling[3] <- ra$posterior_gap_mean[1]

for (i in 4:6) {
  nm <- c("Cemetery Start", "Cemetery End", "Cemetery Duration")[i - 3]
  key$main[i] <- bm$mu[bm$name == nm][1]
  key$alt_sibling[i] <- ba$mu[ba$name == nm][1]
}
key$alt_minus_main <- key$alt_sibling - key$main

write.csv(
  key,
  file.path(base_dir, "south_two_scenarios_key_results.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nSouth Cemetery two-scenario comparison completed.\n")
cat("Created:\n")
cat("  south_two_scenarios_burial_mean_comparison.csv\n")
cat("  south_two_scenarios_relationship_diagnostics.csv\n")
cat("  south_two_scenarios_boundary_duration_comparison.csv\n")
cat("  south_two_scenarios_key_results.csv\n")
