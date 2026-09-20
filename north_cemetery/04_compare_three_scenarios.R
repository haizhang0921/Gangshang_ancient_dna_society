#!/usr/bin/env Rscript

# Summarize the three Northern Cemetery kinship-sensitivity scenarios.

scenarios <- c(main = "main", alt_close = "alt_close", alt_far = "alt_far")

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing result file: ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

# Burial posterior comparison -------------------------------------------------
burial_tables <- lapply(names(scenarios), function(s) {
  d <- read_required(file.path(paste0("north_R_primary_", s), "burial_posterior_summary.csv"))
  d$scenario <- s
  d
})
names(burial_tables) <- names(scenarios)

getcol <- function(d, candidates) {
  hit <- candidates[candidates %in% names(d)]
  if (length(hit) == 0L) return(rep(NA_real_, nrow(d)))
  d[[hit[1L]]]
}

burial_long <- do.call(rbind, lapply(names(burial_tables), function(s) {
  d <- burial_tables[[s]]
  data.frame(
    scenario = s,
    burial = d$burial,
    stage = d$stage,
    posterior_mean = getcol(d, c("posterior_mean", "kinship_mean")),
    posterior_median = getcol(d, c("posterior_median", "kinship_median", "median")),
    rhat = getcol(d, c("rhat")),
    approximate_ess = getcol(d, c("approximate_ess", "ess")),
    stringsAsFactors = FALSE
  )
}))
write.csv(burial_long, "three_scenarios_burial_summary_long.csv", row.names = FALSE)

# Wide table of burial means
wide <- Reduce(function(x, y) merge(x, y, by = c("burial", "stage"), all = TRUE),
  lapply(names(burial_tables), function(s) {
    d <- burial_tables[[s]]
    z <- data.frame(
      burial = d$burial,
      stage = d$stage,
      value = getcol(d, c("posterior_mean", "kinship_mean")),
      stringsAsFactors = FALSE
    )
    names(z)[3] <- paste0("mean_", s)
    z
  })
)
wide$alt_close_minus_main <- wide$mean_alt_close - wide$mean_main
wide$alt_far_minus_main <- wide$mean_alt_far - wide$mean_main
write.csv(wide, "three_scenarios_burial_mean_comparison.csv", row.names = FALSE)

# Relationship diagnostics ---------------------------------------------------
rel_long <- do.call(rbind, lapply(names(scenarios), function(s) {
  d <- read_required(file.path(paste0("north_R_primary_", s), "relationship_posterior_diagnostics.csv"))
  d$scenario <- s
  d
}))
write.csv(rel_long, "three_scenarios_relationship_diagnostics.csv", row.names = FALSE)

# MCMC diagnostics summary ---------------------------------------------------
mcmc_summary <- do.call(rbind, lapply(names(scenarios), function(s) {
  d <- read_required(file.path(paste0("north_R_primary_", s), "mcmc_burial_diagnostics.csv"))
  rhat <- getcol(d, c("rhat", "Rhat"))
  ess <- getcol(d, c("approximate_ess", "ESS", "ess"))
  acc <- getcol(d, c("single_site_acceptance_rate", "acceptance_rate", "acceptance", "single_site_acceptance"))
  data.frame(
    scenario = s,
    maximum_rhat = max(rhat, na.rm = TRUE),
    minimum_approximate_ess = min(ess, na.rm = TRUE),
    median_approximate_ess = median(ess, na.rm = TRUE),
    minimum_acceptance = if (all(is.na(acc))) NA_real_ else min(acc, na.rm = TRUE),
    maximum_acceptance = if (all(is.na(acc))) NA_real_ else max(acc, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.csv(mcmc_summary, "three_scenarios_primary_mcmc_summary.csv", row.names = FALSE)

# Boundary and duration comparisons -----------------------------------------
read_boundary_summary <- function(s) {
  d <- read_required(file.path(paste0("north_R_boundaries_", s), "boundary_summary.csv"))
  d$scenario <- s
  d
}
read_duration_summary <- function(s) {
  d <- read_required(file.path(paste0("north_R_boundaries_", s), "stage_duration_summary.csv"))
  d$scenario <- s
  d
}

boundary_long <- do.call(rbind, lapply(names(scenarios), read_boundary_summary))
duration_long <- do.call(rbind, lapply(names(scenarios), read_duration_summary))
write.csv(boundary_long, "three_scenarios_boundary_summary_long.csv", row.names = FALSE)
write.csv(duration_long, "three_scenarios_duration_summary_long.csv", row.names = FALSE)

cat("Comparison files written successfully.\n")
