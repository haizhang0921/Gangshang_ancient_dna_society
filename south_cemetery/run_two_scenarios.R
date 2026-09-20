#!/usr/bin/env Rscript

# Run the two South Cemetery kinship scenarios from start to finish.
# Main model is used in the manuscript; sibling model is a sensitivity analysis.

run_rscript <- function(script, args = character(0)) {
  rscript <- file.path(R.home("bin"), "Rscript")
  cmd <- c(script, args)
  cat("\n============================================================\n")
  cat("Running: ", paste(cmd, collapse = " "), "\n", sep = "")
  cat("============================================================\n")
  status <- system2(rscript, args = cmd)
  if (!identical(status, 0L)) {
    stop("Rscript failed: ", script, call. = FALSE)
  }
}

# All paths are relative to the current working directory.
required <- c(
  "gangshang_s.csv",
  "genetic_relations_s.csv",
  "genetic_relations_s_alt.csv",
  "01_south_exact_scenario.R",
  "02_reconstruct_south_boundaries_scenario.R",
  "03_compare_two_scenarios.R"
)
missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
  stop(
    "Missing file(s) in the working directory:\n",
    paste(missing, collapse = "\n"),
    "\n\nSet the R working directory to the extracted South Cemetery workflow folder.",
    call. = FALSE
  )
}

seed <- "20260920"
ndraw <- "12000"

# MAIN: SM1-SM2 second-degree relationship, 35 +/- 26 y
run_rscript(
  "01_south_exact_scenario.R",
  c("gangshang_s.csv", "genetic_relations_s.csv", "south_R_exact_main", seed, ndraw)
)
run_rscript(
  "02_reconstruct_south_boundaries_scenario.R",
  c("south_R_exact_main", "gangshang_s.csv", "south_R_boundaries_main", seed)
)

# ALTERNATIVE: SM1-SM2 sibling relationship, 26 +/- 22 y
run_rscript(
  "01_south_exact_scenario.R",
  c("gangshang_s.csv", "genetic_relations_s_alt.csv", "south_R_exact_alt_sibling", seed, ndraw)
)
run_rscript(
  "02_reconstruct_south_boundaries_scenario.R",
  c("south_R_exact_alt_sibling", "gangshang_s.csv", "south_R_boundaries_alt_sibling", seed)
)

# Compare scenarios
run_rscript("03_compare_two_scenarios.R", c("."))

cat("\n============================================================\n")
cat("All South Cemetery scenario analyses completed successfully.\n")
cat("============================================================\n")
