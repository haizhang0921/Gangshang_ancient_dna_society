#!/usr/bin/env Rscript

# Run the Northern Cemetery chronology under three alternative kinship specifications.
# Formal settings are identical across scenarios:
# 4 chains x 80,000 iterations; burn-in 20,000; thin 20; seed 20260806.

seed <- 20260806L
oxcal <- "gangshang_n.csv"
primary_script <- "01_north_primary_scenario.R"
ordered_script <- "02_north_stage_order_scenario.R"
boundary_script <- "03_reconstruct_boundaries_scenario.R"

scenarios <- data.frame(
  scenario = c("main", "alt_close", "alt_far"),
  relationship_file = c(
    "genetic_relations_n.csv",
    "genetic_relations_n_alt.csv",
    "genetic_relations_n_alt2.csv"
  ),
  stringsAsFactors = FALSE
)

required <- c(oxcal, primary_script, ordered_script, boundary_script, scenarios$relationship_file)
missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
  stop("Missing required file(s):\n", paste(missing, collapse = "\n"), call. = FALSE)
}

rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) {
  cand <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  if (file.exists(cand)) rscript <- cand
}
if (!nzchar(rscript)) stop("Rscript executable could not be found.", call. = FALSE)

run_one <- function(script, args) {
  cat("\n============================================================\n")
  cat("Running: ", script, " ", paste(args, collapse = " "), "\n", sep = "")
  cat("============================================================\n")
  status <- system2(rscript, c(script, args), stdout = "", stderr = "")
  if (!identical(status, 0L)) stop("Rscript failed: ", script, call. = FALSE)
}

is_primary_complete <- function(out_dir) {
  all(file.exists(file.path(out_dir, c(
    "burial_posterior_density.csv",
    "joint_posterior_samples.csv",
    "run_metadata.csv",
    "relationships_used.csv"
  ))))
}

is_boundary_complete <- function(out_dir) {
  all(file.exists(file.path(out_dir, c(
    "kinship_updated_oxcal_like_summary.csv",
    "stage_duration_summary.csv",
    "boundary_summary.csv"
  ))))
}

for (i in seq_len(nrow(scenarios))) {
  s <- scenarios$scenario[i]
  rel <- scenarios$relationship_file[i]
  primary_out <- paste0("north_R_primary_", s)
  ordered_out <- paste0("north_R_stage_order_", s)
  boundary_out <- paste0("north_R_boundaries_", s)

  if (is_primary_complete(primary_out)) {
    cat("\nSkipping completed primary model: ", primary_out, "\n", sep = "")
  } else {
    run_one(primary_script, c(oxcal, rel, primary_out, seed))
  }

  if (is_primary_complete(ordered_out)) {
    cat("\nSkipping completed Stage-order model: ", ordered_out, "\n", sep = "")
  } else {
    run_one(ordered_script, c(oxcal, rel, ordered_out, seed))
  }

  if (is_boundary_complete(boundary_out)) {
    cat("\nSkipping completed boundary reconstruction: ", boundary_out, "\n", sep = "")
  } else {
    run_one(boundary_script, c(ordered_out, boundary_out, seed))
  }
}

cat("\nAll three kinship scenarios completed.\n")
cat("Now run: source(\"04_compare_three_scenarios.R\")\n")
