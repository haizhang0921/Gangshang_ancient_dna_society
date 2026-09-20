#!/usr/bin/env Rscript

# ======================================================================
# Northern Cemetery:
# Reconstruct boundaries and Stage durations after scenario-specific kinship update
#
# PURPOSE
# -------
# This script starts from the joint posterior burial-event samples produced
# by a scenario-specific kinship sensitivity model with:
#
#   ENFORCE_GLOBAL_STAGE_ORDER=TRUE
#
# It then reconstructs four ordered phase boundaries under the same basic
# "uniform events within each Phase" assumption used by a simple OxCal Phase.
# From those boundary draws it calculates:
#
#   Stage I Duration   = Boundary(Stage I to II) - Boundary(Cemetery Start)
#   Stage II Duration  = Boundary(Stage II to III) - Boundary(Stage I to II)
#   Stage III Duration = Boundary(Cemetery End) - Boundary(Stage II to III)
#
# It outputs:
#   1. An OxCal-like model summary table in cal BP / years
#   2. Full fixed-grid posterior probability densities
#   3. HPD interval tables
#   4. Boundary and duration posterior samples
#
# IMPORTANT
# ---------
# This is a modular reconstruction:
#
#   updated burial posterior × conditional Phase-boundary model
#
# It is not an exact recovery of OxCal's original joint MCMC because the
# exported OxCal CSV contains marginal distributions rather than joint draws.
#
# COMMAND LINE
# ------------
# Rscript reconstruct_stage_boundaries_durations_corrected.R \
#   north_R_stage_order \
#   north_R_boundaries \
#   20260806
#
# In RStudio, with the sensitivity-model result folder in the working
# directory, simply use:
#
#   source("reconstruct_stage_boundaries_durations_corrected.R")
#
#
# ENVIRONMENT OPTIONS
# -------------------
# BOUNDARY_GIBBS_BURN=20
# BOUNDARY_GIBBS_THIN=4
# BOUNDARY_DRAWS_PER_SAMPLE=1
# BOUNDARY_GRID_STEP=0.5
# MAX_JOINT_SAMPLES=0        # 0 means use all samples
# ======================================================================

args <- commandArgs(trailingOnly = TRUE)

# The first argument is the directory produced by a scenario-specific
# global-Stage-order sensitivity model.
results_dir <- if (length(args) >= 1L) {
  args[1L]
} else {
  "north_R_stage_order"
}

output_dir <- if (length(args) >= 2L) {
  args[2L]
} else {
  "north_R_boundaries"
}

seed <- if (length(args) >= 3L) {
  as.integer(args[3L])
} else {
  20260806L
}

joint_samples_file <- file.path(
  results_dir,
  "joint_posterior_samples.csv"
)

burial_density_file <- file.path(
  results_dir,
  "burial_posterior_density.csv"
)

burials_used_file <- file.path(
  results_dir,
  "burials_used.csv"
)

source_metadata_file <- file.path(
  results_dir,
  "run_metadata.csv"
)

relationships_used_file <- file.path(
  results_dir,
  "relationships_used.csv"
)

get_env_integer <- function(name, default, minimum = 0L) {
  value <- suppressWarnings(
    as.integer(Sys.getenv(name, unset = as.character(default)))
  )

  if (is.na(value) || value < minimum) {
    value <- default
  }

  value
}

get_env_numeric <- function(name, default, minimum = 0) {
  value <- suppressWarnings(
    as.numeric(Sys.getenv(name, unset = as.character(default)))
  )

  if (!is.finite(value) || value <= minimum) {
    value <- default
  }

  value
}

gibbs_burn <- get_env_integer(
  "BOUNDARY_GIBBS_BURN",
  default = 20L,
  minimum = 1L
)

gibbs_thin <- get_env_integer(
  "BOUNDARY_GIBBS_THIN",
  default = 4L,
  minimum = 1L
)

draws_per_sample <- get_env_integer(
  "BOUNDARY_DRAWS_PER_SAMPLE",
  default = 1L,
  minimum = 1L
)

boundary_grid_step <- get_env_numeric(
  "BOUNDARY_GRID_STEP",
  default = 0.5,
  minimum = 0
)

max_joint_samples <- get_env_integer(
  "MAX_JOINT_SAMPLES",
  default = 0L,
  minimum = 0L
)

required_input_files <- c(
  joint_samples_file,
  burial_density_file,
  burials_used_file,
  source_metadata_file,
  relationships_used_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0L) {
  stop(
    "Missing required sensitivity-model file(s):\n",
    paste(missing_input_files, collapse = "\n"),
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

# Validate that the input really is the corrected global-Stage-order model.
source_metadata <- read.csv(
  source_metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

if (!all(c("setting", "value") %in% names(source_metadata))) {
  stop("run_metadata.csv has an unexpected structure.", call. = FALSE)
}

metadata_value <- setNames(
  as.character(source_metadata$value),
  as.character(source_metadata$setting)
)

if (!identical(metadata_value[["enforce_global_stage_order"]], "TRUE")) {
  stop(
    "Boundary reconstruction requires a sensitivity model ",
    "with enforce_global_stage_order = TRUE.",
    call. = FALSE
  )
}

# The boundary reconstruction is intentionally scenario-agnostic.
# The three kinship scenarios differ in relationship degree and/or whether
# NM30-NM27 has an independently established direction.  Those choices have
# already been incorporated into the joint posterior samples produced by the
# Stage-order model.  Boundary reconstruction therefore must NOT hard-code a
# particular number or identity of directed relationships.
#
# We retain only structural validation that is relevant to this step:
#   (i) the source model enforced the global Stage order; and
#   (ii) the relationship audit table is readable and consistent with the
#        metadata count.  The actual relationship specification is inherited
#        from the source model rather than redefined here.

relationships_used <- read.csv(
  relationships_used_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

required_relationship_columns <- c(
  "burial_a", "burial_b", "relationship", "mean", "sd",
  "direction_established", "effective_direction"
)
missing_relationship_columns <- setdiff(
  required_relationship_columns,
  names(relationships_used)
)
if (length(missing_relationship_columns) > 0L) {
  stop(
    "relationships_used.csv is missing required column(s): ",
    paste(missing_relationship_columns, collapse = ", "),
    call. = FALSE
  )
}

metadata_relationship_count <- suppressWarnings(
  as.integer(metadata_value[["relationship_count"]])
)
if (is.finite(metadata_relationship_count) &&
    metadata_relationship_count != nrow(relationships_used)) {
  stop(
    "Relationship count in run_metadata.csv does not match relationships_used.csv.",
    call. = FALSE
  )
}

metadata_direction_count <- suppressWarnings(
  as.integer(metadata_value[["established_direction_count"]])
)

# relationships_used.csv preserves the original text values from the input
# relationship table (typically YES/NO).  Parse these with the same semantics
# as the modelling scripts instead of comparing only with TRUE/1.
parse_direction_flag <- function(x) {
  x_chr <- toupper(trimws(as.character(x)))
  x_chr %in% c("YES", "Y", "TRUE", "T", "1")
}

observed_direction_flags <- parse_direction_flag(
  relationships_used$direction_established
)
observed_direction_count <- sum(observed_direction_flags, na.rm = TRUE)

if (is.finite(metadata_direction_count) &&
    metadata_direction_count != observed_direction_count) {
  stop(
    paste0(
      "Established-direction count mismatch after parsing YES/NO flags: ",
      "run_metadata.csv = ", metadata_direction_count,
      "; relationships_used.csv = ", observed_direction_count, "."
    ),
    call. = FALSE
  )
}

observed_directions_text <- metadata_value[["established_directions"]]
if (is.null(observed_directions_text) || is.na(observed_directions_text) ||
    !nzchar(trimws(observed_directions_text))) {
  observed_directions_text <- "none"
}

cat("Validated scenario-specific global-Stage-order source model:\n")
cat("  global Stage order = TRUE\n")
cat("  relationships = ", nrow(relationships_used), "\n", sep = "")
cat("  independently established directions = ", observed_direction_count, "\n", sep = "")
cat("  directions = ", observed_directions_text, "\n\n", sep = "")

# ----------------------------------------------------------------------
# 1. Stage definitions from the corrected sensitivity-model output
# ----------------------------------------------------------------------

burials_used <- read.csv(
  burials_used_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

required_burial_columns <- c("burial", "stage")
missing_burial_columns <- setdiff(
  required_burial_columns,
  names(burials_used)
)

if (length(missing_burial_columns) > 0L) {
  stop(
    "burials_used.csv is missing: ",
    paste(missing_burial_columns, collapse = ", "),
    call. = FALSE
  )
}

burials_used$burial <- trimws(as.character(burials_used$burial))
burials_used$stage <- trimws(as.character(burials_used$stage))

stage_I <- burials_used$burial[
  burials_used$stage == "Stage I"
]
stage_II <- burials_used$burial[
  burials_used$stage == "Stage II"
]
stage_III <- burials_used$burial[
  burials_used$stage == "Stage III"
]

if (
  length(stage_I) != 22L ||
  length(stage_II) != 5L ||
  length(stage_III) != 7L
) {
  stop(
    "Unexpected Stage counts. Detected Stage I/II/III = ",
    length(stage_I), "/", length(stage_II), "/", length(stage_III),
    "; expected 22/5/7.",
    call. = FALSE
  )
}

stage_map <- setNames(
  burials_used$stage,
  burials_used$burial
)

stage_event_count <- c(
  "Stage I" = length(stage_I),
  "Stage II" = length(stage_II),
  "Stage III" = length(stage_III)
)

parameter_order <- c(
  "Cemetery Start",
  stage_I,
  "Stage I Duration",
  "Stage I to II",
  stage_II,
  "Stage II Duration",
  "Stage II to III",
  stage_III,
  "Stage III Duration",
  "Cemetery End"
)

cat(
  "Stage counts: I=", length(stage_I),
  ", II=", length(stage_II),
  ", III=", length(stage_III), "\n\n",
  sep = ""
)

# ----------------------------------------------------------------------
# 2. Read and validate joint posterior burial samples
# ----------------------------------------------------------------------

joint <- read.csv(
  joint_samples_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

required_burials <- c(stage_I, stage_II, stage_III)
missing_burials <- setdiff(required_burials, names(joint))

if (length(missing_burials) > 0L) {
  stop(
    "联合样本缺少墓葬：",
    paste(missing_burials, collapse = ", "),
    call. = FALSE
  )
}

for (burial in required_burials) {
  joint[[burial]] <- as.numeric(joint[[burial]])
}

bad_joint_rows <- !apply(
  joint[, required_burials, drop = FALSE],
  1L,
  function(x) all(is.finite(x))
)

if (any(bad_joint_rows)) {
  warning(sum(bad_joint_rows), " 行联合样本包含无效值，已删除。")
  joint <- joint[!bad_joint_rows, , drop = FALSE]
}

if (nrow(joint) == 0L) {
  stop("没有可用的联合后验样本。", call. = FALSE)
}

# The boundary reconstruction requires strict global Stage order.
stage_I_latest <- apply(joint[, stage_I, drop = FALSE], 1L, max)
stage_II_earliest <- apply(joint[, stage_II, drop = FALSE], 1L, min)
stage_II_latest <- apply(joint[, stage_II, drop = FALSE], 1L, max)
stage_III_earliest <- apply(joint[, stage_III, drop = FALSE], 1L, min)

global_order_valid <- (
  stage_I_latest < stage_II_earliest &
  stage_II_latest < stage_III_earliest
)

if (!all(global_order_valid)) {
  stop(
    "联合样本中有 ",
    sum(!global_order_valid),
    " 行不满足 Stage I < Stage II < Stage III。",
    "\n请使用 ENFORCE_GLOBAL_STAGE_ORDER=TRUE 的敏感性模型结果。",
    call. = FALSE
  )
}

if (max_joint_samples > 0L && nrow(joint) > max_joint_samples) {
  keep <- sort(
    sample.int(
      nrow(joint),
      size = max_joint_samples,
      replace = FALSE
    )
  )
  joint <- joint[keep, , drop = FALSE]
}

n_joint <- nrow(joint)

cat("Joint posterior samples used: ", n_joint, "\n", sep = "")
cat("Boundary draws per burial sample: ", draws_per_sample, "\n", sep = "")
cat("Boundary Gibbs burn: ", gibbs_burn, "\n", sep = "")
cat("Boundary Gibbs thin: ", gibbs_thin, "\n\n", sep = "")

# ----------------------------------------------------------------------
# 3. Boundary reconstruction under uniform Phase likelihood
# ----------------------------------------------------------------------

# For one Stage with n dated events uniformly distributed between its
# boundaries, the conditional phase likelihood is proportional to:
#
#   duration^(-n)
#
# with the boundaries outside the earliest/latest event.
#
# Shared internal boundaries couple adjacent Stage durations, so we use
# a short Gibbs sampler for each retained burial-event state.

sample_pareto_duration <- function(minimum_duration, shape) {
  if (!is.finite(minimum_duration) || minimum_duration <= 0) {
    stop("Invalid minimum duration in Pareto draw.", call. = FALSE)
  }

  if (!is.finite(shape) || shape <= 0) {
    stop("Invalid Pareto shape.", call. = FALSE)
  }

  u <- runif(1L)
  minimum_duration * u^(-1 / shape)
}

sample_log_density_grid <- function(
  lower,
  upper,
  log_density_function,
  step = 0.5
) {
  if (!is.finite(lower) || !is.finite(upper) || upper <= lower) {
    stop(
      "Invalid conditional boundary interval: ",
      lower,
      " to ",
      upper,
      call. = FALSE
    )
  }

  width <- upper - lower

  if (width <= step) {
    grid <- seq(
      lower + width / 20,
      upper - width / 20,
      length.out = 25L
    )
  } else {
    grid <- seq(
      lower + min(step / 2, width / 20),
      upper - min(step / 2, width / 20),
      by = step
    )

    if (length(grid) < 25L) {
      grid <- seq(
        lower + width / 100,
        upper - width / 100,
        length.out = 50L
      )
    }
  }

  log_weights <- log_density_function(grid)

  if (!all(is.finite(log_weights))) {
    finite <- is.finite(log_weights)

    if (!any(finite)) {
      stop("Conditional boundary density has no finite support.", call. = FALSE)
    }

    grid <- grid[finite]
    log_weights <- log_weights[finite]
  }

  log_weights <- log_weights - max(log_weights)
  weights <- exp(log_weights)
  weights <- weights / sum(weights)

  selected <- sample.int(
    length(grid),
    size = 1L,
    prob = weights
  )

  if (length(grid) == 1L) {
    return(grid[selected])
  }

  local_half_width <- if (selected == 1L) {
    (grid[2L] - grid[1L]) / 2
  } else if (selected == length(grid)) {
    (grid[length(grid)] - grid[length(grid) - 1L]) / 2
  } else {
    min(
      grid[selected] - grid[selected - 1L],
      grid[selected + 1L] - grid[selected]
    ) / 2
  }

  runif(
    1L,
    min = max(lower, grid[selected] - local_half_width),
    max = min(upper, grid[selected] + local_half_width)
  )
}

gibbs_boundary_step <- function(
  boundary_state,
  event_limits,
  event_counts,
  grid_step
) {
  b0 <- boundary_state[1L]
  b1 <- boundary_state[2L]
  b2 <- boundary_state[3L]
  b3 <- boundary_state[4L]

  earliest_I <- event_limits[1L]
  latest_I <- event_limits[2L]
  earliest_II <- event_limits[3L]
  latest_II <- event_limits[4L]
  earliest_III <- event_limits[5L]
  latest_III <- event_limits[6L]

  n_I <- event_counts[1L]
  n_II <- event_counts[2L]
  n_III <- event_counts[3L]

  # External start boundary:
  # p(duration_I | b1, events) proportional to duration_I^(-n_I).
  minimum_duration_I <- b1 - earliest_I
  duration_I <- sample_pareto_duration(
    minimum_duration_I,
    shape = n_I - 1
  )
  b0 <- b1 - duration_I

  # External end boundary.
  minimum_duration_III <- latest_III - b2
  duration_III <- sample_pareto_duration(
    minimum_duration_III,
    shape = n_III - 1
  )
  b3 <- b2 + duration_III

  # Shared boundary between Stage I and Stage II.
  b1 <- sample_log_density_grid(
    lower = latest_I,
    upper = earliest_II,
    log_density_function = function(x) {
      -n_I * log(x - b0) -
        n_II * log(b2 - x)
    },
    step = grid_step
  )

  # Shared boundary between Stage II and Stage III.
  b2 <- sample_log_density_grid(
    lower = latest_II,
    upper = earliest_III,
    log_density_function = function(x) {
      -n_II * log(x - b1) -
        n_III * log(b3 - x)
    },
    step = grid_step
  )

  c(b0, b1, b2, b3)
}

draw_boundaries_for_event_state <- function(
  event_limits,
  event_counts,
  burn,
  thin,
  number_of_draws,
  grid_step
) {
  earliest_I <- event_limits[1L]
  latest_I <- event_limits[2L]
  earliest_II <- event_limits[3L]
  latest_II <- event_limits[4L]
  earliest_III <- event_limits[5L]
  latest_III <- event_limits[6L]

  b1 <- (latest_I + earliest_II) / 2
  b2 <- (latest_II + earliest_III) / 2

  boundary_state <- c(
    earliest_I - 10,
    b1,
    b2,
    latest_III + 10
  )

  for (iteration in seq_len(burn)) {
    boundary_state <- gibbs_boundary_step(
      boundary_state,
      event_limits,
      event_counts,
      grid_step
    )
  }

  output <- matrix(
    NA_real_,
    nrow = number_of_draws,
    ncol = 4L,
    dimnames = list(
      NULL,
      c(
        "Cemetery Start",
        "Stage I to II",
        "Stage II to III",
        "Cemetery End"
      )
    )
  )

  for (draw in seq_len(number_of_draws)) {
    for (iteration in seq_len(thin)) {
      boundary_state <- gibbs_boundary_step(
        boundary_state,
        event_limits,
        event_counts,
        grid_step
      )
    }

    output[draw, ] <- boundary_state
  }

  output
}

event_counts <- c(
  length(stage_I),
  length(stage_II),
  length(stage_III)
)

total_boundary_draws <- n_joint * draws_per_sample

boundary_draws <- matrix(
  NA_real_,
  nrow = total_boundary_draws,
  ncol = 4L,
  dimnames = list(
    NULL,
    c(
      "Cemetery Start",
      "Stage I to II",
      "Stage II to III",
      "Cemetery End"
    )
  )
)

sample_index_output <- integer(total_boundary_draws)
chain_output <- integer(total_boundary_draws)
retained_iteration_output <- integer(total_boundary_draws)

output_row <- 0L

for (sample_index in seq_len(n_joint)) {
  values_I <- as.numeric(
    unlist(joint[sample_index, stage_I, drop = FALSE], use.names = FALSE)
  )
  values_II <- as.numeric(
    unlist(joint[sample_index, stage_II, drop = FALSE], use.names = FALSE)
  )
  values_III <- as.numeric(
    unlist(joint[sample_index, stage_III, drop = FALSE], use.names = FALSE)
  )

  event_limits <- c(
    min(values_I),
    max(values_I),
    min(values_II),
    max(values_II),
    min(values_III),
    max(values_III)
  )

  reconstructed <- draw_boundaries_for_event_state(
    event_limits = event_limits,
    event_counts = event_counts,
    burn = gibbs_burn,
    thin = gibbs_thin,
    number_of_draws = draws_per_sample,
    grid_step = boundary_grid_step
  )

  rows <- output_row + seq_len(draws_per_sample)
  boundary_draws[rows, ] <- reconstructed
  sample_index_output[rows] <- sample_index

  chain_output[rows] <- if ("chain" %in% names(joint)) {
    as.integer(joint$chain[sample_index])
  } else {
    NA_integer_
  }

  retained_iteration_output[rows] <- if (
    "retained_iteration" %in% names(joint)
  ) {
    as.integer(joint$retained_iteration[sample_index])
  } else {
    sample_index
  }

  output_row <- output_row + draws_per_sample

  if (
    sample_index %% max(1L, floor(n_joint / 20L)) == 0L ||
    sample_index == n_joint
  ) {
    cat(
      "Boundary reconstruction: ",
      sample_index,
      "/",
      n_joint,
      "\n",
      sep = ""
    )
  }
}

boundary_sample_table <- data.frame(
  source_sample_index = sample_index_output,
  chain = chain_output,
  retained_iteration = retained_iteration_output,
  cemetery_start_calendar_year = boundary_draws[, "Cemetery Start"],
  stage_I_to_II_calendar_year = boundary_draws[, "Stage I to II"],
  stage_II_to_III_calendar_year = boundary_draws[, "Stage II to III"],
  cemetery_end_calendar_year = boundary_draws[, "Cemetery End"],
  stringsAsFactors = FALSE
)

boundary_sample_table$stage_I_duration_years <- (
  boundary_sample_table$stage_I_to_II_calendar_year -
    boundary_sample_table$cemetery_start_calendar_year
)

boundary_sample_table$stage_II_duration_years <- (
  boundary_sample_table$stage_II_to_III_calendar_year -
    boundary_sample_table$stage_I_to_II_calendar_year
)

boundary_sample_table$stage_III_duration_years <- (
  boundary_sample_table$cemetery_end_calendar_year -
    boundary_sample_table$stage_II_to_III_calendar_year
)

boundary_sample_table$total_cemetery_duration_years <- (
  boundary_sample_table$cemetery_end_calendar_year -
    boundary_sample_table$cemetery_start_calendar_year
)

write.csv(
  boundary_sample_table,
  file.path(output_dir, "boundary_and_stage_duration_samples.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ----------------------------------------------------------------------
# 4. Probability-density helpers
# ----------------------------------------------------------------------

normalize_mass <- function(mass) {
  mass[!is.finite(mass) | mass < 0] <- 0
  total <- sum(mass)

  if (!is.finite(total) || total <= 0) {
    stop("Cannot normalize zero probability mass.", call. = FALSE)
  }

  mass / total
}

weighted_mean <- function(x, mass) {
  mass <- normalize_mass(mass)
  sum(x * mass)
}

weighted_quantile <- function(x, mass, probabilities) {
  order_index <- order(x)
  x <- x[order_index]
  mass <- normalize_mass(mass[order_index])
  cumulative <- cumsum(mass)

  vapply(
    probabilities,
    function(probability) {
      x[which(cumulative >= probability)[1L]]
    },
    numeric(1L)
  )
}

sample_to_fixed_grid <- function(
  samples,
  bin_width = 5,
  center_offset = 0.5
) {
  samples <- samples[is.finite(samples)]

  if (length(samples) == 0L) {
    stop("No finite samples for density grid.", call. = FALSE)
  }

  first_center <- (
    floor((min(samples) - center_offset) / bin_width) * bin_width +
      center_offset
  )

  last_center <- (
    ceiling((max(samples) - center_offset) / bin_width) * bin_width +
      center_offset
  )

  centers <- seq(first_center, last_center, by = bin_width)
  breaks <- c(
    centers - bin_width / 2,
    tail(centers, 1L) + bin_width / 2
  )

  counts <- hist(
    samples,
    breaks = breaks,
    plot = FALSE,
    include.lowest = TRUE,
    right = FALSE
  )$counts

  mass <- normalize_mass(counts)
  density <- mass / bin_width

  data.frame(
    value = centers,
    bin_lower = centers - bin_width / 2,
    bin_upper = centers + bin_width / 2,
    bin_width = bin_width,
    probability_density = density,
    probability_mass = mass,
    stringsAsFactors = FALSE
  )
}

hpd_grid_intervals <- function(
  values,
  density,
  mass,
  bin_width,
  level,
  parameter_name,
  parameter_type,
  scale
) {
  mass <- normalize_mass(mass)
  density[!is.finite(density)] <- 0

  ranking <- order(density, decreasing = TRUE)
  cumulative <- cumsum(mass[ranking])
  selected_count <- which(cumulative >= level)[1L]
  selected <- ranking[seq_len(selected_count)]

  selected_values <- sort(values[selected])

  if (length(selected_values) == 0L) {
    return(data.frame())
  }

  gaps <- c(Inf, diff(selected_values))
  group <- cumsum(gaps > bin_width * 1.5)

  groups <- split(selected_values, group)

  do.call(
    rbind,
    lapply(
      seq_along(groups),
      function(interval_number) {
        interval_values <- groups[[interval_number]]
        included <- values %in% interval_values

        data.frame(
          parameter_type = parameter_type,
          name = parameter_name,
          scale = scale,
          level = level,
          interval = interval_number,
          lower = min(interval_values) - bin_width / 2,
          upper = max(interval_values) + bin_width / 2,
          interval_probability = sum(mass[included]),
          stringsAsFactors = FALSE
        )
      }
    )
  )
}

summarize_grid_distribution <- function(
  density_table,
  parameter_type,
  parameter_name,
  stage,
  scale
) {
  x <- density_table$value
  density <- density_table$probability_density
  mass <- normalize_mass(density_table$probability_mass)
  bin_width <- unique(density_table$bin_width)

  if (length(bin_width) != 1L) {
    stop("A distribution contains more than one bin width.", call. = FALSE)
  }

  quantiles <- weighted_quantile(
    x,
    mass,
    c(0.023, 0.159, 0.5, 0.841, 0.977)
  )

  hpd_683 <- hpd_grid_intervals(
    values = x,
    density = density,
    mass = mass,
    bin_width = bin_width,
    level = 0.683,
    parameter_name = parameter_name,
    parameter_type = parameter_type,
    scale = scale
  )

  hpd_954 <- hpd_grid_intervals(
    values = x,
    density = density,
    mass = mass,
    bin_width = bin_width,
    level = 0.954,
    parameter_name = parameter_name,
    parameter_type = parameter_type,
    scale = scale
  )

  envelope_683 <- c(
    min(hpd_683$lower),
    max(hpd_683$upper)
  )

  envelope_954 <- c(
    min(hpd_954$lower),
    max(hpd_954$upper)
  )

  if (scale == "cal BP") {
    from_683 <- envelope_683[2L]
    to_683 <- envelope_683[1L]
    from_954 <- envelope_954[2L]
    to_954 <- envelope_954[1L]
  } else {
    from_683 <- envelope_683[1L]
    to_683 <- envelope_683[2L]
    from_954 <- envelope_954[1L]
    to_954 <- envelope_954[2L]
  }

  summary_row <- data.frame(
    parameter_type = parameter_type,
    name = parameter_name,
    stage = stage,
    scale = scale,
    mu = weighted_mean(x, mass),
    from_68_3 = from_683,
    to_68_3 = to_683,
    from_95_4 = from_954,
    to_95_4 = to_954,
    q02_3 = quantiles[1L],
    q15_9 = quantiles[2L],
    median = quantiles[3L],
    q84_1 = quantiles[4L],
    q97_7 = quantiles[5L],
    hpd_68_3_number_of_intervals = nrow(hpd_683),
    hpd_95_4_number_of_intervals = nrow(hpd_954),
    stringsAsFactors = FALSE
  )

  list(
    summary = summary_row,
    hpd = rbind(hpd_683, hpd_954)
  )
}

# ----------------------------------------------------------------------
# 5. Read updated burial-event posterior densities
# ----------------------------------------------------------------------

burial_density <- read.csv(
  burial_density_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

required_density_columns <- c(
  "burial",
  "stage",
  "calendar_year",
  "bin_width",
  "kinship_posterior_density",
  "kinship_posterior_probability_mass"
)

missing_density_columns <- setdiff(
  required_density_columns,
  names(burial_density)
)

if (length(missing_density_columns) > 0L) {
  stop(
    "墓葬概率密度文件缺少字段：",
    paste(missing_density_columns, collapse = ", "),
    call. = FALSE
  )
}

burial_density$calendar_year <- as.numeric(
  burial_density$calendar_year
)

burial_density$cal_BP <- 1950 - burial_density$calendar_year
burial_density$bin_width <- as.numeric(burial_density$bin_width)
burial_density$probability_density <- as.numeric(
  burial_density$kinship_posterior_density
)
burial_density$probability_mass <- as.numeric(
  burial_density$kinship_posterior_probability_mass
)

# ----------------------------------------------------------------------
# 6. Build full density table: burials + boundaries + durations
# ----------------------------------------------------------------------

all_density_tables <- list()
summary_rows <- list()
hpd_rows <- list()
density_counter <- 0L
summary_counter <- 0L
hpd_counter <- 0L

append_distribution <- function(
  table,
  parameter_type,
  parameter_name,
  stage,
  scale
) {
  density_counter <<- density_counter + 1L
  summary_counter <<- summary_counter + 1L

  table$parameter_type <- parameter_type
  table$name <- parameter_name
  table$stage <- stage
  table$scale <- scale

  all_density_tables[[density_counter]] <<- table[, c(
    "parameter_type",
    "name",
    "stage",
    "scale",
    "value",
    "bin_lower",
    "bin_upper",
    "bin_width",
    "probability_density",
    "probability_mass"
  )]

  summarized <- summarize_grid_distribution(
    density_table = table,
    parameter_type = parameter_type,
    parameter_name = parameter_name,
    stage = stage,
    scale = scale
  )

  summary_rows[[summary_counter]] <<- summarized$summary

  if (nrow(summarized$hpd) > 0L) {
    hpd_counter <<- hpd_counter + 1L
    summarized$hpd$stage <- stage
    hpd_rows[[hpd_counter]] <<- summarized$hpd
  }
}

# Burial distributions in cal BP.
for (burial in required_burials) {
  d <- burial_density[
    burial_density$burial == burial,
    ,
    drop = FALSE
  ]

  if (nrow(d) == 0L) {
    stop("墓葬概率密度中缺少：", burial, call. = FALSE)
  }

  d <- d[order(d$cal_BP), , drop = FALSE]

  table <- data.frame(
    value = d$cal_BP,
    bin_lower = d$cal_BP - d$bin_width / 2,
    bin_upper = d$cal_BP + d$bin_width / 2,
    bin_width = d$bin_width,
    probability_density = d$probability_density,
    probability_mass = normalize_mass(d$probability_mass),
    stringsAsFactors = FALSE
  )

  append_distribution(
    table,
    parameter_type = "Burial",
    parameter_name = burial,
    stage = unname(stage_map[[burial]]),
    scale = "cal BP"
  )
}

# Boundary distributions, converted from calendar year to cal BP.
boundary_variables <- list(
  "Cemetery Start" =
    boundary_sample_table$cemetery_start_calendar_year,
  "Stage I to II" =
    boundary_sample_table$stage_I_to_II_calendar_year,
  "Stage II to III" =
    boundary_sample_table$stage_II_to_III_calendar_year,
  "Cemetery End" =
    boundary_sample_table$cemetery_end_calendar_year
)

boundary_stage_labels <- c(
  "Cemetery Start" = "Before Stage I",
  "Stage I to II" = "Stage I / Stage II",
  "Stage II to III" = "Stage II / Stage III",
  "Cemetery End" = "After Stage III"
)

for (boundary_name in names(boundary_variables)) {
  cal_BP_samples <- 1950 - boundary_variables[[boundary_name]]

  table <- sample_to_fixed_grid(
    cal_BP_samples,
    bin_width = 5,
    center_offset = 0.5
  )

  append_distribution(
    table,
    parameter_type = "Boundary",
    parameter_name = boundary_name,
    stage = unname(boundary_stage_labels[[boundary_name]]),
    scale = "cal BP"
  )
}

duration_variables <- list(
  "Stage I Duration" =
    boundary_sample_table$stage_I_duration_years,
  "Stage II Duration" =
    boundary_sample_table$stage_II_duration_years,
  "Stage III Duration" =
    boundary_sample_table$stage_III_duration_years
)

for (duration_name in names(duration_variables)) {
  table <- sample_to_fixed_grid(
    duration_variables[[duration_name]],
    bin_width = 5,
    center_offset = 2.5
  )

  stage_name <- sub(" Duration$", "", duration_name)

  append_distribution(
    table,
    parameter_type = "Interval",
    parameter_name = duration_name,
    stage = stage_name,
    scale = "years"
  )
}

full_density <- do.call(rbind, all_density_tables)
model_summary <- do.call(rbind, summary_rows)
hpd_intervals <- do.call(rbind, hpd_rows)

model_summary$order <- match(model_summary$name, parameter_order)
model_summary <- model_summary[
  order(model_summary$order),
  ,
  drop = FALSE
]
model_summary$order <- NULL

full_density$order <- match(full_density$name, parameter_order)
full_density <- full_density[
  order(full_density$order, full_density$value),
  ,
  drop = FALSE
]
full_density$order <- NULL

hpd_intervals$order <- match(hpd_intervals$name, parameter_order)
hpd_intervals <- hpd_intervals[
  order(
    hpd_intervals$order,
    hpd_intervals$level,
    hpd_intervals$interval
  ),
  ,
  drop = FALSE
]
hpd_intervals$order <- NULL

# ----------------------------------------------------------------------
# 7. Comparison table
# ----------------------------------------------------------------------

# This corrected reconstruction is based on the kinship-updated,
# globally ordered joint posterior. The original OxCal boundary posterior
# cannot be recovered exactly from marginal CSV exports, so no artificial
# "original vs updated boundary" comparison is constructed here.
comparison_summary <- model_summary

# ----------------------------------------------------------------------
# 8. Save outputs
# ----------------------------------------------------------------------

# Save an explicit check that every source joint sample obeyed global order.
stage_order_check <- data.frame(
  joint_samples_checked = nrow(joint),
  stage_I_before_II_fraction = mean(stage_I_latest < stage_II_earliest),
  stage_II_before_III_fraction = mean(stage_II_latest < stage_III_earliest),
  complete_global_order_fraction = mean(global_order_valid),
  minimum_stage_I_to_II_gap_years = min(
    stage_II_earliest - stage_I_latest
  ),
  minimum_stage_II_to_III_gap_years = min(
    stage_III_earliest - stage_II_latest
  ),
  stringsAsFactors = FALSE
)

write.csv(
  stage_order_check,
  file.path(output_dir, "source_global_stage_order_check.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  source_metadata,
  file.path(output_dir, "source_model_metadata_copy.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  model_summary,
  file.path(output_dir, "kinship_updated_oxcal_like_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  comparison_summary,
  file.path(output_dir, "kinship_updated_oxcal_like_summary_copy.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  full_density,
  file.path(output_dir, "kinship_updated_all_parameter_density.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  hpd_intervals,
  file.path(output_dir, "kinship_updated_hpd_intervals.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

duration_summary <- model_summary[
  model_summary$parameter_type == "Interval",
  ,
  drop = FALSE
]

write.csv(
  duration_summary,
  file.path(output_dir, "stage_duration_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

duration_density <- full_density[
  full_density$parameter_type == "Interval",
  ,
  drop = FALSE
]

write.csv(
  duration_density,
  file.path(output_dir, "stage_duration_probability_density.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

boundary_summary <- model_summary[
  model_summary$parameter_type == "Boundary",
  ,
  drop = FALSE
]

write.csv(
  boundary_summary,
  file.path(output_dir, "boundary_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Normalization check.
normalization_check <- aggregate(
  probability_mass ~ parameter_type + name + scale,
  data = full_density,
  FUN = sum
)

names(normalization_check)[
  names(normalization_check) == "probability_mass"
] <- "probability_mass_sum"

normalization_check$absolute_error_from_one <- abs(
  normalization_check$probability_mass_sum - 1
)

write.csv(
  normalization_check,
  file.path(output_dir, "density_normalization_check.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

run_metadata <- data.frame(
  item = c(
    "source_results_dir",
    "joint_samples_file",
    "burial_density_file",
    "burials_used_file",
    "source_metadata_file",
    "seed",
    "joint_samples_used",
    "boundary_draws_per_sample",
    "total_boundary_draws",
    "boundary_gibbs_burn",
    "boundary_gibbs_thin",
    "boundary_grid_step_years",
    "boundary_model",
    "calendar_conversion",
    "summary_interval_method"
  ),
  value = c(
    normalizePath(
      results_dir,
      winslash = "/",
      mustWork = FALSE
    ),
    normalizePath(
      joint_samples_file,
      winslash = "/",
      mustWork = FALSE
    ),
    normalizePath(
      burial_density_file,
      winslash = "/",
      mustWork = FALSE
    ),
    normalizePath(
      burials_used_file,
      winslash = "/",
      mustWork = FALSE
    ),
    normalizePath(
      source_metadata_file,
      winslash = "/",
      mustWork = FALSE
    ),
    as.character(seed),
    as.character(n_joint),
    as.character(draws_per_sample),
    as.character(total_boundary_draws),
    as.character(gibbs_burn),
    as.character(gibbs_thin),
    as.character(boundary_grid_step),
    paste0(
      "Conditional uniform-Phase boundary reconstruction; ",
      "phase likelihood proportional to duration^(-n)"
    ),
    "cal BP = 1950 - calendar_year",
    paste0(
      "68.3% and 95.4% fixed-grid HPD; summary table uses ",
      "the envelope when HPD is disjoint"
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  run_metadata,
  file.path(output_dir, "run_metadata.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nCompleted.\n")
cat("Output directory: ", output_dir, "\n", sep = "")
cat(
  "Main summary: ",
  file.path(output_dir, "kinship_updated_oxcal_like_summary.csv"),
  "\n",
  sep = ""
)
cat(
  "Stage durations: ",
  file.path(output_dir, "stage_duration_summary.csv"),
  "\n",
  sep = ""
)
cat(
  "All densities: ",
  file.path(output_dir, "kinship_updated_all_parameter_density.csv"),
  "\n",
  sep = ""
)
