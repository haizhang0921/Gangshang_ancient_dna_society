#!/usr/bin/env Rscript

# ============================================================================
# South Cemetery — conditional reconstruction of Cemetery Start, Cemetery End
# and Cemetery Duration from the corrected kinship-updated modular chronology
#
# Input:
#   south_cemetery_kinship_results/
#       joint_posterior_samples.csv
#       burial_posterior_density.csv
#       burials_used.csv
#       run_metadata.csv
#
#   gangshang_s.csv
#
# Statistical model
# -----------------
# For the South Cemetery burial set with n top-level burial events, the SM14-SM9
# archaeological ordering is already encoded in the OxCal marginal posteriors.
# It is therefore not applied a second time in this modular boundary step.
# Conditional on burial dates, the enclosing Phase-boundary contribution is:
#
# For a single Phase with n burial events:
#
#     B0 < all burial events < B1
#
# and uniform event intensity within the Phase,
#
#     L(B0,B1 | events) proportional to (B1-B0)^(-n).
#
# Conditional on a retained joint burial state, let:
#
#     a = earliest burial date
#     b = latest burial date
#     R = b-a
#
# Define external extensions:
#
#     e0 = a-B0 >= 0
#     e1 = B1-b >= 0
#     S  = e0+e1
#
# Then the exact conditional law is:
#
#     X = S/(R+S) ~ Beta(2, n-2)
#     S = R*X/(1-X)
#     e0 | S ~ Uniform(0,S)
#     e1 = S-e0
#
# Hence:
#
#     B0 = a-e0
#     B1 = b+e1
#     D  = B1-B0 = R+S
#
# No within-sample Gibbs sampler is required.
#
# IMPORTANT
# ---------
# This is a conditional modular reconstruction, not an exact recovery of the
# original OxCal joint Boundary posterior, because the OxCal CSV contains only
# marginal posterior distributions.
#
# Usage:
#   source("reconstruct_south_cemetery_boundaries_duration.R")
#
# ============================================================================


# ---------------------------------------------------------------------------
# 1. Arguments
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

results_dir <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "south_R_exact_main"
}

oxcal_file <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  "gangshang_s.csv"
}

output_dir <- if (length(args) >= 3L) {
  args[[3L]]
} else {
  "south_R_boundaries_main"
}

seed <- if (length(args) >= 4L) {
  as.integer(args[[4L]])
} else {
  20260920L
}

if (!is.finite(seed)) {
  seed <- 20260920L
}

set.seed(seed)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

joint_file <- file.path(
  results_dir,
  "joint_posterior_samples.csv"
)

density_file <- file.path(
  results_dir,
  "burial_posterior_density.csv"
)

burials_file <- file.path(
  results_dir,
  "burials_used.csv"
)

metadata_file <- file.path(
  results_dir,
  "run_metadata.csv"
)

required_files <- c(
  joint_file,
  density_file,
  burials_file,
  metadata_file,
  oxcal_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}


# ---------------------------------------------------------------------------
# 2. Helpers
# ---------------------------------------------------------------------------

normalize_mass <- function(mass) {

  mass[!is.finite(mass) | mass < 0] <- 0

  s <- sum(mass)

  if (!is.finite(s) || s <= 0) {
    stop("Cannot normalize probability mass.", call. = FALSE)
  }

  mass / s
}


weighted_mean <- function(x, mass) {

  mass <- normalize_mass(mass)
  sum(x * mass)
}


weighted_quantile <- function(x, mass, probabilities) {

  o <- order(x)

  x <- x[o]
  mass <- normalize_mass(mass[o])

  cumulative <- cumsum(mass)

  vapply(
    probabilities,
    function(p) {
      x[which(cumulative >= p)[1L]]
    },
    numeric(1L)
  )
}


sample_to_fixed_grid <- function(
  samples,
  bin_width = 5,
  center_offset = 0.5
) {

  samples <- samples[
    is.finite(samples)
  ]

  if (length(samples) == 0L) {
    stop("No finite samples.", call. = FALSE)
  }

  first_center <- (
    floor(
      (min(samples) - center_offset) /
        bin_width
    ) *
      bin_width +
      center_offset
  )

  last_center <- (
    ceiling(
      (max(samples) - center_offset) /
        bin_width
    ) *
      bin_width +
      center_offset
  )

  centers <- seq(
    first_center,
    last_center,
    by = bin_width
  )

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

  data.frame(
    value = centers,
    bin_lower = centers - bin_width / 2,
    bin_upper = centers + bin_width / 2,
    bin_width = bin_width,
    probability_density = mass / bin_width,
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
  parameter_type,
  parameter_name,
  scale
) {

  mass <- normalize_mass(mass)

  density[
    !is.finite(density) |
      density < 0
  ] <- 0

  ranking <- order(
    density,
    decreasing = TRUE
  )

  cumulative <- cumsum(
    mass[ranking]
  )

  selected_count <- which(
    cumulative >= level
  )[1L]

  selected <- ranking[
    seq_len(selected_count)
  ]

  selected_values <- sort(
    values[selected]
  )

  gaps <- c(
    Inf,
    diff(selected_values)
  )

  group <- cumsum(
    gaps >
      bin_width * 1.5
  )

  groups <- split(
    selected_values,
    group
  )

  do.call(
    rbind,
    lapply(
      seq_along(groups),
      function(interval_number) {

        interval_values <- groups[[interval_number]]

        included <- values %in%
          interval_values

        data.frame(
          parameter_type = parameter_type,
          name = parameter_name,
          scale = scale,
          level = level,
          interval = interval_number,
          lower = min(interval_values) -
            bin_width / 2,
          upper = max(interval_values) +
            bin_width / 2,
          interval_probability = sum(
            mass[included]
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )
}


summarize_grid <- function(
  table,
  parameter_type,
  parameter_name,
  scale
) {

  x <- table$value

  mass <- normalize_mass(
    table$probability_mass
  )

  density <- table$probability_density

  bw <- unique(
    table$bin_width
  )

  if (length(bw) != 1L) {
    stop("Multiple bin widths detected.", call. = FALSE)
  }

  q <- weighted_quantile(
    x,
    mass,
    c(
      0.023,
      0.159,
      0.5,
      0.841,
      0.977
    )
  )

  h68 <- hpd_grid_intervals(
    x,
    density,
    mass,
    bw,
    0.683,
    parameter_type,
    parameter_name,
    scale
  )

  h95 <- hpd_grid_intervals(
    x,
    density,
    mass,
    bw,
    0.954,
    parameter_type,
    parameter_name,
    scale
  )

  envelope68 <- c(
    min(h68$lower),
    max(h68$upper)
  )

  envelope95 <- c(
    min(h95$lower),
    max(h95$upper)
  )

  if (scale == "cal BP") {

    from68 <- envelope68[2L]
    to68 <- envelope68[1L]

    from95 <- envelope95[2L]
    to95 <- envelope95[1L]

  } else {

    from68 <- envelope68[1L]
    to68 <- envelope68[2L]

    from95 <- envelope95[1L]
    to95 <- envelope95[2L]
  }

  list(
    summary = data.frame(
      parameter_type = parameter_type,
      name = parameter_name,
      scale = scale,
      mu = weighted_mean(x, mass),
      from_68_3 = from68,
      to_68_3 = to68,
      from_95_4 = from95,
      to_95_4 = to95,
      q02_3 = q[1L],
      q15_9 = q[2L],
      median = q[3L],
      q84_1 = q[4L],
      q97_7 = q[5L],
      hpd_68_3_number_of_intervals = nrow(h68),
      hpd_95_4_number_of_intervals = nrow(h95),
      stringsAsFactors = FALSE
    ),

    hpd = rbind(
      h68,
      h95
    )
  )
}


# ---------------------------------------------------------------------------
# 3. Validate source model
# ---------------------------------------------------------------------------

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

if (!all(
  c("setting", "value") %in%
    names(metadata)
)) {
  stop(
    "Unexpected run_metadata.csv structure.",
    call. = FALSE
  )
}

metadata_value <- setNames(
  as.character(metadata$value),
  as.character(metadata$setting)
)

if (
  metadata_value[["n_top_level_burial_events"]] != "14"
) {
  stop(
    "Expected 14 burial events.",
    call. = FALSE
  )
}

if (
  metadata_value[["n_kinship_relationships"]] != "1"
) {
  stop(
    "Expected one kinship relationship.",
    call. = FALSE
  )
}

if (
  metadata_value[["kinship_pair"]] != "SM1--SM2"
) {
  stop(
    "Expected SM1--SM2 kinship pair.",
    call. = FALSE
  )
}

gap_mean_value <- suppressWarnings(as.numeric(metadata_value[["gap_mean_years"]]))
gap_sd_value <- suppressWarnings(as.numeric(metadata_value[["gap_sd_years"]]))
if (!is.finite(gap_mean_value) || !is.finite(gap_sd_value) || gap_sd_value <= 0) {
  stop("Invalid SM1-SM2 relationship gap parameters in run_metadata.csv.", call. = FALSE)
}

if (
  metadata_value[["direction_established"]] != "FALSE"
) {
  stop(
    "SM1-SM2 must be undirected.",
    call. = FALSE
  )
}


# ---------------------------------------------------------------------------
# 4. Read burial structure and joint draws
# ---------------------------------------------------------------------------

burials <- read.csv(
  burials_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

if (!all(
  c(
    "burial",
    "phase"
  ) %in%
    names(burials)
)) {
  stop(
    "burials_used.csv has an unexpected structure.",
    call. = FALSE
  )
}

burial_names <- as.character(
  burials$burial
)

if (
  length(burial_names) != 14L ||
  anyDuplicated(burial_names)
) {
  stop(
    "Expected 14 unique burial events.",
    call. = FALSE
  )
}


joint <- read.csv(
  joint_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

missing_joint <- setdiff(
  burial_names,
  names(joint)
)

if (length(missing_joint) > 0L) {
  stop(
    "joint_posterior_samples.csv is missing burial(s): ",
    paste(
      missing_joint,
      collapse = ", "
    ),
    call. = FALSE
  )
}

joint[
  burial_names
] <- lapply(
  joint[
    burial_names
  ],
  as.numeric
)

bad_rows <- !apply(
  joint[
    burial_names
  ],
  1L,
  function(x) {
    all(
      is.finite(x)
    )
  }
)

if (any(bad_rows)) {

  warning(
    sum(bad_rows),
    " invalid joint sample row(s) removed."
  )

  joint <- joint[
    !bad_rows,
    ,
    drop = FALSE
  ]
}

if (nrow(joint) == 0L) {
  stop(
    "No valid joint posterior samples.",
    call. = FALSE
  )
}


# ---------------------------------------------------------------------------
# 5. Exact conditional single-Phase boundary reconstruction
# ---------------------------------------------------------------------------

n_events <- length(
  burial_names
)

if (n_events <= 2L) {
  stop(
    "The flat-boundary uniform-Phase posterior requires n > 2.",
    call. = FALSE
  )
}

n_draws <- nrow(
  joint
)

boundary_samples <- data.frame(
  source_sample_index = seq_len(
    n_draws
  ),
  cemetery_start_calendar_year = NA_real_,
  cemetery_end_calendar_year = NA_real_,
  cemetery_duration_years = NA_real_,
  observed_burial_span_years = NA_real_,
  start_extension_years = NA_real_,
  end_extension_years = NA_real_,
  stringsAsFactors = FALSE
)


for (i in seq_len(n_draws)) {

  event_dates <- as.numeric(
    unlist(
      joint[
        i,
        burial_names,
        drop = FALSE
      ],
      use.names = FALSE
    )
  )

  earliest <- min(
    event_dates
  )

  latest <- max(
    event_dates
  )

  R <- latest -
    earliest

  if (
    !is.finite(R) ||
    R <= 0
  ) {
    stop(
      "Invalid burial span in joint sample ",
      i,
      ".",
      call. = FALSE
    )
  }


  # Exact conditional draw.
  x <- stats::rbeta(
    1L,
    shape1 = 2,
    shape2 = n_events - 2
  )

  # Guard against machine-level x = 1.
  x <- min(
    x,
    1 - .Machine$double.eps
  )

  total_external_extension <- (
    R *
      x /
      (1 - x)
  )

  start_extension <- stats::runif(
    1L,
    min = 0,
    max = total_external_extension
  )

  end_extension <- (
    total_external_extension -
      start_extension
  )

  B0 <- earliest -
    start_extension

  B1 <- latest +
    end_extension

  D <- B1 -
    B0


  boundary_samples$cemetery_start_calendar_year[i] <- B0
  boundary_samples$cemetery_end_calendar_year[i] <- B1
  boundary_samples$cemetery_duration_years[i] <- D
  boundary_samples$observed_burial_span_years[i] <- R
  boundary_samples$start_extension_years[i] <- start_extension
  boundary_samples$end_extension_years[i] <- end_extension
}


# ---------------------------------------------------------------------------
# 6. Validate reconstructed boundary states
# ---------------------------------------------------------------------------

earliest_by_sample <- apply(
  joint[
    burial_names
  ],
  1L,
  min
)

latest_by_sample <- apply(
  joint[
    burial_names
  ],
  1L,
  max
)

valid_start <- (
  boundary_samples$cemetery_start_calendar_year <=
    earliest_by_sample
)

valid_end <- (
  boundary_samples$cemetery_end_calendar_year >=
    latest_by_sample
)

valid_duration <- (
  boundary_samples$cemetery_duration_years >
    0
)

if (
  !all(
    valid_start &
      valid_end &
      valid_duration
  )
) {
  stop(
    "At least one reconstructed boundary state is invalid.",
    call. = FALSE
  )
}


write.csv(
  boundary_samples,
  file.path(
    output_dir,
    "boundary_and_cemetery_duration_samples.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ---------------------------------------------------------------------------
# 7. Build boundary/duration density distributions
# ---------------------------------------------------------------------------

start_calBP <- 1950 -
  boundary_samples$cemetery_start_calendar_year

end_calBP <- 1950 -
  boundary_samples$cemetery_end_calendar_year

duration <- boundary_samples$cemetery_duration_years


start_density <- sample_to_fixed_grid(
  start_calBP,
  bin_width = 5,
  center_offset = 0.5
)

end_density <- sample_to_fixed_grid(
  end_calBP,
  bin_width = 5,
  center_offset = 0.5
)

duration_density <- sample_to_fixed_grid(
  duration,
  bin_width = 5,
  center_offset = 0.5
)


start_summary <- summarize_grid(
  start_density,
  parameter_type = "Boundary",
  parameter_name = "Cemetery Start",
  scale = "cal BP"
)

end_summary <- summarize_grid(
  end_density,
  parameter_type = "Boundary",
  parameter_name = "Cemetery End",
  scale = "cal BP"
)

duration_summary <- summarize_grid(
  duration_density,
  parameter_type = "Duration",
  parameter_name = "Cemetery Duration",
  scale = "years"
)


final_summary <- rbind(
  start_summary$summary,
  end_summary$summary,
  duration_summary$summary
)

final_hpd <- rbind(
  start_summary$hpd,
  end_summary$hpd,
  duration_summary$hpd
)


write.csv(
  final_summary,
  file.path(
    output_dir,
    "boundary_duration_summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  final_hpd,
  file.path(
    output_dir,
    "boundary_duration_hpd_intervals.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


add_labels <- function(
  table,
  parameter_type,
  parameter_name,
  scale
) {

  table$parameter_type <- parameter_type
  table$name <- parameter_name
  table$scale <- scale

  table[
    ,
    c(
      "parameter_type",
      "name",
      "scale",
      "value",
      "bin_lower",
      "bin_upper",
      "bin_width",
      "probability_density",
      "probability_mass"
    )
  ]
}


boundary_duration_density <- rbind(
  add_labels(
    start_density,
    "Boundary",
    "Cemetery Start",
    "cal BP"
  ),
  add_labels(
    end_density,
    "Boundary",
    "Cemetery End",
    "cal BP"
  ),
  add_labels(
    duration_density,
    "Duration",
    "Cemetery Duration",
    "years"
  )
)


write.csv(
  boundary_duration_density,
  file.path(
    output_dir,
    "boundary_duration_probability_density.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ---------------------------------------------------------------------------
# 8. Add all 14 burial distributions to form the final 17-parameter chronology
# ---------------------------------------------------------------------------

burial_density <- read.csv(
  density_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

required_density_columns <- c(
  "burial",
  "calendar_year",
  "bin_width",
  "kinship_posterior_density",
  "kinship_posterior_probability_mass"
)

missing_density <- setdiff(
  required_density_columns,
  names(burial_density)
)

if (length(missing_density) > 0L) {
  stop(
    "burial_posterior_density.csv is missing: ",
    paste(
      missing_density,
      collapse = ", "
    ),
    call. = FALSE
  )
}


all_density <- list()
all_summary <- list()
all_hpd <- list()

counter <- 0L

for (b in burial_names) {

  d <- burial_density[
    burial_density$burial == b,
    ,
    drop = FALSE
  ]

  if (nrow(d) == 0L) {
    stop(
      "Missing burial density for ",
      b,
      ".",
      call. = FALSE
    )
  }

  d$cal_BP <- 1950 -
    as.numeric(
      d$calendar_year
    )

  table <- data.frame(
    value = d$cal_BP,
    bin_lower = d$cal_BP -
      d$bin_width / 2,
    bin_upper = d$cal_BP +
      d$bin_width / 2,
    bin_width = d$bin_width,
    probability_density = d$kinship_posterior_density,
    probability_mass = normalize_mass(
      d$kinship_posterior_probability_mass
    ),
    stringsAsFactors = FALSE
  )

  s <- summarize_grid(
    table,
    parameter_type = "Burial",
    parameter_name = b,
    scale = "cal BP"
  )

  counter <- counter + 1L

  all_density[[counter]] <- add_labels(
    table,
    "Burial",
    b,
    "cal BP"
  )

  all_summary[[counter]] <- s$summary
  all_hpd[[counter]] <- s$hpd
}


counter <- counter + 1L
all_density[[counter]] <- add_labels(
  start_density,
  "Boundary",
  "Cemetery Start",
  "cal BP"
)
all_summary[[counter]] <- start_summary$summary
all_hpd[[counter]] <- start_summary$hpd

counter <- counter + 1L
all_density[[counter]] <- add_labels(
  end_density,
  "Boundary",
  "Cemetery End",
  "cal BP"
)
all_summary[[counter]] <- end_summary$summary
all_hpd[[counter]] <- end_summary$hpd

counter <- counter + 1L
all_density[[counter]] <- add_labels(
  duration_density,
  "Duration",
  "Cemetery Duration",
  "years"
)
all_summary[[counter]] <- duration_summary$summary
all_hpd[[counter]] <- duration_summary$hpd


final_all_density <- do.call(
  rbind,
  all_density
)

final_all_summary <- do.call(
  rbind,
  all_summary
)

final_all_hpd <- do.call(
  rbind,
  all_hpd
)


write.csv(
  final_all_density,
  file.path(
    output_dir,
    "south_cemetery_final_17_parameter_density.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  final_all_summary,
  file.path(
    output_dir,
    "south_cemetery_final_17_parameter_summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  final_all_hpd,
  file.path(
    output_dir,
    "south_cemetery_final_17_parameter_hpd_intervals.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ---------------------------------------------------------------------------
# 9. Density normalization
# ---------------------------------------------------------------------------

normalization <- aggregate(
  cbind(
    probability_mass,
    density_integral =
      final_all_density$probability_density *
        final_all_density$bin_width
  ) ~ parameter_type +
    name +
    scale,
  data = final_all_density,
  FUN = sum
)

normalization$mass_error <- abs(
  normalization$probability_mass -
    1
)

normalization$density_error <- abs(
  normalization$density_integral -
    1
)


write.csv(
  normalization,
  file.path(
    output_dir,
    "density_normalization_check.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ---------------------------------------------------------------------------
# 10. Original OxCal archaeology-only boundary distributions for comparison
# ---------------------------------------------------------------------------

ox <- read.csv(
  oxcal_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

ox$name <- trimws(
  as.character(
    ox$name
  )
)

ox$op <- trimws(
  as.character(
    ox$op
  )
)

ox$type <- trimws(
  as.character(
    ox$type
  )
)

ox$value <- as.numeric(
  ox$value
)

ox$probability <- as.numeric(
  ox$probability
)


original_boundary_summary <- list()

for (boundary_name in c(
  "Cemetery Start",
  "Cemetery End"
)) {

  d <- ox[
    ox$op == "Boundary" &
      ox$name == boundary_name &
      ox$type == "posterior",
    ,
    drop = FALSE
  ]

  if (nrow(d) == 0L) {
    next
  }

  calBP <- 1950 -
    d$value

  # OxCal exported probability column behaves as a fixed-grid density.
  mass <- d$probability * 5
  mass <- normalize_mass(
    mass
  )

  table <- data.frame(
    value = calBP,
    bin_lower = calBP - 2.5,
    bin_upper = calBP + 2.5,
    bin_width = 5,
    probability_density = mass / 5,
    probability_mass = mass,
    stringsAsFactors = FALSE
  )

  original_boundary_summary[[boundary_name]] <- summarize_grid(
    table,
    parameter_type = "Original OxCal Boundary",
    parameter_name = boundary_name,
    scale = "cal BP"
  )$summary
}


if (length(original_boundary_summary) > 0L) {

  original_boundary_summary <- do.call(
    rbind,
    original_boundary_summary
  )

  write.csv(
    original_boundary_summary,
    file.path(
      output_dir,
      "original_oxcal_boundary_summary.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}


# ---------------------------------------------------------------------------
# 11. Metadata
# ---------------------------------------------------------------------------

output_metadata <- data.frame(
  setting = c(
    "source_results_dir",
    "source_model_version",
    "oxcal_file",
    "n_burial_events",
    "source_joint_samples",
    "boundary_draws",
    "seed",
    "boundary_model",
    "boundary_sampling_method",
    "statistical_note"
  ),

  value = c(
    normalizePath(
      results_dir,
      winslash = "/",
      mustWork = FALSE
    ),
    metadata_value[["model_version"]],
    normalizePath(
      oxcal_file,
      winslash = "/",
      mustWork = FALSE
    ),
    as.character(
      n_events
    ),
    as.character(
      nrow(joint)
    ),
    as.character(
      nrow(boundary_samples)
    ),
    as.character(
      seed
    ),
    "Single unordered uniform Phase; likelihood proportional to duration^(-n)",
    "Exact conditional Beta transformation; no Gibbs sampler",
    paste(
      "Conditional modular reconstruction from kinship-updated burial chronology.",
      "Not an exact recovery of the original OxCal joint Boundary posterior because",
      "the exported OxCal file contains marginal posterior distributions only."
    )
  ),

  stringsAsFactors = FALSE
)


write.csv(
  output_metadata,
  file.path(
    output_dir,
    "run_metadata.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)




# ---------------------------------------------------------------------------
# 11b. Publication-scale cal BCE summary files
# ---------------------------------------------------------------------------
# Internal calculations above retain the original cal BP convention of the
# validated script. For publication, calendar boundaries are additionally
# exported as cal BCE (= cal BP - 1950). Durations remain unchanged.

if (file.exists(file.path(output_dir, "boundary_duration_summary.csv"))) {
  pub <- read.csv(file.path(output_dir, "boundary_duration_summary.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
  is_boundary <- pub$parameter_type == "Boundary"
  for (cc in c("mu", "from_68_3", "to_68_3", "from_95_4", "to_95_4",
               "q02_3", "q15_9", "median", "q84_1", "q97_7")) {
    if (cc %in% names(pub)) {
      pub[[cc]][is_boundary] <- pub[[cc]][is_boundary] - 1950
    }
  }
  pub$scale[is_boundary] <- "cal BCE"
  write.csv(pub, file.path(output_dir, "boundary_duration_summary_calBCE.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}

if (file.exists(file.path(output_dir, "original_oxcal_boundary_summary.csv"))) {
  pubox <- read.csv(file.path(output_dir, "original_oxcal_boundary_summary.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  for (cc in c("mu", "from_68_3", "to_68_3", "from_95_4", "to_95_4",
               "q02_3", "q15_9", "median", "q84_1", "q97_7")) {
    if (cc %in% names(pubox)) pubox[[cc]] <- pubox[[cc]] - 1950
  }
  pubox$scale <- "cal BCE"
  write.csv(pubox, file.path(output_dir, "original_oxcal_boundary_summary_calBCE.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}


# ---------------------------------------------------------------------------
# 12. Console summary
# ---------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("South Cemetery boundary/duration reconstruction completed\n")
cat("============================================================\n\n")

cat(
  "Burial events: ",
  n_events,
  "\n",
  sep = ""
)

cat(
  "Joint burial states: ",
  nrow(joint),
  "\n",
  sep = ""
)

cat(
  "Boundary states: ",
  nrow(boundary_samples),
  "\n",
  sep = ""
)

cat(
  "Boundary sampler: exact conditional Beta transformation\n\n"
)

cat(
  "Mean observed burial span: ",
  round(
    mean(
      boundary_samples$observed_burial_span_years
    ),
    2
  ),
  " years\n",
  sep = ""
)

cat(
  "Mean reconstructed cemetery duration: ",
  round(
    mean(
      boundary_samples$cemetery_duration_years
    ),
    2
  ),
  " years\n",
  sep = ""
)

cat(
  "Median reconstructed cemetery duration: ",
  round(
    median(
      boundary_samples$cemetery_duration_years
    ),
    2
  ),
  " years\n\n",
  sep = ""
)

cat(
  "Maximum final probability-mass error: ",
  format(
    max(
      normalization$mass_error
    ),
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Maximum final density-integral error: ",
  format(
    max(
      normalization$density_error
    ),
    scientific = TRUE
  ),
  "\n\n",
  sep = ""
)

cat(
  "Output directory: ",
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat("============================================================\n")
