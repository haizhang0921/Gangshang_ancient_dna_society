#!/usr/bin/env Rscript

# ============================================================================
# South Cemetery — corrected modular Bayesian chronology
#
# Archaeology-only OxCal marginal posteriors
#        +
# one ancient-DNA kinship likelihood (SM1-SM2)
#
# Authoritative inputs
# --------------------
#   gangshang_s.csv
#   genetic_relations_s.xlsx
#
# Scenario relationship
# ---------------------
#   SM1 -- SM2, undirected in both scenarios
#   Main: 2nd relationship, 35 +/- 26 years
#   Alternative: sibling, 26 +/- 22 years
#
# Archaeological information already encoded in the OxCal marginals includes
# the SM14-SM9 Sequence and the SM1 top-level Combine event.
#
# IMPORTANT
# ---------
# The statistical logic is the same as for the Northern Cemetery:
#
#   posterior ∝ archaeological OxCal marginal information
#               × kinship temporal-gap likelihood
#
# However, the South Cemetery has only one kinship edge. Therefore the
# two-dimensional posterior for SM1 and SM2 can be evaluated EXACTLY on
# the OxCal calendar grid. No MCMC approximation is needed for the kinship
# update itself.
#
# The script additionally generates 12,000 draws from the exact modular
# joint chronology for later conditional reconstruction of Cemetery Start,
# Cemetery End and Cemetery Duration.
#
# Usage in RStudio
# ----------------
#   setwd("D:/YOUR_PROJECT_FOLDER")
#   source("south_cemetery_kinship_model_exact.R")
#
# Required package
# ----------------
#   install.packages("readxl")
# ============================================================================


# ---------------------------------------------------------------------------
# 1. Arguments and settings
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

oxcal_file <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "gangshang_s.csv"
}

relationship_file <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  "genetic_relations_s.csv"
}

out_dir <- if (length(args) >= 3L) {
  args[[3L]]
} else {
  "south_R_exact_main"
}

master_seed <- if (length(args) >= 4L) {
  as.integer(args[[4L]])
} else {
  20260920L
}

n_joint_draws <- if (length(args) >= 5L) {
  as.integer(args[[5L]])
} else {
  12000L
}

if (!is.finite(master_seed)) {
  master_seed <- 20260920L
}

if (!is.finite(n_joint_draws) || n_joint_draws < 1000L) {
  stop("n_joint_draws must be at least 1000.", call. = FALSE)
}

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

set.seed(master_seed)


# ---------------------------------------------------------------------------
# 2. Package
# ---------------------------------------------------------------------------

# Relationship tables can be supplied as CSV (recommended, no extra package)
# or XLSX. XLSX input requires the readxl package.
read_relationship_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop(
        "XLSX relationship input requires package 'readxl'.\n",
        "Install it with install.packages('readxl'), or use the supplied CSV file.",
        call. = FALSE
      )
    }
    return(as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE))
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}


# ---------------------------------------------------------------------------
# 3. Helper functions
# ---------------------------------------------------------------------------

weighted_mean_safe <- function(x, w) {

  keep <- is.finite(x) & is.finite(w) & w >= 0

  x <- x[keep]
  w <- w[keep]

  if (length(x) == 0L || sum(w) <= 0) {
    return(NA_real_)
  }

  sum(x * w) / sum(w)
}


weighted_sd_safe <- function(x, w) {

  keep <- is.finite(x) & is.finite(w) & w >= 0

  x <- x[keep]
  w <- w[keep]

  if (length(x) == 0L || sum(w) <= 0) {
    return(NA_real_)
  }

  w <- w / sum(w)
  mu <- sum(x * w)

  sqrt(sum(w * (x - mu)^2))
}


weighted_quantile_discrete <- function(x, w, probs) {

  keep <- is.finite(x) & is.finite(w) & w >= 0

  x <- x[keep]
  w <- w[keep]

  if (length(x) == 0L || sum(w) <= 0) {
    return(rep(NA_real_, length(probs)))
  }

  o <- order(x)

  x <- x[o]
  w <- w[o] / sum(w)

  cw <- cumsum(w)

  vapply(
    probs,
    function(p) {
      x[which(cw >= p)[1L]]
    },
    numeric(1L)
  )
}


mode_discrete <- function(x, w) {

  if (length(x) == 0L || sum(w, na.rm = TRUE) <= 0) {
    return(NA_real_)
  }

  x[which.max(w)]
}


infer_bin_width <- function(x) {

  ux <- sort(unique(x[is.finite(x)]))

  if (length(ux) < 2L) {
    return(NA_real_)
  }

  diffs <- diff(ux)
  diffs <- diffs[is.finite(diffs) & diffs > 0]

  if (length(diffs) == 0L) {
    return(NA_real_)
  }

  stats::median(diffs)
}


normalize_density <- function(year, density) {

  keep <- is.finite(year) & is.finite(density) & density >= 0

  year <- year[keep]
  density <- density[keep]

  o <- order(year)
  year <- year[o]
  density <- density[o]

  bw <- infer_bin_width(year)

  if (!is.finite(bw) || bw <= 0) {
    stop("Could not infer a positive calendar-grid width.", call. = FALSE)
  }

  integral <- sum(density * bw)

  if (!is.finite(integral) || integral <= 0) {
    stop("A posterior density has zero or invalid total probability.", call. = FALSE)
  }

  density <- density / integral
  mass <- density * bw

  list(
    year = year,
    bin_width = bw,
    density = density,
    mass = mass
  )
}


hpd_intervals_discrete <- function(x, mass, level, bin_width) {

  keep <- is.finite(x) & is.finite(mass) & mass >= 0

  x <- x[keep]
  mass <- mass[keep]

  if (length(x) == 0L || sum(mass) <= 0) {
    return(
      data.frame(
        lower = numeric(0),
        upper = numeric(0),
        probability = numeric(0)
      )
    )
  }

  mass <- mass / sum(mass)

  ord_prob <- order(
    mass,
    decreasing = TRUE
  )

  cumulative <- cumsum(
    mass[ord_prob]
  )

  k <- which(
    cumulative >= level
  )[1L]

  selected_idx <- ord_prob[
    seq_len(k)
  ]

  sx <- sort(
    x[selected_idx]
  )

  if (length(sx) == 1L) {

    return(
      data.frame(
        lower = sx,
        upper = sx,
        probability = sum(
          mass[x == sx]
        )
      )
    )
  }

  break_after <- which(
    diff(sx) >
      bin_width * 1.01
  )

  group_start <- c(
    1L,
    break_after + 1L
  )

  group_end <- c(
    break_after,
    length(sx)
  )

  out <- vector(
    "list",
    length(group_start)
  )

  for (g in seq_along(group_start)) {

    vals <- sx[
      group_start[g]:
        group_end[g]
    ]

    in_group <- x %in% vals

    out[[g]] <- data.frame(
      lower = min(vals),
      upper = max(vals),
      probability = sum(
        mass[in_group]
      )
    )
  }

  do.call(rbind, out)
}


hpd_envelope <- function(x, mass, level, bin_width) {

  h <- hpd_intervals_discrete(
    x,
    mass,
    level,
    bin_width
  )

  if (nrow(h) == 0L) {
    return(c(NA_real_, NA_real_))
  }

  c(
    min(h$lower),
    max(h$upper)
  )
}


format_hpd_intervals <- function(h) {

  if (nrow(h) == 0L) {
    return("")
  }

  paste(
    paste0(
      format(
        h$lower,
        trim = TRUE,
        scientific = FALSE
      ),
      " to ",
      format(
        h$upper,
        trim = TRUE,
        scientific = FALSE
      )
    ),
    collapse = "; "
  )
}


# ---------------------------------------------------------------------------
# 4. Read and validate archaeology-only OxCal export
# ---------------------------------------------------------------------------

if (!file.exists(oxcal_file)) {
  stop(
    "OxCal file not found: ",
    oxcal_file,
    call. = FALSE
  )
}

ox <- read.csv(
  oxcal_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

required_oxcal_columns <- c(
  "index",
  "op",
  "name",
  "type",
  "value",
  "probability"
)

missing_oxcal_columns <- setdiff(
  required_oxcal_columns,
  names(ox)
)

if (length(missing_oxcal_columns) > 0L) {
  stop(
    "OxCal CSV is missing column(s): ",
    paste(
      missing_oxcal_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

ox$index <- suppressWarnings(
  as.integer(ox$index)
)

ox$op <- trimws(
  as.character(ox$op)
)

ox$name <- trimws(
  as.character(ox$name)
)

ox$type <- trimws(
  as.character(ox$type)
)

ox$value <- suppressWarnings(
  as.numeric(ox$value)
)

ox$probability <- suppressWarnings(
  as.numeric(ox$probability)
)


# This must be an archaeology-only export.
# The current archaeology-only OxCal export may legitimately contain the
# Cemetery Duration Interval. This is an archaeological model output and is
# NOT a kinship constraint. We only reject legacy N() pedigree constraints.
if (any(
  ox$op %in% c("N"),
  na.rm = TRUE
)) {
  stop(
    "The OxCal export still contains legacy N() constraints. ",
    "Use the archaeology-only South Cemetery export.",
    call. = FALSE
  )
}

if (any(
  grepl(
    "GenM2-M1",
    ox$name,
    fixed = TRUE
  ),
  na.rm = TRUE
)) {
  stop(
    "The OxCal export still contains GenM2-M1. ",
    "Use the archaeology-only South Cemetery export.",
    call. = FALSE
  )
}


boundary_parameters <- unique(
  ox[
    ox$type == "posterior" &
      ox$op == "Boundary",
    c("index", "name"),
    drop = FALSE
  ]
)

expected_boundaries <- c(
  "Cemetery Start",
  "Cemetery End"
)

if (
  nrow(boundary_parameters) != 2L ||
  !setequal(
    boundary_parameters$name,
    expected_boundaries
  )
) {
  stop(
    "Expected exactly two Phase boundaries: Cemetery Start and Cemetery End.",
    call. = FALSE
  )
}


# ---------------------------------------------------------------------------
# 5. Identify top-level burial-event parameters
# ---------------------------------------------------------------------------

posterior_parameters <- unique(
  ox[
    ox$type == "posterior" &
      ox$op %in% c(
        "R_Date",
        "Combine"
      ),
    c(
      "index",
      "op",
      "name"
    ),
    drop = FALSE
  ]
)

combine_names <- posterior_parameters$name[
  posterior_parameters$op == "Combine"
]

is_combine_component <- rep(
  FALSE,
  nrow(posterior_parameters)
)

if (length(combine_names) > 0L) {

  # IMPORTANT:
  # Do NOT use startsWith(name, combine_name) here.
  #
  # For Combine("SM1"), startsWith() would incorrectly classify ordinary
  # burials SM10, SM11, SM12, SM13, SM14, SM15 and SM17 as components of
  # SM1. In the South Cemetery export, the true component determinations
  # are named SM1R1, SM1R2, SM1R3 and SM1R4.
  #
  # We therefore require the exact component naming pattern:
  #
  #   <Combine name>R<integer>
  #
  # e.g. SM1R1 ... SM1R4.
  for (cn in combine_names) {

    component_pattern <- paste0(
      "^",
      cn,
      "R[0-9]+$"
    )

    is_combine_component <- (
      is_combine_component |
        (
          posterior_parameters$op == "R_Date" &
          grepl(
            component_pattern,
            posterior_parameters$name
          )
        )
    )
  }
}

burial_parameters <- posterior_parameters[
  !is_combine_component,
  ,
  drop = FALSE
]

burial_parameters <- burial_parameters[
  order(burial_parameters$index),
  ,
  drop = FALSE
]

burial_names <- burial_parameters$name


# Dataset-specific checks.
if (length(burial_names) != 14L) {
  stop(
    "Expected 14 top-level South Cemetery burial events, but detected ",
    length(burial_names),
    ".",
    call. = FALSE
  )
}

if (!all(
  c("SM1", "SM2") %in%
    burial_names
)) {
  stop(
    "SM1 and/or SM2 are missing from the top-level burial-event parameters.",
    call. = FALSE
  )
}

if (!identical(
  burial_parameters$op[
    burial_parameters$name == "SM1"
  ],
  "Combine"
)) {
  stop(
    "SM1 should be the top-level Combine burial event.",
    call. = FALSE
  )
}


cat(
  "Detected top-level South Cemetery burial events (",
  length(burial_names),
  "): ",
  paste(burial_names, collapse = ", "),
  "\n\n",
  sep = ""
)


# Count measurement-level R_Date parameters, including SM1 components.
measurement_parameters <- unique(
  ox[
    ox$type == "posterior" &
      ox$op == "R_Date",
    c(
      "index",
      "name"
    ),
    drop = FALSE
  ]
)

if (nrow(measurement_parameters) != 17L) {
  warning(
    "Detected ",
    nrow(measurement_parameters),
    " R_Date parameters; expected 17 from the current South Cemetery export."
  )
}


# ---------------------------------------------------------------------------
# 6. Extract and normalize burial-level archaeological posteriors
# ---------------------------------------------------------------------------

archaeological <- setNames(
  vector(
    "list",
    length(burial_names)
  ),
  burial_names
)

for (b in burial_names) {

  d <- ox[
    ox$type == "posterior" &
      ox$name == b &
      ox$op ==
        burial_parameters$op[
          burial_parameters$name == b
        ],
    c(
      "value",
      "probability"
    ),
    drop = FALSE
  ]

  archaeological[[b]] <- normalize_density(
    year = d$value,
    density = d$probability
  )
}


all_bin_widths <- vapply(
  archaeological,
  function(x) x$bin_width,
  numeric(1L)
)

if (
  any(
    abs(
      all_bin_widths -
        5
    ) >
      1e-6
  )
) {
  stop(
    "Expected a 5-year OxCal calendar grid for every burial.",
    call. = FALSE
  )
}


# ---------------------------------------------------------------------------
# 7. Read and validate the single kinship relationship
# ---------------------------------------------------------------------------

if (!file.exists(relationship_file)) {
  stop(
    "Relationship file not found: ",
    relationship_file,
    call. = FALSE
  )
}

rel_raw <- read_relationship_table(relationship_file)

required_relationship_columns <- c(
  "Earlier death",
  "Later death",
  "Relationship",
  "Mean",
  "SD",
  "Direction established"
)

missing_relationship_columns <- setdiff(
  required_relationship_columns,
  names(rel_raw)
)

if (length(missing_relationship_columns) > 0L) {
  stop(
    "Relationship workbook is missing column(s): ",
    paste(
      missing_relationship_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (nrow(rel_raw) != 1L) {
  stop(
    "Expected exactly one South Cemetery kinship relationship.",
    call. = FALSE
  )
}

burial_a <- trimws(
  as.character(
    rel_raw[["Earlier death"]][1L]
  )
)

burial_b <- trimws(
  as.character(
    rel_raw[["Later death"]][1L]
  )
)

relationship_type <- trimws(
  as.character(
    rel_raw[["Relationship"]][1L]
  )
)

relationship_mean <- as.numeric(
  rel_raw[["Mean"]][1L]
)

relationship_sd <- as.numeric(
  rel_raw[["SD"]][1L]
)

direction_text <- tolower(
  trimws(
    as.character(
      rel_raw[["Direction established"]][1L]
    )
  )
)

direction_established <- (
  direction_text %in%
    c(
      "yes",
      "y",
      "true",
      "1"
    )
)


if (!setequal(
  c(
    burial_a,
    burial_b
  ),
  c(
    "SM1",
    "SM2"
  )
)) {
  stop(
    "Expected the single kinship relation to connect SM1 and SM2.",
    call. = FALSE
  )
}

if (
  !relationship_type %in% c("2nd relationship", "Sibling") ||
  !is.finite(relationship_mean) ||
  !is.finite(relationship_sd) ||
  relationship_sd <= 0 ||
  direction_established
) {
  stop(
    "Expected one UNDIRECTED SM1-SM2 relationship with a valid Mean and SD. ",
    "Allowed scenario labels are '2nd relationship' or 'Sibling'.",
    call. = FALSE
  )
}

# Scenario-specific parameter audit. This prevents accidental swapping of
# relationship labels and temporal-gap distributions.
if (relationship_type == "2nd relationship" &&
    (relationship_mean != 35 || relationship_sd != 26)) {
  stop("For '2nd relationship', expected Mean=35 and SD=26 years.", call. = FALSE)
}
if (relationship_type == "Sibling" &&
    (relationship_mean != 26 || relationship_sd != 22)) {
  stop("For 'Sibling', expected Mean=26 and SD=22 years.", call. = FALSE)
}


# ---------------------------------------------------------------------------
# 8. Exact two-dimensional Bayesian update for SM1 and SM2
# ---------------------------------------------------------------------------

A <- archaeological[[burial_a]]
B <- archaeological[[burial_b]]

year_a <- A$year
year_b <- B$year

mass_a <- A$mass / sum(A$mass)
mass_b <- B$mass / sum(B$mass)


# Matrix of absolute time gaps.
gap_matrix <- abs(
  outer(
    year_a,
    year_b,
    "-"
  )
)


# Soft relationship likelihood.
relationship_likelihood <- stats::dnorm(
  gap_matrix,
  mean = relationship_mean,
  sd = relationship_sd
)


# Modular archaeological joint information.
archaeological_joint <- outer(
  mass_a,
  mass_b,
  "*"
)


# Exact integrated joint posterior.
joint_unnormalized <- (
  archaeological_joint *
    relationship_likelihood
)

joint_constant <- sum(
  joint_unnormalized
)

if (
  !is.finite(joint_constant) ||
  joint_constant <= 0
) {
  stop(
    "The exact SM1-SM2 joint posterior has zero or invalid normalization.",
    call. = FALSE
  )
}

joint_posterior <- (
  joint_unnormalized /
    joint_constant
)


# Updated marginals.
updated_mass_a <- rowSums(
  joint_posterior
)

updated_mass_b <- colSums(
  joint_posterior
)


# ---------------------------------------------------------------------------
# 9. Build complete burial posterior-density output
# ---------------------------------------------------------------------------

burial_density_list <- vector(
  "list",
  length(burial_names)
)

names(
  burial_density_list
) <- burial_names


for (b in burial_names) {

  arch <- archaeological[[b]]

  kinship_mass <- arch$mass

  related_degree <- 0L

  if (b == burial_a) {

    kinship_mass <- updated_mass_a
    related_degree <- 1L

  } else if (b == burial_b) {

    kinship_mass <- updated_mass_b
    related_degree <- 1L
  }

  kinship_mass <- (
    kinship_mass /
      sum(kinship_mass)
  )

  kinship_density <- (
    kinship_mass /
      arch$bin_width
  )

  burial_density_list[[b]] <- data.frame(
    burial = b,
    phase = "South Cemetery",
    calendar_year = arch$year,
    bin_width = arch$bin_width,
    oxcal_probability_mass = arch$mass,
    oxcal_density = arch$density,
    kinship_posterior_probability_mass = kinship_mass,
    kinship_posterior_density = kinship_density,
    related_degree = related_degree,
    stringsAsFactors = FALSE
  )
}

burial_density <- do.call(
  rbind,
  burial_density_list
)

rownames(
  burial_density
) <- NULL


write.csv(
  burial_density,
  file.path(
    out_dir,
    "burial_posterior_density.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ---------------------------------------------------------------------------
# 10. Burial posterior summaries and HPD intervals
# ---------------------------------------------------------------------------

summary_rows <- vector(
  "list",
  length(burial_names)
)

hpd_rows <- list()
hpd_counter <- 0L


for (i in seq_along(burial_names)) {

  b <- burial_names[i]

  d <- burial_density[
    burial_density$burial == b,
    ,
    drop = FALSE
  ]

  x <- d$calendar_year
  bw <- unique(d$bin_width)[1L]

  ox_mass <- d$oxcal_probability_mass
  kin_mass <- d$kinship_posterior_probability_mass

  ox_mass <- ox_mass / sum(ox_mass)
  kin_mass <- kin_mass / sum(kin_mass)

  ox_q <- weighted_quantile_discrete(
    x,
    ox_mass,
    c(
      0.025,
      0.5,
      0.975
    )
  )

  kin_q <- weighted_quantile_discrete(
    x,
    kin_mass,
    c(
      0.025,
      0.5,
      0.975
    )
  )

  ox_h68 <- hpd_intervals_discrete(
    x,
    ox_mass,
    0.683,
    bw
  )

  ox_h95 <- hpd_intervals_discrete(
    x,
    ox_mass,
    0.954,
    bw
  )

  kin_h68 <- hpd_intervals_discrete(
    x,
    kin_mass,
    0.683,
    bw
  )

  kin_h95 <- hpd_intervals_discrete(
    x,
    kin_mass,
    0.954,
    bw
  )

  ox_h95_env <- hpd_envelope(
    x,
    ox_mass,
    0.954,
    bw
  )

  kin_h95_env <- hpd_envelope(
    x,
    kin_mass,
    0.954,
    bw
  )

  ox_mean <- weighted_mean_safe(
    x,
    ox_mass
  )

  kin_mean <- weighted_mean_safe(
    x,
    kin_mass
  )

  summary_rows[[i]] <- data.frame(
    burial = b,
    phase = "South Cemetery",
    related_degree = d$related_degree[1L],

    oxcal_mean = ox_mean,
    oxcal_median = ox_q[2L],
    oxcal_mode = mode_discrete(
      x,
      ox_mass
    ),
    oxcal_q025 = ox_q[1L],
    oxcal_q975 = ox_q[3L],
    oxcal_hpd95_lower_envelope = ox_h95_env[1L],
    oxcal_hpd95_upper_envelope = ox_h95_env[2L],
    oxcal_hpd95_width = diff(ox_h95_env),

    posterior_mean = kin_mean,
    posterior_median = kin_q[2L],
    posterior_mode = mode_discrete(
      x,
      kin_mass
    ),
    posterior_q025 = kin_q[1L],
    posterior_q975 = kin_q[3L],
    posterior_hpd95_lower_envelope = kin_h95_env[1L],
    posterior_hpd95_upper_envelope = kin_h95_env[2L],
    posterior_hpd95_width = diff(kin_h95_env),

    posterior_mean_change_years = (
      kin_mean -
        ox_mean
    ),

    hpd95_width_change_years = (
      diff(kin_h95_env) -
        diff(ox_h95_env)
    ),

    oxcal_hpd68_intervals = format_hpd_intervals(
      ox_h68
    ),

    oxcal_hpd95_intervals = format_hpd_intervals(
      ox_h95
    ),

    posterior_hpd68_intervals = format_hpd_intervals(
      kin_h68
    ),

    posterior_hpd95_intervals = format_hpd_intervals(
      kin_h95
    ),

    stringsAsFactors = FALSE
  )


  add_hpd_rows <- function(
    h,
    source,
    level
  ) {

    if (nrow(h) == 0L) {
      return(NULL)
    }

    data.frame(
      burial = b,
      source = source,
      level = level,
      interval_id = seq_len(
        nrow(h)
      ),
      lower_calendar_year = h$lower,
      upper_calendar_year = h$upper,
      included_probability = h$probability,
      stringsAsFactors = FALSE
    )
  }


  temp_hpd <- list(
    add_hpd_rows(
      ox_h68,
      "OxCal archaeological posterior",
      0.683
    ),
    add_hpd_rows(
      ox_h95,
      "OxCal archaeological posterior",
      0.954
    ),
    add_hpd_rows(
      kin_h68,
      "Kinship-updated posterior",
      0.683
    ),
    add_hpd_rows(
      kin_h95,
      "Kinship-updated posterior",
      0.954
    )
  )

  temp_hpd <- temp_hpd[
    !vapply(
      temp_hpd,
      is.null,
      logical(1L)
    )
  ]

  if (length(temp_hpd) > 0L) {

    for (hh in temp_hpd) {

      hpd_counter <- hpd_counter + 1L
      hpd_rows[[hpd_counter]] <- hh
    }
  }
}


burial_summary <- do.call(
  rbind,
  summary_rows
)

burial_hpd <- do.call(
  rbind,
  hpd_rows
)


write.csv(
  burial_summary,
  file.path(
    out_dir,
    "burial_posterior_summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  burial_hpd,
  file.path(
    out_dir,
    "burial_hpd_intervals.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ---------------------------------------------------------------------------
# 11. Exact relationship-gap posterior
# ---------------------------------------------------------------------------

gap_values <- sort(
  unique(
    as.vector(
      gap_matrix
    )
  )
)

gap_mass <- vapply(
  gap_values,
  function(g) {
    sum(
      joint_posterior[
        abs(
          gap_matrix -
            g
        ) <
          1e-8
      ]
    )
  },
  numeric(1L)
)

gap_mass <- gap_mass / sum(
  gap_mass
)

gap_mean <- weighted_mean_safe(
  gap_values,
  gap_mass
)

gap_sd <- weighted_sd_safe(
  gap_values,
  gap_mass
)

gap_quantiles <- weighted_quantile_discrete(
  gap_values,
  gap_mass,
  c(
    0.025,
    0.5,
    0.975
  )
)

gap_h68 <- hpd_intervals_discrete(
  gap_values,
  gap_mass,
  0.683,
  5
)

gap_h95 <- hpd_intervals_discrete(
  gap_values,
  gap_mass,
  0.954,
  5
)


# Calendar year: numerically smaller = earlier.
sm1_earlier_probability <- sum(
  joint_posterior[
    outer(
      year_a,
      year_b,
      "<"
    )
  ]
)

sm2_earlier_probability <- sum(
  joint_posterior[
    outer(
      year_a,
      year_b,
      ">"
    )
  ]
)

same_grid_probability <- sum(
  joint_posterior[
    outer(
      year_a,
      year_b,
      "=="
    )
  ]
)


relationship_diagnostics <- data.frame(
  burial_a = burial_a,
  burial_b = burial_b,
  relationship = relationship_type,
  direction_established = FALSE,
  input_mean = relationship_mean,
  input_sd = relationship_sd,

  posterior_gap_mean = gap_mean,
  posterior_gap_sd = gap_sd,
  posterior_gap_median = gap_quantiles[2L],
  posterior_gap_q025 = gap_quantiles[1L],
  posterior_gap_q975 = gap_quantiles[3L],

  posterior_gap_hpd68 = format_hpd_intervals(
    gap_h68
  ),

  posterior_gap_hpd95 = format_hpd_intervals(
    gap_h95
  ),

  posterior_mean_standardized_difference = (
    gap_mean -
      relationship_mean
  ) /
    relationship_sd,

  probability_SM1_earlier_than_SM2 = sm1_earlier_probability,
  probability_SM2_earlier_than_SM1 = sm2_earlier_probability,
  probability_same_5yr_grid_cell = same_grid_probability,

  stringsAsFactors = FALSE
)


write.csv(
  relationship_diagnostics,
  file.path(
    out_dir,
    "relationship_posterior_diagnostics.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


gap_density_output <- data.frame(
  gap_years = gap_values,
  probability_mass = gap_mass,
  probability_density = gap_mass / 5,
  stringsAsFactors = FALSE
)

write.csv(
  gap_density_output,
  file.path(
    out_dir,
    "relationship_gap_posterior_density.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ---------------------------------------------------------------------------
# 12. Save exact SM1-SM2 two-dimensional joint posterior
# ---------------------------------------------------------------------------

joint_pair_output <- expand.grid(
  SM1_calendar_year = year_a,
  SM2_calendar_year = year_b,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# expand.grid varies the first factor fastest, which is not the same ordering
# as as.vector(matrix). Build the posterior explicitly by matching coordinates.
joint_pair_output$absolute_gap_years <- abs(
  joint_pair_output$SM1_calendar_year -
    joint_pair_output$SM2_calendar_year
)

joint_pair_output$posterior_probability_mass <- mapply(
  function(a, b) {

    ia <- match(
      a,
      year_a
    )

    ib <- match(
      b,
      year_b
    )

    joint_posterior[
      ia,
      ib
    ]
  },
  joint_pair_output$SM1_calendar_year,
  joint_pair_output$SM2_calendar_year
)

write.csv(
  joint_pair_output,
  file.path(
    out_dir,
    "SM1_SM2_exact_joint_posterior.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ---------------------------------------------------------------------------
# 13. Draw a modular joint burial chronology for later boundary reconstruction
# ---------------------------------------------------------------------------
#
# SM1 and SM2 are drawn jointly from their exact two-dimensional posterior.
# The remaining 12 burials are drawn independently from their archaeology-only
# OxCal marginal posteriors. This is the South-Cemetery analogue of the
# modular approximation used for the Northern Cemetery.
# ---------------------------------------------------------------------------

sampled_pair_linear <- sample.int(
  length(
    joint_posterior
  ),
  size = n_joint_draws,
  replace = TRUE,
  prob = as.vector(
    joint_posterior
  )
)

sampled_pair_subscripts <- arrayInd(
  sampled_pair_linear,
  .dim = dim(
    joint_posterior
  )
)

joint_draws <- data.frame(
  sample_id = seq_len(
    n_joint_draws
  ),
  stringsAsFactors = FALSE
)


for (b in burial_names) {

  if (b == burial_a) {

    joint_draws[[b]] <- year_a[
      sampled_pair_subscripts[
        ,
        1L
      ]
    ]

  } else if (b == burial_b) {

    joint_draws[[b]] <- year_b[
      sampled_pair_subscripts[
        ,
        2L
      ]
    ]

  } else {

    arch <- archaeological[[b]]

    joint_draws[[b]] <- sample(
      arch$year,
      size = n_joint_draws,
      replace = TRUE,
      prob = arch$mass
    )
  }
}


write.csv(
  joint_draws,
  file.path(
    out_dir,
    "joint_posterior_samples.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ---------------------------------------------------------------------------
# 14. Probability-normalization checks
# ---------------------------------------------------------------------------

normalization_rows <- vector(
  "list",
  length(burial_names)
)

for (i in seq_along(burial_names)) {

  b <- burial_names[i]

  d <- burial_density[
    burial_density$burial == b,
    ,
    drop = FALSE
  ]

  normalization_rows[[i]] <- data.frame(
    burial = b,

    oxcal_probability_mass_sum = sum(
      d$oxcal_probability_mass
    ),

    oxcal_density_integral = sum(
      d$oxcal_density *
        d$bin_width
    ),

    kinship_probability_mass_sum = sum(
      d$kinship_posterior_probability_mass
    ),

    kinship_density_integral = sum(
      d$kinship_posterior_density *
        d$bin_width
    ),

    stringsAsFactors = FALSE
  )
}

normalization_check <- do.call(
  rbind,
  normalization_rows
)

normalization_check$max_absolute_error_from_one <- apply(
  abs(
    normalization_check[
      ,
      c(
        "oxcal_probability_mass_sum",
        "oxcal_density_integral",
        "kinship_probability_mass_sum",
        "kinship_density_integral"
      )
    ] -
      1
  ),
  1L,
  max
)


write.csv(
  normalization_check,
  file.path(
    out_dir,
    "density_normalization_check.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ---------------------------------------------------------------------------
# 15. Burial and relationship input records
# ---------------------------------------------------------------------------

burials_used <- burial_parameters

names(
  burials_used
)[
  names(
    burials_used
  ) == "index"
] <- "oxcal_index"

burials_used$phase <- "South Cemetery"

burials_used$related_degree <- ifelse(
  burials_used$name %in%
    c(
      burial_a,
      burial_b
    ),
  1L,
  0L
)

names(
  burials_used
)[
  names(
    burials_used
  ) == "name"
] <- "burial"


write.csv(
  burials_used,
  file.path(
    out_dir,
    "burials_used.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


relationships_used <- data.frame(
  burial_a = burial_a,
  burial_b = burial_b,
  relationship = relationship_type,
  mean = relationship_mean,
  sd = relationship_sd,
  direction_established = FALSE,
  direction_source = "none; absolute temporal difference used",
  stringsAsFactors = FALSE
)

write.csv(
  relationships_used,
  file.path(
    out_dir,
    "relationships_used.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ---------------------------------------------------------------------------
# 16. Run metadata
# ---------------------------------------------------------------------------

metadata <- data.frame(
  setting = c(
    "model",
    "model_version",
    "oxcal_file",
    "relationship_file",
    "n_top_level_burial_events",
    "n_radiocarbon_R_Date_parameters",
    "n_kinship_relationships",
    "kinship_pair",
    "relationship_type",
    "gap_mean_years",
    "gap_sd_years",
    "direction_established",
    "calendar_grid_years",
    "inference_method",
    "exact_joint_grid_cells",
    "joint_posterior_draws",
    "master_seed",
    "statistical_note"
  ),

  value = c(
    "OxCal marginal posteriors x kinship gap likelihood",
    "South Cemetery two-scenario exact single-edge modular model; SM14-SM9 archaeology encoded in OxCal marginals; SM1 Combine retained",
    normalizePath(
      oxcal_file,
      winslash = "/",
      mustWork = FALSE
    ),
    normalizePath(
      relationship_file,
      winslash = "/",
      mustWork = FALSE
    ),
    as.character(
      length(
        burial_names
      )
    ),
    as.character(
      nrow(
        measurement_parameters
      )
    ),
    "1",
    paste0(
      burial_a,
      "--",
      burial_b
    ),
    relationship_type,
    as.character(
      relationship_mean
    ),
    as.character(
      relationship_sd
    ),
    "FALSE",
    "5",
    "Exact two-dimensional grid enumeration for SM1-SM2; independent modular marginals for unconnected burials",
    as.character(
      length(
        joint_posterior
      )
    ),
    as.character(
      n_joint_draws
    ),
    as.character(
      master_seed
    ),
    paste(
      "Modular approximation using burial-level OxCal marginal posteriors.",
      "SM1-SM2 is evaluated exactly under an undirected absolute-gap likelihood.",
      "The SM14-SM9 archaeological Sequence and SM1 Combine are already encoded in the OxCal marginals and are not multiplied again in R. Original full OxCal joint boundary correlations are unavailable in the CSV."
    )
  ),

  stringsAsFactors = FALSE
)


write.csv(
  metadata,
  file.path(
    out_dir,
    "run_metadata.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ---------------------------------------------------------------------------
# 17. Compact model summary
# ---------------------------------------------------------------------------

posterior_change_summary <- burial_summary[
  order(
    -abs(
      burial_summary$posterior_mean_change_years
    )
  ),
  c(
    "burial",
    "related_degree",
    "oxcal_mean",
    "posterior_mean",
    "posterior_mean_change_years",
    "oxcal_hpd95_width",
    "posterior_hpd95_width",
    "hpd95_width_change_years"
  ),
  drop = FALSE
]

write.csv(
  posterior_change_summary,
  file.path(
    out_dir,
    "posterior_change_summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ---------------------------------------------------------------------------
# 18. Console report
# ---------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("South Cemetery exact kinship chronology completed\n")
cat("============================================================\n\n")

cat(
  "Archaeology-only OxCal input: ",
  oxcal_file,
  "\n",
  sep = ""
)

cat(
  "Top-level burial events: ",
  length(
    burial_names
  ),
  "\n",
  sep = ""
)

cat(
  "R_Date parameters: ",
  nrow(
    measurement_parameters
  ),
  "\n",
  sep = ""
)

cat(
  "Kinship relationship: ",
  burial_a,
  " -- ",
  burial_b,
  " (",
  relationship_type,
  ", undirected, ",
  relationship_mean,
  " +/- ",
  relationship_sd,
  " years)\n",
  sep = ""
)

cat(
  "Exact SM1-SM2 joint grid cells: ",
  length(
    joint_posterior
  ),
  "\n",
  sep = ""
)

cat(
  "Joint draws for later boundary reconstruction: ",
  n_joint_draws,
  "\n\n",
  sep = ""
)

cat(
  "Posterior absolute gap mean: ",
  round(
    gap_mean,
    2
  ),
  " years\n",
  sep = ""
)

cat(
  "Posterior absolute gap median: ",
  round(
    gap_quantiles[2L],
    2
  ),
  " years\n",
  sep = ""
)

cat(
  "Standardized difference from input mean: ",
  round(
    relationship_diagnostics$posterior_mean_standardized_difference,
    3
  ),
  " SD\n",
  sep = ""
)

cat(
  "P(SM1 earlier than SM2): ",
  round(
    sm1_earlier_probability,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "P(SM2 earlier than SM1): ",
  round(
    sm2_earlier_probability,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "P(same 5-year grid cell): ",
  round(
    same_grid_probability,
    4
  ),
  "\n\n",
  sep = ""
)

cat(
  "Maximum density-normalization error: ",
  format(
    max(
      normalization_check$max_absolute_error_from_one
    ),
    scientific = TRUE
  ),
  "\n\n",
  sep = ""
)

cat(
  "Output directory: ",
  normalizePath(
    out_dir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)

cat("Key output files:\n")
cat("  burial_posterior_density.csv\n")
cat("  burial_posterior_summary.csv\n")
cat("  burial_hpd_intervals.csv\n")
cat("  relationship_posterior_diagnostics.csv\n")
cat("  relationship_gap_posterior_density.csv\n")
cat("  SM1_SM2_exact_joint_posterior.csv\n")
cat("  joint_posterior_samples.csv\n")
cat("  density_normalization_check.csv\n")
cat("  posterior_change_summary.csv\n")
cat("  run_metadata.csv\n")
cat("============================================================\n")
