#!/usr/bin/env Rscript

# =============================================================================
# Northern Cemetery: two-stage kinship-constrained chronology
# =============================================================================
# Stage 1 (already completed in OxCal 4.4):
#   radiocarbon + archaeological Sequence/Phase/Combine model
#   exported as marginal posterior densities in gangshang_n.csv
#
# THREE-SCENARIO INPUT VERSION (2026-09-20)
# -----------------------
# Authoritative inputs: gangshang_n.csv plus one scenario-specific relationship table.
# The relationship table itself is authoritative for degree and direction.
# Main, closer-kin and more-distant-kin scenarios are run with the same chronology code.
#
# New OxCal export includes three additional human radiocarbon determinations
# (48 human determinations in the northern cemetery); burial-event structure remains 34 events.
#
# Stage 2 (this script):
#   use each burial's OxCal marginal posterior as its archaeological prior,
#   then update all burial dates jointly with the 20 kinship gap constraints.
#
# IMPORTANT STATISTICAL LIMITATION
# --------------------------------
# The OxCal CSV contains marginal posterior densities, not the original joint
# OxCal MCMC sample. Therefore the R model uses the product of burial-level
# OxCal marginal posteriors as a modular approximation to the archaeological
# joint posterior. It does not reconstruct the correlations among dates and
# boundaries in the original OxCal Sequence.
#
# Main model:
#   pi(t_1,...,t_n) proportional to
#       product_i p_OxCal,i(t_i) * product_r L_kinship,r(t_a, t_b)
#
# A burial has exactly one date variable, even when it appears in several
# kinship relationships. The resulting graph may contain branches and cycles.
#
# Outputs include posterior probability density on the original OxCal calendar
# grid, posterior mass, HPD sets, relationship diagnostics, MCMC diagnostics,
# joint samples, and a multi-page comparison PDF.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

oxcal_file <- if (length(args) >= 1L) args[[1L]] else "gangshang_n.csv"
relationship_file <- if (length(args) >= 2L) args[[2L]] else "genetic_relations_n.csv"
out_dir <- if (length(args) >= 3L) args[[3L]] else "north_R_primary"
master_seed <- if (length(args) >= 4L) as.integer(args[[4L]]) else 20260806L

# ------------------------------ Runtime settings -----------------------------
# These can be changed without editing the script, for example:
# N_CHAINS=4 N_ITER=100000 BURNIN=25000 THIN=20 Rscript ...

get_env_int <- function(name, default) {
  x <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(x)) stop("Environment variable ", name, " must be an integer.")
  x
}

get_env_num <- function(name, default) {
  x <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(x)) stop("Environment variable ", name, " must be numeric.")
  x
}

as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y", "是")
}

get_env_flag <- function(name, default) {
  raw <- Sys.getenv(name, unset = if (default) "TRUE" else "FALSE")
  as_flag(raw)
}

n_chains <- get_env_int("N_CHAINS", 4L)
n_iter <- get_env_int("N_ITER", 80000L)
burnin <- get_env_int("BURNIN", 20000L)
thin <- get_env_int("THIN", 20L)
rb_max_samples <- get_env_int("RB_MAX_SAMPLES", 12000L)
component_redraw_prob <- get_env_num("COMPONENT_REDRAW_PROB", 0.03)
component_shift_prob <- get_env_num("COMPONENT_SHIFT_PROB", 0.20)
max_shift_years <- get_env_int("MAX_SHIFT_YEARS", 50L)
infer_stage_direction <- get_env_flag("INFER_STAGE_DIRECTION", TRUE)
enforce_global_stage_order <- get_env_flag("ENFORCE_GLOBAL_STAGE_ORDER", FALSE)
save_joint_samples <- get_env_flag("SAVE_JOINT_SAMPLES", TRUE)
make_plots <- get_env_flag("MAKE_PLOTS", TRUE)

if (n_chains < 1L) stop("N_CHAINS must be at least 1.")
if (n_iter <= burnin || burnin < 0L || thin < 1L) {
  stop("Require N_ITER > BURNIN >= 0 and THIN >= 1.")
}
if (rb_max_samples < 100L) stop("RB_MAX_SAMPLES must be at least 100.")
if (component_redraw_prob < 0 || component_redraw_prob > 1 ||
    component_shift_prob < 0 || component_shift_prob > 1) {
  stop("Component move probabilities must be between 0 and 1.")
}
if (max_shift_years < 5L) stop("MAX_SHIFT_YEARS must be at least 5.")

n_keep_per_chain <- floor((n_iter - burnin) / thin)
if (n_keep_per_chain < 200L) {
  warning("Fewer than 200 retained samples per chain; increase N_ITER or reduce THIN.")
}

if (is.na(master_seed)) stop("Seed must be a valid integer.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Northern Cemetery kinship chronology\n")
cat("OxCal file: ", oxcal_file, "\n", sep = "")
cat("Relationship file: ", relationship_file, "\n", sep = "")
cat("Output directory: ", out_dir, "\n", sep = "")
cat(sprintf("MCMC: %d chains; %d iterations; burn-in %d; thin %d\n",
            n_chains, n_iter, burnin, thin))
cat("Global Stage-order indicator: ", enforce_global_stage_order, "\n", sep = "")
cat("Infer direction for cross-Stage relationships: ", infer_stage_direction, "\n\n", sep = "")

# ------------------------------- Helper functions ----------------------------

normalize_names <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y <- gsub("[^a-z0-9]+", "_", y)
  y <- gsub("^_+|_+$", "", y)
  y
}

format_cal_year <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    ifelse(x < 0,
           paste0(format(round(abs(x)), scientific = FALSE, trim = TRUE), " cal BCE"),
           paste0(format(round(x), scientific = FALSE, trim = TRUE), " cal CE"))
  )
}

cell_widths <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n == 1L) return(1)
  if (any(diff(x) <= 0)) stop("Calendar grid must be strictly increasing.")
  edges <- c(
    x[[1L]] - (x[[2L]] - x[[1L]]) / 2,
    (x[-n] + x[-1L]) / 2,
    x[[n]] + (x[[n]] - x[[n - 1L]]) / 2
  )
  diff(edges)
}

log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(x - m)))
}

sample_value <- function(values, probabilities = NULL) {
  if (length(values) < 1L) stop("Cannot sample from an empty vector.")
  j <- sample.int(length(values), size = 1L, prob = probabilities)
  values[[j]]
}

weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  x <- x[ok]
  w <- w[ok]
  if (length(x) == 0L || sum(w) <= 0) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1L]], numeric(1L))
}

weighted_summary <- function(grid, mass) {
  mass <- pmax(as.numeric(mass), 0)
  mass <- mass / sum(mass)
  qs <- weighted_quantile(grid, mass, c(0.025, 0.16, 0.5, 0.84, 0.975))
  c(
    mean = sum(grid * mass),
    mode = grid[which.max(mass)],
    q025 = qs[[1L]], q16 = qs[[2L]], median = qs[[3L]],
    q84 = qs[[4L]], q975 = qs[[5L]]
  )
}

hpd_intervals <- function(grid, density, mass, level, burial, stage_label) {
  density <- pmax(as.numeric(density), 0)
  mass <- pmax(as.numeric(mass), 0)
  mass <- mass / sum(mass)
  ord <- order(density, decreasing = TRUE)
  k <- which(cumsum(mass[ord]) >= level)[1L]
  selected <- rep(FALSE, length(grid))
  selected[ord[seq_len(k)]] <- TRUE

  widths <- cell_widths(grid)
  lower <- grid - widths / 2
  upper <- grid + widths / 2
  idx <- which(selected)
  if (length(idx) == 0L) return(NULL)

  groups <- cumsum(c(TRUE, diff(idx) > 1L))
  pieces <- split(idx, groups)
  do.call(rbind, lapply(seq_along(pieces), function(j) {
    ii <- pieces[[j]]
    data.frame(
      burial = burial,
      stage = stage_label,
      level = level,
      interval = j,
      lower_cal_year = min(lower[ii]),
      upper_cal_year = max(upper[ii]),
      interval_probability = sum(mass[ii]),
      lower_label = format_cal_year(min(lower[ii])),
      upper_label = format_cal_year(max(upper[ii])),
      stringsAsFactors = FALSE
    )
  }))
}

read_relationship_table <- function(path) {
  if (!file.exists(path)) stop("Relationship file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop(
        "Reading an Excel relationship file requires the R package 'readxl'.\n",
        "Install it once with: install.packages(\"readxl\")\n",
        "Alternatively save the sheet as CSV and pass that CSV to the script."
      )
    }
    z <- as.data.frame(readxl::read_excel(path, sheet = 1), stringsAsFactors = FALSE)
  } else if (ext %in% c("csv", "txt")) {
    z <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported relationship file type: ", ext, ". Use XLSX, XLS, or CSV.")
  }
  z
}

# ------------------------- 1. Read OxCal posterior CSV -----------------------

if (!file.exists(oxcal_file)) stop("OxCal CSV not found: ", oxcal_file)
ox <- read.csv(oxcal_file, stringsAsFactors = FALSE, check.names = FALSE)
required_ox <- c("index", "op", "name", "type", "value", "probability")
missing_ox <- setdiff(required_ox, names(ox))
if (length(missing_ox) > 0L) {
  stop("OxCal CSV is missing columns: ", paste(missing_ox, collapse = ", "))
}

post <- ox[
  ox$type == "posterior" & ox$op %in% c("Combine", "R_Combine", "R_Date"),
  required_ox,
  drop = FALSE
]
if (nrow(post) == 0L) stop("No posterior Combine/R_Combine/R_Date rows found.")

# Keep top-level burial-event distributions. When a Combine/R_Combine named NM24
# exists, remove component R_Date rows such as NM24-1, NM24-3, etc.
combined_names <- unique(post$name[post$op %in% c("Combine", "R_Combine")])
component_base <- sub("-[0-9]+$", "", post$name)
is_component <- post$op == "R_Date" & grepl("-[0-9]+$", post$name) &
  component_base %in% combined_names
post <- post[!is_component, , drop = FALSE]

post$value <- suppressWarnings(as.numeric(post$value))
post$probability <- suppressWarnings(as.numeric(post$probability))
post <- post[is.finite(post$value) & is.finite(post$probability), , drop = FALSE]
post$probability <- pmax(post$probability, 0)

# Collapse accidental duplicate rows on the same grid point.
post <- aggregate(
  probability ~ index + op + name + type + value,
  data = post,
  FUN = sum
)
post <- post[order(post$name, post$value), , drop = FALSE]

burial_names <- sort(unique(post$name))
if (length(burial_names) < 2L) stop("Fewer than two burial distributions were found.")

base_grid <- base_width <- base_density <- base_mass <- base_log_mass <-
  vector("list", length(burial_names))
names(base_grid) <- names(base_width) <- names(base_density) <-
  names(base_mass) <- names(base_log_mass) <- burial_names
base_op <- setNames(character(length(burial_names)), burial_names)

for (p in burial_names) {
  d <- post[post$name == p, c("op", "value", "probability"), drop = FALSE]
  d <- d[order(d$value), , drop = FALSE]
  grid <- as.numeric(d$value)
  dens <- pmax(as.numeric(d$probability), 0)
  if (anyDuplicated(grid)) stop("Duplicate grid values remain for burial ", p, ".")
  widths <- cell_widths(grid)
  mass <- dens * widths
  if (sum(mass) <= 0) stop("OxCal posterior has zero probability for burial ", p, ".")
  mass <- mass / sum(mass)
  dens <- mass / widths

  base_grid[[p]] <- grid
  base_width[[p]] <- widths
  base_density[[p]] <- dens
  base_mass[[p]] <- mass
  base_log_mass[[p]] <- ifelse(mass > 0, log(mass), -Inf)
  base_op[[p]] <- paste(sort(unique(d$op)), collapse = "+")
}

# Stage membership is reconstructed directly from the OxCal export.
# The top-level burial indices must lie between the named OxCal boundaries:
#   Cemetery Start < Stage I burials < Stage I to II
#   Stage I to II  < Stage II burials < Stage II to III
#   Stage II to III < Stage III burials < Cemetery End
boundary_rows <- ox[
  ox$type == "posterior" & ox$op == "Boundary",
  c("index", "name"),
  drop = FALSE
]
boundary_rows$index <- suppressWarnings(as.integer(boundary_rows$index))
boundary_index <- tapply(boundary_rows$index, boundary_rows$name, function(z) unique(z)[1L])

required_boundaries <- c(
  "Cemetery Start", "Stage I to II", "Stage II to III", "Cemetery End"
)
missing_boundaries <- setdiff(required_boundaries, names(boundary_index))
if (length(missing_boundaries) > 0L) {
  stop(
    "OxCal export is missing required boundary row(s): ",
    paste(missing_boundaries, collapse = ", ")
  )
}

b_start <- as.integer(boundary_index[["Cemetery Start"]])
b_I_II <- as.integer(boundary_index[["Stage I to II"]])
b_II_III <- as.integer(boundary_index[["Stage II to III"]])
b_end <- as.integer(boundary_index[["Cemetery End"]])

if (!(b_start < b_I_II && b_I_II < b_II_III && b_II_III < b_end)) {
  stop("OxCal boundary indices are not in the expected chronological code order.")
}

burial_index <- tapply(
  suppressWarnings(as.integer(post$index)),
  post$name,
  function(z) unique(z)[1L]
)
burial_index <- burial_index[burial_names]

stage_number <- setNames(rep(NA_integer_, length(burial_names)), burial_names)
stage_number[burial_index > b_start & burial_index < b_I_II] <- 1L
stage_number[burial_index > b_I_II & burial_index < b_II_III] <- 2L
stage_number[burial_index > b_II_III & burial_index < b_end] <- 3L

if (anyNA(stage_number)) {
  stop(
    "Could not infer Stage for burial(s): ",
    paste(names(stage_number)[is.na(stage_number)], collapse = ", ")
  )
}

stage_label <- setNames(
  c("Stage I", "Stage II", "Stage III")[stage_number],
  burial_names
)

stage_I <- names(stage_number)[stage_number == 1L]
stage_II <- names(stage_number)[stage_number == 2L]
stage_III <- names(stage_number)[stage_number == 3L]

# Dataset-level integrity checks for the corrected Northern Cemetery input.
if (length(burial_names) != 34L ||
    length(stage_I) != 22L ||
    length(stage_II) != 5L ||
    length(stage_III) != 7L) {
  stop(
    "Unexpected burial/Stage counts. Detected total=", length(burial_names),
    ", Stage I=", length(stage_I),
    ", Stage II=", length(stage_II),
    ", Stage III=", length(stage_III),
    ". Expected 34 / 22 / 5 / 7."
  )
}

cat("Detected ", length(burial_names), " burial-event posterior distributions.\n", sep = "")
cat("Stage I: ", sum(stage_number == 1L), "; Stage II: ",
    sum(stage_number == 2L), "; Stage III: ", sum(stage_number == 3L), "\n\n", sep = "")

# ------------------------- 2. Read relationship data -------------------------

rel_raw <- read_relationship_table(relationship_file)
names(rel_raw) <- normalize_names(names(rel_raw))

required_rel <- c(
  "earlier_death", "later_death", "relationship", "mean", "sd",
  "direction_established"
)
missing_rel <- setdiff(required_rel, names(rel_raw))
if (length(missing_rel) > 0L) {
  stop("Relationship table is missing columns: ", paste(missing_rel, collapse = ", "))
}

rel <- rel_raw[, required_rel, drop = FALSE]
rel <- rel[!is.na(rel$earlier_death) & !is.na(rel$later_death), , drop = FALSE]
rel$earlier_death <- trimws(as.character(rel$earlier_death))
rel$later_death <- trimws(as.character(rel$later_death))
rel$relationship <- trimws(as.character(rel$relationship))
rel$direction_established <- trimws(as.character(rel$direction_established))
rel$mean <- suppressWarnings(as.numeric(rel$mean))
rel$sd <- suppressWarnings(as.numeric(rel$sd))

if (nrow(rel) == 0L) stop("No relationship rows were found.")
if (any(!is.finite(rel$mean)) || any(!is.finite(rel$sd)) ||
    any(rel$mean < 0) || any(rel$sd <= 0)) {
  stop("Every relationship requires Mean >= 0 and SD > 0.")
}
if (any(rel$earlier_death == rel$later_death)) {
  stop("A relationship cannot link a burial to itself.")
}
unknown_burials <- setdiff(
  unique(c(rel$earlier_death, rel$later_death)), burial_names
)
if (length(unknown_burials) > 0L) {
  stop("Relationship burial(s) absent from OxCal output: ",
       paste(unknown_burials, collapse = ", "))
}

rel$id <- sprintf("K%02d", seq_len(nrow(rel)))
rel$direction_yes <- as_flag(rel$direction_established)
rel$burial_a <- rel$earlier_death
rel$burial_b <- rel$later_death
rel$effective_direction <- "unknown"
rel$earlier_burial <- NA_character_
rel$later_burial <- NA_character_
rel$direction_source <- "none"

for (k in seq_len(nrow(rel))) {
  a <- rel$burial_a[[k]]
  b <- rel$burial_b[[k]]
  sa <- stage_number[[a]]
  sb <- stage_number[[b]]

  if (rel$direction_yes[[k]]) {
    rel$effective_direction[[k]] <- "A_before_B"
    rel$earlier_burial[[k]] <- a
    rel$later_burial[[k]] <- b
    rel$direction_source[[k]] <- "Direction established = YES"
    if (sa > sb) {
      stop("Relationship ", rel$id[[k]],
           " has a confirmed direction that conflicts with the OxCal Stage order.")
    }
  } else if (infer_stage_direction && sa != sb) {
    if (sa < sb) {
      rel$effective_direction[[k]] <- "A_before_B"
      rel$earlier_burial[[k]] <- a
      rel$later_burial[[k]] <- b
    } else {
      rel$effective_direction[[k]] <- "B_before_A"
      rel$earlier_burial[[k]] <- b
      rel$later_burial[[k]] <- a
    }
    rel$direction_source[[k]] <- "inferred from OxCal Stage order"
  }
}

# Record the explicitly directed pairs exactly as supplied by the scenario table.
# This object is used only for diagnostics/metadata; it must be defined for
# every scenario, including alt2 where NM30-NM27 is undirected.
observed_directed_pairs <- if (any(rel$direction_yes)) {
  paste(
    rel$burial_a[rel$direction_yes],
    rel$burial_b[rel$direction_yes],
    sep = "->"
  )
} else {
  character(0)
}

# Scenario-level integrity checks.
# All three uncertainty scenarios contain 20 kinship edges. Relationship class,
# mean, SD and explicit direction are read directly from the supplied table.
if (nrow(rel) != 20L) {
  stop("Expected 20 kinship relationships, but detected ", nrow(rel), ".")
}

# Validate the empirical parameterization used throughout the manuscript.
allowed_parameters <- data.frame(
  relationship = c(
    "Parent-Offspring",
    "Sibling",
    "Grandparent-Grandchild",
    "2nd relationship"
  ),
  mean = c(29, 26, 35, 35),
  sd = c(19, 22, 32, 26),
  stringsAsFactors = FALSE
)

for (k in seq_len(nrow(rel))) {
  hit <- which(
    allowed_parameters$relationship == rel$relationship[[k]] &
    allowed_parameters$mean == rel$mean[[k]] &
    allowed_parameters$sd == rel$sd[[k]]
  )
  if (length(hit) != 1L) {
    stop(
      "Unexpected relationship parameterization in row ", k, ": ",
      rel$relationship[[k]], ", Mean=", rel$mean[[k]],
      ", SD=", rel$sd[[k]], "."
    )
  }
}

# NM27-NM21 remains undirected Grandparent-Grandchild in all three scenarios.
nm27_nm21 <- which(
  (rel$burial_a == "NM27" & rel$burial_b == "NM21") |
  (rel$burial_a == "NM21" & rel$burial_b == "NM27")
)
if (length(nm27_nm21) != 1L) {
  stop("Expected exactly one NM27-NM21 relationship.")
}
if (rel$direction_yes[nm27_nm21] ||
    rel$relationship[nm27_nm21] != "Grandparent-Grandchild" ||
    rel$mean[nm27_nm21] != 35 ||
    rel$sd[nm27_nm21] != 32) {
  stop(
    "NM27-NM21 must be an UNDIRECTED Grandparent-Grandchild relationship ",
    "with Mean=35 and SD=32."
  )
}

# NM30-NM27 is Parent-Offspring (29 +/- 19 y); whether its death order is
# established is deliberately scenario-specific and is therefore not hard-coded.
nm30_nm27 <- which(
  (rel$burial_a == "NM30" & rel$burial_b == "NM27") |
  (rel$burial_a == "NM27" & rel$burial_b == "NM30")
)
if (length(nm30_nm27) != 1L) {
  stop("Expected exactly one NM30-NM27 relationship.")
}
if (rel$relationship[nm30_nm27] != "Parent-Offspring" ||
    rel$mean[nm30_nm27] != 29 ||
    rel$sd[nm30_nm27] != 19) {
  stop("NM30-NM27 must be Parent-Offspring with Mean=29 and SD=19.")
}

cat("Loaded ", nrow(rel), " kinship relationships.\n", sep = "")
cat("Direction fixed by relationship table: ", sum(rel$direction_yes), "\n", sep = "")
if (length(observed_directed_pairs) > 0L) {
  cat("Established directions: ", paste(sort(observed_directed_pairs), collapse = ", "), "\n", sep = "")
} else {
  cat("Established directions: none\n")
}
cat(
  "NM30-NM27 explicit direction: ",
  if (rel$direction_yes[nm30_nm27]) "NM30 before NM27" else "undirected",
  "\n", sep = ""
)
cat("NM27-NM21: undirected Grandparent-Grandchild, 35 +/- 32 years.\n")
cat("Additional cross-Stage directions inferred: ",
    sum(rel$direction_source == "inferred from OxCal Stage order"), "\n\n", sep = "")

# --------------------------- 3. Graph and likelihood -------------------------

edge_index_by_burial <- setNames(vector("list", length(burial_names)), burial_names)
for (p in burial_names) {
  edge_index_by_burial[[p]] <- which(rel$burial_a == p | rel$burial_b == p)
}

edge_loglik <- function(k, state) {
  mu <- rel$mean[[k]]
  sig <- rel$sd[[k]]
  if (rel$effective_direction[[k]] == "unknown") {
    gap <- abs(state[[rel$burial_b[[k]]]] - state[[rel$burial_a[[k]]]])
  } else {
    gap <- state[[rel$later_burial[[k]]]] - state[[rel$earlier_burial[[k]]]]
    if (!is.finite(gap) || gap < 0) return(-Inf)
  }
  dnorm(gap, mean = mu, sd = sig, log = TRUE)
}

local_loglik <- function(p, state) {
  idx <- edge_index_by_burial[[p]]
  if (length(idx) == 0L) return(0)
  vals <- vapply(idx, edge_loglik, numeric(1L), state = state)
  sum(vals)
}

total_loglik <- function(state) {
  vals <- vapply(seq_len(nrow(rel)), edge_loglik, numeric(1L), state = state)
  sum(vals)
}

global_stage_valid <- function(state) {
  if (!enforce_global_stage_order) return(TRUE)
  max(state[stage_I]) < min(state[stage_II]) &&
    max(state[stage_II]) < min(state[stage_III])
}

local_stage_valid <- function(p, proposed, state) {
  if (!enforce_global_stage_order) return(TRUE)
  s <- stage_number[[p]]
  if (s == 1L) {
    proposed < min(state[stage_II])
  } else if (s == 2L) {
    proposed > max(state[stage_I]) && proposed < min(state[stage_III])
  } else {
    proposed > max(state[stage_II])
  }
}

# Connected components of the kinship graph. These support occasional block
# proposals that improve movement through correlated posterior regions.
related_nodes <- unique(c(rel$burial_a, rel$burial_b))
adjacency <- setNames(vector("list", length(related_nodes)), related_nodes)
for (k in seq_len(nrow(rel))) {
  a <- rel$burial_a[[k]]
  b <- rel$burial_b[[k]]
  adjacency[[a]] <- unique(c(adjacency[[a]], b))
  adjacency[[b]] <- unique(c(adjacency[[b]], a))
}

components <- list()
unseen <- related_nodes
while (length(unseen) > 0L) {
  seed_node <- unseen[[1L]]
  queue <- seed_node
  comp <- character(0)
  while (length(queue) > 0L) {
    x <- queue[[1L]]
    queue <- queue[-1L]
    if (x %in% comp) next
    comp <- c(comp, x)
    queue <- unique(c(queue, adjacency[[x]]))
  }
  components[[length(components) + 1L]] <- sort(unique(comp))
  unseen <- setdiff(unseen, comp)
}

component_edge_indices <- lapply(components, function(comp) {
  which(rel$burial_a %in% comp | rel$burial_b %in% comp)
})

# -------------------------- 4. Initial-state functions -----------------------

draw_from_prior <- function(p) {
  sample_value(base_grid[[p]], base_mass[[p]])
}

draw_independent_state <- function() {
  setNames(vapply(burial_names, draw_from_prior, numeric(1L)), burial_names)
}

# Optional strict Stage-order mode needs a valid starting state. This function
# finds two cut points with high joint marginal probability and draws each Stage
# from its appropriately truncated OxCal marginal.
draw_stage_ordered_state <- function() {
  all_grid <- sort(unique(unlist(base_grid, use.names = FALSE)))
  cuts <- (all_grid[-length(all_grid)] + all_grid[-1L]) / 2
  if (length(cuts) < 2L) stop("Insufficient calendar grid for Stage initialization.")

  log_m1 <- vapply(cuts, function(cut) {
    sum(vapply(stage_I, function(p) {
      z <- sum(base_mass[[p]][base_grid[[p]] < cut])
      if (z > 0) log(z) else -Inf
    }, numeric(1L)))
  }, numeric(1L))

  log_m3 <- vapply(cuts, function(cut) {
    sum(vapply(stage_III, function(p) {
      z <- sum(base_mass[[p]][base_grid[[p]] > cut])
      if (z > 0) log(z) else -Inf
    }, numeric(1L)))
  }, numeric(1L))

  best <- -Inf
  best_i <- best_j <- NA_integer_
  for (i in seq_len(length(cuts) - 1L)) {
    if (!is.finite(log_m1[[i]])) next
    for (j in (i + 1L):length(cuts)) {
      if (!is.finite(log_m3[[j]])) next
      lm2 <- sum(vapply(stage_II, function(p) {
        z <- sum(base_mass[[p]][base_grid[[p]] > cuts[[i]] &
                                  base_grid[[p]] < cuts[[j]]])
        if (z > 0) log(z) else -Inf
      }, numeric(1L)))
      score <- log_m1[[i]] + lm2 + log_m3[[j]]
      if (is.finite(score) && score > best) {
        best <- score
        best_i <- i
        best_j <- j
      }
    }
  }

  if (!is.finite(best)) stop("Could not construct a Stage-ordered initial state.")
  c12 <- cuts[[best_i]]
  c23 <- cuts[[best_j]]
  state <- setNames(numeric(length(burial_names)), burial_names)

  for (p in burial_names) {
    g <- base_grid[[p]]
    w <- base_mass[[p]]
    s <- stage_number[[p]]
    keep <- if (s == 1L) g < c12 else if (s == 2L) g > c12 & g < c23 else g > c23
    if (!any(keep) || sum(w[keep]) <= 0) stop("No valid initial support for ", p, ".")
    state[[p]] <- sample_value(g[keep], w[keep])
  }
  state
}

initialize_state <- function() {
  for (attempt in seq_len(100000L)) {
    state <- if (enforce_global_stage_order) draw_stage_ordered_state() else draw_independent_state()
    if (global_stage_valid(state) && is.finite(total_loglik(state))) return(state)
  }
  stop(
    "Could not find an initial state satisfying all directed constraints. ",
    "Check relationship directions, Mean/SD values, and Stage assignments."
  )
}

# ------------------------------- 5. MCMC sampler -----------------------------

run_chain <- function(chain_id, chain_seed) {
  set.seed(chain_seed)
  state <- initialize_state()
  samples <- matrix(
    NA_real_, nrow = n_keep_per_chain, ncol = length(burial_names),
    dimnames = list(NULL, burial_names)
  )

  site_proposed <- setNames(integer(length(burial_names)), burial_names)
  site_accepted <- setNames(integer(length(burial_names)), burial_names)
  redraw_proposed <- redraw_accepted <- integer(length(components))
  shift_proposed <- shift_accepted <- integer(length(components))
  kept <- 0L
  progress_points <- unique(pmax(1L, round(c(0.25, 0.5, 0.75, 1) * n_iter)))
  shift_values <- seq(-max_shift_years, max_shift_years, by = 5L)
  shift_values <- shift_values[shift_values != 0]

  for (iter in seq_len(n_iter)) {
    # Independence Metropolis-within-Gibbs. Proposal = OxCal prior, so the
    # burial's prior term and proposal density cancel in the acceptance ratio.
    for (p in sample(burial_names, length(burial_names), replace = FALSE)) {
      site_proposed[[p]] <- site_proposed[[p]] + 1L
      old_value <- state[[p]]
      proposed <- draw_from_prior(p)
      if (identical(proposed, old_value)) {
        site_accepted[[p]] <- site_accepted[[p]] + 1L
        next
      }
      if (!local_stage_valid(p, proposed, state)) next

      old_score <- local_loglik(p, state)
      state[[p]] <- proposed
      new_score <- local_loglik(p, state)
      log_alpha <- new_score - old_score

      accept <- is.finite(new_score) &&
        (log_alpha >= 0 || log(runif(1L)) < log_alpha)
      if (accept) {
        site_accepted[[p]] <- site_accepted[[p]] + 1L
      } else {
        state[[p]] <- old_value
      }
    }

    # Occasional component redraw from the product of component priors.
    # All component prior/proposal terms cancel; only relationship likelihood
    # and optional Stage support determine acceptance.
    for (cidx in seq_along(components)) {
      comp <- components[[cidx]]
      if (runif(1L) < component_redraw_prob) {
        redraw_proposed[[cidx]] <- redraw_proposed[[cidx]] + 1L
        old_values <- state[comp]
        old_score <- total_loglik(state)
        state[comp] <- vapply(comp, draw_from_prior, numeric(1L))
        new_score <- if (global_stage_valid(state)) total_loglik(state) else -Inf
        log_alpha <- new_score - old_score
        accept <- is.finite(new_score) &&
          (log_alpha >= 0 || log(runif(1L)) < log_alpha)
        if (accept) {
          redraw_accepted[[cidx]] <- redraw_accepted[[cidx]] + 1L
        } else {
          state[comp] <- old_values
        }
      }

      # Common-shift block move. This preserves all internal time gaps and can
      # move a connected group together across calibration-curve modes.
      if (runif(1L) < component_shift_prob) {
        shift_proposed[[cidx]] <- shift_proposed[[cidx]] + 1L
        shift <- sample(shift_values, size = 1L)
        old_values <- state[comp]
        new_values <- old_values + shift
        valid_grid <- vapply(seq_along(comp), function(j) {
          new_values[[j]] %in% base_grid[[comp[[j]]]]
        }, logical(1L))

        if (all(valid_grid)) {
          old_score <- total_loglik(state)
          log_prior_ratio <- 0
          for (j in seq_along(comp)) {
            p <- comp[[j]]
            old_idx <- match(old_values[[j]], base_grid[[p]])
            new_idx <- match(new_values[[j]], base_grid[[p]])
            log_prior_ratio <- log_prior_ratio +
              base_log_mass[[p]][[new_idx]] - base_log_mass[[p]][[old_idx]]
          }
          state[comp] <- new_values
          new_score <- if (global_stage_valid(state)) total_loglik(state) else -Inf
          log_alpha <- log_prior_ratio + new_score - old_score
          accept <- is.finite(new_score) && is.finite(log_alpha) &&
            (log_alpha >= 0 || log(runif(1L)) < log_alpha)
          if (accept) {
            shift_accepted[[cidx]] <- shift_accepted[[cidx]] + 1L
          } else {
            state[comp] <- old_values
          }
        }
      }
    }

    if (iter > burnin && ((iter - burnin) %% thin == 0L)) {
      kept <- kept + 1L
      samples[kept, ] <- state
    }

    if (iter %in% progress_points) {
      cat(sprintf("Chain %d progress: %d%%\n", chain_id, round(100 * iter / n_iter)))
    }
  }

  if (kept != n_keep_per_chain) samples <- samples[seq_len(kept), , drop = FALSE]

  list(
    samples = samples,
    site_proposed = site_proposed,
    site_accepted = site_accepted,
    redraw_proposed = redraw_proposed,
    redraw_accepted = redraw_accepted,
    shift_proposed = shift_proposed,
    shift_accepted = shift_accepted
  )
}

chain_results <- vector("list", n_chains)
for (ch in seq_len(n_chains)) {
  cat("Starting chain ", ch, " of ", n_chains, "...\n", sep = "")
  chain_results[[ch]] <- run_chain(ch, master_seed + 10007L * ch)
}

chain_samples <- lapply(chain_results, function(z) z$samples)
all_samples <- do.call(rbind, chain_samples)
chain_id_vector <- rep(seq_len(n_chains), vapply(chain_samples, nrow, integer(1L)))
cat("\nRetained ", nrow(all_samples), " joint posterior samples.\n\n", sep = "")

# -------------------------- 6. MCMC convergence checks -----------------------

basic_rhat <- function(chain_vectors) {
  m <- length(chain_vectors)
  if (m < 2L) return(NA_real_)
  n <- min(vapply(chain_vectors, length, integer(1L)))
  if (n < 2L) return(NA_real_)
  z <- vapply(chain_vectors, function(x) x[seq_len(n)], numeric(n))
  chain_means <- colMeans(z)
  chain_vars <- apply(z, 2L, var)
  W <- mean(chain_vars)
  B <- n * var(chain_means)
  if (!is.finite(W) || W <= 0) return(ifelse(B == 0, 1, NA_real_))
  var_hat <- ((n - 1) / n) * W + B / n
  sqrt(var_hat / W)
}

single_chain_ess <- function(x) {
  n <- length(x)
  if (n < 4L || var(x) == 0) return(as.numeric(n))
  lag_max <- min(1000L, n - 1L)
  ac <- as.numeric(acf(x, lag.max = lag_max, plot = FALSE, demean = TRUE)$acf)[-1L]
  if (length(ac) == 0L) return(as.numeric(n))
  pair_sums <- ac[seq(1L, length(ac), by = 2L)] +
    c(ac[seq(2L, length(ac), by = 2L)], 0)[seq_along(ac[seq(1L, length(ac), by = 2L)])]
  first_negative <- which(pair_sums < 0)[1L]
  if (!is.na(first_negative)) pair_sums <- pair_sums[seq_len(first_negative - 1L)]
  tau <- 1 + 2 * sum(pair_sums)
  if (!is.finite(tau) || tau < 1) tau <- 1
  min(n, n / tau)
}

mcmc_diag <- do.call(rbind, lapply(burial_names, function(p) {
  vectors <- lapply(chain_samples, function(z) z[, p])
  site_prop <- sum(vapply(chain_results, function(z) z$site_proposed[[p]], integer(1L)))
  site_acc <- sum(vapply(chain_results, function(z) z$site_accepted[[p]], integer(1L)))
  data.frame(
    burial = p,
    stage = stage_label[[p]],
    related_degree = length(edge_index_by_burial[[p]]),
    rhat = basic_rhat(vectors),
    approximate_ess = sum(vapply(vectors, single_chain_ess, numeric(1L))),
    single_site_acceptance_rate = if (site_prop > 0) site_acc / site_prop else NA_real_,
    stringsAsFactors = FALSE
  )
}))

component_diag <- do.call(rbind, lapply(seq_along(components), function(cidx) {
  rp <- sum(vapply(chain_results, function(z) z$redraw_proposed[[cidx]], integer(1L)))
  ra <- sum(vapply(chain_results, function(z) z$redraw_accepted[[cidx]], integer(1L)))
  sp <- sum(vapply(chain_results, function(z) z$shift_proposed[[cidx]], integer(1L)))
  sa <- sum(vapply(chain_results, function(z) z$shift_accepted[[cidx]], integer(1L)))
  data.frame(
    component = cidx,
    burials = paste(components[[cidx]], collapse = ";"),
    redraw_proposals = rp,
    redraw_acceptance_rate = if (rp > 0) ra / rp else NA_real_,
    shift_proposals = sp,
    shift_acceptance_rate = if (sp > 0) sa / sp else NA_real_,
    stringsAsFactors = FALSE
  )
}))

write.csv(mcmc_diag, file.path(out_dir, "mcmc_burial_diagnostics.csv"), row.names = FALSE)
write.csv(component_diag, file.path(out_dir, "mcmc_component_diagnostics.csv"), row.names = FALSE)

# --------------------- 7. Rao-Blackwell posterior densities ------------------
# For burial i, average p(t_i | all other sampled dates, kinship constraints)
# over retained MCMC states. This is smoother and lower-noise than a histogram
# of sampled dates while retaining the original OxCal grid and multimodality.

set.seed(master_seed + 909091L)
if (nrow(all_samples) > rb_max_samples) {
  rb_rows <- sort(sample(seq_len(nrow(all_samples)), rb_max_samples))
} else {
  rb_rows <- seq_len(nrow(all_samples))
}
rb_samples <- all_samples[rb_rows, , drop = FALSE]
cat("Rao-Blackwell density reconstruction uses ", nrow(rb_samples), " states.\n", sep = "")

conditional_loglik_grid <- function(p, grid, sampled_state) {
  out <- rep(0, length(grid))
  idx <- edge_index_by_burial[[p]]
  if (length(idx) > 0L) {
    for (k in idx) {
      mu <- rel$mean[[k]]
      sig <- rel$sd[[k]]
      a <- rel$burial_a[[k]]
      b <- rel$burial_b[[k]]

      if (rel$effective_direction[[k]] == "unknown") {
        other <- if (p == a) sampled_state[[b]] else sampled_state[[a]]
        gap <- abs(grid - other)
        out <- out + dnorm(gap, mean = mu, sd = sig, log = TRUE)
      } else {
        earlier <- rel$earlier_burial[[k]]
        later <- rel$later_burial[[k]]
        if (p == earlier) {
          gap <- sampled_state[[later]] - grid
        } else {
          gap <- grid - sampled_state[[earlier]]
        }
        ll <- dnorm(gap, mean = mu, sd = sig, log = TRUE)
        ll[gap < 0] <- -Inf
        out <- out + ll
      }
    }
  }

  if (enforce_global_stage_order) {
    s <- stage_number[[p]]
    valid <- if (s == 1L) {
      grid < min(sampled_state[stage_II])
    } else if (s == 2L) {
      grid > max(sampled_state[stage_I]) & grid < min(sampled_state[stage_III])
    } else {
      grid > max(sampled_state[stage_II])
    }
    out[!valid] <- -Inf
  }
  out
}

posterior_mass <- posterior_density <- empirical_mass <-
  vector("list", length(burial_names))
names(posterior_mass) <- names(posterior_density) <- names(empirical_mass) <- burial_names

for (pi in seq_along(burial_names)) {
  p <- burial_names[[pi]]
  grid <- base_grid[[p]]

  # Empirical sample mass is retained as a diagnostic.
  empirical_counts <- tabulate(match(all_samples[, p], grid), nbins = length(grid))
  empirical_mass[[p]] <- empirical_counts / sum(empirical_counts)

  if (length(edge_index_by_burial[[p]]) == 0L && !enforce_global_stage_order) {
    # An unrelated burial factorizes from the kinship likelihood, so its exact
    # second-stage marginal is its original OxCal marginal.
    pmass <- base_mass[[p]]
  } else {
    acc <- numeric(length(grid))
    for (r in seq_len(nrow(rb_samples))) {
      sampled_state <- rb_samples[r, ]
      logw <- base_log_mass[[p]] + conditional_loglik_grid(p, grid, sampled_state)
      lse <- log_sum_exp(logw)
      if (!is.finite(lse)) {
        stop("No finite conditional support for burial ", p,
             " during Rao-Blackwell reconstruction.")
      }
      acc <- acc + exp(logw - lse)
    }
    pmass <- acc / nrow(rb_samples)
    pmass <- pmass / sum(pmass)
  }

  posterior_mass[[p]] <- pmass
  posterior_density[[p]] <- pmass / base_width[[p]]
  cat(sprintf("Density reconstruction: %d/%d %s\n", pi, length(burial_names), p))
}

# ------------------------------- 8. Main outputs -----------------------------

density_rows <- vector("list", length(burial_names))
summary_rows <- vector("list", length(burial_names))
hpd_rows <- list()
hpd_count <- 0L

for (i in seq_along(burial_names)) {
  p <- burial_names[[i]]
  grid <- base_grid[[p]]
  widths <- base_width[[p]]
  lower <- grid - widths / 2
  upper <- grid + widths / 2
  prior_mass <- base_mass[[p]]
  post_mass <- posterior_mass[[p]]
  prior_density <- base_density[[p]]
  post_density <- posterior_density[[p]]

  density_rows[[i]] <- data.frame(
    burial = p,
    stage = stage_label[[p]],
    calendar_year = grid,
    calendar_label = format_cal_year(grid),
    bin_lower = lower,
    bin_upper = upper,
    bin_width = widths,
    oxcal_density = prior_density,
    oxcal_probability_mass = prior_mass,
    kinship_posterior_density = post_density,
    kinship_posterior_probability_mass = post_mass,
    empirical_mcmc_probability_mass = empirical_mass[[p]],
    stringsAsFactors = FALSE
  )

  ps <- weighted_summary(grid, prior_mass)
  ks <- weighted_summary(grid, post_mass)
  diag_row <- mcmc_diag[mcmc_diag$burial == p, , drop = FALSE]
  summary_rows[[i]] <- data.frame(
    burial = p,
    stage = stage_label[[p]],
    related_degree = length(edge_index_by_burial[[p]]),
    oxcal_mean = ps[["mean"]],
    oxcal_mode = ps[["mode"]],
    oxcal_q025 = ps[["q025"]],
    oxcal_q16 = ps[["q16"]],
    oxcal_median = ps[["median"]],
    oxcal_q84 = ps[["q84"]],
    oxcal_q975 = ps[["q975"]],
    posterior_mean = ks[["mean"]],
    posterior_mode = ks[["mode"]],
    posterior_q025 = ks[["q025"]],
    posterior_q16 = ks[["q16"]],
    posterior_median = ks[["median"]],
    posterior_q84 = ks[["q84"]],
    posterior_q975 = ks[["q975"]],
    posterior_mean_shift_years = ks[["mean"]] - ps[["mean"]],
    rhat = diag_row$rhat,
    approximate_ess = diag_row$approximate_ess,
    single_site_acceptance_rate = diag_row$single_site_acceptance_rate,
    stringsAsFactors = FALSE
  )

  for (lev in c(0.682, 0.954)) {
    hpd_count <- hpd_count + 1L
    hpd_rows[[hpd_count]] <- hpd_intervals(
      grid, post_density, post_mass, lev, p, stage_label[[p]]
    )
  }
}

density_out <- do.call(rbind, density_rows)
summary_out <- do.call(rbind, summary_rows)
hpd_out <- do.call(rbind, hpd_rows)

# Human-readable calendar labels in the summary, while keeping numeric columns.
for (nm in c("oxcal_mode", "oxcal_q025", "oxcal_median", "oxcal_q975",
             "posterior_mode", "posterior_q025", "posterior_median", "posterior_q975")) {
  summary_out[[paste0(nm, "_label")]] <- format_cal_year(summary_out[[nm]])
}

write.csv(
  density_out,
  file.path(out_dir, "burial_posterior_density.csv"),
  row.names = FALSE
)
write.csv(
  summary_out,
  file.path(out_dir, "burial_posterior_summary.csv"),
  row.names = FALSE
)
write.csv(
  hpd_out,
  file.path(out_dir, "burial_posterior_hpd_intervals.csv"),
  row.names = FALSE
)

# Input audit: burials and relationships actually used.
burial_audit <- data.frame(
  burial = burial_names,
  oxcal_operation = unname(base_op[burial_names]),
  stage = unname(stage_label[burial_names]),
  relationship_degree = vapply(burial_names, function(p) length(edge_index_by_burial[[p]]), integer(1L)),
  grid_points = vapply(base_grid, length, integer(1L)),
  median_grid_step = vapply(base_grid, function(x) if (length(x) > 1L) median(diff(x)) else NA_real_, numeric(1L)),
  oxcal_density_integral = vapply(burial_names, function(p) {
    sum(base_density[[p]] * base_width[[p]])
  }, numeric(1L)),
  stringsAsFactors = FALSE
)
write.csv(burial_audit, file.path(out_dir, "burials_used.csv"), row.names = FALSE)

rel_out <- rel[, c(
  "id", "burial_a", "burial_b", "relationship", "mean", "sd",
  "direction_established", "effective_direction", "earlier_burial",
  "later_burial", "direction_source"
), drop = FALSE]
write.csv(rel_out, file.path(out_dir, "relationships_used.csv"), row.names = FALSE)

# ----------------------- 9. Relationship posterior checks --------------------

relationship_diag <- do.call(rbind, lapply(seq_len(nrow(rel)), function(k) {
  a <- rel$burial_a[[k]]
  b <- rel$burial_b[[k]]
  signed_ab <- all_samples[, b] - all_samples[, a]
  if (rel$effective_direction[[k]] == "unknown") {
    gap <- abs(signed_ab)
    direction_probability <- NA_real_
  } else {
    earlier <- rel$earlier_burial[[k]]
    later <- rel$later_burial[[k]]
    gap <- all_samples[, later] - all_samples[, earlier]
    direction_probability <- mean(gap >= 0)
  }
  qs <- as.numeric(quantile(gap, probs = c(0.025, 0.16, 0.5, 0.84, 0.975),
                            names = FALSE, type = 1))
  data.frame(
    id = rel$id[[k]],
    burial_a = a,
    burial_b = b,
    relationship = rel$relationship[[k]],
    input_mean = rel$mean[[k]],
    input_sd = rel$sd[[k]],
    effective_direction = rel$effective_direction[[k]],
    earlier_burial = rel$earlier_burial[[k]],
    later_burial = rel$later_burial[[k]],
    direction_source = rel$direction_source[[k]],
    posterior_gap_mean = mean(gap),
    posterior_gap_sd = sd(gap),
    posterior_gap_q025 = qs[[1L]],
    posterior_gap_q16 = qs[[2L]],
    posterior_gap_median = qs[[3L]],
    posterior_gap_q84 = qs[[4L]],
    posterior_gap_q975 = qs[[5L]],
    posterior_probability_direction_satisfied = direction_probability,
    posterior_mean_standardized_difference =
      (mean(gap) - rel$mean[[k]]) / rel$sd[[k]],
    stringsAsFactors = FALSE
  )
}))
write.csv(
  relationship_diag,
  file.path(out_dir, "relationship_posterior_diagnostics.csv"),
  row.names = FALSE
)

# ---------------------------- 10. Joint samples ------------------------------

if (save_joint_samples) {
  sample_out <- data.frame(
    chain = chain_id_vector,
    retained_iteration = unlist(lapply(chain_samples, function(z) seq_len(nrow(z)))),
    all_samples,
    check.names = FALSE
  )
  write.csv(
    sample_out,
    file.path(out_dir, "joint_posterior_samples.csv"),
    row.names = FALSE
  )
}

# ----------------------------- 11. Comparison PDF ----------------------------

if (make_plots) {
  pdf_file <- file.path(out_dir, "burial_posterior_density_comparison.pdf")
  grDevices::pdf(pdf_file, width = 11.7, height = 8.3, onefile = TRUE)
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(3, 3), mar = c(3.3, 3.5, 2.2, 0.8), oma = c(0, 0, 1.5, 0))
  for (p in burial_names) {
    z <- density_out[density_out$burial == p, , drop = FALSE]
    ymax <- max(c(z$oxcal_density, z$kinship_posterior_density), na.rm = TRUE)
    plot(
      z$calendar_year, z$kinship_posterior_density,
      type = "n", xlab = "Calendar year (cal BCE)",
      ylab = "Probability density", main = paste0(p, " — ", stage_label[[p]]),
      ylim = c(0, ymax * 1.08), xaxt = "n"
    )
    polygon(
      c(z$calendar_year, rev(z$calendar_year)),
      c(z$kinship_posterior_density, rep(0, nrow(z))),
      col = "grey80", border = NA
    )
    lines(z$calendar_year, z$kinship_posterior_density, lwd = 1.5)
    lines(z$calendar_year, z$oxcal_density, lty = 2, lwd = 1.1, col = "grey35")
    ticks <- pretty(range(z$calendar_year), n = 5)
    axis(1, at = ticks, labels = abs(round(ticks)))
    box()
    if (p == burial_names[[1L]]) {
      legend(
        "topright", legend = c("Kinship-updated", "OxCal archaeological"),
        lty = c(1, 2), lwd = c(1.5, 1.1), col = c("black", "grey35"),
        bty = "n", cex = 0.75
      )
    }
  }
  mtext("Northern Cemetery: archaeological vs kinship-updated posterior densities",
        outer = TRUE, cex = 1.1)
  par(old_par)
  grDevices::dev.off()
}

# ----------------------------- 12. Run metadata ------------------------------

metadata <- data.frame(
  setting = c(
    "model", "model_version", "oxcal_file", "relationship_file", "master_seed", "n_chains",
    "n_iter", "burnin", "thin", "retained_samples", "rb_samples",
    "relationship_count", "established_direction_count", "established_directions",
    "nm27_nm21_treatment",
    "infer_stage_direction", "enforce_global_stage_order",
    "component_redraw_prob", "component_shift_prob", "max_shift_years",
    "statistical_note"
  ),
  value = c(
    "OxCal marginal posteriors x kinship gap likelihoods",
    "Three-scenario framework 2026-09-20; relationship degree/direction read from scenario input",
    normalizePath(oxcal_file, winslash = "/", mustWork = FALSE),
    normalizePath(relationship_file, winslash = "/", mustWork = FALSE),
    master_seed, n_chains, n_iter, burnin, thin, nrow(all_samples), nrow(rb_samples),
    nrow(rel), sum(rel$direction_yes),
    paste(sort(observed_directed_pairs), collapse = "; "),
    "Undirected Grandparent-Grandchild; absolute gap N(35,32)",
    infer_stage_direction, enforce_global_stage_order,
    component_redraw_prob, component_shift_prob, max_shift_years,
    paste(
      "Modular approximation using the product of burial-level OxCal marginal",
      "posteriors; original OxCal joint correlations are unavailable in the CSV."
    )
  ),
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(out_dir, "run_metadata.csv"), row.names = FALSE)

cat("\nCompleted successfully. Main density output:\n")
cat(file.path(out_dir, "burial_posterior_density.csv"), "\n")
cat("For every burial, sum(posterior probability mass) = 1 and\n")
cat("sum(posterior density x bin width) = 1.\n")
