#!/usr/bin/env Rscript

# =============================================================================
# Figure 4A — updated site-wide chronology (MAIN models only)
# Gangshang Northern Cemetery + moat + Southern Cemetery
#
# Updated 2026-09-20 after completion of:
#   1) Northern Cemetery: 48 human 14C determinations, revised MAIN kinship model
#      -> north_R_primary_main/burial_posterior_density.csv
#   2) Southern Cemetery: archaeology-only OxCal model + MAIN SM1-SM2
#      second-degree relationship
#      -> south_R_exact_main/burial_posterior_density.csv
#   3) Moat: latest OxCal posterior export
#      -> gangshang_moat.csv, Combine("Moat date") posterior
#
# IMPORTANT
#   - Figure 4A uses ONLY the MAIN kinship reconstructions used in the main text.
#   - Northern alternative kinship scenarios and the Southern sibling alternative
#     belong in the Supplementary sensitivity analyses, not in Figure 4A.
#   - All modelled calendar ages are plotted as cal BCE.
#   - Raw conventional radiocarbon ages are not plotted here.
#
# Figure style follows the previous Figure 4A:
#   * one row per top-level burial event
#   * mean + 68.3%, 95.4%, 99.7% HPD intervals
#   * pedigree-specific colours A-I
#   * grey for burials not assigned a pedigree colour in the original figure
#   * black for the moat
#   * North displayed as archaeological Phases I-IV
#
# -----------------------------------------------------------------------------
# EASY USE IN RSTUDIO
# -----------------------------------------------------------------------------
# Put this script in the project folder that contains either the result folders
# or GS_north_cemetery.zip / GS_south_cemetery.zip and gangshang_moat.csv.
# Then run:
#
#   source("06_plot_Figure4A_updated_main_models.R")
#
# Optional command line:
#   Rscript 06_plot_Figure4A_updated_main_models.R PROJECT_ROOT OUTPUT_STEM OUTPUT_DIR
#
# Example:
#   Rscript 06_plot_Figure4A_updated_main_models.R . Figure4A_updated Figure4A_output
# =============================================================================


# =============================================================================
# 1. Arguments and paths
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

resolve_script_dir <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  hit <- full_args[grepl("^--file=", full_args)]
  if (length(hit) == 0L) return(getwd())
  normalizePath(dirname(sub("^--file=", "", hit[1L])), winslash = "/", mustWork = FALSE)
}

script_dir <- resolve_script_dir()

project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
} else {
  script_dir
}

output_stem <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  "Figure4A_updated_main_models"
}

output_dir <- if (length(args) >= 3L) {
  args[[3L]]
} else {
  file.path(project_root, "Figure4A_output")
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(output_dir)) {
  stop("Cannot create output directory: ", output_dir, call. = FALSE)
}


# =============================================================================
# 2. Input discovery helpers
# =============================================================================

norm_slash <- function(x) gsub("\\\\", "/", x)

find_result_source <- function(root, exact_suffix, preferred_zip_name = NULL) {

  # 2.1 Search extracted folders first.
  hits <- list.files(
    root,
    pattern = "burial_posterior_density\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(hits) > 0L) {
    hits_norm <- norm_slash(normalizePath(hits, winslash = "/", mustWork = FALSE))
    keep <- endsWith(hits_norm, exact_suffix)
    hits <- hits[keep]
  }

  if (length(hits) == 1L) {
    return(list(kind = "file", path = hits[[1L]], inner = NA_character_))
  }

  if (length(hits) > 1L) {
    stop(
      "More than one extracted result matched:\n  ",
      paste(hits, collapse = "\n  "),
      "\nPlease remove duplicate old result folders or edit the input path in the script.",
      call. = FALSE
    )
  }

  # 2.2 If no extracted folder exists, search ZIP archives.
  zip_hits <- list.files(root, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)

  if (!is.null(preferred_zip_name)) {
    preferred <- zip_hits[basename(zip_hits) == preferred_zip_name]
    if (length(preferred) > 0L) {
      zip_hits <- c(preferred, setdiff(zip_hits, preferred))
    }
  }

  found <- list()

  for (zz in zip_hits) {
    zlist <- tryCatch(unzip(zz, list = TRUE), error = function(e) NULL)
    if (is.null(zlist)) next

    inner_hits <- zlist$Name[endsWith(norm_slash(zlist$Name), exact_suffix)]

    if (length(inner_hits) == 1L) {
      found[[length(found) + 1L]] <- list(kind = "zip", path = zz, inner = inner_hits[[1L]])
    }
  }

  if (length(found) == 1L) return(found[[1L]])

  if (length(found) > 1L) {
    desc <- vapply(found, function(z) paste0(z$path, " :: ", z$inner), character(1))
    stop(
      "More than one ZIP result matched:\n  ",
      paste(desc, collapse = "\n  "),
      "\nPlease keep only the current result archive in PROJECT_ROOT.",
      call. = FALSE
    )
  }

  stop(
    "Could not find result ending in: ", exact_suffix,
    "\nPROJECT_ROOT = ", root,
    call. = FALSE
  )
}

read_result_source <- function(src) {
  if (identical(src$kind, "file")) {
    return(read.csv(
      src$path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8"
    ))
  }

  read.csv(
    unz(src$path, src$inner),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

find_unique_named_file <- function(root, filename) {
  direct <- file.path(root, filename)
  if (file.exists(direct)) return(direct)

  hits <- list.files(root, pattern = paste0("^", gsub("\\.", "\\\\.", filename), "$"),
                     recursive = TRUE, full.names = TRUE)

  if (length(hits) == 1L) return(hits[[1L]])
  if (length(hits) > 1L) {
    stop(
      "More than one file named ", filename, " was found:\n  ",
      paste(hits, collapse = "\n  "),
      call. = FALSE
    )
  }

  stop("Cannot find ", filename, " under ", root, call. = FALSE)
}

north_source <- find_result_source(
  project_root,
  "/north_R_primary_main/burial_posterior_density.csv",
  preferred_zip_name = "GS_north_cemetery.zip"
)

south_source <- find_result_source(
  project_root,
  "/south_R_exact_main/burial_posterior_density.csv",
  preferred_zip_name = "GS_south_cemetery.zip"
)

moat_file <- find_unique_named_file(project_root, "gangshang_moat.csv")


# =============================================================================
# 3. Read inputs
# =============================================================================

north <- read_result_source(north_source)
south <- read_result_source(south_source)
moat <- read.csv(
  moat_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

required_north <- c("burial", "calendar_year", "bin_width", "kinship_posterior_density")
required_south <- c("burial", "calendar_year", "bin_width", "kinship_posterior_density")
required_moat <- c("op", "name", "type", "value", "probability")

if (length(setdiff(required_north, names(north))) > 0L) {
  stop("Northern MAIN density file is missing required columns.", call. = FALSE)
}
if (length(setdiff(required_south, names(south))) > 0L) {
  stop("Southern MAIN density file is missing required columns.", call. = FALSE)
}
if (length(setdiff(required_moat, names(moat))) > 0L) {
  stop("gangshang_moat.csv is missing required columns.", call. = FALSE)
}

north$burial <- trimws(as.character(north$burial))
north$calendar_year <- as.numeric(north$calendar_year)
north$bin_width <- as.numeric(north$bin_width)
north$kinship_posterior_density <- as.numeric(north$kinship_posterior_density)

south$burial <- trimws(as.character(south$burial))
south$calendar_year <- as.numeric(south$calendar_year)
south$bin_width <- as.numeric(south$bin_width)
south$kinship_posterior_density <- as.numeric(south$kinship_posterior_density)

moat$value <- as.numeric(moat$value)
moat$probability <- as.numeric(moat$probability)

moat_subset <- moat[
  moat$op == "Combine" &
    moat$name == "Moat date" &
    moat$type == "posterior",
  , drop = FALSE
]

if (nrow(moat_subset) == 0L) {
  stop("Could not find Combine('Moat date') posterior in gangshang_moat.csv.", call. = FALSE)
}


# =============================================================================
# 4. Probability / HPD helper functions
# =============================================================================

normalize_mass <- function(density, bin_width) {
  mass <- density * bin_width
  s <- sum(mass, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) stop("Invalid probability mass.", call. = FALSE)
  mass / s
}

weighted_mean_safe <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w >= 0
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0L || sum(w) <= 0) return(NA_real_)
  sum(x * w) / sum(w)
}

# Highest-posterior-density intervals on a regular grid; multimodal intervals retained.
get_hpd_intervals <- function(x, density, bin_width, level) {

  x <- as.numeric(x)
  d <- as.numeric(density)
  bw <- as.numeric(bin_width[1L])

  keep <- is.finite(x) & is.finite(d) & d >= 0
  x <- x[keep]
  d <- d[keep]

  if (length(x) == 0L) {
    return(data.frame(lower = numeric(0), upper = numeric(0)))
  }

  mass <- normalize_mass(d, bw)
  ord <- order(d, decreasing = TRUE)
  k <- which(cumsum(mass[ord]) >= level)[1L]
  chosen <- sort(x[ord[seq_len(k)]])

  if (length(chosen) == 0L) {
    return(data.frame(lower = numeric(0), upper = numeric(0)))
  }

  groups <- list()
  start_idx <- 1L

  if (length(chosen) == 1L) {
    groups[[1L]] <- chosen
  } else {
    for (ii in 2:length(chosen)) {
      if ((chosen[ii] - chosen[ii - 1L]) > bw * 1.5) {
        groups[[length(groups) + 1L]] <- chosen[start_idx:(ii - 1L)]
        start_idx <- ii
      }
    }
    groups[[length(groups) + 1L]] <- chosen[start_idx:length(chosen)]
  }

  out <- do.call(rbind, lapply(groups, function(g) {
    data.frame(
      lower = min(g) - bw / 2,
      upper = max(g) + bw / 2,
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

summarise_distribution <- function(cal_BCE, density, bin_width) {
  x <- as.numeric(cal_BCE)
  d <- as.numeric(density)
  bw <- as.numeric(bin_width[1L])

  keep <- is.finite(x) & is.finite(d) & d >= 0
  x <- x[keep]
  d <- d[keep]

  mass <- normalize_mass(d, bw)

  list(
    mean = weighted_mean_safe(x, mass),
    hpd_683 = get_hpd_intervals(x, d, bw, 0.683),
    hpd_954 = get_hpd_intervals(x, d, bw, 0.954),
    hpd_997 = get_hpd_intervals(x, d, bw, 0.997)
  )
}

intervals_to_string <- function(x) {
  if (is.null(x) || nrow(x) == 0L) return(NA_character_)
  paste(
    apply(x, 1L, function(z) paste0(round(z[["lower"]]), "–", round(z[["upper"]]))),
    collapse = "; "
  )
}


# =============================================================================
# 5. Figure order and pedigree colour assignments
# =============================================================================
# These preserve the pedigree colour assignments used in the previous Figure 4A.
# They represent the aDNA pedigree display, not merely the subset of relationships
# used as temporal likelihoods in the chronology model.

event_spec <- data.frame(
  event_id = c(
    # North archaeological Phase I
    "NM24", "NM28", "NM25", "NM30", "NM27", "NM21", "NM22", "NM12",
    "NM23", "NM18", "NM26", "NM20",

    # North archaeological Phase II
    "NM29", "NM6", "NM5", "NM11", "NM8", "NM9", "NM10", "NM2", "NM3", "NM4",

    # North archaeological Phase III (urn burials)
    "NW13", "NW10", "NW7", "NW6", "NW14",

    # North archaeological Phase IV
    "NM34", "NM31", "NM14", "NM32", "NM15", "NM16", "NM13",

    # Moat
    "Moat date",

    # South
    "SM13", "SM10", "SM9", "SM17", "SM15", "SM3", "SM14",
    "SM12", "SM6", "SM2", "SM1", "SM8", "SM11", "SM5"
  ),
  source = c(
    rep("north", 12), rep("north", 10), rep("north", 5), rep("north", 7),
    "moat", rep("south", 14)
  ),
  block = c(
    rep("North Cemetery\nPhase I", 12),
    rep("North Cemetery\nPhase II", 10),
    rep("North Cemetery\nPhase III", 5),
    rep("North Cemetery\nPhase IV", 7),
    "Moat",
    rep("South Cemetery", 14)
  ),
  display_label = c(
    "NM24", "NM28", "NM25", "NM30", "NM27", "NM21", "NM22", "NM12",
    "NM23", "NM18", "NM26", "NM20",

    "NM29", "NM6-1", "NM5-2", "NM11", "NM8", "NM9", "NM10", "NM2", "NM3", "NM4-1",

    "NW13", "NW10", "NW7", "NW6", "NW14",

    "NM34", "NM31", "NM14", "NM32", "NM15", "NM16", "NM13",

    "Moat",

    "SM13", "SM10", "SM9", "SM17", "SM15", "SM3", "SM14",
    "SM12", "SM6", "SM2", "SM1", "SM8", "SM11", "SM5"
  ),
  pedigree = c(
    # North Phase I
    "A", "A", "A", "B", "B", "B", "C", "C", "C", NA, NA, "E",

    # North Phase II
    "E", "F", "F", "G", "G", NA, "D", "D", "D", "D",

    # North Phase III
    NA, "D", NA, NA, NA,

    # North Phase IV
    "D", "D", "D", "D", "H", "H", NA,

    # Moat
    "Moat",

    # South
    NA, NA, NA, NA, NA, NA, NA, NA, NA, "I", "I", NA, NA, NA
  ),
  stringsAsFactors = FALSE
)

# Verify current authoritative event sets.
missing_north <- setdiff(event_spec$event_id[event_spec$source == "north"], unique(north$burial))
missing_south <- setdiff(event_spec$event_id[event_spec$source == "south"], unique(south$burial))

if (length(missing_north) > 0L) {
  stop("Missing Northern burial event(s): ", paste(missing_north, collapse = ", "), call. = FALSE)
}
if (length(missing_south) > 0L) {
  stop("Missing Southern burial event(s): ", paste(missing_south, collapse = ", "), call. = FALSE)
}

if (length(unique(north$burial)) != 34L) {
  warning("Northern input does not contain exactly 34 top-level burial events.")
}
if (length(unique(south$burial)) != 14L) {
  warning("Southern input does not contain exactly 14 top-level burial events.")
}


# =============================================================================
# 6. Summarise all displayed posterior distributions in cal BCE
# =============================================================================

summary_list <- vector("list", nrow(event_spec))
names(summary_list) <- event_spec$event_id

summary_table <- data.frame(
  order_id = seq_len(nrow(event_spec)),
  source = event_spec$source,
  block = event_spec$block,
  event_id = event_spec$event_id,
  display_label = event_spec$display_label,
  pedigree = event_spec$pedigree,
  mean_cal_BCE = NA_real_,
  hpd_68_3_cal_BCE = NA_character_,
  hpd_95_4_cal_BCE = NA_character_,
  hpd_99_7_cal_BCE = NA_character_,
  stringsAsFactors = FALSE
)

for (ii in seq_len(nrow(event_spec))) {

  id <- event_spec$event_id[ii]
  src <- event_spec$source[ii]

  if (src == "north") {

    d <- north[north$burial == id,
               c("calendar_year", "bin_width", "kinship_posterior_density"),
               drop = FALSE]

    tmp <- summarise_distribution(
      cal_BCE = -d$calendar_year,
      density = d$kinship_posterior_density,
      bin_width = d$bin_width
    )

  } else if (src == "south") {

    d <- south[south$burial == id,
               c("calendar_year", "bin_width", "kinship_posterior_density"),
               drop = FALSE]

    tmp <- summarise_distribution(
      cal_BCE = -d$calendar_year,
      density = d$kinship_posterior_density,
      bin_width = d$bin_width
    )

  } else if (src == "moat") {

    # In the OxCal export, probability is a density on a 5-year grid.
    d <- moat_subset[, c("value", "probability"), drop = FALSE]

    tmp <- summarise_distribution(
      cal_BCE = -d$value,
      density = d$probability,
      bin_width = rep(5, nrow(d))
    )

  } else {
    stop("Unknown source: ", src, call. = FALSE)
  }

  summary_list[[id]] <- tmp
  summary_table$mean_cal_BCE[ii] <- tmp$mean
  summary_table$hpd_68_3_cal_BCE[ii] <- intervals_to_string(tmp$hpd_683)
  summary_table$hpd_95_4_cal_BCE[ii] <- intervals_to_string(tmp$hpd_954)
  summary_table$hpd_99_7_cal_BCE[ii] <- intervals_to_string(tmp$hpd_997)
}


# =============================================================================
# 7. Drawing helpers
# =============================================================================

pedigree_cols <- c(
  "A" = "#FF0000",
  "B" = "#1E2BFF",
  "C" = "#1B9E54",
  "D" = "#8E199B",
  "E" = "#A05252",
  "F" = "#20D7E6",
  "G" = "#FF7F2A",
  "H" = "#FFB000",
  "I" = "#FF00FF",
  "Moat" = "#000000"
)

default_col <- "#A9A9A9"
panel_bg <- "#FFFFFF"

# Rounded interval polygon, matching previous Figure 4A visual grammar.
draw_rounded_interval <- function(x0, x1, y, col, fill_col, height,
                                  lwd = 1.2, corner_x = 7) {

  width <- abs(x1 - x0)
  if (!is.finite(width) || width <= 0) return(invisible(NULL))

  rx <- min(corner_x, width / 2)
  ry <- min(height * 0.45, height / 2)

  theta1 <- seq(pi, pi / 2, length.out = 20)
  theta2 <- seq(pi / 2, 0, length.out = 20)
  theta3 <- seq(0, -pi / 2, length.out = 20)
  theta4 <- seq(-pi / 2, -pi, length.out = 20)

  x_left <- min(x0, x1)
  x_right <- max(x0, x1)

  cx_tl <- x_left + rx;  cy_tl <- y + height / 2 - ry
  cx_tr <- x_right - rx; cy_tr <- y + height / 2 - ry
  cx_br <- x_right - rx; cy_br <- y - height / 2 + ry
  cx_bl <- x_left + rx;  cy_bl <- y - height / 2 + ry

  xs <- c(
    cx_tl + rx * cos(theta1),
    cx_tr + rx * cos(theta2),
    cx_br + rx * cos(theta3),
    cx_bl + rx * cos(theta4)
  )
  ys <- c(
    cy_tl + ry * sin(theta1),
    cy_tr + ry * sin(theta2),
    cy_br + ry * sin(theta3),
    cy_bl + ry * sin(theta4)
  )

  polygon(xs, ys, col = fill_col, border = col, lwd = lwd, xpd = NA)
}

draw_event <- function(s, y, col) {

  fill_997 <- adjustcolor(col, alpha.f = 0.12)
  fill_954 <- adjustcolor(col, alpha.f = 0.22)
  fill_683 <- adjustcolor(col, alpha.f = 0.38)

  if (nrow(s$hpd_997) > 0L) {
    apply(s$hpd_997, 1L, function(z) {
      draw_rounded_interval(z[["lower"]], z[["upper"]], y,
                            col, fill_997, 0.68, 1.0)
    })
  }

  if (nrow(s$hpd_954) > 0L) {
    apply(s$hpd_954, 1L, function(z) {
      draw_rounded_interval(z[["lower"]], z[["upper"]], y,
                            col, fill_954, 0.52, 1.0)
    })
  }

  if (nrow(s$hpd_683) > 0L) {
    apply(s$hpd_683, 1L, function(z) {
      draw_rounded_interval(z[["lower"]], z[["upper"]], y,
                            col, fill_683, 0.36, 1.1)
    })
  }

  points(
    x = s$mean,
    y = y,
    pch = 21,
    bg = "white",
    col = col,
    cex = 0.85,
    lwd = 1.1,
    xpd = NA
  )
}

safe_output_path <- function(directory, filename) {
  target <- file.path(directory, filename)
  if (!file.exists(target)) return(target)

  removed <- suppressWarnings(file.remove(target))
  if (isTRUE(removed) && !file.exists(target)) return(target)

  ext <- tools::file_ext(filename)
  stem <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", filename) else filename
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  new_name <- if (nzchar(ext)) paste0(stem, "_", stamp, ".", ext) else paste0(stem, "_", stamp)

  warning("Output file appears to be open/locked; saving as: ", new_name)
  file.path(directory, new_name)
}


# =============================================================================
# 8. Layout
# =============================================================================

n_events <- nrow(event_spec)
event_spec$y <- rev(seq_len(n_events))

block_levels <- unique(event_spec$block)
block_ends <- cumsum(rle(event_spec$block)$lengths)
separator_y <- event_spec$y[block_ends[-length(block_ends)]] - 0.5
block_midpoints <- sapply(block_levels, function(b) mean(event_spec$y[event_spec$block == b]))

# Preserve the broad temporal framing of the old Figure 4A, now in cal BCE.
xlim <- c(3500, 2200)
major_ticks <- seq(3400, 2200, by = -200)
minor_ticks <- seq(3500, 2200, by = -20)

# Left-side label columns.
event_label_x <- 3485
block_label_x <- 3375


# =============================================================================
# 9. Plot function
# =============================================================================

make_plot <- function() {

  op <- par(
    mar = c(4.5, 15.8, 2.0, 2.2),
    bg = "white",
    lend = "round",
    xaxs = "i",
    yaxs = "i",
    family = "serif"
  )
  on.exit(par(op), add = TRUE)

  plot(
    NA, NA,
    xlim = xlim,
    ylim = c(0.2, n_events + 1.0),
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    bty = "n"
  )

  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = panel_bg, border = NA)

  # Major chronological guides.
  for (xx in major_ticks) {
    abline(v = xx, col = "grey35", lwd = 0.8)
  }

  # Archaeological-area/phase separators.
  for (yy in separator_y) {
    abline(h = yy, col = "grey20", lwd = 1.0, lty = 2)
  }

  # Event intervals.
  for (ii in seq_len(nrow(event_spec))) {
    id <- event_spec$event_id[ii]
    ped <- event_spec$pedigree[ii]
    col_i <- if (!is.na(ped) && ped %in% names(pedigree_cols)) pedigree_cols[[ped]] else default_col
    draw_event(summary_list[[id]], event_spec$y[ii], col_i)
  }

  # Burial labels.
  text(
    x = event_label_x,
    y = event_spec$y,
    labels = event_spec$display_label,
    adj = c(0, 0.5),
    cex = 0.80,
    xpd = NA
  )

  # Group labels.
  for (b in block_levels) {
    text(
      x = block_label_x,
      y = block_midpoints[[b]],
      labels = b,
      adj = c(0, 0.5),
      cex = if (b == "South Cemetery") 1.10 else 1.05,
      xpd = NA
    )
  }

  # X axis: older BCE values on left.
  axis(
    side = 1,
    at = major_ticks,
    labels = major_ticks,
    cex.axis = 1.15,
    lwd = 1,
    lwd.ticks = 1,
    tcl = -0.7
  )

  axis(
    side = 1,
    at = minor_ticks,
    labels = FALSE,
    lwd = 0,
    lwd.ticks = 0.8,
    tcl = -0.35
  )

  mtext(
    "Modeled date (cal BCE)",
    side = 1,
    line = 3.0,
    cex = 1.30,
    font = 2
  )

  box(lwd = 1.2)

  # Panel label A.
  mtext("A", side = 3, line = 0.25, adj = -0.035, cex = 1.65, font = 2)

  # Probability-interval legend.
  legend(
    x = 2495,
    y = n_events + 0.45,
    legend = c(
      "Mean",
      "68.3% HPD",
      "95.4% HPD",
      "99.7% HPD"
    ),
    pch = c(21, 22, 22, 22),
    pt.bg = c(
      "white",
      adjustcolor("grey45", alpha.f = 0.38),
      adjustcolor("grey45", alpha.f = 0.22),
      adjustcolor("grey45", alpha.f = 0.12)
    ),
    pt.cex = c(0.95, 1.35, 1.35, 1.35),
    col = rep("grey45", 4),
    pt.lwd = c(1.1, 1.0, 1.0, 1.0),
    bty = "o",
    bg = "#FFFFFF",
    cex = 0.95,
    xjust = 0,
    yjust = 1
  )

  # Pedigree-colour legend.
  legend(
    x = 2480,
    y = n_events - 6.2,
    legend = c(
      "Pedigree A", "Pedigree B", "Pedigree C", "Pedigree D", "Pedigree E",
      "Pedigree F", "Pedigree G", "Pedigree H", "Pedigree I", "Moat"
    ),
    pch = rep(22, 10),
    pt.bg = c(
      adjustcolor(pedigree_cols["A"], alpha.f = 0.22),
      adjustcolor(pedigree_cols["B"], alpha.f = 0.22),
      adjustcolor(pedigree_cols["C"], alpha.f = 0.22),
      adjustcolor(pedigree_cols["D"], alpha.f = 0.22),
      adjustcolor(pedigree_cols["E"], alpha.f = 0.22),
      adjustcolor(pedigree_cols["F"], alpha.f = 0.22),
      adjustcolor(pedigree_cols["G"], alpha.f = 0.22),
      adjustcolor(pedigree_cols["H"], alpha.f = 0.22),
      adjustcolor(pedigree_cols["I"], alpha.f = 0.22),
      adjustcolor(pedigree_cols["Moat"], alpha.f = 0.22)
    ),
    pt.cex = rep(1.4, 10),
    pt.lwd = rep(1.0, 10),
    col = c(
      pedigree_cols["A"], pedigree_cols["B"], pedigree_cols["C"], pedigree_cols["D"],
      pedigree_cols["E"], pedigree_cols["F"], pedigree_cols["G"], pedigree_cols["H"],
      pedigree_cols["I"], pedigree_cols["Moat"]
    ),
    bty = "o",
    bg = "#FFFFFF",
    cex = 0.92,
    xjust = 0,
    yjust = 1
  )
}


# =============================================================================
# 10. Save files
# =============================================================================

graphics.off()

pdf_file <- safe_output_path(output_dir, paste0(output_stem, ".pdf"))
tiff_file <- safe_output_path(output_dir, paste0(output_stem, "_600dpi.tiff"))
png_file <- safe_output_path(output_dir, paste0(output_stem, "_preview.png"))
csv_file <- safe_output_path(output_dir, paste0(output_stem, "_summary.csv"))

write.csv(summary_table, csv_file, row.names = FALSE, fileEncoding = "UTF-8")

# PDF
pdf_open <- FALSE
tryCatch({
  grDevices::pdf(
    pdf_file,
    width = 10.5,
    height = 11.0,
    family = "serif",
    onefile = TRUE,
    useDingbats = FALSE
  )
  pdf_open <- TRUE
  make_plot()
  grDevices::dev.off()
  pdf_open <- FALSE
}, error = function(e) {
  if (pdf_open && grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
  stop("PDF output failed: ", conditionMessage(e), call. = FALSE)
})

# TIFF 600 dpi
tryCatch({
  grDevices::tiff(
    tiff_file,
    width = 10.5,
    height = 11.0,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  make_plot()
  grDevices::dev.off()
}, error = function(e) {
  if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
  warning("TIFF output failed: ", conditionMessage(e))
})

# PNG preview
tryCatch({
  grDevices::png(
    png_file,
    width = 10.5,
    height = 11.0,
    units = "in",
    res = 220
  )
  make_plot()
  grDevices::dev.off()
}, error = function(e) {
  if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
  warning("PNG output failed: ", conditionMessage(e))
})


# =============================================================================
# 11. Console summary
# =============================================================================

source_description <- function(src) {
  if (src$kind == "file") return(normalizePath(src$path, winslash = "/", mustWork = FALSE))
  paste0(normalizePath(src$path, winslash = "/", mustWork = FALSE), " :: ", src$inner)
}

cat("\n============================================================\n")
cat("Figure 4A completed successfully\n")
cat("============================================================\n")
cat("Northern MAIN posterior:\n  ", source_description(north_source), "\n", sep = "")
cat("Southern MAIN posterior:\n  ", source_description(south_source), "\n", sep = "")
cat("Moat posterior:\n  ", normalizePath(moat_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("\nDisplayed events:\n")
cat("  North = 34 top-level burial events\n")
cat("  Moat  = Combine('Moat date') posterior\n")
cat("  South = 14 top-level burial events\n")
cat("\nCalendar scale: cal BCE\n")
cat("Intervals: 68.3%, 95.4%, 99.7% HPD\n")
cat("\nOutput:\n")
cat("  PDF  : ", pdf_file, "\n", sep = "")
cat("  TIFF : ", tiff_file, "\n", sep = "")
cat("  PNG  : ", png_file, "\n", sep = "")
cat("  CSV  : ", csv_file, "\n", sep = "")
cat("============================================================\n\n")
