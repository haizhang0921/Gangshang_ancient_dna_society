#!/usr/bin/env Rscript

# ============================================================================
# Figure 4B — Updated summed posterior densities (SPD), revised layout
# Gangshang Northern Cemetery + Southern Cemetery + moat
# Updated chronology: final MAIN kinship models, 2026-09-20
#
# IMPORTANT
# ---------
# Main-text chronology only:
#   North: north_R_primary_main/burial_posterior_density.csv
#   South: south_R_exact_main/burial_posterior_density.csv
#   Moat : gangshang_moat.csv -> Combine("Moat date") posterior
#
# Alternative kinship scenarios are NOT used in Figure 4B; they remain
# supplementary sensitivity analyses.
#
# SPD definition used here:
#   - Each top-level burial-event posterior integrates to 1.
#   - Northern Stage SPD = raw sum of burial-event posterior densities
#     within that Stage (22 / 5 / 7 events).
#   - Southern Cemetery SPD = raw sum of the 14 burial-event posterior
#     densities.
#   - Moat curve = ONE Combine("Moat date") posterior, integrating to 1.
#   - No normalization by the number of events and no kernel smoothing.
#
# Therefore SPD height reflects BOTH the number of represented burial events
# and their temporal concentration. It is not a direct estimate of population
# size, demographic growth, mortality, or cemetery-use intensity.
#
# x-axis is cal BCE, in accordance with the revised manuscript convention.
# ============================================================================

# ----------------------------------------------------------------------------
# 0. User settings
# ----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

get_script_dir <- function() {
  z <- commandArgs(trailingOnly = FALSE)
  hit <- z[grepl("^--file=", z)]
  if (length(hit) > 0L) {
    return(normalizePath(dirname(sub("^--file=", "", hit[1L])),
                         winslash = "/", mustWork = FALSE))
  }
  getwd()
}

script_dir <- get_script_dir()
project_dir <- if (length(args) >= 1L) {
  normalizePath(args[1L], winslash = "/", mustWork = FALSE)
} else {
  script_dir
}

output_dir <- file.path(project_dir, "Figure4B_output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_stem <- if (length(args) >= 2L) args[2L] else "Figure4B_updated_main_models_v3"

# Plot window: reproduces approximately the chronological window of the old
# Figure 4B (formerly ~5200–4300 cal BP), now expressed in cal BCE.
x_limits_BCE <- c(3250, 2350)  # older -> younger
x_major_break <- 100
x_minor_break <- 50

# Main figure size, suitable for later A+B+C assembly.
figure_width_in  <- 7.25
figure_height_in <- 3.55
tiff_dpi <- 600
base_family <- "Arial"

# Optional interpretive labels. Set FALSE for a purely chronological SPD panel.
show_interpretive_labels <- TRUE
show_social_transition_arrows <- TRUE

# These labels are deliberately editable here rather than hard-coded below.
# Current wording is aligned with the revised manuscript's more cautious
# interpretation of kinship organization.
label_yellow_river <- "Yangshao-related gene inflow"
label_high_ranked  <- "High-ranked tombs"
label_steppe       <- "Steppe-related ancestry"
label_power        <- "Social power"

left_social_label  <- "Maternal anchoring"
right_social_label <- "Increasing paternal-line prominence"

# To reproduce the older wording exactly, replace the two strings above with:
# left_social_label  <- "Matrilineal Clan"
# right_social_label <- "Patrilineal Clan"


# ----------------------------------------------------------------------------
# 1. Packages
# ----------------------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Install with install.packages('ggplot2').",
       call. = FALSE)
}
library(ggplot2)


# ----------------------------------------------------------------------------
# 2. Input helpers: read either extracted result folders or ZIP archives
# ----------------------------------------------------------------------------

find_extracted_file <- function(base_dir, suffix_pattern) {
  ff <- list.files(base_dir, recursive = TRUE, full.names = TRUE)
  ff_norm <- gsub("\\\\", "/", ff)
  hit <- ff[grepl(suffix_pattern, ff_norm)]
  if (length(hit) == 0L) return(NA_character_)
  hit[1L]
}

read_csv_from_zip_pattern <- function(zip_file, member_pattern) {
  entries <- unzip(zip_file, list = TRUE)$Name
  hit <- entries[grepl(member_pattern, gsub("\\\\", "/", entries))]
  if (length(hit) != 1L) {
    stop(
      "Expected exactly one ZIP member matching:\n  ", member_pattern,
      "\nin:\n  ", zip_file,
      "\nFound: ", length(hit),
      if (length(hit) > 0L) paste0("\n", paste(hit, collapse = "\n")) else "",
      call. = FALSE
    )
  }
  read.csv(unz(zip_file, hit), stringsAsFactors = FALSE,
           check.names = FALSE, fileEncoding = "UTF-8")
}

read_north_main <- function(base_dir) {
  # Prefer extracted formal-R MAIN result.
  f <- find_extracted_file(
    base_dir,
    "north_R_primary_main/burial_posterior_density\\.csv$"
  )
  if (!is.na(f)) {
    message("North MAIN input: ", f)
    return(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE,
                    fileEncoding = "UTF-8"))
  }

  zip_candidates <- c(
    file.path(base_dir, "GS_north_cemetery.zip"),
    file.path(base_dir, "Gangshang_North_3_kinship_scenarios_20260920_FIXED3.zip")
  )
  zip_candidates <- zip_candidates[file.exists(zip_candidates)]
  if (length(zip_candidates) == 0L) {
    stop("Cannot find Northern Cemetery MAIN result or GS_north_cemetery.zip.",
         call. = FALSE)
  }
  message("North MAIN input ZIP: ", zip_candidates[1L])
  read_csv_from_zip_pattern(
    zip_candidates[1L],
    "north_R_primary_main/burial_posterior_density\\.csv$"
  )
}

read_south_main <- function(base_dir) {
  f <- find_extracted_file(
    base_dir,
    "south_R_exact_main/burial_posterior_density\\.csv$"
  )
  if (!is.na(f)) {
    message("South MAIN input: ", f)
    return(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE,
                    fileEncoding = "UTF-8"))
  }

  zip_candidates <- c(
    file.path(base_dir, "GS_south_cemetery.zip"),
    file.path(base_dir, "Gangshang_South_2_kinship_scenarios_20260920.zip")
  )
  zip_candidates <- zip_candidates[file.exists(zip_candidates)]
  if (length(zip_candidates) == 0L) {
    stop("Cannot find Southern Cemetery MAIN result or GS_south_cemetery.zip.",
         call. = FALSE)
  }
  message("South MAIN input ZIP: ", zip_candidates[1L])
  read_csv_from_zip_pattern(
    zip_candidates[1L],
    "south_R_exact_main/burial_posterior_density\\.csv$"
  )
}

find_moat_file <- function(base_dir) {
  candidates <- c(
    file.path(base_dir, "gangshang_moat.csv"),
    file.path(script_dir, "gangshang_moat.csv"),
    file.path(getwd(), "gangshang_moat.csv")
  )
  candidates <- unique(candidates[file.exists(candidates)])
  if (length(candidates) == 0L) {
    # final fallback: recursive search
    ff <- list.files(base_dir, pattern = "^gangshang_moat\\.csv$",
                     recursive = TRUE, full.names = TRUE)
    if (length(ff) == 0L) stop("Cannot find gangshang_moat.csv.", call. = FALSE)
    return(ff[1L])
  }
  candidates[1L]
}

north <- read_north_main(project_dir)
south <- read_south_main(project_dir)
moat_file <- find_moat_file(project_dir)
message("Moat input: ", moat_file)
moat <- read.csv(moat_file, stringsAsFactors = FALSE, check.names = FALSE,
                 fileEncoding = "UTF-8")


# ----------------------------------------------------------------------------
# 3. Validate authoritative datasets
# ----------------------------------------------------------------------------

required_north <- c("burial", "stage", "calendar_year", "bin_width",
                    "kinship_posterior_density")
required_south <- c("burial", "calendar_year", "bin_width",
                    "kinship_posterior_density")
required_moat <- c("op", "name", "type", "value", "probability")

for (x in setdiff(required_north, names(north))) {
  stop("Northern result is missing column: ", x, call. = FALSE)
}
for (x in setdiff(required_south, names(south))) {
  stop("Southern result is missing column: ", x, call. = FALSE)
}
for (x in setdiff(required_moat, names(moat))) {
  stop("Moat export is missing column: ", x, call. = FALSE)
}

north$burial <- trimws(as.character(north$burial))
north$stage <- trimws(as.character(north$stage))
north$calendar_year <- as.numeric(north$calendar_year)
north$bin_width <- as.numeric(north$bin_width)
north$kinship_posterior_density <- as.numeric(north$kinship_posterior_density)

south$burial <- trimws(as.character(south$burial))
south$calendar_year <- as.numeric(south$calendar_year)
south$bin_width <- as.numeric(south$bin_width)
south$kinship_posterior_density <- as.numeric(south$kinship_posterior_density)

north_events <- unique(north[, c("burial", "stage")])
stage_counts <- table(factor(north_events$stage,
                             levels = c("Stage I", "Stage II", "Stage III")))
if (nrow(north_events) != 34L ||
    any(as.integer(stage_counts) != c(22L, 5L, 7L))) {
  stop("Northern MAIN input does not contain 34 events in Stage I/II/III = 22/5/7.",
       call. = FALSE)
}
if (length(unique(south$burial)) != 14L) {
  stop("Southern MAIN input does not contain 14 top-level burial events.",
       call. = FALSE)
}

# Check that every burial-event posterior integrates to approximately 1.
check_event_normalization <- function(dat, id_col, density_col) {
  ids <- unique(dat[[id_col]])
  vals <- vapply(ids, function(id) {
    d <- dat[dat[[id_col]] == id, , drop = FALSE]
    sum(d[[density_col]] * d$bin_width, na.rm = TRUE)
  }, numeric(1L))
  max(abs(vals - 1))
}

north_norm_error <- check_event_normalization(
  north, "burial", "kinship_posterior_density"
)
south_norm_error <- check_event_normalization(
  south, "burial", "kinship_posterior_density"
)

if (north_norm_error > 1e-5) warning("North normalization deviation > 1e-5: ", north_norm_error)
if (south_norm_error > 1e-5) warning("South normalization deviation > 1e-5: ", south_norm_error)


# ----------------------------------------------------------------------------
# 4. Build raw SPDs
# ----------------------------------------------------------------------------

sum_density_by_year <- function(dat, group_name) {
  ans <- aggregate(
    dat$kinship_posterior_density,
    by = list(calendar_year = dat$calendar_year),
    FUN = sum,
    na.rm = TRUE
  )
  names(ans)[2L] <- "summed_density"
  ans$group <- group_name
  ans$cal_BCE <- -ans$calendar_year
  ans
}

spd_list <- list(
  sum_density_by_year(north[north$stage == "Stage I", , drop = FALSE],
                      "Northern Cemetery Stage I"),
  sum_density_by_year(north[north$stage == "Stage II", , drop = FALSE],
                      "Northern Cemetery Stage II"),
  sum_density_by_year(north[north$stage == "Stage III", , drop = FALSE],
                      "Northern Cemetery Stage III"),
  sum_density_by_year(south, "Southern Cemetery")
)

# Moat: use the latest formal OxCal Combine("Moat date") posterior.
moat_date <- moat[
  moat$op == "Combine" &
    moat$name == "Moat date" &
    moat$type == "posterior",
  , drop = FALSE
]
if (nrow(moat_date) == 0L) {
  stop("Could not find Combine('Moat date') posterior in gangshang_moat.csv.",
       call. = FALSE)
}
moat_date$value <- as.numeric(moat_date$value)
moat_date$probability <- as.numeric(moat_date$probability)

moat_spd <- data.frame(
  calendar_year = moat_date$value,
  summed_density = moat_date$probability,
  group = "Wall and Moat",
  cal_BCE = -moat_date$value,
  stringsAsFactors = FALSE
)

spd <- do.call(rbind, c(spd_list, list(moat_spd)))

# Display order matches the previous Figure 4B legend.
group_levels <- c(
  "Northern Cemetery Stage I",
  "Northern Cemetery Stage II",
  "Northern Cemetery Stage III",
  "Southern Cemetery",
  "Wall and Moat"
)
spd$group <- factor(spd$group, levels = group_levels)


# ----------------------------------------------------------------------------
# 5. Diagnostic summaries: area and peak of each curve
# ----------------------------------------------------------------------------

integrate_group <- function(d) {
  d <- d[order(d$calendar_year), , drop = FALSE]
  if (nrow(d) < 2L) return(NA_real_)
  bw <- median(diff(sort(unique(d$calendar_year))))
  sum(d$summed_density * bw, na.rm = TRUE)
}

peak_group <- function(d) {
  d[which.max(d$summed_density), c("cal_BCE", "summed_density")]
}

group_diagnostics <- do.call(rbind, lapply(group_levels, function(g) {
  d <- spd[spd$group == g, , drop = FALSE]
  pk <- peak_group(d)
  data.frame(
    group = g,
    integrated_area = integrate_group(d),
    peak_cal_BCE = pk$cal_BCE,
    peak_summed_density = pk$summed_density,
    stringsAsFactors = FALSE
  )
}))

write.csv(
  spd,
  file.path(output_dir, paste0(output_stem, "_SPD_plot_data.csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  group_diagnostics,
  file.path(output_dir, paste0(output_stem, "_SPD_diagnostics.csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nSPD integrated areas (expected approximately 22, 5, 7, 14, 1):\n")
print(group_diagnostics)


# ----------------------------------------------------------------------------
# 6. Colours — retain Figure 4B visual grammar
# ----------------------------------------------------------------------------

figure4b_cols <- c(
  "Northern Cemetery Stage I"   = "#8BBE22",  # green
  "Northern Cemetery Stage II"  = "#22B9B2",  # turquoise
  "Northern Cemetery Stage III" = "#B78AE2",  # purple
  "Southern Cemetery"           = "#EF7E78",  # salmon
  "Wall and Moat"               = "#163A70"   # deep blue
)

# Slight transparency allows overlap to remain legible.
area_alpha <- 0.72
line_width <- 0.55


# ----------------------------------------------------------------------------
# 7. Annotation helpers
# ----------------------------------------------------------------------------

# Calculate y values from actual new curves so labels remain stable after reruns.
curve_peak <- setNames(group_diagnostics$peak_summed_density,
                       group_diagnostics$group)
curve_peak_x <- setNames(group_diagnostics$peak_cal_BCE,
                         group_diagnostics$group)

y_data_max <- max(spd$summed_density, na.rm = TRUE)
y_top <- y_data_max * 1.32

# Function returning interpolated SPD height at a requested BCE date.
get_y_at <- function(group_name, x_bce) {
  d <- spd[spd$group == group_name, c("cal_BCE", "summed_density"), drop = FALSE]
  d <- d[order(d$cal_BCE), , drop = FALSE]
  approx(d$cal_BCE, d$summed_density, xout = x_bce,
         rule = 2, ties = mean)$y
}

# Revised annotation positions. These are interpretive labels only; they do
# not enter the chronology model. Coordinates can be edited here if needed.
# "Urban society" has been removed. "Steppe-related ancestry" and
# "Social power" are vertically stacked above the Southern Cemetery curve.
southern_label_x <- unname(curve_peak_x[["Southern Cemetery"]])
southern_label_base_y <- unname(curve_peak[["Southern Cemetery"]])

annotation_table <- data.frame(
  label = c(
    label_yellow_river,
    label_high_ranked,
    label_steppe,
    label_power
  ),
  x = c(
    3075,                                 # within early Northern use
    curve_peak_x[["Northern Cemetery Stage III"]],
    southern_label_x,
    southern_label_x
  ),
  y = c(
    get_y_at("Northern Cemetery Stage I", 3075) + 0.045 * y_top,
    curve_peak[["Northern Cemetery Stage III"]] + 0.055 * y_top,
    southern_label_base_y + 0.125 * y_top,  # upper label
    southern_label_base_y + 0.060 * y_top   # lower label
  ),
  stringsAsFactors = FALSE
)


# ----------------------------------------------------------------------------
# 8. Build Figure 4B
# ----------------------------------------------------------------------------

p <- ggplot(spd, aes(x = cal_BCE, y = summed_density,
                     fill = group, colour = group, group = group)) +
  geom_area(alpha = area_alpha, position = "identity", linewidth = 0) +
  geom_line(linewidth = line_width) +
  scale_fill_manual(values = figure4b_cols, breaks = group_levels, name = NULL) +
  scale_colour_manual(values = figure4b_cols, breaks = group_levels, guide = "none") +
  scale_x_reverse(
    limits = x_limits_BCE,
    breaks = seq(2400, 3200, by = x_major_break),
    minor_breaks = seq(2350, 3250, by = x_minor_break),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = pretty(c(0, y_data_max), n = 4),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Modeled date (cal BCE)",
    y = "Summed Probability Density"
  ) +
  coord_cartesian(ylim = c(0, y_top), clip = "off") +
  theme_classic(base_size = 9, base_family = base_family) +
  theme(
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25),
    panel.grid.minor.x = element_line(colour = "grey95", linewidth = 0.18),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.30, colour = "black"),
    axis.text = element_text(colour = "black", size = 8),
    axis.title.x = element_text(size = 9, margin = margin(t = 6)),
    axis.title.y = element_text(size = 9, margin = margin(r = 6)),
    legend.position = c(0.81, 0.80),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "white", colour = "grey65", linewidth = 0.35),
    legend.key.height = grid::unit(3.4, "mm"),
    legend.key.width = grid::unit(4.2, "mm"),
    legend.text = element_text(size = 7.3),
    plot.margin = margin(t = 6, r = 7, b = 46, l = 5)
  )

# Panel label B.
p <- p + annotate(
  "text", x = x_limits_BCE[1] - 6, y = y_top * 1.015,
  label = "B", hjust = 0, vjust = 1,
  size = 7.2, family = base_family
)

# Interpretive labels, using updated chronology and wording.
if (isTRUE(show_interpretive_labels)) {
  p <- p + geom_text(
    data = annotation_table,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    family = base_family,
    size = 3.05,
    colour = "black",
    fontface = "bold"
  )
}

# Social-organization arrows are drawn OUTSIDE the plotting panel. coord_cartesian
# uses y >= 0, while clip = "off" allows these negative-y annotations to appear
# below the x-axis without changing the probability-density axis.
if (isTRUE(show_social_transition_arrows)) {
  # Move arrows and labels farther below the x-axis to avoid overlap with
  # the x-axis title. Two arrows only, with a small gap centred on 2900 BCE.
  arrow_y <- -0.205 * y_top
  text_y  <- -0.285 * y_top
  transition_bce <- 2900
  arrow_gap_half_width <- 10  # 20-year visual gap: 2910-2890 BCE
  left_arrow_end  <- transition_bce + arrow_gap_half_width
  right_arrow_start <- transition_bce - arrow_gap_half_width

  p <- p +
    annotate(
      "segment",
      x = 3230, xend = left_arrow_end,
      y = arrow_y, yend = arrow_y,
      linewidth = 0.38,
      colour = "grey25",
      arrow = grid::arrow(length = grid::unit(2.0, "mm"), type = "closed")
    ) +
    annotate(
      "text",
      x = (3230 + left_arrow_end) / 2, y = text_y,
      label = left_social_label,
      size = 3.15,
      family = base_family,
      fontface = "bold",
      colour = "grey15"
    ) +
    annotate(
      "segment",
      x = right_arrow_start, xend = 2410,
      y = arrow_y, yend = arrow_y,
      linewidth = 0.38,
      colour = "grey25",
      arrow = grid::arrow(length = grid::unit(2.0, "mm"), type = "closed")
    ) +
    annotate(
      "text",
      x = (right_arrow_start + 2410) / 2, y = text_y,
      label = right_social_label,
      size = 3.00,
      family = base_family,
      fontface = "bold",
      colour = "grey15"
    )
}


# ----------------------------------------------------------------------------
# 9. Save publication files
# ----------------------------------------------------------------------------

pdf_file <- file.path(output_dir, paste0(output_stem, ".pdf"))
png_file <- file.path(output_dir, paste0(output_stem, "_preview.png"))
tiff_file <- file.path(output_dir, paste0(output_stem, "_600dpi.tiff"))

pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"

ggsave(
  pdf_file, p,
  width = figure_width_in, height = figure_height_in,
  units = "in", device = pdf_device, limitsize = FALSE
)

ggsave(
  png_file, p,
  width = figure_width_in, height = figure_height_in,
  units = "in", dpi = 220, limitsize = FALSE
)

ggsave(
  tiff_file, p,
  width = figure_width_in, height = figure_height_in,
  units = "in", dpi = tiff_dpi,
  compression = "lzw", limitsize = FALSE
)

cat("\nFigure 4B completed.\n")
cat("North: MAIN kinship model only\n")
cat("South: MAIN second-degree SM1-SM2 model only\n")
cat("Moat : latest Combine('Moat date') posterior\n")
cat("Calendar scale: cal BCE\n")
cat("\nOutputs:\n")
cat("  ", pdf_file, "\n", sep = "")
cat("  ", png_file, "\n", sep = "")
cat("  ", tiff_file, "\n", sep = "")
cat("  ", file.path(output_dir, paste0(output_stem, "_SPD_plot_data.csv")), "\n", sep = "")
cat("  ", file.path(output_dir, paste0(output_stem, "_SPD_diagnostics.csv")), "\n", sep = "")
