#!/usr/bin/env Rscript

# ============================================================================
# Figure 4C — Updated burial statistics
# Input: gs_burials.csv
#
# The pedigree colours are EXACTLY the same as Figure 4A:
#   A = #FF0000
#   B = #1E2BFF
#   C = #1B9E54
#   D = #8E199B
#   E = #A05252
#   F = #20D7E6
#   G = #FF7F2A
#   H = #FFB000
#   I = #FF00FF
#   Unassigned = #A9A9A9
#
# Figure structure:
#   1. Area (m²)
#   2. Pottery
#   3. Jadeware
#   4. Bone/horn/ivory artefacts
#   5. Animal bone
#   6. Grave value
#
# Tomb order is taken directly from the first appearance in gs_burials.csv,
# so the archaeological order supplied in the input file is preserved.
# Northern and Southern cemeteries are separated by a light dashed line.
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Input / output
# ----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)

  if (length(file_arg) > 0L) {
    p <- sub("^--file=", "", file_arg[1L])
    return(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(e) ""
    )
    if (nzchar(p)) {
      return(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()

input_file <- if (length(args) >= 1L) {
  args[1L]
} else {
  candidates <- unique(c(
    file.path(script_dir, "gs_burials.csv"),
    file.path(getwd(), "gs_burials.csv"),
    "/mnt/data/gs_burials.csv"
  ))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    message("gs_burials.csv was not found automatically; please select it.")
    file.choose()
  } else {
    hit[1L]
  }
}

input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)

output_dir <- if (length(args) >= 2L) {
  args[2L]
} else {
  file.path(dirname(input_file), "Figure4C_output")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_stem <- "Figure4C_updated_burial_statistics"

pdf_file  <- file.path(output_dir, paste0(output_stem, ".pdf"))
tiff_file <- file.path(output_dir, paste0(output_stem, "_600dpi.tiff"))
png_file  <- file.path(output_dir, paste0(output_stem, "_preview.png"))
csv_file  <- file.path(output_dir, paste0(output_stem, "_plot_data.csv"))

# Stand-alone panel size for later assembly with Figure 4A/B.
fig_width_mm  <- 96
fig_height_mm <- 118
base_font     <- 7.2

# Set TRUE only if Figure 4C is to be shown independently.
# In the combined Figure 4, Figure 4A already contains the pedigree legend.
show_pedigree_legend <- FALSE


# ----------------------------------------------------------------------------
# 1. Packages
# ----------------------------------------------------------------------------

required_pkgs <- c("ggplot2", "scales")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Please install the following R package(s) first: ",
    paste(missing_pkgs, collapse = ", "),
    "\nExample: install.packages(c(\"ggplot2\", \"scales\"))",
    call. = FALSE
  )
}

library(ggplot2)
library(scales)


# ----------------------------------------------------------------------------
# 2. Read and validate input
# ----------------------------------------------------------------------------

burials <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

required_columns <- c("tomb", "pedigree", "metric", "value")
missing_columns <- setdiff(required_columns, names(burials))

if (length(missing_columns) > 0L) {
  stop(
    "gs_burials.csv is missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

burials$tomb     <- trimws(as.character(burials$tomb))
burials$pedigree <- trimws(as.character(burials$pedigree))
burials$metric   <- trimws(as.character(burials$metric))
burials$value    <- suppressWarnings(as.numeric(burials$value))

bad_rows <- (
  !nzchar(burials$tomb) |
  !nzchar(burials$pedigree) |
  !nzchar(burials$metric) |
  !is.finite(burials$value)
)

if (any(bad_rows)) {
  stop(
    "The input contains ", sum(bad_rows),
    " invalid row(s) in tomb/pedigree/metric/value.",
    call. = FALSE
  )
}

# Expected metric set and order for Figure 4C.
metric_levels <- c(
  "Area (m²)",
  "Pottery",
  "Jadeware",
  "Bone/horn/ivory artefacts",
  "Animal bone",
  "Grave value"
)

unexpected_metrics <- setdiff(unique(burials$metric), metric_levels)
missing_metrics <- setdiff(metric_levels, unique(burials$metric))

if (length(unexpected_metrics) > 0L) {
  stop(
    "Unexpected metric(s): ",
    paste(unexpected_metrics, collapse = ", "),
    call. = FALSE
  )
}

if (length(missing_metrics) > 0L) {
  stop(
    "Missing metric(s): ",
    paste(missing_metrics, collapse = ", "),
    call. = FALSE
  )
}

# Preserve the archaeological tomb order exactly as supplied by the CSV.
tomb_order <- unique(burials$tomb)

# Check that each tomb has one pedigree assignment only.
pedigree_per_tomb <- split(burials$pedigree, burials$tomb)
non_unique_pedigree <- names(pedigree_per_tomb)[
  vapply(pedigree_per_tomb, function(x) length(unique(x)) != 1L, logical(1))
]

if (length(non_unique_pedigree) > 0L) {
  stop(
    "These tombs have more than one pedigree assignment: ",
    paste(non_unique_pedigree, collapse = ", "),
    call. = FALSE
  )
}

# Check that every tomb has all six metrics exactly once.
cell_counts <- table(burials$tomb, burials$metric)
if (any(cell_counts != 1L)) {
  stop(
    "Each tomb must occur exactly once for each of the six metrics.",
    call. = FALSE
  )
}

# Figure 4A pedigree categories.
allowed_pedigrees <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "Unassigned")
unknown_pedigrees <- setdiff(unique(burials$pedigree), allowed_pedigrees)

if (length(unknown_pedigrees) > 0L) {
  stop(
    "Unknown pedigree label(s): ",
    paste(unknown_pedigrees, collapse = ", "),
    call. = FALSE
  )
}

burials$tomb <- factor(burials$tomb, levels = tomb_order)
burials$metric <- factor(burials$metric, levels = metric_levels, ordered = TRUE)
burials$pedigree <- factor(
  burials$pedigree,
  levels = allowed_pedigrees
)

# Identify the Northern/Southern split from tomb prefixes.
north_tombs <- tomb_order[grepl("^N[MW]", tomb_order)]
south_tombs <- tomb_order[grepl("^SM", tomb_order)]

if ((length(north_tombs) + length(south_tombs)) != length(tomb_order)) {
  warning(
    "Some tomb labels are not recognised as Northern (NM/NW) or Southern (SM). " ,
    "The cemetery separator may need manual checking."
  )
}

# This assumes all northern tombs precede all southern tombs, as in the supplied CSV.
separator_x <- length(north_tombs) + 0.5


# ----------------------------------------------------------------------------
# 3. Colours — EXACTLY synchronized with Figure 4A
# ----------------------------------------------------------------------------

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
  "Unassigned" = "#A9A9A9"
)

# Export the exact plot data used for the figure.
write.csv(
  burials,
  csv_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ----------------------------------------------------------------------------
# 4. Plot
# ----------------------------------------------------------------------------

p <- ggplot(
  burials,
  aes(x = tomb, y = value, fill = pedigree)
) +
  geom_col(
    width = 0.72,
    colour = NA
  ) +
  facet_wrap(
    ~ metric,
    ncol = 1,
    scales = "free_y",
    strip.position = "top"
  ) +
  geom_vline(
    xintercept = separator_x,
    linewidth = 0.28,
    linetype = "dashed",
    colour = "#A8A8A8"
  ) +
  scale_fill_manual(
    values = pedigree_cols,
    breaks = allowed_pedigrees,
    drop = FALSE,
    name = NULL
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.07)),
    breaks = breaks_pretty(n = 4)
  ) +
  labs(
    x = NULL,
    y = NULL,
    tag = "C"
  ) +
  theme_classic(
    base_family = "Helvetica",
    base_size = base_font
  ) +
  theme(
    # Compact vertical structure for integration with panels A and B.
    panel.spacing.y = grid::unit(0.55, "mm"),

    # Restrained horizontal guides, matching the visual weight of A/B.
    panel.grid.major.y = element_line(
      colour = "#E5E5E5",
      linewidth = 0.26
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),

    # Fine axes.
    axis.line = element_line(
      colour = "#777777",
      linewidth = 0.28
    ),
    axis.ticks.y = element_line(
      colour = "#777777",
      linewidth = 0.28
    ),
    axis.ticks.x = element_blank(),
    axis.ticks.length = grid::unit(1.0, "mm"),

    axis.text.y = element_text(
      colour = "#333333",
      size = 5.7,
      margin = margin(r = 1.4)
    ),
    axis.text.x = element_text(
      colour = "#333333",
      size = 4.7,
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      margin = margin(t = 1.0)
    ),

    strip.background = element_blank(),
    strip.text = element_text(
      colour = "#222222",
      size = 6.4,
      face = "bold",
      margin = margin(t = 1.0, b = 0.5)
    ),

    plot.tag = element_text(
      family = "Helvetica",
      size = 14,
      face = "plain",
      colour = "#111111"
    ),
    plot.tag.position = c(0.005, 0.995),

    plot.margin = margin(
      t = 3.5,
      r = 1.5,
      b = 1.5,
      l = 2.0,
      unit = "mm"
    )
  )

if (isTRUE(show_pedigree_legend)) {
  p <- p +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(colour = NA)
      )
    ) +
    theme(
      legend.position = "top",
      legend.justification = "left",
      legend.text = element_text(size = 5.8),
      legend.key.size = grid::unit(3.0, "mm"),
      legend.spacing.x = grid::unit(1.0, "mm")
    )
} else {
  p <- p + guides(fill = "none")
}


# ----------------------------------------------------------------------------
# 5. Export
# ----------------------------------------------------------------------------

pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"

ggsave(
  filename = pdf_file,
  plot = p,
  width = fig_width_mm,
  height = fig_height_mm,
  units = "mm",
  device = pdf_device,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = tiff_file,
  plot = p,
  width = fig_width_mm,
  height = fig_height_mm,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = png_file,
  plot = p,
  width = fig_width_mm,
  height = fig_height_mm,
  units = "mm",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

cat("\nFigure 4C completed.\n")
cat("Input : ", input_file, "\n", sep = "")
cat("Tombs : ", length(tomb_order), " (Northern = ", length(north_tombs),
    "; Southern = ", length(south_tombs), ")\n", sep = "")
cat("Metrics: ", paste(metric_levels, collapse = "; "), "\n", sep = "")
cat("\nPedigree colours are identical to Figure 4A.\n")
cat("Unassigned burials use the same Figure 4A default grey: #A9A9A9.\n")
cat("\nOutputs:\n")
cat("  PDF  : ", pdf_file, "\n", sep = "")
cat("  TIFF : ", tiff_file, "\n", sep = "")
cat("  PNG  : ", png_file, "\n", sep = "")
cat("  Data : ", csv_file, "\n", sep = "")
