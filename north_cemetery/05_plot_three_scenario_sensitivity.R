#!/usr/bin/env Rscript

# ======================================================================
# Gangshang Northern Cemetery — kinship-scenario sensitivity figure
#
# Panel a: reconstructed boundaries and Stage durations under three
#          alternative kinship reconstructions.
# Panel b: posterior distributions for the four burial events that best
#          illustrate scenario sensitivity (NM30, NM27, NM22, NM23).
#
# Expected inputs (in working directory):
#   three_scenarios_boundary_summary_long.csv
#   three_scenarios_duration_summary_long.csv
#   north_R_boundaries_main/kinship_updated_hpd_intervals.csv
#   north_R_boundaries_alt_close/kinship_updated_hpd_intervals.csv
#   north_R_boundaries_alt_far/kinship_updated_hpd_intervals.csv
#   north_R_primary_main/burial_posterior_density.csv
#   north_R_primary_alt_close/burial_posterior_density.csv
#   north_R_primary_alt_far/burial_posterior_density.csv
#
# Output:
#   Gangshang_North_kinship_sensitivity_figure.pdf
#   Gangshang_North_kinship_sensitivity_figure_600dpi.tiff
#   Gangshang_North_kinship_sensitivity_figure_preview.png
# ======================================================================

pkgs <- c("ggplot2", "patchwork")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Package '", p, "' is required. Install with install.packages('", p, "').", call. = FALSE)
  }
}
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) args[[1L]] else getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

scenario_ids <- c("main", "alt_close", "alt_far")
scenario_labels <- c(
  main = "Main",
  alt_close = "Closer-kin alternative",
  alt_far = "More-distant alternative"
)
scenario_cols <- c(
  main = "#8C3B3B",
  alt_close = "#3B75A2",
  alt_far = "#C47A20"
)
scenario_shapes <- c(main = 16, alt_close = 15, alt_far = 17)
scenario_offsets <- c(main = 0.18, alt_close = 0, alt_far = -0.18)

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
}

# ----------------------------------------------------------------------
# 1. Summary means
# ----------------------------------------------------------------------
boundary_summary <- read_required(file.path(project_dir, "three_scenarios_boundary_summary_long.csv"))
duration_summary <- read_required(file.path(project_dir, "three_scenarios_duration_summary_long.csv"))
boundary_summary$mean_BCE <- boundary_summary$mu - 1950

# ----------------------------------------------------------------------
# 2. Exact HPD segments (retain multimodal/disjoint intervals)
# ----------------------------------------------------------------------
hpd_all <- do.call(rbind, lapply(scenario_ids, function(sc) {
  z <- read_required(file.path(project_dir, paste0("north_R_boundaries_", sc), "kinship_updated_hpd_intervals.csv"))
  z$scenario <- sc
  z
}))

# Boundary HPD segments converted from cal BP to cal BCE.
b_hpd <- subset(hpd_all, parameter_type == "Boundary" & name %in% c(
  "Cemetery Start", "Stage I to II", "Stage II to III", "Cemetery End"
))
b_hpd$xmin <- pmin(b_hpd$lower - 1950, b_hpd$upper - 1950)
b_hpd$xmax <- pmax(b_hpd$lower - 1950, b_hpd$upper - 1950)

boundary_levels <- c("Cemetery End", "Stage II to III", "Stage I to II", "Cemetery Start")
boundary_labels <- c(
  "Cemetery End",
  "Stage II → III",
  "Stage I → II",
  "Cemetery Start"
)
b_hpd$y_base <- match(b_hpd$name, boundary_levels)
b_hpd$y <- b_hpd$y_base + scenario_offsets[b_hpd$scenario]

b_mean <- subset(boundary_summary, name %in% boundary_levels)
b_mean$y_base <- match(b_mean$name, boundary_levels)
b_mean$y <- b_mean$y_base + scenario_offsets[b_mean$scenario]

p_boundaries <- ggplot() +
  geom_segment(
    data = subset(b_hpd, abs(level - 0.954) < 1e-6),
    aes(x = xmin, xend = xmax, y = y, yend = y, colour = scenario),
    linewidth = 0.45, alpha = 0.60
  ) +
  geom_segment(
    data = subset(b_hpd, abs(level - 0.683) < 1e-6),
    aes(x = xmin, xend = xmax, y = y, yend = y, colour = scenario),
    linewidth = 1.25, alpha = 0.95
  ) +
  geom_point(
    data = b_mean,
    aes(x = mean_BCE, y = y, colour = scenario, shape = scenario),
    size = 2.2, stroke = 0.25
  ) +
  scale_colour_manual(values = scenario_cols, labels = scenario_labels) +
  scale_shape_manual(values = scenario_shapes, labels = scenario_labels) +
  scale_y_continuous(breaks = seq_along(boundary_levels), labels = boundary_labels) +
  scale_x_reverse() +
  labs(title = "Boundaries", x = "Calendar age (cal BCE)", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 7.5),
    axis.title.x = element_text(size = 8.5),
    plot.title = element_text(size = 9.5, face = "bold"),
    legend.position = "none"
  )

# Duration HPD segments.
d_hpd <- subset(hpd_all, parameter_type == "Interval" & name %in% c(
  "Stage I Duration", "Stage II Duration", "Stage III Duration"
))
duration_levels <- c("Stage III Duration", "Stage II Duration", "Stage I Duration")
duration_labels <- c("Stage III", "Stage II", "Stage I")
d_hpd$y_base <- match(d_hpd$name, duration_levels)
d_hpd$y <- d_hpd$y_base + scenario_offsets[d_hpd$scenario]

d_mean <- subset(duration_summary, name %in% duration_levels)
d_mean$y_base <- match(d_mean$name, duration_levels)
d_mean$y <- d_mean$y_base + scenario_offsets[d_mean$scenario]

p_durations <- ggplot() +
  geom_segment(
    data = subset(d_hpd, abs(level - 0.954) < 1e-6),
    aes(x = lower, xend = upper, y = y, yend = y, colour = scenario),
    linewidth = 0.45, alpha = 0.60
  ) +
  geom_segment(
    data = subset(d_hpd, abs(level - 0.683) < 1e-6),
    aes(x = lower, xend = upper, y = y, yend = y, colour = scenario),
    linewidth = 1.25, alpha = 0.95
  ) +
  geom_point(
    data = d_mean,
    aes(x = mu, y = y, colour = scenario, shape = scenario),
    size = 2.2, stroke = 0.25
  ) +
  scale_colour_manual(values = scenario_cols, labels = scenario_labels) +
  scale_shape_manual(values = scenario_shapes, labels = scenario_labels) +
  scale_y_continuous(breaks = seq_along(duration_levels), labels = duration_labels) +
  labs(title = "Stage durations", x = "Duration (years)", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 7.5),
    axis.title.x = element_text(size = 8.5),
    plot.title = element_text(size = 9.5, face = "bold"),
    legend.position = "none"
  )

# ----------------------------------------------------------------------
# 3. Selected burial posterior distributions
# ----------------------------------------------------------------------
burials <- c("NM30", "NM27", "NM22", "NM23")
posterior_all <- do.call(rbind, lapply(scenario_ids, function(sc) {
  z <- read_required(file.path(project_dir, paste0("north_R_primary_", sc), "burial_posterior_density.csv"))
  z <- z[z$burial %in% burials, , drop = FALSE]
  z$scenario <- sc
  z$cal_BCE <- -as.numeric(z$calendar_year)
  z
}))
posterior_all$burial <- factor(posterior_all$burial, levels = burials)

p_post <- ggplot(
  posterior_all,
  aes(x = cal_BCE, y = kinship_posterior_density, colour = scenario, fill = scenario, group = scenario)
) +
  geom_area(alpha = 0.07, position = "identity") +
  geom_line(linewidth = 0.65) +
  facet_grid(burial ~ ., scales = "free_y", switch = "y") +
  scale_x_reverse(limits = c(3145, 2945)) +
  scale_colour_manual(values = scenario_cols, labels = scenario_labels, name = NULL) +
  scale_fill_manual(values = scenario_cols, labels = scenario_labels, name = NULL) +
  labs(title = "Selected burial-event posteriors", x = "Calendar age (cal BCE)", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 7.5),
    axis.title.x = element_text(size = 8.5),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 8.5),
    plot.title = element_text(size = 9.5, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 8),
    legend.key.width = grid::unit(8, "mm")
  )

left_panel <- p_boundaries / p_durations + plot_layout(heights = c(1.35, 1))
final_plot <- left_panel | p_post
final_plot <- final_plot +
  plot_layout(widths = c(1.05, 1.35)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12),
      plot.caption = element_text(size = 7, colour = "grey35", hjust = 0)
    ),
    caption = "Thin segments: 95.4% HPD; thick segments: 68.3% HPD; symbols: posterior means."
  )

out_pdf <- file.path(project_dir, "Gangshang_North_kinship_sensitivity_figure.pdf")
out_tif <- file.path(project_dir, "Gangshang_North_kinship_sensitivity_figure_600dpi.tiff")
out_png <- file.path(project_dir, "Gangshang_North_kinship_sensitivity_figure_preview.png")

ggsave(out_pdf, final_plot, width = 7.2, height = 7.8, units = "in", device = cairo_pdf)
ggsave(out_tif, final_plot, width = 7.2, height = 7.8, units = "in", dpi = 600, compression = "lzw")
ggsave(out_png, final_plot, width = 7.2, height = 7.8, units = "in", dpi = 250)

cat("Sensitivity figure written to:\n")
cat(out_pdf, "\n", out_tif, "\n", out_png, "\n", sep = "")
