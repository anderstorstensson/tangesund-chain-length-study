# Shared figure style: ggplot theme, season palette and order, raster resolution.
#
# Reads nothing, writes nothing. Deps: ggplot2.
# Use: source("R/lib_plotting_theme.R")
# Figure width and height are set by 28_render_submission_figures.R only.

suppressPackageStartupMessages(library(ggplot2))

# ---- Season -----------------------------------------------------------------
SEASON_LEVELS <- c("Winter", "Spring", "Summer", "Autumn")

# Okabe-Ito colors, chosen for the largest minimum CIE76 distance between seasons
# under simulated protanopia, deuteranopia and tritanopia (Vienot et al. 1999).
SEASON_COLOURS <- c(Winter = "#0072B2",   # blue
                    Spring = "#009E73",   # bluish green
                    Summer = "#E69F00",   # orange
                    Autumn = "#CC79A7")   # reddish purple

season_factor <- function(x) factor(as.character(x), levels = SEASON_LEVELS)

# limits + drop = FALSE keep each season on its color when a panel lacks a season
scale_colour_season <- function(name = "Season", ...) {
  scale_colour_manual(values = SEASON_COLOURS, limits = SEASON_LEVELS,
                      drop = FALSE, name = name, ...)
}
scale_color_season <- scale_colour_season
scale_fill_season <- function(name = "Season", ...) {
  scale_fill_manual(values = SEASON_COLOURS, limits = SEASON_LEVELS,
                    drop = FALSE, name = name, ...)
}

# ---- Theme ------------------------------------------------------------------
# 11 pt for every figure, composites included
THEME_BASE_SIZE <- 11

theme_paper <- function(base_size = THEME_BASE_SIZE, ...) {
  theme_bw(base_size = base_size, ...) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text       = element_text(size = rel(0.9)),
      legend.key       = element_blank(),
      plot.title       = element_text(size = rel(1.0), face = "plain"),
      plot.tag         = element_text(size = rel(1.1), face = "bold")
    )
}

# ---- Output -----------------------------------------------------------------
# Working-figure resolution; 600-dpi rasters come from 28_render_submission_figures.R
FIG_DPI <- 150

# Nominal column widths (in), not enforced
FIG_W_SINGLE <- 3.5
FIG_W_DOUBLE <- 7.2

save_fig <- function(filename, plot = ggplot2::last_plot(), width, height,
                     dpi = FIG_DPI, ...) {
  ggsave(filename, plot, width = width, height = height, dpi = dpi, ...)
}
