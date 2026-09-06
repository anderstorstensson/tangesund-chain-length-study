#!/usr/bin/env Rscript
# =============================================================================
# Vector PDF and 600-dpi PNG of each manuscript figure. Each panel script is re-run
# in an R subprocess with a ggsave() shim that drops panel titles and writes the
# figure files, sized so axis text lands at about 6 pt on the page.
#
# Reads:  the panel scripts in the table below and their inputs.
# Writes: outputs/tangesund/figures/Fig1-Fig7 and FigS1, S2, S4-S8 (.pdf + .png).
#         Fig. S3 is written there by 15b_single_cell_gallery.py.
# Env:    FIG_ONLY="Fig4,FigS5" renders a subset.
# Deps:   ggplot2 plus whatever the panel scripts load.
# Run:    Rscript R/28_render_submission_figures.R   (after all panel scripts)
# =============================================================================

figdir <- "outputs/tangesund/figures"
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

# script, panel PNG it saves, manuscript label (Fig. S3 comes from 15b)
fig <- rbind(
  c("R/02_coverage.R",             "02_coverage_panel.png",     "Fig1"),
  c("R/03_chain_length.R",         "03_chain_panel.png",        "Fig2"),
  c("R/17_plasticity_turnover.R", "17_plasticity_panel.png",  "Fig3"),
  c("R/16_percell_size.R",        "16_percell_panel.png",     "Fig4"),
  c("R/05_composition.R",          "05_composition_panel.png",  "Fig5"),
  c("R/10_env_mechanism.R",        "10_mechanism_panel.png",    "Fig6"),
  c("R/08_diel_chain.R",           "08_diel_panel.png",         "Fig7"),
  c("R/06_env_buoy.R",             "06_buoy_panel.png",         "FigS1"),
  c("R/27c_day_bootstrap.R",      "27c_contrasts_panel.png",   "FigS2"),
  c("R/10b_october_halocline.R",  "10b_october_panel.png",     "FigS4"),
  c("R/09b_regional_salinity.R",  "09b_svea_salinity_panel.png", "FigS5"),
  c("R/16c_chaeto_size_density.R", "16c_chaeto_size_density.png", "FigS6"),
  c("R/23i_validation_figure.R",  "23i_validation_panel.png",  "FigS7"),
  c("R/27d_si_robustness.R",      "27d_robustness_panel.png",  "FigS8")
)
colnames(fig) <- c("script", "png", "label")

shim <- '
library(ggplot2)
.rg <- ggplot2::ggsave
.PNG <- Sys.getenv("FIG_PNG"); .LAB <- Sys.getenv("FIG_LABEL"); .DIR <- Sys.getenv("FIG_DIR")
# drop panel titles in every panel of a (nested) patchwork
.small_titles <- function(p) {
  if (inherits(p, "patchwork"))
    p$patches$plots <- lapply(p$patches$plots, .small_titles)
  p + ggplot2::theme(plot.title = ggplot2::element_blank())
}
# Deposit canvas scaled so axis text lands at .TARGET_PT when the figure is placed in a
# 6.5 x 8.5 in text block; the working PNG keeps the size set in the panel script.
.DOCX_W <- 6.5; .DOCX_H <- 8.5
.BASE <- 11; .AXIS_REL <- 0.8      # THEME_BASE_SIZE in R/lib_plotting_theme.R; axis text rel(0.8)
.TARGET_PT <- 6
# capped at 1: figures already at page width (S2, S6) are never enlarged
.deposit_scale <- function(w, h) {
  placed <- min(.DOCX_W, .DOCX_H * w / h)                 # width the manuscript places it at
  min(1, placed * (.BASE * .AXIS_REL / .TARGET_PT) / w)   # ggsave multiplies both dims by this
}
# A smaller canvas also enlarges legends; shrunk figures get compacted legends.
# Resize elements in place: ggplot2 will not merge element_text over element_markdown.
.set_size <- function(el, r) {
  if (is.null(el)) return(ggplot2::element_text(size = ggplot2::rel(r)))
  # element_blank has no size property
  if (inherits(el, "element_blank")) return(el)
  el$size <- ggplot2::rel(r)
  el
}
.compact_legends <- function(p) {
  if (inherits(p, "patchwork"))
    p$patches$plots <- lapply(p$patches$plots, .compact_legends)
  p + ggplot2::theme(
    axis.title         = .set_size(p$theme$axis.title, 0.85),
    legend.key.size    = grid::unit(0.5, "lines"),
    legend.text        = .set_size(p$theme$legend.text,  0.60),
    legend.title       = .set_size(p$theme$legend.title, 0.65),
    legend.margin      = ggplot2::margin(1, 1, 1, 1),
    legend.box.spacing = grid::unit(3, "pt"))
}
ggsave <- function(filename, plot = ggplot2::last_plot(), ...) {
  plot <- .small_titles(plot)
  .rg(filename, plot, ...)
  if (basename(filename) == .PNG) {
    a <- list(...)
    if (is.null(a$width) || is.null(a$height))
      stop(sprintf("%s: the panel is saved without an explicit width and height, so the
 deposit cannot be sized for the page; add both to its ggsave() call", .LAB))
    sc <- .deposit_scale(a$width, a$height)
    # compact legends only for figures being shrunk onto the page
    dep <- if (sc < 1) .compact_legends(plot) else plot
    dims <- list(plot = dep, width = a$width, height = a$height, scale = sc)
    do.call(.rg, c(list(file.path(.DIR, paste0(.LAB, ".pdf")), device = grDevices::cairo_pdf), dims))
    do.call(.rg, c(list(file.path(.DIR, paste0(.LAB, ".png")), dpi = 600), dims))
    cat(sprintf("  -> %s.pdf / %s.png  (%.2f x %.2f in, scale %.3f)\\n",
                .LAB, .LAB, a$width * sc, a$height * sc, sc))
  }
}
source(Sys.getenv("FIG_SCRIPT"))
'
shimfile <- tempfile(fileext = ".R"); writeLines(shim, shimfile)

# FIG_ONLY: comma-separated labels; default is all
only <- Sys.getenv("FIG_ONLY")
if (nzchar(only)) fig <- fig[fig[, "label"] %in% trimws(strsplit(only, ",")[[1]]), , drop = FALSE]

ok <- 0L
for (i in seq_len(nrow(fig))) {
  cat(sprintf("[%2d/%2d] %-28s -> %s\n", i, nrow(fig), fig[i, "script"], fig[i, "label"]))
  st <- system2("Rscript", shimfile,
    env = c(sprintf("FIG_PNG=%s",    fig[i, "png"]),
            sprintf("FIG_LABEL=%s",  fig[i, "label"]),
            sprintf("FIG_DIR=%s",    figdir),
            sprintf("FIG_SCRIPT=%s", fig[i, "script"]),
            "SVEA_FIG_ONLY=1",    # 09b redraws Fig. S5 from its CSVs instead of refitting (30-60 min)
            "HALOCLINE_FIG_ONLY=1",  # 10b redraws Fig. S4 from its CSVs instead of bootstrapping
            "BOOT_FIG_ONLY=1"),   # 27c redraws Fig. S2 from its CSV instead of bootstrapping
    stdout = FALSE, stderr = FALSE)
  if (identical(st, 0L)) ok <- ok + 1L else cat("   ! FAILED (exit", st, ")\n")
}
cat(sprintf("\nDone: %d/%d figures rendered to %s\n", ok, nrow(fig), figdir))
