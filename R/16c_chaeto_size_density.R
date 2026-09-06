# =============================================================================
# Per-cell size distributions of the Chaetoceros chain class at 3 and 16 m, all ROIs
# and single-cell ROIs: density modes, quantiles, bounding-box width modes and a
# two-component mixture (Fig. S6).
#
# Reads:  data/roi_level/roi_2016.txt.gz, outputs/tangesund/02_sample_meta.csv.
# Writes: outputs/tangesund/16c_chaeto_size_density.csv, 16c_chaeto_size_density.{png,pdf}.
# Deps:   data.table, ggplot2, patchwork.
# Run:    Rscript R/16c_chaeto_size_density.R
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("R/lib_plotting_theme.R")

cfg <- list(
  roi_file = "data/roi_level/roi_2016.txt.gz", meta = "outputs/tangesund/02_sample_meta.csv",
  out_dir = "outputs/tangesund", depths = c(3, 16), um_per_px = 2.77,
  esd_cut = 8, esd_lim = c(3, 25), dens_adjust = 1.5, plot_adjust = 1.2,
  probs = c(.05, .25, .5, .75, .95)
)

meta <- fread(cfg$meta)[qc_pass == TRUE, .(sample, depth, date)]
d <- fread(cfg$roi_file,
           select = c("sample", "class_name", "biovolume_um3", "numBlobs", "cell_count_resolved",
                      "BoundingBox_xwidth", "BoundingBox_ywidth"),
           showProgress = FALSE)
d <- merge(d[grepl("^Chaetoceros", class_name)], meta, by = "sample")
d <- d[depth %in% cfg$depths & numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 &
         is.finite(cell_count_resolved) & cell_count_resolved >= 1]
chain_cls <- d[, .N, by = class_name][order(-N)][1, class_name]
cat("Chaetoceros chain class:", chain_cls, "\n")
x <- d[class_name == chain_cls]
x[, `:=`(esd = (6 * biovolume_um3 / cell_count_resolved / pi)^(1/3),
         lbv = log10(biovolume_um3 / cell_count_resolved),
         width_um = pmin(BoundingBox_xwidth, BoundingBox_ywidth) / cfg$um_per_px,
         depth_f = factor(paste(depth, "m"), levels = paste(cfg$depths, "m")),
         single = cell_count_resolved == 1)]
n <- x[, .(n_all = .N, n_single = sum(single), n_multi = sum(!single), n_days = uniqueN(date)), by = depth][order(depth)]
cat("ROIs by depth:\n"); print(n)

# ---- density modes and quantiles ------------------------------------------------
modes <- function(v, adj = cfg$dens_adjust, to_esd = TRUE) {
  dd <- density(v, adjust = adj, n = 1024)
  i <- which(diff(sign(diff(dd$y))) == -2) + 1
  i <- i[dd$y[i] > 0.05 * max(dd$y)]
  m <- if (to_esd) (6 * 10^dd$x[i] / pi)^(1/3) else dd$x[i]
  list(modes = m, rel = dd$y[i] / max(dd$y))
}
mode_rows <- rbindlist(lapply(cfg$depths, function(dp) rbindlist(list(
  { m <- modes(x[depth == dp, lbv]);            data.table(depth = dp, subset = "all", quantity = "mode_esd_um", value = m$modes, extra = m$rel) },
  { m <- modes(x[depth == dp & single, lbv]);   data.table(depth = dp, subset = "single", quantity = "mode_esd_um", value = m$modes, extra = m$rel) },
  { m <- modes(x[depth == dp & !single, lbv]);  data.table(depth = dp, subset = "multi", quantity = "mode_esd_um", value = m$modes, extra = m$rel) },
  { m <- modes(x[depth == dp, width_um], to_esd = FALSE)
    data.table(depth = dp, subset = "all", quantity = "mode_bbox_width_um", value = m$modes, extra = m$rel) },
  { m <- modes(x[depth == dp & !single, width_um], to_esd = FALSE)
    data.table(depth = dp, subset = "multi", quantity = "mode_bbox_width_um", value = m$modes, extra = m$rel) }))))
q_rows <- rbindlist(lapply(cfg$depths, function(dp) rbind(
  data.table(depth = dp, subset = "all", quantity = paste0("esd_q", 100 * cfg$probs), value = quantile(x[depth == dp, esd], cfg$probs), extra = NA_real_),
  data.table(depth = dp, subset = "single", quantity = paste0("esd_q", 100 * cfg$probs), value = quantile(x[depth == dp & single, esd], cfg$probs), extra = NA_real_))))
qs <- x[single == TRUE, .(q = quantile(esd, cfg$probs), p = cfg$probs), by = depth]
ratio_rows <- merge(qs[depth == 3], qs[depth == 16], by = "p", suffixes = c("_3", "_16"))[
  , .(depth = NA_real_, subset = "single", quantity = paste0("ratio_3m_over_16m_q", 100 * p), value = q_3 / q_16, extra = NA_real_)]
share_rows <- x[, .(subset = "all", quantity = paste0("share_esd_over_", cfg$esd_cut, "um"), value = mean(esd > cfg$esd_cut), extra = NA_real_), by = depth]
n_rows <- melt(n, id.vars = "depth", variable.name = "quantity", value.name = "value")[, `:=`(subset = "all", extra = NA_real_)]

# ---- two-component mixture on log10 per-cell biovolume -----------------------------
em2 <- function(v, iter = 200) {
  mu <- as.numeric(quantile(v, c(.25, .75))); s <- rep(sd(v) / 2, 2); p <- c(.5, .5)
  for (k in seq_len(iter)) {
    r <- cbind(p[1] * dnorm(v, mu[1], s[1]), p[2] * dnorm(v, mu[2], s[2])); r <- r / rowSums(r)
    p <- colMeans(r); mu <- colSums(r * v) / colSums(r)
    s <- sqrt(colSums(r * (v - rep(mu, each = length(v)))^2) / colSums(r))
  }
  ll <- sum(log(p[1] * dnorm(v, mu[1], s[1]) + p[2] * dnorm(v, mu[2], s[2])))
  o <- order(mu)
  list(p = p[o], mu = mu[o], s = s[o], mode_esd = (6 * 10^mu[o] / pi)^(1/3), bic = -2 * ll + 5 * log(length(v)))
}
mix_rows <- rbindlist(lapply(cfg$depths, function(dp) {
  v <- x[depth == dp, lbv]; f <- em2(v)
  bic1 <- -2 * sum(dnorm(v, mean(v), sd(v), log = TRUE)) + 2 * log(length(v))
  data.table(depth = dp, subset = "all",
             quantity = c("mix_weight_small", "mix_weight_large", "mix_mode_small_esd_um", "mix_mode_large_esd_um", "mix_dBIC_1_minus_2"),
             value = c(f$p, f$mode_esd, bic1 - f$bic), extra = NA_real_)
}))
fp <- em2(x$lbv)
r <- cbind(fp$p[1] * dnorm(x$lbv, fp$mu[1], fp$s[1]), fp$p[2] * dnorm(x$lbv, fp$mu[2], fp$s[2]))
x[, resp_small := r[, 1] / rowSums(r)]
pooled_rows <- rbind(
  data.table(depth = NA_real_, subset = "all", quantity = c("pooled_mix_mode_small_esd_um", "pooled_mix_mode_large_esd_um"), value = fp$mode_esd, extra = NA_real_),
  x[, .(subset = "all", quantity = "pooled_mix_share_small_component", value = mean(resp_small), extra = NA_real_), by = depth])

out <- rbind(n_rows, mode_rows, q_rows, ratio_rows, share_rows, mix_rows, pooled_rows, use.names = TRUE)
out[, `:=`(class = chain_cls, value = round(value, 4), extra = round(extra, 3))]
setcolorder(out, c("class", "depth", "subset", "quantity", "value", "extra"))
fwrite(out, file.path(cfg$out_dir, "16c_chaeto_size_density.csv"))
cat("\n-- 16c_chaeto_size_density.csv --\n"); print(out[, .(depth, subset, quantity, value, extra)], nrows = 200)

# ---- figure (Fig. S6) --------------------------------------------------------------
pal <- c("3 m" = "firebrick", "16 m" = "steelblue")
lab_n <- function(sub) {
  nn <- sub[, .N, by = depth_f][order(depth_f)]
  setNames(sprintf("%s (n = %s)", nn$depth_f, format(nn$N, big.mark = ",")), nn$depth_f)
}
dens_panel <- function(sub, xlab) {
  ggplot(sub, aes(esd, colour = depth_f)) +
    geom_vline(xintercept = cfg$esd_cut, linetype = 2, colour = "grey50") +
    geom_density(adjust = cfg$plot_adjust, linewidth = 0.8, key_glyph = "path") +
    scale_x_log10(breaks = c(3, 5, 8, 12, 20)) + coord_cartesian(xlim = cfg$esd_lim) +
    scale_colour_manual(values = pal, labels = lab_n(sub), name = NULL) +
    theme_paper() + theme(legend.position = "top", legend.key.width = unit(1.2, "lines")) +
    labs(x = xlab, y = "Density")
}
pa <- dens_panel(x, expression(paste("Per-cell ESD (", mu, "m), all chain-class ROIs")))
pb <- dens_panel(x[single == TRUE], expression(paste("Per-cell ESD (", mu, "m), single-cell ROIs")))
ggsave(file.path(cfg$out_dir, "16c_chaeto_size_density.png"),
       ((pa + labs(tag = "a")) | (pb + labs(tag = "b"))) &
         theme(plot.tag = element_text(face = "bold", size = 13)),
       width = 6.5, height = 3.2, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "16c_chaeto_size_density.pdf"),
       ((pa + labs(tag = "a")) | (pb + labs(tag = "b"))) &
         theme(plot.tag = element_text(face = "bold", size = 13)),
       width = 6.5, height = 3.2, device = cairo_pdf)
cat(sprintf("\nfigure -> %s/16c_chaeto_size_density.png (+ .pdf)\n", cfg$out_dir))
