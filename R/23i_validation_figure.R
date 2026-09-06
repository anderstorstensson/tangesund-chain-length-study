# =============================================================================
# The two blind manual counts in one figure (Fig. S7) and the rows added to Table S5.
# Nothing is refitted except a day-clustered bootstrap of the mean count error per
# stratum, which 23e and 23h do not tabulate.
#  (a) detector against manual count per chain, both taxa, by depth
#  (b) mean count error (detector minus manual) by size bin and depth
#  (c) the weighted differential over-count next to the surface-minus-deep contrast
#  (d) unclassified non-junk ROIs and chain-eligible ROIs with no detected cell, by depth
# Frames flagged "not the taxon" or "cells not resolvable" are excluded.
#
# Reads:  data/chaetoceros_validation/chaetoceros_manual_counts.csv,
#         data/asterionellopsis_validation/asterionellopsis_manual_counts.csv,
#         outputs/tangesund/chaetoceros_validation/23e_differential_bias.csv, 23e_class_purity.csv,
#         outputs/tangesund/asterionellopsis_validation/23h_differential_bias.csv, 23h_class_purity.csv,
#         outputs/tangesund/23b_detection_by_depth.csv.
# Writes: outputs/tangesund/23i_bias_by_stratum_boot.csv, 23i_pooled_accuracy.csv,
#         23i_validation_rows.md (Table S5 rows), 23i_validation_panel.png/.pdf (Fig. S7).
# Deps:   data.table, ggplot2, patchwork, ggtext.
# Run:    Rscript R/23i_validation_figure.R   (after 23b, 23e and 23h)
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(ggtext) })
source("R/lib_plotting_theme.R")
cfg <- list(out_dir = "outputs/tangesund", B = 2000, seed = 41,
            taxa = list(
              list(key = "chaeto", label = "Chaetoceros", counts = "data/chaetoceros_validation/chaetoceros_manual_counts.csv",
                   bin = "width_bin", bin_label = "ROI width (µm)", bins = c("<=8", "8-12", "12-18", ">18"),
                   diff = "outputs/tangesund/chaetoceros_validation/23e_differential_bias.csv",
                   purity = "outputs/tangesund/chaetoceros_validation/23e_class_purity.csv"),
              list(key = "ag", label = "A. glacialis", counts = "data/asterionellopsis_validation/asterionellopsis_manual_counts.csv",
                   bin = "feret_bin", bin_label = "ROI maximum Feret diameter (µm)", bins = c("<=18", "18-34", "34-98", ">98"),
                   diff = "outputs/tangesund/asterionellopsis_validation/23h_differential_bias.csv",
                   purity = "outputs/tangesund/asterionellopsis_validation/23h_class_purity.csv")))
pal_depth <- c(`3` = "firebrick", `16` = "steelblue4")

kept <- rbindlist(lapply(cfg$taxa, function(tx) {
  d <- fread(tx$counts)
  d <- d[flags == ""]
  d[, .(taxon = tx$label, anon_id, depth, date = as.Date(date), bin = factor(get(tx$bin), levels = tx$bins),
        detector = detector_count, manual = manual_count, err = detector_count - manual_count)]
}))
# the lumped class is "Chaetoceros spp." as in the other figures; the species stays a binomial
strip_lab <- c(Chaetoceros = "*Chaetoceros* spp.", `A. glacialis` = "*A. glacialis*")
kept[, taxon := factor(taxon, levels = names(strip_lab))]
kept[, strip := factor(strip_lab[as.character(taxon)], levels = strip_lab)]
kept[, depth_f := factor(depth, levels = c(3, 16))]

# ---- (b) day-clustered bootstrap of the mean error per stratum ------------------
boot_days <- function(err, day, B = cfg$B) {
  days <- unique(day); idx <- split(seq_along(err), day)
  draws <- vapply(seq_len(B), function(b) mean(err[unlist(idx[as.character(sample(days, length(days), TRUE))], use.names = FALSE)]), numeric(1))
  c(lo = unname(quantile(draws, 0.025)), hi = unname(quantile(draws, 0.975)))
}
set.seed(cfg$seed)
strata <- kept[, {
  q <- boot_days(err, as.character(date))
  .(n = .N, days = uniqueN(date), bias = mean(err), lo = q[["lo"]], hi = q[["hi"]],
    exact_pct = 100 * mean(err == 0), within1_pct = 100 * mean(abs(err) <= 1), mae = mean(abs(err)))
}, by = .(taxon, strip, depth_f, bin)][order(taxon, depth_f, bin)]
fwrite(strata[, !"strip"], file.path(cfg$out_dir, "23i_bias_by_stratum_boot.csv"))
cat("-- count error by stratum with day-clustered bootstrap intervals (23i_bias_by_stratum_boot.csv) --\n")
print(strata[, .(taxon, depth = depth_f, bin, n, bias = round(bias, 2), CI = sprintf("[%.2f, %.2f]", lo, hi))])

# ---- pooled accuracy on the kept chains, per taxon (Table S5 rows) ---------------
set.seed(cfg$seed + 1)
pooled <- kept[, {
  q <- boot_days(err, as.character(date))
  .(n = .N, days = uniqueN(date), exact_pct = 100 * mean(err == 0), within1_pct = 100 * mean(abs(err) <= 1),
    mae = mean(abs(err)), bias = mean(err), bias_lo = q[["lo"]], bias_hi = q[["hi"]])
}, by = taxon]
fwrite(pooled, file.path(cfg$out_dir, "23i_pooled_accuracy.csv"))
cat("\n-- pooled accuracy on the kept chains (23i_pooled_accuracy.csv) --\n"); print(pooled)

pur <- rbindlist(lapply(cfg$taxa, function(tx) {
  p <- fread(tx$purity)
  p[, .(taxon = tx$label, depth, not_pct = get(grep("^not_.*_pct$", names(p), value = TRUE)[1]),
        contrast_raw, contrast_purified, n_flagged = get(grep("^n_flagged", names(p), value = TRUE)[1]), n_drawn)]
}))
minus <- function(s) gsub("-", "−", s, fixed = TRUE)
mrow <- function(...) paste0("| ", paste(c(...), collapse = " | "), " |")
g <- function(tx) pooled[taxon == tx]
p3 <- function(tx) pur[taxon == tx & depth == 3]; p16 <- function(tx) pur[taxon == tx & depth == 16]
rows <- c(
  "## Supplementary Table 5: rows appended", "",
  mrow("Quantity", "*Chaetoceros*", "95% CI", "*A. glacialis*", "95% CI"), "|:--|--:|:--:|--:|:--:|",
  mrow("Frames flagged not the taxon, 3 m (% of the imaged class, weighted by size bin)", sprintf("%.0f", p3("Chaetoceros")$not_pct), "", sprintf("%.0f", p3("A. glacialis")$not_pct), ""),
  mrow("Frames flagged not the taxon, 16 m (%)", sprintf("%.0f", p16("Chaetoceros")$not_pct), "", sprintf("%.0f", p16("A. glacialis")$not_pct), ""),
  mrow("Surface-minus-deep contrast, all imaged ROIs, after removing such frames at their detector counts (cells)",
       sprintf("+%.2f", p3("Chaetoceros")$contrast_purified), "", sprintf("+%.2f", p3("A. glacialis")$contrast_purified), ""),
  mrow("Kept chains counted exactly (%)", sprintf("%.0f", g("Chaetoceros")$exact_pct), "", sprintf("%.0f", g("A. glacialis")$exact_pct), ""),
  mrow("Kept chains counted within one cell (%)", sprintf("%.0f", g("Chaetoceros")$within1_pct), "", sprintf("%.0f", g("A. glacialis")$within1_pct), ""),
  mrow("Mean absolute error on the kept chains (cells)", sprintf("%.2f", g("Chaetoceros")$mae), "", sprintf("%.2f", g("A. glacialis")$mae), ""),
  mrow("Net bias on the kept chains, unweighted (cells)", minus(sprintf("%.2f", g("Chaetoceros")$bias)), minus(sprintf("[%.2f, %.2f]", g("Chaetoceros")$bias_lo, g("Chaetoceros")$bias_hi)),
       minus(sprintf("%+.2f", g("A. glacialis")$bias)), minus(sprintf("[%.2f, %.2f]", g("A. glacialis")$bias_lo, g("A. glacialis")$bias_hi))),
  "")
writeLines(rows, file.path(cfg$out_dir, "23i_validation_rows.md"))
cat(sprintf("\nrows -> %s/23i_validation_rows.md\n", cfg$out_dir))

# ---- figure (Fig. S7) ----------------------------------------------------------
pa <- ggplot(kept, aes(manual, detector, colour = depth_f)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_jitter(width = 0.15, height = 0.15, alpha = 0.6, size = 1.4) +
  facet_wrap(~ strip, scales = "free") +
  scale_colour_manual(values = pal_depth, name = "Depth (m)") +
  theme_paper() + theme(legend.position = "top", strip.text = element_markdown()) +
  labs(x = "Manual count (cells per chain)", y = "Detector count (cells per chain)", title = "Detector against manual count")
pb <- ggplot(strata, aes(bias, bin, colour = depth_f)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.25, position = position_dodge(width = 0.5)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.5)) +
  facet_wrap(~ strip, scales = "free") +
  scale_colour_manual(values = pal_depth, name = "Depth (m)") +
  theme_paper() + theme(legend.position = "top", strip.text = element_markdown()) +
  labs(x = "Count error, detector − manual (cells; 95% CI)", y = "Size bin (µm)", title = "Count error by size bin and depth")
dif <- rbindlist(lapply(cfg$taxa, function(tx) fread(tx$diff)[, .(taxon = tx$label, quantity, estimate, lo95, hi95)]))
qlab <- c(diff_bias = "Differential\nover-count", contrast_manual = "Contrast,\nmanual counts", contrast_detector_same_chains = "Contrast, detector\non the same chains",
          contrast_raw = "Contrast, all ROIs,\nuncorrected", contrast_corrected = "Contrast, all ROIs,\nbias-corrected")
dif <- dif[quantity %in% names(qlab)]
dif[, quantity := factor(qlab[quantity], levels = rev(qlab))]
dif[, strip := factor(strip_lab[taxon], levels = strip_lab)]
pc <- ggplot(dif, aes(estimate, quantity)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo95, xmax = hi95), orientation = "y", width = 0.25, colour = "grey30") +
  geom_point(size = 2.2, colour = "firebrick") +
  facet_wrap(~ strip, scales = "free_x") +
  theme_paper() + theme(strip.text = element_markdown()) +
  labs(x = "Cells per chain (day-clustered bootstrap 95% CI)", y = NULL, title = "The differential against the contrast")
cen <- fread(file.path(cfg$out_dir, "23b_detection_by_depth.csv"))[!is.na(depth)]
cen <- melt(cen[, .(depth, `Unclassified non-junk ROIs` = unclassified_pct, `Chain-eligible ROIs with no detected cell` = zero_detection_pct)], id.vars = "depth")
pd <- ggplot(cen, aes(value, depth, colour = variable, shape = variable)) +
  geom_path(orientation = "y") + geom_point(size = 2) +
  scale_y_reverse(breaks = c(3, 6, 8, 11, 13, 16)) +
  scale_colour_manual(values = c("grey30", "steelblue4"), name = NULL) + scale_shape_manual(values = c(16, 17), name = NULL) +
  theme_paper() + theme(legend.position = "top", legend.direction = "vertical") +
  labs(x = "Share of ROIs (%)", y = "Depth (m)", title = "Censoring by depth")
fig <- ((pa + labs(tag = "a")) | (pb + labs(tag = "b"))) / (free(pc + labs(tag = "c")) | (pd + labs(tag = "d"))) +
  plot_layout(widths = c(1, 1), heights = c(1, 1)) & theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(file.path(cfg$out_dir, "23i_validation_panel.png"), fig, width = 12, height = 9, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "23i_validation_panel.pdf"), fig, width = 12, height = 9, device = cairo_pdf)
cat(sprintf("figure -> %s/23i_validation_panel.{png,pdf}\n", cfg$out_dir))
