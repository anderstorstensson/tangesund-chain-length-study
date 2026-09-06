# =============================================================================
# Chaetoceros validation set, stage 3 of 3: depth-differential count bias from the
# manual counts (manual count = "cell" boxes per task), post-stratified to each depth's
# width composition, with day-clustered bootstrap intervals (Tables S4, S5); class
# purity and its effect on the chain-length and per-cell size contrasts.
# Frames flagged "not Chaetoceros" or "cells not resolvable" are excluded.
#
# Reads:  the Label Studio export (default data/chaetoceros_validation/label_studio_export_project24.json),
#         the bundled data/chaetoceros_validation/ key and population bins (copied from 23c's
#         outputs if absent), data/roi_level/roi_2016.txt.gz, outputs/tangesund/02_sample_meta.csv
# Writes: data/chaetoceros_validation/ chaetoceros_manual_counts.csv, chaetoceros_cell_boxes.csv;
#         outputs/tangesund/chaetoceros_validation/ 23e_bias_by_stratum.csv, 23e_width_quartiles.csv,
#         23e_class_purity.csv, 23e_class_purity_size.csv, 23e_differential_bias.csv
# Deps:   data.table, jsonlite
# Run:    Rscript R/23e_chaetoceros_analysis.R [export.json]   (no Label Studio server needed)
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

cfg <- list(
  work_dir = "outputs/tangesund/chaetoceros_validation",
  data_dir = "data/chaetoceros_validation",            # tracked: the hand-made counts cannot be regenerated
  roi_file = "data/roi_level/roi_2016.txt.gz",
  meta_csv = "outputs/tangesund/02_sample_meta.csv",
  micron   = 1 / 2.77,
  width_breaks = c(0, 8, 12, 18, Inf),        # um, as R/23c_chaetoceros_select.R
  n_boot   = 2000, seed = 1
)
export <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(export)) export <- file.path(cfg$data_dir, "label_studio_export_project24.json")
if (!file.exists(export)) stop("usage: Rscript R/23e_chaetoceros_analysis.R [label-studio export.json]")
for (p in c(cfg$data_dir, cfg$work_dir)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# ---- manual counts from the export -------------------------------------------
tasks <- fromJSON(export, simplifyVector = FALSE)
manual <- rbindlist(lapply(tasks, function(t) {
  img <- t$data$image
  if (is.null(img)) return(NULL)
  id <- sub("\\.png$", "", basename(sub(".*[?&]d=", "", utils::URLdecode(img))))
  ann <- Filter(function(a) !isTRUE(a$was_cancelled), t$annotations)
  if (!length(ann)) return(data.table(anon_id = id, manual = NA_integer_, flags = ""))
  a <- ann[[length(ann)]]                                    # latest completed annotation
  res <- a$result
  n_box <- sum(vapply(res, function(r) identical(r$type, "rectanglelabels"), logical(1)))
  fl <- unlist(lapply(res, function(r) if (identical(r$type, "choices")) r$value$choices))
  data.table(anon_id = id, manual = as.integer(n_box), flags = paste(fl, collapse = "; "))
}))
boxes <- rbindlist(lapply(tasks, function(t) {
  img <- t$data$image
  ann <- Filter(function(a) !isTRUE(a$was_cancelled), t$annotations)
  if (is.null(img) || !length(ann)) return(NULL)
  id <- sub("\\.png$", "", basename(sub(".*[?&]d=", "", utils::URLdecode(img))))
  a <- ann[[length(ann)]]
  rb <- Filter(function(r) identical(r$type, "rectanglelabels"), a$result)
  if (!length(rb)) return(NULL)
  rbindlist(lapply(rb, function(r) data.table(
    anon_id = id, image_width_px = r$original_width, image_height_px = r$original_height,
    x_pct = r$value$x, y_pct = r$value$y, width_pct = r$value$width, height_pct = r$value$height)))
}))
key_file <- file.path(cfg$data_dir, "chaetoceros_selection_key.csv")
if (!file.exists(key_file)) file.copy(file.path(cfg$work_dir, "23c_selection_key.csv"), key_file)
key <- fread(key_file)
key[, date := as.Date(date)]
v <- merge(key, manual, by = "anon_id", all.x = TRUE)
setnames(v, c("det_count", "manual"), c("detector_count", "manual_count"))
fwrite(v[order(anon_id)], file.path(cfg$data_dir, "chaetoceros_manual_counts.csv"))      # one row per blinded chain
fwrite(boxes[order(anon_id)], file.path(cfg$data_dir, "chaetoceros_cell_boxes.csv"))       # one row per manually drawn cell
setnames(v, c("detector_count", "manual_count"), c("det_count", "manual"))
cat(sprintf("tasks in export: %d; annotated: %d of %d drawn\n", nrow(manual), sum(!is.na(v$manual)), nrow(v)))
v_all <- copy(v)                                                   # every drawn frame, flags included
excl <- v[grepl("not Chaetoceros|not resolvable", flags)]
cat(sprintf("excluded by flag: %d (%s)\n", nrow(excl), paste(table(excl$depth), collapse = " / ")))
v <- v[!is.na(manual) & !grepl("not Chaetoceros|not resolvable", flags)]
v[, err := det_count - manual]
v[, width_bin := factor(width_bin, levels = c("<=8", "8-12", "12-18", ">18"))]

# ---- (a) per stratum ------------------------------------------------------------
cell <- v[, .(n = .N, days = uniqueN(date), bias = mean(err), mae = mean(abs(err)),
              exact_pct = 100 * mean(err == 0), within1_pct = 100 * mean(abs(err) <= 1),
              mean_manual = mean(manual), mean_det = mean(det_count)), by = .(depth, width_bin)]
setorder(cell, depth, width_bin)
fwrite(cell, file.path(cfg$work_dir, "23e_bias_by_stratum.csv"))
cat("\n-- count error (detector - manual) by depth x short-axis width --\n")
print(cell[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

# ---- (b)/(c) post-stratified bias and the differential -------------------------
pop_file <- file.path(cfg$data_dir, "chaetoceros_population_bins.csv")
if (!file.exists(pop_file)) file.copy(file.path(cfg$work_dir, "23c_population_bins.csv"), pop_file)
pop <- fread(pop_file)
pop[, width_bin := factor(width_bin, levels = levels(v$width_bin))]
meta <- fread(cfg$meta_csv)[qc_pass == TRUE, .(sample, depth, date = as.Date(date))]
d <- fread(cfg$roi_file, select = c("sample", "roi_number", "class_name", "numBlobs", "biovolume_um3", "cell_count",
                                    "cell_count_resolved", "BoundingBox_xwidth", "BoundingBox_ywidth"), showProgress = FALSE)
d <- d[grepl("^Chaetoceros", class_name) & !is.na(cell_count) & cell_count != -1 &
       numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0]
d[, det_count := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1,
                         as.integer(round(cell_count_resolved)), 1L)]
d[, width_um := pmin(BoundingBox_xwidth, BoundingBox_ywidth) * cfg$micron]
d[, width_bin := cut(width_um, cfg$width_breaks, labels = levels(v$width_bin), right = TRUE)]
d <- merge(d, meta, by = "sample")
# pooled short-axis width quartiles
wq_of <- function(x, lab) data.table(population = lab, q25_um = round(quantile(x, .25), 2), q50_um = round(quantile(x, .5), 2),
                                     q75_um = round(quantile(x, .75), 2), n_roi = length(x))
wq <- rbind(wq_of(d$width_um, "all six depths"), wq_of(d[depth %in% c(3, 16), width_um], "3 m and 16 m"))
fwrite(wq, file.path(cfg$work_dir, "23e_width_quartiles.csv"))
cat("\n-- Chaetoceros short-axis width quartiles (um) --\n"); print(wq)
d <- d[depth %in% c(3, 16)]
day_mean <- d[, .(m = mean(det_count), n = .N), by = .(depth, date)]   # day-level means + weights (clustered)

post_strat <- function(vv, dd) {
  b <- vv[, .(bias = mean(err), man = mean(manual), det = mean(det_count)), by = .(depth, width_bin)]
  b <- merge(pop[depth %in% c(3, 16), .(depth, width_bin, w)], b, by = c("depth", "width_bin"), all.x = TRUE)
  b[is.na(bias), `:=`(bias = 0, man = 0, det = 0)]                 # an empty stratum contributes nothing
  ps <- b[, .(bias = sum(w * bias), man = sum(w * man), det = sum(w * det)), by = depth]
  raw <- dd[, .(m = weighted.mean(m, n)), by = depth]                # ROI-pooled mean
  c(bias_3 = ps[depth == 3, bias], bias_16 = ps[depth == 16, bias],
    diff_bias = ps[depth == 3, bias] - ps[depth == 16, bias],
    manual_3 = ps[depth == 3, man], manual_16 = ps[depth == 16, man],
    contrast_manual = ps[depth == 3, man] - ps[depth == 16, man],
    contrast_detector_same_chains = ps[depth == 3, det] - ps[depth == 16, det],
    contrast_raw = raw[depth == 3, m] - raw[depth == 16, m],
    contrast_corrected = (raw[depth == 3, m] - raw[depth == 16, m]) - (ps[depth == 3, bias] - ps[depth == 16, bias]))
}
est <- post_strat(v, day_mean)

# ---- (d) class purity on the blinded draw ------------------------------------------
# "not Chaetoceros" rate weighted by population width bins; the purified contrast
# removes those frames, at their detector counts, from the ROI-pooled depth means
pur <- v_all[, .(rate = mean(grepl("not Chaetoceros", flags)),
                 det_nc = mean(det_count[grepl("not Chaetoceros", flags)])), by = .(depth, width_bin)]
pur <- merge(pop[depth %in% c(3, 16), .(depth, width_bin, w)], pur, by = c("depth", "width_bin"))
pur[is.na(det_nc), det_nc := 0]
raw_pooled <- day_mean[, .(pooled_mean = weighted.mean(m, n)), by = depth]
purity <- pur[, .(not_chaeto_pct = 100 * sum(w * rate),
                  mean_det_not_chaeto = sum(w * rate * det_nc) / sum(w * rate)), by = depth]
purity <- merge(purity, raw_pooled, by = "depth")
purity[, purified_mean := (pooled_mean - (not_chaeto_pct / 100) * mean_det_not_chaeto) / (1 - not_chaeto_pct / 100)]
purity[, `:=`(contrast_raw = pooled_mean[depth == 3] - pooled_mean[depth == 16],
              contrast_purified = purified_mean[depth == 3] - purified_mean[depth == 16])]
# the same quantities unweighted over the drawn frames (equal allocation)
pur_eq <- v_all[, .(not_chaeto_pct_equal_alloc = 100 * mean(grepl("not Chaetoceros", flags)),
                    mean_det_not_chaeto_equal_alloc = mean(det_count[grepl("not Chaetoceros", flags)]),
                    n_flagged_not_chaeto = sum(grepl("not Chaetoceros", flags)), n_drawn = .N), by = depth]
purity <- merge(purity, pur_eq, by = "depth")
fwrite(purity, file.path(cfg$work_dir, "23e_class_purity.csv"))
cat("\n-- class purity on the blinded draw (23e_class_purity.csv) --\n")
print(purity[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

# ---- (e) per-cell size of the flagged frames ---------------------------------
# log10 per-cell biovolume, weighted like the purity rates, with and without them
fr <- merge(v_all[, .(sample, roi_number, depth, width_bin, flags)],
            d[, .(sample, roi_number, lbv = log10(biovolume_um3 / det_count))], by = c("sample", "roi_number"))
fr[, not_ch := grepl("not Chaetoceros", flags)]
popsz <- d[, .(lbv_pop = mean(log10(biovolume_um3 / det_count))), by = .(depth, width_bin)]
flg <- fr[, .(rate = mean(not_ch), lbv_nc = mean(lbv[not_ch]), lbv_kept = mean(lbv[!not_ch])), by = .(depth, width_bin)]
sz <- merge(merge(pop[depth %in% c(3, 16), .(depth, width_bin, w)], popsz, by = c("depth", "width_bin")),
            flg, by = c("depth", "width_bin"))
sz[is.na(lbv_nc), lbv_nc := 0]
sz[, lbv_pur := (lbv_pop - rate * lbv_nc) / (1 - rate)]
esd_of <- function(l) (6 * 10^l / pi)^(1/3)
size <- sz[, .(log10_percell_bv_raw = sum(w * lbv_pop),
               log10_percell_bv_flagged = sum(w * rate * lbv_nc) / sum(w * rate),
               log10_percell_bv_kept_drawn = sum(w * lbv_kept),
               log10_percell_bv_purified = sum(w * lbv_pur)), by = depth]
size[, `:=`(esd_raw_um = esd_of(log10_percell_bv_raw), esd_flagged_um = esd_of(log10_percell_bv_flagged),
            esd_purified_um = esd_of(log10_percell_bv_purified))]
size[, `:=`(contrast_raw_log10 = log10_percell_bv_raw[depth == 3] - log10_percell_bv_raw[depth == 16],
            contrast_purified_log10 = log10_percell_bv_purified[depth == 3] - log10_percell_bv_purified[depth == 16])]
fwrite(size, file.path(cfg$work_dir, "23e_class_purity_size.csv"))
cat("\n-- per-cell size of the flagged frames and the purified size contrast (23e_class_purity_size.csv) --\n")
print(size[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

set.seed(cfg$seed)
days <- unique(rbind(v[, .(depth, date)], day_mean[, .(depth, date)]))
bt <- replicate(cfg$n_boot, {
  dr <- days[, .(date = sample(date, .N, replace = TRUE)), by = depth]
  vv <- merge(dr, v, by = c("depth", "date"), allow.cartesian = TRUE)
  dd <- merge(dr, day_mean, by = c("depth", "date"), allow.cartesian = TRUE)
  post_strat(vv, dd)
})
ci <- apply(bt, 1, quantile, c(.025, .975), na.rm = TRUE)
out <- data.table(quantity = names(est), estimate = round(est, 3), lo95 = round(ci[1, ], 3), hi95 = round(ci[2, ], 3),
                  note = c("post-stratified mean over-count at 3 m", "post-stratified mean over-count at 16 m",
                           "3 m minus 16 m over-count (the part that can move the contrast)",
                           "post-stratified mean MANUAL count at 3 m", "post-stratified mean MANUAL count at 16 m",
                           "surface-minus-deep contrast from the manual counts alone (post-stratified)",
                           "the detector's contrast on the same validated chains (post-stratified)",
                           "surface-minus-deep Chaetoceros contrast of ROI-pooled detector counts (all Chaetoceros classes)",
                           "the same contrast after subtracting the differential bias"))
fwrite(out, file.path(cfg$work_dir, "23e_differential_bias.csv"))
cat("unrounded:\n"); print(signif(rbind(estimate = est, lo95 = ci[1, ], hi95 = ci[2, ]), 6))
cat(sprintf("\n-- post-stratified bias (day-clustered bootstrap, %d resamples) --\n", cfg$n_boot))
print(out)
cat(sprintf("\noverall: n = %d chains on %d days; exact %.1f%%, within 1 cell %.1f%%, MAE %.2f, bias %+.2f\n",
            nrow(v), uniqueN(v$date), 100 * mean(v$err == 0), 100 * mean(abs(v$err) <= 1), mean(abs(v$err)), mean(v$err)))
