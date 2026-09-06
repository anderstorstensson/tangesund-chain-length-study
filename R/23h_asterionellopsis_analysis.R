# =============================================================================
# Asterionellopsis glacialis validation set, stage 3 of 3: depth-differential count bias
# from the manual counts (manual count = "cell" boxes per task), post-stratified to each
# depth's Feret composition, with day-clustered bootstrap intervals (Tables S4, S5); class
# purity and its effect on the chain-length contrast.
# Frames flagged "not Asterionellopsis glacialis" or "cells not resolvable" are excluded.
#
# Reads:  the Label Studio export (default data/asterionellopsis_validation/label_studio_export_project25.json),
#         the bundled data/asterionellopsis_validation/ key and population bins (copied from 23f's
#         outputs if absent), data/roi_level/roi_2016.txt.gz, outputs/tangesund/02_sample_meta.csv
# Writes: data/asterionellopsis_validation/ asterionellopsis_manual_counts.csv, asterionellopsis_cell_boxes.csv;
#         outputs/tangesund/asterionellopsis_validation/ 23h_bias_by_stratum.csv, 23h_feret_quartiles.csv,
#         23h_class_purity.csv, 23h_differential_bias.csv
# Deps:   data.table, jsonlite
# Run:    Rscript R/23h_asterionellopsis_analysis.R [export.json]   (no Label Studio server needed)
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

cfg <- list(
  work_dir = "outputs/tangesund/asterionellopsis_validation",
  data_dir = "data/asterionellopsis_validation",            # tracked once counted: the hand-made counts cannot be regenerated
  roi_file = "data/roi_level/roi_2016.txt.gz",
  meta_csv = "outputs/tangesund/02_sample_meta.csv",
  micron   = 1 / 2.77,
  feret_breaks = c(0, 18, 34, 98, Inf),       # um, as R/23f_asterionellopsis_select.R
  feret_labels = c("<=18", "18-34", "34-98", ">98"),
  flag_not_taxon = "not Asterionellopsis glacialis",
  flag_exclude   = "not Asterionellopsis glacialis|not resolvable",
  n_boot   = 2000, seed = 1
)
export <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(export)) export <- file.path(cfg$data_dir, "label_studio_export_project25.json")
if (!file.exists(export)) stop("usage: Rscript R/23h_asterionellopsis_analysis.R [label-studio export.json]")
for (p in c(cfg$data_dir, cfg$work_dir)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# ---- manual counts and boxes from the export ------------------------------------
tasks <- fromJSON(export, simplifyVector = FALSE)
task_id <- function(t) sub("\\.png$", "", basename(sub(".*[?&]d=", "", utils::URLdecode(t$data$image))))
latest <- function(t) {
  ann <- Filter(function(a) !isTRUE(a$was_cancelled), t$annotations)
  if (length(ann)) ann[[length(ann)]] else NULL                   # latest completed annotation
}
manual <- rbindlist(lapply(tasks, function(t) {
  if (is.null(t$data$image)) return(NULL)
  a <- latest(t)
  if (is.null(a)) return(data.table(anon_id = task_id(t), manual = NA_integer_, flags = ""))
  n_box <- sum(vapply(a$result, function(r) identical(r$type, "rectanglelabels"), logical(1)))
  fl <- unlist(lapply(a$result, function(r) if (identical(r$type, "choices")) r$value$choices))
  data.table(anon_id = task_id(t), manual = as.integer(n_box), flags = paste(fl, collapse = "; "))
}))
boxes <- rbindlist(lapply(tasks, function(t) {
  a <- if (is.null(t$data$image)) NULL else latest(t)
  rb <- if (is.null(a)) list() else Filter(function(r) identical(r$type, "rectanglelabels"), a$result)
  if (!length(rb)) return(NULL)
  rbindlist(lapply(rb, function(r) data.table(
    anon_id = task_id(t), image_width_px = r$original_width, image_height_px = r$original_height,
    x_pct = r$value$x, y_pct = r$value$y, width_pct = r$value$width, height_pct = r$value$height)))
}))
key_file <- file.path(cfg$data_dir, "asterionellopsis_selection_key.csv")
if (!file.exists(key_file)) invisible(file.copy(file.path(cfg$work_dir, "23f_selection_key.csv"), key_file))
key <- fread(key_file)
key[, date := as.Date(date)]
v <- merge(key, manual, by = "anon_id", all.x = TRUE)
setnames(v, c("det_count", "manual"), c("detector_count", "manual_count"))
fwrite(v[order(anon_id)], file.path(cfg$data_dir, "asterionellopsis_manual_counts.csv"))      # one row per blinded colony
fwrite(boxes[order(anon_id)], file.path(cfg$data_dir, "asterionellopsis_cell_boxes.csv"))       # one row per manually drawn cell
setnames(v, c("detector_count", "manual_count"), c("det_count", "manual"))
cat(sprintf("tasks in export: %d; annotated: %d of %d drawn\n", nrow(manual), sum(!is.na(v$manual)), nrow(v)))
v[, feret_bin := factor(feret_bin, levels = cfg$feret_labels)]
v_all <- copy(v)                                                   # every drawn frame, flags included
excl <- v[grepl(cfg$flag_exclude, flags)]
cat(sprintf("excluded by flag: %d (%s at 3 / 16 m)\n", nrow(excl),
            paste(excl[, .N, keyby = factor(depth, c(3, 16))]$N, collapse = " / ")))
v <- v[!is.na(manual) & !grepl(cfg$flag_exclude, flags)]
v[, err := det_count - manual]

# ---- (a) per stratum ------------------------------------------------------------
drawn <- v_all[, .(n_drawn = .N, n_unannotated = sum(is.na(manual)),
                   n_not_taxon = sum(grepl(cfg$flag_not_taxon, flags)),
                   n_unresolvable = sum(grepl("not resolvable", flags))), by = .(depth, feret_bin)]
cell <- v[, .(n = .N, days = uniqueN(date), bias = mean(err), mae = mean(abs(err)),
              exact_pct = 100 * mean(err == 0), within1_pct = 100 * mean(abs(err) <= 1),
              mean_manual = mean(manual), mean_det = mean(det_count)), by = .(depth, feret_bin)]
cell <- merge(drawn, cell, by = c("depth", "feret_bin"), all.x = TRUE)
setorder(cell, depth, feret_bin)
fwrite(cell, file.path(cfg$work_dir, "23h_bias_by_stratum.csv"))
cat("\n-- count error (detector - manual) by depth x maximum Feret diameter --\n")
print(cell[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

# ---- (b)/(c) post-stratified bias and the differential -------------------------
pop_file <- file.path(cfg$data_dir, "asterionellopsis_population_bins.csv")
if (!file.exists(pop_file)) invisible(file.copy(file.path(cfg$work_dir, "23f_population_bins.csv"), pop_file))
pop <- fread(pop_file)
pop[, feret_bin := factor(feret_bin, levels = cfg$feret_labels)]
meta <- fread(cfg$meta_csv)[qc_pass == TRUE, .(sample, depth, date = as.Date(date))]
d <- fread(cfg$roi_file, select = c("sample", "roi_number", "class_name", "numBlobs", "biovolume_um3", "cell_count",
                                    "cell_count_resolved", "maxFeretDiameter"), showProgress = FALSE)
d <- d[grepl("^Asterionellopsis_glacialis", class_name) & !is.na(cell_count) & cell_count != -1 &
       numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0]
d[, det_count := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1,
                         as.integer(round(cell_count_resolved)), 1L)]
d <- merge(d, meta, by = "sample")
# the pooled maximum-Feret quartiles the draw's bins were placed at (Supplementary Methods)
d[, feret_um := maxFeretDiameter * cfg$micron]
fq_of <- function(x, lab) data.table(population = lab, q25_um = round(quantile(x, .25), 2), q50_um = round(quantile(x, .5), 2),
                                     q75_um = round(quantile(x, .75), 2), n_roi = length(x))
fq <- rbind(fq_of(d$feret_um, "all six depths"), fq_of(d[depth %in% c(3, 16), feret_um], "3 m and 16 m"))
fwrite(fq, file.path(cfg$work_dir, "23h_feret_quartiles.csv"))
cat("\n-- A. glacialis maximum Feret diameter quartiles (um) --\n"); print(fq)
d <- d[depth %in% c(3, 16)]
day_mean <- d[, .(m = mean(det_count), n = .N), by = .(depth, date)]   # day-level means + weights (clustered)

post_strat <- function(vv, dd) {
  b <- vv[, .(bias = mean(err), man = mean(manual), det = mean(det_count)), by = .(depth, feret_bin)]
  b <- merge(pop[depth %in% c(3, 16), .(depth, feret_bin, w)], b, by = c("depth", "feret_bin"), all.x = TRUE)
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

# ---- (d) class purity on the blinded draw ---------------------------------------
# The frames judged not A. glacialis, weighted by each Feret bin's share of the
# imaged A. glacialis at that depth, give the population-weighted misclassification
# rate of the class on this deployment. The purified contrast removes those
# frames, at their detector counts, from the ROI-pooled depth means (as 23e).
pur <- v_all[, .(rate = mean(grepl(cfg$flag_not_taxon, flags)),
                 det_nt = mean(det_count[grepl(cfg$flag_not_taxon, flags)])), by = .(depth, feret_bin)]
pur <- merge(pop[depth %in% c(3, 16), .(depth, feret_bin, w)], pur, by = c("depth", "feret_bin"))
pur[is.na(det_nt), det_nt := 0]
raw_pooled <- day_mean[, .(pooled_mean = weighted.mean(m, n)), by = depth]
purity <- pur[, .(not_taxon_pct = 100 * sum(w * rate),
                  mean_det_not_taxon = sum(w * rate * det_nt) / sum(w * rate)), by = depth]
purity <- merge(purity, raw_pooled, by = "depth")
purity[, purified_mean := (pooled_mean - (not_taxon_pct / 100) * mean_det_not_taxon) / (1 - not_taxon_pct / 100)]
purity[, `:=`(contrast_raw = pooled_mean[depth == 3] - pooled_mean[depth == 16],
              contrast_purified = purified_mean[depth == 3] - purified_mean[depth == 16])]
# the same at the equal allocation of the draw (40 frames per Feret bin)
pur_eq <- v_all[, .(not_taxon_pct_equal_alloc = 100 * mean(grepl(cfg$flag_not_taxon, flags)),
                    mean_det_not_taxon_equal_alloc = mean(det_count[grepl(cfg$flag_not_taxon, flags)]),
                    n_flagged_not_taxon = sum(grepl(cfg$flag_not_taxon, flags)),
                    n_unresolvable = sum(grepl("not resolvable", flags)), n_drawn = .N), by = depth]
purity <- merge(purity, pur_eq, by = "depth")
fwrite(purity, file.path(cfg$work_dir, "23h_class_purity.csv"))
cat("\n-- class purity on the blinded draw (23h_class_purity.csv) --\n")
print(purity[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

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
                           "the detector's contrast on the same validated colonies (post-stratified)",
                           "surface-minus-deep A. glacialis contrast of ROI-pooled detector counts",
                           "the same contrast after subtracting the differential bias"))
fwrite(out, file.path(cfg$work_dir, "23h_differential_bias.csv"))
cat("unrounded:\n"); print(signif(rbind(estimate = est, lo95 = ci[1, ], hi95 = ci[2, ]), 6))
cat(sprintf("\n-- post-stratified bias (day-clustered bootstrap, %d resamples) --\n", cfg$n_boot))
print(out)
cat(sprintf("\noverall: n = %d colonies on %d days; exact %.1f%%, within 1 cell %.1f%%, MAE %.2f, bias %+.2f\n",
            nrow(v), uniqueN(v$date), 100 * mean(v$err == 0), 100 * mean(abs(v$err) <= 1), mean(abs(v$err)), mean(v$err)))
