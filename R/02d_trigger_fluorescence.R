# =============================================================================
# Fluorescence-trigger check (Note S1): PMT B fluorescence per cell of single-cell and
# multi-cell chain-former ROIs by depth and month, and the share of ROIs whose PMT B
# peak lies within 20% above the sample's trigger threshold.
#
# Reads:  raw .adc files under data/raw_bins/2016/, data/roi_level/roi_2016.txt.gz,
#         outputs/tangesund/02_sample_meta.csv, 02c_ifcb_settings.csv
# Writes: outputs/tangesund/02d_fluorescence_by_depth.csv, 02d_near_threshold.csv,
#         02d_join_check.csv
# Deps:   data.table
# Run:    Rscript R/02d_trigger_fluorescence.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))
cfg <- list(raw_root = "data/raw_bins/2016", out_dir = "outputs/tangesund",
            roi_file = "data/roi_level/roi_2016.txt.gz", meta_csv = "outputs/tangesund/02_sample_meta.csv",
            settings_csv = "outputs/tangesund/02c_ifcb_settings.csv",
            junk = c("Debris", "Beads", "Air_bubbles", "Zooplankton"),
            depths = c(3, 6, 8, 11, 13, 16), near = 1.2)

meta <- fread(cfg$meta_csv)[qc_pass == TRUE & depth %in% cfg$depths, .(sample, depth, date = as.Date(date))]
meta[, period := factor(month.abb[as.integer(format(date, "%m"))], levels = c("Aug", "Sep", "Oct"))]
sett <- fread(cfg$settings_csv)[, .(sample, pmtB_threshold = as.numeric(pmtB_threshold), trigger_selection)]
meta <- merge(meta, sett, by = "sample")
cat(sprintf("QC-pass samples at the six depths with settings: %d; trigger selection values: %s; PMT B thresholds: %s\n",
            nrow(meta), paste(unique(meta$trigger_selection), collapse = ","),
            paste(sort(unique(meta$pmtB_threshold)), collapse = ", ")))

# ---- per-ROI PMT signals from the ADC files ------------------------------------
# ADC columns: 1 trigger, 3-4 PMT A-B integrated, 7-8 PMT A-B peak, 16-17 ROI width and
# height; ROI number = ADC row
adc_files <- list.files(cfg$raw_root, pattern = "^D2016.*_IFCB110\\.adc$", recursive = TRUE, full.names = TRUE)
adc_files <- adc_files[sub("\\.adc$", "", basename(adc_files)) %in% meta$sample]
read_adc <- function(path) {
  a <- fread(path, header = FALSE, select = c(1, 3, 4, 7, 8, 16, 17), showProgress = FALSE)
  setnames(a, c("trigger", "pmtA_int", "pmtB_int", "pmtA_peak", "pmtB_peak", "roi_w", "roi_h"))
  a[, `:=`(sample = sub("\\.adc$", "", basename(path)), roi_number = seq_len(.N))]
  a[roi_w > 0 & roi_h > 0]
}
adc <- rbindlist(lapply(adc_files, read_adc))
cat(sprintf("ADC rows with a ROI: %d over %d samples\n", nrow(adc), uniqueN(adc$sample)))

# ---- chain-counted ROIs (same rule as 20_within_taxon_chain.R) -----------------
d <- fread(cfg$roi_file, select = c("sample", "roi_number", "class_name", "biovolume_um3", "numBlobs",
                                    "cell_count", "cell_count_resolved"), showProgress = FALSE)
d <- d[numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 & !class_name %in% cfg$junk]
d <- d[!is.na(cell_count) & cell_count != -1]
d[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1, round(cell_count_resolved), 1L)]
d <- merge(d, meta, by = "sample")
n_roi <- nrow(d)
d <- merge(d, adc[, .(sample, roi_number, pmtA_int, pmtB_int, pmtA_peak, pmtB_peak)], by = c("sample", "roi_number"))
d[, `:=`(peak_per_cell = pmtB_peak / nchain, int_per_cell = pmtB_int / nchain,
         int_per_biovolume = 1e3 * pmtB_int / biovolume_um3,      # V s per 1000 um3
         above_threshold = pmtB_peak >= pmtB_threshold,
         near_threshold = pmtB_peak < cfg$near * pmtB_threshold,
         group = fifelse(nchain == 1, "single cell", "multi-cell"))]
jc <- data.table(chain_counted_roi = n_roi, joined_to_adc = nrow(d), pct_joined = 100 * nrow(d) / n_roi,
                 pct_pmtB_peak_at_or_above_threshold = 100 * mean(d$above_threshold),
                 pct_pmtA_peak_at_or_above_pmtB_threshold = 100 * mean(d$pmtA_peak >= d$pmtB_threshold))
fwrite(jc, file.path(cfg$out_dir, "02d_join_check.csv"))
cat(sprintf("joined %d of %d chain-counted ROIs (%.1f%%); PMT B peak at or above the sample threshold in %.1f%% (PMT A: %.1f%%)\n",
            jc$joined_to_adc, jc$chain_counted_roi, jc$pct_joined, jc$pct_pmtB_peak_at_or_above_threshold,
            jc$pct_pmtA_peak_at_or_above_pmtB_threshold))

# ---- (1) fluorescence per cell by depth (and month) -----------------------------
qs <- function(x) as.list(setNames(quantile(x, c(.1, .25, .5, .75, .9)), c("q10", "q25", "q50", "q75", "q90")))
summ <- function(dd, by) rbindlist(list(
  dd[, c(list(quantity = "PMT B peak per cell (V)", n = .N, n_days = uniqueN(date)), qs(peak_per_cell)), by = by],
  dd[, c(list(quantity = "PMT B integrated per cell (V s)", n = .N, n_days = uniqueN(date)), qs(int_per_cell)), by = by],
  dd[, c(list(quantity = "PMT B integrated per 1000 um3 biovolume", n = .N, n_days = uniqueN(date)), qs(int_per_biovolume)), by = by],
  dd[, c(list(quantity = "PMT B peak of the ROI (V)", n = .N, n_days = uniqueN(date)), qs(pmtB_peak)), by = by]))
fl <- rbind(summ(d, c("depth", "group"))[, period := "all"],
            summ(d, c("period", "depth", "group")), fill = TRUE)
fl <- rbind(fl, summ(d, "depth")[, `:=`(group = "all chain-counted", period = "all")], fill = TRUE)
setcolorder(fl, c("period", "depth", "group", "quantity"))
setorder(fl, period, quantity, group, depth)
fwrite(fl, file.path(cfg$out_dir, "02d_fluorescence_by_depth.csv"))

# ---- (2) near-threshold fraction ----------------------------------------------
nt <- rbind(d[, .(n = .N, n_days = uniqueN(date), pct_near_threshold = 100 * mean(near_threshold),
                  pct_below_threshold = 100 * mean(!above_threshold)), by = .(depth, group)][, period := "all"],
            d[, .(n = .N, n_days = uniqueN(date), pct_near_threshold = 100 * mean(near_threshold),
                  pct_below_threshold = 100 * mean(!above_threshold)), by = .(period, depth, group)],
            d[, .(n = .N, n_days = uniqueN(date), pct_near_threshold = 100 * mean(near_threshold),
                  pct_below_threshold = 100 * mean(!above_threshold)), by = .(depth)][, `:=`(group = "all chain-counted", period = "all")],
            fill = TRUE)
nt[, near_factor := cfg$near]
setcolorder(nt, c("period", "depth", "group"))
setorder(nt, period, group, depth)
fwrite(nt, file.path(cfg$out_dir, "02d_near_threshold.csv"))

cat("\n-- median PMT B peak per cell (V) by depth: single cells | multi-cell ROIs --\n")
print(dcast(fl[period == "all" & quantity == "PMT B peak per cell (V)" & group != "all chain-counted",
               .(depth, group, q50 = round(q50, 4))], depth ~ group, value.var = "q50"))
cat("\n-- median PMT B integrated per 1000 um3 biovolume by depth (light-acclimation index) --\n")
print(dcast(fl[period == "all" & quantity == "PMT B integrated per 1000 um3 biovolume" & group != "all chain-counted",
               .(depth, group, q50 = round(q50, 4))], depth ~ group, value.var = "q50"))
cat(sprintf("\n-- ROIs with PMT B peak below %.1f x the sample threshold (%%) --\n", cfg$near))
print(dcast(nt[period == "all", .(depth, group, pct = round(pct_near_threshold, 2))], depth ~ group, value.var = "pct"))
s3 <- nt[period == "all" & group == "single cell" & depth == 3, pct_near_threshold]
s16 <- nt[period == "all" & group == "single cell" & depth == 16, pct_near_threshold]
cat(sprintf("\nsingle-cell chain-former ROIs within 20%% of the trigger: %.2f%% at 3 m, %.2f%% at 16 m\n", s3, s16))
cat(sprintf("tables -> %s/02d_fluorescence_by_depth.csv, 02d_near_threshold.csv, 02d_join_check.csv\n", cfg$out_dir))
