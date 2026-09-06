# =============================================================================
# IFCB acquisition settings (.hdr) and analyzed volume per sample, 2016 deployment.
#
# Reads:  raw .hdr files under data/raw_bins/2016; outputs/tangesund/02_sample_meta.csv;
#         data/tangesund_samples.csv (analyzed volume per sample, from the SMHI delivery).
# Writes: outputs/tangesund/02c_ifcb_settings.csv (one row per recorded sample);
#         prints runs of identical settings and analyzed-volume quantiles.
# Deps:   data.table.
# Run:    Rscript R/02c_ifcb_settings.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))

raw_root <- "data/raw_bins/2016"
keys <- c(pmtA_threshold = "PMTAtriggerThreshold_DAQ_MCConly",
          pmtB_threshold = "PMTBtriggerThreshold_DAQ_MCConly",
          pmtA_hv = "PMTAhighVoltage", pmtB_hv = "PMTBhighVoltage",
          trigger_selection = "PMTtriggerSelection_DAQ_MCConly",
          syringe_volume_ml = "SyringeSampleVolume",
          binarize_threshold = "binarizeThreshold", minimum_blob_area = "minimumBlobArea",
          software = "softwareVersion")

read_hdr <- function(path) {
  l <- readLines(path, warn = FALSE)
  v <- vapply(keys, function(k) {
    x <- l[startsWith(l, paste0(k, ":"))]
    if (length(x)) trimws(sub("^[^:]*:", "", x[1])) else NA_character_
  }, character(1))
  as.data.table(c(list(sample = sub("\\.hdr$", "", basename(path))), as.list(v)))
}
hdr <- rbindlist(lapply(list.files(raw_root, "^D2016.*_IFCB110\\.hdr$", recursive = TRUE, full.names = TRUE), read_hdr))
hdr[, datetime := as.POSIXct(substr(sample, 2, 16), format = "%Y%m%dT%H%M%S", tz = "UTC")]

meta <- fread("outputs/tangesund/02_sample_meta.csv")[, .(sample, depth)]
vol <- fread("data/tangesund_samples.csv", select = c("sample", "ml_analyzed"))
out <- merge(merge(hdr, meta, by = "sample", all.x = TRUE), vol, by = "sample", all.x = TRUE)
out[, analysed := !is.na(depth)]
setorder(out, datetime)
fwrite(out, "outputs/tangesund/02c_ifcb_settings.csv")

sig <- c("pmtA_threshold", "pmtB_threshold", "pmtA_hv", "pmtB_hv", "trigger_selection", "syringe_volume_ml")
out[, run := rleid(do.call(paste, .SD)), .SDcols = sig]
runs <- out[, c(list(from = format(min(datetime)), to = format(max(datetime)), n = .N, n_analysed = sum(analysed)),
                lapply(.SD, first)), by = run, .SDcols = sig]
print(runs)
cat("\nAnalyzed volume (mL) over the analyzed samples:\n")
print(out[analysed == TRUE, as.list(round(quantile(ml_analyzed, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), 3))])
cat("software:", unique(out$software), "\n")
