# =============================================================================
# Triggers, ROIs and ROIs with a segmentable object per sample, 2016 deployment.
#
# Reads:  raw .adc files under data/raw_bins/2016; outputs/tangesund/02_sample_meta.csv.
# Writes: outputs/tangesund/02b_roi_trigger_counts.csv (one row per recorded sample);
#         prints totals for analyzed and all recorded samples.
# Deps:   data.table.
# Run:    Rscript R/02b_roi_trigger_counts.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))

raw_root <- "data/raw_bins/2016"
meta <- fread("outputs/tangesund/02_sample_meta.csv")[, .(sample, n_roi)]

adc_files <- list.files(raw_root, pattern = "^D2016.*_IFCB110\\.adc$", recursive = TRUE, full.names = TRUE)
count_one <- function(path) {
  a <- fread(path, header = FALSE, select = c(1, 16, 17))
  setnames(a, c("trigger", "roi_w", "roi_h"))
  data.table(sample = sub("\\.adc$", "", basename(path)),
             adc_rows = nrow(a),
             triggers = uniqueN(a$trigger),
             rois = sum(a$roi_w > 0 & a$roi_h > 0))
}
cnt <- rbindlist(lapply(adc_files, count_one))
cnt <- merge(cnt, meta, by = "sample", all.x = TRUE)
setnames(cnt, "n_roi", "rois_with_object")
cnt[, analysed := !is.na(rois_with_object)]
setorder(cnt, sample)
fwrite(cnt, "outputs/tangesund/02b_roi_trigger_counts.csv")

tot <- cnt[, .(samples = .N, triggers = sum(triggers), rois = sum(rois),
               rois_with_object = sum(rois_with_object, na.rm = TRUE)), by = analysed]
print(tot)
cat(sprintf("\nAnalysed samples: %d triggers, %d ROIs, %d ROIs with a segmentable object (the rows that enter the analysis).\nAll recorded samples: %d ROIs (the dashboard total).\n",
            tot[analysed == TRUE, triggers], tot[analysed == TRUE, rois], tot[analysed == TRUE, rois_with_object],
            sum(tot$rois)))
