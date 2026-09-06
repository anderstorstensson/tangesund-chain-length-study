# =============================================================================
# Svea chain-length contrast over salinity 25-33 on autumn cruises only (Sep-Nov and
# Aug-Oct), with the cruise-cluster bootstrap and both colony weightings of 09b.
#
# Reads:  data/SHARK_PlanktonImaging_<year>_SMHI/processed_data/data.txt (2022-2026),
#         data/svea/svea_cruise_corrections.csv, outputs/tangesund/09b_svea_chain_by_sample_taxon.rds;
#         falls back to data/svea/svea_samples.csv and svea_chain_by_sample_taxon.csv
#         when the deliveries and the cache are absent. Missing years are skipped.
# Writes: outputs/tangesund/09c_svea_autumn.csv
# Deps:   data.table, mgcv, parallel
# Run:    Rscript R/09c_svea_autumn.R   (after 09b)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(mgcv); library(parallel) })
cfg <- list(data_glob = "data/SHARK_PlanktonImaging_%d_SMHI/processed_data/data.txt", years = 2022:2026,
            cache = "outputs/tangesund/09b_svea_chain_by_sample_taxon.rds", out_dir = "outputs/tangesund",
            deposit_meta = "data/svea/svea_samples.csv", deposit_agg = "data/svea/svea_chain_by_sample_taxon.csv",
            min_chain_roi = 20, qc_exclude = list(c("2024-06-02", "2025-07-14"), c("2022-01-01", "2022-09-30")),
            marine_lo = 25, marine_hi = 33, B = 1000, seed = 42, cores = 4)
# ---- sample metadata from the SHARK deliveries, as in 09b -------------------------
meta <- rbindlist(lapply(cfg$years, function(yr) {
  f <- sprintf(cfg$data_glob, yr); if (!file.exists(f)) return(NULL)
  e <- fread(f, select = c("SMPNO", "CRUISE_NO", "IMAGE_VERIFICATION", "PSAL_FB", "SDATE"), na.strings = c("", "NA"), showProgress = FALSE)
  u <- unique(e[IMAGE_VERIFICATION == "PredictedByMachine", .(sample = SMPNO, cruise_no = as.character(CRUISE_NO),
                                                             salinity = as.numeric(PSAL_FB), date = as.Date(SDATE))])
  u[, year := yr][]
}), use.names = TRUE)
# cruise numbers missing from the 2026 delivery, as in 09b
fix_csv <- "data/svea/svea_cruise_corrections.csv"
if (nrow(meta) > 0 && file.exists(fix_csv)) {
  fx <- fread(fix_csv, colClasses = list(character = c("sample", "cruise_no")))
  meta[fx, on = "sample", cruise_no := i.cruise_no]
}
from_deposit <- nrow(meta) == 0 && file.exists(cfg$deposit_meta)
if (from_deposit) {                                   # bundled QC-pass subset
  meta <- fread(cfg$deposit_meta, colClasses = list(character = c("sample", "cruise_no")))
  meta[, date := as.Date(date)]
}
meta <- meta[is.finite(salinity) & !is.na(cruise_no) & !is.na(date)]
meta[, `:=`(cruise = paste(year, cruise_no, sep = "_"), month = as.integer(format(date, "%m")), qc_pass = TRUE)]
for (p in cfg$qc_exclude) meta[date >= as.Date(p[1]) & date <= as.Date(p[2]), qc_pass := FALSE]
if (!file.exists(cfg$cache) && !file.exists(cfg$deposit_agg))
  stop("missing ", cfg$cache, " and ", cfg$deposit_agg, ": run R/09b_regional_salinity.R first")
agg <- if (file.exists(cfg$cache)) readRDS(cfg$cache) else fread(cfg$deposit_agg)
setorder(agg, sample, taxon)
cm <- agg[, .(chain_equal = sum(chain_len * n_roi) / sum(n_roi), chain_bv = sum(chain_len * biovol) / sum(biovol), n_roi = sum(n_roi)), by = sample]
cm <- merge(cm, meta, by = "sample")[qc_pass == TRUE & n_roi >= cfg$min_chain_roi]
mar <- cm[salinity >= cfg$marine_lo]
cat(sprintf("marine span: %d samples, %d cruises | by month: %s\n", nrow(mar), uniqueN(mar$cruise),
            paste(sprintf("%02d:%d", sort(unique(mar$month)), table(mar$month)), collapse = " ")))

# ---- cruise-cluster bootstrap of the salinity contrast, as in 09b -----------------
RNGkind("L'Ecuyer-CMRG"); set.seed(cfg$seed)
resample_cruises <- function(dat) {
  cr <- unique(dat$cruise); pick <- sample(cr, length(cr), replace = TRUE)
  rbindlist(lapply(seq_along(pick), function(i) { x <- dat[cruise == pick[i]]; x[, cruise := paste0("b", i)]; x }))
}
# no season term: each set is one season
fit_contrast <- function(dd, resp) {
  dd <- copy(dd); dd[, `:=`(y = log(get(resp)), cruise = factor(cruise))]
  f <- bam(y ~ s(salinity, k = 5) + s(cruise, bs = "re"), data = dd, discrete = TRUE, nthreads = 1)
  nd <- data.frame(salinity = c(cfg$marine_lo, cfg$marine_hi), cruise = dd$cruise[1])
  p <- as.numeric(predict(f, nd, type = "terms", terms = "s(salinity)")[, 1]); p[2] - p[1]
}
sets <- list(`September to November` = 9:11, `August to October` = 8:10)
out <- rbindlist(lapply(names(sets), function(s) rbindlist(lapply(c("chain_equal", "chain_bv"), function(resp) {
  dat <- mar[month %in% sets[[s]]]
  obs <- fit_contrast(dat, resp)
  draws <- unlist(mclapply(seq_len(cfg$B), function(b) {
    r <- try(fit_contrast(resample_cruises(dat), resp), silent = TRUE); if (inherits(r, "try-error")) NA_real_ else r
  }, mc.cores = cfg$cores, mc.set.seed = TRUE))
  q <- quantile(draws, c(.025, .975), na.rm = TRUE)
  data.table(months = s, weighting = ifelse(resp == "chain_equal", "colonies weighted equally", "biovolume-weighted"),
             lo_sal = cfg$marine_lo, hi_sal = cfg$marine_hi, n = nrow(dat), cruises = uniqueN(dat$cruise),
             salinity_range = sprintf("%.1f-%.1f", min(dat$salinity), max(dat$salinity)),
             contrast = obs, ci_lo = q[1], ci_hi = q[2], n_draws = sum(is.finite(draws)))
}))))
fwrite(out, file.path(cfg$out_dir, "09c_svea_autumn.csv"))
cat("\n-- autumn-only marine-span contrasts, salinity 25 to 33 (09c_svea_autumn.csv) --\n")
print(out[, .(months, weighting, n, cruises, salinity_range, contrast = round(contrast, 3), CI = sprintf("[%.3f, %.3f]", ci_lo, ci_hi))])
