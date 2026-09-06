# =============================================================================
# Samples and sampling days per month in each salinity bin of Fig. 6 (Note S2).
#
# Reads:  outputs/tangesund/06_env_matched.csv, 03_chain_metrics.csv.
# Writes: outputs/tangesund/10d_salinity_bin_months.csv.
# Deps:   data.table.
# Run:    Rscript R/10d_salinity_bin_months.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))
cfg <- list(out_dir = "outputs/tangesund")
env <- fread(file.path(cfg$out_dir, "06_env_matched.csv"))
chn <- fread(file.path(cfg$out_dir, "03_chain_metrics.csv"))[, .(sample, qc_pass, reliable_chain, mean_chain)]
ana <- merge(env, chn, by = "sample")[qc_pass == TRUE & reliable_chain == TRUE & is.finite(insitu_temp) & is.finite(mld)]
ana[, month := factor(month.abb[as.integer(format(as.Date(date), "%m"))], levels = c("Aug", "Sep", "Oct"))]
cat(sprintf("mechanism sample set reproduced: %d samples\n", nrow(ana)))
# same binning schemes as sal_slopes() in R/10_env_mechanism.R
bins <- function(dat, scheme) {
  dat <- copy(dat)
  if (scheme == "Fixed salinity boundaries") {
    dat[, sbin := cut(insitu_sal, c(-Inf, 26, 28, 29, 30, 31, Inf), labels = c("≤ 26","26–28","28–29","29–30","30–31","> 31"))]
  } else {
    q <- unique(quantile(dat$insitu_sal, seq(0, 1, length = 7)))
    dat[, sbin := cut(insitu_sal, q, include.lowest = TRUE, labels = sprintf("%.1f–%.1f", head(q, -1), tail(q, -1)))]
  }
  dat[, scheme := scheme][]
}
b <- rbind(bins(ana, "Fixed salinity boundaries"), bins(ana, "Equal-n salinity bins"))
out <- b[, .(n = .N, n_days = uniqueN(date),
             n_Aug = sum(month == "Aug"), n_Sep = sum(month == "Sep"), n_Oct = sum(month == "Oct"),
             pct_Aug = round(100 * mean(month == "Aug"), 1), pct_Sep = round(100 * mean(month == "Sep"), 1),
             pct_Oct = round(100 * mean(month == "Oct"), 1),
             days_Aug = uniqueN(date[month == "Aug"]), days_Sep = uniqueN(date[month == "Sep"]),
             days_Oct = uniqueN(date[month == "Oct"]),
             n_depths = uniqueN(depth), sal_lo = min(insitu_sal)), by = .(scheme, bin = sbin)]
setorder(out, scheme, -sal_lo)
fwrite(out, file.path(cfg$out_dir, "10d_salinity_bin_months.csv"))
cat("-- month composition of the salinity bins (10d_salinity_bin_months.csv) --\n")
print(out[, .(scheme, bin, n, n_days, pct_Aug, pct_Sep, pct_Oct, days_Aug, days_Sep, days_Oct)])
