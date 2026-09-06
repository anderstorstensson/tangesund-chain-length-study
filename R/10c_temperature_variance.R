# =============================================================================
# Between-day, within-day vertical and within-day residual shares of the matched
# in-situ temperature variance, per month and overall.
#
# Reads:  outputs/tangesund/06_env_matched.csv, 03_chain_metrics.csv, 10_aic_comparison.csv.
# Writes: outputs/tangesund/10c_temperature_variance.csv.
# Deps:   data.table.
# Run:    Rscript R/10c_temperature_variance.R   (after 10)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
cfg <- list(out_dir = "outputs/tangesund")
env <- fread(file.path(cfg$out_dir, "06_env_matched.csv"))
chn <- fread(file.path(cfg$out_dir, "03_chain_metrics.csv"))[, .(sample, qc_pass, reliable_chain)]
ana <- merge(env, chn, by = "sample")[qc_pass == TRUE & reliable_chain == TRUE & is.finite(insitu_temp) & is.finite(mld)]
ana[, `:=`(date = as.Date(date), month = factor(month.abb[as.integer(format(as.Date(date), "%m"))],
                                                levels = c("Aug", "Sep", "Oct")))]
decomp <- function(x) {
  x[, tday := mean(insitu_temp), by = date]
  x[, tdd  := mean(insitu_temp), by = .(date, depth)]
  tot <- var(x$insitu_temp)
  data.table(n = nrow(x), n_days = uniqueN(x$date), sd_total = sqrt(tot),
             pct_between_day = 100 * var(x$tday) / tot,
             pct_within_day_vertical = 100 * var(x$tdd - x$tday) / tot,
             pct_within_day_residual = 100 * var(x$insitu_temp - x$tdd) / tot,
             within_day_sd = sqrt(var(x$insitu_temp - x$tday)),
             median_3m_minus_16m = {
               w <- dcast(x[depth %in% c(3, 16), .(t = mean(insitu_temp)), by = .(date, depth)], date ~ depth, value.var = "t")
               if (all(c("3", "16") %in% names(w))) median(w$`3` - w$`16`, na.rm = TRUE) else NA_real_ })
}
out <- rbind(ana[, decomp(copy(.SD)), by = month], cbind(month = "All", decomp(copy(ana))))
fwrite(out, file.path(cfg$out_dir, "10c_temperature_variance.csv"))
cat("-- in-situ temperature variance decomposition on the mechanism-model samples (10c_temperature_variance.csv) --\n")
print(out[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
aic <- fread(file.path(cfg$out_dir, "10_aic_comparison.csv"))
a <- setNames(aic$aic, aic$model)
cat(sprintf("\nAIC: depth only %.1f | depth + temperature %.1f | salinity + temperature %.1f | full hydrography %.1f\n",
            a["depth"], a["depth + temperature"], a["salinity + temperature"], a["full hydrography"]))
cat(sprintf("  depth-only -> depth+temperature: %.0f points | depth+temperature -> salinity+temperature: %.0f points\n",
            a["depth"] - a["depth + temperature"], a["depth + temperature"] - a["salinity + temperature"]))
