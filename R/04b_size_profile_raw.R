# =============================================================================
# Per-depth size-spectrum exponent b, size diversity and median ESD: raw mean with a
# day-clustered bootstrap SE, raw median, mean of daily means, and GAMM day-adjusted value.
#
# Reads:  outputs/tangesund/04_size_spectra.csv
# Writes: outputs/tangesund/04b_size_profile_raw.csv (read by 04 for its panel a)
# Deps:   data.table, mgcv
# Run:    Rscript R/04b_size_profile_raw.R   (after 04; then rerun 04)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(mgcv) })
cfg <- list(csv = "outputs/tangesund/04_size_spectra.csv", out_dir = "outputs/tangesund",
            depths = c(3, 6, 8, 11, 13, 16), B = 2000)
ana <- fread(cfg$csv)[qc_pass == TRUE & reliable == TRUE]
ana[, `:=`(day = factor(date))]

# day-cluster bootstrap SE of the sample-weighted mean (as in 03 and 05)
boot_se_day <- function(value, day, B = cfg$B, seed = 1) {
  days <- unique(day); nd <- length(days); idx <- split(seq_along(value), day)
  set.seed(seed)
  sd(vapply(seq_len(B), function(b) mean(value[unlist(idx[as.character(sample(days, nd, TRUE))], use.names = FALSE)]),
            numeric(1)))
}
vars <- c(median_esd = "Median community ESD (um)", size_div = "Size diversity", b = "Size-spectrum exponent b")
out <- rbindlist(lapply(names(vars), function(v) {
  # same GAMM as 04
  m <- bam(reformulate(c("s(depth,k=5)","s(yday,k=8)","ti(depth,yday,k=c(5,6))","s(day,bs=\"re\")"), v),
           data = ana, method = "fREML")
  g <- data.frame(depth = cfg$depths, yday = median(ana$yday), day = ana$day[1])
  adj <- as.numeric(predict(m, g, exclude = "s(day)"))
  daily <- ana[, .(v = mean(get(v))), by = .(depth, day)][, .(mean_of_daily_means = mean(v)), by = depth]
  r <- ana[, .(n = .N, n_days = uniqueN(day), raw_mean = mean(get(v)), raw_median = median(get(v)),
               se_day = boot_se_day(get(v), as.character(day))), by = depth][order(depth)]
  r <- merge(r, daily, by = "depth")
  r[, `:=`(variable = v, label = vars[v], gamm_day_adjusted = adj[match(depth, cfg$depths)])]
  r[]
}))
setcolorder(out, c("variable", "label", "depth"))
fwrite(out, file.path(cfg$out_dir, "04b_size_profile_raw.csv"))
cat("-- community size profiles: raw mean (+ day-clustered SE), raw median, mean of daily means, GAMM day-adjusted --\n")
print(out[, .(variable, depth, n, raw_mean = round(raw_mean, 3), se_day = round(se_day, 3),
              raw_median = round(raw_median, 3), daily = round(mean_of_daily_means, 3),
              adj = round(gamm_day_adjusted, 3))])
for (v in names(vars)) {
  r <- out[variable == v]
  cat(sprintf("%s: raw mean 3 m %.3f vs 16 m %.3f | extreme at %d m (%.3f)\n", v,
              r[depth == 3, raw_mean], r[depth == 16, raw_mean],
              r[if (v == "size_div") which.max(raw_mean) else which.min(raw_mean), depth],
              r[if (v == "size_div") which.max(raw_mean) else which.min(raw_mean), raw_mean]))
}
cat(sprintf("table -> %s/04b_size_profile_raw.csv\n", cfg$out_dir))
