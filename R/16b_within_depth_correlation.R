# =============================================================================
# Sample-level correlation of chain length and per-cell ESD across samples, within
# depth, within depth and day, and per depth, with day-cluster bootstrap intervals.
#
# Reads:  outputs/tangesund/16_percell_metrics.csv.
# Writes: outputs/tangesund/16b_within_depth_correlation.csv (read by 16 for Fig. 4b).
# Deps:   data.table.
# Run:    Rscript R/16b_within_depth_correlation.R   (after 16; then rerun 16 or let 28 do it)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
cfg <- list(csv = "outputs/tangesund/16_percell_metrics.csv", out_dir = "outputs/tangesund", B = 2000)
s <- fread(cfg$csv)
s[, `:=`(day = factor(date), depth_f = factor(depth))]
resid_on <- function(y, f) resid(lm(y ~ f))
s[, `:=`(chain_d = resid_on(mean_chain, depth_f), esd_d = resid_on(percell_esd, depth_f))]
s[, `:=`(chain_dd = resid(lm(mean_chain ~ depth_f + day)), esd_dd = resid(lm(percell_esd ~ depth_f + day)))]

# day-clustered bootstrap: resample whole sampling days, recompute r each time
boot_r <- function(x, y, day, B = cfg$B, seed = 1) {
  days <- unique(day); nd <- length(days); idx <- split(seq_along(x), day)
  set.seed(seed)
  r <- vapply(seq_len(B), function(b) {
    i <- unlist(idx[as.character(sample(days, nd, TRUE))], use.names = FALSE)
    cor(x[i], y[i]) }, numeric(1))
  c(lo = unname(quantile(r, .025, na.rm = TRUE)), hi = unname(quantile(r, .975, na.rm = TRUE)),
    p = 2 * min(mean(r <= 0, na.rm = TRUE), mean(r >= 0, na.rm = TRUE)))
}
row <- function(label, x, y, day) {
  b <- boot_r(x, y, as.character(day))
  data.table(correlation = label, n = length(x), n_days = uniqueN(day), r = cor(x, y),
             ci_lo = b["lo"], ci_hi = b["hi"], p_day_bootstrap = b["p"])
}
out <- rbind(
  row("across samples (as in Fig. 4b)", s$mean_chain, s$percell_esd, s$day),
  row("within depth (depth means removed)", s$chain_d, s$esd_d, s$day),
  row("within depth and day (depth and day means removed)", s$chain_dd, s$esd_dd, s$day),
  rbindlist(lapply(sort(unique(s$depth)), function(dp) {
    x <- s[depth == dp]; row(sprintf("within %d m", dp), x$mean_chain, x$percell_esd, x$day) })))
out[, p_day_bootstrap := pmax(p_day_bootstrap, 1 / cfg$B)]   # a bootstrap p cannot be below 1/B
fwrite(out, file.path(cfg$out_dir, "16b_within_depth_correlation.csv"))
cat("-- chain length x per-cell ESD correlations (16b_within_depth_correlation.csv) --\n")
print(out[, .(correlation, n, n_days, r = round(r, 3), CI = sprintf("[%.3f, %.3f]", ci_lo, ci_hi), p = signif(p_day_bootstrap, 2))])
