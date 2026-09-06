# =============================================================================
# Sensitivity to the chain-trait sample threshold (>= 20 chain ROIs in 03): retention
# by depth and period (Table 1), the depth profile and the 3 m - 16 m GAMM contrast
# at thresholds 0-50, and an unthresholded ROI-weighted refit.
#
# Reads:  outputs/tangesund/02_sample_meta.csv, 03_chain_metrics.csv
# Writes: outputs/tangesund/24_threshold_retention.csv, 24_threshold_retention_by_period.csv,
#         24_threshold_profile_sweep.csv, 24_threshold_contrast_sweep.csv
# Deps:   data.table, mgcv
# Run:    Rscript R/24_trait_threshold_sensitivity.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(mgcv) })

cfg <- list(out_dir = "outputs/tangesund", meta = "outputs/tangesund/02_sample_meta.csv",
            thresholds = c(0, 5, 10, 20, 50), depths = c(3, 6, 8, 11, 13, 16))

meta <- fread(cfg$meta)[, .(sample, depth, date, yday, qc_pass)]
chn  <- fread(file.path(cfg$out_dir, "03_chain_metrics.csv"))[
  , .(sample, n_chain_roi, mean_chain, cells_chain)]
D <- merge(meta, chn, by = "sample", all.x = TRUE)[qc_pass == TRUE & is.finite(mean_chain)]
D[, `:=`(period = factor(month.abb[as.integer(format(as.Date(date), "%m"))],
                         levels = c("Aug", "Sep", "Oct")), day = factor(date))]

# ---- (A) retention by depth, at the reported threshold of 20 chain ROIs -------
ret <- D[, .(samples = .N, pass20 = sum(n_chain_roi >= 20)), by = depth][order(depth)]
ret[, retained_pct := round(100 * pass20 / samples, 1)]
cat("== (A) chain-trait threshold retention by depth ==\n"); print(ret)
fwrite(ret, file.path(cfg$out_dir, "24_threshold_retention.csv"))

# retention by period (quoted in the Results)
ret_p <- D[, .(samples = .N, pass20 = sum(n_chain_roi >= 20)), by = period][order(period)]
ret_p[, retained_pct := round(100 * pass20 / samples, 1)]
cat("== (A2) chain-trait threshold retention by period ==\n"); print(ret_p)
fwrite(ret_p, file.path(cfg$out_dir, "24_threshold_retention_by_period.csv"))

# ---- (B) depth profile by threshold: mean of sample means and pooled ------------
prof <- rbindlist(lapply(cfg$thresholds, function(th) {
  d <- D[n_chain_roi >= th]
  d[, .(threshold = th, n = .N,
        mean_of_means = mean(mean_chain),
        pooled_weighted = sum(cells_chain) / sum(n_chain_roi)), by = depth][order(depth)]
}))
cat("\n== (B) chain-length profile by sample threshold ==\n")
print(dcast(prof, threshold ~ depth, value.var = "mean_of_means")[, lapply(.SD, function(x) round(x, 2))])
cat("\n-- same, abundance-weighted (pooled cells / chains) --\n")
print(dcast(prof, threshold ~ depth, value.var = "pooled_weighted")[, lapply(.SD, function(x) round(x, 2))])

# monotonicity and raw surface-deep difference per threshold
smry <- prof[, .(n = sum(n),
                 surf_mm = mean_of_means[depth == 3], deep_mm = mean_of_means[depth == 16],
                 surf_pw = pooled_weighted[depth == 3], deep_pw = pooled_weighted[depth == 16],
                 monotonic = all(diff(mean_of_means[order(depth)]) < 0)), by = threshold]
smry[, `:=`(diff_mm = surf_mm - deep_mm, diff_pw = surf_pw - deep_pw)]
cat("\n-- surface-deep difference and monotonicity by threshold --\n")
print(smry[, .(threshold, n, surf_mm = round(surf_mm,2), deep_mm = round(deep_mm,2),
               diff_mm = round(diff_mm,2), diff_pw = round(diff_pw,2), monotonic)])
fwrite(prof, file.path(cfg$out_dir, "24_threshold_profile_sweep.csv"))

# ---- (C) formal GAMM surface-deep contrast per period, per threshold ---------
# same specification as contrast_depth() in 13
contrast_at <- function(th) {
  d2 <- D[n_chain_roi >= th]
  m <- bam(mean_chain ~ period + s(depth, by = period, k = 5) + s(yday, k = 8) + s(day, bs = "re"),
           data = d2, method = "fREML")
  rbindlist(lapply(levels(D$period), function(pp) {
    sub <- d2[period == pp]; if (nrow(sub) < 40) return(NULL)
    g <- data.frame(depth = c(3, 16), period = factor(pp, levels = levels(D$period)),
                    yday = median(sub$yday), day = sub$day[1])
    Xp <- predict(m, g, type = "lpmatrix")
    cvec <- Xp[1, ] - Xp[2, ]
    est <- sum(cvec * coef(m)); se <- sqrt(drop(cvec %*% vcov(m) %*% cvec))
    data.table(threshold = th, period = pp, est = est, lo = est - 2 * se, hi = est + 2 * se,
               n = nrow(sub))
  }))
}
contr <- rbindlist(lapply(cfg$thresholds, contrast_at))

# (C2) no threshold, samples weighted by chain ROI count (reported as threshold = -1)
contrast_weighted <- function() {
  d2 <- D
  m <- bam(mean_chain ~ period + s(depth, by = period, k = 5) + s(yday, k = 8) + s(day, bs = "re"),
           data = d2, weights = d2$n_chain_roi, method = "fREML")
  rbindlist(lapply(levels(D$period), function(pp) {
    sub <- d2[period == pp]
    g <- data.frame(depth = c(3, 16), period = factor(pp, levels = levels(D$period)),
                    yday = median(sub$yday), day = sub$day[1])
    Xp <- predict(m, g, type = "lpmatrix")
    cvec <- Xp[1, ] - Xp[2, ]
    est <- sum(cvec * coef(m)); se <- sqrt(drop(cvec %*% vcov(m) %*% cvec))
    data.table(threshold = -1, period = pp, est = est, lo = est - 2 * se, hi = est + 2 * se,
               n = nrow(sub))
  }))
}
contr <- rbind(contr, contrast_weighted())
contr[, sig := ifelse(lo > 0 | hi < 0, "*", "")]
cat("\n(threshold = -1 denotes the unthresholded, chain-ROI-weighted refit)\n")
cat("\n== (C) GAMM surface(3 m)-deep(16 m) chain-length contrast by sample threshold ==\n")
print(contr[, .(threshold, period, est = round(est, 2), lo = round(lo, 2), hi = round(hi, 2), sig, n)])
fwrite(contr, file.path(cfg$out_dir, "24_threshold_contrast_sweep.csv"))
