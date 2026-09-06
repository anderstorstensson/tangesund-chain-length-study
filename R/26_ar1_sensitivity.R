# =============================================================================
# AR1 sensitivity: the 3 m - 16 m chain-length contrast (day RE, within-day AR1, and
# day x depth aggregate), the hydrography smooths and mixed-layer contrast, and the
# day-night contrast, each naive vs AR1.
#
# Reads:  outputs/tangesund/02_sample_meta.csv, 03_chain_metrics.csv, 06_env_matched.csv
# Writes: outputs/tangesund/26_ar1_contrasts.csv, 26_ar1_hydrography_terms.csv,
#         26_ar1_mixed_layer.csv, 26_ar1_rho.csv
# Deps:   data.table, mgcv
# Run:    Rscript R/26_ar1_sensitivity.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(mgcv) })

cfg <- list(out_dir = "outputs/tangesund", lat = 58.075, lon = 11.493, surf = 3, deep = 16)

meta <- fread(file.path(cfg$out_dir, "02_sample_meta.csv"))[, .(sample, depth, date, yday, qc_pass)]
chn  <- fread(file.path(cfg$out_dir, "03_chain_metrics.csv"))[, .(sample, datetime, reliable_chain,
             mean_chain, n_chain_roi)]
env  <- fread(file.path(cfg$out_dir, "06_env_matched.csv"))[, .(sample, insitu_sal, insitu_temp, dsigma, mld)]
D <- Reduce(function(a, b) merge(a, b, by = "sample", all.x = TRUE), list(meta, chn, env))
D[, datetime := as.POSIXct(datetime, tz = "UTC")]
D[, `:=`(period = factor(month.abb[as.integer(format(as.Date(date), "%m"))], levels = c("Aug","Sep","Oct")),
         day    = factor(date))]
D <- D[qc_pass == TRUE & reliable_chain == TRUE & is.finite(mean_chain)]

# ---- within-day lag-1 rho (median residual ACF) and AR.start ------------------
prep_ar <- function(dd, m0) {
  setorder(dd, day, datetime)
  dd[, ar_start := !duplicated(day)]
  res <- residuals(m0)
  rho <- dd[, .(r = if (.N > 6) acf(res[.I], plot = FALSE, lag.max = 1)$acf[2] else NA_real_),
            by = day][, median(r, na.rm = TRUE)]
  list(rho = rho, dd = dd)
}
contr_1 <- function(m, dd, pp) {                          # surface-deep contrast in one period
  sub <- dd[period == pp]
  g <- data.frame(depth = c(cfg$surf, cfg$deep), period = factor(pp, levels = levels(dd$period)),
                  yday = median(sub$yday), day = sub$day[1])
  Xp <- predict(m, g, type = "lpmatrix"); cv <- Xp[1, ] - Xp[2, ]
  est <- sum(cv * coef(m)); se <- sqrt(drop(cv %*% vcov(m) %*% cv))
  data.table(period = pp, est = est, se = se, lo = est - 2*se, hi = est + 2*se,
             sig = (est - 2*se > 0) | (est + 2*se < 0), n_days = uniqueN(sub$day))
}

# ---- (A) chain-length contrast: naive, AR1, day x depth aggregate --------------
f_chain <- mean_chain ~ period + s(depth, by = period, k = 5) + s(yday, k = 8) + s(day, bs = "re")

dA <- copy(D); setorder(dA, day, datetime)
m0 <- bam(f_chain, data = dA, method = "fREML")
ar <- prep_ar(dA, m0); dA <- ar$dd
m1 <- bam(f_chain, data = dA, method = "fREML", rho = ar$rho, AR.start = dA$ar_start)

# one abundance-weighted value per day x depth
agg <- D[, .(mean_chain = weighted.mean(mean_chain, n_chain_roi), yday = yday[1], period = period[1]),
         by = .(date, depth)][, day := factor(date)]
mAg <- bam(f_chain, data = agg, method = "fREML")

levs <- levels(D$period)
cA <- rbindlist(list(
  cbind(model = "naive (day-RE)",  rbindlist(lapply(levs, contr_1, m = m0,  dd = dA))),
  cbind(model = "AR1-corrected",   rbindlist(lapply(levs, contr_1, m = m1,  dd = dA))),
  cbind(model = "day x depth aggregate", rbindlist(lapply(levs, contr_1, m = mAg, dd = agg)))))
fwrite(cA, file.path(cfg$out_dir, "26_ar1_contrasts.csv"))

# lag-1 ACF between successive revisits to the same depth
setorder(D, depth, day, datetime)
r_same <- D[, .(r = if (.N > 6) acf(mean_chain - mean(mean_chain), plot = FALSE, lag.max = 1)$acf[2] else NA_real_),
            by = .(depth, day)][, median(r, na.rm = TRUE)]

cat(sprintf("\n== within-day autocorrelation: lag-1 rho (time-ordered) = %.2f | same-depth revisit ACF = %.2f ==\n",
            ar$rho, r_same))
cat("\n== (A) surface(3m)-deep(16m) chain-length contrast, three rigour levels ==\n")
print(cA[, .(model, period, est = round(est, 3), se = round(se, 3),
             CI = sprintf("[%.2f, %.2f]", lo, hi), sig, n_days)])

# ---- (B) hydrography smooths and mixed-layer contrast, naive vs AR1 ----------
dB <- D[is.finite(insitu_temp) & is.finite(mld) & is.finite(insitu_sal) & is.finite(dsigma)]
dB[, in_ml := factor(fifelse(depth <= mld, "In mixed layer", "Below pycnocline"),
                     levels = c("Below pycnocline", "In mixed layer"))]
setorder(dB, day, datetime)
f_hyd <- mean_chain ~ s(insitu_temp, k=6) + s(insitu_sal, k=6) + s(dsigma, k=6) + s(mld, k=6) +
          s(yday, k=8) + ti(insitu_sal, yday, k=c(5,6)) + s(day, bs="re")   # as 10; no shrinkage
mB0 <- bam(f_hyd, data = dB, method = "fREML")
arB <- prep_ar(dB, mB0); dB <- arB$dd
mB1 <- bam(f_hyd, data = dB, method = "fREML", rho = arB$rho, AR.start = dB$ar_start)
Fof <- function(m) { st <- summary(m)$s.table
  sapply(c("s(insitu_sal)","s(insitu_temp)","s(dsigma)","s(mld)","s(yday)"), function(t) st[t,"F"]) }
Pof <- function(m) { st <- summary(m)$s.table
  sapply(c("s(insitu_sal)","s(insitu_temp)","s(dsigma)","s(mld)","s(yday)"), function(t) st[t,"p-value"]) }
Ftab <- data.table(term = c("s(insitu_sal)","s(insitu_temp)","s(dsigma)","s(mld)","s(yday)"),
                   F_naive = round(Fof(mB0), 1), F_AR1 = round(Fof(mB1), 1),
                   p_naive = signif(Pof(mB0), 2), p_AR1 = signif(Pof(mB1), 2))
fwrite(Ftab, file.path(cfg$out_dir, "26_ar1_hydrography_terms.csv"))
f_ml <- mean_chain ~ in_ml + s(insitu_temp, k=6) + s(yday, k=8) + s(day, bs="re")
mlc <- function(use_ar) {
  m <- if (!use_ar) bam(f_ml, data = dB, method = "fREML")
       else         bam(f_ml, data = dB, method = "fREML", rho = arB$rho, AR.start = dB$ar_start)
  summary(m)$p.table["in_mlIn mixed layer", c("Estimate","Std. Error","Pr(>|t|)")]
}
c0 <- mlc(FALSE); c1 <- mlc(TRUE)
cat(sprintf("\n== (B) mechanism smooth F: naive vs AR1 (rho=%.2f) ==\n", arB$rho)); print(Ftab)
cat(sprintf("   mixed-layer vs below:  naive %+.2f (SE %.2f, p=%.3g) | AR1 %+.2f (SE %.2f, p=%.3g)\n",
            c0[1], c0[2], c0[3], c1[1], c1[2], c1[3]))
fwrite(data.table(model = c("naive (day-RE)", "AR1-corrected"), estimate = c(c0[1], c1[1]), se = c(c0[2], c1[2]),
                  lo = c(c0[1], c1[1]) - 1.96 * c(c0[2], c1[2]), hi = c(c0[1], c1[1]) + 1.96 * c(c0[2], c1[2]),
                  p = c(c0[3], c1[3]), rho = c(NA, arB$rho), n = nrow(dB), n_days = uniqueN(dB$day)),
       file.path(cfg$out_dir, "26_ar1_mixed_layer.csv"))

# ---- (C) day-night contrast, naive vs AR1 ------------------------------------
solar_elev <- function(dt) {
  doy <- as.integer(format(dt, "%j")); hutc <- as.integer(format(dt, "%H")) + as.integer(format(dt, "%M"))/60
  decl <- 23.45 * sin(2*pi*(284 + doy)/365) * pi/180
  H <- (15 * (hutc + cfg$lon/15 - 12)) * pi/180
  asin(sin(cfg$lat*pi/180)*sin(decl) + cos(cfg$lat*pi/180)*cos(decl)*cos(H)) * 180/pi
}
dC <- copy(D)[, daynight := factor(fifelse(solar_elev(datetime) > 0, "day", "night"),
                                   levels = c("night", "day"))]
setorder(dC, day, datetime)
f_dn <- mean_chain ~ daynight + s(depth, k=5) + s(yday, k=6) + s(day, bs="re")
mC0 <- bam(f_dn, data = dC, method = "fREML")
arC <- prep_ar(dC, mC0); dC <- arC$dd
mC1 <- bam(f_dn, data = dC, method = "fREML", rho = arC$rho, AR.start = dC$ar_start)
d0 <- summary(mC0)$p.table["daynightday", c("Estimate","Std. Error","Pr(>|t|)")]
d1 <- summary(mC1)$p.table["daynightday", c("Estimate","Std. Error","Pr(>|t|)")]
cat(sprintf("\n== (C) diel day-night contrast (rho=%.2f) ==\n", arC$rho))
cat(sprintf("   naive %+.3f (SE %.3f, p=%.3g) | AR1 %+.3f (SE %.3f, p=%.3g)\n",
            d0[1], d0[2], d0[3], d1[1], d1[2], d1[3]))

fwrite(data.table(model = c("chain-length contrasts (A)", "hydrography (B)", "diel day-night (C)"),
                  rho_lag1_within_day = c(ar$rho, arB$rho, arC$rho),
                  same_depth_revisit_acf = c(r_same, NA, NA), n = c(nrow(dA), nrow(dB), nrow(dC))),
       file.path(cfg$out_dir, "26_ar1_rho.csv"))
cat(sprintf("table -> %s/26_ar1_contrasts.csv\n", cfg$out_dir))
