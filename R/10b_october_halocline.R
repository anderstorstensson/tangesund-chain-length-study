# =============================================================================
# Monthly chain-length profiles, 6-16 m depth and salinity slopes, 3 m minus 6 m
# contrasts and buoy salinity structure around the October halocline shoaling (Fig. S4).
#
# Reads:  outputs/tangesund/03_chain_metrics.csv, 06_buoy_profile_at_ifcb_depths.csv,
#         06_env_matched.csv.
# Writes: outputs/tangesund/10b_monthly_profile.csv, 10b_monthly_tests.csv,
#         10b_halocline_by_month.csv, 10b_salinity_denominators.csv,
#         10b_slope_per_salinity.csv, 10b_slope_boot.csv, 10b_shoaling_date.csv (only if a
#         shoaling date exists), 10b_october_panel.png.
# Env:    HALOCLINE_FIG_ONLY=1 redraws Fig. S4 from 10b_slope_boot.csv instead of
#         bootstrapping (used by 28).
# Deps:   data.table, mgcv, ggplot2, patchwork, parallel.
# Run:    Rscript R/10b_october_halocline.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(patchwork); library(parallel) })
source("R/lib_plotting_theme.R")
cfg <- list(out_dir = "outputs/tangesund", depths = c(3, 6, 8, 11, 13, 16), B = 2000,
            cores = min(6L, max(1L, detectCores() - 2L)),
            fig_only = nzchar(Sys.getenv("HALOCLINE_FIG_ONLY")))
chn <- fread(file.path(cfg$out_dir, "03_chain_metrics.csv"))[qc_pass == TRUE & reliable_chain == TRUE]
chn <- merge(chn, fread(file.path(cfg$out_dir, "06_env_matched.csv"))[, .(sample, insitu_sal)], by = "sample", all.x = TRUE)
chn[, `:=`(date = as.Date(date), day = factor(date))]
chn[, month := factor(month.abb[as.integer(format(date, "%m"))], levels = c("Aug", "Sep", "Oct"))]

boot_se_day <- function(value, day, B = cfg$B, seed = 1) {
  days <- unique(day); nd <- length(days); if (nd < 2) return(NA_real_)
  idx <- split(seq_along(value), day); set.seed(seed)
  sd(vapply(seq_len(B), function(b) mean(value[unlist(idx[as.character(sample(days, nd, TRUE))], use.names = FALSE)]),
            numeric(1)))
}
prof <- chn[, .(n = .N, n_days = uniqueN(day), mean_chain = mean(mean_chain),
                se_day = boot_se_day(mean_chain, as.character(day))), by = .(month, depth)][order(month, depth)]
fwrite(prof, file.path(cfg$out_dir, "10b_monthly_profile.csv"))
cat("-- chain length by month and depth (10b_monthly_profile.csv) --\n")
print(dcast(prof[, .(month, depth, v = sprintf("%.2f ± %.2f (n=%d)", mean_chain, se_day, n))], depth ~ month))

# ---- tests per month --------------------------------------------------------
tests <- rbindlist(lapply(levels(chn$month), function(mo) {
  x <- chn[month == mo]
  deep <- x[depth >= 6]
  m1 <- bam(mean_chain ~ depth + s(day, bs = "re"), data = deep, method = "fREML")
  s1 <- summary(m1)$p.table["depth", ]
  # log scale
  m1l <- bam(log(mean_chain) ~ depth + s(day, bs = "re"), data = deep, method = "fREML")
  s1l <- summary(m1l)$p.table["depth", ]
  # log chain length on matched per-sample salinity over 6-16 m
  deep_s <- deep[is.finite(insitu_sal)]
  m1s <- bam(log(mean_chain) ~ insitu_sal + s(day, bs = "re"), data = deep_s, method = "fREML")
  s1s <- summary(m1s)$p.table["insitu_sal", ]
  # depth as a factor over 6-16 m, F test vs day-only
  deep[, depth_f := factor(depth)]
  m0 <- bam(mean_chain ~ s(day, bs = "re"), data = deep, method = "fREML")
  mf <- bam(mean_chain ~ depth_f + s(day, bs = "re"), data = deep, method = "fREML")
  av <- anova(m0, mf, test = "F")
  sh <- x[depth <= 6]; sh[, is3 := factor(depth == 3, levels = c(FALSE, TRUE))]
  m3 <- bam(mean_chain ~ is3 + s(day, bs = "re"), data = sh, method = "fREML")
  s3 <- summary(m3)$p.table["is3TRUE", ]
  rbind(
    data.table(month = mo, test = "depth slope over 6-16 m (cells per m)", n = nrow(deep), n_days = uniqueN(deep$day),
               estimate = s1["Estimate"], se = s1["Std. Error"], lo = s1["Estimate"] - 1.96 * s1["Std. Error"],
               hi = s1["Estimate"] + 1.96 * s1["Std. Error"], statistic = s1["t value"], df = NA_real_, p = s1["Pr(>|t|)"]),
    data.table(month = mo, test = "log depth slope over 6-16 m (per m)", n = nrow(deep), n_days = uniqueN(deep$day),
               estimate = s1l["Estimate"], se = s1l["Std. Error"], lo = s1l["Estimate"] - 1.96 * s1l["Std. Error"],
               hi = s1l["Estimate"] + 1.96 * s1l["Std. Error"], statistic = s1l["t value"], df = NA_real_, p = s1l["Pr(>|t|)"]),
    data.table(month = mo, test = "log slope on matched salinity over 6-16 m (per salinity unit)", n = nrow(deep_s), n_days = uniqueN(deep_s$day),
               estimate = s1s["Estimate"], se = s1s["Std. Error"], lo = s1s["Estimate"] - 1.96 * s1s["Std. Error"],
               hi = s1s["Estimate"] + 1.96 * s1s["Std. Error"], statistic = s1s["t value"], df = NA_real_, p = s1s["Pr(>|t|)"]),
    data.table(month = mo, test = "depth as a factor over 6-16 m (F test vs day-only)", n = nrow(deep), n_days = uniqueN(deep$day),
               estimate = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_,
               statistic = av$F[2], df = av$Df[2], p = av$`Pr(>F)`[2]),
    data.table(month = mo, test = "3 m minus 6 m (cells)", n = nrow(sh), n_days = uniqueN(sh$day),
               estimate = s3["Estimate"], se = s3["Std. Error"], lo = s3["Estimate"] - 1.96 * s3["Std. Error"],
               hi = s3["Estimate"] + 1.96 * s3["Std. Error"], statistic = s3["t value"], df = NA_real_, p = s3["Pr(>|t|)"]))
}))
fwrite(tests, file.path(cfg$out_dir, "10b_monthly_tests.csv"))
cat("\n-- per-month tests (10b_monthly_tests.csv) --\n")
print(tests[, .(month, test, n, n_days, estimate = round(estimate, 3), CI = ifelse(is.na(lo), "", sprintf("[%.3f, %.3f]", lo, hi)),
                statistic = round(statistic, 2), p = signif(p, 3))])

# ---- the buoy's salinity structure per month ----------------------------------
bp <- fread(file.path(cfg$out_dir, "06_buoy_profile_at_ifcb_depths.csv"))
bp[, month := factor(format(time, "%b"), levels = c("Aug", "Sep", "Oct"))]
w <- dcast(bp, time + month ~ depth, value.var = "sal")
setnames(w, as.character(cfg$depths), paste0("s", cfg$depths))
w[, `:=`(d3_6 = s3 - s6, d6_16 = s6 - s16, d3_16 = s3 - s16)]
# steepest salinity step between consecutive target depths
steps <- bp[order(time, depth), .(step_top = depth[which.max(diff(sal) / diff(depth))]), by = .(time, month)]
hal <- merge(w[, .(hours = .N, sal_3m = median(s3), sal_6m = median(s6), sal_16m = median(s16),
                   diff_3_6 = median(d3_6), diff_6_16 = median(d6_16),
                   share_of_gradient_above_6m = median(d3_6 / d3_16, na.rm = TRUE)), by = month],
             steps[, .(pct_steepest_step_3_6 = 100 * mean(step_top == 3)), by = month], by = "month")
fwrite(hal, file.path(cfg$out_dir, "10b_halocline_by_month.csv"))
# candidate 6-16 m salinity denominators and within-day depth-salinity correlation
salmed_m <- w[, .(diff_of_monthly_medians_6_16 = median(s6) - median(s16),
                  diff_of_monthly_medians_3_6 = median(s3) - median(s6)), by = month]
matched_m <- chn[is.finite(insitu_sal), .(matched_median_diff_6_16 = median(insitu_sal[depth == 16]) - median(insitu_sal[depth == 6]),
                                          matched_median_diff_3_6 = median(insitu_sal[depth == 6]) - median(insitu_sal[depth == 3]),
                                          n_samples = .N), by = month]
cor_m <- chn[depth >= 6 & is.finite(insitu_sal), .(r = if (.N > 5) cor(depth, insitu_sal) else NA_real_), by = .(month, day)][
  , .(median_within_day_cor_depth_matched_sal_6_16 = median(r, na.rm = TRUE)), by = month]
den <- Reduce(function(a, b) merge(a, b, by = "month"),
              list(hal[, .(month, buoy_median_hourly_diff_6_16 = -diff_6_16, buoy_median_hourly_diff_3_6 = -diff_3_6)],
                   salmed_m[, .(month, diff_of_monthly_medians_6_16 = -diff_of_monthly_medians_6_16,
                                diff_of_monthly_medians_3_6 = -diff_of_monthly_medians_3_6)],
                   matched_m, cor_m))
den[, month := factor(month, levels = c("Aug", "Sep", "Oct"))]; setorder(den, month)
fwrite(den, file.path(cfg$out_dir, "10b_salinity_denominators.csv"))
cat("\n-- 6-16 m salinity denominators by month, positive = saltier at depth (10b_salinity_denominators.csv) --\n")
print(den[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
# slope per salinity unit
sps <- merge(tests[test == "depth slope over 6-16 m (cells per m)", .(month, slope_cells_per_m = estimate)],
             tests[test == "log depth slope over 6-16 m (per m)", .(month, log_slope_per_m = estimate)], by = "month")
sps <- merge(sps, hal[, .(month, sal_diff_6_16 = -diff_6_16)], by = "month")
sps <- merge(sps, tests[test == "log slope on matched salinity over 6-16 m (per salinity unit)",
                        .(month, log_slope_per_matched_salinity_unit = estimate,
                          log_slope_per_matched_salinity_unit_lo = lo, log_slope_per_matched_salinity_unit_hi = hi)], by = "month")
sps[, `:=`(cells_per_m_per_salinity_unit = slope_cells_per_m / sal_diff_6_16,      # buoy-median denominator
           pct_decline_6_to_16 = 100 * (1 - exp(10 * log_slope_per_m)))]
sps[, month := factor(month, levels = c("Aug", "Sep", "Oct"))]; setorder(sps, month)
fwrite(sps, file.path(cfg$out_dir, "10b_slope_per_salinity.csv"))
cat("\n-- slope per salinity unit (10b_slope_per_salinity.csv) --\n"); print(sps)
cat("\n-- buoy salinity structure by month (10b_halocline_by_month.csv; medians over hourly profiles) --\n")
print(hal[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

# ---- day-cluster bootstrap of the slopes and derived ratios -------------------
# Days resampled within period and relabeled so the day RE stays identifiable; the
# salinity difference is fixed. Slopes are per m; "per salinity unit" = change over
# the 10 m layer divided by its salinity difference.
slope_stats <- function(x, sal6_16) {
  deep <- x[depth >= 6]
  b  <- coef(bam(mean_chain ~ depth + s(day, bs = "re"), data = deep, method = "fREML"))[["depth"]]
  bl <- coef(bam(log(mean_chain) ~ depth + s(day, bs = "re"), data = deep, method = "fREML"))[["depth"]]
  deep_s <- deep[is.finite(insitu_sal)]
  bs <- coef(bam(log(mean_chain) ~ insitu_sal + s(day, bs = "re"), data = deep_s, method = "fREML"))[["insitu_sal"]]
  bsc <- coef(bam(mean_chain ~ insitu_sal + s(day, bs = "re"), data = deep_s, method = "fREML"))[["insitu_sal"]]
  pm <- x[, .(m = mean(mean_chain)), by = depth]
  mof <- function(dp) if (dp %in% pm$depth) pm[depth == dp, m] else NA_real_
  d3 <- mof(3); d6 <- mof(6); d16 <- mof(16)
  c(slope_cells_per_m = b, log_slope_per_m = bl,
    cells_per_m_per_salinity_unit = b / sal6_16,
    dchain_per_salinity_unit_6_16 = 10 * b / sal6_16,
    dlog_per_salinity_unit_6_16 = 10 * bl / sal6_16,          # buoy-median denominator
    dlog_per_matched_salinity_unit_6_16 = bs,                 # regression on matched per-sample salinity
    dchain_per_matched_salinity_unit_6_16 = bsc,
    pct_decline_6_to_16 = 100 * (1 - exp(10 * bl)),
    share_of_chain_difference_above_6m = (d3 - d6) / (d3 - d16))
}
if (!cfg$fig_only) {
RNGkind("L'Ecuyer-CMRG")
boot_period <- function(x, lab, sal6_16, B = cfg$B, seed = 1) {
  x <- copy(x)[, day := factor(as.character(day))]
  days <- levels(x$day); obs <- slope_stats(x, sal6_16)
  set.seed(seed)
  draws <- mclapply(seq_len(B), function(b) {
    pick <- sample(days, length(days), replace = TRUE)
    xb <- rbindlist(lapply(seq_along(pick), function(i) copy(x[day == pick[i]])[, day := paste0("b", i)]))
    xb[, day := factor(day)]
    r <- try(slope_stats(xb, sal6_16), silent = TRUE)
    if (inherits(r, "try-error")) rep(NA_real_, length(obs)) else r
  }, mc.cores = cfg$cores, mc.set.seed = TRUE)
  dr <- do.call(rbind, draws)
  data.table(period = lab, quantity = names(obs), estimate = obs,
             lo = apply(dr, 2, quantile, .025, na.rm = TRUE), hi = apply(dr, 2, quantile, .975, na.rm = TRUE),
             n = nrow(x), n_days = length(days), sal_diff_6_16 = sal6_16, n_draws = sum(complete.cases(dr)))
}
sb <- rbindlist(lapply(levels(chn$month), function(mo)
  boot_period(chn[month == mo], mo, hal[month == mo, -diff_6_16])))

# ---- October split at the shoaling date ----------------------------------------
# first day from which the steepest salinity step is at 3-6 m in > 50% of hourly
# profiles on that day and every later day
dstep <- steps[, .(frac_3_6 = mean(step_top == 3), hours = .N), by = .(date = as.Date(time))][order(date)]
dstep[, all_later := rev(cumsum(!rev(frac_3_6 > 0.5)) == 0)]
shoal_date <- dstep[all_later == TRUE, min(date)]
if (is.finite(shoal_date)) {
  oct <- chn[month == "Oct"]
  wsub <- w[month == "Oct"]
  pre  <- oct[date <  shoal_date]; post <- oct[date >= shoal_date]
  cat(sprintf("\nhalocline shoaling date (steepest step at 3-6 m in > 50%% of hours from then on): %s; October days before %d, from %d\n",
              format(shoal_date), uniqueN(pre$day), uniqueN(post$day)))
  if (uniqueN(pre$day) >= 4)  sb <- rbind(sb, boot_period(pre,  "Oct before shoaling", wsub[as.Date(time) <  shoal_date, -median(d6_16)]))
  if (uniqueN(post$day) >= 4) sb <- rbind(sb, boot_period(post, "Oct from shoaling",   wsub[as.Date(time) >= shoal_date, -median(d6_16)]))
  fwrite(data.table(shoaling_date = shoal_date, october_days_before = uniqueN(pre$day), october_days_from = uniqueN(post$day),
                    samples_before = nrow(pre), samples_from = nrow(post)),
         file.path(cfg$out_dir, "10b_shoaling_date.csv"))
}
fwrite(sb, file.path(cfg$out_dir, "10b_slope_boot.csv"))
cat(sprintf("\n-- slopes and derived quantities with day-cluster bootstrap intervals (10b_slope_boot.csv, %d draws) --\n", cfg$B))
print(sb[, .(period, quantity, estimate = round(estimate, 3), CI = sprintf("[%.3f, %.3f]", lo, hi), n_days)])
} else {
  sb <- fread(file.path(cfg$out_dir, "10b_slope_boot.csv"))
  cat("\nHALOCLINE_FIG_ONLY set: the bootstrap is read from 10b_slope_boot.csv\n")
}

# ---- figure (Fig. S4) ---------------------------------------------------------
pal <- c(Aug = "#E69F00", Sep = "#0072B2", Oct = "#009E73")
pa <- ggplot(prof, aes(mean_chain, depth, colour = month, shape = month)) +
  geom_errorbar(aes(xmin = mean_chain - se_day, xmax = mean_chain + se_day), orientation = "y", width = 0.4, show.legend = FALSE) +
  geom_path(orientation = "y") + geom_point(size = 2.2) +
  scale_y_reverse(breaks = cfg$depths) + scale_colour_manual(values = pal, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  theme_paper() + theme(legend.position = "top") +
  labs(x = "Mean chain length (cells; ± SE clustered by day)", y = "Depth (m)", title = "Chain length by month")
salm <- melt(w[, .(month, time, s3, s6, s8, s11, s13, s16)], id.vars = c("month", "time"))
salm[, depth := as.numeric(sub("s", "", variable))]
salmed <- salm[, .(sal = median(value)), by = .(month, depth)][order(month, depth)]
pb <- ggplot(salmed, aes(sal, depth, colour = month, shape = month)) +
  geom_path(orientation = "y") + geom_point(size = 2.2) +
  scale_y_reverse(breaks = cfg$depths) + scale_colour_manual(values = pal, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  theme_paper() + theme(legend.position = "top") +
  labs(x = "Salinity (median of hourly buoy profiles)", y = "Depth (m)", title = "Salinity by month")
sl <- tests[test == "depth slope over 6-16 m (cells per m)"]
sl[, month := factor(month, levels = c("Aug", "Sep", "Oct"))]
pc <- ggplot(sl, aes(estimate, month, colour = month, shape = month)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.25) + geom_point(size = 2.4) +
  scale_colour_manual(values = pal, guide = "none") + scale_shape_manual(values = c(16, 17, 15), guide = "none") +
  scale_y_discrete(limits = rev) +
  theme_paper() +
  labs(x = "Depth slope, 6–16 m (cells per m, 95% CI)", y = NULL, title = "Depth slope below 6 m")
# (d) log slope, slope per matched salinity unit, 6-to-16 m decline and the share of the
# 3-16 m difference above 6 m, with their day-cluster bootstrap intervals
qd <- c(log_slope_per_m = "Log slope per m, 6\u201316 m",
        dlog_per_matched_salinity_unit_6_16 = "Log slope per matched salinity unit",
        pct_decline_6_to_16 = "Decline from 6 to 16 m (%)",
        share_of_chain_difference_above_6m = "Share above 6 m")
sd_ <- sb[period %in% levels(chn$month) & quantity %in% names(qd)]
sd_[, `:=`(month = factor(period, levels = c("Aug", "Sep", "Oct")), quantity = factor(qd[quantity], levels = qd))]
pd <- ggplot(sd_, aes(estimate, month, colour = month, shape = month)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.25) + geom_point(size = 2.4) +
  facet_wrap(~ quantity, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = pal, guide = "none") + scale_shape_manual(values = c(16, 17, 15), guide = "none") +
  scale_y_discrete(limits = rev) +
  theme_paper() +
  labs(x = "Estimate with day-cluster bootstrap 95% CI", y = NULL, title = "Slopes on the log scale and per salinity unit")
ggsave(file.path(cfg$out_dir, "10b_october_panel.png"),
       (((pa + labs(tag = "a")) | (pb + labs(tag = "b")) | (pc + labs(tag = "c"))) +
          plot_layout(widths = c(1.2, 1.2, 1), guides = "collect")) / (pd + labs(tag = "d")) +
         plot_layout(heights = c(1.3, 1)) &
         theme(legend.position = "top", plot.tag = element_text(face = "bold", size = 13)),
       width = 12, height = 8.2, dpi = FIG_DPI)
cat(sprintf("\nfigure -> %s/10b_october_panel.png\n", cfg$out_dir))
