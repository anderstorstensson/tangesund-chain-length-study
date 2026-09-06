# =============================================================================
# The within/between-taxon partition of 17 on days with both 3 m and 16 m sampled, and
# the 3 m - 16 m GAMM contrasts of 13 without week 5 and without weeks 5-7, plus the
# weekly 3 m and 16 m means behind Fig. 2c.
#
# Reads:  data/roi_level/roi_2016.txt.gz, outputs/tangesund/02_sample_meta.csv,
#         03_chain_metrics.csv
# Writes: outputs/tangesund/17b_matched_day_partition.csv, 17b_sep_without_week5.csv,
#         17b_weekly_surface_deep.csv
# Deps:   data.table, mgcv
# Run:    Rscript R/17b_matched_day_partition.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(mgcv) })
cfg <- list(roi_file = "data/roi_level/roi_2016.txt.gz", meta_csv = "outputs/tangesund/02_sample_meta.csv",
            chain_csv = "outputs/tangesund/03_chain_metrics.csv",
            out_dir = "outputs/tangesund", junk = c("Debris", "Beads", "Air_bubbles", "Zooplankton"),
            surf = 3, deep = 16, min_taxon_roi = 5, nboot = 2000, seed = 1)   # as in 17_plasticity_turnover.R
meta <- fread(cfg$meta_csv); meta[, date := as.Date(date)]
meta[, period := factor(month.abb[as.integer(format(date, "%m"))], levels = c("Aug", "Sep", "Oct"))]

# ---- (A) partition on matched days ------------------------------------------------
d <- fread(cfg$roi_file, select = c("sample", "class_name", "biovolume_um3", "numBlobs", "cell_count", "cell_count_resolved"),
           showProgress = FALSE)
d <- d[numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 & !class_name %in% cfg$junk]
d <- d[!is.na(cell_count) & cell_count != -1]                          # chain-eligible only
d[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1, round(cell_count_resolved), 1L)]
d <- merge(d, meta[qc_pass == TRUE, .(sample, depth, date, period)], by = "sample")[depth %in% c(cfg$surf, cfg$deep)]
d[, day := as.integer(date)]
days_both <- d[, .(both = uniqueN(depth) == 2), by = day][both == TRUE, day]
cat(sprintf("days with chain-eligible ROIs at 3 m: %d | at 16 m: %d | at both: %d\n",
            uniqueN(d[depth == cfg$surf, day]), uniqueN(d[depth == cfg$deep, day]), length(days_both)))
dm <- d[day %in% days_both]
# same additive (Price-style) decomposition as 17
decompose <- function(dat) {
  agg <- dat[, .(L = mean(nchain), nn = .N), by = .(depth, class_name)][nn >= cfg$min_taxon_roi]
  agg[, f := nn / sum(nn), by = depth]
  w <- merge(agg[depth == cfg$surf, .(class_name, L_s = L, f_s = f)], agg[depth == cfg$deep, .(class_name, L_d = L, f_d = f)],
             by = "class_name", all = TRUE)
  for (col in c("L_s", "L_d", "f_s", "f_d")) w[is.na(get(col)), (col) := 0]
  w[, `:=`(fbar = (f_s + f_d)/2, Lbar = (L_s + L_d)/2)]
  within <- sum(w$fbar * (w$L_s - w$L_d)); between <- sum(w$Lbar * (w$f_s - w$f_d))
  list(within = within, between = between, total = within + between)
}
boot_ci <- function(dat, nboot = cfg$nboot) {
  days <- unique(dat$day)
  out <- data.table(within = numeric(nboot), between = numeric(nboot), frac_within = numeric(nboot))
  for (b in seq_len(nboot)) {
    sub <- dat[.(sample(days, length(days), replace = TRUE)), on = "day", allow.cartesian = TRUE, nomatch = NULL]
    r <- decompose(sub); out[b, `:=`(within = r$within, between = r$between, frac_within = r$within / (r$within + r$between))]
  }
  out
}
set.seed(cfg$seed)
qs <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
res <- rbindlist(lapply(c("All", "Aug", "Sep", "Oct"), function(per) {
  dat <- if (per == "All") dm else dm[period == per]
  if (uniqueN(dat$day) < 3) return(NULL)
  r <- decompose(dat); bb <- boot_ci(dat)
  data.table(set = "matched days only", period = per, n_days = uniqueN(dat$day), n_roi = nrow(dat),
             dCWM = r$total, within = r$within, between = r$between, frac_within = 100 * r$within / r$total,
             frac_within_lo = 100 * qs(bb$frac_within)[1], frac_within_hi = 100 * qs(bb$frac_within)[2])
}))
fwrite(res, file.path(cfg$out_dir, "17b_matched_day_partition.csv"))
cat("\n-- partition on days with both 3 m and 16 m sampled (17b_matched_day_partition.csv) --\n")
print(res[, .(period, n_days, n_roi, dCWM = round(dCWM, 2), within = round(within, 2), between = round(between, 2),
              frac_within = round(frac_within), CI = sprintf("[%.0f, %.0f]", frac_within_lo, frac_within_hi))])
rm(d, dm); invisible(gc())

# ---- (B) September contrast without week 5 ---------------------------------------------
# same specification as 13, on reduced week sets
chn <- fread(cfg$chain_csv)[qc_pass == TRUE & reliable_chain == TRUE]
chn[, `:=`(date = as.Date(date), day = factor(date))]
chn[, period := factor(month.abb[as.integer(format(date, "%m"))], levels = c("Aug", "Sep", "Oct"))]
contrast <- function(dd, label) {
  m <- bam(mean_chain ~ period + s(depth, by = period, k = 5) + s(yday, k = 8) + s(day, bs = "re"), data = dd, method = "fREML")
  rbindlist(lapply(levels(dd$period), function(pp) {
    sub <- dd[period == pp]; if (nrow(sub) < 40) return(NULL)
    g <- data.frame(depth = c(3, 16), period = factor(pp, levels = levels(dd$period)), yday = median(sub$yday), day = sub$day[1])
    Xp <- predict(m, g, type = "lpmatrix"); cv <- Xp[1, ] - Xp[2, ]
    est <- sum(cv * coef(m)); se <- sqrt(drop(cv %*% vcov(m) %*% cv))
    data.table(set = label, period = pp, n = nrow(sub), n_days = uniqueN(sub$day), est = est, lo = est - 2 * se, hi = est + 2 * se)
  }))
}
sens <- rbind(contrast(chn, "all weeks"), contrast(chn[week != 5], "without week 5"),
              contrast(chn[!week %in% 5:7], "without weeks 5-7"))
fwrite(sens, file.path(cfg$out_dir, "17b_sep_without_week5.csv"))
cat("\n-- surface-minus-deep contrasts by period, with weeks removed (17b_sep_without_week5.csv) --\n")
print(sens[, .(set, period, n, n_days, est = round(est, 2), CI = sprintf("[%.2f, %.2f]", lo, hi))])
wk <- chn[depth %in% c(3, 16), .(n = .N, n_days = uniqueN(day), mean_chain = mean(mean_chain)), by = .(week, depth)][order(week, depth)]
fwrite(wk, file.path(cfg$out_dir, "17b_weekly_surface_deep.csv"))
cat("\n-- weekly 3 m and 16 m means (17b_weekly_surface_deep.csv) --\n")
print(dcast(wk[, .(week, depth, v = sprintf("%.2f (n=%d, %d d)", mean_chain, n, n_days))], week ~ depth))
