# =============================================================================
# Hydrography from SMHI buoy 33027 (TANGESUND BOJ, at the IFCB station): hourly
# stratification metrics, column descriptives, and in-situ salinity, temperature,
# oxygen and current matched to each IFCB sample (Fig. S1).
#
# Reads:  outputs/tangesund/02_sample_meta.csv; the four buoy CSVs in data/tångesund_bouy.
# Writes: outputs/tangesund/06_buoy_profile_at_ifcb_depths.csv, 06_env_matched.csv,
#         06_column_descriptives_by_depth.csv, 06_column_descriptives.csv, the
#         time-depth sections and stratification time series (06_*.png),
#         06_hydrographic_profile.png and 06_buoy_panel.png.
# Deps:   data.table, oce, ggplot2, patchwork.
# Run:    Rscript R/06_env_buoy.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")  # English month names on date axes
suppressPackageStartupMessages({ library(data.table); library(oce); library(ggplot2); library(patchwork) })
source("R/lib_plotting_theme.R")

cfg <- list(
  buoy_dir = "data/tångesund_bouy",
  meta_csv = "outputs/tangesund/02_sample_meta.csv",
  out_dir  = "outputs/tangesund",
  files = list(salinity = "smhi-opendata_4_33027_20260702_055347.csv",
               temp     = "smhi-opendata_5_33027_20260702_055238.csv",
               oxygen   = "smhi-opendata_15_33027_20260702_055403.csv",
               current  = "smhi-opendata_2_3_33027_20260702_055356.csv"),
  lat = 58.0747, match_tol_min = 90,        # accept buoy hour within +-90 min of a sample
  mld_crit = 0.125                           # sigma-t change (kg/m3) defining the mixed layer
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- parse one SMHI buoy CSV (7 header lines; time;value;qual;depth) -> long ---
read_buoy1 <- function(file, varname) {
  d <- fread(file.path(cfg$buoy_dir, file), skip = 7, sep = ";", header = FALSE,
             fill = TRUE, select = 1:4, na.strings = c("",""), showProgress = FALSE)
  setnames(d, c("time","val","qual","depth"))
  d <- d[!is.na(val) & !is.na(depth)]
  d[, `:=`(time = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
           val = as.numeric(val), depth = as.numeric(depth), variable = varname)]
  d[is.finite(val) & is.finite(depth) & !is.na(time), .(time, depth, variable, val, qual)]
}
# currents: two params per row (dir col2, speed col5; depth col4)
read_current <- function(file) {
  d <- fread(file.path(cfg$buoy_dir, file), skip = 7, sep = ";", header = FALSE,
             fill = TRUE, select = c(1,2,4,5), na.strings = c("",""), showProgress = FALSE)
  setnames(d, c("time","cur_dir","depth","cur_speed"))
  d[, `:=`(time = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"), depth = as.numeric(depth),
           cur_dir = as.numeric(cur_dir), cur_speed = as.numeric(cur_speed))]
  rbind(d[is.finite(cur_speed), .(time, depth, variable = "cur_speed", val = cur_speed, qual = "G")],
        d[is.finite(cur_dir),   .(time, depth, variable = "cur_dir",   val = cur_dir,   qual = "G")])
}

buoy <- rbindlist(list(
  read_buoy1(cfg$files$salinity, "sal"),
  read_buoy1(cfg$files$temp,     "temp"),
  read_buoy1(cfg$files$oxygen,   "oxy"),
  read_current(cfg$files$current)), use.names = TRUE)
buoy <- buoy[qual %in% c("G","Y")]                          # good and aggregated
cat(sprintf("buoy rows: %s | %s to %s | depths %s\n", format(nrow(buoy), big.mark=","),
            as.Date(min(buoy$time)), as.Date(max(buoy$time)),
            paste(sort(unique(buoy$depth)), collapse=",")))

# ---- per-hour column stratification metrics (from T & S profiles) ------------
ts <- dcast(buoy[variable %in% c("sal","temp")], time + depth ~ variable, value.var = "val",
            fun.aggregate = mean)
ts <- ts[is.finite(sal) & is.finite(temp)]
ts[, sigma_t := swSigmaT(salinity = sal, temperature = temp, pressure = depth)]
setorder(ts, time, depth)
strat <- ts[, {
  o <- order(depth); dz <- depth[o]; sg <- sigma_t[o]; sa <- sal[o]; te <- temp[o]
  surf <- which.min(dz); deep <- which.max(dz)
  # MLD: shallowest depth where sigma_t exceeds surface sigma_t by mld_crit (interp)
  dsig <- sg - sg[surf]
  mld <- if (any(dsig >= cfg$mld_crit)) {
    j <- which(dsig >= cfg$mld_crit)[1]
    if (j == 1) dz[1] else approx(dsig[(j-1):j], dz[(j-1):j], xout = cfg$mld_crit)$y
  } else max(dz)
  # max buoyancy frequency proxy: largest sigma_t gradient between adjacent depths
  n2 <- if (length(dz) >= 2) max(diff(sg) / diff(dz)) else NA_real_
  .(n_depth = length(dz), surf_temp = te[surf], surf_sal = sa[surf], surf_sigma = sg[surf],
    dsigma = sg[deep] - sg[surf], mld = mld, dsigma_dz_max = n2)
}, by = time]
cat(sprintf("stratification hours: %d | median MLD %.1f m, median surf-bottom d-sigma-t %.2f kg/m3\n",
            nrow(strat), median(strat$mld, na.rm=TRUE), median(strat$dsigma, na.rm=TRUE)))

# ---- the column itself, at the IFCB target depths (Fig. 6a) ------------------
# every hourly buoy profile interpolated to all six depths; each depth covers the same
# hours (medians over IFCB-matched rows mix in sampling composition)
ifcb_depths <- c(3, 6, 8, 11, 13, 16)
prof_hourly <- ts[time >= as.POSIXct("2016-08-10", tz = "UTC") &
                  time <= as.POSIXct("2016-10-17", tz = "UTC")][, {
  o <- order(depth)
  if (uniqueN(depth) >= 4)
    .(depth = ifcb_depths,
      sal  = approx(depth[o], sal[o],  xout = ifcb_depths, rule = 2)$y,
      temp = approx(depth[o], temp[o], xout = ifcb_depths, rule = 2)$y)
  else NULL
}, by = time]
prof_hourly[, sigt := swSigmaT(salinity = sal, temperature = temp, pressure = depth)]
setorder(prof_hourly, time, depth)
fwrite(prof_hourly, file.path(cfg$out_dir, "06_buoy_profile_at_ifcb_depths.csv"))
cat(sprintf("column profile at IFCB depths: %d hourly profiles | 3 m S=%.1f T=%.1f -> 16 m S=%.1f T=%.1f\n",
            uniqueN(prof_hourly$time),
            prof_hourly[depth == 3,  median(sal)],  prof_hourly[depth == 3,  median(temp)],
            prof_hourly[depth == 16, median(sal)],  prof_hourly[depth == 16, median(temp)]))

# ---- column descriptives (10 Aug - 17 Oct 2016, medians over hourly profiles) --
win0 <- as.POSIXct("2016-08-10", tz = "UTC"); win1 <- as.POSIXct("2016-10-17", tz = "UTC")
# (a) per target depth: median and range of the interpolated hourly profile
desc_depth <- prof_hourly[, .(n_hours = .N,
  sal_med = median(sal), sal_min = min(sal), sal_max = max(sal),
  temp_med = median(temp), temp_min = min(temp), temp_max = max(temp),
  sigt_med = median(sigt)), by = depth][order(depth)]
# per-month temperature contrast 3 m - 16 m (medians of the hourly profiles)
tm <- dcast(prof_hourly[depth %in% c(3, 16)][, month := format(time, "%b")],
            month + time ~ depth, value.var = "temp")
setnames(tm, c("3", "16"), c("t3", "t16"))
temp_month <- tm[, .(temp_3m_minus_16m = median(t3 - t16)), by = month]
# vertical temperature range over the six target depths per hourly profile, by month
trange <- prof_hourly[, .(rng = max(temp) - min(temp), d3_8 = temp[depth == 3] - temp[depth == 8]),
                      by = .(time, month = format(time, "%b"))]
trange_month <- trange[, .(range_med = median(rng), range_p90 = quantile(rng, .9), d3_8_med = median(d3_8)), by = month]
# (b) per-hour column stratification within the window, from the raw buoy profiles
tsw <- ts[time >= win0 & time <= win1]
strat_w <- tsw[, if (uniqueN(depth) < 4) NULL else {          # same profile rule as prof_hourly
  o <- order(depth); dz <- depth[o]; sg <- sigma_t[o]
  dsig <- sg - sg[1]
  mld <- if (any(dsig >= cfg$mld_crit)) {
    j <- which(dsig >= cfg$mld_crit)[1]
    if (j == 1) dz[1] else approx(dsig[(j-1):j], dz[(j-1):j], xout = cfg$mld_crit)$y
  } else max(dz)
  g <- diff(sg) / diff(dz); k <- which.max(g)
  .(dsigma = sg[length(dz)] - sg[1], mld = mld, step_max = max(g),
    step_bin = paste0(dz[k], "-", dz[k + 1]))
}, by = time]
strat_w[, month := format(time, "%b")]
step_bins <- strat_w[, .N, by = step_bin][order(-N)]
dsig_month <- strat_w[, .(dsigma_med = median(dsigma), mld_med = median(mld)), by = month]
# (c) weekly column temperature (all buoy depths; week 1 starts 10 Aug)
tsw[, week := as.integer(floor(as.numeric(as.Date(time) - as.Date("2016-08-10")) / 7)) + 1L]
temp_week <- tsw[, .(temp_med = median(temp)), by = week][order(week)]
# (d) current speed over the sampled depths (3-16 m) within the window
cw <- buoy[variable == "cur_speed" & time >= win0 & time <= win1 & depth >= 3 & depth <= 16]
cur_q <- cw[, .(n = .N, median = median(val), p90 = quantile(val, 0.9), max = max(val))]
long <- rbindlist(list(
  data.table(block = "column_hourly", stat = c("n_hours", "dsigma_median", "mld_median",
                                               "step_max_median", "step_max_max"),
             value = c(nrow(strat_w), median(strat_w$dsigma), median(strat_w$mld),
                       median(strat_w$step_max), max(strat_w$step_max))),
  data.table(block = "steepest_step_bin_hours", stat = step_bins$step_bin, value = step_bins$N),
  data.table(block = "dsigma_median_by_month", stat = dsig_month$month, value = dsig_month$dsigma_med),
  data.table(block = "mld_median_by_month",    stat = dsig_month$month, value = dsig_month$mld_med),
  data.table(block = "temp_3m_minus_16m_median_by_month", stat = temp_month$month,
             value = temp_month$temp_3m_minus_16m),
  data.table(block = "column_temp_range_3_16m_median_by_month", stat = trange_month$month, value = trange_month$range_med),
  data.table(block = "column_temp_range_3_16m_p90_by_month", stat = trange_month$month, value = trange_month$range_p90),
  data.table(block = "temp_3m_minus_8m_median_by_month", stat = trange_month$month, value = trange_month$d3_8_med),
  data.table(block = "column_temp_median_by_week", stat = paste0("week", temp_week$week),
             value = temp_week$temp_med),
  data.table(block = "current_3_16m", stat = names(cur_q), value = unlist(cur_q))))
fwrite(desc_depth, file.path(cfg$out_dir, "06_column_descriptives_by_depth.csv"))
fwrite(long,       file.path(cfg$out_dir, "06_column_descriptives.csv"))
cat("\n-- column descriptives (06_column_descriptives_by_depth.csv / 06_column_descriptives.csv) --\n")
print(desc_depth[, lapply(.SD, function(x) round(x, 2))])
print(long[, .(block, stat, value = round(value, 3))], nrows = 200)

# ---- match IFCB samples to nearest buoy hour + interpolate to sample depth ----
# each variable is matched on its own timestamps (the time bases differ)
meta <- fread(cfg$meta_csv); meta[, datetime := as.POSIXct(datetime, tz = "UTC")]

# nearest-time (within tol) roll join of samples onto one variable's timestamps
match_times <- function(times) {
  vt <- data.table(time = sort(unique(times)))[, tkey := time]
  s  <- meta[, .(sample, depth, datetime)][, tkey := datetime]
  setkey(vt, tkey); setkey(s, tkey)
  s <- vt[s, roll = "nearest"]
  s[abs(as.numeric(difftime(datetime, time, units = "mins"))) <= cfg$match_tol_min,
    .(sample, depth, datetime, time)]
}
# interpolate a variable's profile at each sample's matched time to the sample depth
interp_var <- function(vname) {
  bp <- buoy[variable == vname, .(val = mean(val)), by = .(time, bdepth = depth)]
  s  <- match_times(bp$time)
  x  <- merge(s, bp, by = "time", allow.cartesian = TRUE)
  x[, .(value = if (uniqueN(bdepth) >= 2) approx(bdepth, val, xout = depth[1], rule = 2)$y
                else val[which.min(abs(bdepth - depth[1]))]),
    by = .(sample, depth)][, .(sample, variable = vname, value)]
}
envlong <- rbindlist(lapply(c("sal","temp","oxy","cur_speed"), interp_var))
envwide <- dcast(envlong, sample ~ variable, value.var = "value")
setnames(envwide, c("sal","temp","oxy","cur_speed"),
         c("insitu_sal","insitu_temp","insitu_oxy","insitu_cur"), skip_absent = TRUE)
for (v in c("insitu_sal","insitu_temp","insitu_oxy","insitu_cur"))
  cat(sprintf("  %-12s matched %d/%d samples\n", v, sum(is.finite(envwide[[v]])), nrow(meta)))

# attach column stratification metrics: nearest strat hour (within tol) per sample
sm <- match_times(strat$time)[, .(sample, buoy_time = time)][!duplicated(sample)]
sm[, dt_min := NA_real_]
env <- merge(sm, envwide, by = "sample", all = TRUE)
env <- merge(env, strat[, .(time, surf_temp, surf_sal, dsigma, mld, dsigma_dz_max)],
             by.x = "buoy_time", by.y = "time", all.x = TRUE)
env <- merge(meta[, .(sample, depth, datetime, date, week, yday)], env, by = "sample")
fwrite(env, file.path(cfg$out_dir, "06_env_matched.csv"))
cat(sprintf("env-matched table: %d samples | in-situ T/S/O2/current + stratification\n", nrow(env)))

# ---- descriptive figures: buoy time-depth sections (Aug-Oct 2016) ------------
win <- ts[time >= as.POSIXct("2016-08-08", tz="UTC") & time <= as.POSIXct("2016-10-19", tz="UTC")]
sec <- function(var, lab, opt = "D") ggplot(win, aes(as.Date(time), depth, fill = get(var))) +
  # y ticks at the buoy sensor depths
  geom_tile() + scale_y_reverse(breaks = sort(unique(win$depth))) +
  scale_fill_viridis_c(option = opt, name = lab) +
  theme_paper() + labs(x = NULL, y = "Depth (m)", title = lab)
p1 <- sec("temp", "Temperature (°C)", "B")
p2 <- sec("sal", "Salinity", "D")
p3 <- sec("sigma_t", "Density anomaly σt (kg m⁻³)", "G")
sw <- strat[time >= as.POSIXct("2016-08-08", tz="UTC") & time <= as.POSIXct("2016-10-19", tz="UTC")]
# daily means (hourly values on a date axis plot as jagged noise)
swd <- sw[, .(mld = mean(mld, na.rm = TRUE), dsigma = mean(dsigma, na.rm = TRUE)),
          by = .(date = as.Date(time))][order(date)]
p4 <- ggplot(swd, aes(date)) +
  geom_line(aes(y = mld, colour = "MLD (m)"), linewidth = 0.8) +
  geom_line(aes(y = dsigma, colour = "Surface-to-bottom Δσt"), linewidth = 0.8) +
  scale_colour_manual(values = c("MLD (m)" = "steelblue", "Surface-to-bottom Δσt" = "firebrick"), name = NULL) +
  theme_paper() +
  labs(x = NULL, y = "Daily mean", title = "Stratification through the season (daily mean)")

ggsave(file.path(cfg$out_dir, "06_buoy_temp_section.png"), p1, width = 9, height = 3.6, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "06_buoy_sal_section.png"),  p2, width = 9, height = 3.6, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "06_buoy_density_section.png"), p3, width = 9, height = 3.6, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "06_stratification_timeseries.png"), p4, width = 9, height = 3.6, dpi = FIG_DPI)
# ---- (e) the sampled column by month --------------------------------------------
# monthly medians over prof_hourly; hourly profiles drawn faintly because October has
# two hydrographic states that a quantile band would blur
pm <- copy(prof_hourly)[, month := factor(format(time, "%b"), levels = c("Aug", "Sep", "Oct"))]
pm_long <- melt(pm[, .(time, depth, month, `Salinity` = sal, `Temperature (°C)` = temp,
                       `Density anomaly σt (kg m⁻³)` = sigt)],
                id.vars = c("time", "depth", "month"), variable.name = "var", value.name = "val")
pm_med <- pm_long[, .(med = median(val, na.rm = TRUE)), by = .(var, month, depth)]
setorder(pm_med, var, month, depth); setorder(pm_long, var, month, time, depth)
p5 <- ggplot(pm_med, aes(x = med, y = depth, colour = month, shape = month)) +
  geom_line(data = pm_long, aes(x = val, y = depth, group = interaction(time, month)),
            alpha = 0.02, linewidth = 0.3, orientation = "y") +
  geom_line(orientation = "y", linewidth = 0.9) + geom_point(size = 1.6) +
  facet_wrap(~ var, scales = "free_x") +
  scale_y_reverse(breaks = ifcb_depths) +
  scale_colour_manual(values = c(Aug = "#E69F00", Sep = "#0072B2", Oct = "#009E73"), name = NULL) +
  scale_shape_manual(values = c(Aug = 16, Sep = 17, Oct = 15), name = NULL) +
  theme_paper() + theme(legend.position = "top") +
  labs(x = NULL, y = "Depth (m)", title = "The sampled column by month")
ggsave(file.path(cfg$out_dir, "06_hydrographic_profile.png"), p5, width = 11, height = 4.2, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "06_buoy_panel.png"),
       ((p1 + labs(tag = "a")) / (p2 + labs(tag = "b")) / (p3 + labs(tag = "c")) / (p4 + labs(tag = "d")) /
        (p5 + labs(tag = "e"))) + plot_layout(heights = c(1, 1, 1, 1, 1.5)) &
         theme(plot.tag = element_text(face = "bold", size = 15)),
       width = 10, height = 16, dpi = FIG_DPI)

cat("\n========================= ENV-BUOY SUMMARY =========================\n")
cat(sprintf("buoy 33027: %d/%d samples matched (<=%d min).\n",
            nrow(env), nrow(meta), cfg$match_tol_min))
cat(sprintf("table -> 06_env_matched.csv | sections -> %s/06_buoy_panel.png (plus the individual 06_*.png sections)\n", cfg$out_dir))
cat("======================================================================\n")
