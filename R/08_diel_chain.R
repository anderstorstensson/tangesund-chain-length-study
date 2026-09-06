# =============================================================================
# Diel signal in chain length: cyclic smooth of local solar time (primary, Fig. 7),
# solar-elevation smooth, day-night contrast and depth-by-time interaction.
#
# Reads:  outputs/tangesund/03_chain_metrics.csv.
# Writes: outputs/tangesund/08_diel_summary.csv, 08_diel_terms.csv, 08_diel_bin_depth.csv,
#         the 08_*.png panels and 08_diel_panel.png.
# Deps:   data.table, ggplot2, mgcv, patchwork.
# Run:    Rscript R/08_diel_chain.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")  # English month names on date axes
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(mgcv); library(patchwork) })
source("R/lib_plotting_theme.R")

cfg <- list(
  chain_csv = "outputs/tangesund/03_chain_metrics.csv",
  out_dir   = "outputs/tangesund",
  lat = 58.075, lon = 11.493, depths = c(3, 6, 8, 11, 13, 16),
  boot_B = 2000
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# approximate local solar time (h) and solar elevation (deg) from UTC datetime
solar_time <- function(dt, lon) {
  hutc <- as.integer(format(dt, "%H")) + as.integer(format(dt, "%M"))/60
  (hutc + lon/15) %% 24                                   # E positive; no equation of time
}
solar_elev <- function(dt, lat, lon) {
  doy  <- as.integer(format(dt, "%j"))
  decl <- 23.45 * sin(2*pi*(284 + doy)/365) * pi/180
  H    <- (15*(solar_time(dt, lon) - 12)) * pi/180        # hour angle
  asin(sin(lat*pi/180)*sin(decl) + cos(lat*pi/180)*cos(decl)*cos(H)) * 180/pi
}

dat <- fread(cfg$chain_csv)
ts  <- regmatches(dat$sample, regexpr("D[0-9]{8}T[0-9]{6}", dat$sample))
dat[, datetime := as.POSIXct(ts, format = "D%Y%m%dT%H%M%S", tz = "UTC")]
dat[, `:=`(stime = solar_time(datetime, cfg$lon),
           sun   = solar_elev(datetime, cfg$lat, cfg$lon),
           hour  = as.integer(format(datetime, "%H")),
           date  = as.Date(datetime))]
dat[, daynight := factor(fifelse(sun > 0, "day", "night"), levels = c("night","day"))]
ana <- dat[qc_pass == TRUE & reliable_chain == TRUE]
ana[, `:=`(depth_f = factor(depth, levels = cfg$depths), day = factor(date))]
cat(sprintf("diel analysis: %d samples on %d days | solar time %.1f-%.1f h | sun elevation %.0f to %.0f deg\n",
            nrow(ana), uniqueN(ana$day), min(ana$stime), max(ana$stime), min(ana$sun), max(ana$sun)))

# ---- primary: cyclic smooth of solar time -------------------------------------
m_cc <- bam(mean_chain ~ s(stime, bs = "cc", k = 6) + s(depth, k = 5) + s(yday, k = 6) + s(day, bs = "re"),
            knots = list(stime = c(0, 24)), data = ana, method = "fREML")
# depth-by-time interaction
m_int <- bam(mean_chain ~ s(stime, bs = "cc", k = 6) + s(depth, k = 5) +
               ti(stime, depth, bs = c("cc", "tp"), k = c(6, 5)) + s(yday, k = 6) + s(day, bs = "re"),
             knots = list(stime = c(0, 24)), data = ana, method = "fREML")
# secondary: thin-plate smooth of solar elevation
m_el <- bam(mean_chain ~ s(sun, k = 6) + s(depth, k = 5) + s(yday, k = 6) + s(day, bs = "re"),
            data = ana, method = "fREML")

term_row <- function(m, label, term) {
  st <- summary(m)$s.table
  data.table(model = label, term = term, edf = st[term, "edf"], ref_df = st[term, "Ref.df"],
             F = st[term, "F"], p = st[term, "p-value"], n = nrow(ana), n_days = uniqueN(ana$day),
             dev_expl = summary(m)$dev.expl)
}
terms <- rbindlist(list(
  term_row(m_cc,  "cyclic solar time", "s(stime)"),
  term_row(m_cc,  "cyclic solar time", "s(depth)"),
  term_row(m_cc,  "cyclic solar time", "s(yday)"),
  term_row(m_int, "cyclic solar time x depth", "ti(stime,depth)"),
  term_row(m_el,  "solar elevation", "s(sun)"),
  term_row(m_el,  "solar elevation", "s(depth)")))
cat("\n-- smooth terms --\n"); print(terms[, .(model, term, edf = round(edf,2), ref_df = round(ref_df,2),
                                             F = round(F,1), p = signif(p,3))])

# day vs night contrast (mixed, day RE)
m_dn <- bam(mean_chain ~ daynight + s(depth, k = 5) + s(yday, k = 6) + s(day, bs = "re"),
            data = ana, method = "fREML")
co <- summary(m_dn)$p.table
dn_est <- co["daynightday","Estimate"]; dn_se <- co["daynightday","Std. Error"]
dn_lo <- dn_est - 1.96*dn_se; dn_hi <- dn_est + 1.96*dn_se
cat(sprintf("day-night contrast (day - night): %+.3f cells (SE %.3f, t = %.2f, 95%% CI [%+.2f, %+.2f], p=%.3g)\n",
            dn_est, dn_se, co["daynightday","t value"], dn_lo, dn_hi, co["daynightday","Pr(>|t|)"]))
terms <- rbind(terms, data.table(model = "day-night factor", term = "daynight (day - night)",
                                 edf = NA_real_, ref_df = m_dn$df.residual, F = NA_real_,
                                 p = co["daynightday","Pr(>|t|)"], n = nrow(ana), n_days = uniqueN(ana$day),
                                 dev_expl = summary(m_dn)$dev.expl))
terms[, t := c(rep(NA_real_, .N - 1), co["daynightday","t value"])]
fwrite(terms, file.path(cfg$out_dir, "08_diel_terms.csv"))

# ---- (1) diel curve over solar time: amplitude, phase and intervals ----------
g <- data.frame(stime = seq(0, 24, length = 97), depth = median(ana$depth),
                yday = median(ana$yday), day = ana$day[1])
p <- predict(m_cc, g, se.fit = TRUE, exclude = "s(day)")
pr <- data.table(stime = g$stime, fit = p$fit, se = p$se.fit)
cc_amp <- diff(range(pr$fit))
cat(sprintf("cyclic smooth spans %.2f cells over the 24-h cycle (shortest chains at solar time %.1f h, longest at %.1f h)\n",
            cc_amp, pr[which.min(fit), stime], pr[which.max(fit), stime]))
# amplitude interval by simulation from the coefficient posterior
amp_ci <- function(m, g, n_sim = 1000, seed = 1) {
  Xp <- predict(m, g, type = "lpmatrix", exclude = "s(day)")
  set.seed(seed); br <- rmvn(n_sim, coef(m), vcov(m))
  amps <- apply(Xp %*% t(br), 2, function(z) diff(range(z)))
  c(point = diff(range(Xp %*% coef(m))), lo = unname(quantile(amps, .025)), hi = unname(quantile(amps, .975)))
}
cc_ci <- amp_ci(m_cc, g)
cat(sprintf("  amplitude 95%% interval by posterior simulation: [%.2f, %.2f] cells\n", cc_ci["lo"], cc_ci["hi"]))
# null range: s(stime) coefficients set to zero, same covariance; bias-corrected
# amplitude = observed range minus null median
amp_null <- function(m, g, n_sim = 2000, seed = 1) {
  Xp <- predict(m, g, type = "lpmatrix", exclude = "s(day)")
  b0 <- coef(m); b0[grep("s\\(stime\\)", names(b0))] <- 0
  set.seed(seed); br <- rmvn(n_sim, b0, vcov(m))
  amps <- apply(Xp %*% t(br), 2, function(z) diff(range(z)))
  c(median = unname(median(amps)), lo = unname(quantile(amps, .025)), hi = unname(quantile(amps, .975)),
    p_ge_obs = mean(amps >= diff(range(Xp %*% coef(m)))))
}
cc_null <- amp_null(m_cc, g)
cat(sprintf("  range under a flat diel truth (same covariance, 2000 draws): median %.2f, 95%% [%.2f, %.2f]; bias-corrected amplitude %.2f cells\n",
            cc_null["median"], cc_null["lo"], cc_null["hi"], cc_amp - cc_null["median"]))
# sensitivity: primary seasonal terms s(yday, k = 8) + ti(depth, yday)
m_cc8 <- bam(mean_chain ~ s(stime, bs = "cc", k = 6) + s(depth, k = 5) + s(yday, k = 8) +
               ti(depth, yday, k = c(5, 6)) + s(day, bs = "re"),
             knots = list(stime = c(0, 24)), data = ana, method = "fREML")
st8 <- summary(m_cc8)$s.table["s(stime)", ]
cc8_ci <- amp_ci(m_cc8, g)
cat(sprintf("  with s(yday, k = 8) + ti(depth, yday): s(stime) edf %.2f F %.1f p %.3g; amplitude %.2f [%.2f, %.2f]\n",
            st8["edf"], st8["F"], st8["p-value"], cc8_ci["point"], cc8_ci["lo"], cc8_ci["hi"]))
cc8_null <- amp_null(m_cc8, g)
cat(sprintf("  its range under a flat diel truth: median %.2f, 95%% [%.2f, %.2f]; bias-corrected amplitude %.2f cells\n",
            cc8_null["median"], cc8_null["lo"], cc8_null["hi"], cc8_ci["point"] - cc8_null["median"]))
terms <- rbind(terms, data.table(model = "cyclic solar time, primary seasonal terms (k = 8, ti(depth, yday))", term = "s(stime)",
                                 edf = st8["edf"], ref_df = st8["Ref.df"], F = st8["F"], p = st8["p-value"],
                                 n = nrow(ana), n_days = uniqueN(ana$day), dev_expl = summary(m_cc8)$dev.expl, t = NA_real_))
fwrite(terms, file.path(cfg$out_dir, "08_diel_terms.csv"))
# elevation curve
ge <- data.frame(sun = seq(min(ana$sun), max(ana$sun), length = 80),
                 depth = median(ana$depth), yday = median(ana$yday), day = ana$day[1])
pe <- predict(m_el, ge, se.fit = TRUE, exclude = "s(day)")
pre <- data.table(sun = ge$sun, fit = pe$fit, se = pe$se.fit)
el_amp <- diff(range(pre$fit))
cat(sprintf("elevation smooth spans %.2f cells from %.0f to %.0f deg (shortest chains at %.1f deg)\n",
            el_amp, min(pre$sun), max(pre$sun), pre[which.min(fit), sun]))

# night band: solar-time bins with > 50% of samples taken with the sun below the horizon
nb <- ana[, .(night = mean(sun < 0)), by = .(bin = floor(stime * 2) / 2)][order(bin)]
night_lo <- nb[night > 0.5 & bin < 12, max(bin) + 0.5]   # mean sunrise
night_hi <- nb[night > 0.5 & bin > 12, min(bin)]          # mean sunset
cat(sprintf("mean night on the solar-time axis: before %.1f h and after %.1f h\n", night_lo, night_hi))

fwrite(data.table(statistic = c("cyclic smooth amplitude over 24 h (cells)",
                                "cyclic smooth amplitude, primary seasonal terms (k = 8, ti(depth, yday)) (cells)",
                                "solar time of minimum (h)", "solar time of maximum (h)", "day-night contrast (cells)",
                                "elevation smooth range (cells)", "mean night starts after (solar h)", "mean night ends before (solar h)",
                                "range of the fitted diel curve under a flat diel truth (cells)",
                                "cyclic smooth amplitude, bias-corrected (observed minus null median) (cells)",
                                "range under a flat diel truth, primary seasonal terms (cells)",
                                "cyclic smooth amplitude, primary seasonal terms, bias-corrected (cells)"),
                  estimate = c(cc_amp, cc8_ci["point"], pr[which.min(fit), stime], pr[which.max(fit), stime], dn_est, el_amp,
                               night_hi, night_lo,
                               cc_null["median"], cc_amp - cc_null["median"], cc8_null["median"], cc8_ci["point"] - cc8_null["median"]),
                  lo = c(cc_ci["lo"], cc8_ci["lo"], NA, NA, dn_lo, NA, NA, NA, cc_null["lo"], NA, cc8_null["lo"], NA),
                  hi = c(cc_ci["hi"], cc8_ci["hi"], NA, NA, dn_hi, NA, NA, NA, cc_null["hi"], NA, cc8_null["hi"], NA),
                  p  = c(terms[model == "cyclic solar time" & term == "s(stime)", p], st8["p-value"], NA, NA,
                         co["daynightday","Pr(>|t|)"], terms[term == "s(sun)", p], NA, NA,
                         cc_null["p_ge_obs"], NA, cc8_null["p_ge_obs"], NA),
                  note = c("posterior-simulation 95% interval (1000 draws)", "posterior-simulation 95% interval (1000 draws)",
                           "", "", "±1.96 SE", "", "solar-time bins with > 50% of samples at night", "",
                           "s(stime) coefficients set to zero, same covariance, 2000 draws; median and 95% range; p = share of null draws at or above the observed range",
                           "", "as above, for the k = 8 + ti(depth, yday) model", "")),
       file.path(cfg$out_dir, "08_diel_summary.csv"))

p1 <- ggplot(pr, aes(stime, fit)) +
  annotate("rect", xmin = -Inf, xmax = night_lo, ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "grey55") +
  annotate("rect", xmin = night_hi, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "grey55") +
  annotate("text", x = 0.3, y = Inf, hjust = 0, vjust = 1.6,
           label = "night", fontface = "italic", colour = "grey35", size = 3.2) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se), alpha = 0.25, fill = "steelblue") +
  geom_line(colour = "steelblue", linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24), expand = c(0, 0)) +
  theme_paper() +
  labs(x = "Local solar time (h)", y = "Mean chain length (cells; day-adjusted)",
       title = "Diel cycle of diatom chain length")

# ---- (2) raw means by solar-time bin and depth, day-cluster bootstrap SE -----
boot_se_day <- function(value, day, B = cfg$boot_B, seed = 1) {
  days <- unique(day); nd <- length(days)
  if (nd < 2) return(NA_real_)
  idx_by_day <- split(seq_along(value), day)
  set.seed(seed)
  sd(vapply(seq_len(B), function(b) {
    d_samp <- sample(days, nd, replace = TRUE)
    mean(value[unlist(idx_by_day[as.character(d_samp)], use.names = FALSE)])
  }, numeric(1)))
}
ana[, tbin := floor(stime / 3) * 3 + 1.5]                 # eight 3-h bins of solar time
hh <- ana[, .(mean = mean(mean_chain), n = .N, n_day = uniqueN(day),
              se = boot_se_day(mean_chain, as.character(day))), by = .(tbin, depth_f)][order(depth_f, tbin)]
fwrite(hh, file.path(cfg$out_dir, "08_diel_bin_depth.csv"))
p2 <- ggplot(hh, aes(tbin, mean, colour = depth_f, shape = depth_f, linetype = depth_f)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.6, linetype = 1, alpha = 0.6) +
  geom_line() + geom_point(size = 2) +
  scale_colour_viridis_d(name = "Depth (m)") +
  scale_shape_manual(values = c(16, 17, 15, 18, 8, 3), name = "Depth (m)") +
  scale_linetype_manual(values = c(1, 2, 3, 4, 5, 6), name = "Depth (m)") +
  scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24)) + theme_paper() +
  labs(x = "Local solar time (h; 3-h bins)", y = "Mean chain length (cells; ± SE clustered by day)",
       title = "Chain length by time of day and depth")

# ---- (3) diel curve by depth (interaction) -----------------------------------
gg <- CJ(stime = seq(0, 24, length = 97), depth = cfg$depths)
gg[, `:=`(yday = median(ana$yday), day = ana$day[1])]
gg[, fit := predict(m_int, gg, exclude = "s(day)")]
gg[, depth_f := factor(depth, levels = cfg$depths)]
p3 <- ggplot(gg, aes(stime, fit, colour = depth_f)) +
  geom_line(linewidth = 0.8) + scale_colour_viridis_d(name = "Depth (m)") +
  scale_x_continuous(breaks = seq(0, 24, 6)) + theme_paper() +
  labs(x = "Local solar time (h)", y = "Mean chain length (fitted)",
       title = "Diel chain-length cycle by depth")
# ---- (4) chain length vs solar elevation -------------------------------------
p4 <- ggplot(pre, aes(sun, fit)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "grey55") +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se), alpha = 0.25, fill = "steelblue") +
  geom_line(colour = "steelblue", linewidth = 1) + theme_paper() +
  labs(x = "Solar elevation (deg;  <0 = night)", y = "Mean chain length (day-adjusted)",
       title = "Chain length against solar elevation")

ggsave(file.path(cfg$out_dir, "08_diel_chain_curve.png"), p1, width = 7, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "08_chain_by_hour_depth.png"), p2, width = 8, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "08_diel_chain_by_depth.png"), p3, width = 8, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "08_chain_vs_elevation.png"), p4, width = 7, height = 4.5, dpi = FIG_DPI)
# 10 x 4.5 in canvas = 180 mm double-column width
ggsave(file.path(cfg$out_dir, "08_diel_panel.png"),
       ((p1 + labs(tag = "a")) |
        (p2 + labs(tag = "b") + theme(legend.position = "right"))) &
         theme(plot.tag = element_text(face = "bold", size = 13)),
       width = 10, height = 4.5, dpi = FIG_DPI)

cat("\n========================= DIEL SUMMARY =========================\n")
cat(sprintf("tables -> 08_diel_summary.csv, 08_diel_terms.csv, 08_diel_bin_depth.csv | panel -> %s/08_diel_panel.png (plus the individual 08_*.png panels)\n", cfg$out_dir))
cat("==================================================================\n")
