# =============================================================================
# Chain length vs hydrography: GAMMs of salinity, stratification, temperature and
# depth, the mixed-layer contrast, depth slopes within salinity bins, and Fig. 6.
#
# Reads:  outputs/tangesund/06_env_matched.csv, 03_chain_metrics.csv, 04_size_spectra.csv,
#         05_composition.csv, 06_buoy_profile_at_ifcb_depths.csv; optional, for Fig. 6d:
#         09b_svea_contrasts.csv, 09b_svea_pertaxon_marine.csv, 25_log_contrast.csv,
#         27c_log_contrast_boot.csv, 06_column_descriptives_by_depth.csv
# Writes: outputs/tangesund/10_aic_comparison.csv, 10_depth_in_hydro_models.csv,
#         10_joint_model_terms.csv, 10_depth_salinity_medians.csv, 10_mixed_layer_contrast.csv,
#         10_mixed_layer_by_depth.csv, 10_salinity_bin_depth_slopes.csv,
#         10_freshest_bin_depth_means.csv, 10_conditioned_slope.csv, the 10_*.png panels
#         and 10_mechanism_panel.png (Fig. 6)
# Deps:   data.table, mgcv, ggplot2, patchwork, ggtext (Fig. 6d)
# Run:    Rscript R/10_env_mechanism.R   (rerun after 09b and 27c to draw Fig. 6d)
# =============================================================================
Sys.setlocale("LC_TIME", "C")  # English month names on date axes
suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(patchwork) })
source("R/lib_plotting_theme.R")

cfg <- list(out_dir = "outputs/tangesund", depths = c(3,6,8,11,13,16))
env  <- fread("outputs/tangesund/06_env_matched.csv")
chn  <- fread("outputs/tangesund/03_chain_metrics.csv")[, .(sample, qc_pass, reliable_chain,
             mean_chain, frac_multi, frac_cells_chain)]
spc  <- fread("outputs/tangesund/04_size_spectra.csv")[, .(sample, reliable, b, size_div, median_esd)]
cmp  <- fread("outputs/tangesund/05_composition.csv")[, .(sample, richness, shannon,
             f_chain_diatom = `f_chain diatom`)]

dat <- Reduce(function(a,b) merge(a,b,by="sample",all.x=TRUE), list(env, chn, spc, cmp))
dat[, `:=`(day = factor(date), depth_f = factor(depth, levels = cfg$depths),
           in_ml = factor(fifelse(depth <= mld, "In mixed layer", "Below pycnocline"),
                          levels = c("Below pycnocline","In mixed layer")))]
ana <- dat[qc_pass == TRUE & reliable_chain == TRUE & is.finite(insitu_temp) & is.finite(mld)]
cat(sprintf("mechanism samples: %d | in-situ temp %.1f-%.1f C, sal %.1f-%.1f, MLD median %.1f m\n",
            nrow(ana), min(ana$insitu_temp), max(ana$insitu_temp),
            min(ana$insitu_sal), max(ana$insitu_sal), median(ana$mld)))

# ---- axis decomposition -----------------------------------------------------
cat(sprintf("axis decomposition: cor(depth,sal)=%.2f  cor(depth,temp)=%.2f  cor(yday,temp)=%.2f\n",
            cor(ana$depth, ana$insitu_sal), cor(ana$depth, ana$insitu_temp), cor(ana$yday, ana$insitu_temp)))

# ---- (A) hydrography GAMM ----------------------------------------------------
# every model has a seasonal interaction with its vertical predictor; no shrinkage (select = FALSE)
mA <- bam(mean_chain ~ s(insitu_temp, k=6) + s(insitu_sal, k=6) + s(dsigma, k=6) +
            s(mld, k=6) + s(yday, k=8) + ti(insitu_sal, yday, k=c(5,6)) + s(day, bs="re"),
          data = ana, method = "fREML")
ssA <- summary(mA); stA <- ssA$s.table
cat(sprintf("\n(A) chain ~ hydrography: dev.expl %.1f%%\n", 100*ssA$dev.expl))
for (tm in c("s(insitu_temp)","s(insitu_sal)","s(dsigma)","s(mld)","s(yday)"))
  cat(sprintf("  %-16s edf=%.2f F=%.1f p=%.3g\n", tm, stA[tm,"edf"], stA[tm,"F"], stA[tm,"p-value"]))

# AIC: depth vs hydrography
mD <- bam(mean_chain ~ s(depth, k=5) + s(yday, k=8) + ti(depth, yday, k=c(5,6)) + s(day, bs="re"),
          data = ana, method = "fREML")
mT <- bam(mean_chain ~ s(insitu_temp, k=6) + s(insitu_sal, k=6) + s(yday, k=8) +
            ti(insitu_sal, yday, k=c(5,6)) + s(day, bs="re"), data = ana, method = "fREML")
cat(sprintf("model comparison (AIC): depth %.1f | T+S %.1f | full hydrography %.1f\n",
            AIC(mD), AIC(mT), AIC(mA)))

# temperature in both models: mDT (depth + T) vs mT (salinity + T)
mDT <- bam(mean_chain ~ s(depth, k=5) + s(insitu_temp, k=6) + s(yday, k=8) +
             ti(depth, yday, k=c(5,6)) + s(day, bs="re"), data = ana, method = "fREML")
cat(sprintf("  temperature held in both: depth+T %.1f | salinity+T %.1f\n", AIC(mDT), AIC(mT)))
cat(sprintf("  -> adding temperature to depth: %.0f AIC points; swapping depth for salinity: %.0f; full hydrography vs salinity+T: %.0f\n",
            AIC(mD) - AIC(mDT), AIC(mDT) - AIC(mT), AIC(mT) - AIC(mA)))
# joint depth + salinity + temperature model (concurvity in 27_model_diagnostics.R)
mDST <- bam(mean_chain ~ s(depth, k=5) + s(insitu_sal, k=6) + s(insitu_temp, k=6) + s(yday, k=8) +
              ti(insitu_sal, yday, k=c(5,6)) + s(day, bs="re"), data = ana, method = "fREML")
stJ <- summary(mDST)$s.table
cat(sprintf("  joint depth + salinity + temperature: AIC %.1f (salinity+T %.1f); s(depth) edf=%.2f F=%.1f p=%.3g | s(insitu_sal) edf=%.2f F=%.1f p=%.3g\n",
            AIC(mDST), AIC(mT), stJ["s(depth)","edf"], stJ["s(depth)","F"], stJ["s(depth)","p-value"],
            stJ["s(insitu_sal)","edf"], stJ["s(insitu_sal)","F"], stJ["s(insitu_sal)","p-value"]))
# full hydrography model + s(depth), with and without ti(depth, yday)
mAD <- bam(mean_chain ~ s(depth, k=5) + s(insitu_temp, k=6) + s(insitu_sal, k=6) + s(dsigma, k=6) + s(mld, k=6) +
             s(yday, k=8) + ti(insitu_sal, yday, k=c(5,6)) + s(day, bs="re"), data = ana, method = "fREML")
mAD2 <- bam(mean_chain ~ s(depth, k=5) + s(insitu_temp, k=6) + s(insitu_sal, k=6) + s(dsigma, k=6) + s(mld, k=6) +
              s(yday, k=8) + ti(insitu_sal, yday, k=c(5,6)) + ti(depth, yday, k=c(5,6)) + s(day, bs="re"),
            data = ana, method = "fREML")
fwrite(data.table(model = c("depth","salinity + temperature","full hydrography","depth + temperature",
                            "depth + salinity + temperature (joint)", "full hydrography + depth",
                            "full hydrography + depth + ti(depth, yday)"),
                  aic = c(AIC(mD), AIC(mT), AIC(mA), AIC(mDT), AIC(mDST), AIC(mAD), AIC(mAD2)),
                  dev_expl = c(summary(mD)$dev.expl, summary(mT)$dev.expl, summary(mA)$dev.expl, summary(mDT)$dev.expl,
                               summary(mDST)$dev.expl, summary(mAD)$dev.expl, summary(mAD2)$dev.expl)),
       file.path(cfg$out_dir, "10_aic_comparison.csv"))
# s(depth) 3 m - 16 m contrast at the median of the other predictors
depth_in <- rbindlist(lapply(list(list("depth + salinity + temperature (joint)", mDST),
                                  list("full hydrography + depth", mAD),
                                  list("full hydrography + depth + ti(depth, yday)", mAD2)), function(x) {
  m <- x[[2]]; st <- summary(m)$s.table
  g <- data.frame(depth = c(3, 16), insitu_sal = median(ana$insitu_sal), insitu_temp = median(ana$insitu_temp),
                  dsigma = median(ana$dsigma), mld = median(ana$mld), yday = median(ana$yday), day = ana$day[1])
  Xp <- predict(m, g, type = "lpmatrix"); cv <- Xp[1, ] - Xp[2, ]
  est <- sum(cv * coef(m)); se <- sqrt(drop(cv %*% vcov(m) %*% cv))
  data.table(model = x[[1]], aic = AIC(m), dev_expl = summary(m)$dev.expl, n = nrow(ana), n_days = uniqueN(ana$day),
             term = "s(depth)", edf = st["s(depth)", "edf"], ref_df = st["s(depth)", "Ref.df"],
             F = st["s(depth)", "F"], p = st["s(depth)", "p-value"],
             contrast_3_16_at_median_hydrography = est, se = se, lo = est - 2 * se, hi = est + 2 * se)
}))
fwrite(depth_in, file.path(cfg$out_dir, "10_depth_in_hydro_models.csv"))
cat("  s(depth) inside the hydrography models (10_depth_in_hydro_models.csv):\n")
print(depth_in[, .(model, aic = round(aic, 1), edf = round(edf, 2), F = round(F, 1), p = signif(p, 2),
                   contrast = sprintf("%+.3f [%+.3f, %+.3f]", contrast_3_16_at_median_hydrography, lo, hi))])
fwrite(data.table(term = rownames(stJ), edf = stJ[, "edf"], ref_df = stJ[, "Ref.df"], F = stJ[, "F"], p = stJ[, "p-value"],
                  n = nrow(ana), n_days = uniqueN(ana$day), dev_expl = summary(mDST)$dev.expl, aic = AIC(mDST)),
       file.path(cfg$out_dir, "10_joint_model_terms.csv"))
# median matched salinity per depth (the 3 m to 16 m span)
sal_depth <- ana[, .(n = .N, sal_median = median(insitu_sal), sal_q25 = quantile(insitu_sal, .25),
                     sal_q75 = quantile(insitu_sal, .75)), by = depth][order(depth)]
sal_depth[, span_3_to_16 := sal_median[depth == 16] - sal_median[depth == 3]]
fwrite(sal_depth, file.path(cfg$out_dir, "10_depth_salinity_medians.csv"))
cat(sprintf("  median matched salinity 3 m %.2f, 16 m %.2f: span %.2f\n", sal_depth[depth == 3, sal_median],
            sal_depth[depth == 16, sal_median], sal_depth[1, span_3_to_16]))

# ---- (B) mixed-layer contrast (day RE) ---------------------------------------
mB <- bam(mean_chain ~ in_ml + s(insitu_temp, k=6) + s(yday, k=8) + s(day, bs="re"),
          data = ana, method = "fREML")
cB <- summary(mB)$p.table
mlb_est <- cB["in_mlIn mixed layer","Estimate"]; mlb_se <- cB["in_mlIn mixed layer","Std. Error"]
mlb_p <- cB["in_mlIn mixed layer","Pr(>|t|)"]
cat(sprintf("(B) chain length IN mixed layer vs below: %+.2f cells (SE %.2f, 95%% CI [%.2f, %.2f], p=%.3g)\n",
            mlb_est, mlb_se, mlb_est - 1.96*mlb_se, mlb_est + 1.96*mlb_se, mlb_p))

# refits with s(depth), and within 3 m samples only (in_ml is nearly a depth indicator)
mB_d  <- bam(mean_chain ~ in_ml + s(depth, k=5) + s(insitu_temp, k=6) + s(yday, k=8) +
               s(day, bs="re"), data = ana, method = "fREML")
cBd <- summary(mB_d)$p.table["in_mlIn mixed layer", ]
mB_s <- bam(mean_chain ~ in_ml + s(insitu_temp, k=6) + s(yday, k=8) + s(day, bs="re"),
            data = ana[depth == 3], method = "fREML")
cBs <- summary(mB_s)$p.table["in_mlIn mixed layer", ]
cat(sprintf("    + s(depth): %+.2f cells (SE %.2f, p=%.2g) | within 3 m only: %+.2f (SE %.2f, p=%.2g)\n",
            cBd["Estimate"], cBd["Std. Error"], cBd["Pr(>|t|)"],
            cBs["Estimate"], cBs["Std. Error"], cBs["Pr(>|t|)"]))
fwrite(data.table(model = c("no depth term (as plotted)", "+ s(depth)", "within 3 m only"),
                  estimate = c(mlb_est, cBd["Estimate"], cBs["Estimate"]),
                  se = c(mlb_se, cBd["Std. Error"], cBs["Std. Error"]),
                  p  = c(mlb_p,  cBd["Pr(>|t|)"],    cBs["Pr(>|t|)"])),
       file.path(cfg$out_dir, "10_mixed_layer_contrast.csv"))

# mixed-layer sample counts per depth
ml_depth <- dcast(ana[, .N, by = .(depth, in_ml)], depth ~ in_ml, value.var = "N", fill = 0)
setnames(ml_depth, c("depth","below_pycnocline","in_mixed_layer"))
ml_depth[, pct_in_ml := 100 * in_mixed_layer / (in_mixed_layer + below_pycnocline)]
fwrite(ml_depth, file.path(cfg$out_dir, "10_mixed_layer_by_depth.csv"))
cat(sprintf("    mixed-layer samples: %d of %d, %d of them at 3 m, %d below 11 m, on %d of %d days\n",
            sum(ml_depth$in_mixed_layer), nrow(ana), ml_depth[depth == 3, in_mixed_layer],
            ml_depth[depth > 11, sum(in_mixed_layer)],
            uniqueN(ana[in_ml == "In mixed layer", day]), uniqueN(ana$day)))

# ---- partial effects for figures --------------------------------------------
pterm <- function(m, var, rng) {
  g <- data.table(insitu_temp = median(ana$insitu_temp), insitu_sal = median(ana$insitu_sal),
                  dsigma = median(ana$dsigma), mld = median(ana$mld), depth = median(ana$depth),
                  yday = median(ana$yday), day = ana$day[1])
  g <- g[rep(1, 80)]; g[[var]] <- seq(rng[1], rng[2], length = 80)
  p <- predict(m, g, se.fit = TRUE, exclude = "s(day)")
  data.table(x = g[[var]], fit = p$fit, se = p$se.fit, term = var)
}
q <- function(v) quantile(ana[[v]], c(.02,.98), na.rm = TRUE)
pe <- rbindlist(list(pterm(mA,"insitu_temp",q("insitu_temp")), pterm(mA,"dsigma",q("dsigma")),
                     pterm(mA,"mld",q("mld"))))
lab <- c(insitu_temp="In-situ temperature (°C)", dsigma="Surface-to-bottom Δσt (kg/m³)", mld="Mixed-layer depth (m)")
pe[, term := factor(lab[term], levels = lab)]

p1_partial <- ggplot(pe, aes(x, fit)) + geom_ribbon(aes(ymin=fit-se, ymax=fit+se), alpha=0.2) +
  geom_line(colour="steelblue", linewidth=1) + facet_wrap(~term, scales="free_x") +
  theme_paper() + labs(x=NULL, y="Mean chain length (partial)",
    title="Chain length vs hydrography (day/season-adjusted)")

# ---- monthly column profiles (10_hydrographic_profile.png) --------------------
# From the hourly buoy profiles (06), not the IFCB-matched rows, which sample depths unequally.
# Keep the setorder calls: geom_line follows row order.
prof <- fread(file.path(cfg$out_dir, "06_buoy_profile_at_ifcb_depths.csv"))
prof[, month := factor(format(time, "%b"), levels = c("Aug","Sep","Oct"))]
prof_long <- melt(prof[, .(time, depth, month, `Salinity` = sal,
                           `Temperature (°C)` = temp, `σt (kg m⁻³)` = sigt)],
                  id.vars = c("time","depth","month"), variable.name = "var", value.name = "val")
prof_med <- prof_long[, .(med = median(val, na.rm = TRUE)), by = .(var, month, depth)]
setorder(prof_med, var, month, depth)

# spread drawn as faint hourly profiles (October temperature is bimodal)
setorder(prof_long, var, month, time, depth)

p1 <- ggplot(prof_med, aes(x = med, y = depth, colour = month)) +
  geom_line(data = prof_long, aes(x = val, y = depth, group = interaction(time, month)),
            alpha = 0.02, linewidth = 0.3, orientation = "y") +
  geom_line(orientation = "y", linewidth = 0.9) + geom_point(size = 1.4) +
  facet_wrap(~ var, scales = "free_x") +
  scale_y_reverse(breaks = cfg$depths) +
  # Okabe-Ito
  scale_colour_manual(values = c(Aug = "#E69F00", Sep = "#0072B2", Oct = "#009E73"), name = NULL) +
  scale_fill_manual(values   = c(Aug = "#E69F00", Sep = "#0072B2", Oct = "#009E73"), guide = "none") +
  theme_paper() + theme(legend.position = "top") +
  labs(x = NULL, y = "Depth (m)",
       title = "The sampled column: haline-stratified, near-isothermal")
# Fig. 6a: chain length vs salinity; pooled GAM smooth (±SE) and per-month smooths
ana[, month := factor(month.abb[as.integer(format(as.Date(date), "%m"))], levels = c("Aug", "Sep", "Oct"))]
p2 <- ggplot(ana, aes(insitu_sal, mean_chain)) +
  geom_point(aes(colour = month), size=0.5, alpha=0.4) +
  geom_smooth(aes(colour = month), method="gam", formula=y~s(x, k = 5), se=FALSE, linewidth = 0.6) +
  geom_smooth(method="gam", formula=y~s(x), colour="firebrick", fill="firebrick", se=TRUE, level=0.68, alpha=0.2) +
  scale_colour_manual(values = c(Aug = "#E69F00", Sep = "#0072B2", Oct = "#009E73"), name = "Month") +
  theme_paper() + theme(legend.position = "top") +
  labs(x="In situ salinity", y="Mean chain length (cells)",
    title="Fresher upper column has longer chains")
# Fig. 6b: chain length vs temperature
p2b <- ggplot(ana, aes(insitu_temp, mean_chain)) +
  geom_point(aes(colour=yday), size=0.5, alpha=0.35) + scale_colour_viridis_c(option="C", name="Day of year") +
  geom_smooth(method="gam", formula=y~s(x), colour="firebrick", fill="firebrick", se=TRUE, level=0.68, alpha=0.2) +
  theme_paper() + labs(x="In situ temperature (°C)", y="Mean chain length (cells)",
    title="Chains shorten as the season cools")
# mixed-layer panel, annotated with the unadjusted and s(depth)-adjusted contrasts
ml_n <- ana[, .N, by = in_ml]
ml_lab <- sprintf("Unadjusted %+.2f cells (p = %.3f)\nWith s(depth) %+.2f (p = %.2f)",
                  mlb_est, mlb_p, cBd["Estimate"], cBd["Pr(>|t|)"])
p3 <- ggplot(ana, aes(in_ml, mean_chain, fill=in_ml)) +
  geom_violin(alpha=0.5, colour=NA) + geom_boxplot(width=0.15, outlier.size=0.3) +
  geom_text(data = ml_n, aes(x = in_ml, y = -Inf, label = paste0("n = ", N)),
            inherit.aes = FALSE, vjust = -0.6, size = 3, colour = "grey30") +
  annotate("text", x = 1.5, y = Inf, label = ml_lab, vjust = 1.2, size = 3.1, colour = "grey20") +
  scale_fill_manual(values=c("Below pycnocline"="grey70","In mixed layer"="steelblue"), guide="none") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
  theme_paper() + labs(x=NULL, y="Mean chain length", title="Mixed layer vs below pycnocline")
p4 <- ggplot(ana[is.finite(insitu_oxy)], aes(insitu_oxy, f_chain_diatom, colour=depth)) +
  geom_point(size=0.5, alpha=0.35) + scale_colour_viridis_c(name="Depth (m)") +
  geom_smooth(method="gam", formula=y~s(x), colour="firebrick", se=FALSE) +
  theme_paper() + labs(x="In-situ oxygen (ml/l)", y="Chain-diatom share of cells",
    title="Chain-diatom dominance vs oxygen")

# ---- (C) depth slopes within salinity bins (Fig. 6c) --------------------------
# fixed-boundary and equal-n bins; linear depth slope with a day RE per bin
sal_slopes <- function(dat, scheme) {
  dat <- copy(dat)
  if (scheme == "Fixed salinity boundaries") {
    dat[, sbin := cut(insitu_sal, c(-Inf, 26, 28, 29, 30, 31, Inf),
                      labels = c("≤ 26","26–28","28–29","29–30","30–31","> 31"))]
  } else {
    q <- unique(quantile(dat$insitu_sal, seq(0, 1, length = 7)))
    dat[, sbin := cut(insitu_sal, q, include.lowest = TRUE,
                      labels = sprintf("%.1f–%.1f", head(q,-1), tail(q,-1)))]
  }
  # not `b`: dat has a column b (size-spectrum slope) that data.table would match instead
  rbindlist(lapply(levels(dat$sbin), function(bin_lab) {
    s <- dat[sbin == bin_lab]
    if (nrow(s) < 60 || uniqueN(s$depth) < 3) return(NULL)
    m <- tryCatch(bam(mean_chain ~ depth + s(day, bs = "re"), data = s, method = "fREML"),
                  error = function(e) { message("salinity bin ", bin_lab, ": ", conditionMessage(e)); NULL })
    if (is.null(m)) return(NULL)
    pt <- summary(m)$p.table["depth", ]
    data.table(scheme = scheme, bin = bin_lab, sal_lo = min(s$insitu_sal),
               n = nrow(s), n_days = uniqueN(s$day),
               depth_sd = sd(s$depth), n_depths = uniqueN(s$depth),
               # samples at the deep end of the bin
               n_deep = sum(s$depth >= 13),
               slope = pt[1], se = pt[2], p = pt[4],
               lo = pt[1] - 1.96*pt[2], hi = pt[1] + 1.96*pt[2])
  }))
}
m_uncond <- bam(mean_chain ~ depth + s(day, bs = "re"), data = ana, method = "fREML")
slope_uncond <- summary(m_uncond)$p.table["depth", "Estimate"]
salb <- rbindlist(lapply(c("Fixed salinity boundaries", "Equal-n salinity bins"),
                         function(sc) sal_slopes(ana, sc)))
salb[, `:=`(scheme = factor(scheme, levels = c("Fixed salinity boundaries","Equal-n salinity bins")),
            retains = p < 0.05)]
setorder(salb, scheme, -sal_lo)                      # freshest bin at the top of each facet
salb[, bin := factor(bin, levels = unique(bin))]
fwrite(salb, file.path(cfg$out_dir, "10_salinity_bin_depth_slopes.csv"))
cat(sprintf("\n(C) depth slope on chain length, unconditioned: %+.3f cells/m\n", slope_uncond))
cat("    within salinity bins (day RE), both binning schemes:\n")
print(salb[, .(scheme, bin, n, n_days, n_deep, depth_sd = round(depth_sd,1),
               slope = round(slope,3), CI = sprintf("[%+.3f, %+.3f]", lo, hi),
               p = signif(p,3), retains)])
# depth SD (full vs within bin) and per-depth means in the freshest bin
cat(sprintf("    depth SD: full data %.1f m | within-bin %.1f-%.1f m\n",
            sd(ana$depth), min(salb$depth_sd), max(salb$depth_sd)))
fresh_bin <- salb[scheme == "Fixed salinity boundaries"][which.min(sal_lo), bin]
fresh <- ana[insitu_sal <= 26, .(n = .N, mean_chain = mean(mean_chain)), by = depth][order(depth)]
fwrite(fresh, file.path(cfg$out_dir, "10_freshest_bin_depth_means.csv"))
cat(sprintf("    freshest bin (%s) mean chain length by depth: %s\n", fresh_bin,
            paste(sprintf("%dm %.2f (n=%d)", fresh$depth, fresh$mean_chain, fresh$n), collapse = ", ")))
cat(sprintf("    slope attenuates from %+.3f to a within-bin range of %+.3f to %+.3f cells/m\n",
            slope_uncond, min(salb$slope), max(salb$slope)))

# linear depth slope with s(salinity) held, all months and per month (Supplementary Note 2)
cond_row <- function(m, label, subset_lab, dat) {
  pt <- summary(m)$p.table["depth", ]
  data.table(subset = subset_lab, model = label, n = nrow(dat), n_days = uniqueN(dat$day),
             slope = pt[1], se = pt[2], lo = pt[1] - 1.96 * pt[2], hi = pt[1] + 1.96 * pt[2], t = pt[3], p = pt[4])
}
cond <- rbindlist(list(
  cond_row(m_uncond, "depth + day RE (unconditioned)", "all months", ana),
  cond_row(bam(mean_chain ~ depth + s(insitu_sal, k=6) + s(day, bs="re"), data = ana, method = "fREML"),
           "depth + s(salinity) + day RE", "all months", ana),
  cond_row(bam(mean_chain ~ depth + s(insitu_sal, k=6) + s(yday, k=8) + s(day, bs="re"), data = ana, method = "fREML"),
           "depth + s(salinity) + s(yday) + day RE", "all months", ana),
  cond_row(bam(mean_chain ~ depth + s(insitu_sal, k=6) + s(yday, k=8) + ti(insitu_sal, yday, k=c(5,6)) + s(day, bs="re"),
               data = ana, method = "fREML"),
           "depth + s(salinity) + s(yday) + ti(salinity, yday) + day RE", "all months", ana)))
ana[, month := factor(month.abb[as.integer(format(as.Date(date), "%m"))], levels = c("Aug", "Sep", "Oct"))]
cond <- rbind(cond, rbindlist(lapply(levels(ana$month), function(mo) {
  s <- ana[month == mo]
  rbind(cond_row(bam(mean_chain ~ depth + s(day, bs="re"), data = s, method = "fREML"), "depth + day RE (unconditioned)", mo, s),
        cond_row(bam(mean_chain ~ depth + s(insitu_sal, k=6) + s(day, bs="re"), data = s, method = "fREML"),
                 "depth + s(salinity) + day RE", mo, s))
})))
cond[, share_of_unconditioned := slope / slope[model == "depth + day RE (unconditioned)"][1], by = subset]
fwrite(cond, file.path(cfg$out_dir, "10_conditioned_slope.csv"))
cat("    depth slope with matched salinity held as a smooth (10_conditioned_slope.csv):\n")
print(cond[, .(subset, model, n, n_days, slope = round(slope, 4), CI = sprintf("[%+.4f, %+.4f]", lo, hi),
               p = signif(p, 3), share = round(share_of_unconditioned, 2))])

p5 <- ggplot(salb, aes(slope, bin, colour = retains, shape = retains)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_vline(xintercept = slope_uncond, linetype = 1, colour = "grey25") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.25) +
  geom_point(size = 2.4, fill = "white") +
  facet_wrap(~ scheme, ncol = 1, scales = "free_y") +
  # filled orange = CI excludes zero, open blue = includes it
  scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  theme_paper() +
  labs(x = "Depth slope within salinity bin (cells per m, ±95% CI)",
       y = "In-situ salinity bin",
       title = "Depth effect with salinity held")

ggsave(file.path(cfg$out_dir,"10_chain_vs_hydrography.png"), p1_partial, width=11, height=3.8, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir,"10_hydrographic_profile.png"), p1, width=11, height=4.2, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir,"10_depth_within_salinity_bins.png"), p5, width=6, height=6, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir,"10_chain_vs_salinity_vertical.png"), p2,  width=7, height=4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir,"10_chain_vs_temp_seasonal.png"),    p2b, width=7, height=4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir,"10_chain_mixedlayer.png"),     p3, width=6,  height=4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir,"10_chaindiatom_vs_oxygen.png"),p4, width=7,  height=4.5, dpi = FIG_DPI)
# ---- (E) regional surface record (Fig. 6d, same as Fig. S5b) --------------------
# Read from 09b's CSVs, not refitted; skipped if absent. Tangesund reference from
# 25_log_contrast.csv (fallback log(2.68/4.27)), interval band from 27c.
svea_csv <- file.path(cfg$out_dir, "09b_svea_contrasts.csv")
tang_csv <- file.path(cfg$out_dir, "25_log_contrast.csv")
tang_boot <- file.path(cfg$out_dir, "27c_log_contrast_boot.csv")
if (file.exists(tang_csv)) {
  tr <- fread(tang_csv)[set == "season-averaged"]
  tang_log <- -tr$contrast; tang_ci <- -c(tr$hi, tr$lo)
  if (file.exists(tang_boot)) { tb <- fread(tang_boot)[set == "season-averaged"]; tang_ci <- -c(tb$hi, tb$lo) }
  else warning("27c_log_contrast_boot.csv missing; the Tangesund band uses the model SE")
} else { tang_log <- log(2.68 / 4.27); tang_ci <- c(NA, NA); warning("25_log_contrast.csv missing; falling back to log(2.68/4.27)") }
p6 <- NULL
if (file.exists(svea_csv) && requireNamespace("ggtext", quietly = TRUE)) {
  sv <- fread(svea_csv)
  # add the Chaetoceros and Pseudo-nitzschia within-taxon rows above Skeletonema
  pt_csv <- file.path(cfg$out_dir, "09b_svea_pertaxon_marine.csv")
  if (file.exists(pt_csv)) {
    ch <- fread(pt_csv)[taxon == "Chaetoceros_spp_chain"]
    ch_row <- data.table(set = "Chaetoceros, within taxon, marine span", lo_sal = 25, hi_sal = 31,
                         contrast = ch$contrast, ci_lo = ch$ci_lo, ci_hi = ch$ci_hi, n = ch$n, cruises = ch$cruises,
                         n_draws = ch$n_draws, flat = NA)
    pn <- fread(pt_csv)[taxon == "Pseudo-nitzschia_spp"]
    pn_row <- data.table(set = "Pseudo-nitzschia, within taxon, marine span", lo_sal = 25, hi_sal = 31,
                         contrast = pn$contrast, ci_lo = pn$ci_lo, ci_hi = pn$ci_hi, n = pn$n, cruises = pn$cruises,
                         n_draws = pn$n_draws, flat = NA)
    sv <- rbind(sv[!grepl("^Skeletonema", set)], ch_row, pn_row, sv[grepl("^Skeletonema", set)], fill = TRUE)
  }
  sv[, label := factor(set, levels = rev(set))]
  sv[, marine := grepl("marine", set)]
  # label span from the buoy median profile (06), else the matched medians
  buoy_csv <- file.path(cfg$out_dir, "06_column_descriptives_by_depth.csv")
  span_lab <- if (file.exists(buoy_csv)) { b <- fread(buoy_csv); c(b[depth == 3, sal_med], b[depth == 16, sal_med]) } else
    c(sal_depth[depth == 3, sal_median], sal_depth[depth == 16, sal_median])
  tang_lab <- sprintf("T\u00e5ngesund 3 m to 16 m\n(salinity %.1f to %.1f)", span_lab[1], span_lab[2])
  p6 <- ggplot(sv, aes(contrast, label, colour = marine)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
    { if (all(is.finite(tang_ci))) annotate("rect", xmin = tang_ci[1], xmax = tang_ci[2], ymin = -Inf, ymax = Inf,
                                            fill = "#d7301f", alpha = 0.08) } +
    geom_vline(xintercept = tang_log, linetype = "dashed", colour = "#d7301f", linewidth = 0.6) +
    annotate("text", x = tang_log + 0.01, y = nrow(sv) + 0.45, label = tang_lab,
             colour = "#d7301f", size = 3, vjust = 1, hjust = 0) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0, linewidth = 0.5) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = c(`TRUE` = "#2c7fb8", `FALSE` = "grey30"), guide = "none") +
    scale_y_discrete(labels = function(x) {
      x <- paste0(toupper(substr(x, 1, 1)), substring(x, 2))
      # label spans by salinity range
      x <- sub("^Full gradient", "Salinity 6 to 18", x)
      x <- sub("^Marine span", "Salinity 25 to 33", x)
      x <- sub("^Skeletonema", "*Skeletonema*", x)
      x <- sub("^Chaetoceros", "*Chaetoceros*", x)
      x <- sub("^Pseudo-nitzschia", "*Pseudo-nitzschia*", x)
      x <- sub("within taxon, marine span", "within taxon, salinity 25 to 31", x)   # state every span
      sub(", ", ",<br>", x)
    }, expand = expansion(add = c(0.6, 1.2))) +
    labs(x = "Change in log chain length over the stated salinity span (\u00b195% CI)", y = NULL,
         title = "Regional surface record: salinity alone") +
    theme_paper() +
    theme(axis.text.y = ggtext::element_markdown())
  cat("panel (d): Svea contrasts read from", svea_csv, "\n")
} else cat("panel (d) skipped: 09b_svea_contrasts.csv missing or ggtext not installed\n")

# Fig. 6: 2 x 2 of p2, p2b, p5, p6
row2 <- if (is.null(p6)) {
  (plot_spacer() | (p5 + labs(tag = "c")) | plot_spacer()) + plot_layout(widths = c(1, 2, 1))
} else {
  (p5 + labs(tag = "c")) | (p6 + labs(tag = "d"))
}
ggsave(file.path(cfg$out_dir,"10_mechanism_panel.png"),
       # p3 and p4 are written standalone only
       ((p2 + labs(tag = "a") | p2b + labs(tag = "b")) / row2) &
         theme(plot.tag = element_text(face = "bold", size = 15)),
       width=15, height=9, dpi = FIG_DPI)

cat("\n=======================================================================\n")
cat(sprintf("table joins 06+03+04+05 | panel -> %s/10_mechanism_panel.png (plus the individual 10_*.png panels)\n", cfg$out_dir))
cat("=======================================================================\n")

# Fig. 6a/b pooled smooths refitted as drawn (geom_smooth gam uses mgcv::gam defaults)
tryCatch(local({
  fs <- gam(mean_chain ~ s(insitu_sal), data = ana)
  gs <- seq(min(ana$insitu_sal), max(ana$insitu_sal), by = 0.1)
  ps <- predict(fs, data.frame(insitu_sal = gs))
  cat(sprintf("\nFig. 6a pooled salinity smooth: peak %.2f cells at salinity %.1f (data %.1f-%.1f)\n",
              max(ps), gs[which.max(ps)], min(gs), max(gs)))
  ft <- gam(mean_chain ~ s(insitu_temp), data = ana)
  gt <- seq(ceiling(min(ana$insitu_temp)), floor(max(ana$insitu_temp)), by = 1)
  pt <- predict(ft, data.frame(insitu_temp = gt), se.fit = TRUE)
  cat("Fig. 6b pooled temperature smooth (cells +- SE) at each whole degree:\n")
  cat(sprintf("  %2d C: %.2f +- %.2f (samples within 0.5 C: %d)\n", gt, pt$fit, pt$se.fit,
              sapply(gt, function(t) sum(abs(ana$insitu_temp - t) < 0.5))), sep = "")
}), error = function(e) message("smooth printout failed: ", conditionMessage(e)))
