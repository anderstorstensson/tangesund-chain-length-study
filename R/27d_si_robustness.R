# =============================================================================
# The robustness display items of the SI. Nothing is refitted: every number is read
# from CSVs that earlier scripts wrote.
#  Fig. S8: (a) the 3 m minus 16 m chain-length contrast per period under each
#      sensitivity analysis (trait threshold and weighting, 24; within-day AR1 error
#      and day-by-depth means, 26; without study week 5 and on matched days, 17b; the
#      primary contrast with its day-resampled interval, 27c), and season-averaged
#      under the model choices of 25 and on the accurate-to-moderate detectors (19);
#      (b) the within-taxon share under each partition variant (17, 20, 17b);
#      (c) the share of chain-former ROIs whose fluorescence peak lies within 20% above
#      the trigger threshold, by depth, cell number and month (02d).
#  27d_si_tables.md: markdown rows for Table S3 (sensitivity rows and days with both
#      depths), Data S1 (observed concurvity, AR1 hydrography terms, residual
#      autocorrelation), Table S2 with a day-varying-slope column, Table S7 (depth slope
#      within salinity bins with their month composition) and Table S8 (the Svea
#      marine-span contrast under each definition and its share of the Tangesund decline).
#
# Reads:  outputs/tangesund/ 02d_near_threshold, 03_chain_metrics, 09b_svea_contrasts,
#         09b_svea_contrasts_endpoint, 09b_svea_pertaxon_marine, 09c_svea_autumn,
#         10_salinity_bin_depth_slopes, 10d_salinity_bin_months, 17_plasticity_turnover,
#         17_partition_variants, 17b_matched_day_partition, 17b_sep_without_week5,
#         19_counter_robustness, 20_partition_variants, 20_sample_log_contrasts,
#         20_within_taxon_pertaxon, 20_within_taxon_pertaxon_model,
#         24_threshold_contrast_sweep, 25_contrast_conventions, 26_ar1_contrasts,
#         26_ar1_hydrography_terms, 26_ar1_rho, 27_concurvity_full, 27_concurvity_joint,
#         27c_log_contrast_boot, 27c_period_contrast_boot, 27c_pertaxon_model_boot,
#         27c_primary_contrast_boot, 27c_svea_share_of_reference (all .csv).
# Writes: outputs/tangesund/27d_robustness_panel.png/.pdf (Fig. S8), 27d_si_tables.md.
# Deps:   data.table, ggplot2, patchwork.
# Run:    Rscript R/27d_si_robustness.R   (last before 28; needs 27c run with
#         BOOT_PARTS=E and the BOOT_PARTS=B rerun after 09b)
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("R/lib_plotting_theme.R")
cfg <- list(out_dir = "outputs/tangesund", surf = 3, deep = 16,
            span_tangesund = 7,      # salinity units between 3 m and 16 m (Fig. 6 legend)
            fig = "27d_robustness_panel", tables = "27d_si_tables.md")
rd <- function(f) fread(file.path(cfg$out_dir, f))
periods <- c("Aug", "Sep", "Oct")
pname <- c(Aug = "August", Sep = "September", Oct = "October")
f2  <- function(x) sprintf("%.2f", x)
ci2 <- function(e, lo, hi) sprintf("%.2f [%.2f, %.2f]", e, lo, hi)
mrow <- function(...) paste0("| ", paste(c(...), collapse = " | "), " |")
minus <- function(s) gsub("-", "−", s, fixed = TRUE)     # typographic minus for table cells
out <- character(0)
say <- function(...) out <<- c(out, ...)

# ---- (1) Table S3 rows and the per-period contrasts -----------------------------
prim <- rd("27c_period_contrast_boot.csv")[response == "Chain length (cells)",
                                            .(period, est = estimate, lo = boot_lo, hi = boot_hi)]
prim[, `:=`(analysis = "Primary model (day-cluster bootstrap interval)", interval = "bootstrap")]
sweep <- rd("24_threshold_contrast_sweep.csv")[threshold %in% c(0, 5, 10, 50, -1)]
sweep[, analysis := fifelse(threshold == -1, "Unthresholded, weighted by chain-counted ROIs",
                            sprintf("Threshold %s chain-counted ROIs", ifelse(threshold == 0, "0", paste0("≥ ", threshold))))]
sweep <- sweep[, .(period, est, lo, hi, analysis, interval = "model")]
ar1 <- rd("26_ar1_contrasts.csv")[model != "naive (day-RE)"]
ar1[, analysis := fifelse(model == "AR1-corrected", "Within-day AR1 error", "Day-by-depth means, chain-ROI-weighted")]
ar1 <- ar1[, .(period, est, lo, hi, analysis, interval = "model")]
wk <- rd("17b_sep_without_week5.csv")[set != "all weeks"]
wk[, analysis := fifelse(set == "without week 5", "Without study week 5", "Without study weeks 5–7")]
wk <- wk[, .(period, est, lo, hi, analysis, interval = "model")]
md <- rd("17b_matched_day_partition.csv")[period != "All", .(period, est = dCWM, lo = NA_real_, hi = NA_real_,
                                                              analysis = "Days with both depths only (community-weighted mean)", interval = "none")]
fore <- rbind(prim, sweep, ar1, wk, md)
fore[, period := factor(period, levels = periods)]
lev <- c("Primary model (day-cluster bootstrap interval)", "Threshold 0 chain-counted ROIs", "Threshold ≥ 5 chain-counted ROIs",
         "Threshold ≥ 10 chain-counted ROIs", "Threshold ≥ 50 chain-counted ROIs", "Unthresholded, weighted by chain-counted ROIs",
         "Within-day AR1 error", "Day-by-depth means, chain-ROI-weighted", "Without study week 5", "Without study weeks 5–7",
         "Days with both depths only (community-weighted mean)")
fore[, analysis := factor(analysis, levels = lev)]
stopifnot(!anyNA(fore$analysis))
# days on which both a 3 m and a 16 m sample pass the chain-trait threshold
chn <- rd("03_chain_metrics.csv")[qc_pass == TRUE & reliable_chain == TRUE]
chn[, date := as.Date(date)][, period := month.abb[as.integer(format(date, "%m"))]]
days_both <- chn[, .(days_both = sum(tapply(depth, date, function(d) all(c(cfg$surf, cfg$deep) %in% d)))), by = period]
days_both <- days_both[match(periods, period)]

say("## Supplementary Table 3: rows appended below the two primary rows", "")
say(mrow("Response", "August", "September", "October"), "|:--|:--:|:--:|:--:|")
say(mrow("Sampling days with both a 3 m and a 16 m sample above the chain-trait threshold (n)", days_both$days_both))
for (a in lev[-1]) {
  x <- fore[analysis == a][match(periods, period)]
  cell <- if (x$interval[1] == "none") f2(x$est) else ci2(x$est, x$lo, x$hi)
  say(mrow(paste0("Mean chain length, ", tolower(substr(a, 1, 1)), substr(a, 2, nchar(a))), minus(cell)))
}
say("")

# season-averaged strip: model choices (25) and the detector restriction (19)
conv <- rd("25_contrast_conventions.csv")[convention == "averaged over observed days of year"]
pb <- rd("27c_primary_contrast_boot.csv")[quantity == "season_averaged"]
season <- rbind(
  data.table(analysis = "Primary model (day-cluster bootstrap interval)", est = pb$estimate, lo = pb$boot_lo, hi = pb$boot_hi),
  conv[model == "precision-weighted", .(analysis = "Weighted by chain-counted ROIs", est = contrast, lo, hi)],
  conv[grepl("high k", model), .(analysis = "Seasonal smooth k = 20, depth k = 6", est = contrast, lo, hi)],
  conv[model == "Gamma(log)", .(analysis = "Gamma errors, log link", est = contrast, lo, hi)],
  rd("19_counter_robustness.csv")[set == "accurate/moderate only", .(analysis = "Accurate-to-moderate detectors only (difference of means)", est = contrast, lo = NA_real_, hi = NA_real_)])
season[, period := "Season"]
fore_all <- rbind(fore[, .(period = as.character(period), analysis = as.character(analysis), est, lo, hi)], season)
fore_all[, period := factor(period, levels = c(periods, "Season"), labels = c(pname, "Season-averaged"))]
fore_all[, analysis := factor(analysis, levels = rev(c(lev, setdiff(season$analysis, lev))))]
pa <- ggplot(fore_all, aes(est, analysis)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.3, colour = "grey30", na.rm = TRUE) +
  geom_point(aes(shape = is.na(lo)), size = 2.2, colour = "firebrick") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  facet_grid(. ~ period, scales = "free_x") +
  theme_paper() +
  labs(x = "Surface (3 m) minus deep (16 m) chain-length contrast (cells; 95% interval)", y = NULL,
       title = "The vertical contrast under every sensitivity analysis")

# ---- (2) the within-taxon share under each partition variant -------------------
sh <- rd("17_plasticity_turnover.csv")
sh[, variant := fifelse(period == "All", "All detectors, ROI-weighted (primary)", paste0("Primary, ", pname[period]))]
share <- rbind(sh[, .(variant, share = frac_within, lo = frac_within_lo, hi = frac_within_hi)],
               rd("17_partition_variants.csv")[!grepl("primary", variant),
                 .(variant = fifelse(grepl("accurate", variant), "Accurate-to-moderate detectors only", "All detectors, cell-weighted"),
                   share = frac_within, lo = NA_real_, hi = NA_real_)],
               rd("20_partition_variants.csv")[grepl("purified", variant),
                 .(variant = "Chaetoceros class purified (blind count)", share = frac_within, lo = NA_real_, hi = NA_real_)],
               rd("17b_matched_day_partition.csv")[period %in% c("All", "Sep"),
                 .(variant = fifelse(period == "All", "Days with both depths only", "Days with both depths only, September"),
                   share = frac_within, lo = frac_within_lo, hi = frac_within_hi)])
share[, variant := factor(variant, levels = rev(unique(variant)))]
pb_ <- ggplot(share, aes(share, variant)) +
  geom_vline(xintercept = 50, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.3, colour = "grey30", na.rm = TRUE) +
  geom_point(aes(shape = is.na(lo)), size = 2.2, colour = "steelblue4") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  scale_x_continuous(limits = c(0, 100)) +
  theme_paper() +
  labs(x = "Within-taxon share of the contrast (%)", y = NULL,
       title = "The within-taxon share under each partition variant")

# ---- (3) fluorescence near-threshold fraction ---------------------------------
nt <- rd("02d_near_threshold.csv")[group %in% c("single cell", "multi-cell")]
nt[, period := factor(period, levels = c("all", periods), labels = c("All", "Aug", "Sep", "Oct"))]
nt[, group := factor(group, levels = c("single cell", "multi-cell"), labels = c("Single-cell ROIs", "Multi-cell ROIs"))]
pc <- ggplot(nt, aes(pct_near_threshold, depth, colour = group, shape = group)) +
  geom_path(orientation = "y") + geom_point(size = 2) +
  scale_y_reverse(breaks = c(3, 6, 8, 11, 13, 16)) + scale_x_continuous(breaks = c(20, 40)) +
  scale_colour_manual(values = c("firebrick", "steelblue4"), name = NULL) + scale_shape_manual(values = c(16, 17), name = NULL) +
  facet_grid(. ~ period) +
  theme_paper() + theme(legend.position = "top") +
  labs(x = "ROIs within 20% above the threshold (%)", y = "Depth (m)",
       title = "Near-threshold fraction by depth")

bottom <- ((pb_ + labs(tag = "b")) | (pc + labs(tag = "c"))) + plot_layout(widths = c(1, 1.4))
fig <- (pa + labs(tag = "a")) / bottom +
  plot_layout(heights = c(1.15, 1)) & theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(file.path(cfg$out_dir, paste0(cfg$fig, ".png")), fig, width = 13, height = 9.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, paste0(cfg$fig, ".pdf")), fig, width = 13, height = 9.5, device = cairo_pdf)
cat(sprintf("figure -> %s/%s.{png,pdf}\n", cfg$out_dir, cfg$fig))

# ---- (4) Data S1 rows ----------------------------------------------------------
cf <- rd("27_concurvity_full.csv")[measure == "observed"]
cj <- rd("27_concurvity_joint.csv")[measure == "observed"]
ar1t <- rd("26_ar1_hydrography_terms.csv")
rho <- rd("26_ar1_rho.csv")
say("## Supplementary Data 1: rows appended", "")
say(mrow("Block", "Model", "Response", "Unit", "n", "Days", "Dev. expl. (%)", "Term", "edf", "Ref. df", "F", "t", "p", "Estimate [95% CI]"))
say("|:--|:--|:--|:--|--:|--:|--:|:--|--:|--:|--:|--:|--:|:--|")
terms_full <- c("s(insitu_temp)", "s(insitu_sal)", "s(dsigma)", "s(mld)", "ti(insitu_sal,yday)")
first <- TRUE
for (tm in terms_full) {
  say(mrow(if (first) "Halocline as the vertical axis" else "", if (first) "Hydrography GAMM, observed concurvity of each term with the other terms" else "",
           if (first) "mean chain length (cells)" else "", if (first) "sample" else "", if (first) "1,892" else "", if (first) "56" else "", "",
           tm, "", "", "", "", "", sprintf("%.3f", cf[[tm]])))
  first <- FALSE
}
say(mrow("", "Joint GAMM (depth + salinity + temperature), observed concurvity of the depth smooth with the other terms", "mean chain length (cells)", "sample", "1,892", "56", "",
         "s(depth)", "", "", "", "", "", sprintf("%.2f", cj[["s(depth)"]])))
first <- TRUE
fmtp <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
for (i in seq_len(nrow(ar1t))) {
  say(mrow(if (first) "Halocline as the vertical axis" else "", if (first) "Hydrography GAMM refitted with a within-day AR1 error (Supplementary Note 1)" else "",
           if (first) "mean chain length (cells)" else "", if (first) "sample" else "", if (first) "1,892" else "", if (first) "56" else "", "",
           ar1t$term[i], "", "", sprintf("%.1f", ar1t$F_AR1[i]), "", fmtp(ar1t$p_AR1[i]), ""))
  first <- FALSE
}
r1 <- rho[grepl("chain-length", model)]
say(mrow("Vertical decline in chain length", "Chain-length period model, residual autocorrelation within a sampling day", "mean chain length (cells)", "sample", "1,893", "56", "",
         "lag-1 ρ between consecutive samples; between successive same-depth revisits", "", "", "", "", "",
         sprintf("%.2f; %.2f", r1$rho_lag1_within_day, r1$same_depth_revisit_acf)))
say("")

# ---- (5) Table S2 with the day-varying-slope column -----------------------------
pt <- rd("20_within_taxon_pertaxon.csv")[, row := .I]        # the CSV order is the table's order (by contribution)
pm <- rd("20_within_taxon_pertaxon_model.csv")[, .(class_name, n_surf, pct, pct_lo, pct_hi)]
rs <- rd("27c_pertaxon_model_boot.csv")[, .(class_name, pct_rs, pct_lo_rs, pct_hi_rs)]
pt <- merge(pt, pm, by = "n_surf", sort = FALSE); pt <- merge(pt, rs, by = "class_name", sort = FALSE)
setorder(pt, row)
label <- c(Chaetoceros_spp_chain = "*Chaetoceros* spp. (chain class)", Asterionellopsis_glacialis = "*Asterionellopsis glacialis*",
           Cerataulina_pelagica = "*Cerataulina pelagica*", Dinobryon_spp = "*Dinobryon* spp.", `Pseudo-nitzschia_spp` = "*Pseudo-nitzschia* spp.",
           Guinardia_delicatula = "*Guinardia delicatula*", Dactyliosolen_fragilissimus = "*Dactyliosolen fragilissimus*",
           `Chaetoceros_tenuissimus-like` = "*Chaetoceros* cf. *tenuissimus*", Uroglenopsis_sp = "*Uroglenopsis* sp.",
           Leptocylindrus_danicus_Leptocylindrus_minimus = "*Leptocylindrus danicus/minimus*", Chaetoceros_spp_single = "*Chaetoceros* spp. (single-cell class)",
           Thalassiosira_levanderi = "*Thalassiosira levanderi*", `Thalassiosira_gravida-like` = "*Thalassiosira* cf. *gravida*", Skeletonema_spp = "*Skeletonema* spp.")
stopifnot(all(pt$class_name %in% names(label)))
pctci <- function(e, lo, hi) sprintf("%+.0f%% [%+.0f, %+.0f]", e, lo, hi)
sgn <- function(x, d) minus(ifelse(x > 0, sprintf(paste0("+%.", d, "f"), x), ifelse(x < 0, sprintf(paste0("%.", d, "f"), x), sprintf(paste0("%.", d, "f"), x))))
say("## Supplementary Table 2: regenerated (new last-but-two column, day-varying slope)", "")
say(mrow("Taxon", "ROIs 3 m", "ROIs 16 m", "Chain 3 m", "Chain 16 m", "Δ", "Contribution", "% of within-taxon term", "ROI model, 3 m vs 16 m", "ROI model, slope varying by day", "Exact %", "Flag"))
say("|:--|--:|--:|--:|--:|--:|--:|--:|:--:|:--:|--:|:--|")
for (i in seq_len(nrow(pt))) {
  x <- pt[i]
  say(mrow(label[[x$class_name]], x$n_surf, x$n_deep, f2(x$L_surf), f2(x$L_deep), sgn(x$delta, 2),
           minus(sprintf("%.3f", x$within_i)), minus(sprintf("%.1f", x$pct_of_within)),
           pctci(x$pct, x$pct_lo, x$pct_hi), pctci(x$pct_rs, x$pct_lo_rs, x$pct_hi_rs),
           sprintf("%.1f", x$exact_pct), fifelse(x$flag == "ok", "accurate", x$flag)))
}
say("")

# ---- (6) Table S7: depth slope within salinity bins -----------------------------
sl <- rd("10_salinity_bin_depth_slopes.csv")
mo <- rd("10d_salinity_bin_months.csv")
bins <- merge(sl, mo[, .(scheme, bin, pct_Aug, pct_Sep, pct_Oct)], by = c("scheme", "bin"), sort = FALSE)
bins[, scheme := factor(scheme, levels = c("Fixed salinity boundaries", "Equal-n salinity bins"))]
setorder(bins, scheme, -sal_lo)
say("## Supplementary Table 7: depth slope within salinity bins (new table)", "")
say(mrow("Binning", "Salinity bin", "Samples", "Days", "Depths sampled", "Samples at 13 and 16 m", "August (%)", "September (%)", "October (%)", "Depth slope (cells m^−1^) [95% CI]", "Excludes zero"))
say("|:--|:--|--:|--:|--:|--:|--:|--:|--:|:--|:--:|")
for (i in seq_len(nrow(bins))) {
  x <- bins[i]
  say(mrow(as.character(x$scheme), x$bin, x$n, x$n_days, x$n_depths, x$n_deep, sprintf("%.0f", x$pct_Aug), sprintf("%.0f", x$pct_Sep), sprintf("%.0f", x$pct_Oct),
           minus(sprintf("%.3f [%.3f, %.3f]", x$slope, x$lo, x$hi)), fifelse(x$retains, "yes", "no")))
}
say("")

# ---- (7) Table S8: the Svea marine-span contrast under every definition ---------
ref <- rd("27c_log_contrast_boot.csv")[quantity == "season_averaged", .(est = estimate, lo = boot_lo, hi = boot_hi)]
tax_ref <- rd("20_sample_log_contrasts.csv")[, .(class_name, est, lo, hi)]
svea <- rbind(
  rd("09b_svea_contrasts.csv")[grepl("^marine span", set), .(definition = fifelse(grepl("equally", set), "Marine span 25–33, colonies weighted equally", "Marine span 25–33, biovolume-weighted"),
                                                             lo_sal, hi_sal, n, cruises, contrast, ci_lo, ci_hi)],
  rd("09b_svea_contrasts_endpoint.csv")[, .(definition = fifelse(grepl("equally", set), "Marine span ending at 32, colonies weighted equally", "Marine span ending at 32, biovolume-weighted"),
                                            lo_sal, hi_sal, n, cruises, contrast, ci_lo, ci_hi)],
  rd("09c_svea_autumn.csv")[order(months != "August to October"), .(definition = sprintf("%s cruises only, %s", months, weighting), lo_sal, hi_sal, n, cruises, contrast, ci_lo, ci_hi)],
  rd("09b_svea_pertaxon_marine.csv")[taxon %in% c("Chaetoceros_spp_chain", "Pseudo-nitzschia_spp", "Skeletonema_spp")][order(match(taxon, c("Chaetoceros_spp_chain", "Pseudo-nitzschia_spp", "Skeletonema_spp"))),
    .(definition = c(Chaetoceros_spp_chain = "*Chaetoceros* spp., within taxon", `Pseudo-nitzschia_spp` = "*Pseudo-nitzschia* spp., within taxon",
                     Skeletonema_spp = "*Skeletonema* spp., within taxon")[taxon], lo_sal = 25, hi_sal = 31, n, cruises, contrast, ci_lo, ci_hi, class_name = taxon)],
  fill = TRUE)
svea[, span := hi_sal - lo_sal]
svea[, scaled := contrast * cfg$span_tangesund / span]
svea[, `:=`(ref_est = ref$est, ref_lo = ref$lo, ref_hi = ref$hi, reference = "Community, 3 m to 16 m")]
for (i in which(!is.na(svea$class_name))) {
  r <- tax_ref[class_name == svea$class_name[i]]
  svea[i, `:=`(ref_est = r$est, ref_lo = r$lo, ref_hi = r$hi, reference = paste0(sub(", within taxon", "", definition), ", 3 m to 16 m"))]
}
# share of the reference decline: the Svea contrast is negative for a decline, the reference is stored as a positive decline
svea[, share := 100 * (-scaled) / ref_est]
svea[, share_lo := 100 * pmax(0, -ci_hi * cfg$span_tangesund / span) / ref_est]
svea[, share_hi := 100 * (-ci_lo * cfg$span_tangesund / span) / ref_est]
# the two 25-33 shares must reproduce 27c's own
chk <- rd("27c_svea_share_of_reference.csv")
stopifnot(all(abs(svea[1:2, share] - chk$share_point) < 0.05))
say("## Supplementary Table 8: the Svea marine-span contrast under every definition (new table)", "")
say(mrow("Definition", "Salinity span", "Samples", "Cruises", "Contrast (log units) [95% CI]", "Scaled to 7 salinity units", "Tångesund reference decline [95% CI]", "Share of the reference (%)", "Share over the *Svea* interval (%)"))
say("|:--|:--:|--:|--:|:--|--:|:--|--:|:--|")
for (i in seq_len(nrow(svea))) {
  x <- svea[i]
  say(mrow(x$definition, sprintf("%g–%g", x$lo_sal, x$hi_sal), x$n, x$cruises, minus(ci2(x$contrast, x$ci_lo, x$ci_hi)), minus(f2(x$scaled)),
           paste0(x$reference, ": ", minus(ci2(x$ref_est, x$ref_lo, x$ref_hi))), sprintf("%.0f", if (abs(x$share) < 0.5) 0 else x$share),
           sprintf("%.0f to %.0f", x$share_lo, x$share_hi)))
}
say("")
writeLines(out, file.path(cfg$out_dir, cfg$tables))
cat(sprintf("tables -> %s/%s\n", cfg$out_dir, cfg$tables))
