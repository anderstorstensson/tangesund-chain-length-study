# =============================================================================
# Day-resampled (percentile) intervals for the 3 m - 16 m contrasts: (A) primary
# GAMM, (B) its log-scale reference, (C) per-period contrasts (Fig. S2,
# Supplementary Table 3), and the diagnostics (D) Table 2 ROI models, (E) per-taxon ROI models.
#
# Reads:  outputs/tangesund/02_sample_meta.csv, 03_chain_metrics.csv, 05_composition.csv,
#         13_depth_contrasts.csv, 09b_svea_contrasts.csv (optional, part B share rows);
#         data/roi_level/roi_2016.txt.gz for (D) and (E).
# Writes: outputs/tangesund/27c_primary_contrast_boot.csv, 27c_log_contrast_boot.csv,
#         27c_svea_share_of_reference.csv, 27c_period_contrast_boot.csv,
#         27c_roi_model_boot.csv, 27c_pertaxon_model_boot.csv, 27c_contrasts_panel.png (Fig. S2).
# Env:    BOOT_PARTS (default "ABC"; D and E are diagnostics); BOOT_FIG_ONLY=1 redraws
#         Fig. S2 from the CSVs (used by 28); BOOT_THREADS, threads for the (D)/(E) fits.
# Deps:   data.table, mgcv, ggplot2, parallel.
# Run:    Rscript R/27c_day_bootstrap.R
#         Run after 13 and before 09b and 10; rerun with BOOT_PARTS=B after 09b to
#         write 27c_svea_share_of_reference.csv.
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(parallel) })
source("R/lib_plotting_theme.R")

cfg <- list(
  out_dir  = "outputs/tangesund",
  meta     = "outputs/tangesund/02_sample_meta.csv",
  chain    = "outputs/tangesund/03_chain_metrics.csv",
  comp     = "outputs/tangesund/05_composition.csv",
  contrasts = "outputs/tangesund/13_depth_contrasts.csv",
  roi_file = "data/roi_level/roi_2016.txt.gz",
  junk     = c("Debris", "Beads", "Air_bubbles", "Zooplankton"),
  weak_taxa = c("Chaetoceros", "Asterionellopsis", "Pauliella"),   # as 20
  min_taxon_model = 200,                                           # as 20
  surf = 3, deep = 16, min_chain_roi = 20,
  B_sample = 1000, B_roi = 200,
  cores = min(6L, max(1L, detectCores() - 2L)),
  parts = strsplit(Sys.getenv("BOOT_PARTS", "ABC"), "")[[1]],
  fig_only = nzchar(Sys.getenv("BOOT_FIG_ONLY"))
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
RNGkind("L'Ecuyer-CMRG")

# ---- shared machinery ---------------------------------------------------------
# Resample days with replacement (within strata if given) and refit. Days are
# relabeled b1..bn so the day RE stays identifiable. stat_fun(dat) returns a named
# numeric vector; the result is a draws x statistics matrix.
day_boot <- function(dat, stat_fun, B, seed, strata = NULL, cores = cfg$cores) {
  dat <- copy(dat)[, day := as.character(day)]
  days <- if (is.null(strata)) list(unique(dat$day)) else split(unique(dat[, .(day, s = get(strata))]), by = "s")
  idx <- split(seq_len(nrow(dat)), dat$day)
  one <- function(b) {
    pick <- if (is.null(strata)) sample(days[[1]], length(days[[1]]), replace = TRUE)
            else unlist(lapply(days, function(g) sample(g$day, nrow(g), replace = TRUE)), use.names = FALSE)
    xb <- rbindlist(lapply(seq_along(pick), function(i) {
      x <- dat[idx[[pick[i]]]]; x[, day := paste0("b", i)]
      if ("samp" %in% names(x)) x[, samp := paste0(samp, "_b", i)]
      x }))
    xb[, day := factor(day)]; if ("samp" %in% names(xb)) xb[, samp := factor(samp)]
    r <- try(stat_fun(xb), silent = TRUE)
    if (inherits(r, "try-error")) NULL else r
  }
  set.seed(seed)
  draws <- mclapply(seq_len(B), one, mc.cores = cores, mc.set.seed = TRUE)
  ok <- !vapply(draws, is.null, logical(1))
  do.call(rbind, draws[ok])
}
summarise_boot <- function(obs, se, dr, B, label) {
  # obs, se: model estimate and SE (named vectors); dr: draws matrix
  rbindlist(lapply(names(obs), function(q) data.table(
    quantity = q, label = label, estimate = obs[[q]], model_se = se[[q]],
    model_lo = obs[[q]] - 2 * se[[q]], model_hi = obs[[q]] + 2 * se[[q]],
    boot_median = median(dr[, q]), boot_sd = sd(dr[, q]),
    boot_lo = unname(quantile(dr[, q], .025)), boot_hi = unname(quantile(dr[, q], .975)),
    width_ratio = (quantile(dr[, q], .975) - quantile(dr[, q], .025)) / (4 * se[[q]]),
    B = B, n_draws = nrow(dr))))
}

meta <- fread(cfg$meta)[, .(sample, depth, date, yday, qc_pass)]
chn  <- fread(cfg$chain)[, .(sample, n_chain_roi, mean_chain)]
D <- merge(meta, chn, by = "sample")[qc_pass == TRUE & is.finite(mean_chain)]
D[, `:=`(day = factor(date),
         period = factor(month.abb[as.integer(format(as.Date(date), "%m"))], levels = c("Aug", "Sep", "Oct")))]
ana <- D[n_chain_roi >= cfg$min_chain_roi]
cat(sprintf("chain-thresholded set: %d samples on %d days\n", nrow(ana), uniqueN(ana$day)))

# season-averaged and median-day 3 m - 16 m contrasts from the lpmatrix (as in 25)
form_primary <- mean_chain ~ s(depth, k = 5) + s(yday, k = 8) + ti(depth, yday, k = c(5, 6)) + s(day, bs = "re")
med_yday <- median(ana$yday)
conv <- function(m, dat, with_se = FALSE) {
  Xa <- function(dp) colMeans(predict(m, data.frame(depth = dp, yday = dat$yday, day = dat$day[1]), type = "lpmatrix"))
  Xm <- function(dp) predict(m, data.frame(depth = dp, yday = med_yday, day = dat$day[1]), type = "lpmatrix")[1, ]
  ca <- Xa(cfg$surf) - Xa(cfg$deep); cm <- Xm(cfg$surf) - Xm(cfg$deep)
  est <- c(season_averaged = sum(ca * coef(m)), median_day = sum(cm * coef(m)))
  if (!with_se) return(est)
  list(est = est, se = c(season_averaged = sqrt(drop(ca %*% vcov(m) %*% ca)),
                         median_day = sqrt(drop(cm %*% vcov(m) %*% cm))))
}

# ---- (A) primary contrast -----------------------------------------------------
if ("A" %in% cfg$parts && !cfg$fig_only) {
  m <- bam(form_primary, data = ana, method = "fREML")
  obs <- conv(m, ana, with_se = TRUE)
  t0 <- Sys.time()
  dr <- day_boot(ana, function(xb) conv(bam(form_primary, data = xb, method = "fREML"), xb), cfg$B_sample, seed = 1)
  outA <- summarise_boot(obs$est, obs$se, dr, cfg$B_sample, "primary GAMM, mean chain length (cells), 3 m - 16 m")
  outA[, `:=`(n = nrow(ana), n_days = uniqueN(ana$day))]
  fwrite(outA, file.path(cfg$out_dir, "27c_primary_contrast_boot.csv"))
  cat(sprintf("\n== (A) primary contrast, day-cluster bootstrap (B = %d, %.1f min) ==\n", cfg$B_sample,
              as.numeric(Sys.time() - t0, units = "mins")))
  print(outA[, .(quantity, estimate = round(estimate, 3), model_CI = sprintf("[%.2f, %.2f]", model_lo, model_hi),
                 boot_CI = sprintf("[%.2f, %.2f]", boot_lo, boot_hi), boot_sd = round(boot_sd, 3),
                 width_ratio = round(width_ratio, 2), n_draws)])
}

# ---- (B) log-scale reference ----------------------------------------------------
if ("B" %in% cfg$parts && !cfg$fig_only) {
  form_log <- update(form_primary, log(mean_chain) ~ .)
  ml <- bam(form_log, data = ana, method = "fREML")
  obs <- conv(ml, ana, with_se = TRUE)
  t0 <- Sys.time()
  dr <- day_boot(ana, function(xb) conv(bam(form_log, data = xb, method = "fREML"), xb), cfg$B_sample, seed = 2)
  outB <- summarise_boot(obs$est, obs$se, dr, cfg$B_sample, "primary GAMM on log chain length, 3 m - 16 m")
  outB[, `:=`(n = nrow(ana), n_days = uniqueN(ana$day))]
  # same columns as 25_log_contrast.csv (set, contrast, lo, hi); 09b and 10 read the interval
  outB[, set := fifelse(quantity == "season_averaged", "season-averaged", "median-day")]
  outB[, `:=`(contrast = estimate, lo = boot_lo, hi = boot_hi)]
  fwrite(outB, file.path(cfg$out_dir, "27c_log_contrast_boot.csv"))
  cat(sprintf("\n== (B) log-scale reference, day-cluster bootstrap (B = %d, %.1f min) ==\n", cfg$B_sample,
              as.numeric(Sys.time() - t0, units = "mins")))
  print(outB[, .(quantity, estimate = round(estimate, 3), model_CI = sprintf("[%.3f, %.3f]", model_lo, model_hi),
                 boot_CI = sprintf("[%.3f, %.3f]", boot_lo, boot_hi), boot_sd = round(boot_sd, 4),
                 width_ratio = round(width_ratio, 2), n_draws)])
  # Svea marine-span contrasts as shares of the reference (scaled 8 -> 7 salinity units)
  sv_csv <- file.path(cfg$out_dir, "09b_svea_contrasts.csv")
  if (file.exists(sv_csv)) {
    sv <- fread(sv_csv)[grepl("^marine span", set)]
    ref <- outB[set == "season-averaged"]
    span <- 7 / 8
    sh <- sv[, .(set, contrast_scaled = -contrast * span,
                 share_point = 100 * (-contrast * span) / ref$estimate,
                 share_at_ref_hi = 100 * (-contrast * span) / ref$boot_hi,    # reference upper bound -> smallest share
                 share_at_ref_lo = 100 * (-contrast * span) / ref$boot_lo,    # reference lower bound -> largest share
                 share_svea_lo_at_ref_point = 100 * (-ci_lo * span) / ref$estimate)]
    fwrite(sh, file.path(cfg$out_dir, "27c_svea_share_of_reference.csv"))
    cat("\nSvea marine-span contrasts as a share of the Tangesund reference (27c_svea_share_of_reference.csv):\n")
    print(sh[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)])
  }
}

# ---- (C) per-period contrasts ---------------------------------------------------
form_period <- mean_chain ~ period + s(depth, by = period, k = 5) + s(yday, k = 8) + s(day, bs = "re")
med_by_period <- ana[, .(yday = median(yday)), by = period]      # fixed evaluation points, as in 13
period_contrasts <- function(m, dat, with_se = FALSE) {
  r <- lapply(levels(dat$period), function(pp) {
    g <- data.frame(depth = c(cfg$surf, cfg$deep), period = factor(pp, levels = levels(dat$period)),
                    yday = med_by_period[period == pp, yday], day = dat$day[1])
    Xp <- predict(m, g, type = "lpmatrix"); cv <- Xp[1, ] - Xp[2, ]
    c(est = sum(cv * coef(m)), se = sqrt(drop(cv %*% vcov(m) %*% cv)))
  })
  est <- setNames(vapply(r, `[[`, numeric(1), "est"), levels(dat$period))
  if (!with_se) return(est)
  list(est = est, se = setNames(vapply(r, `[[`, numeric(1), "se"), levels(dat$period)))
}
if ("C" %in% cfg$parts && !cfg$fig_only) {
  mp <- bam(form_period, data = ana, method = "fREML")
  obs <- period_contrasts(mp, ana, with_se = TRUE)
  t0 <- Sys.time()
  dr <- day_boot(ana, function(xb) period_contrasts(bam(form_period, data = xb, method = "fREML"), xb),
                 cfg$B_sample, seed = 3, strata = "period")
  outC <- summarise_boot(obs$est, obs$se, dr, cfg$B_sample, "chain length (cells), 3 m - 16 m, period model (13)")
  outC[, response := "Chain length (cells)"]
  # chain-diatom share (Supplementary Table 3), on every QC-pass sample
  cmp <- fread(cfg$comp); setnames(cmp, "f_chain diatom", "f_chain_diatom")
  dc <- merge(D[, .(sample, depth, yday, day, period)], cmp[, .(sample, f_chain_diatom)], by = "sample")[is.finite(f_chain_diatom)]
  form_share <- f_chain_diatom ~ period + s(depth, by = period, k = 5) + s(yday, k = 8) + s(day, bs = "re")
  ms <- bam(form_share, data = dc, method = "fREML")
  obs_s <- period_contrasts(ms, dc, with_se = TRUE)
  dr_s <- day_boot(dc, function(xb) period_contrasts(bam(form_share, data = xb, method = "fREML"), xb),
                   cfg$B_sample, seed = 4, strata = "period")
  outS <- summarise_boot(obs_s$est, obs_s$se, dr_s, cfg$B_sample, "chain-diatom cell share, 3 m - 16 m, period model (13)")
  outS[, response := "Chain-diatom share"]
  outC <- rbind(outC, outS)
  outC[, `:=`(period = quantity, n = fifelse(response == "Chain length (cells)", nrow(ana), nrow(dc)),
              n_days = fifelse(response == "Chain length (cells)", uniqueN(ana$day), uniqueN(dc$day)))]
  cnt <- rbind(ana[, .(response = "Chain length (cells)", n_period = .N, days_period = uniqueN(day)), by = period],
               dc[, .(response = "Chain-diatom share", n_period = .N, days_period = uniqueN(day)), by = period])
  cnt[, period := as.character(period)]
  outC <- merge(outC, cnt, by = c("response", "period"), all.x = TRUE, sort = FALSE)
  fwrite(outC, file.path(cfg$out_dir, "27c_period_contrast_boot.csv"))
  cat(sprintf("\n== (C) per-period contrasts, days resampled within period (B = %d, %.1f min) ==\n", cfg$B_sample,
              as.numeric(Sys.time() - t0, units = "mins")))
  print(outC[, .(response, period, estimate = round(estimate, 3), model_CI = sprintf("[%.2f, %.2f]", model_lo, model_hi),
                 boot_CI = sprintf("[%.2f, %.2f]", boot_lo, boot_hi), width_ratio = round(width_ratio, 2), n_draws, days_period)])
}

# ---- Supplementary Fig. 2 from the bootstrap intervals -------------------------
if ("C" %in% cfg$parts || cfg$fig_only) {
  pc <- fread(file.path(cfg$out_dir, "27c_period_contrast_boot.csv"))
  pc[, response := factor(response, levels = c("Chain length (cells)", "Chain-diatom share"))]
  pc[, period := factor(period, levels = c("Oct", "Sep", "Aug"))]
  pc[, sig := boot_lo > 0 | boot_hi < 0]
  pS2 <- ggplot(pc, aes(estimate, period, colour = sig)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbar(aes(xmin = boot_lo, xmax = boot_hi), orientation = "y", width = 0.25) + geom_point(size = 2.4) +
    facet_wrap(~ response, scales = "free_x") +
    scale_colour_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey60"), guide = "none") +
    theme_paper() +
    labs(x = "Surface (3 m) - deep (16 m) contrast (day-cluster bootstrap 95% CI)", y = NULL,
         title = "Surface-deep contrasts by period")
  ggsave(file.path(cfg$out_dir, "27c_contrasts_panel.png"), pS2, width = 6.5, height = 3.4, dpi = FIG_DPI)
  cat(sprintf("Supplementary Fig. 2 -> %s/27c_contrasts_panel.png\n", cfg$out_dir))
}

# ---- shared by (D) and (E): the ROI table of R/20_within_taxon_chain.R ---------
load_roi_model_table <- function() {
  metaq <- fread(cfg$meta)[qc_pass == TRUE, .(sample, depth, yday, date)]
  metaq[, date := as.Date(date)]
  d <- fread(cfg$roi_file, select = c("sample", "class_name", "biovolume_um3", "numBlobs", "cell_count", "cell_count_resolved"),
             showProgress = FALSE)
  d <- d[numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 & !class_name %in% cfg$junk]
  d <- d[!is.na(cell_count) & cell_count != -1]
  d[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1, round(cell_count_resolved), 1L)]
  d <- merge(d, metaq, by = "sample")
  d[, `:=`(day = factor(date), samp = factor(sample), is_weak = grepl(paste(cfg$weak_taxa, collapse = "|"), class_name))]
  keep <- d[, .N, by = class_name][N >= cfg$min_taxon_model, class_name]
  dm <- d[class_name %in% keep][, .(sample, class_name, nchain, depth, yday, day, samp, is_weak)]
  dm[, taxon := factor(class_name)]
  rm(d); invisible(gc())
  dm
}
# 3 m - 16 m contrast of the fixed part of a fitted ROI model; random-effect
# lpmatrix columns are zeroed (same as predict(..., exclude =))
fixed_contrast <- function(m, g, re = c("s(day)", "s(day,depth)", "s(samp)")) {
  Xp <- predict(m, g, type = "lpmatrix"); cv <- Xp[1, ] - Xp[2, ]
  for (sm in m$smooth) if (sm$label %in% re) cv[sm$first.para:sm$last.para] <- 0
  c(est = sum(cv * coef(m)), se = sqrt(drop(cv %*% vcov(m) %*% cv)))
}

# ---- (D) ROI-level within-taxon models (Table 2) ----------------------------------
if ("D" %in% cfg$parts && !cfg$fig_only) {
  dm <- load_roi_model_table()
  cat(sprintf("\n== (D) ROI-level models: %d ROIs, %d classes, %d days ==\n", nrow(dm), uniqueN(dm$taxon), uniqueN(dm$day)))

  roi_fit <- function(dd, taxon_term, nthreads = 1L) {
    dd <- copy(dd)[, taxon := droplevels(factor(taxon))]
    f <- if (taxon_term) log(nchain) ~ taxon + s(depth, k = 5) + s(yday, k = 8) + s(day, bs = "re") + s(samp, bs = "re")
         else            log(nchain) ~         s(depth, k = 5) + s(yday, k = 8) + s(day, bs = "re") + s(samp, bs = "re")
    m <- bam(f, data = dd, method = "fREML", discrete = TRUE, nthreads = nthreads)
    g <- data.frame(depth = c(cfg$surf, cfg$deep), taxon = dd$taxon[1], yday = median(dd$yday), day = dd$day[1], samp = dd$samp[1])
    Xp <- predict(m, g, type = "lpmatrix"); cv <- Xp[1, ] - Xp[2, ]
    c(est = sum(cv * coef(m)), se = sqrt(drop(cv %*% vcov(m) %*% cv)))
  }
  sets <- list(
    list(label = "all detectors",                        sub = quote(TRUE),                          taxon = TRUE),
    list(label = "accurate/moderate detectors only",     sub = quote(is_weak == FALSE),              taxon = TRUE),
    list(label = "excluding Chaetoceros",                sub = quote(!grepl("Chaetoceros", class_name)), taxon = TRUE),
    list(label = "all detectors, no taxon factor (total contrast)",                     sub = quote(TRUE),             taxon = FALSE),
    list(label = "accurate/moderate detectors, no taxon factor (total contrast)",       sub = quote(is_weak == FALSE), taxon = FALSE))
  # time one fit to decide the bootstrap scope: all rows if a fit takes under ~90 s,
  # the two main rows under ~5 min, otherwise a random depth slope per day instead
  t0 <- Sys.time(); obs1 <- roi_fit(dm, TRUE, nthreads = 4L)
  fit_min <- as.numeric(Sys.time() - t0, units = "mins")
  cat(sprintf("one all-detector fit: %.1f min (4 threads); est %.3f se %.4f\n", fit_min, obs1["est"], obs1["se"]))
  boot_all <- fit_min * 4 < 1.5
  boot_headline <- fit_min * 4 < 5
  nthr <- as.integer(Sys.getenv("BOOT_THREADS", "8"))
  draws <- list()
  outD <- rbindlist(lapply(sets, function(s) {
    dd <- dm[eval(s$sub)]
    obs <- if (identical(s$label, "all detectors")) obs1 else roi_fit(dd, s$taxon, nthreads = nthr)
    headline <- s$label %in% c("all detectors", "accurate/moderate detectors only")
    if (boot_all || (boot_headline && headline)) {
      # one seed per ROI subset: within-taxon and total contrasts share resamples (paired ratio)
      seed <- if (grepl("accurate", s$label)) 11L else if (grepl("excluding", s$label)) 12L else 10L
      t1 <- Sys.time()
      dr <- day_boot(dd, function(xb) roi_fit(xb, s$taxon), cfg$B_roi, seed = seed)
      cat(sprintf("  %-62s bootstrap B = %d done in %.1f min (%d draws)\n", s$label, cfg$B_roi,
                  as.numeric(Sys.time() - t1, units = "mins"), nrow(dr)))
      draws[[s$label]] <<- dr[, "est"]
      r <- summarise_boot(c(contrast = obs[["est"]]), c(contrast = obs[["se"]]),
                          matrix(dr[, "est"], ncol = 1, dimnames = list(NULL, "contrast")), cfg$B_roi, s$label)
      r[, interval := "day-cluster bootstrap"]
      r
    } else {
      # random depth slope per day on top of the day and sample intercepts
      dd2 <- copy(dd)[, taxon := droplevels(factor(taxon))]
      f <- if (s$taxon) log(nchain) ~ taxon + s(depth, k = 5) + s(yday, k = 8) + s(day, bs = "re") + s(day, depth, bs = "re") + s(samp, bs = "re")
           else         log(nchain) ~         s(depth, k = 5) + s(yday, k = 8) + s(day, bs = "re") + s(day, depth, bs = "re") + s(samp, bs = "re")
      t1 <- Sys.time()
      m <- bam(f, data = dd2, method = "fREML", discrete = TRUE, nthreads = nthr)
      cat(sprintf("  %-62s random-depth-slope refit in %.1f min\n", s$label, as.numeric(Sys.time() - t1, units = "mins")))
      g <- data.frame(depth = c(cfg$surf, cfg$deep), taxon = dd2$taxon[1], yday = median(dd2$yday), day = dd2$day[1], samp = dd2$samp[1])
      # zero the RE columns: a day's random slope does not cancel between 3 m and 16 m
      rs <- fixed_contrast(m, g); est_rs <- rs[["est"]]; se_rs <- rs[["se"]]
      data.table(quantity = "contrast", label = s$label, estimate = obs[["est"]], model_se = obs[["se"]],
                 model_lo = obs[["est"]] - 2 * obs[["se"]], model_hi = obs[["est"]] + 2 * obs[["se"]],
                 boot_median = est_rs, boot_sd = se_rs, boot_lo = est_rs - 2 * se_rs, boot_hi = est_rs + 2 * se_rs,
                 width_ratio = se_rs / obs[["se"]], B = 0L, n_draws = 0L,
                 interval = "random depth slope per day (model-based), estimate in boot_median")
    }
  }))
  # within-taxon share = within / total on the same ROIs, paired by draw where both were bootstrapped
  share_rows <- rbindlist(lapply(list(c("all detectors", "all detectors, no taxon factor (total contrast)"),
                                      c("accurate/moderate detectors only", "accurate/moderate detectors, no taxon factor (total contrast)")),
    function(pr) {
      w <- outD[label == pr[1]]; tt <- outD[label == pr[2]]
      pt <- 100 * w$estimate / tt$estimate
      if (all(pr %in% names(draws))) {
        n <- min(length(draws[[pr[1]]]), length(draws[[pr[2]]]))
        rat <- 100 * draws[[pr[1]]][seq_len(n)] / draws[[pr[2]]][seq_len(n)]
        data.table(quantity = "within-taxon share (%)", label = pr[1], estimate = pt, model_se = NA_real_, model_lo = NA_real_, model_hi = NA_real_,
                   boot_median = median(rat), boot_sd = sd(rat), boot_lo = unname(quantile(rat, .025)), boot_hi = unname(quantile(rat, .975)),
                   width_ratio = NA_real_, B = cfg$B_roi, n_draws = n, interval = "day-cluster bootstrap, paired draws")
      } else data.table(quantity = "within-taxon share (%)", label = pr[1], estimate = pt, model_se = NA_real_, model_lo = NA_real_, model_hi = NA_real_,
                        boot_median = NA_real_, boot_sd = NA_real_, boot_lo = NA_real_, boot_hi = NA_real_,
                        width_ratio = NA_real_, B = 0L, n_draws = 0L, interval = "total not bootstrapped")
    }))
  outD <- rbind(outD, share_rows)
  outD[quantity == "contrast", `:=`(pct = 100 * (exp(estimate) - 1), pct_lo = 100 * (exp(boot_lo) - 1), pct_hi = 100 * (exp(boot_hi) - 1),
                                    model_pct_lo = 100 * (exp(model_lo) - 1), model_pct_hi = 100 * (exp(model_hi) - 1))]
  nroi <- setNames(vapply(sets, function(s) dm[eval(s$sub), .N], integer(1)), vapply(sets, `[[`, "", "label"))
  outD[, n_roi := nroi[label]]
  outD[, n_days := uniqueN(dm$day)]
  fwrite(outD, file.path(cfg$out_dir, "27c_roi_model_boot.csv"))
  cat("\n-- Table 2 rows with day-cluster intervals (27c_roi_model_boot.csv) --\n")
  print(outD[, .(label, estimate = round(estimate, 3), model_CI = sprintf("[%.3f, %.3f]", model_lo, model_hi),
                 boot_CI = sprintf("[%.3f, %.3f]", boot_lo, boot_hi), pct = round(pct, 1),
                 pct_CI = sprintf("[%+.0f%%, %+.0f%%]", pct_lo, pct_hi), width_ratio = round(width_ratio, 2), interval)])
}

# ---- (E) per-taxon ROI models (Supplementary Table 2) ---------------------------
# 20's per-class model, refitted with a random depth slope per day. The first
# contrast (no slope) reproduces 20_within_taxon_pertaxon_model.csv as a check.
if ("E" %in% cfg$parts && !cfg$fig_only) {
  if (!exists("dm")) dm <- load_roi_model_table()
  nthr <- as.integer(Sys.getenv("BOOT_THREADS", "8"))
  fit_rs_one <- function(dd) {
    dd <- copy(dd)[, `:=`(day = droplevels(day), samp = droplevels(samp))]
    base <- log(nchain) ~ s(depth, k = 5) + s(yday, k = 8) + s(day, bs = "re") + s(samp, bs = "re")
    fit <- function(f) tryCatch(bam(f, data = dd, method = "fREML", discrete = TRUE, nthreads = nthr),
                                error = function(e) { message("  fit failed: ", conditionMessage(e)); NULL })
    m0 <- fit(base); m1 <- fit(update(base, . ~ . + s(day, depth, bs = "re")))
    g <- data.frame(depth = c(cfg$surf, cfg$deep), yday = median(dd$yday), day = dd$day[1], samp = dd$samp[1])
    a <- if (is.null(m0)) c(est = NA_real_, se = NA_real_) else fixed_contrast(m0, g)
    b <- if (is.null(m1)) c(est = NA_real_, se = NA_real_) else fixed_contrast(m1, g)
    data.table(est = a[["est"]], se = a[["se"]], est_rs = b[["est"]], se_rs = b[["se"]])
  }
  t0 <- Sys.time()
  outE <- dm[, c(list(n_roi = .N, n_days = uniqueN(day), n_surf = sum(depth == cfg$surf), n_deep = sum(depth == cfg$deep)),
                 fit_rs_one(.SD)), by = class_name]
  outE[, `:=`(lo = est - 2 * se, hi = est + 2 * se, lo_rs = est_rs - 2 * se_rs, hi_rs = est_rs + 2 * se_rs,
              width_ratio = se_rs / se, pct = 100 * (exp(est) - 1), pct_rs = 100 * (exp(est_rs) - 1),
              pct_lo_rs = 100 * (exp(est_rs - 2 * se_rs) - 1), pct_hi_rs = 100 * (exp(est_rs + 2 * se_rs) - 1))]
  outE[, `:=`(excludes_zero_model = lo > 0 | hi < 0, excludes_zero_rs = lo_rs > 0 | hi_rs < 0)]
  setorder(outE, -n_roi)
  fwrite(outE, file.path(cfg$out_dir, "27c_pertaxon_model_boot.csv"))
  cat(sprintf("\n== (E) per-taxon ROI models with a random depth slope per day (%.1f min; 27c_pertaxon_model_boot.csv) ==\n",
              as.numeric(Sys.time() - t0, units = "mins")))
  print(outE[, .(class_name, n_roi, model = sprintf("%+.0f%% [%+.0f, %+.0f]", pct, 100 * (exp(lo) - 1), 100 * (exp(hi) - 1)),
                 random_slope = sprintf("%+.0f%% [%+.0f, %+.0f]", pct_rs, pct_lo_rs, pct_hi_rs),
                 width_ratio = round(width_ratio, 1), excludes_zero_model, excludes_zero_rs)])
}

cat("\ndone.\n")
