# =============================================================================
# Sensitivity of the primary chain-length GAMMs (03, 13) on the same thresholded
# set: prior weights = n_chain_roi, k.check and a larger basis, a Gamma(log) refit,
# season-averaged vs median-day contrasts, and the log-scale 3 m - 16 m reference.
#
# Reads:  outputs/tangesund/02_sample_meta.csv, 03_chain_metrics.csv.
# Writes: outputs/tangesund/25_precision_by_depth.csv, 25_weighted_profile.csv,
#         25_weighted_contrasts.csv, 25_kcheck.csv, 25_highk_profile.csv,
#         25_gamma_profile.csv, 25_contrast_conventions.csv, 25_log_contrast.csv
#         (read by 09b and 10).
# Deps:   data.table, mgcv.
# Run:    Rscript R/25_weight_sensitivity.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(mgcv) })

cfg <- list(out_dir = "outputs/tangesund", meta = "outputs/tangesund/02_sample_meta.csv",
            depths = c(3, 6, 8, 11, 13, 16), min_chain_roi = 20)

meta <- fread(cfg$meta)[, .(sample, depth, date, yday, qc_pass)]
chn  <- fread(file.path(cfg$out_dir, "03_chain_metrics.csv"))[
  , .(sample, n_chain_roi, mean_chain)]
D <- merge(meta, chn, by = "sample", all.x = TRUE)[qc_pass == TRUE & is.finite(mean_chain)]
D[, `:=`(period = factor(month.abb[as.integer(format(as.Date(date), "%m"))],
                         levels = c("Aug", "Sep", "Oct")),
         day = factor(date))]
ana <- D[n_chain_roi >= cfg$min_chain_roi]                 # chain-trait thresholded set
cat(sprintf("analysis set: %d samples (>= %d chain ROIs)\n", nrow(ana), cfg$min_chain_roi))

# ---- (A) chain-ROI count by depth --------------------------------------------
spread <- ana[, .(n = .N, med_roi = as.numeric(median(n_chain_roi)),
                  q25 = as.numeric(quantile(n_chain_roi, .25)),
                  q75 = as.numeric(quantile(n_chain_roi, .75))), by = depth][order(depth)]
cat("\n== (A) chain-ROI count per sample, by depth (the precision the weights encode) ==\n")
print(spread)
cat(sprintf("surface:deep median precision ratio = %.2f\n",
            spread[depth == 3, med_roi] / spread[depth == 16, med_roi]))
fwrite(spread, file.path(cfg$out_dir, "25_precision_by_depth.csv"))

# ---- (B) depth profile, unweighted vs weighted ------------------------------
# 03 specification; day RE excluded from predictions
form <- mean_chain ~ s(depth, k = 5) + s(yday, k = 8) + ti(depth, yday, k = c(5, 6)) +
  s(day, bs = "re")
m_un <- bam(form, data = ana, method = "fREML")
m_wt <- bam(form, data = ana, method = "fREML", weights = ana$n_chain_roi)

# averaged over the observed days of year, not predicted at the median yday
fitted_profile <- function(m, lab) {
  rbindlist(lapply(cfg$depths, function(dd) {
    g <- data.frame(depth = dd, yday = ana$yday, day = ana$day[1])
    p <- predict(m, g, type = "response", exclude = "s(day)")
    data.table(set = lab, depth = dd, fit = mean(as.numeric(p)))
  }))
}
prof <- rbind(fitted_profile(m_un, "unweighted"), fitted_profile(m_wt, "weighted"))
raw <- ana[, .(raw_mean = mean(mean_chain),
               raw_wtd  = sum(n_chain_roi * mean_chain) / sum(n_chain_roi)), by = depth][order(depth)]
prof <- merge(prof, raw, by = "depth"); setorder(prof, set, depth)
cat("\n== (B) season-marginal fitted depth profile, unweighted vs precision-weighted ==\n")
cat("   raw_mean = unweighted mean of per-sample means (the reported Fig. 2a profile)\n")
cat("   raw_wtd  = the same means weighted by chain-ROI count, which is arithmetically\n")
cat("              the abundance-weighted (pooled) profile the manuscript already reports\n")
print(prof[, .(set, depth, fit = round(fit, 3), raw_mean = round(raw_mean, 3),
               raw_wtd = round(raw_wtd, 3))])
wide <- dcast(prof, depth ~ set, value.var = "fit")
cat("\nfitted surface-minus-deep, unweighted: ", round(wide[depth == 3, unweighted] - wide[depth == 16, unweighted], 3),
    " | weighted: ", round(wide[depth == 3, weighted] - wide[depth == 16, weighted], 3), "\n", sep = "")
cat("monotonic, unweighted fit: ", all(diff(wide$unweighted) < 0),
    " | weighted fit: ", all(diff(wide$weighted) < 0),
    " | raw: ", all(diff(raw$raw_mean) < 0),
    " | raw weighted: ", all(diff(raw$raw_wtd) < 0), "\n", sep = "")
fwrite(prof, file.path(cfg$out_dir, "25_weighted_profile.csv"))

# ---- (C) monthly contrasts, unweighted vs weighted --------------------------
# same specification and lpmatrix contrast as R/13_contrasts.R
contrasts_for <- function(use_weights, lab) {
  m <- if (use_weights)
    bam(mean_chain ~ period + s(depth, by = period, k = 5) + s(yday, k = 8) + s(day, bs = "re"),
        data = ana, method = "fREML", weights = ana$n_chain_roi)
  else
    bam(mean_chain ~ period + s(depth, by = period, k = 5) + s(yday, k = 8) + s(day, bs = "re"),
        data = ana, method = "fREML")
  rbindlist(lapply(levels(ana$period), function(pp) {
    sub <- ana[period == pp]; if (nrow(sub) < 40) return(NULL)
    g <- data.frame(depth = c(3, 16), period = factor(pp, levels = levels(ana$period)),
                    yday = median(sub$yday), day = sub$day[1])
    Xp <- predict(m, g, type = "lpmatrix")
    cvec <- Xp[1, ] - Xp[2, ]
    est <- sum(cvec * coef(m)); se <- sqrt(drop(cvec %*% vcov(m) %*% cvec))
    data.table(set = lab, period = pp, est = est, lo = est - 2 * se, hi = est + 2 * se,
               se = se, n = nrow(sub))
  }))
}
ctr <- rbind(contrasts_for(FALSE, "unweighted"), contrasts_for(TRUE, "weighted"))
cat("\n== (C) surface(3 m) - deep(16 m) contrast per period, same thresholded set ==\n")
print(ctr[, .(set, period, est = round(est, 3), lo = round(lo, 3), hi = round(hi, 3),
              sig = ifelse(lo > 0 | hi < 0, "*", ""), n)])
fwrite(ctr, file.path(cfg$out_dir, "25_weighted_contrasts.csv"))

# ---- (D) basis adequacy on the primary model --------------------------------
cat("\n== (D) k.check on the primary (unweighted) chain-length model ==\n")
set.seed(1)  # k.check simulates its p-value; seed so the CSV reproduces
kc <- k.check(m_un); print(round(kc, 4))
set.seed(1)
kcw <- k.check(m_wt)
fwrite(data.table(set = "unweighted", term = rownames(kc), as.data.table(kc)),
       file.path(cfg$out_dir, "25_kcheck.csv"))
fwrite(data.table(set = "weighted", term = rownames(kcw), as.data.table(kcw)),
       file.path(cfg$out_dir, "25_kcheck.csv"), append = TRUE)

# larger basis: if the fit does not move, a low k-index is residual structure, not k
form_k <- mean_chain ~ s(depth, k = 6) + s(yday, k = 20) + ti(depth, yday, k = c(6, 10)) +
  s(day, bs = "re")
m_k <- bam(form_k, data = ana, method = "fREML")
pk <- fitted_profile(m_k, "high_k")
cat("\n== (D2) refit with a larger basis (yday k = 20, depth k = 6) ==\n")
print(pk[, .(depth, fit = round(fit, 3))])
cat("high-k surface-minus-deep: ", round(pk[depth == 3, fit] - pk[depth == 16, fit], 3), "\n", sep = "")
fwrite(pk, file.path(cfg$out_dir, "25_highk_profile.csv"))     # 3 columns, unlike 25_weighted_profile.csv

# ---- (E) family check: Gamma(log) -------------------------------------------
m_gam <- tryCatch(bam(form, data = ana, family = Gamma(link = "log"), method = "fREML"),
                  error = function(e) { message("Gamma refit failed: ", conditionMessage(e)); NULL })
if (!is.null(m_gam)) {
  pg <- fitted_profile(m_gam, "gamma_log")
  cat("\n== (E) Gamma(log) season-marginal profile ==\n")
  print(pg[, .(depth, fit = round(fit, 3))])
  cat("Gamma surface-minus-deep: ", round(pg[depth == 3, fit] - pg[depth == 16, fit], 3),
      " | monotonic: ", all(diff(pg$fit) < 0), "\n", sep = "")
  fwrite(pg, file.path(cfg$out_dir, "25_gamma_profile.csv"))
}
# ---- (F) contrast conventions and the log-scale reference -------------------
# season-averaged vs median-day contrast for each model above
med_yday <- median(ana$yday)
conv_one <- function(m, lab, link_identity = TRUE) {
  Xa <- function(dp) colMeans(predict(m, data.frame(depth = dp, yday = ana$yday, day = ana$day[1]), type = "lpmatrix"))
  Xm <- function(dp) predict(m, data.frame(depth = dp, yday = med_yday, day = ana$day[1]), type = "lpmatrix")[1, ]
  fa <- function(dp) mean(as.numeric(predict(m, data.frame(depth = dp, yday = ana$yday, day = ana$day[1]),
                                             type = "response", exclude = "s(day)")))
  fm <- function(dp) as.numeric(predict(m, data.frame(depth = dp, yday = med_yday, day = ana$day[1]),
                                        type = "response", exclude = "s(day)"))
  se_of <- function(cv) if (link_identity) sqrt(drop(cv %*% vcov(m) %*% cv)) else NA_real_
  rbind(data.table(model = lab, convention = "averaged over observed days of year",
                   surface = fa(3), deep = fa(16), contrast = fa(3) - fa(16), se = se_of(Xa(3) - Xa(16))),
        data.table(model = lab, convention = "median day of year",
                   surface = fm(3), deep = fm(16), contrast = fm(3) - fm(16), se = se_of(Xm(3) - Xm(16))))
}
conv <- rbind(conv_one(m_un, "primary, unweighted"), conv_one(m_wt, "precision-weighted"),
              conv_one(m_k, "high k (yday k = 20, depth k = 6)"),
              if (!is.null(m_gam)) conv_one(m_gam, "Gamma(log)", link_identity = FALSE))
conv[, `:=`(lo = contrast - 2 * se, hi = contrast + 2 * se, median_yday = med_yday, n = nrow(ana))]
fwrite(conv, file.path(cfg$out_dir, "25_contrast_conventions.csv"))
cat("\n== (F) surface-minus-deep contrast by convention (25_contrast_conventions.csv) ==\n")
print(conv[, .(model, convention, surface = round(surface, 2), deep = round(deep, 2), contrast = round(contrast, 3),
               CI = ifelse(is.na(se), "", sprintf("[%.2f, %.2f]", lo, hi)))])

# log-scale refit: difference of mean logs, the scale 09b models
m_log <- bam(update(form, log(mean_chain) ~ .), data = ana, method = "fREML")
lg <- conv_one(m_log, "log(mean_chain), unweighted")
lg[, `:=`(lo = contrast - 2 * se, hi = contrast + 2 * se, n = nrow(ana), n_days = uniqueN(ana$day),
          raw_mean_log_surface = ana[depth == 3, mean(log(mean_chain))],
          raw_mean_log_deep = ana[depth == 16, mean(log(mean_chain))])]
lg[, set := fifelse(grepl("averaged", convention), "season-averaged", "median-day")]
setcolorder(lg, c("set", "model", "convention"))
fwrite(lg, file.path(cfg$out_dir, "25_log_contrast.csv"))
cat("\n== (F2) 3 m - 16 m contrast of log per-sample chain length (25_log_contrast.csv; the Fig. 6d reference) ==\n")
print(lg[, .(set, surface = round(surface, 3), deep = round(deep, 3), contrast = round(contrast, 3),
             CI = sprintf("[%.3f, %.3f]", lo, hi), raw_diff = round(raw_mean_log_surface - raw_mean_log_deep, 3))])
cat("\ndone.\n")
