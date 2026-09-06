# =============================================================================
# Within-taxon chain length by depth: per-taxon partition terms (A), ROI-level
# models with and without the taxon factor (B; Table 2), per-taxon models (C),
# partition variants (D), per-taxon contrasts on 09b's scale (E), and
# negative-binomial refits (B2).
#
# Reads:  data/roi_level/roi_2016.txt.gz, outputs/tangesund/02_sample_meta.csv,
#         data/functional_group_map.csv, outputs/tangesund/19_counter_validation.csv;
#         outputs/tangesund/chaetoceros_validation/23e_class_purity.csv (optional; run 23e first,
#         or the purified-Chaetoceros variant is dropped).
# Writes: outputs/tangesund/20_within_taxon_pertaxon.csv, 20_within_taxon_chain.csv,
#         20_class_counts.csv, 20_within_taxon_pertaxon_model.csv,
#         20_sample_log_contrasts.csv, 20_partition_variants.csv, 20_within_taxon_nb.csv.
# Deps:   data.table, mgcv.
# Run:    Rscript R/20_within_taxon_chain.R   (the negative-binomial refits take about 2 h)
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(mgcv) })

cfg <- list(
  roi_file  = "data/roi_level/roi_2016.txt.gz",
  meta_csv  = "outputs/tangesund/02_sample_meta.csv",
  fg_map    = "data/functional_group_map.csv",
  val_csv   = "outputs/tangesund/19_counter_validation.csv",
  out_dir   = "outputs/tangesund",
  junk      = c("Debris", "Beads", "Air_bubbles", "Zooplankton"),
  surf = 3, deep = 16,
  min_taxon_roi   = 5,     # ROIs per taxon and depth in the partition, as 17
  min_taxon_model = 200,   # ROIs per taxon in the ROI-level models
  min_taxon_show  = 100,   # ROIs per taxon at each depth in the printed table
  # detectors flagged weak in Table S6 (same set as 19_chain_counter_validation.R)
  weak_taxa = c("Chaetoceros", "Asterionellopsis", "Pauliella")
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

meta <- fread(cfg$meta_csv)[qc_pass == TRUE, .(sample, depth, yday, date)]
meta[, date := as.Date(date)]

d <- fread(cfg$roi_file, select = c("sample","class_name","biovolume_um3","numBlobs",
                                    "cell_count","cell_count_resolved"), showProgress = FALSE)
d <- d[numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 & !class_name %in% cfg$junk]
d <- d[!is.na(cell_count) & cell_count != -1]                          # chain-eligible only
d[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1,
                      round(cell_count_resolved), 1L)]
d <- merge(d, meta, by = "sample")
d[, `:=`(day = factor(date), samp = factor(sample),
         is_weak = grepl(paste(cfg$weak_taxa, collapse = "|"), class_name))]
cat(sprintf("chain-eligible ROIs: %d over %d samples, %d days | weak-detector share %.1f%%\n",
            nrow(d), uniqueN(d$sample), uniqueN(d$day), 100 * mean(d$is_weak)))
n_class_any <- uniqueN(d$class_name)        # classes with at least one chain-counted ROI

# =====================  (A) which taxa contribute the within-taxon term  ======
# 17's partition on the 3 m and 16 m ROIs, per taxon: within_i = fbar_i * (L_i,surf - L_i,deep)
sd_pool <- d[depth %in% c(cfg$surf, cfg$deep)]
agg <- sd_pool[, .(L = mean(nchain), nn = .N), by = .(depth, class_name)][nn >= cfg$min_taxon_roi]
agg[, f := nn / sum(nn), by = depth]
w <- merge(agg[depth == cfg$surf, .(class_name, L_s = L, f_s = f, n_s = nn)],
           agg[depth == cfg$deep, .(class_name, L_d = L, f_d = f, n_d = nn)],
           by = "class_name", all = TRUE)
for (col in c("L_s","L_d","f_s","f_d","n_s","n_d")) w[is.na(get(col)), (col) := 0]
w[, `:=`(fbar = (f_s + f_d)/2, Lbar = (L_s + L_d)/2)]
w[, `:=`(within_i = fbar * (L_s - L_d), between_i = Lbar * (f_s - f_d))]
within_tot <- sum(w$within_i); between_tot <- sum(w$between_i)

disp <- fread(cfg$fg_map)[, setNames(display_name, class)]
val  <- fread(cfg$val_csv)[, .(detector, exact_pct, mae, flag)]
# class -> genus detector (as in 19_chain_counter_validation.R); join validation on the detector
map_det <- function(cn) fcase(
  grepl("^Chaetoceros", cn), "chaetoceros",
  grepl("^Skeletonema", cn), "skeletonema",
  grepl("^Pseudo-nitzschia", cn), "pseudo-nitzschia",
  grepl("^Leptocylindrus", cn), "leptocylindrus",
  grepl("^Asterionellopsis", cn), "asterionellopsis",
  grepl("^Cerataulina", cn), "cerataulina_pelagica",
  grepl("^Thalassiosira|^Cyclotella", cn), "thalassiosira",
  grepl("^Guinardia", cn), "guinardia_delicatula",
  grepl("^Thalassionema|^Diatoma", cn), "thalassionema_nitzschioides",
  grepl("^Melosira", cn), "melosira",
  grepl("^Dactyliosolen", cn), "dactyliosolen_fragilissimus",
  grepl("^Peridiniella", cn), "peridiniella_catenata",
  grepl("^Pauliella", cn), "pauliella",
  grepl("^Dinobryon", cn), "dinobryon",
  grepl("^Uroglenopsis", cn), "uroglenopsis",
  grepl("^Eucampia", cn), "eucampia",
  grepl("^Lauderia|^Detonula", cn), "lauderia_detonula",
  default = NA_character_)
w[, taxon := fifelse(class_name %in% names(disp), disp[class_name], class_name)]
w[, detector := map_det(class_name)]
if (anyNA(w$detector))
  stop("chain-eligible classes with no detector: ", paste(unique(w[is.na(detector), class_name]), collapse = ", "))
w <- merge(w, val, by = "detector", all.x = TRUE)[, detector := NULL]
w[, pct_of_within := 100 * within_i / within_tot]
setorder(w, -within_i)

out_a <- w[n_s >= cfg$min_taxon_show & n_d >= cfg$min_taxon_show,
           .(taxon, n_surf = n_s, n_deep = n_d,
             L_surf = round(L_s, 2), L_deep = round(L_d, 2), delta = round(L_s - L_d, 2),
             f_surf = round(f_s, 3), f_deep = round(f_d, 3),
             within_i = round(within_i, 3), pct_of_within = round(pct_of_within, 1),
             exact_pct, mae, flag)]
fwrite(out_a, file.path(cfg$out_dir, "20_within_taxon_pertaxon.csv"))
cat(sprintf("\n== (A) partition: dCWM %.2f = within %.2f (%.0f%%) + between %.2f ==\n",
            within_tot + between_tot, within_tot,
            100 * within_tot / (within_tot + between_tot), between_tot))
print(out_a)
weak_share <- 100 * sum(w[is.na(flag) | flag == "weak", within_i]) / within_tot
cat(sprintf("\nweak-detector taxa contribute %.0f%% of the within-taxon term (Chaetoceros %.0f%%, Asterionellopsis %.0f%%)\n",
            weak_share,
            100 * sum(w[grepl("^Chaetoceros", taxon), within_i]) / within_tot,
            100 * sum(w[grepl("^Asterionellopsis", taxon), within_i]) / within_tot))

# =====================  (B) ROI-level within-taxon depth response  ============
keep <- d[, .N, by = class_name][N >= cfg$min_taxon_model, class_name]
dm <- d[class_name %in% keep][, taxon := factor(class_name)]
fwrite(data.table(classes_with_chain_counted_roi = n_class_any, classes_in_roi_models = length(keep),
                  min_roi_per_class = cfg$min_taxon_model, n_roi_in_models = nrow(dm)),
       file.path(cfg$out_dir, "20_class_counts.csv"))
cat(sprintf("classes with >= 1 chain-counted ROI: %d; with >= %d ROIs (ROI models): %d\n",
            n_class_any, cfg$min_taxon_model, length(keep)))

fit_within <- function(dd, label, taxon_term = TRUE, count_model = FALSE) {
  dd <- copy(dd)[, taxon := droplevels(factor(taxon))]
  if (nlevels(dd$taxon) < 2) return(NULL)
  # taxon_term = FALSE: total contrast (within plus between taxa) on the same scale
  f <- if (taxon_term) log(nchain) ~ taxon + s(depth, k = 5) + s(yday, k = 8) + s(day, bs = "re") + s(samp, bs = "re")
       else             log(nchain) ~         s(depth, k = 5) + s(yday, k = 8) + s(day, bs = "re") + s(samp, bs = "re")
  # count_model = TRUE: negative binomial on nchain - 1 (contrast on its log-mean scale)
  if (count_model) f <- update(f, I(nchain - 1) ~ .)
  m <- if (count_model) tryCatch(bam(f, data = dd, family = nb(), method = "fREML", discrete = TRUE, nthreads = 4),
                                 error = function(e) { message("negative-binomial refit failed: ", conditionMessage(e)); NULL })
       else bam(f, data = dd, method = "fREML", discrete = TRUE, nthreads = 4)
  if (is.null(m)) return(NULL)
  st <- summary(m)$s.table
  # 3 m - 16 m on the depth smooth; the other terms cancel between the two rows
  g <- data.frame(depth = c(cfg$surf, cfg$deep), taxon = dd$taxon[1],
                  yday = median(dd$yday), day = dd$day[1], samp = dd$samp[1])
  Xp <- predict(m, g, type = "lpmatrix"); cv <- Xp[1, ] - Xp[2, ]
  est <- sum(cv * coef(m)); se <- sqrt(drop(cv %*% vcov(m) %*% cv))
  data.table(set = label, family = if (count_model) "negative binomial, nchain - 1" else "Gaussian, log nchain",
             n_roi = nrow(dd), n_taxa = nlevels(dd$taxon), n_days = uniqueN(dd$day),
             dev_expl = summary(m)$dev.expl,
             edf = st["s(depth)","edf"], ref_df = st["s(depth)","Ref.df"],
             F = st["s(depth)", intersect(c("F", "Chi.sq"), colnames(st))[1]], p = st["s(depth)","p-value"],
             est = est, se = se, lo = est - 2*se, hi = est + 2*se,
             ratio = exp(est), pct = 100 * (exp(est) - 1),
             pct_lo = 100 * (exp(est - 2*se) - 1), pct_hi = 100 * (exp(est + 2*se) - 1))
}

res <- rbindlist(list(
  fit_within(dm,                    "all detectors"),
  fit_within(dm[is_weak == FALSE],  "accurate/moderate detectors only"),
  fit_within(dm[!grepl("Chaetoceros", class_name)], "excluding Chaetoceros"),
  fit_within(dm, "all detectors, no taxon factor (total contrast)", taxon_term = FALSE),
  fit_within(dm[is_weak == FALSE], "accurate/moderate detectors, no taxon factor (total contrast)", taxon_term = FALSE)))
fwrite(res, file.path(cfg$out_dir, "20_within_taxon_chain.csv"))
cat("\n== (B) ROI-level within-taxon depth response, log chain length (day + sample RE) ==\n")
print(res[, .(set, n_roi, n_taxa, edf = round(edf,2), F = round(F,1),
              p = signif(p,3), surf_minus_deep = round(est,3),
              CI = sprintf("[%.3f, %.3f]", lo, hi),
              pct = sprintf("%+.0f%% [%+.0f, %+.0f]", pct, pct_lo, pct_hi))])

# =====================  (C) per-taxon ROI-level contrasts  ===================
# the ROI model fitted to one class at a time
fit_one <- function(dd) {
  dd <- copy(dd)[, `:=`(day = droplevels(day), samp = droplevels(samp))]
  m <- tryCatch(bam(log(nchain) ~ s(depth, k = 5) + s(yday, k = 8) + s(day, bs = "re") + s(samp, bs = "re"),
                    data = dd, method = "fREML", discrete = TRUE, nthreads = 4),
                error = function(e) NULL)
  if (is.null(m)) return(data.table(est = NA_real_, se = NA_real_))
  g <- data.frame(depth = c(cfg$surf, cfg$deep), yday = median(dd$yday), day = dd$day[1], samp = dd$samp[1])
  Xp <- predict(m, g, type = "lpmatrix"); cv <- Xp[1, ] - Xp[2, ]
  data.table(est = sum(cv * coef(m)), se = sqrt(drop(cv %*% vcov(m) %*% cv)))
}
pt <- dm[, c(list(n_roi = .N, n_days = uniqueN(day), n_surf = sum(depth == cfg$surf), n_deep = sum(depth == cfg$deep)),
             fit_one(.SD)), by = class_name]
pt[, `:=`(lo = est - 2*se, hi = est + 2*se, pct = 100 * (exp(est) - 1),
          pct_lo = 100 * (exp(est - 2*se) - 1), pct_hi = 100 * (exp(est + 2*se) - 1))]
pt[, taxon := fifelse(class_name %in% names(disp), disp[class_name], class_name)]
pt[, detector := map_det(class_name)]
pt <- merge(pt, val, by = "detector", all.x = TRUE)[, detector := NULL]
setorder(pt, -n_roi)
fwrite(pt, file.path(cfg$out_dir, "20_within_taxon_pertaxon_model.csv"))
cat("\n== (C) per-taxon ROI-level surface-minus-deep contrasts (log chain length) ==\n")
print(pt[, .(taxon, n_roi, flag, est = round(est, 3), CI = sprintf("[%.3f, %.3f]", lo, hi),
             pct = sprintf("%+.0f%%", pct))])

# =====================  (E) per-taxon contrasts on the regional scale  =======
# log per-sample mean chain length for the six most abundant chain-formers (>= 5 ROIs
# per sample x taxon, as 09b), primary GAMM; contrast averaged over observed days as in 25
top6 <- pt[order(-n_roi)][1:6, class_name]
sa <- d[class_name %in% top6, .(lmean = log(mean(nchain)), nn = .N, yday = yday[1], day = day[1], depth = depth[1]),
        by = .(class_name, sample)][nn >= cfg$min_taxon_roi]
fit_sample <- function(dd) {
  dd <- copy(dd)[, day := droplevels(day)]
  m <- bam(lmean ~ s(depth, k = 5) + s(yday, k = 8) + ti(depth, yday, k = c(5, 6)) + s(day, bs = "re"),
           data = dd, method = "fREML")
  Xa <- function(dp) colMeans(predict(m, data.frame(depth = dp, yday = dd$yday, day = dd$day[1]), type = "lpmatrix"))
  Xm <- function(dp) predict(m, data.frame(depth = dp, yday = median(dd$yday), day = dd$day[1]), type = "lpmatrix")[1, ]
  cv <- Xa(cfg$surf) - Xa(cfg$deep); cm <- Xm(cfg$surf) - Xm(cfg$deep)   # the day RE cancels
  data.table(n_samples = nrow(dd), n_days = uniqueN(dd$day), n_surf = sum(dd$depth == cfg$surf), n_deep = sum(dd$depth == cfg$deep),
             est = sum(cv * coef(m)), se = sqrt(drop(cv %*% vcov(m) %*% cv)),
             est_median_day = sum(cm * coef(m)), se_median_day = sqrt(drop(cm %*% vcov(m) %*% cm)),
             raw_mean_log_surf = mean(dd$lmean[dd$depth == cfg$surf]), raw_mean_log_deep = mean(dd$lmean[dd$depth == cfg$deep]))
}
se_tab <- sa[, fit_sample(.SD), by = class_name]
se_tab[, `:=`(lo = est - 2*se, hi = est + 2*se)]
se_tab[, taxon := fifelse(class_name %in% names(disp), disp[class_name], class_name)]
setcolorder(se_tab, c("class_name", "taxon"))
fwrite(se_tab, file.path(cfg$out_dir, "20_sample_log_contrasts.csv"))
cat("\n== (E) per-taxon 3 m - 16 m contrasts of log per-sample chain length, day-of-year averaged (20_sample_log_contrasts.csv) ==\n")
print(se_tab[, .(taxon, n_samples, n_days, est = round(est, 3), CI = sprintf("[%.3f, %.3f]", lo, hi),
                 median_day = round(est_median_day, 3))])

# =====================  (D) partition variants  ==============================
# (i) only taxa with >= 5 ROIs at both depths; (ii) Chaetoceros chain class corrected
# for the non-Chaetoceros frames in 23e_class_purity.csv
partition <- function(ag) {
  ag <- copy(ag)[, f := nn / sum(nn), by = depth]
  ww <- merge(ag[depth == cfg$surf, .(class_name, L_s = L, f_s = f)],
              ag[depth == cfg$deep, .(class_name, L_d = L, f_d = f)], by = "class_name", all = TRUE)
  for (col in c("L_s","L_d","f_s","f_d")) ww[is.na(get(col)), (col) := 0]
  wi <- sum((ww$f_s + ww$f_d)/2 * (ww$L_s - ww$L_d)); be <- sum((ww$L_s + ww$L_d)/2 * (ww$f_s - ww$f_d))
  data.table(dCWM = wi + be, within = wi, between = be, frac_within = 100 * wi / (wi + be))
}
agg0 <- sd_pool[, .(L = mean(nchain), nn = .N), by = .(depth, class_name)][nn >= cfg$min_taxon_roi]
both <- agg0[, .N, by = class_name][N == 2, class_name]
pur_file <- file.path(cfg$out_dir, "chaetoceros_validation", "23e_class_purity.csv")
variants <- list(`primary (taxa zero-filled at a depth with < 5 ROIs)` = partition(agg0),
                 `taxa with >= 5 ROIs at both depths only` = partition(agg0[class_name %in% both]))
if (file.exists(pur_file)) {
  pur <- fread(pur_file)
  aggp <- copy(agg0)[, nn := as.numeric(nn)]                       # fractional counts after purification
  for (dp in c(cfg$surf, cfg$deep)) {
    q <- pur[depth == dp, not_chaeto_pct] / 100; mnc <- pur[depth == dp, mean_det_not_chaeto]
    aggp[depth == dp & class_name == "Chaetoceros_spp_chain", `:=`(L = (L - q * mnc) / (1 - q), nn = nn * (1 - q))]
  }
  variants[["Chaetoceros chain class purified (23e misclassification removed)"]] <- partition(aggp)
}
pv <- rbindlist(variants, idcol = "variant")
fwrite(pv, file.path(cfg$out_dir, "20_partition_variants.csv"))
cat("\n== (D) partition variants (20_partition_variants.csv) ==\n"); print(pv)

# =====================  (B2) negative-binomial refits  =======================
res_nb <- rbindlist(list(
  fit_within(dm,                   "all detectors",                   count_model = TRUE),
  fit_within(dm[is_weak == FALSE], "accurate/moderate detectors only", count_model = TRUE)))
if (nrow(res_nb)) {
  res_nb <- merge(res_nb, res[, .(set, est_gauss = est, se_gauss = se, halfwidth_gauss_pct = (pct_hi - pct_lo) / 2)], by = "set")
  res_nb[, halfwidth_nb_pct := (pct_hi - pct_lo) / 2]
  res_nb[, halfwidth_ratio_nb_over_gauss := halfwidth_nb_pct / halfwidth_gauss_pct]
  fwrite(res_nb, file.path(cfg$out_dir, "20_within_taxon_nb.csv"))
  cat("\n== (B2) negative-binomial refits on nchain - 1 (20_within_taxon_nb.csv) ==\n")
  print(res_nb[, .(set, est = round(est, 3), CI = sprintf("[%.3f, %.3f]", lo, hi), pct = sprintf("%+.0f%% [%+.0f, %+.0f]", pct, pct_lo, pct_hi),
                   gauss_pct_halfwidth = round(halfwidth_gauss_pct, 2), nb_pct_halfwidth = round(halfwidth_nb_pct, 2))])
}

cat("\n===================== WITHIN-TAXON CHAIN SUMMARY =====================\n")
cat(sprintf("(A) weak-detector share of the within-taxon term: %.0f%%\n", weak_share))
cat(sprintf("(B) within-taxon contrast: all detectors %+.0f%%; accurate/moderate %+.0f%% [%+.0f, %+.0f]\n",
            res[set == "all detectors", pct],
            res[set == "accurate/moderate detectors only", pct],
            res[set == "accurate/moderate detectors only", pct_lo],
            res[set == "accurate/moderate detectors only", pct_hi]))
cat(sprintf("tables -> %s/20_within_taxon_chain.csv, 20_within_taxon_pertaxon.csv\n", cfg$out_dir))
cat("=========================================================================\n")
