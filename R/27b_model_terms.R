# =============================================================================
# Model-term table (Data S1): edf, reference df, F or t, p, n, days and deviance
# explained for every cited model. Terms other scripts already write are read from
# their CSVs.
#
# Reads:  outputs/tangesund/02_sample_meta.csv, 03_chain_metrics.csv, 04_size_spectra.csv,
#         05_composition.csv, 06_env_matched.csv, 08_diel_terms.csv, 10_aic_comparison.csv,
#         10_joint_model_terms.csv, 10_depth_in_hydro_models.csv, 10_conditioned_slope.csv,
#         10_mixed_layer_contrast.csv, 10_salinity_bin_depth_slopes.csv,
#         10b_monthly_tests.csv, 16_percell_metrics.csv, 16_roi_correlation.csv,
#         16b_within_depth_correlation.csv, 20_within_taxon_chain.csv;
#         data/roi_level/roi_2016.txt.gz, data/functional_group_map.csv
# Writes: outputs/tangesund/27b_model_terms.csv
# Deps:   data.table, mgcv
# Run:    Rscript R/27b_model_terms.R   (after 03-06, 08, 10, 10b, 13, 16, 16b, 20)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(mgcv) })
cfg <- list(out_dir = "outputs/tangesund", roi_file = "data/roi_level/roi_2016.txt.gz",
            meta_csv = "outputs/tangesund/02_sample_meta.csv",
            fg_map = "data/functional_group_map.csv",
            junk = c("Debris", "Beads", "Air_bubbles", "Zooplankton"), depths = c(3, 6, 8, 11, 13, 16))
o <- function(f) file.path(cfg$out_dir, f)
rows <- list()
add_smooths <- function(block, model, response, unit, m, terms = NULL) {
  ss <- summary(m); st <- ss$s.table
  if (is.null(terms)) terms <- rownames(st)
  terms <- terms[terms %in% rownames(st)]
  n_days <- tryCatch(length(levels(m$model$day)), error = function(e) NA_integer_)
  rows[[length(rows) + 1]] <<- data.table(
    block = block, model = model, response = response, unit = unit, n = nrow(m$model), n_days = n_days,
    dev_expl = ss$dev.expl, term = terms, edf = st[terms, "edf"], ref_df = st[terms, "Ref.df"],
    F = st[terms, "F"], t = NA_real_, p = st[terms, "p-value"], estimate = NA_real_, lo = NA_real_, hi = NA_real_)
}
add_row <- function(...) rows[[length(rows) + 1]] <<- data.table(...)
# depth-and-season GAMM shared by 03, 04, 05 and 16
gamm_form <- function(yv) reformulate(c("s(depth,k=5)", "s(yday,k=8)", "ti(depth,yday,k=c(5,6))", "s(day,bs=\"re\")"), yv)

# ---- 1. sample-level depth-and-season models --------------------------------------
chn <- fread(o("03_chain_metrics.csv"))[qc_pass == TRUE & reliable_chain == TRUE]; chn[, day := factor(date)]
m1 <- bam(gamm_form("mean_chain"), data = chn, method = "fREML")
add_smooths("Vertical decline in chain length", "Chain-length GAMM (Fig. 2; Methods, Depth and season models)",
            "mean chain length (cells)", "sample", m1)
pc <- fread(o("16_percell_metrics.csv")); pc[, day := factor(date)]
m2 <- bam(gamm_form("log_percell"), data = pc, method = "fREML")
add_smooths("Parallel vertical size gradient", "Per-cell size GAMM (Fig. 4a)", "log10 per-cell biovolume", "sample", m2)
sp <- fread(o("04_size_spectra.csv"))[qc_pass == TRUE & reliable == TRUE]; sp[, day := factor(date)]
for (v in c("b", "size_div", "median_esd")) {
  m <- bam(gamm_form(v), data = sp, method = "fREML")
  add_smooths("Community size spectrum", "Size-spectrum GAMM",
              c(b = "size-spectrum exponent b", size_div = "size diversity", median_esd = "median community ESD (um)")[v],
              "sample", m, "s(depth)")
}
cm <- fread(o("05_composition.csv"))[qc_pass == TRUE]; cm[, `:=`(day = factor(date), cd_share = `f_chain diatom`)]
for (v in c("rich_rare100", "shannon", "cd_share")) {
  m <- bam(gamm_form(v), data = cm[is.finite(get(v))], method = "fREML")
  add_smooths("Layered community (Fig. 5)", "Composition GAMM (Fig. 5)",
              c(rich_rare100 = "rarefied richness (100 cells)", shannon = "Shannon diversity", cd_share = "chain-diatom cell share")[v],
              "sample", m, "s(depth)")
}

# ---- 2. ROI-level within-taxon per-cell size (16 specification) ------------------
meta <- fread(cfg$meta_csv); meta[, date := as.Date(date)]
d <- fread(cfg$roi_file, select = c("sample", "class_name", "biovolume_um3", "numBlobs", "cell_count", "cell_count_resolved"),
           showProgress = FALSE)
d <- d[numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 & !class_name %in% cfg$junk]
d <- d[!is.na(cell_count) & cell_count != -1]
d[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1, round(cell_count_resolved), 1L)]
d[, percell_bv := biovolume_um3 / nchain]
d <- merge(d, meta[, .(sample, depth, yday, date, qc_pass)], by = "sample")[qc_pass == TRUE & depth %in% cfg$depths & percell_bv > 0]
top <- d[, .N, by = class_name][order(-N)][1:6, class_name]
cet <- d[class_name %in% top][, `:=`(taxon = factor(class_name), day = factor(date), samp = factor(sample))]
m3 <- bam(log10(percell_bv) ~ taxon + s(depth, k = 5) + s(yday, k = 6) + s(day, bs = "re") + s(samp, bs = "re"),
          data = cet, method = "fREML", discrete = TRUE, nthreads = 4)
add_smooths("Parallel vertical size gradient", "Within-taxon per-cell size, ROI-level model with taxon factor (Fig. 4c)",
            "log10 per-cell biovolume", "ROI", m3, "s(depth)")
rm(d, cet); invisible(gc())

# ---- 3. within-taxon chain length, ROI level (20 CSV) -----------------------------
wt <- fread(o("20_within_taxon_chain.csv"))
for (i in seq_len(nrow(wt))) add_row(
  block = "Within-taxon change", model = sprintf("Within-taxon chain length, ROI-level model, %s (Table 2)", wt$set[i]),
  response = "log chain length", unit = "ROI", n = wt$n_roi[i], n_days = wt$n_days[i], dev_expl = wt$dev_expl[i],
  term = "s(depth)", edf = wt$edf[i], ref_df = if ("ref_df" %in% names(wt)) wt$ref_df[i] else NA_real_,
  F = wt$F[i], t = NA_real_, p = wt$p[i],
  estimate = wt$est[i], lo = wt$lo[i], hi = wt$hi[i])

# ---- 4. hydrography and fragmentation models (10 / 13 specifications) -------------
env <- fread(o("06_env_matched.csv"))
hy <- merge(env, fread(o("03_chain_metrics.csv"))[, .(sample, qc_pass, reliable_chain, mean_chain)], by = "sample")
hy <- hy[qc_pass == TRUE & reliable_chain == TRUE]; hy[, day := factor(date)]
ha <- hy[is.finite(insitu_temp) & is.finite(mld)]
mA <- bam(mean_chain ~ s(insitu_temp, k = 6) + s(insitu_sal, k = 6) + s(dsigma, k = 6) + s(mld, k = 6) + s(yday, k = 8) +
            ti(insitu_sal, yday, k = c(5, 6)) + s(day, bs = "re"), data = ha, method = "fREML")   # as 10; no shrinkage
add_smooths("Halocline as the vertical axis", "Hydrography GAMM, depth replaced by matched in situ variables (Fig. 6)",
            "mean chain length (cells)", "sample", mA, c("s(insitu_temp)", "s(insitu_sal)", "s(dsigma)", "s(mld)", "s(yday)", "ti(insitu_sal,yday)"))
db <- hy[is.finite(insitu_cur)]
mB <- bam(mean_chain ~ s(insitu_cur, k = 6) + s(insitu_sal, k = 6) + s(yday, k = 8) + s(day, bs = "re"),
          data = db, method = "fREML")   # as 13; no shrinkage
add_smooths("Fragmentation test", "Current-speed GAMM (Results, bulk flow)", "mean chain length (cells)", "sample", mB,
            c("s(insitu_cur)", "s(insitu_sal)", "s(yday)"))
jt <- fread(o("10_joint_model_terms.csv"))       # depth + salinity + temperature together
for (i in seq_len(nrow(jt))) if (jt$term[i] != "s(day)") add_row(block = "Halocline as the vertical axis",
  model = "Joint GAMM, depth and matched salinity and temperature together (concurve smooths; Results)",
  response = "mean chain length (cells)", unit = "sample", n = jt$n[i], n_days = jt$n_days[i], dev_expl = jt$dev_expl[i],
  term = jt$term[i], edf = jt$edf[i], ref_df = jt$ref_df[i], F = jt$F[i], t = NA_real_, p = jt$p[i],
  estimate = NA_real_, lo = NA_real_, hi = NA_real_)
# depth inside the hydrography models, and the unbinned conditioned slope
di <- fread(o("10_depth_in_hydro_models.csv"))
for (i in seq_len(nrow(di))) if (!grepl("joint", di$model[i])) add_row(block = "Halocline as the vertical axis",
  model = sprintf("GAMM, %s (s(depth) contrast 3 m minus 16 m at median hydrography; Supplementary Note 2)", di$model[i]),
  response = "mean chain length (cells)", unit = "sample", n = di$n[i], n_days = di$n_days[i], dev_expl = di$dev_expl[i],
  term = di$term[i], edf = di$edf[i], ref_df = di$ref_df[i], F = di$F[i], t = NA_real_, p = di$p[i],
  estimate = di$contrast_3_16_at_median_hydrography[i], lo = di$lo[i], hi = di$hi[i])
cs <- fread(o("10_conditioned_slope.csv"))
for (i in seq_len(nrow(cs))) add_row(block = "Depth slope conditioned on matched salinity (Supplementary Note 2)",
  model = sprintf("%s, %s", cs$model[i], cs$subset[i]), response = "mean chain length (cells per m)", unit = "sample",
  n = cs$n[i], n_days = cs$n_days[i], dev_expl = NA_real_, term = "depth (parametric slope)", edf = NA_real_, ref_df = NA_real_,
  F = NA_real_, t = cs$t[i], p = cs$p[i], estimate = cs$slope[i], lo = cs$lo[i], hi = cs$hi[i])
aic <- fread(o("10_aic_comparison.csv"))
for (i in seq_len(nrow(aic))) add_row(block = "Halocline as the vertical axis", model = sprintf("AIC, %s", aic$model[i]),
  response = "mean chain length (cells)", unit = "sample", n = nrow(ha), n_days = uniqueN(ha$day), dev_expl = NA_real_,
  term = "AIC", edf = NA_real_, ref_df = NA_real_, F = NA_real_, t = NA_real_, p = NA_real_, estimate = aic$aic[i], lo = NA_real_, hi = NA_real_)
ml <- fread(o("10_mixed_layer_contrast.csv"))
for (i in seq_len(nrow(ml))) add_row(block = "Mixed-layer contrast",
  model = sprintf("Mixed-layer versus below-pycnocline contrast, %s", ml$model[i]), response = "mean chain length (cells)",
  unit = "sample", n = NA_integer_, n_days = NA_integer_, dev_expl = NA_real_, term = "in mixed layer (parametric)",
  edf = NA_real_, ref_df = NA_real_, F = NA_real_, t = ml$estimate[i] / ml$se[i], p = ml$p[i],
  estimate = ml$estimate[i], lo = ml$estimate[i] - 1.96 * ml$se[i], hi = ml$estimate[i] + 1.96 * ml$se[i])
sb <- fread(o("10_salinity_bin_depth_slopes.csv"))
for (i in seq_len(nrow(sb))) add_row(block = "Depth slope within salinity bins (Fig. 6c, Note S2)",
  model = sprintf("%s, bin %s", sb$scheme[i], sb$bin[i]), response = "mean chain length (cells per m)", unit = "sample",
  n = sb$n[i], n_days = sb$n_days[i], dev_expl = NA_real_, term = "depth (parametric slope)", edf = NA_real_, ref_df = NA_real_,
  F = NA_real_, t = sb$slope[i] / sb$se[i], p = sb$p[i], estimate = sb$slope[i], lo = sb$lo[i], hi = sb$hi[i])

# ---- 5. diel (08), October (10b), correlations (16b) --------------------------------
dl <- fread(o("08_diel_terms.csv"))
for (i in seq_len(nrow(dl))) add_row(block = "Diel cycle (Fig. 7)", model = sprintf("Diel GAMM, %s", dl$model[i]),
  response = "mean chain length (cells)", unit = "sample", n = dl$n[i], n_days = dl$n_days[i], dev_expl = dl$dev_expl[i],
  term = dl$term[i], edf = dl$edf[i], ref_df = dl$ref_df[i], F = dl$F[i], t = if ("t" %in% names(dl)) dl$t[i] else NA_real_,
  p = dl$p[i], estimate = NA_real_, lo = NA_real_, hi = NA_real_)
oc <- fread(o("10b_monthly_tests.csv"))
for (i in seq_len(nrow(oc))) add_row(block = "October halocline shoaling (Fig. S4c, Note S2)",
  model = sprintf("%s, %s", oc$month[i], oc$test[i]), response = "mean chain length (cells)", unit = "sample",
  n = oc$n[i], n_days = oc$n_days[i], dev_expl = NA_real_, term = oc$test[i], edf = NA_real_, ref_df = oc$df[i],
  F = ifelse(grepl("F test", oc$test[i]), oc$statistic[i], NA_real_),
  t = ifelse(grepl("F test", oc$test[i]), NA_real_, oc$statistic[i]), p = oc$p[i],
  estimate = oc$estimate[i], lo = oc$lo[i], hi = oc$hi[i])
rc <- fread(o("16_roi_correlation.csv"))          # ROI-level r
for (i in seq_len(nrow(rc))) add_row(block = "Chain length and cell size (Fig. 4b)",
  model = sprintf("%s correlation, ROI level, all chain-counted ROIs", rc$method[i]), response = "log10 per-cell biovolume vs cells per chain",
  unit = "ROI", n = rc$n[i], n_days = rc$n_days[i], dev_expl = NA_real_, term = "r (day-clustered bootstrap)", edf = NA_real_,
  ref_df = NA_real_, F = NA_real_, t = NA_real_, p = NA_real_, estimate = rc$r[i], lo = rc$ci_lo[i], hi = rc$ci_hi[i])
wd <- fread(o("16b_within_depth_correlation.csv"))
for (i in seq_len(nrow(wd))) add_row(block = "Chain length and cell size (Fig. 4b)",
  model = sprintf("Pearson correlation, %s", wd$correlation[i]), response = "per-cell ESD vs mean chain length", unit = "sample",
  n = wd$n[i], n_days = wd$n_days[i], dev_expl = NA_real_, term = "r (day-clustered bootstrap)", edf = NA_real_,
  ref_df = NA_real_, F = NA_real_, t = NA_real_, p = wd$p_day_bootstrap[i], estimate = wd$r[i], lo = wd$ci_lo[i], hi = wd$ci_hi[i])

out <- rbindlist(rows, use.names = TRUE)
fwrite(out, o("27b_model_terms.csv"))
cat("-- model terms (27b_model_terms.csv) --\n")
print(out[, .(block, model = substr(model, 1, 60), term, n, n_days, dev = round(100 * dev_expl, 1), edf = round(edf, 2),
              ref_df = round(ref_df, 2), F = round(F, 1), t = round(t, 2), p = signif(p, 2),
              est = round(estimate, 3))], nrows = 200)
