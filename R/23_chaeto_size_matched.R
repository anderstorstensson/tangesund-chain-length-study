# =============================================================================
# Chaetoceros size coverage relative to the validation chains of 21, and the
# 3 m - 16 m contrast corrected for the width-dependent count bias inside that range.
#
# Reads:  data/roi_level/roi_2016.txt.gz, outputs/tangesund/02_sample_meta.csv,
#         outputs/tangesund/21_chaeto_size_validation_raw.csv
# Writes: outputs/tangesund/ 23_chaeto_coverage_by_depth.csv, 23_chaeto_bias_in_support.csv
# Deps:   data.table, mgcv
# Run:    Rscript R/23_chaeto_size_matched.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(mgcv) })

cfg <- list(
  roi_file = "data/roi_level/roi_2016.txt.gz",
  meta_csv = "outputs/tangesund/02_sample_meta.csv",
  val_raw  = "outputs/tangesund/21_chaeto_size_validation_raw.csv",
  out_dir  = "outputs/tangesund",
  micron   = 1 / 2.77,
  n_boot   = 2000
)

meta <- fread(cfg$meta_csv)[qc_pass == TRUE, .(sample, depth, date)]
d <- fread(cfg$roi_file, select = c("sample","class_name","numBlobs","biovolume_um3","cell_count",
             "cell_count_resolved","BoundingBox_xwidth","BoundingBox_ywidth"), showProgress = FALSE)
d <- d[grepl("^Chaetoceros", class_name) & !is.na(cell_count) & cell_count != -1 &
       numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0]
d[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1,
                      round(cell_count_resolved), 1L)]
d <- merge(d, meta, by = "sample")
d[, mindim := pmin(BoundingBox_xwidth, BoundingBox_ywidth) * cfg$micron]   # count-independent
d[, fp_esd := 2 * sqrt((BoundingBox_xwidth * BoundingBox_ywidth) / nchain / pi) * cfg$micron]
d <- d[is.finite(mindim) & is.finite(fp_esd)]

v <- fread(cfg$val_raw)                       # manual, pred, err, percell_esd (footprint), mindim_um
md_lo <- min(v$mindim_um);    md_hi <- max(v$mindim_um)
fp_lo <- min(v$percell_esd);  fp_hi <- max(v$percell_esd)
cat(sprintf("validation support: short axis %.1f-%.1f um | footprint per-cell ESD %.1f-%.1f um (n=%d)\n",
            md_lo, md_hi, fp_lo, fp_hi, nrow(v)))

# ---- (A) coverage by depth, both size axes -----------------------------------
cov <- d[, .(n = .N,
             mindim_med = round(median(mindim), 1),
             pct_below_mindim = round(100 * mean(mindim < md_lo), 1),
             fp_esd_med = round(median(fp_esd), 1),
             pct_below_fp = round(100 * mean(fp_esd < fp_lo), 1),
             mean_chain = round(mean(nchain), 2)), by = depth][order(depth)]
cat("\n== (A) Chaetoceros size coverage against the validation support, by depth ==\n")
print(cov)
fwrite(cov, file.path(cfg$out_dir, "23_chaeto_coverage_by_depth.csv"))

# ---- (B) measured bias inside the support ------------------------------------
# bias model on the short axis (not count-derived); footprint ESD version as a check
mb_md <- gam(err ~ s(mindim_um, k = 6), data = v, method = "REML")
mb_fp <- gam(err ~ s(percell_esd, k = 6), data = v, method = "REML")

sd2 <- d[depth %in% c(3, 16)]
insup <- sd2[mindim >= md_lo & mindim <= md_hi]
insup[, bias_md := as.numeric(predict(mb_md, newdata = data.frame(mindim_um = mindim)))]
insup[, bias_fp := as.numeric(predict(mb_fp, newdata = data.frame(percell_esd = fp_esd)))]
insup[, nchain_corr := nchain - bias_md]

bias_surf <- insup[depth == 3,  mean(bias_md)]; bias_deep <- insup[depth == 16, mean(bias_md)]
raw_in    <- insup[depth == 3,  mean(nchain)]      - insup[depth == 16, mean(nchain)]
corr_in   <- insup[depth == 3,  mean(nchain_corr)] - insup[depth == 16, mean(nchain_corr)]

set.seed(43)
bt <- replicate(cfg$n_boot, {
  dd <- insup[date %in% sample(unique(insup$date), uniqueN(insup$date), replace = TRUE)]
  c(dd[depth == 3, mean(nchain)]      - dd[depth == 16, mean(nchain)],
    dd[depth == 3, mean(nchain_corr)] - dd[depth == 16, mean(nchain_corr)],
    dd[depth == 3, mean(bias_md)]     - dd[depth == 16, mean(bias_md)])
})
ci <- apply(bt, 1, quantile, c(0.025, 0.975), na.rm = TRUE)

cat(sprintf("\n== (B) inside the validation support (n=%d: %d surface, %d deep) ==\n",
            nrow(insup), insup[depth == 3, .N], insup[depth == 16, .N]))
cat(sprintf("  mean predicted bias, surface: %+.3f cells | deep: %+.3f cells\n", bias_surf, bias_deep))
cat(sprintf("  DIFFERENTIAL bias (surface - deep): %+.3f cells [%.3f, %.3f]  <- the quantity of interest\n",
            bias_surf - bias_deep, ci[1,3], ci[2,3]))
cat(sprintf("  contrast uncorrected: %+.2f cells [%.2f, %.2f]\n", raw_in, ci[1,1], ci[2,1]))
cat(sprintf("  contrast bias-corrected: %+.2f cells [%.2f, %.2f]\n", corr_in, ci[1,2], ci[2,2]))
cat(sprintf("  check, footprint-based bias model: surface %+.3f, deep %+.3f, differential %+.3f\n",
            insup[depth == 3, mean(bias_fp)], insup[depth == 16, mean(bias_fp)],
            insup[depth == 3, mean(bias_fp)] - insup[depth == 16, mean(bias_fp)]))

out <- data.table(
  quantity = c("mean predicted bias, surface", "mean predicted bias, deep",
               "differential bias (surface - deep)",
               "contrast, uncorrected", "contrast, bias-corrected",
               "footprint model: mean predicted bias, surface", "footprint model: mean predicted bias, deep",
               "footprint model: differential bias (surface - deep)"),
  est = c(bias_surf, bias_deep, bias_surf - bias_deep, raw_in, corr_in,
          insup[depth == 3, mean(bias_fp)], insup[depth == 16, mean(bias_fp)],
          insup[depth == 3, mean(bias_fp)] - insup[depth == 16, mean(bias_fp)]),
  lo  = c(NA, NA, ci[1,3], ci[1,1], ci[1,2], NA, NA, NA),
  hi  = c(NA, NA, ci[2,3], ci[2,1], ci[2,2], NA, NA, NA))
out[, (c("est","lo","hi")) := lapply(.SD, round, 3), .SDcols = c("est","lo","hi")]
# the n behind every row
out[, `:=`(n = nrow(insup), n_surface = insup[depth == 3, .N], n_deep = insup[depth == 16, .N],
           n_days = uniqueN(insup$date))]
fwrite(out, file.path(cfg$out_dir, "23_chaeto_bias_in_support.csv"))

print(out)
