# =============================================================================
# Chaetoceros at 3 and 16 m on the footprint size metric of the count validation
# (bounding-box area per cell): shares below the validation's size range and in its
# mid and large tertiles.
#
# Reads:  data/roi_level/roi_2016.txt.gz, outputs/tangesund/02_sample_meta.csv,
#         outputs/tangesund/21_chaeto_size_validation_raw.csv.
# Writes: outputs/tangesund/22_chaeto_size_reconcile.csv.
# Deps:   data.table.
# Run:    Rscript R/22_chaeto_size_reconcile.R   (after 21)
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table) })

cfg <- list(
  roi_file  = "data/roi_level/roi_2016.txt.gz",
  meta_csv  = "outputs/tangesund/02_sample_meta.csv",
  val_raw   = "outputs/tangesund/21_chaeto_size_validation_raw.csv",
  out_dir   = "outputs/tangesund",
  micron    = 1 / 2.77,                 # um per pixel
  depths    = c(3, 16)
)

meta <- fread(cfg$meta_csv)[qc_pass == TRUE, .(sample, depth)]
d <- fread(cfg$roi_file, select = c("sample","class_name","numBlobs","biovolume_um3","cell_count",
             "cell_count_resolved","BoundingBox_xwidth","BoundingBox_ywidth"), showProgress = FALSE)
d <- d[grepl("^Chaetoceros", class_name) & !is.na(cell_count) & cell_count != -1 &
       numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0]
d[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1, round(cell_count_resolved), 1L)]
d <- merge(d, meta, by = "sample")[depth %in% cfg$depths]

# the two per-cell size metrics on the same ROIs
d[, bv_esd := 2 * (3 / (4*pi) * (biovolume_um3 / nchain))^(1/3)]                       # biovolume ESD, as in script 16
d[, fp_esd := 2 * sqrt((BoundingBox_xwidth * BoundingBox_ywidth) / nchain / pi) * cfg$micron] # validation footprint ESD

# validation tertile thresholds, as 21's bins
v <- fread(cfg$val_raw)                                  # cols: manual, pred, err, percell_esd (footprint), ...
cat(sprintf("Spearman rho of count error vs short-axis width on the %d validation chains: %.3f (vs footprint ESD: %.3f)\n",
            nrow(v), cor(v$err, v$mindim_um, method = "spearman"), cor(v$err, v$percell_esd, method = "spearman")))
q <- quantile(v$percell_esd, c(1/3, 2/3))               # mid / large thresholds
val_support_min <- min(v$percell_esd)
big_bias <- v[percell_esd >= q[2], mean(err)]

prof <- d[, .(n = .N,
              bv_esd_med = round(median(bv_esd), 1),
              fp_esd_med = round(median(fp_esd), 1),
              mean_chain = round(mean(nchain), 2),
              pct_below_val_support = round(100 * mean(fp_esd < val_support_min), 1),
              pct_ge_mid   = round(100 * mean(fp_esd >= q[1]), 2),
              pct_ge_large = round(100 * mean(fp_esd >= q[2]), 2)), by = depth][order(depth)]
fwrite(prof, file.path(cfg$out_dir, "22_chaeto_size_reconcile.csv"))

cat(sprintf("validation footprint per-cell ESD: %.1f-%.1f um (median %.1f, n=%d) | tertiles mid>=%.1f large>=%.1f (large-bin bias %+.2f cells)\n",
            min(v$percell_esd), max(v$percell_esd), median(v$percell_esd), nrow(v), q[1], q[2], big_bias))
cat("\n== study Chaetoceros vs validation size support, by depth ==\n")
print(prof)

cat(sprintf("\ntable -> %s/22_chaeto_size_reconcile.csv\n", cfg$out_dir))
