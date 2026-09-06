# =============================================================================
# Instrument size limits and surface light: (A) chain ROIs longer than the intake
# screen or spanning the image frame, and the depth profile with only ROIs of
# 150 um or less; (B) concentration of chains by length class and depth; (C) the
# 3 m minus 16 m difference in per-cell ESD of the six dominant chain-formers;
# (D) daily modeled surface PAR (STRANG, SMHI) added to the depth and season model
# of mean chain length (Supplementary Note 1, Methods).
#
# Reads:  data/roi_level/roi_2016.txt.gz, data/functional_group_map.csv,
#         outputs/tangesund/02_sample_meta.csv, 02c_ifcb_settings.csv, 03_chain_metrics.csv,
#         data/strang/strang_par_daily_2016.json (downloaded from the SMHI open-data
#         API and cached on the first run).
# Writes: outputs/tangesund/27e_truncation_by_depth.csv, 27e_chain_concentration_by_depth.csv,
#         27e_percell_esd_contrast_by_taxon.csv, 27e_light_model_terms.csv,
#         27e_light_day_level.csv.
# Deps:   data.table, mgcv, jsonlite.
# Run:    Rscript R/27e_size_limits_light.R
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(mgcv); library(jsonlite) })
cfg <- list(
  roi_file = "data/roi_level/roi_2016.txt.gz",
  meta_csv = "outputs/tangesund/02_sample_meta.csv",
  vol_csv  = "outputs/tangesund/02c_ifcb_settings.csv",
  chain_csv = "outputs/tangesund/03_chain_metrics.csv",
  fg_map   = "data/functional_group_map.csv",
  strang_par = "data/strang/strang_par_daily_2016.json",
  # STRANG point series at the station: parameter 120 is photosynthetically active
  # radiation, daily sums in Wh m-2
  strang_url = paste0("https://opendata-download-metanalys.smhi.se/api/category/strang1g/version/1/",
                      "geotype/point/lon/11.493/lat/58.075/parameter/120/data.json",
                      "?from=2016-08-01&to=2016-10-20&interval=daily"),
  out_dir  = "outputs/tangesund",
  junk     = c("Debris", "Beads", "Air_bubbles", "Zooplankton"),
  depths   = c(3, 6, 8, 11, 13, 16),
  px_per_um = 2.77, screen_um = 150, min_chain_roi = 20, n_boot = 2000
)
if (!file.exists(cfg$strang_par)) {
  dir.create(dirname(cfg$strang_par), recursive = TRUE, showWarnings = FALSE)
  tryCatch(download.file(cfg$strang_url, cfg$strang_par, quiet = TRUE),
           error = function(e) stop("could not download the STRANG series from SMHI: ", conditionMessage(e)))
}
meta <- fread(cfg$meta_csv); meta[, date := as.Date(date)]
vol  <- fread(cfg$vol_csv)[, .(sample, ml_analyzed)]
d <- fread(cfg$roi_file, select = c("sample", "class_name", "biovolume_um3", "numBlobs", "cell_count",
                                    "cell_count_resolved", "maxFeretDiameter", "BoundingBox_xwidth",
                                    "BoundingBox_ywidth"), showProgress = FALSE)
# image frame = the largest bounding box in the deployment, in each direction
frame_px <- c(x = max(d$BoundingBox_xwidth, na.rm = TRUE), y = max(d$BoundingBox_ywidth, na.rm = TRUE))
d <- d[numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 & !class_name %in% cfg$junk]
d[, chain_elig := !is.na(cell_count) & cell_count != -1]
d[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1, round(cell_count_resolved), 1)]
ce <- merge(d[chain_elig == TRUE], meta[, .(sample, depth, yday, date, qc_pass)], by = "sample")
ce <- merge(ce[qc_pass == TRUE & depth %in% cfg$depths], vol, by = "sample")
ce[, feret_um := maxFeretDiameter / cfg$px_per_um]
ce[, full_frame := BoundingBox_xwidth >= 0.95 * frame_px["x"] | BoundingBox_ywidth >= 0.95 * frame_px["y"]]

# ---- (A) truncation ----------------------------------------------------------
ps <- ce[, .(n = .N, mc = mean(nchain), mc_le = mean(nchain[feret_um <= cfg$screen_um])),
         by = .(sample, depth)][n >= cfg$min_chain_roi]
A <- merge(
  ce[, .(n_roi = .N, feret_median_um = median(feret_um), feret_p99_um = quantile(feret_um, .99),
         feret_max_um = max(feret_um), pct_gt_screen = 100 * mean(feret_um > cfg$screen_um),
         pct_full_frame = 100 * mean(full_frame)), by = depth],
  ps[, .(mean_chain = mean(mc), mean_chain_le_screen = mean(mc_le)), by = depth], by = "depth")[order(depth)]
A[, `:=`(frame_x_um = frame_px["x"] / cfg$px_per_um, frame_y_um = frame_px["y"] / cfg$px_per_um)]
fwrite(A, file.path(cfg$out_dir, "27e_truncation_by_depth.csv"))
cat("\n== (A) truncation ==\n"); print(A)

# ---- (B) concentration by chain-length class ---------------------------------
allqc <- merge(meta[qc_pass == TRUE & depth %in% cfg$depths, .(sample, depth)], vol, by = "sample")
ce[, len_class := cut(nchain, c(0, 1, 3, 7, Inf), labels = c("1", "2-3", "4-7", "8+"))]
cc <- ce[, .(n = .N, cells = sum(nchain), bv = sum(biovolume_um3)), by = .(sample, len_class)]
cc <- merge(CJ(sample = allqc$sample, len_class = levels(ce$len_class)), cc, by = c("sample", "len_class"), all.x = TRUE)
cc[is.na(n), `:=`(n = 0L, cells = 0, bv = 0)]
cc <- merge(cc, allqc, by = "sample")
B <- cc[, .(roi_per_ml = mean(n / ml_analyzed), cells_per_ml = mean(cells / ml_analyzed),
            biovolume_um3_per_ml = mean(bv / ml_analyzed), n_samples = .N), by = .(depth, len_class)][order(depth, len_class)]
fwrite(B, file.path(cfg$out_dir, "27e_chain_concentration_by_depth.csv"))
cat("\n== (B) chain ROIs per mL by length class ==\n"); print(dcast(B, depth ~ len_class, value.var = "roi_per_ml"))

# ---- (C) per-taxon per-cell ESD contrast --------------------------------------
esd <- function(bv) 2 * (3 / (4 * pi) * bv)^(1/3)
disp <- fread(cfg$fg_map)[, setNames(display_name, class)]
ce[, `:=`(taxon = disp[class_name], pc = biovolume_um3 / nchain)]
top <- ce[, .N, by = taxon][order(-N)][1:6, taxon]
set.seed(1)
C <- rbindlist(lapply(top, function(tx) {
  x <- ce[taxon == tx & depth %in% c(3, 16)]; days <- unique(x$date); idx <- split(seq_len(nrow(x)), x$date)
  f <- function(z) esd(mean(z[depth == 3, pc])) - esd(mean(z[depth == 16, pc]))
  bt <- vapply(seq_len(cfg$n_boot), function(b) {
    z <- x[unlist(idx[as.character(sample(days, length(days), replace = TRUE))], use.names = FALSE)]
    if (z[depth == 3, .N] == 0 || z[depth == 16, .N] == 0) NA_real_ else f(z)
  }, numeric(1))
  data.table(taxon = tx, esd_3m = esd(mean(x[depth == 3, pc])), esd_16m = esd(mean(x[depth == 16, pc])),
             diff = f(x), ci_lo = quantile(bt, .025, na.rm = TRUE), ci_hi = quantile(bt, .975, na.rm = TRUE),
             n_days = length(days))
}))
fwrite(C, file.path(cfg$out_dir, "27e_percell_esd_contrast_by_taxon.csv"))
cat("\n== (C) per-cell ESD, 3 m minus 16 m ==\n"); print(C)

# ---- (D) light ---------------------------------------------------------------
# daily PAR in kWh m-2: same day, day before, 3- and 7-day trailing means
L <- as.data.table(fromJSON(cfg$strang_par))[, .(date = as.Date(substr(date_time, 1, 10)), par = value / 1000)][order(date)]
L[, `:=`(par_prev = shift(par), par3 = frollmean(par, 3, align = "right"), par7 = frollmean(par, 7, align = "right"))]
s <- fread(cfg$chain_csv)[reliable_chain == TRUE & depth %in% cfg$depths]; s[, date := as.Date(date)]
s <- merge(s, L, by = "date"); s[, day := factor(date)]
base <- mean_chain ~ s(depth, k = 5) + s(yday, k = 8) + ti(depth, yday, k = c(5, 6)) + s(day, bs = "re")
m0 <- bam(base, data = s, method = "ML")
D <- rbindlist(lapply(c("par", "par_prev", "par3", "par7"), function(v) {
  s[, x := get(v)]
  m1 <- bam(update(base, . ~ . + x), data = s, method = "ML")
  m2 <- bam(update(base, . ~ . + x + x:depth), data = s, method = "ML")
  p1 <- summary(m1)$p.table; p2 <- summary(m2)$p.table
  data.table(window = v, slope_cells_per_kWh = p1["x", 1], se = p1["x", 2], p = p1["x", 4], dAIC = AIC(m1) - AIC(m0),
             interaction_per_m = p2["x:depth", 1], p_interaction = p2["x:depth", 4], dAIC_interaction = AIC(m2) - AIC(m0))
}))
fwrite(D, file.path(cfg$out_dir, "27e_light_model_terms.csv"))
cat("\n== (D) daily PAR (kWh m-2) added to the depth x season model ==\n"); print(D)
# day-level correlations with the seasonal trend removed from both series
dm <- s[, .(mc = mean(mean_chain)), by = .(date, yday, depth, par, par3)]
w <- dcast(dm[depth %in% c(3, 16)], date + yday + par + par3 ~ depth, value.var = "mc")
setnames(w, c("3", "16"), c("chain_3m", "chain_16m")); w[, contrast := chain_3m - chain_16m]
D2 <- rbindlist(lapply(c("chain_3m", "chain_16m", "contrast"), function(y) rbindlist(lapply(c("par", "par3"), function(v) {
  z <- w[is.finite(get(y)) & is.finite(get(v))]
  ct <- cor.test(resid(gam(as.formula(paste(y, "~ s(yday, k = 6)")), data = z)),
                 resid(gam(as.formula(paste(v, "~ s(yday, k = 6)")), data = z)))
  data.table(response = y, window = v, r = ct$estimate, ci_lo = ct$conf.int[1], ci_hi = ct$conf.int[2], p = ct$p.value, n_days = nrow(z))
}))))
fwrite(D2, file.path(cfg$out_dir, "27e_light_day_level.csv"))
cat("\n== (D) day-level correlations, seasonal trend removed ==\n"); print(D2)
