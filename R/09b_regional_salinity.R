# =============================================================================
# Chain length along the R/V Svea surface salinity gradient, 2022-2026, with
# cruise-cluster bootstrap intervals (Fig. S5; Fig. 6d reads its contrasts).
#
# Reads:  data/SHARK_PlanktonImaging_<year>_SMHI/processed_data/data.txt and
#         data/roi_level/roi_<year>.txt.gz for 2022-2026 (missing years are dropped);
#         when neither is present, the bundled data/svea/ subset. Also
#         outputs/tangesund/25_log_contrast.csv and 27c_log_contrast_boot.csv
#         (optional; without it the Tangesund band uses the model SE).
# Writes: outputs/tangesund/09b_svea_chain_by_sample_taxon.rds (extraction cache),
#         09b_svea_samples.csv, 09b_svea_marine_bins.csv, 09b_svea_marine_cruises.csv,
#         09b_svea_curve.csv, 09b_svea_contrasts.csv, 09b_svea_contrasts_endpoint.csv,
#         09b_svea_pertaxon_marine.csv, 09b_svea_pooled.csv, 09b_svea_taxon_filter.csv,
#         09b_svea_summary.txt, 09b_svea_salinity_panel.{png,pdf} (Fig. S5). A run from
#         the SHARK deliveries also rewrites data/svea/.
# Env:    SVEA_FIG_ONLY=1 redraws Fig. S5 from the CSVs without refitting (used by 28).
# Deps:   data.table, mgcv, ggplot2, patchwork, ggtext, parallel.
# Run:    Rscript R/09b_regional_salinity.R   (30-60 min; the extraction is cached)
#         Run after 25 and 27c.
# =============================================================================
Sys.setlocale("LC_TIME", "C")  # English month names on date axes
suppressPackageStartupMessages({
  library(data.table); library(mgcv); library(ggplot2); library(patchwork); library(parallel); library(ggtext)
})

cfg <- list(
  roi_glob   = "data/roi_level/roi_%d.txt.gz",
  data_glob  = "data/SHARK_PlanktonImaging_%d_SMHI/processed_data/data.txt",
  deposit_dir  = "data/svea",                                  # bundled QC-pass subset
  cruise_fix   = "data/svea/svea_cruise_corrections.csv",       # sample -> resolved cruise_no
  deposit_meta = "data/svea/svea_samples.csv",
  deposit_agg  = "data/svea/svea_chain_by_sample_taxon.csv",
  years      = 2022:2026,
  out_dir    = "outputs/tangesund",
  junk       = c("Debris", "Beads", "Air_bubbles", "Zooplankton"),
  min_chain_roi = 20,   # chain-eligible ROIs per sample, as at Tangesund
  min_roi_taxon = 5,    # ROIs per sample x taxon row in the within-taxon models
  min_n_taxon   = 200,  # sample x taxon rows to keep a taxon
  min_cr_taxon  = 8,    # cruises to keep a taxon
  qc_exclude = list(c("2024-06-02", "2025-07-14"), c("2022-01-01", "2022-09-30")),  # IFCB134 PMT fault; commissioning
  full_lo = 6,  full_hi = 18,
  marine_min = 25, marine_lo = 25, marine_hi = 33,
  taxon_hi = 31, marine_mid = 28,   # within-taxon span and its coverage midpoint
  grid       = seq(3, 35, by = 0.5),
  tang_ref   = "outputs/tangesund/25_log_contrast.csv",        # Tangesund 3 m - 16 m log contrast
  tang_boot  = "outputs/tangesund/27c_log_contrast_boot.csv",  # its day-resampled interval
  B          = 1000,
  seed       = 42,
  cores      = min(8L, max(1L, detectCores() - 2L))
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
season_of <- function(month) factor(fcase(month %in% c(12, 1, 2), "Winter", month %in% 3:5, "Spring",
                                          month %in% 6:8, "Summer", month %in% 9:11, "Autumn"),
                                    levels = c("Winter", "Spring", "Summer", "Autumn"))
log_lines <- character(0)
say <- function(...) { s <- sprintf(...); cat(s, "\n"); log_lines <<- c(log_lines, s); invisible(s) }
# Tangesund reference on the log scale; point estimate from 25, interval from 27c
tang_ref <- fread(cfg$tang_ref)[set == "season-averaged"]
tang_log <- -tang_ref$contrast; tang_ci <- -c(tang_ref$hi, tang_ref$lo)
if (file.exists(cfg$tang_boot)) {
  tb <- fread(cfg$tang_boot)[set == "season-averaged"]
  tang_ci <- -c(tb$hi, tb$lo)
} else warning("27c_log_contrast_boot.csv missing; the Tangesund band uses the model SE")

# ---- sample metadata from the SHARK deliveries -------------------------------
meta <- rbindlist(lapply(cfg$years, function(yr) {
  f <- sprintf(cfg$data_glob, yr); if (!file.exists(f)) return(NULL)
  e <- fread(f, select = c("SMPNO", "CRUISE_NO", "IMAGE_VERIFICATION", "PSAL_FB", "SDATE"),
             na.strings = c("", "NA"), showProgress = FALSE)
  u <- unique(e[IMAGE_VERIFICATION == "PredictedByMachine",
                .(sample = SMPNO, cruise_no = as.character(CRUISE_NO),
                  salinity = as.numeric(PSAL_FB), date = as.Date(SDATE))])
  u[, year := yr][]
}), use.names = TRUE)
# Placeholder and missing 2026 cruise numbers, resolved from the 2026-09-07 delivery
if (nrow(meta) > 0 && file.exists(cfg$cruise_fix)) {
  fx <- fread(cfg$cruise_fix, colClasses = list(character = c("sample", "cruise_no")))
  meta[fx, on = "sample", cruise_no := i.cruise_no]
}
from_deposit <- nrow(meta) == 0 && file.exists(cfg$deposit_meta)
if (from_deposit) {                                   # fallback: bundled QC-pass subset
  meta <- fread(cfg$deposit_meta, colClasses = list(character = c("sample", "cruise_no")))
  meta[, date := as.Date(date)]
}
meta <- meta[is.finite(salinity) & !is.na(cruise_no) & !is.na(date)]
meta[, cruise := paste(year, cruise_no, sep = "_")]
meta[, season := season_of(as.integer(format(date, "%m")))]
meta[, qc_pass := TRUE]
for (p in cfg$qc_exclude) meta[date >= as.Date(p[1]) & date <= as.Date(p[2]), qc_pass := FALSE]
say("metadata: %d machine-classified samples with ferrybox salinity, %d cruises, %d QC-pass",
    nrow(meta), uniqueN(meta$cruise), sum(meta$qc_pass))
if (from_deposit) say("inputs: the bundled QC-pass subset (data/svea); pre-QC sample counts are not recoverable from it")

# ---- chain-eligible ROIs -> sample x taxon (same rule as R/03_chain_length.R) --
cache <- file.path(cfg$out_dir, "09b_svea_chain_by_sample_taxon.rds")
if (!file.exists(cache) && !any(file.exists(sprintf(cfg$roi_glob, cfg$years))) && file.exists(cfg$deposit_agg)) {
  agg <- fread(cfg$deposit_agg)                       # fallback: bundled sample x taxon table
  saveRDS(agg, cache)
} else if (!file.exists(cache)) {
  agg <- rbindlist(lapply(cfg$years, function(yr) {
    f <- sprintf(cfg$roi_glob, yr); if (!file.exists(f)) return(NULL)
    cat("  extracting", f, "\n")
    d <- fread(f, select = c("sample", "class", "cell_count", "cell_count_resolved", "biovolume_um3", "numBlobs"),
               showProgress = FALSE)
    d <- d[numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 & !class %in% cfg$junk]
    d <- d[!is.na(cell_count) & cell_count != -1]                                   # counter ran
    d[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1,
                          round(cell_count_resolved), 1)]
    d[, percell_bv := biovolume_um3 / nchain]
    d[, .(chain_len = mean(nchain), percell_bv = median(percell_bv),
          biovol = sum(biovolume_um3), n_roi = .N), by = .(sample, taxon = class)]
  }))
  saveRDS(agg, cache)
} else agg <- readRDS(cache)
setorder(agg, sample, taxon)                          # same row order from either source; pooled fits sum in row order
if (!from_deposit) {                                  # rewrite the bundled subset
  dir.create(cfg$deposit_dir, showWarnings = FALSE, recursive = TRUE)
  dm <- meta[qc_pass == TRUE, .(sample, cruise_no, date, salinity, year)][order(sample)]
  da <- agg[sample %in% dm$sample][order(sample, taxon)]
  g17 <- function(x) sprintf("%.17g", x)             # 17 digits so the fallback run is byte-identical
  fwrite(dm[, .(sample, cruise_no, date, salinity = g17(salinity), year)], cfg$deposit_meta)
  fwrite(da[, .(sample, taxon, chain_len = g17(chain_len), percell_bv = g17(percell_bv),
                biovol = g17(biovol), n_roi)], cfg$deposit_agg)
}

# ---- per-sample community chain length ---------------------------------------
cm <- agg[, .(chain_equal = sum(chain_len * n_roi) / sum(n_roi),      # = mean over ROIs
              chain_bv    = sum(chain_len * biovol) / sum(biovol),
              n_roi = sum(n_roi)), by = sample]
cm <- merge(cm, meta, by = "sample")
n_before_qc <- cm[n_roi >= cfg$min_chain_roi, .N]
cm <- cm[qc_pass == TRUE & n_roi >= cfg$min_chain_roi]
cm[, `:=`(cruise = factor(cruise), season = factor(season))]
fwrite(cm[, .(sample, cruise, date, year, season, salinity, chain_equal, chain_bv, n_roi)],
       file.path(cfg$out_dir, "09b_svea_samples.csv"))
say("community set: %d samples (%d before QC), %d cruises, %s chain-eligible ROIs, salinity %.1f-%.1f, years %d-%d",
    nrow(cm), n_before_qc, uniqueN(cm$cruise), format(sum(cm$n_roi), big.mark = ","),
    min(cm$salinity), max(cm$salinity), min(cm$year), max(cm$year))
say("community set by season: %s", paste(sprintf("%s %d", names(table(cm$season)), table(cm$season)), collapse = ", "))
say("community mean chain length: %.2f cells (equal weighting), %.2f (biovolume weighting)",
    mean(cm$chain_equal), mean(cm$chain_bv))
mar <- cm[salinity >= cfg$marine_min]
say("marine span (salinity >= %g): %d samples, %d cruises, %s ROIs, salinity %.1f-%.1f, mean chain length %.2f cells",
    cfg$marine_min, nrow(mar), uniqueN(mar$cruise), format(sum(mar$n_roi), big.mark = ","),
    min(mar$salinity), max(mar$salinity), mean(mar$chain_equal))
say("marine span by season: %s", paste(sprintf("%s %d", names(table(droplevels(mar)$season)), table(droplevels(mar)$season)), collapse = ", "))
# raw means by salinity bin and samples per cruise in the marine span
mar[, sbin := cut(salinity, c(25, 27, 29, 31, 33, 36), right = FALSE)]
fwrite(mar[, .(n = .N, cruises = uniqueN(cruise), mean_log_equal = mean(log(chain_equal)),
               mean_log_bv = mean(log(chain_bv)), mean_chain_equal = mean(chain_equal)), by = sbin][order(sbin)],
       file.path(cfg$out_dir, "09b_svea_marine_bins.csv"))
fwrite(mar[, .(n = .N, salinity_lo = min(salinity), salinity_hi = max(salinity), date_lo = min(date), date_hi = max(date)),
           by = cruise][order(-n)], file.path(cfg$out_dir, "09b_svea_marine_cruises.csv"))
mar[, sbin := NULL]

fig_only <- nzchar(Sys.getenv("SVEA_FIG_ONLY"))
if (fig_only) {
  contrasts <- fread(file.path(cfg$out_dir, "09b_svea_contrasts.csv"))
  pertaxon  <- fread(file.path(cfg$out_dir, "09b_svea_pertaxon_marine.csv"))
  curve     <- fread(file.path(cfg$out_dir, "09b_svea_curve.csv"))
} else {
# ---- bootstrap machinery -----------------------------------------------------
RNGkind("L'Ecuyer-CMRG"); set.seed(cfg$seed)
resample_cruises <- function(dat) {
  cr <- unique(dat$cruise); pick <- sample(cr, length(cr), replace = TRUE)
  rbindlist(lapply(seq_along(pick), function(i) {
    x <- dat[cruise == pick[i]]; x[, cruise := paste0("b", i)]; x }))
}
cluster_boot <- function(dat, fun, B = cfg$B) {
  # fun(dat) -> numeric vector; returns the observed vector and a B x length matrix of refits
  dat <- droplevels(as.data.table(dat))
  obs <- fun(dat)
  draws <- mclapply(seq_len(B), function(b) {
    dd <- resample_cruises(dat); dd[, cruise := factor(cruise)]
    r <- try(fun(dd), silent = TRUE)
    if (inherits(r, "try-error")) rep(NA_real_, length(obs)) else r
  }, mc.cores = cfg$cores, mc.set.seed = TRUE)
  list(obs = obs, draws = do.call(rbind, draws))
}
ci <- function(v) quantile(v, c(.025, .975), na.rm = TRUE)

# ---- (1) community: log chain length ~ s(salinity) + season + cruise RE ------
# returns the contrast followed by the partial on the grid (both from the same draws)
comm_fun <- function(resp, grid, lo, hi) function(dd) {
  dd <- copy(droplevels(dd)); dd[, y := log(get(resp))]
  f <- bam(y ~ s(salinity, k = 5) + season + s(cruise, bs = "re"), data = dd, discrete = TRUE, nthreads = 1)
  nd <- data.frame(salinity = c(lo, hi, grid), season = dd$season[1], cruise = dd$cruise[1])
  p <- as.numeric(predict(f, nd, type = "terms", terms = "s(salinity)")[, 1])
  c(p[2] - p[1], p[-(1:2)])
}
comm_sets <- list(
  list(label = "full gradient, colonies weighted equally",  resp = "chain_equal", dat = cm,  lo = cfg$full_lo,   hi = cfg$full_hi,   grid = cfg$grid),
  list(label = "full gradient, biovolume-weighted",         resp = "chain_bv",    dat = cm,  lo = cfg$full_lo,   hi = cfg$full_hi,   grid = numeric(0)),
  list(label = "marine span, colonies weighted equally",     resp = "chain_equal", dat = mar, lo = cfg$marine_lo, hi = cfg$marine_hi, grid = numeric(0)),
  list(label = "marine span, biovolume-weighted",            resp = "chain_bv",    dat = mar, lo = cfg$marine_lo, hi = cfg$marine_hi, grid = numeric(0)))
say("\n== (1) community mean chain length vs salinity (%d draws, %d cores) ==", cfg$B, cfg$cores)
contrasts <- list(); curve <- NULL
for (s in comm_sets) {
  t0 <- Sys.time()
  r <- cluster_boot(s$dat, comm_fun(s$resp, s$grid, s$lo, s$hi))
  q <- ci(r$draws[, 1])
  contrasts[[s$label]] <- data.table(set = s$label, lo_sal = s$lo, hi_sal = s$hi,
                                     contrast = r$obs[1], ci_lo = q[1], ci_hi = q[2],
                                     n = nrow(s$dat), cruises = uniqueN(s$dat$cruise),
                                     n_draws = sum(is.finite(r$draws[, 1])))
  say("%-42s n=%-5d cruises=%-3d S%g->S%g %+.3f log [%+.3f, %+.3f]  (%.1f min)", s$label, nrow(s$dat),
      uniqueN(s$dat$cruise), s$lo, s$hi, r$obs[1], q[1], q[2], as.numeric(Sys.time() - t0, units = "mins"))
  if (length(s$grid)) {
    g <- r$draws[, -1, drop = FALSE]; o <- r$obs[-1]
    ref <- which.min(abs(s$grid - cfg$full_lo))               # change from salinity 6
    curve <- data.table(salinity = s$grid, fit = o - o[ref],
                        lo = apply(sweep(g, 1, g[, ref]), 2, quantile, .025, na.rm = TRUE),
                        hi = apply(sweep(g, 1, g[, ref]), 2, quantile, .975, na.rm = TRUE))
    # trim to the 1st-99th percentile of salinity
    curve <- curve[salinity >= quantile(cm$salinity, .01) & salinity <= quantile(cm$salinity, .99)]
  }
}
fwrite(curve, file.path(cfg$out_dir, "09b_svea_curve.csv"))
# Endpoint sensitivity: the marine-span contrast to salinity 32. Kept in its own file
# because Figs. 6d and S5b plot every row of 09b_svea_contrasts.csv.
end_sets <- list(
  list(label = "marine span to 32, colonies weighted equally", resp = "chain_equal", dat = mar, lo = cfg$marine_lo, hi = 32, grid = numeric(0)),
  list(label = "marine span to 32, biovolume-weighted",        resp = "chain_bv",    dat = mar, lo = cfg$marine_lo, hi = 32, grid = numeric(0)))
endpoint <- rbindlist(lapply(end_sets, function(s) {
  t0 <- Sys.time()
  r <- cluster_boot(s$dat, comm_fun(s$resp, s$grid, s$lo, s$hi)); q <- ci(r$draws[, 1])
  say("%-42s n=%-5d cruises=%-3d S%g->S%g %+.3f log [%+.3f, %+.3f]  (%.1f min)", s$label, nrow(s$dat),
      uniqueN(s$dat$cruise), s$lo, s$hi, r$obs[1], q[1], q[2], as.numeric(Sys.time() - t0, units = "mins"))
  data.table(set = s$label, lo_sal = s$lo, hi_sal = s$hi, contrast = r$obs[1], ci_lo = q[1], ci_hi = q[2],
             n = nrow(s$dat), cruises = uniqueN(s$dat$cruise), n_draws = sum(is.finite(r$draws[, 1])))
}))
endpoint[, flat := ci_lo <= 0 & ci_hi >= 0]
fwrite(endpoint, file.path(cfg$out_dir, "09b_svea_contrasts_endpoint.csv"))

# ---- within-taxon set --------------------------------------------------------
d <- merge(agg[n_roi >= cfg$min_roi_taxon], meta[qc_pass == TRUE], by = "sample")
d <- d[is.finite(chain_len) & chain_len > 0]
d[, `:=`(lchain = log(chain_len), lbv = log(percell_bv))]
tf <- d[, .(rows_ge_min_roi = .N, cruises = uniqueN(cruise)), by = taxon]
tf <- merge(tf, merge(agg, meta[qc_pass == TRUE, .(sample)], by = "sample")[, .(rows_any = .N, n_roi_total = sum(n_roi)), by = taxon],
            by = "taxon", all = TRUE)
tf[is.na(rows_ge_min_roi), `:=`(rows_ge_min_roi = 0L, cruises = 0L)]
tf[, kept := rows_ge_min_roi >= cfg$min_n_taxon & cruises >= cfg$min_cr_taxon]
setorder(tf, -rows_ge_min_roi)
fwrite(tf, file.path(cfg$out_dir, "09b_svea_taxon_filter.csv"))      # inclusion rule per taxon
ag <- tf[taxon == "Asterionellopsis_glacialis"]
if (nrow(ag)) say("Asterionellopsis glacialis: %d sample x taxon rows (%d ROIs) in the QC-pass record, %d with >= %d ROIs, below the %d-row rule",
                  ag$rows_any, ag$n_roi_total, ag$rows_ge_min_roi, cfg$min_roi_taxon, cfg$min_n_taxon)
tk <- tf[kept == TRUE, taxon]
d  <- d[taxon %in% tk]
d[, `:=`(taxon = factor(taxon), cruise = factor(cruise), season = factor(season))]
say("\nwithin-taxon set: %d sample x taxon rows, %d taxa, %d cruises", nrow(d), length(tk), uniqueN(d$cruise))

# taxon and season fixed, cruise random; size_term adds s(log per-cell biovolume).
# size_term = FALSE matches the Tangesund per-taxon models, which have no size term.
taxon_fun <- function(pooled, lo, hi, size_term = TRUE) function(dd) {
  dd <- droplevels(dd)
  fm <- if (pooled && size_term)  lchain ~ s(salinity, k = 6) + s(lbv, k = 5) + taxon + season + s(cruise, bs = "re")
        else if (pooled)          lchain ~ s(salinity, k = 6) + taxon + season + s(cruise, bs = "re")
        else if (size_term)       lchain ~ s(salinity, k = 6) + s(lbv, k = 5) + season + s(cruise, bs = "re")
        else                      lchain ~ s(salinity, k = 6) + season + s(cruise, bs = "re")
  f <- bam(fm, data = dd, discrete = TRUE, nthreads = 1)
  nd <- data.frame(salinity = c(lo, hi), lbv = median(dd$lbv), season = dd$season[1], cruise = dd$cruise[1])
  if (pooled) nd$taxon <- dd$taxon[1]
  p <- as.numeric(predict(f, nd, type = "terms", terms = "s(salinity)")[, 1])
  p[2] - p[1]
}
# ---- (2) within-taxon contrasts within the marine span -----------------------
say("\n== (2) within-taxon salinity partials within the marine span, S%g -> S%g ==", cfg$marine_lo, cfg$taxon_hi)
dm <- droplevels(d[salinity >= cfg$marine_min])
dm[, `:=`(taxon = factor(taxon), cruise = factor(cruise), season = factor(season))]
say("marine-span within-taxon set: %d sample x taxon rows, %d taxa, %d cruises", nrow(dm), uniqueN(dm$taxon), uniqueN(dm$cruise))
covered_marine <- function(x) uniqueN(x$cruise) >= 6 && nrow(x) >= 100 &&
  x[salinity < cfg$marine_mid, .N] >= 20 && x[salinity > cfg$marine_mid, .N] >= 20
pertaxon <- rbindlist(lapply(levels(dm$taxon), function(tx) {
  x <- droplevels(dm[taxon == tx])
  if (!covered_marine(x)) return(NULL)
  r  <- cluster_boot(x, taxon_fun(FALSE, cfg$marine_lo, cfg$taxon_hi, size_term = FALSE))
  rf <- cluster_boot(x, taxon_fun(FALSE, cfg$marine_lo, cfg$taxon_hi, size_term = TRUE))    # at fixed cell size
  q <- ci(r$draws[, 1]); qf <- ci(rf$draws[, 1])
  say("  %-45s n=%-5d cruises=%-3d %+.3f [%+.3f, %+.3f]  %.1f%% of draws positive | fixed size %+.3f [%+.3f, %+.3f]",
      tx, nrow(x), uniqueN(x$cruise), r$obs, q[1], q[2], 100 * mean(r$draws[, 1] > 0, na.rm = TRUE), rf$obs, qf[1], qf[2])
  data.table(taxon = tx, n = nrow(x), cruises = uniqueN(x$cruise), contrast = r$obs, ci_lo = q[1], ci_hi = q[2],
             frac_pos = mean(r$draws[, 1] > 0, na.rm = TRUE), n_draws = sum(is.finite(r$draws[, 1])),
             contrast_fixed_size = rf$obs, ci_lo_fs = qf[1], ci_hi_fs = qf[2], n_draws_fs = sum(is.finite(rf$draws[, 1])))
}))
setorder(pertaxon, contrast)
fwrite(pertaxon, file.path(cfg$out_dir, "09b_svea_pertaxon_marine.csv"))
sk <- pertaxon[taxon == "Skeletonema_spp"]
if (nrow(sk)) contrasts[["Skeletonema, within taxon, marine span"]] <- data.table(
  set = "Skeletonema, within taxon, marine span", lo_sal = cfg$marine_lo, hi_sal = cfg$taxon_hi,
  contrast = sk$contrast, ci_lo = sk$ci_lo, ci_hi = sk$ci_hi, n = sk$n, cruises = sk$cruises, n_draws = sk$n_draws)

# ---- (3) full-gradient positive control --------------------------------------
say("\n== (3) full-gradient positive control, S%g -> S%g ==", cfg$full_lo, cfg$full_hi)
pooled <- rbindlist(lapply(list(
  list(label = "Skeletonema, within taxon", dat = d[taxon == "Skeletonema_spp"], pooled = FALSE),
  list(label = "all chain-formers, within taxon", dat = d, pooled = TRUE),
  list(label = "all chain-formers except Skeletonema, within taxon", dat = d[taxon != "Skeletonema_spp"], pooled = TRUE)),
  function(s) {
    r <- cluster_boot(s$dat, taxon_fun(s$pooled, cfg$full_lo, cfg$full_hi)); q <- ci(r$draws[, 1])
    say("  %-52s n=%-5d cruises=%-3d %+.3f [%+.3f, %+.3f]", s$label, nrow(s$dat), uniqueN(s$dat$cruise), r$obs, q[1], q[2])
    data.table(set = s$label, lo_sal = cfg$full_lo, hi_sal = cfg$full_hi, contrast = r$obs, ci_lo = q[1], ci_hi = q[2],
               n = nrow(s$dat), cruises = uniqueN(s$dat$cruise), n_draws = sum(is.finite(r$draws[, 1])))
  }))
fwrite(pooled, file.path(cfg$out_dir, "09b_svea_pooled.csv"))

contrasts <- rbindlist(contrasts)
contrasts[, flat := ci_lo <= 0 & ci_hi >= 0]
fwrite(contrasts, file.path(cfg$out_dir, "09b_svea_contrasts.csv"))

say("\nTangesund 3 m -> 16 m: %+.3f log units [%+.3f, %+.3f] (mean log per-sample chain length, day-of-year averaged; 25)",
    tang_log, tang_ci[1], tang_ci[2])
writeLines(log_lines, file.path(cfg$out_dir, "09b_svea_summary.txt"))
}  # end of the analysis block

# ---- Fig. S5 -----------------------------------------------------------------
pretty_taxon <- function(x) {
  x <- gsub("_", " ", x)
  x <- sub(" spp chain$", " spp. (chains)", x); x <- sub(" spp single$", " spp. (single)", x)
  x <- sub(" spp$", " spp.", x); x <- sub("^Chaetoceros cf convolutus$", "Chaetoceros cf. convolutus", x)
  # "Genus epithet-like" -> "Genus cf. epithet", as in data/functional_group_map.csv
  x <- sub("^([A-Z][a-z]+) ([a-z]+)-like", "\\1 cf. \\2", x)
  sub("^Leptocylindrus danicus Leptocylindrus minimus$", "Leptocylindrus danicus/minimus", x)
}
red <- "#d7301f"; blue <- "#2c7fb8"

pa <- ggplot(curve, aes(salinity, fit)) +
  annotate("rect", xmin = cfg$marine_min, xmax = max(curve$salinity), ymin = -Inf, ymax = Inf,
           fill = blue, alpha = 0.08) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey40", alpha = 0.2) +
  geom_line(colour = "grey20", linewidth = 0.9) +
  geom_rug(data = cm, aes(x = salinity, y = NULL), sides = "b", alpha = 0.04, length = unit(0.02, "npc")) +
  annotate("text", x = cfg$marine_min + 0.3, y = max(curve$hi), hjust = 0, vjust = 1, size = 3,
           label = "Marine span") +
  labs(x = "Salinity (ferrybox)", y = "Change in log community chain length\n(from salinity 6)",
       tag = "a", title = "Community chain length against surface salinity") +
  theme_bw(base_size = 11)

fl <- copy(contrasts)
fl[, label := factor(set, levels = rev(set))]
fl[, marine := grepl("marine", set)]
pb <- ggplot(fl, aes(contrast, label, colour = marine)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
  annotate("rect", xmin = tang_ci[1], xmax = tang_ci[2], ymin = -Inf, ymax = Inf, fill = red, alpha = 0.08) +
  geom_vline(xintercept = tang_log, linetype = "dashed", colour = red, linewidth = 0.6) +
  annotate("text", x = tang_log + 0.01, y = nrow(fl) + 0.45, label = "Tångesund\n3 m to 16 m", colour = red,
           size = 3, vjust = 1, hjust = 0) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0, linewidth = 0.5) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(`TRUE` = blue, `FALSE` = "grey30"), guide = "none") +
  scale_y_discrete(labels = function(x) {
    x <- paste0(toupper(substr(x, 1, 1)), substring(x, 2))
    # label spans by their salinity limits (the record spans salinity 3 to 35)
    x <- sub("^Full gradient", "Salinity 6 to 18", x)
    x <- sub("^Marine span", "Salinity 25 to 33", x)
    x <- sub("within taxon, marine span", "within taxon, salinity 25 to 31", x)   # as in Fig. 6d
    x <- sub("^Skeletonema", "*Skeletonema*", x)
    sub(", ", ",<br>", x)
  }, expand = expansion(add = c(0.6, 1.2))) +
  labs(x = "Change in log chain length\nover the stated span", y = NULL,
       tag = "b", title = "Contrasts (cruise-cluster 95% CI)") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_markdown())

# italicize a pretty_taxon() label, leaving spp., cf., /, (chains), (single) and -like upright
ital_taxon <- function(x) vapply(x, function(lab) {
  w <- strsplit(lab, " ", fixed = TRUE)[[1]]
  upright <- c("spp.", "sp.", "cf.", "/", "(chains)", "(single)")
  w <- ifelse(w %in% upright, w,
              ifelse(grepl("-like$", w), sub("^(.*)-like$", "*\\1*-like", w), paste0("*", w, "*")))
  paste(w, collapse = " ")
}, character(1), USE.NAMES = FALSE)

# Thalassiosira nordenskioeldii stays in the CSV but is not drawn (interval about -0.5 to +0.37)
pt <- copy(pertaxon)[!grepl("nordenskioeldii", taxon)]
pt[, taxon_lab := factor(pretty_taxon(taxon), levels = pretty_taxon(taxon))]
pt[, focal := taxon == "Skeletonema_spp"]
pc <- ggplot(pt, aes(contrast, taxon_lab, colour = focal)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0, linewidth = 0.5) +
  geom_point(aes(size = n)) +
  geom_point(aes(x = contrast_fixed_size), shape = 1, size = 2.2, colour = "grey30") +   # open: at fixed cell size
  scale_colour_manual(values = c(`TRUE` = red, `FALSE` = "grey45"), guide = "none") +
  scale_size_area(max_size = 3.2, guide = "none") +
  scale_y_discrete(labels = ital_taxon) +
  labs(x = sprintf("Change in log chain length within taxon, salinity %g to %g", cfg$marine_lo, cfg$taxon_hi), y = NULL,
       tag = "c", title = "Within-taxon change within the marine span") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_markdown())

fig <- (free(pa, type = "label", side = "l") | pb) / pc + plot_layout(heights = c(1, 1.15)) &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(file.path(cfg$out_dir, "09b_svea_salinity_panel.png"), fig, width = 11, height = 9, dpi = 200)
ggsave(file.path(cfg$out_dir, "09b_svea_salinity_panel.pdf"), fig, width = 11, height = 9)
cat("\nfigure -> outputs/tangesund/09b_svea_salinity_panel.{png,pdf}\n")

# Salinity at the minimum of the equal-weighted marine-span smooth (Supplementary Note 3).
# Last in the script so it does not shift the bootstrap RNG streams.
tryCatch(local({
  dd <- copy(droplevels(mar)); dd[, y := log(chain_equal)]; dd[, cruise := factor(cruise)]
  f  <- bam(y ~ s(salinity, k = 5) + season + s(cruise, bs = "re"), data = dd, discrete = TRUE, nthreads = 1)
  gr <- seq(cfg$marine_min, max(dd$salinity), by = 0.1)
  p  <- as.numeric(predict(f, data.frame(salinity = gr, season = dd$season[1], cruise = dd$cruise[1]),
                           type = "terms", terms = "s(salinity)")[, 1])
  say("marine-span equal-weighted smooth: minimum at salinity %.1f (grid %.1f-%.1f)", gr[which.min(p)], min(gr), max(gr))
}), error = function(e) message("marine-span minimum printout failed: ", conditionMessage(e)))
