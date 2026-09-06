# =============================================================================
# Chaetoceros depth gallery, stage 1: draw representative 3 m and 16 m ROIs (random
# within each depth's per-cell ESD interquartile range, one per day) and extract PNGs.
#
# Reads:  data/roi_level/roi_2016.txt.gz, outputs/tangesund/02_sample_meta.csv,
#         raw bins under data/raw_bins/2016/
# Writes: outputs/tangesund/chaeto_gallery_work/selection.csv and pngs/
# Env:    GALLERY_SEED (default 43, reproduces Fig. 4d), SINGLE_SEED (default 3)
# Deps:   data.table, iRfcb
# Run:    Rscript R/14_chaeto_gallery_select.R, then python3 R/15_chaeto_gallery.py
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(iRfcb) })
seed <- as.integer(Sys.getenv("GALLERY_SEED", "43"))
# single-cell strip (unused; Fig. S3 comes from 14b/15b) has its own seed
single_seed <- as.integer(Sys.getenv("SINGLE_SEED", "3"))
cat(sprintf("gallery selection seed = %d | single-cell seed = %d\n", seed, single_seed))

cfg <- list(
  roi_file  = "data/roi_level/roi_2016.txt.gz",
  meta_csv  = "outputs/tangesund/02_sample_meta.csv",
  raw_root  = "data/raw_bins/2016",   # raw .roi/.adc/.hdr bins, in per-day folders
  work_dir  = "outputs/tangesund/chaeto_gallery_work",     # selection.csv + extracted PNGs
  micron    = 1/2.77,     # project micron_factor (matches biovolume_um3); for the scale bar
  min_score = 0.9,        # confident Chaetoceros only
  k_chain   = 6,          # chains per layer
  k_single  = 4,          # single cells per layer
  max_chain = 8L          # display cap so one long chain does not shrink the common scale
)
dir.create(file.path(cfg$work_dir, "pngs"), recursive = TRUE, showWarnings = FALSE)

meta <- fread(cfg$meta_csv); meta[, date := as.Date(date)]
d <- fread(cfg$roi_file, select = c("sample","roi_number","class_name","biovolume_um3",
                                    "numBlobs","cell_count","cell_count_resolved","score"))
d <- d[grepl("^Chaetoceros_spp_chain", class_name) &
       is.finite(biovolume_um3) & biovolume_um3 > 0 & score >= cfg$min_score]
d[, ncell := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1,
                     round(cell_count_resolved), 1L)]
d <- merge(d, meta[qc_pass == TRUE, .(sample, depth, date)], by = "sample")
esd <- function(bv) 2 * (3 / (4 * pi) * bv)^(1/3)
d[, pc_esd := esd(biovolume_um3 / ncell)]

# random draw within the depth's per-cell-ESD IQR, one ROI per day
draw <- function(dep, single = FALSE, k = cfg$k_chain) {
  x <- d[depth == dep]
  x <- if (single) x[ncell == 1] else x[ncell >= 2 & ncell <= cfg$max_chain]
  qq <- quantile(x$pc_esd, c(.25, .75))
  b <- x[pc_esd >= qq[1] & pc_esd <= qq[2]]
  b <- b[b[, .I[sample(.N, 1)], by = date]$V1]     # one random ROI per sampling day
  b[sample(.N, min(k, .N))]
}
# chain seed reproduces Fig. 4d
set.seed(seed)
surf_chain <- draw(3)
deep_chain <- draw(16)
set.seed(single_seed)
single_surf <- draw(3,  TRUE, cfg$k_single)
single_deep <- draw(16, TRUE, cfg$k_single)
lays <- list(surface_3m        = surf_chain,
             deep_16m          = deep_chain,
             single_surface_3m = single_surf,
             single_deep_16m   = single_deep)
sel <- rbindlist(lapply(names(lays), function(n) cbind(layer = n, lays[[n]])))
sel <- sel[, .(layer, sample, roi_number, depth, date, ncell,
               pc_esd = round(pc_esd, 2), biovolume_um3 = round(biovolume_um3),
               score = round(score, 3))]
fwrite(sel, file.path(cfg$work_dir, "selection.csv"))

cat("== chosen set vs full population (representativeness) ==\n")
pop <- function(dep, single) { x <- d[depth == dep]
  x <- if (single) x[ncell == 1] else x[ncell >= 2]
  data.table(pop_ncell = round(mean(x$ncell), 2), pop_pcesd = round(mean(x$pc_esd), 2), pop_n = nrow(x)) }
for (n in names(lays)) { s <- lays[[n]]; dep <- s$depth[1]; sg <- grepl("single", n); p <- pop(dep, sg)
  cat(sprintf("  %-18s chosen ncell=%.2f pc_esd=%.2f (n=%d) | pop ncell=%.2f pc_esd=%.2f (n=%d)\n",
              n, mean(s$ncell), mean(s$pc_esd), nrow(s), p$pop_ncell, p$pop_pcesd, p$pop_n)) }

# no burnt-in scale bar: it corrupts the YOLO count (15 draws the bar)
outdir <- file.path(cfg$work_dir, "pngs")
for (s in unique(sel$sample)) {
  rf <- file.path(cfg$raw_root, substr(s, 1, 9), paste0(s, ".roi"))
  ifcb_extract_pngs(rf, out_folder = outdir, ROInumbers = sel[sample == s, roi_number],
                    scale_bar_um = NULL, gamma = 1, normalize = FALSE,
                    overwrite = TRUE, verbose = FALSE)
}
cat(sprintf("\nselection + %d PNGs -> %s\n",
            length(list.files(outdir, recursive = TRUE, pattern = "png$")), cfg$work_dir))
cat("next: python3 R/15_chaeto_gallery.py\n")
