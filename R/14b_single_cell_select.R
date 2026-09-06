# =============================================================================
# Single-cell Chaetoceros exemplars (Fig. S3), stage 1: candidate selection and PNGs.
# Rule at 3 m and 16 m: Chaetoceros chain or single class, score >= 0.9, count 1,
# numBlobs == 1, Area fills >= 25% of the bounding box, per-cell ESD inside the depth's
# IQR, one per day, 40 candidates in random order (seed 429). 15b keeps the first four
# per depth with exactly one detector box.
#
# Reads:  data/roi_level/roi_2016.txt.gz, outputs/tangesund/02_sample_meta.csv,
#         raw bins under data/raw_bins/2016/
# Writes: outputs/tangesund/single_cell_work/selection.csv,
#         14b_single_cell_population.csv, pngs/
# Deps:   data.table, iRfcb
# Run:    Rscript R/14b_single_cell_select.R, then python3 R/15b_single_cell_gallery.py
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(iRfcb) })
cfg <- list(
  roi_file = "data/roi_level/roi_2016.txt.gz",
  meta_csv = "outputs/tangesund/02_sample_meta.csv",
  raw_root = "data/raw_bins/2016",
  work_dir = "outputs/tangesund/single_cell_work",
  min_score = 0.9, min_fill = 0.25, k_candidates = 40, seed = 429)
dir.create(file.path(cfg$work_dir, "pngs"), recursive = TRUE, showWarnings = FALSE)
meta <- fread(cfg$meta_csv); meta[, date := as.Date(date)]
d <- fread(cfg$roi_file, select = c("sample", "roi_number", "class_name", "biovolume_um3", "numBlobs", "Area",
                                    "BoundingBox_xwidth", "BoundingBox_ywidth", "cell_count", "cell_count_resolved", "score"),
           showProgress = FALSE)
d <- d[class_name %in% c("Chaetoceros_spp_chain", "Chaetoceros_spp_single") & is.finite(biovolume_um3) & biovolume_um3 > 0 &
       !is.na(cell_count) & cell_count != -1]
d[, ncell := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1, round(cell_count_resolved), 1L)]
d <- merge(d, meta[qc_pass == TRUE, .(sample, depth, date)], by = "sample")
esd <- function(bv) 2 * (3 / (4 * pi) * bv)^(1/3)
d[, `:=`(pc_esd = esd(biovolume_um3 / ncell), fill = Area / (BoundingBox_xwidth * BoundingBox_ywidth))]

# single-cell population summaries (quoted in the Results)
pop <- rbind(
  d[ncell == 1 & class_name == "Chaetoceros_spp_chain" & score >= cfg$min_score & depth %in% c(3, 16),
    .(definition = "chain class, count 1, score >= 0.9 (previous strip)", n = .N, mean_esd = mean(pc_esd), median_esd = median(pc_esd)), by = depth],
  d[ncell == 1 & score >= cfg$min_score & depth %in% c(3, 16),
    .(definition = "both Chaetoceros classes, count 1, score >= 0.9", n = .N, mean_esd = mean(pc_esd), median_esd = median(pc_esd)), by = depth],
  d[ncell == 1 & score >= cfg$min_score & numBlobs == 1 & fill >= cfg$min_fill & depth %in% c(3, 16),
    .(definition = "selection rule (single object, fill >= 0.25)", n = .N, mean_esd = mean(pc_esd), median_esd = median(pc_esd)), by = depth])
setorder(pop, definition, depth)
fwrite(pop, file.path(cfg$work_dir, "14b_single_cell_population.csv"))
cat("-- single-cell Chaetoceros per-cell ESD, 3 m vs 16 m --\n"); print(pop[, .(definition, depth, n, mean_esd = round(mean_esd, 2), median_esd = round(median_esd, 2))])

draw <- function(dep) {
  x <- d[depth == dep & ncell == 1 & score >= cfg$min_score & numBlobs == 1 & fill >= cfg$min_fill]
  qq <- quantile(x$pc_esd, c(.25, .75))
  x <- x[pc_esd >= qq[1] & pc_esd <= qq[2]]
  x <- x[x[, .I[sample(.N, 1)], by = date]$V1]           # one candidate per sampling day
  x[sample(.N, min(cfg$k_candidates, .N))][, rank := .I][]
}
set.seed(cfg$seed)
sel <- rbind(cbind(layer = "surface_3m", draw(3)), cbind(layer = "deep_16m", draw(16)))
sel <- sel[, .(layer, rank, sample, roi_number, depth, date, class_name, ncell, pc_esd = round(pc_esd, 2),
               fill = round(fill, 2), score = round(score, 3))]
fwrite(sel, file.path(cfg$work_dir, "selection.csv"))
cat(sprintf("\ncandidates: %d surface, %d deep (IQR-bounded, one per day, seed %d)\n",
            sel[layer == "surface_3m", .N], sel[layer == "deep_16m", .N], cfg$seed))
outdir <- file.path(cfg$work_dir, "pngs")
for (s in unique(sel$sample)) {
  rf <- file.path(cfg$raw_root, substr(s, 1, 9), paste0(s, ".roi"))
  ifcb_extract_pngs(rf, out_folder = outdir, ROInumbers = sel[sample == s, roi_number],
                    scale_bar_um = NULL, gamma = 1, normalize = FALSE, overwrite = TRUE, verbose = FALSE)
}
cat(sprintf("selection + PNGs -> %s\nnext: python3 R/15b_single_cell_gallery.py\n", cfg$work_dir))
