# =============================================================================
# Asterionellopsis glacialis validation set, stage 1 of 3: draw 320 A. glacialis ROIs at
# 3 m and 16 m, 40 per depth x maximum-Feret-diameter bin (at most 3 per day), and write
# them as blinded PNGs. The draw is fixed by cfg$seed. Feret diameter tracks the cell count
# of a stellate colony (Spearman 0.91), so the bins are close to count bins.
#
# Stage 2: python3 R/23g_asterionellopsis_ls.py (Label Studio counting project)
# Stage 3: Rscript R/23h_asterionellopsis_analysis.R (runs from the bundled data/asterionellopsis_validation)
#
# Reads:  data/roi_level/roi_2016.txt.gz, outputs/tangesund/02_sample_meta.csv,
#         data/raw_bins/2016
# Writes: outputs/tangesund/asterionellopsis_validation/ images/ (blinded PNGs), 23f_selection_key.csv
#         (blinding key, kept outside images/), 23f_population_bins.csv (stage-3 weights)
# Deps:   data.table, iRfcb
# Run:    Rscript R/23f_asterionellopsis_select.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(iRfcb) })

cfg <- list(
  roi_file    = "data/roi_level/roi_2016.txt.gz",
  meta_csv    = "outputs/tangesund/02_sample_meta.csv",
  raw_root    = "data/raw_bins/2016",
  work_dir    = "outputs/tangesund/asterionellopsis_validation",
  micron      = 1 / 2.77,
  depths      = c(3, 16),
  feret_breaks = c(0, 18, 34, 98, Inf),         # um, ROI maximum Feret diameter
  feret_labels = c("<=18", "18-34", "34-98", ">98"),
  per_cell    = 40,                             # chains per depth x Feret stratum
  max_per_day = 3,                              # cap per sampling day within a stratum
  seed        = 20260913
)
img_dir <- file.path(cfg$work_dir, "images")
tmp_dir <- file.path(cfg$work_dir, "png_raw")
if (length(list.files(img_dir, pattern = "\\.png$")))
  stop(img_dir, " already holds PNGs; the seed pins the blinded ids to the Label Studio tasks, ",
       "so delete the folder deliberately before re-drawing")
for (p in c(img_dir, tmp_dir)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# ---- population -------------------------------------------------------------
meta <- fread(cfg$meta_csv)[qc_pass == TRUE, .(sample, depth, date = as.Date(date))]
d <- fread(cfg$roi_file, select = c("sample", "roi_number", "class_name", "numBlobs", "biovolume_um3",
                                    "cell_count", "cell_count_resolved", "maxFeretDiameter", "score"),
           showProgress = FALSE)
d <- d[grepl("^Asterionellopsis_glacialis", class_name) & !is.na(cell_count) & cell_count != -1 &
       numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0]
d[, det_count := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1,
                         as.integer(round(cell_count_resolved)), 1L)]
d <- merge(d, meta, by = "sample")
d[, feret_um := maxFeretDiameter * cfg$micron]
d[, feret_bin := cut(feret_um, cfg$feret_breaks, labels = cfg$feret_labels, right = TRUE)]
if (anyNA(d$feret_bin)) stop(sum(is.na(d$feret_bin)), " ROIs without a Feret bin")

q <- d[depth %in% cfg$depths, quantile(feret_um, c(0.25, 0.5, 0.75))]
cat(sprintf("pooled 3 m + 16 m Feret quartiles: %s um\n", paste(round(q, 1), collapse = ", ")))
cat(sprintf("Spearman rho, Feret diameter vs detector count: %.2f\n",
            d[depth %in% cfg$depths, cor(feret_um, det_count, method = "spearman")]))

# Feret-bin composition of the whole population, all depths (stage-3 weights)
pop <- d[, .(n = .N, mean_det_count = mean(det_count)), by = .(depth, feret_bin)]
pop[, w := n / sum(n), by = depth]
setorder(pop, depth, feret_bin)
fwrite(pop, file.path(cfg$work_dir, "23f_population_bins.csv"))
cat("-- A. glacialis population by depth x Feret bin --\n")
print(dcast(pop, depth ~ feret_bin, value.var = "w")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

# ---- stratified draw ----------------------------------------------------------
set.seed(cfg$seed)
draw_cell <- function(x, k, cap) {
  x <- x[sample(.N)]                                   # random order
  x[, rank_in_day := seq_len(.N), by = date]
  y <- x[rank_in_day <= cap]
  y <- y[order(rank_in_day, sample(.N))]               # one per day, then a second, then a third
  y[seq_len(min(k, .N))]
}
sel <- d[depth %in% cfg$depths, draw_cell(.SD, cfg$per_cell, cfg$max_per_day), by = .(depth, feret_bin)]
short <- sel[, .N, by = .(depth, feret_bin)][N < cfg$per_cell]
if (nrow(short)) stop("strata below ", cfg$per_cell, " chains:\n", paste(capture.output(print(short)), collapse = "\n"))
sel <- sel[sample(.N)]                                 # shuffle before assigning ids
sel[, anon_id := sprintf("ag_%04d", .I)]
sel <- sel[, .(anon_id, sample, roi_number, depth, date, class_name, feret_um = round(feret_um, 2),
               feret_bin, det_count, score = round(score, 3), rank_in_day)]
cat(sprintf("\ndrawn %d chains: %s\n", nrow(sel),
            paste(capture.output(print(table(sel$depth, sel$feret_bin))), collapse = "\n")))
cat("sampling days and chains beyond the second per day, by stratum:\n")
print(sel[, .(days = uniqueN(date), third_per_day = sum(rank_in_day == 3), detector_cells = sum(det_count)),
          by = .(depth, feret_bin)][order(depth, feret_bin)])
cat(sprintf("detector cells to box in total: %d\n", sel[, sum(det_count)]))

# ---- extract PNGs, then copy to blinded ids ----------------------------------
for (s in unique(sel$sample)) {
  rf <- file.path(cfg$raw_root, substr(s, 1, 9), paste0(s, ".roi"))
  if (!file.exists(rf)) stop("missing raw file: ", rf)
  ifcb_extract_pngs(rf, out_folder = tmp_dir, ROInumbers = sel[sample == s, roi_number],
                    scale_bar_um = NULL, gamma = 1, normalize = FALSE, overwrite = TRUE, verbose = FALSE)
}
raw_pngs <- list.files(tmp_dir, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
key <- data.table(path = raw_pngs,
                  sample = sub("_(\\d+)\\.png$", "", basename(raw_pngs)),
                  roi_number = as.integer(sub(".*_(\\d+)\\.png$", "\\1", basename(raw_pngs))))
sel <- merge(sel, key, by = c("sample", "roi_number"), all.x = TRUE)
if (anyNA(sel$path)) stop("PNG not extracted for ", sum(is.na(sel$path)), " ROIs")
ok <- file.copy(sel$path, file.path(img_dir, paste0(sel$anon_id, ".png")), overwrite = TRUE)
stopifnot(all(ok))
unlink(tmp_dir, recursive = TRUE)
sel[, path := NULL]
setorder(sel, anon_id)
fwrite(sel, file.path(cfg$work_dir, "23f_selection_key.csv"))   # the blinding key: NOT in images/
cat(sprintf("\n%d blinded PNGs -> %s\nkey -> %s\nnext: python3 R/23g_asterionellopsis_ls.py\n",
            nrow(sel), img_dir, file.path(cfg$work_dir, "23f_selection_key.csv")))
