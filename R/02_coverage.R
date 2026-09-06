# =============================================================================
# Per-sample metadata (depth, UTC time), first-pass QC, sampling coverage and the
# station map (Fig. S1). Run first: every later script reads 02_sample_meta.csv.
#
# Reads:  data/roi_level/roi_2016.txt.gz, data/tangesund_samples.csv (bundled; rebuilt
#         from the SMHI 2016 SHARK delivery when that file is present),
#         data/tangesund_coastline_eea.gpkg (fetched from the EEA service if absent).
# Writes: outputs/tangesund/ 02_sample_meta.csv, 02_coverage_by_depth.csv, and the
#         02_*.png figures.
# Deps:   data.table, ggplot2, patchwork, sf, rnaturalearth
# Run:    Rscript R/02_coverage.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")  # English month names on date axes
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(sf); library(rnaturalearth) })
source("R/lib_plotting_theme.R")

cfg <- list(
  roi_file  = "data/roi_level/roi_2016.txt.gz",
  shark_txt = "data/SHARK_PlanktonImaging_2016_Tangesund_SMHI/received_data/shark_data_2016_Tangesund_2026-08-10.txt",
  samples   = "data/tangesund_samples.csv",
  out_dir   = "outputs/tangesund",
  meta_out  = "outputs/tangesund/02_sample_meta.csv",
  junk      = c("Debris", "Beads", "Air_bubbles", "Zooplankton"),  # not-a-cell classes
  depths    = c(3, 6, 8, 11, 13, 16),
  min_cells = 100,   # per-sample classified-cell floor (same as 04)
  min_taxa  = 5      # per-sample taxon floor
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- per-sample metadata: depth + UTC datetime ------------------------------
# One row per sample; depth = MNDEP (== MXDEP); datetime from the sample id (UTC).
# The sample list is the SMHI delivery of the 2016 data to SHARK, made with
# tools/ifcb_delivery_data_to_shark_tangesund.R, which drops blacklisted and flagged
# samples (86 of the extracted samples). That delivery is not public, so its sample
# list is bundled as data/tangesund_samples.csv; with the delivery present the
# bundled file is rebuilt from it, with the analyzed volume (SMVOL) that 02c reads.
# PredictedByMachine rows only: ValidatedByHuman rows would double-count.
if (file.exists(cfg$shark_txt)) {
  sk <- fread(cfg$shark_txt, select = c("SMPNO","STATN","MNDEP","MXDEP","WADEP","SDATE","SMVOL","IFCBNO","IMAGE_VERIFICATION"),
              na.strings = c("","NA"), showProgress = FALSE)
  sk <- sk[IMAGE_VERIFICATION == "PredictedByMachine"]
  samples <- unique(sk[, .(sample = SMPNO, station = STATN, depth = as.numeric(MNDEP),
                           depth_max = as.numeric(MXDEP), water_depth = as.numeric(WADEP),
                           sdate = as.Date(SDATE), ml_analyzed = as.numeric(SMVOL))])
  fwrite(samples, cfg$samples)
  meta <- samples[, !"ml_analyzed"]
} else {
  meta <- fread(cfg$samples, colClasses = list(character = c("sample", "station"),
                                                numeric = c("depth", "depth_max", "water_depth")))
  meta[, `:=`(sdate = as.Date(sdate), ml_analyzed = NULL)]
}
stopifnot(!anyDuplicated(meta$sample))                       # one depth per sample

# datetime from id: D20160810T104734_IFCB110  -> UTC POSIXct
ts <- regmatches(meta$sample, regexpr("D[0-9]{8}T[0-9]{6}", meta$sample))
meta[, datetime := as.POSIXct(ts, format = "D%Y%m%dT%H%M%S", tz = "UTC")]
meta[, `:=`(date  = as.Date(datetime),
            hour  = as.integer(format(datetime, "%H")),
            yday  = as.integer(format(datetime, "%j")))]
meta[, week := as.integer(floor((yday - min(yday)) / 7)) + 1L]           # study-week 1..
meta[, depth_f := factor(depth, levels = cfg$depths)]
setorder(meta, datetime)

cat(sprintf("Tangesund 2016: %d samples | %s to %s | depths %s m | station(s): %s\n",
            nrow(meta), min(meta$date), max(meta$date),
            paste(sort(unique(meta$depth)), collapse = ","), paste(unique(meta$station), collapse = ",")))

# ---- light per-sample ROI/QC summary (stream the one 2016 file) --------------
d <- fread(cfg$roi_file, select = c("sample","class_name","biovolume_um3","numBlobs",
                                    "cell_count","cell_count_resolved"), showProgress = FALSE)
d <- d[numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0]
d[, ncell := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1, round(cell_count_resolved), 1L)]
d[, is_junk := class_name %in% cfg$junk]
qc <- d[, .(n_roi        = .N,
            n_bead       = sum(class_name == "Beads"),
            n_debris     = sum(class_name == "Debris"),
            n_junk       = sum(is_junk),
            n_cells_cls  = sum(ncell[!is_junk]),                 # classified (non-junk) cells
            n_taxa       = uniqueN(class_name[!is_junk]),
            frac_junk    = mean(is_junk)), by = sample]
meta <- merge(meta, qc, by = "sample", all.x = TRUE)
for (col in c("n_roi","n_bead","n_debris","n_junk","n_cells_cls","n_taxa"))
  meta[is.na(get(col)), (col) := 0L]
meta[is.na(frac_junk), frac_junk := 1]

# ---- first-pass QC flags (no samples dropped) --------------------------------
# high_junk: bead/debris-dominated samples (bubbles, blanks)
meta[, `:=`(reliable_spec = n_cells_cls >= cfg$min_cells,
            ok_taxa       = n_taxa      >= cfg$min_taxa,
            high_junk     = frac_junk   >  0.8)]
meta[, qc_pass := reliable_spec & !high_junk]
fwrite(meta, cfg$meta_out)

# ---- coverage tables --------------------------------------------------------
cov_depth <- meta[, .(n = .N, n_days = uniqueN(date), n_reliable = sum(reliable_spec),
                      med_cells = as.integer(round(median(n_cells_cls))),
                      med_taxa  = as.integer(round(median(n_taxa)))), by = depth][order(depth)]
cat("\n-- coverage & reliability per depth --\n"); print(cov_depth)
fwrite(cov_depth, file.path(cfg$out_dir, "02_coverage_by_depth.csv"))

cat(sprintf("\nQC: %d/%d samples pass (>=%d classified cells & <80%% junk); median cells/sample %d, taxa %d\n",
            sum(meta$qc_pass), nrow(meta), cfg$min_cells,
            as.integer(median(meta$n_cells_cls)), as.integer(median(meta$n_taxa))))

# ---- figures ----------------------------------------------------------------
# (1) depth x study-week sampling intensity
wk <- meta[, .(n = .N, date0 = min(date)), by = .(week, depth_f)]
# white labels on the dark (low-count) viridis cells, dark labels on the bright ones
wk[, txt_col := fifelse(sqrt(n) < 0.6 * sqrt(max(n)), "white", "grey10")]
p1 <- ggplot(wk, aes(factor(week), depth_f, fill = n)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n, colour = txt_col), size = 2.6) +
  scale_colour_identity() +
  scale_fill_viridis_c(trans = "sqrt") +
  scale_y_discrete(limits = rev(levels(meta$depth_f))) +
  theme_paper() +
  labs(x = "Study week (from 10 Aug 2016)", y = "Depth (m)", fill = "Samples",
       title = "Depth × time sampling intensity")

# (2) samples through time by depth (the depth-cycling design)
p2 <- ggplot(meta, aes(datetime, depth, colour = depth_f)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_y_reverse(breaks = cfg$depths) + scale_colour_viridis_d(guide = "none") +
  theme_paper() +
  labs(x = NULL, y = "Depth (m)", title = "Sampling through the season (depth-cycling)")

# (3) diel coverage: hour x depth (does every depth get sampled across the diel cycle?)
hr <- meta[, .N, by = .(hour, depth_f)]
p3 <- ggplot(hr, aes(hour, depth_f, fill = N)) +
  geom_tile(colour = "white") +
  scale_fill_viridis_c() + scale_y_discrete(limits = rev(levels(meta$depth_f))) +
  scale_x_continuous(breaks = seq(0, 23, 3)) +
  theme_paper() +
  labs(x = "Hour (UTC)", y = "Depth (m)", fill = "Samples",
       title = "Diel × depth coverage")

# (4) per-sample classified-cell counts (QC context)
p4 <- ggplot(meta, aes(n_cells_cls + 1, fill = reliable_spec)) +
  geom_histogram(bins = 40, colour = "white", linewidth = 0.1) +
  scale_x_log10() + geom_vline(xintercept = cfg$min_cells, linetype = 2) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey70"),
                    breaks = c("TRUE", "FALSE"),
                    labels = c("TRUE" = sprintf("Retained (\u2265 %d cells)", cfg$min_cells),
                               "FALSE" = sprintf("Excluded (< %d cells)", cfg$min_cells)),
                    name = NULL) +
  theme_paper() +
  labs(x = "Classified cells per sample (log)", y = "Samples",
       title = "Per-sample cell counts (QC)")

# (0) station map on the EEA coastline (cached after the first fetch), with a locator inset
sf::sf_use_s2(FALSE)
coast_cache <- "data/tangesund_coastline_eea.gpkg"
if (!file.exists(coast_cache)) {
  eea <- paste0("https://marine.discomap.eea.europa.eu/arcgis/rest/services/Marine/EEA_coastline_2017/MapServer/0/query",
    "?where=1%3D1&geometry=11.0,57.8,11.9,58.35&geometryType=esriGeometryEnvelope&inSR=4326",
    "&spatialRel=esriSpatialRelIntersects&outFields=*&f=geojson")
  st_write(st_read(eea, quiet = TRUE), coast_cache, quiet = TRUE)
}
coast_m <- st_transform(st_make_valid(st_read(coast_cache, quiet = TRUE)), 32632)
pxy <- st_coordinates(st_transform(st_sfc(st_point(c(11.493, 58.075)), crs = 4326), 32632))
win <- st_coordinates(st_transform(st_sfc(st_point(c(11.25, 58.00)), st_point(c(11.75, 58.15)), crs = 4326), 32632))
mx0 <- win[1,1]; mx1 <- win[2,1]; my0 <- win[1,2]; my1 <- win[2,2]
sbx <- mx0 + 0.05*(mx1-mx0); sby <- my0 + 0.08*(my1-my0)
loc <- ne_countries(scale = "medium", returnclass = "sf", continent = "Europe")
p0_inset <- ggplot(loc) + geom_sf(fill = "grey85", colour = "grey60", linewidth = 0.15) +
  # heavy outline keeps the small map extent visible at print size
  annotate("rect", xmin = 11.25, xmax = 11.75, ymin = 58.0, ymax = 58.15, colour = "firebrick", fill = "firebrick", linewidth = 1.6) +
  coord_sf(xlim = c(3, 20), ylim = c(54.5, 62), expand = FALSE) +
  theme_void() + theme(panel.background = element_rect(fill = "#eaf3f7", colour = "grey40", linewidth = 0.5))
p0 <- ggplot() +
  geom_sf(data = coast_m, fill = "grey80", colour = "grey55", linewidth = 0.2) +
  annotation_custom(ggplotGrob(p0_inset), xmin = mx0 + 0.70*(mx1-mx0), xmax = mx1, ymin = my0 + 0.55*(my1-my0), ymax = my1) +
  annotate("point", x = pxy[1], y = pxy[2], shape = 21, size = 3.4, stroke = 1.1, fill = "firebrick", colour = "black") +
  annotate("text", x = pxy[1], y = pxy[2] + 0.07*(my1-my0), label = "IFCB + buoy 33027", size = 3.2, fontface = "bold") +
  annotate("segment", x = sbx, xend = sbx + 5000, y = sby, yend = sby, linewidth = 1.1) +
  annotate("text", x = sbx + 2500, y = sby + 0.045*(my1-my0), label = "5 km", size = 3) +
  coord_sf(xlim = c(mx0, mx1), ylim = c(my0, my1), expand = FALSE, crs = 32632, datum = 4326) +
  theme_paper() +
  theme(panel.background = element_rect(fill = "#eaf3f7"), axis.text = element_text(size = 7),
        panel.grid = element_line(colour = "grey90", linewidth = 0.2)) +
  labs(x = NULL, y = NULL, title = "Study site: Tångesund, Swedish Skagerrak coast")

ggsave(file.path(cfg$out_dir, "02_depth_week_intensity.png"), p1, width = 10, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "02_season_depth_cycling.png"), p2, width = 10, height = 4.0, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "02_diel_depth_coverage.png"), p3, width = 8,  height = 4.0, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "02_sample_cell_counts.png"),  p4, width = 8,  height = 4.0, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "02_coverage_panel.png"),
       ((p0 + labs(tag = "a") | p1 + labs(tag = "b")) /
        (p3 + labs(tag = "c") | p4 + labs(tag = "d"))) &
         theme(plot.tag = element_text(face = "bold", size = 15)),
       width = 16, height = 9, dpi = FIG_DPI)

cat(sprintf("meta cache -> %s | figures -> %s/02_coverage_panel.png (plus the individual 02_*.png panels)\n", cfg$meta_out, cfg$out_dir))
