# =============================================================================
# Build the per-ROI table (biovolume, carbon, chain counts, v4 morphometric features,
# class and score) from the raw IFCB archive. Edit the paths in cfg first.
#
# Reads:  data/raw_bins/<year>/ (.adc/.hdr/.roi), data/features_v4/<year>/ (_fea_v4.csv),
#         data/classified/class<year>*/ (*_class.h5 with a cell_count dataset)
# Writes: data/roi_level/roi_<year>.txt.gz, sample_summary.txt, extraction_metadata.txt
# Deps:   iRfcb (feature-chain-counts branch), data.table, dplyr, hdf5r; arrow optional
# Run:    Rscript R/00b_extract_roi_table.R
#
# The *_class.h5 files must come from a classifier run with chain counting enabled;
# without the cell_count dataset the script aborts. Check one first:
#   Rscript -e 'library(hdf5r); h <- H5File$new("<a>_class.h5"); print(names(h))'
# Only samples with a feature file are extracted: an incomplete features folder
# gives a short table without an error. Check the file count before a long run.
# Sample volume is not included; join data.txt SMVOL by SMPNO.
# =============================================================================

suppressPackageStartupMessages({
  library(iRfcb); library(data.table); library(dplyr); library(hdf5r)
})

# ---- config (edit the paths for the raw-data machine) -----------------------
cfg <- list(
  # add 2023:2025 (and point raw_root at their archive) to rebuild the roi tables 09 reads
  years = 2016,

  # raw .adc/.hdr/.roi per year, in per-day folders
  raw_root      = "data/raw_bins/%d",
  # extracted v4 slim feature files (_fea_v4.csv); created by ifcb_extract_features if absent
  features_root = "data/features_v4/%d",
  # blob files, used only when run_feature_extraction = TRUE
  blobs_root    = "data/blobs/%d",
  # classifier output, one D..._class.h5 per sample; folder found by glob per year
  # (class2016_v3 = ResNet50 V6.1 with 17 chain counters). Errors if the glob is
  # ambiguous; set class_dir to choose.
  classified_root = "data/classified",
  class_glob      = "class%d*",                 # %d -> year
  class_dir       = list(),                      # optional per-year override, e.g.
                                                 #   list(`2016` = ".../class2016_v3")

  out_dir       = "data/roi_level",          # gitignored; per-year per-ROI .txt
  out_format    = "txt",                      # "txt" (tab, gz) or "parquet"

  micron_factor   = 1/2.77,                   # current SMHI factor (older data used 1/3.4)
  feature_version = 4,                         # v4 slim features, not v2
  use_cell_counts = TRUE,
  drop_zero_volume = TRUE,                    # degenerate blobs removed at extraction
  diatom_equation = "auto",                   # small/large Menden-Deuer/Lessard chosen per cell
  carbon_conversion = "cell",                 # per-CELL carbon x resolved count, not whole-ROI
  diatom_include = c("Actinocyclus_spp", "Actinocyclus", "Navicula-like", "Navicula", "cf_Proboscia_rhizosolenia"),
  run_feature_extraction = FALSE,             # TRUE only if the feature files are not built yet (slow)
  parallel = TRUE, n_cores = NULL,

  # NULL keeps every v4 feature column:
  #   Area Biovolume BoundingBox_xwidth BoundingBox_ywidth ConvexArea ConvexPerimeter
  #   Eccentricity EquivDiameter Extent MajorAxisLength MinorAxisLength Orientation
  #   Perimeter RepresentativeWidth Solidity SurfaceArea maxFeretDiameter
  #   minFeretDiameter numBlobs summed{Area,Biovolume,ConvexArea,ConvexPerimeter,
  #   MajorAxisLength,MinorAxisLength,Perimeter,SurfaceArea} Area_over_Perimeter
  #   Area_over_PerimeterSquared summedConvexPerimeter_over_Perimeter
  # Biovolume and SurfaceArea are in pixels; biovolume_um3 is the micron-scaled size.
  morpho_cols = NULL
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- helpers ----------------------------------------------------------------
strip_sample <- function(x) sub("_(fea.*|class|multiblob).*$", "", tools::file_path_sans_ext(basename(x)))

# resolve the class .h5 folder for a year: an explicit override takes precedence; else glob under
# classified_root. Errors on ambiguity so an old dev folder can't be picked silently.
resolve_class_dir <- function(year) {
  ov <- cfg$class_dir[[as.character(year)]]
  if (!is.null(ov)) return(ov)
  hits <- Sys.glob(file.path(cfg$classified_root, sprintf(cfg$class_glob, year)))
  hits <- hits[dir.exists(hits)]
  if (length(hits) == 0) stop("no class folder for ", year, " under ", cfg$classified_root)
  if (length(hits) > 1)
    stop("ambiguous class folders for ", year, " (", paste(basename(hits), collapse = ", "),
         ") -> set cfg$class_dir[['", year, "']] to the production one")
  hits
}

# per-ROI class_name (thresholded, incl. "unclassified"), class_name_auto (argmax),
# score and model_name from a _class.h5
read_h5_scores <- function(h5_path) {
  h <- H5File$new(h5_path, mode = "r")
  on.exit(h$close_all())
  # first existing dataset among candidate names (older h5 files use other names)
  rd <- function(names) {
    for (nm in names) { v <- tryCatch(h[[nm]]$read(), error = function(e) NULL)
                        if (!is.null(v)) return(v) }
    NULL
  }
  roi <- rd("roi_numbers")
  if (is.null(roi) || length(roi) == 0) return(NULL)          # empty sample -> skip
  n    <- length(roi)
  thr  <- rd(c("class_name", "class_labels_above_threshold")) # canonical: threshold-adjusted, and what built data.txt
  auto <- rd(c("class_name_auto", "class_labels_auto"))       # argmax, no threshold
  scr  <- rd("output_scores")
  model <- rd(c("model_name", "classifier_name", "classifierName"))
  # winning probability per ROI; handles 1-ROI samples (vector) and either matrix orientation
  if (is.null(dim(scr))) {
    score <- max(scr)
  } else {
    roi_dim <- if (ncol(scr) == n) 2L else 1L
    score <- apply(scr, roi_dim, max)
  }
  data.table(
    sample          = strip_sample(h5_path),
    roi_number      = as.integer(roi),
    class_name      = thr,
    class_name_auto = auto,
    score           = score,
    model_name      = if (is.null(model)) NA_character_ else model[1]
  )
}

# ---- per-year extraction ----------------------------------------------------
process_year <- function(year) {
  raw    <- sprintf(cfg$raw_root, year)
  feats  <- sprintf(cfg$features_root, year)
  blobs  <- sprintf(cfg$blobs_root, year)
  classd <- resolve_class_dir(year)
  message("\n==== year ", year, " ====")
  message("  class folder: ", classd, "  (check that this is the production classifier)")

  # 0) feature extraction (slow, only if needed)
  if (cfg$run_feature_extraction) {
    message("  extracting features (slow) ...")
    # feature_tag = "fea" matches 00a; mixed name styles in one folder are read twice
    iRfcb::ifcb_extract_features(data_folder = raw, features_folder = feats,
                                 blobs_folder = blobs, feature_tag = "fea",
                                 parallel = cfg$parallel,
                                 n_cores = cfg$n_cores, overwrite = FALSE, verbose = TRUE)
  }
  stopifnot(dir.exists(feats))

  # 1) biovolume + carbon + chain counts (one row per ROI)
  message("  ifcb_extract_biovolumes() ...")
  bv <- iRfcb::ifcb_extract_biovolumes(
    feature_files = feats, class_files = classd,
    micron_factor = cfg$micron_factor, use_cell_counts = cfg$use_cell_counts,
    drop_zero_volume = cfg$drop_zero_volume,
    diatom_equation = cfg$diatom_equation,
    carbon_conversion = cfg$carbon_conversion,
    diatom_include = cfg$diatom_include) %>%
    as.data.table()
  bv[, roi_number := as.integer(roi_number)]

  # 2) curated morphometric features
  message("  ifcb_read_features() + select morpho subset ...")
  fl <- iRfcb::ifcb_read_features(feats, feature_version = cfg$feature_version)
  feat <- rbindlist(lapply(names(fl), function(nm) {
    d <- as.data.table(fl[[nm]]); d[, sample := strip_sample(nm)]
    keep <- if (is.null(cfg$morpho_cols)) setdiff(names(d), c("sample", "roi_number"))
            else intersect(cfg$morpho_cols, names(d))
    d[, c("sample", "roi_number", keep), with = FALSE]
  }), use.names = TRUE, fill = TRUE)
  feat[, roi_number := as.integer(roi_number)]

  # 3) class confidence score from production .h5
  message("  reading *_class.h5 scores ...")
  h5s <- list.files(classd, pattern = "_class\\.h5$", full.names = TRUE)
  scores <- rbindlist(lapply(h5s, read_h5_scores), use.names = TRUE, fill = TRUE)

  # 4) left-join onto the biovolume table; every ROI kept. class_name (from the h5)
  # is the canonical class and matches data.txt.
  roi <- bv %>%
    merge(feat,   by = c("sample", "roi_number"), all.x = TRUE) %>%
    merge(scores[, .(sample, roi_number, class_name, class_name_auto, score, model_name)],
          by = c("sample", "roi_number"), all.x = TRUE)
  setorder(roi, sample, roi_number)

  # 5) write per-ROI table
  out <- file.path(cfg$out_dir, sprintf("roi_%d.%s", year,
                   ifelse(cfg$out_format == "parquet", "parquet", "txt.gz")))
  if (cfg$out_format == "parquet" && requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(roi, out)
  } else {
    fwrite(roi, out, sep = "\t")   # .gz inferred from extension
  }
  message(sprintf("  wrote %s  (%s ROIs, %s samples)", out,
                  format(nrow(roi), big.mark = ","), uniqueN(roi$sample)))

  # 6) per-sample summary for QC (n_roi, n_cells, total biovolume; volume joined later)
  smry <- roi[, .(n_roi = .N,
                  n_cells = sum(cell_count_resolved %||% cell_count, na.rm = TRUE),
                  total_biovol_um3 = sum(biovolume_um3, na.rm = TRUE),
                  median_score = median(score, na.rm = TRUE)), by = sample]
  list(roi_meta = data.table(year = year, n_roi = nrow(roi), n_samples = uniqueN(roi$sample),
                             class_dir = basename(classd),
                             model_name = scores$model_name[1] %||% NA), summary = smry)
}

`%||%` <- function(a, b) if (is.null(a) || all(is.na(a))) b else a

# ---- run --------------------------------------------------------------------
results <- lapply(cfg$years, function(y) tryCatch(process_year(y),
                  error = function(e) { message("  YEAR ", y, " FAILED: ", conditionMessage(e)); NULL }))
results <- Filter(Negate(is.null), results)

# stop if every year failed, so the existing summary and sidecar are not overwritten
if (length(results) == 0)
  stop("no year produced a per-ROI table, so sample_summary.txt and ",
       "extraction_metadata.txt were left as they were. See the messages above ",
       "for the per-year failure.")

# per-sample summary across years
smry_all <- rbindlist(lapply(results, `[[`, "summary"), use.names = TRUE, fill = TRUE)
fwrite(smry_all, file.path(cfg$out_dir, "sample_summary.txt"), sep = "\t")

# ---- reproducibility sidecar ------------------------------------------------
desc <- tryCatch(packageDescription("iRfcb"), error = function(e) NULL)
meta <- c(
  sprintf("extracted_on: %s", Sys.Date()),
  sprintf("iRfcb_version: %s", as.character(utils::packageVersion("iRfcb"))),
  sprintf("iRfcb_remote_ref: %s", desc$RemoteRef %||% "NA"),
  sprintf("iRfcb_remote_sha: %s", desc$RemoteSha %||% "NA"),
  sprintf("micron_factor: %s (= %.5f)", "1/2.77", cfg$micron_factor),
  sprintf("feature_version: %s", cfg$feature_version),
  sprintf("use_cell_counts: %s", cfg$use_cell_counts),
  sprintf("biovolume_units: um3 (per-ROI/per-chain; divide by cell_count for per-cell)"),
  sprintf("morpho_cols: %s", paste(cfg$morpho_cols, collapse = ", ")),
  sprintf("years_done: %s", paste(sapply(results, function(r) r$roi_meta$year), collapse = ", ")),
  sprintf("classifiers (VERIFY = production, not dev): %s",
          paste(sapply(results, function(r) sprintf("%s -> %s [%s]",
                r$roi_meta$year, r$roi_meta$model_name, r$roi_meta$class_dir)), collapse = "; ")),
  sprintf("note: sample volume (ml analyzed) NOT here -> join data.txt SMVOL by SMPNO")
)
writeLines(meta, file.path(cfg$out_dir, "extraction_metadata.txt"))

cat("\n================ EXTRACTION DONE ================\n")
print(rbindlist(lapply(results, `[[`, "roi_meta"), fill = TRUE))
cat("\nper-ROI tables + sample_summary.txt + extraction_metadata.txt in", cfg$out_dir, "\n")
cat("join key to delivered data: sample == SMPNO ; volume from data.txt SMVOL\n")
cat("================================================\n")
