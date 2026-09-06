# =============================================================================
# Download the Tangesund raw IFCB bins from the HABON-IFCB Dashboard (dataset
# "tangosund", 2,829 bins) and extract v4 feature files from them.
#
# Reads:  nothing local; needs network access, and Python for the extraction stage.
# Writes: data/raw_bins/2016/<Dyyyymmdd>/<sample>.{roi,adc,hdr}   (~18 GB)
#         data/features_v4/2016/<sample>_fea_v4.csv               (~0.7 GB)
#         data/blobs/2016/                                        (~1.4 GB, not read downstream)
#         data/dashboard_tangosund_bins.csv                       (per-bin manifest)
# Env:    STAGE=download|features|all (default all); SAMPLES=<file of sample ids, one per line>
# Deps:   iRfcb (>= 0.9.0), reticulate, and WHOI ifcb-features pinned to v1.0.0:
#           Rscript -e 'iRfcb::ifcb_py_install(features = TRUE, features_ref = "v1.0.0")'
#         or use an existing Python: RETICULATE_PYTHON=/path/to/python
# Run:    Rscript R/00a_download_dashboard_data.R
#
# Resumable: files already on disk are skipped. Expect hours, mostly extraction.
# The dashboard serves no classifier output. 00b_extract_roi_table.R also needs
# the *_class.h5 files with the per-ROI cell_count dataset (see the README).
# Dashboard v4 feature files are not used: they come from a different
# ifcb-features build and differ slightly from the paper's.
# =============================================================================
suppressPackageStartupMessages(library(iRfcb))

cfg <- list(
  dashboard_url = "https://habon-ifcb.whoi.edu",
  dataset       = "tangosund",          # dashboard spelling
  year          = 2016,

  raw_dir       = "data/raw_bins/2016",   # iRfcb creates the per-day Dyyyymmdd/ folders
  features_dir  = "data/features_v4/2016",
  blobs_dir     = "data/blobs/2016",      # not read downstream
  manifest_csv  = "data/dashboard_tangosund_bins.csv",

  # which stages to run; STAGE=download|features|all overrides
  stage         = Sys.getenv("STAGE", "all"),
  # optional file of sample ids, one per line; empty means every bin in the dataset
  samples_file  = Sys.getenv("SAMPLES", ""),

  chunk_size    = 100L,   # samples per download step
  parallel_downloads = 5L,
  sleep_time    = 2,
  max_retries   = 3L,
  multi_timeout = 300,

  feature_parallel = TRUE,
  n_cores       = NULL    # NULL -> all but one
)

# Expected v4 feature columns, compared with the first extracted file
v4_cols <- c("roi_number", "Area", "Biovolume", "BoundingBox_xwidth", "BoundingBox_ywidth",
             "ConvexArea", "ConvexPerimeter", "Eccentricity", "EquivDiameter", "Extent",
             "MajorAxisLength", "MinorAxisLength", "Orientation", "Perimeter",
             "RepresentativeWidth", "Solidity", "SurfaceArea", "maxFeretDiameter",
             "minFeretDiameter", "numBlobs", "summedArea", "summedBiovolume",
             "summedConvexArea", "summedConvexPerimeter", "summedMajorAxisLength",
             "summedMinorAxisLength", "summedPerimeter", "summedSurfaceArea",
             "Area_over_PerimeterSquared", "Area_over_Perimeter",
             "summedConvexPerimeter_over_Perimeter")

# Pinned ifcb-features release. Accepted commits: v1.0.0 and the docs-only commit after it.
features_ref <- "v1.0.0"
features_ok_commits <- c("7da48de", "68550ea")
install_call <- sprintf("iRfcb::ifcb_py_install(features = TRUE, features_ref = \"%s\")", features_ref)

stopifnot(cfg$stage %in% c("download", "features", "all"))
do_download <- cfg$stage %in% c("download", "all")
do_features <- cfg$stage %in% c("features", "all")
dataset_url <- paste0(sub("/+$", "", cfg$dashboard_url), "/", cfg$dataset, "/")

# ---- helpers ----------------------------------------------------------------

in_chunks <- function(x, n) split(x, ceiling(seq_along(x) / n))

# a raw bin counts as present only with all three files (extraction needs the .adc)
have_raw <- function() {
  ext <- function(e) sub(paste0("\\.", e, "$"), "", basename(Sys.glob(file.path(cfg$raw_dir, "*", paste0("*.", e)))))
  Reduce(intersect, list(ext("roi"), ext("adc"), ext("hdr")))
}

have_features <- function()
  sub("_fea_v4\\.csv$", "", list.files(cfg$features_dir, pattern = "_fea_v4\\.csv$"))

# Fails early if Python or ifcb-features is missing
check_python <- function() {
  if (!reticulate::py_available(initialize = TRUE))
    stop("the feature extraction needs Python, and reticulate found none.\n",
         "  Rscript -e '", install_call, "'\n",
         "  or run with STAGE=download to fetch the raw bins only.")
  missing <- Filter(function(m) inherits(try(reticulate::import(m, delay_load = FALSE), silent = TRUE), "try-error"),
                    c("ifcb", "ifcb_features"))
  if (length(missing))
    stop("the WHOI ifcb-features Python package is not importable from ",
         reticulate::py_exe(), "\n  missing module(s): ", paste(missing, collapse = ", "), "\n",
         "  Rscript -e '", install_call, "'\n",
         "  or point R at a Python that already has it:\n",
         "  RETICULATE_PYTHON=/path/to/python Rscript R/00a_download_dashboard_data.R")
  message("  python: ", reticulate::py_exe())
  report_features_version()
}

# Reports the ifcb-features build from the git commit pip recorded (the package
# version string does not track release tags). A mismatch warns but does not stop.
report_features_version <- function() {
  info <- tryCatch(reticulate::py_run_string(paste(
    "import importlib.metadata as _md, json as _json",
    "_d = _md.distribution('ifcb-features')",
    "_v = _d.version",
    "_t = _d.read_text('direct_url.json')",
    "_c = _json.loads(_t)['vcs_info']['commit_id'] if _t else None",
    sep = "\n"), local = TRUE, convert = TRUE), error = function(e) NULL)
  if (is.null(info)) { message("  ifcb-features: version could not be determined"); return(invisible(NULL)) }
  commit <- info$`_c`
  message(sprintf("  ifcb-features: %s (commit %s)", info$`_v`,
                  if (is.null(commit)) "unknown" else substr(commit, 1, 12)))
  if (!is.null(commit) && !substr(commit, 1, 7) %in% features_ok_commits)
    message("  NOTE: this is not the build the paper used (", features_ref, ").\n",
            "        Feature values may differ slightly. Reinstall with:\n",
            "        Rscript -e '", install_call, "'")
  invisible(NULL)
}

check_v4 <- function() {
  f <- list.files(cfg$features_dir, pattern = "_fea_v4\\.csv$", full.names = TRUE)
  if (length(f) == 0) return(invisible(NULL))
  got <- names(utils::read.csv(f[1], nrows = 1))
  if (!identical(got, v4_cols))
    stop("the extracted feature files are not the expected v4 set.\n",
         "  first file: ", basename(f[1]), "\n",
         "  missing:    ", paste(setdiff(v4_cols, got), collapse = ", "), "\n",
         "  unexpected: ", paste(setdiff(got, v4_cols), collapse = ", "), "\n",
         "00b_extract_roi_table.R assumes v4; do not run it on these.")
  invisible(NULL)
}

download_bins <- function(samples) {
  dir.create(cfg$raw_dir, recursive = TRUE, showWarnings = FALSE)
  chunks <- in_chunks(samples, cfg$chunk_size)
  for (i in seq_along(chunks)) {
    message(sprintf("  [download] chunk %d/%d (%d samples)", i, length(chunks), length(chunks[[i]])))
    tryCatch(
      ifcb_download_dashboard_data(
        dashboard_url = dataset_url, samples = chunks[[i]], file_types = c("roi", "adc", "hdr"),
        dest_dir = cfg$raw_dir, parallel_downloads = cfg$parallel_downloads,
        sleep_time = cfg$sleep_time, multi_timeout = cfg$multi_timeout,
        max_retries = cfg$max_retries, quiet = TRUE),
      error = function(e)
        # keep going; the summary reports shortfalls and a rerun resumes
        message("    chunk failed, continuing: ", conditionMessage(e)))
  }
}

# ---- 1. per-bin manifest ----------------------------------------------------
# Dashboard metadata (sample_time, depth, ml_analyzed, position, n_images per bin);
# ifcb_list_dashboard_bins() no longer works with the dashboard.
message("fetching bin metadata for dataset '", cfg$dataset, "' ...")
meta <- as.data.frame(ifcb_download_dashboard_metadata(cfg$dashboard_url,
                                                       dataset_name = cfg$dataset, quiet = TRUE))
all_samples <- meta$pid
dir.create(dirname(cfg$manifest_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(meta, cfg$manifest_csv, row.names = FALSE)
message(sprintf("  %d bins, %s to %s, %s ROIs -> %s",
                length(all_samples), min(meta$sample_time), max(meta$sample_time),
                format(sum(meta$n_images, na.rm = TRUE), big.mark = ","), cfg$manifest_csv))

if (nzchar(cfg$samples_file)) {
  want <- trimws(readLines(cfg$samples_file, warn = FALSE))
  want <- want[nzchar(want)]
  unknown <- setdiff(want, all_samples)
  if (length(unknown))
    stop(length(unknown), " sample id(s) in ", cfg$samples_file,
         " are not in this dataset, e.g. ", paste(utils::head(unknown, 3), collapse = ", "))
  message(sprintf("  restricted to %d samples from %s", length(want), cfg$samples_file))
} else {
  want <- all_samples
}

if (do_features) { message("checking the Python feature-extraction environment ..."); check_python() }

# ---- 2. raw bins ------------------------------------------------------------
if (do_download) {
  todo <- setdiff(want, have_raw())
  message(sprintf("\nraw bins: %d of %d already present, downloading %d (~%.1f GB) -> %s",
                  length(want) - length(todo), length(want), length(todo),
                  length(todo) * 0.0064, cfg$raw_dir))
  if (length(todo)) download_bins(todo)
}

# ---- 3. feature extraction --------------------------------------------------
# feature_tag = "fea" writes <sample>_fea_v4.csv; do not mix with the default
# "features" names in one folder (ifcb_extract_biovolumes() would read both).
# overwrite = FALSE skips bins that already have a feature file.
if (do_features) {
  dir.create(cfg$features_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$blobs_dir, recursive = TRUE, showWarnings = FALSE)
  present <- have_raw()
  todo <- setdiff(intersect(want, present), have_features())
  message(sprintf("\nfeatures: %d of %d raw bins already extracted, extracting %d -> %s",
                  length(intersect(want, present)) - length(todo), length(intersect(want, present)),
                  length(todo), cfg$features_dir))
  if (length(setdiff(want, present)))
    message(sprintf("  (%d requested sample(s) have no raw bin on disk and are skipped)",
                    length(setdiff(want, present))))
  if (length(todo)) {
    ifcb_extract_features(
      data_folder = cfg$raw_dir, features_folder = cfg$features_dir,
      blobs_folder = cfg$blobs_dir, bins = todo, feature_tag = "fea",
      parallel = cfg$feature_parallel, n_cores = cfg$n_cores,
      overwrite = FALSE, verbose = TRUE)
  }
  check_v4()
}

# ---- 4. tally ---------------------------------------------------------------
n_raw  <- length(have_raw())
n_feat <- length(have_features())
short  <- function(n, requested, ran) {
  if (!ran) "(stage not run)" else if (n >= requested) "complete" else "SHORT -> re-run to resume"
}

cat("\n================ SUMMARY ================\n")
cat(sprintf("dataset               : %s (%s)\n", cfg$dataset, dataset_url))
cat(sprintf("bins in dataset       : %d\n", length(all_samples)))
cat(sprintf("requested             : %d\n", length(want)))
cat(sprintf("raw bins (roi+adc+hdr): %d / %d  %s\n", n_raw, length(want), short(n_raw, length(want), do_download)))
cat(sprintf("feature files         : %d / %d  %s\n", n_feat, length(want), short(n_feat, length(want), do_features)))

cat("\n---- NOT AVAILABLE FROM THE DASHBOARD ----\n")
cat("The 'tangosund' dataset serves no classifier output: the autoclass file type\n")
cat("returns 404. That means this download alone cannot rebuild the per-ROI table.\n")
cat("00b_extract_roi_table.R additionally needs, and will abort without:\n")
cat("  data/classified/class2016_v3/<sample>_class.h5   one per sample, from a\n")
cat("  classifier run with chain counting enabled, carrying the per-ROI `cell_count`\n")
cat("  dataset that every chain-length result in the paper is built on.\n")
cat("Obtain those separately (see 'Data required' in the README). To skip stage 0\n")
cat("entirely, get the finished data/roi_level/roi_2016.txt.gz the same way.\n")

cat("\nnext: Rscript R/00b_extract_roi_table.R   (once the *_class.h5 files are in place)\n")
cat("=========================================\n")
