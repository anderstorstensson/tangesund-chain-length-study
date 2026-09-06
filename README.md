# Tångesund diatom chain-length study

Analysis code for:

> Torstensson, A., M. L. Brosnahan, A. Cembella, and B. Karlson. *Decline in diatom
> chain length and cell size with depth across a coastal halocline.* (manuscript
> submitted)

An Imaging FlowCytobot (IFCB) sampled six depths (3–16 m) at Tångesund on the Swedish
Skagerrak coast from 10 August to 17 October 2016 (2,743 samples). YOLO detectors count
the cells in each imaged diatom chain. The scripts take the per-ROI tables to every figure
and statistic in the paper. Supplementary items are cited by S-number (Fig. S5 is
Supplementary Fig. 5, Data S1 is Supplementary Data 1).

## Repository layout

The repository holds code and five small bundled datasets (`data/chaetoceros_validation/`,
`data/asterionellopsis_validation/`, `data/svea/`, `data/tangesund_samples.csv` and
`data/functional_group_map.csv`). All other inputs are obtained separately (see [Data required](#data-required))
and placed under `data/`. Outputs go to `outputs/tangesund/`. Apart from the bundled files,
both directories are git-ignored. `tools/` holds a record of the SMHI delivery script that
made the sample list; it is not part of the pipeline.

To check the R setup without any data:

```bash
Rscript tests/test_size_spectra_synthetic.R      # prints ALL TESTS PASSED
```

## Running the pipeline

Run every script from the repository root, for example `Rscript R/03_chain_length.R`. Each
script starts with a header listing what it reads, writes and needs. Output file names start
with the number of the script that wrote them (`06_env_matched.csv` comes from `06_env_buoy.R`).

Run the scripts in numerical order, with lettered scripts (`02b`, `09b`, ...) directly after
their parent. The exceptions:

- Run `23e` before `20`, which reads `23e_class_purity.csv` and drops that variant when the
  file is absent.
- `09b` reads outputs of `25` and `27c`. Run it and `09c` after `27c`, then run `10` again.
- After `09b`, run `27c` again with `BOOT_PARTS=B` to write `27c_svea_share_of_reference.csv`.
- `27d` reads outputs from across the pipeline, including `27c_pertaxon_model_boot.csv`
  (run `27c` with `BOOT_PARTS=E`). Run it after the `BOOT_PARTS=B` rerun.
- `04` and `16` read `04b` and `16b` when present. Run the pair, then the parent again, or
  let `28` re-run it.
- Run `28` last.

The negative-binomial refits in `20` take about two hours.

```
R/
  lib_size_spectra.R                bounded power-law MLE for size spectra
  lib_plotting_theme.R              ggplot theme, season palette, raster resolution

  # Stage 0: inputs (needs the raw archive; most users skip it)
  00a_download_dashboard_data.R     raw bins and v4 features from the HABON-IFCB Dashboard
  (ifcb-classify infer)             class files with cell counts (see Classifying the bins)
  00b_extract_roi_table.R           per-ROI tables from raw bins, features and class files
  01_build_functional_groups.R      data/functional_group_map.csv (WoRMS via SHARK4R; bundled)

  # Stage 1: analysis
  02_coverage.R                     sample metadata cache, QC, station map, coverage
  02b_roi_trigger_counts.R          triggers and ROIs per sample (.adc)
  02c_ifcb_settings.R               acquisition settings and analyzed volume (.hdr)
  02d_trigger_fluorescence.R        fluorescence per cell and near-trigger ROIs by depth
  03_chain_length.R                 chain length by depth and season
  04_size_spectra.R                 size-spectrum slope, size diversity, median ESD
  04b_size_profile_raw.R            per-depth means of the size metrics
  05_composition.R                  functional-group shares, rarefied richness and Shannon
  06_env_buoy.R                     buoy hydrography matched to each sample (Fig. S1e)
  08_diel_chain.R                   diel chain-length models
  09_regional_context.R             regional envelope from the Svea transect
  09b_regional_salinity.R           chain length along the Svea salinity gradient (Figs. S5, 6d)
  09c_svea_autumn.R                 the 09b contrast on autumn cruises
  10_env_mechanism.R                chain length vs hydrography GAMMs (Fig. 6)
  10b_october_halocline.R           October halocline shoaling (Fig. S4)
  10c_temperature_variance.R        variance decomposition of matched temperature
  10d_salinity_bin_months.R         month composition of the Fig. 6 salinity bins
  11_nutrients.R                    bottle nutrients and size-fractionated chlorophyll
  13_contrasts.R                    surface-to-deep contrasts and sensitivities
  14_chaeto_gallery_select.R        select Chaetoceros ROIs for the gallery
  14b_single_cell_select.R          select single-cell Chaetoceros exemplars
  15_chaeto_gallery.py              detector counts and gallery panel (Fig. 4d)
  15b_single_cell_gallery.py        single-cell plate (Fig. S3)
  16_percell_size.R                 per-cell size by depth (Fig. 4)
  16b_within_depth_correlation.R    chain length x cell size correlation within depth
  16c_chaeto_size_density.R         Chaetoceros per-cell size densities (Fig. S6)
  17_plasticity_turnover.R          within- vs between-taxon partition
  17b_matched_day_partition.R       the partition on matched days
  18_chain_distribution.R           chain-length quantiles by depth

  # Stage 2: chain-counter validation
  19_chain_counter_validation.R     detector accuracy
  20_within_taxon_chain.R           within-taxon models (Table 2)
  21_chaeto_size_validation.py      size-stratified Chaetoceros count validation
  22_chaeto_size_reconcile.R        Chaetoceros size-bias reconciliation
  23_chaeto_size_matched.R          size coverage and bias-corrected contrast
  23b_detection_by_depth.R          unclassified and zero-detection rates by depth
  23c_chaetoceros_select.R          draw the blinded Chaetoceros validation set
  23d_chaetoceros_ls.py             Label Studio counting project
  23e_chaetoceros_analysis.R        depth-differential count bias (Tables S4, S5)
  23f_asterionellopsis_select.R     draw the blinded Asterionellopsis glacialis validation set
  23g_asterionellopsis_ls.py        Label Studio counting project
  23h_asterionellopsis_analysis.R   depth-differential count bias (Tables S4, S5)
  23i_validation_figure.R           both blind counts in one figure (Fig. S7, Table S5)

  # Stage 3: sensitivity and model checks
  24_trait_threshold_sensitivity.R  chain-trait sample threshold
  25_weight_sensitivity.R           precision-weighted refits; log-scale reference
  26_ar1_sensitivity.R              AR1 residual autocorrelation
  27_model_diagnostics.R            depth x hour coverage, concurvity
  27b_model_terms.R                 term statistics for every cited model (Data S1)
  27c_day_bootstrap.R               day-resampled intervals (Fig. S2)
  27d_si_robustness.R               sensitivity analyses in one figure (Fig. S8, Tables S2, S3, S7, S8, Data S1)
  27e_size_limits_light.R           intake-screen and frame limits; daily surface PAR (STRÅNG)

  28_render_submission_figures.R    all figures as vector PDF and 600-dpi PNG
tests/
  test_size_spectra_synthetic.R     unit test for lib_size_spectra.R
tools/
  ifcb_delivery_data_to_shark_tangesund.R  SMHI delivery to SHARK (record only; see Data required)
```

`28_render_submission_figures.R` re-runs each figure script and writes to
`outputs/tangesund/figures/`. `FIG_ONLY="Fig4,FigS5"` renders a subset.

`23c`, `23d` and `23e` run in that order and need a Label Studio server, as do `23f`, `23g` and
`23h` for *Asterionellopsis glacialis*. The manual counts and blinded images are bundled in
`data/chaetoceros_validation/` and `data/asterionellopsis_validation/`, and `23e` and `23h` run from these without a server.

## Data required

Paths are relative to the repository root. Two of them contain `å`.

### Analysis inputs

| Path | Read by | Source |
|---|---|---|
| `data/roi_level/roi_2016.txt.gz` | most scripts | built by `00b` from the `.h5` files |
| `data/roi_level/roi_2022.txt.gz` to `roi_2026.txt.gz` | `09` (2023–2025), `09b` | not public; `09b` and `09c` fall back to `data/svea/` |
| `data/tangesund_samples.csv` | `02`, `02c` | bundled; the delivered sample list with depths and analyzed volumes (see below) |
| `data/SHARK_PlanktonImaging_<year>_SMHI/processed_data/data.txt` and `.../received_data/class_names.txt` | `01`, `09`, `09b`, `09c` | SHARK (2022–2024 published); `01` needs every year, including the unpublished 2016 delivery |
| `data/functional_group_map.csv` | `03`, `05` and others | bundled; built by `01` |
| `data/tångesund_bouy/smhi-opendata_*_33027_*.csv` | `06`, `10` | SMHI open data, buoy 33027 |
| `data/tångesund_water_samples/skagerrak_chemical_observations.csv` (or the `Results …xlsx` original) | `11` | Zenodo, doi:10.5281/zenodo.18378122 |
| `data/chains/chains_<detector>_yolo11s/results.csv` and `data/chains_counter_accuracy.csv` | `19` (and `20` through it) | not public |
| `data/chains/chains_chaetoceros_yolo11s/weights/best.pt` | `15`, `15b`, `21`, `23d` | SciLifeLab Data Repository, doi:10.17044/scilifelab.33565219 |
| `data/chains/chains_asterionellopsis_yolo11s/weights/best.pt` | `23g` (with `--preannotate`) | SciLifeLab Data Repository, doi:10.17044/scilifelab.33565219 |
| `data/chains_validation/chaetoceros_test/` | `21` | [IFCBChainCounts](https://github.com/SAMS-UK/IFCBChainCounts) |
| `data/tangesund_coastline_eea.gpkg` | `02` | downloaded and cached by `02` |
| `data/strang/strang_par_daily_2016.json` | `27e` | SMHI STRÅNG (https://strang.smhi.se); downloaded and cached by `27e` |
| `data/svea/` | `09b`, `09c` | bundled |
| `data/chaetoceros_validation/` | `23e` | bundled |
| `data/asterionellopsis_validation/` | `23h` | bundled |

`data/tangesund_samples.csv` lists the 2,743 samples of the SMHI delivery of the 2016 data
to SHARK, with their depths and analyzed volumes. That delivery is not public. It was made with
`tools/ifcb_delivery_data_to_shark_tangesund.R` (a record of the SMHI pipeline, not runnable
outside SMHI), which leaves out blacklisted samples and those flagged for bubbles, incomplete
runs or beads. With the delivery file present, `02` rebuilds the sample list from it.
`data/functional_group_map.csv` is bundled because `01` needs that same delivery and queries
WoRMS, whose answers can change.

`data/svea/` contains the *Svea* samples, the sample × taxon chain table and the cruise-number
corrections. `09b` rewrites it whenever it runs from the SHARK deliveries and ROI tables.
The ferrybox salinity is also public as platform `SveaFB` in Copernicus Marine product
INSITU_BAL_PHYBGCWAV_DISCRETE_MYNRT_013_032 (doi:10.48670/moi-00032).

### Raw inputs

Read by stage 0, `02b` to `02d`, `14`, `14b`, `23c` and `23f`.

| Path | Size | Notes |
|---|---|---|
| `data/raw_bins/2016/<Dyyyymmdd>/` | 18 GB | `.adc`, `.hdr` and `.roi` bins, fetched by `00a` |
| `data/features_v4/2016/` | 0.7 GB | 2,829 `_fea_v4.csv` files, made by `00a` |
| `data/classified/class2016_v3/*_class.h5` | 2.7 GB | classifier output with a `cell_count` dataset, made by `ifcb-classify` (see below) |

`00b` extracts only the samples it finds a feature file for. Check the count first
(`ls data/features_v4/2016 | wc -l` should give 2,829). It stops if the `.h5` files have
no `cell_count` dataset:

```bash
Rscript -e 'library(hdf5r); h <- H5File$new("<sample>_class.h5"); print(names(h))'
```

The dashboard holds 2,829 bins, and `02` keeps the 2,743 that are in the SHARK delivery
(`data/tangesund_samples.csv`).

### Downloading from the dashboard

```bash
Rscript -e 'iRfcb::ifcb_py_install(features = TRUE, features_ref = "v1.0.0")'
Rscript R/00a_download_dashboard_data.R      # 18 GB, several hours, resumable
```

`00a` downloads dataset `tangosund` from <https://habon-ifcb.whoi.edu/timeline?dataset=tangosund>,
writes `data/dashboard_tangosund_bins.csv` and extracts the features locally. The dashboard's
own v4 feature files come from a different `ifcb-features` build and are not used.
`STAGE=download` or `STAGE=features` runs one stage, and `SAMPLES=ids.txt` limits the run
to the listed samples. `00b` also needs the `.h5` files, which the dashboard does not serve.

### Classifying the bins

The class files come from [`ifcb-classify`](https://doi.org/10.5281/zenodo.21875810) 0.4.0
with the ResNet-50 V6.1 classifier, its class thresholds and the per-taxon YOLO11s cell
counters, all in the SciLifeLab Data Repository (doi:10.17044/scilifelab.33565219) together
with the inference configuration `inference_config_v6_1.yaml`. Unpack the deposit, then

```bash
pip install "ifcb-classify[chains] @ git+https://github.com/nodc-sweden/ifcb-pytorch-classify@v0.4.0"
ifcb-classify infer --config inference_config_v6_1.yaml \
  --input data/raw_bins/2016 --output data/classified/class2016_v3
```

The configuration enables chain counting; without it the `.h5` files have no `cell_count`
dataset and `00b` stops. The paper's class files were made with 17 detectors. The deposit
also holds a detector for *Nitzschia paleacea*, trained later; if the configuration maps
`Nitzschia_paleacea`, remove that entry to reproduce the paper's files. The paper's files were made on a CUDA GPU. A re-run of five bins
on a CPU with v0.4.0 gave the same class labels and cell counts for every ROI, with class
scores differing by less than 2 × 10⁻⁵.

## Software

R (≥ 4.2):

```r
install.packages(c("data.table","dplyr","ggplot2","ggtext","patchwork","mgcv",
                   "vegan","oce","sf","rnaturalearth","readxl","png","jsonlite"))
```

- [`iRfcb`](https://cran.r-project.org/package=iRfcb) 0.10.0. `00b` needs its
  `feature-chain-counts` branch and `hdf5r`.
- [`SHARK4R`](https://github.com/sharksmhi/SHARK4R) for `01`.
- `arrow` (optional) to write the per-ROI tables as Parquet.

Python (≥ 3.9), for `15`, `15b`, `21`, `23d` and `23g`:

```bash
pip install ultralytics numpy pandas matplotlib pillow scipy requests
```

`23d` and `23g` also need a Label Studio server and `ls_preannotate_api.py` from IFCBChainCounts
(set `IFCBCHAINCOUNTS_TOOLS` to its `tools/` directory).

Feature extraction (`00a`, or `00b` with `run_feature_extraction = TRUE`) uses WHOI
[`ifcb-features`](https://github.com/WHOIGit/ifcb-features) through `reticulate`, pinned
to v1.0.0 (install command above). `00a` prints the installed commit and warns if it is
not v1.0.0. To use another Python installation, set `RETICULATE_PYTHON=/path/to/python`.

## Citation

Please cite the paper above, together with `iRfcb` (Torstensson 2025), `SHARK4R` (Lindh
and Torstensson 2026), WHOI `ifcb-features` and the YOLO colony-cell counters of Groves et
al. (2026).

## License

MIT (see [`LICENSE`](LICENSE)).
