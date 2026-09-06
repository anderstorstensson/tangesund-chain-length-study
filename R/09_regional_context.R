# =============================================================================
# Regional context: Tangesund 3 m chain length, cell size and chain-diatom share
# compared with Skagerrak (salinity >= 25) R/V Svea samples from Aug-Oct 2023-2025.
#
# Reads:  data/SHARK_PlanktonImaging_<year>_SMHI/processed_data/data.txt and
#         data/roi_level/roi_<year>.txt.gz (2023-2025; missing years are skipped),
#         outputs/tangesund/03_chain_metrics.csv, 04_size_spectra.csv
# Writes: outputs/tangesund/09_regional_skagerrak_autumn.csv,
#         09_regional_context_chain.png, 09_regional_context_panel.png
# Deps:   data.table, ggplot2, patchwork; sources R/lib_size_spectra.R
# Run:    Rscript R/09_regional_context.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")  # English month names on date axes
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("R/lib_plotting_theme.R")
source("R/lib_size_spectra.R")

cfg <- list(
  roi_dir   = "data/roi_level",
  years     = c(2023, 2024, 2025),                 # Svea years with Skagerrak autumn coverage
  data_glob = "data/SHARK_PlanktonImaging_%d_SMHI/processed_data/data.txt",
  tang_chain= "outputs/tangesund/03_chain_metrics.csv",
  tang_spec = "outputs/tangesund/04_size_spectra.csv",
  out_dir   = "outputs/tangesund",
  junk      = c("Debris","Beads","Air_bubbles","Zooplankton"),
  sal_skagerrak = 25, months = 8:10,
  esd_min = 7, esd_max = 200, min_cells = 100, min_chain_roi = 20
)
v_lo <- (pi/6)*cfg$esd_min^3; v_hi <- (pi/6)*cfg$esd_max^3

# ---- target Svea samples: Skagerrak, autumn, machine-classified ---------------
targ <- rbindlist(lapply(cfg$years, function(yr) {
  f <- sprintf(cfg$data_glob, yr); if (!file.exists(f)) return(NULL)
  e <- fread(f, select = c("SMPNO","IMAGE_VERIFICATION","PSAL_FB","SDATE"), na.strings = c("","NA"))
  u <- unique(e[IMAGE_VERIFICATION == "PredictedByMachine",
                .(sample = SMPNO, salinity = as.numeric(PSAL_FB), date = as.Date(SDATE))])
  u[is.finite(salinity) & salinity >= cfg$sal_skagerrak &
    as.integer(format(date, "%m")) %in% cfg$months][, year := yr][]
}), use.names = TRUE)
cat(sprintf("Svea Skagerrak autumn target samples: %d (years %s)\n",
            nrow(targ), paste(cfg$years, collapse=",")))

# ---- stream roi per year, filter to targets, compute metrics ------------------
svea <- rbindlist(lapply(cfg$years, function(yr) {
  f <- file.path(cfg$roi_dir, sprintf("roi_%d.txt.gz", yr)); if (!file.exists(f)) return(NULL)
  ids <- targ[year == yr, sample]
  d <- fread(f, select = c("sample","class_name","biovolume_um3","numBlobs","cell_count","cell_count_resolved"),
             showProgress = FALSE)
  d <- d[sample %in% ids & numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 & !class_name %in% cfg$junk]
  d[, chain_elig := !is.na(cell_count) & cell_count != -1]
  d[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1, round(cell_count_resolved), 1L)]
  d[(chain_elig) == FALSE, nchain := 1L]
  d[, per_cell_biovol := biovolume_um3 / nchain]
  chn <- d[chain_elig == TRUE, .(n_chain_roi = .N, mean_chain = mean(nchain), frac_multi = mean(nchain >= 2),
                                 cells_chain = sum(nchain)), by = sample]
  tot <- d[, .(cells_tot = sum(nchain),
               n_cell_fit = sum(nchain[per_cell_biovol >= v_lo & per_cell_biovol <= v_hi]),
               b_part = fit_plb(biovolume_um3, xmin = v_lo, xmax = v_hi)$b,
               median_esd = median(esd_from_biovol(per_cell_biovol[per_cell_biovol >= v_lo & per_cell_biovol <= v_hi]))),
           by = sample]
  m <- merge(tot, chn, by = "sample", all.x = TRUE)
  m[, year := yr][]
}), use.names = TRUE, fill = TRUE)
svea <- merge(svea, targ[, .(sample, salinity)], by = "sample")
svea[, frac_cells_chain := cells_chain / cells_tot]
svea <- svea[n_cell_fit >= cfg$min_cells & n_chain_roi >= cfg$min_chain_roi]
fwrite(svea, file.path(cfg$out_dir, "09_regional_skagerrak_autumn.csv"))

# ---- Tangesund surface (3 m) reference ---------------------------------------
tc <- fread(cfg$tang_chain)[qc_pass == TRUE & reliable_chain == TRUE & depth == 3]
ts <- fread(cfg$tang_spec)[qc_pass == TRUE & reliable == TRUE & depth == 3]

cmp <- rbindlist(list(
  data.table(src = "Tangesund 2016 (3 m)", mean_chain = tc$mean_chain, frac_multi = tc$frac_multi,
             frac_cells_chain = tc$frac_cells_chain),
  data.table(src = "Svea Skagerrak autumn", mean_chain = svea$mean_chain, frac_multi = svea$frac_multi,
             frac_cells_chain = svea$frac_cells_chain)), fill = TRUE)
cmp_esd <- rbindlist(list(
  data.table(src = "Tangesund 2016 (3 m)", median_esd = ts$median_esd, b = ts$b_part),
  data.table(src = "Svea Skagerrak autumn", median_esd = svea$median_esd, b = svea$b_part)))

cat("\n-- surface chain length: Tangesund 3 m vs Svea Skagerrak autumn --\n")
print(cmp[, .(n = .N, mean_chain = round(mean(mean_chain, na.rm=TRUE),2),
              frac_multi = round(mean(frac_multi, na.rm=TRUE),2),
              frac_cells_chain = round(mean(frac_cells_chain, na.rm=TRUE),2)), by = src])
print(cmp_esd[, .(median_esd = round(median(median_esd, na.rm=TRUE),1),
                  b_part = round(median(b, na.rm=TRUE),2)), by = src])

# ---- figures: distribution overlays ------------------------------------------
p1 <- ggplot(cmp, aes(mean_chain, fill = src)) +
  geom_density(alpha = 0.45, colour = NA) + scale_fill_manual(values = c("firebrick","steelblue"), guide = "none") +
  theme_paper() + labs(x = "Mean chain length (cells / chain)", y = "Density",
    title = "Surface chain length")
p2 <- ggplot(cmp_esd, aes(median_esd, fill = src)) +
  geom_density(alpha = 0.45, colour = NA) + scale_fill_manual(values = c("firebrick","steelblue"), guide = "none") +
  theme_paper() + labs(x = "Median cell ESD (µm)", y = "Density", title = "Surface cell size")
p3 <- ggplot(cmp, aes(frac_cells_chain, fill = src)) +
  geom_density(alpha = 0.45, colour = NA) + scale_fill_manual(values = c("firebrick","steelblue"), name = NULL) +
  theme_paper() + labs(x = "Chain-diatom share of cells", y = "Density", title = "Chain-diatom dominance")
ggsave(file.path(cfg$out_dir, "09_regional_context_chain.png"), p1, width = 7.5, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "09_regional_context_panel.png"),
       (p1 + labs(tag = "a") | p2 + labs(tag = "b") | p3 + labs(tag = "c")) &
         theme(plot.tag = element_text(face = "bold", size = 15)),
       width = 15, height = 4.5, dpi = FIG_DPI)

cat(sprintf("tables -> 09_regional_skagerrak_autumn.csv | panel -> %s/09_regional_context_panel.png\n", cfg$out_dir))
