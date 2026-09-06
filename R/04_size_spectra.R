# =============================================================================
# Per-sample size-spectrum exponent b (bounded power-law MLE on per-cell biovolume),
# size diversity and median ESD, with depth and season GAMMs.
#
# Reads:  outputs/tangesund/02_sample_meta.csv, data/roi_level/roi_2016.txt.gz,
#         outputs/tangesund/04b_size_profile_raw.csv (optional; panel a falls back to the
#         GAMM profile without it)
# Writes: outputs/tangesund/ 04_size_spectra.csv and the 04_*.png figures
# Deps:   data.table, ggplot2, mgcv, patchwork
# Run:    Rscript R/04_size_spectra.R, then 04b, then this script again for panel a
# =============================================================================
Sys.setlocale("LC_TIME", "C")  # English month names on date axes
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(mgcv); library(patchwork) })
source("R/lib_plotting_theme.R")
source("R/lib_size_spectra.R")

cfg <- list(
  roi_file  = "data/roi_level/roi_2016.txt.gz",
  meta_csv  = "outputs/tangesund/02_sample_meta.csv",
  out_dir   = "outputs/tangesund",
  junk      = c("Debris", "Beads", "Air_bubbles", "Zooplankton"),
  depths    = c(3, 6, 8, 11, 13, 16),
  esd_min   = 7, esd_max = 200, min_cells = 100, boot_B = 400
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(2016)  # bootstrap intervals b_lo/b_hi
v_lo <- (pi/6) * cfg$esd_min^3; v_hi <- (pi/6) * cfg$esd_max^3
meta <- fread(cfg$meta_csv); meta[, date := as.Date(date)]

# one sample -> per-cell (primary) + per-particle spectrum metrics; the bootstrap is
# seeded from the sample id so results do not depend on processing order
seed_from_id <- function(id) { u <- utf8ToInt(id); as.integer(sum(u * seq_along(u)) %% 2147483647L) }
spectrum_metrics <- function(biovol, pcb, ncells, sid) {
  set.seed(seed_from_id(sid))
  inw_c <- pcb >= v_lo & pcb <= v_hi
  n_cell_fit <- sum(ncells[inw_c])
  out <- list(b=NA_real_, b_lo=NA_real_, b_hi=NA_real_, size_div=NA_real_,
              b_part=NA_real_, size_div_part=NA_real_, total_biovol=sum(biovol),
              n_cell_fit=n_cell_fit, n_all=length(biovol), median_esd=NA_real_,
              nbss_r2=NA_real_, reliable=FALSE)
  if (n_cell_fit < cfg$min_cells) return(out)
  ci <- plb_bootstrap(pcb, w = ncells, xmin = v_lo, xmax = v_hi, B = cfg$boot_B)
  out$b <- unname(ci["b"]); out$b_lo <- unname(ci["lo"]); out$b_hi <- unname(ci["hi"])
  out$size_div <- size_diversity(pcb, w = ncells, xmin = v_lo, xmax = v_hi)
  out$b_part <- fit_plb(biovol, xmin = v_lo, xmax = v_hi)$b
  out$size_div_part <- size_diversity(biovol, xmin = v_lo, xmax = v_hi)
  out$median_esd <- median(esd_from_biovol(pcb[inw_c]))
  out$nbss_r2 <- nbss_binned(pcb[inw_c])$r2
  out$reliable <- TRUE
  out
}

# ---- stream the single 2016 file -> per-cell biovolume -----------------------
d <- fread(cfg$roi_file, select = c("sample","biovolume_um3","numBlobs","class_name",
                                    "cell_count","cell_count_resolved"), showProgress = FALSE)
d <- d[numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 & !class_name %in% cfg$junk]
yolo <- !is.na(d$cell_count) & d$cell_count != -1 & !is.na(d$cell_count_resolved)
d[, n_cells := 1L]
d[yolo, n_cells := pmax(1L, round(cell_count_resolved))]
d[, per_cell_biovol := biovolume_um3 / n_cells]
spec <- d[, spectrum_metrics(biovolume_um3, per_cell_biovol, n_cells, sample[1]), by = sample]

dat <- merge(meta, spec, by = "sample", all.x = TRUE)
fwrite(dat, file.path(cfg$out_dir, "04_size_spectra.csv"))
ana <- dat[qc_pass == TRUE & reliable == TRUE]
ana[, `:=`(depth_f = factor(depth, levels = cfg$depths), day = factor(date))]
cat(sprintf("reliable spectra: %d / %d QC-pass | per-cell b median %.2f | median ESD %.1f um\n",
            nrow(ana), sum(dat$qc_pass), median(ana$b), median(ana$median_esd)))

# ---- GAMMs: b, size_div, median_esd ~ s(depth)+s(yday)+ti + day-RE -----------
fit_gamm <- function(yv) tryCatch(
  bam(reformulate(c("s(depth,k=5)","s(yday,k=8)","ti(depth,yday,k=c(5,6))","s(day,bs=\"re\")"), yv),
      data = ana, method = "fREML"),
  error = function(e) { message(yv, " GAMM failed: ", conditionMessage(e)); NULL })
m_b   <- fit_gamm("b")
m_sd  <- fit_gamm("size_div")
m_esd <- fit_gamm("median_esd")

pdepth <- function(m, yv) {
  if (is.null(m)) return(NULL)
  g <- data.frame(depth = seq(3, 16, length = 60), yday = median(ana$yday), day = ana$day[1])
  p <- predict(m, g, se.fit = TRUE, exclude = "s(day)")
  data.table(var = yv, depth = g$depth, fit = p$fit, se = p$se.fit)
}
prof <- rbindlist(list(pdepth(m_b,"Size-spectrum exponent b"), pdepth(m_sd,"Size diversity"),
                       pdepth(m_esd,"Median community ESD (µm)")))

# ---- (1) depth profiles: raw means with a day-clustered SE ----------------------
# values from 04b_size_profile_raw.csv; the GAMM profile (p1_adj) goes to
# 04_spectrum_depth_profiles.png only
raw_csv <- file.path(cfg$out_dir, "04b_size_profile_raw.csv")
if (file.exists(raw_csv)) {
  rawp <- fread(raw_csv)
  rawp[, var := factor(c(median_esd = "Median community ESD (µm)", size_div = "Size diversity",
                         b = "Size-spectrum exponent b")[variable],
                       levels = c("Median community ESD (µm)", "Size diversity", "Size-spectrum exponent b"))]
  setorder(rawp, var, depth)
  p1 <- ggplot(rawp, aes(raw_mean, depth)) +
    geom_errorbar(aes(xmin = raw_mean - se_day, xmax = raw_mean + se_day), orientation = "y", width = 0.4, colour = "steelblue") +
    geom_path(colour = "steelblue") + geom_point(size = 2.2, colour = "steelblue") +
    scale_y_reverse(breaks = cfg$depths) + facet_wrap(~ var, scales = "free_x") +
    theme_paper() +
    labs(x = "Mean of per-sample values (± SE clustered by day)", y = "Depth (m)",
         title = "Vertical structure of the size spectrum")
} else {
  message("04b_size_profile_raw.csv missing: panel (a) falls back to the day-adjusted GAMM profile")
  p1 <- ggplot(prof, aes(fit, depth)) +
    geom_ribbon(aes(xmin = fit - se, xmax = fit + se), alpha = 0.2, orientation = "y") +
    geom_path(colour = "steelblue", linewidth = 0.9) +
    scale_y_reverse(breaks = cfg$depths) + facet_wrap(~ var, scales = "free_x") +
    theme_paper() +
    labs(x = "Day-adjusted value", y = "Depth (m)", title = "Vertical structure of the size spectrum")
}
p1_adj <- ggplot(prof, aes(fit, depth)) +
  geom_ribbon(aes(xmin = fit - se, xmax = fit + se), alpha = 0.2, orientation = "y") +
  geom_path(colour = "steelblue", linewidth = 0.9) +
  scale_y_reverse(breaks = cfg$depths) + facet_wrap(~ var, scales = "free_x") +
  theme_paper() + labs(x = "Day-adjusted value", y = "Depth (m)")

# ---- (2) phenology of slope b by depth --------------------------------------
p2 <- ggplot(ana, aes(date, b, colour = depth_f, linetype = depth_f)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6), se = FALSE, linewidth = 0.8) +
  scale_colour_viridis_d(name = "Depth (m)") +
  scale_linetype_manual(values = c(1, 2, 3, 4, 5, 6), name = "Depth (m)") + theme_paper() +
  labs(x = NULL, y = "Size-spectrum exponent b", title = "Phenology of size-spectrum slope by depth")

# ---- (3) depth x week heatmap of median cell size ----------------------------
hw <- ana[, .(median_esd = mean(median_esd), n = .N), by = .(week, depth_f)]
p3 <- ggplot(hw, aes(factor(week), depth_f, fill = median_esd)) +
  geom_tile(colour = "white") + scale_fill_viridis_c(option = "D", name = "Median\nESD (µm)") +
  scale_y_discrete(limits = rev(levels(ana$depth_f))) + theme_paper() +
  labs(x = "Study week (from 10 Aug 2016)", y = "Depth (m)",
       title = "Median cell size across depth × time")

# ---- (4) b vs median_esd, colored by depth -----------------------------------
p4 <- ggplot(ana, aes(median_esd, b, colour = depth)) +
  geom_point(size = 0.5, alpha = 0.4) + scale_colour_viridis_c(name = "Depth (m)") +
  geom_smooth(method = "gam", formula = y ~ s(x), colour = "firebrick", se = FALSE) +
  theme_paper() +
  labs(x = "Median cell ESD (µm)", y = "Slope b (per-cell)",
       title = "Steeper spectra have smaller median cells")

ggsave(file.path(cfg$out_dir, "04_spectrum_depth_profiles.png"), p1_adj, width = 10, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "04_slope_phenology_depth.png"),   p2, width = 9,  height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "04_medesd_depth_week_heat.png"),  p3, width = 10, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "04_b_vs_medesd.png"),             p4, width = 7,  height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "04_spectra_panel.png"),
       ((p1 + labs(tag = "a")) /
        (p2 + labs(tag = "b") | p3 + labs(tag = "c"))) &
         theme(plot.tag = element_text(face = "bold", size = 15)),
       width = 15, height = 9, dpi = FIG_DPI)

# ---- summary ----------------------------------------------------------------
cat("\n========================= SIZE-SPECTRUM SUMMARY =========================\n")
for (nm in c("b","size_div","median_esd")) {
  m <- get(paste0("m_", ifelse(nm=="median_esd","esd", ifelse(nm=="size_div","sd","b"))))
  if (is.null(m)) next
  ss <- summary(m); st <- ss$s.table
  cat(sprintf("%-11s dev.expl %.1f%% | s(depth) F=%.1f p=%.2g | s(yday) F=%.1f p=%.2g | ti F=%.1f p=%.2g\n",
              nm, 100*ss$dev.expl, st["s(depth)","F"], st["s(depth)","p-value"],
              st["s(yday)","F"], st["s(yday)","p-value"], st["ti(depth,yday)","F"], st["ti(depth,yday)","p-value"]))
}
bp <- ana[, .(b = mean(b), esd = mean(median_esd)), by = depth][order(depth)]
cat(sprintf("slope b surface(3m)=%.2f vs deep(16m)=%.2f | median ESD %0.1f vs %0.1f um\n",
            bp[depth==3,b], bp[depth==16,b], bp[depth==3,esd], bp[depth==16,esd]))
cat(sprintf("table -> 04_size_spectra.csv | panel -> %s/04_spectra_panel.png (plus the individual 04_*.png panels)\n", cfg$out_dir))
cat("===========================================================================\n")
