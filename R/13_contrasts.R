# =============================================================================
# Surface (3 m) minus deep (16 m) contrasts by period with +-2 SE, the current-speed
# test, and the biovolume-weighting and ESD-window sensitivities.
#
# Reads:  data/roi_level/roi_2016.txt.gz; outputs/tangesund/02_sample_meta.csv,
#         03_chain_metrics.csv, 04_size_spectra.csv, 05_composition.csv, 06_env_matched.csv.
# Writes: outputs/tangesund/13_depth_contrasts.csv, 13_bvw_profile.csv, and
#         13_forest_contrasts.png, 13_chain_vs_current.png,
#         13_slopeb_window_sensitivity.png, 13_contrasts_panel.png.
# Deps:   data.table, mgcv, ggplot2, patchwork; sources R/lib_size_spectra.R.
# Run:    Rscript R/13_contrasts.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")  # English month names on date axes
suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(patchwork) })
source("R/lib_plotting_theme.R")
source("R/lib_size_spectra.R")

cfg <- list(
  roi_file = "data/roi_level/roi_2016.txt.gz", meta = "outputs/tangesund/02_sample_meta.csv",
  out_dir = "outputs/tangesund", junk = c("Debris","Beads","Air_bubbles","Zooplankton"),
  depths = c(3,6,8,11,13,16), esd_windows = c(5,7,10,12,15), esd_max = 200, min_cells = 100
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# assemble the sample-level table (03/04/05 + env)
meta <- fread(cfg$meta)[, .(sample, depth, date, week, yday, qc_pass)]
chn  <- fread("outputs/tangesund/03_chain_metrics.csv")[, .(sample, reliable_chain, mean_chain, frac_cells_chain)]
spc  <- fread("outputs/tangesund/04_size_spectra.csv")[, .(sample, reliable, b, median_esd)]
cmp  <- fread("outputs/tangesund/05_composition.csv"); setnames(cmp, "f_chain diatom", "f_chain_diatom")
env  <- fread("outputs/tangesund/06_env_matched.csv")[, .(sample, insitu_cur, insitu_sal)]
# f_chain_diatom (05, Fig. 5) is not the same quantity as 03's frac_cells_chain
D <- Reduce(function(a,b) merge(a,b,by="sample",all.x=TRUE),
            list(meta, chn, spc, cmp[, .(sample, richness, f_chain_diatom)], env))
D[, `:=`(period = factor(month.abb[as.integer(format(as.Date(date),"%m"))], levels = c("Aug","Sep","Oct")),
         day = factor(date))]
D <- D[qc_pass == TRUE]

# ---- (A) surface(3m)-deep(16m) contrast per period, with CIs -----------------
contrast_depth <- function(yv, dsub) {
  # sample sets: chain-thresholded (chain length), reliable spectra (b, ESD), all QC (share)
  d2 <- D[reliable_chain == TRUE & is.finite(get(yv))]
  if (yv %in% c("b","median_esd")) d2 <- D[reliable == TRUE & is.finite(get(yv))]
  if (yv == "f_chain_diatom")      d2 <- D[is.finite(get(yv))]
  m <- bam(reformulate(c("period","s(depth,by=period,k=5)","s(yday,k=8)","s(day,bs=\"re\")"), yv),
           data = d2, method = "fREML")
  rbindlist(lapply(levels(D$period), function(pp) {
    sub <- d2[period == pp]; if (nrow(sub) < 40) return(NULL)
    g <- data.frame(depth = c(3, 16), period = factor(pp, levels = levels(D$period)),
                    yday = median(sub$yday), day = sub$day[1])
    Xp <- predict(m, g, type = "lpmatrix")
    cvec <- Xp[1, ] - Xp[2, ]                      # surface - deep; shared yday/day terms cancel
    est <- sum(cvec * coef(m)); se <- sqrt(drop(cvec %*% vcov(m) %*% cvec))
    data.table(response = yv, period = pp, est = est, lo = est - 2*se, hi = est + 2*se, n = nrow(sub))
  }))
}
resp <- c(mean_chain = "Chain length (cells)", b = "Slope b", median_esd = "Median ESD (µm)",
          f_chain_diatom = "Chain-diatom share")
contr <- rbindlist(lapply(names(resp), contrast_depth))
contr[, response := factor(resp[response], levels = resp)]
# Aug at the top of the y-axis
contr[, period := factor(period, levels = c("Oct", "Sep", "Aug"))]
cat("== (A) surface(3m) - deep(16m) contrasts (+-2SE), by period ==\n")
print(contr[, .(response, period, est = round(est,3), lo = round(lo,3), hi = round(hi,3),
                sig = ifelse(lo > 0 | hi < 0, "*", ""), n)])
fwrite(contr, file.path(cfg$out_dir, "13_depth_contrasts.csv"))
pA <- ggplot(contr, aes(est, period, colour = lo > 0 | hi < 0)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.25) + geom_point(size = 2.4) +
  facet_wrap(~ response, scales = "free_x") +
  scale_colour_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey60"), guide = "none") +
  theme_paper() + labs(x = "Surface (3 m) - deep (16 m) contrast (+-2SE)", y = NULL,
    title = "Surface−deep contrasts by period")

# ---- (B) current-speed fragmentation test ------------------------------------
db <- D[reliable_chain == TRUE & is.finite(mean_chain) & is.finite(insitu_cur)]
mB <- bam(mean_chain ~ s(insitu_cur, k=6) + s(insitu_sal, k=6) + s(yday, k=8) + s(day, bs="re"),
          data = db, method = "fREML")   # no shrinkage (select = FALSE)
sB <- summary(mB); stB <- sB$s.table
cat(sprintf("\n== (B) fragmentation: chain ~ s(current)+s(sal)+s(yday)+day-RE (n=%d) ==\n", nrow(db)))
cat(sprintf("  s(current) edf=%.2f F=%.1f p=%.3g | s(salinity) edf=%.2f F=%.1f p=%.3g\n",
            stB["s(insitu_cur)","edf"], stB["s(insitu_cur)","F"], stB["s(insitu_cur)","p-value"],
            stB["s(insitu_sal)","edf"], stB["s(insitu_sal)","F"], stB["s(insitu_sal)","p-value"]))
gcur <- data.frame(insitu_cur = seq(quantile(db$insitu_cur,.02), quantile(db$insitu_cur,.98), length=80),
                   insitu_sal = median(db$insitu_sal), yday = median(db$yday), day = db$day[1])
pc <- predict(mB, gcur, se.fit = TRUE, exclude = "s(day)")
pB <- ggplot(data.table(cur=gcur$insitu_cur, fit=pc$fit, se=pc$se.fit), aes(cur, fit)) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se), alpha=0.2) + geom_line(colour="steelblue", linewidth=1) +
  theme_paper() + labs(x="Current speed (cm/s)", y="Mean chain length (partial)",
    title="Chain length vs current speed (fragmentation test)")

# ---- (C) sensitivities: stream roi_2016 --------------------------------------
d <- fread(cfg$roi_file, select = c("sample","biovolume_um3","numBlobs","class_name","cell_count","cell_count_resolved"),
           showProgress = FALSE)
d <- d[numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 & !class_name %in% cfg$junk]
d[, chain_elig := !is.na(cell_count) & cell_count != -1]
d[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1, round(cell_count_resolved), 1L)]
d[, per_cell := biovolume_um3 / nchain]

# (C1) biovolume-weighted chain length per sample
bvw <- d[chain_elig == TRUE, .(mean_chain_bvw = sum(nchain * biovolume_um3) / sum(biovolume_um3),
                               n_chain_roi = .N), by = sample]
bvw <- merge(meta, bvw, by = "sample")[qc_pass == TRUE & n_chain_roi >= 20]
bvw[, day := factor(date)]
mC1 <- bam(mean_chain_bvw ~ s(depth,k=5) + s(yday,k=8) + s(day,bs="re"), data = bvw, method = "fREML")
prof_bvw <- bvw[, .(m = mean(mean_chain_bvw), n = .N, n_days = uniqueN(day)), by = depth][order(depth)]
cat(sprintf("\n== (C1) biovolume-weighted chain length: surface(3m)=%.1f deep(16m)=%.1f cells (s(depth) F=%.1f) ==\n",
            prof_bvw[depth==3,m], prof_bvw[depth==16,m], summary(mC1)$s.table["s(depth)","F"]))
stC1 <- summary(mC1)$s.table
prof_bvw[, `:=`(s_depth_edf = stC1["s(depth)","edf"], s_depth_ref_df = stC1["s(depth)","Ref.df"],
                s_depth_F = stC1["s(depth)","F"], s_depth_p = stC1["s(depth)","p-value"])]
fwrite(prof_bvw, file.path(cfg$out_dir, "13_bvw_profile.csv"))   # 3 m and 16 m values quoted in Note S1

# (C2) slope-b depth gradient across ESD windows
v_hi <- (pi/6) * cfg$esd_max^3
wins <- rbindlist(lapply(cfg$esd_windows, function(emin) {
  v_lo <- (pi/6) * emin^3
  s <- d[per_cell >= v_lo & per_cell <= v_hi,
         .(b = fit_plb(per_cell, w = rep(1,.N), xmin = v_lo, xmax = v_hi)$b, ncell = .N), by = sample]
  s <- merge(meta, s[ncell >= cfg$min_cells], by = "sample")[qc_pass == TRUE]
  s[, `:=`(day = factor(date), emin = emin)]
  m <- bam(b ~ s(depth,k=5) + s(yday,k=8) + s(day,bs="re"), data = s, method = "fREML")
  g <- data.frame(depth = c(3,16), yday = median(s$yday), day = s$day[1])
  Xp <- predict(m, g, type="lpmatrix"); cvec <- Xp[1,]-Xp[2,]
  est <- sum(cvec*coef(m)); se <- sqrt(drop(cvec %*% vcov(m) %*% cvec))
  data.table(emin = emin, est = est, lo = est-2*se, hi = est+2*se, n = nrow(s))
}))
cat("== (C2) slope-b surface-deep contrast across ESD windows ==\n"); print(wins[, .(emin, est=round(est,3), lo=round(lo,3), hi=round(hi,3), n)])
pC <- ggplot(wins, aes(factor(emin), est)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, colour = "steelblue") + geom_point(size = 2.4, colour = "steelblue") +
  theme_paper() + labs(x = "ESD fit window minimum (µm)", y = "Slope b: surface - deep (+-2SE)",
    title = "Slope b across size windows")

ggsave(file.path(cfg$out_dir, "13_forest_contrasts.png"), pA, width = 10, height = 5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "13_chain_vs_current.png"),   pB, width = 7, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "13_slopeb_window_sensitivity.png"), pC, width = 6.5, height = 4.5, dpi = FIG_DPI)
# combined panel: chain length and chain-diatom share only
pA_si <- pA %+% contr[response %in% resp[c("mean_chain", "f_chain_diatom")]]
ggsave(file.path(cfg$out_dir, "13_contrasts_panel.png"), pA_si,
       width = 6.5, height = 3.4, dpi = FIG_DPI)

cat(sprintf("tables -> 13_depth_contrasts.csv | panels -> %s/13_contrasts_panel.png (plus the individual 13_*.png panels)\n", cfg$out_dir))
