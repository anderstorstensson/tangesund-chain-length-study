# =============================================================================
# Per-sample chain-length metrics and the chain-length depth profile, phenology,
# per-taxon and distribution figures. Chain-eligible ROI: cell_count != -1.
#
# Reads:  outputs/tangesund/02_sample_meta.csv, data/roi_level/roi_2016.txt.gz,
#         data/functional_group_map.csv.
# Writes: outputs/tangesund/03_chain_metrics.csv, 03_chain_depth_profile.csv,
#         03_colonial_sensitivity.csv (non-diatom colonial classes removed), the 03_*.png
#         panels and 03_chain_panel.png.
# Deps:   data.table, ggplot2, mgcv, patchwork, ggtext.
# Run:    Rscript R/03_chain_length.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")  # English month names on date axes
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(mgcv); library(patchwork); library(ggtext) })
source("R/lib_plotting_theme.R")

# markdown-italicize a Latin taxon name, leaving spp./sp./cf. upright (ggtext)
ital_taxon <- function(x) vapply(as.character(x), function(s) {
  toks <- strsplit(s, " ", fixed = TRUE)[[1]]
  paste(ifelse(toks %in% c("spp.", "sp.", "cf."), toks, paste0("*", toks, "*")), collapse = " ")
}, character(1))

# Okabe-Ito colors keyed by taxon, shared with 16_percell_size.R
oi_taxa <- c(
  "Chaetoceros spp."               = "#009E73",  # bluish green
  "Pseudo-nitzschia spp."          = "#E69F00",  # orange
  "Asterionellopsis glacialis"     = "#0072B2",  # blue
  "Leptocylindrus danicus/minimus" = "#CC79A7",  # reddish purple
  "Skeletonema spp."               = "#56B4E9",  # sky blue
  "Cerataulina pelagica"           = "#D55E00"   # vermillion
)

cfg <- list(
  roi_file  = "data/roi_level/roi_2016.txt.gz",
  meta_csv  = "outputs/tangesund/02_sample_meta.csv",
  fg_map    = "data/functional_group_map.csv",
  out_dir   = "outputs/tangesund",
  junk      = c("Debris", "Beads", "Air_bubbles", "Zooplankton"),
  depths    = c(3, 6, 8, 11, 13, 16),
  min_chain_roi = 20,     # min chain-eligible ROIs for a reliable per-sample chain metric
  n_taxa_focus  = 6       # per-taxon phenology: this many top chain-formers
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
meta <- fread(cfg$meta_csv)
meta[, date := as.Date(date)]; meta[, depth_f := factor(depth, levels = cfg$depths)]

# ---- ROIs: chain-eligible / solitary / junk -----------------------------------
d <- fread(cfg$roi_file, select = c("sample","class_name","biovolume_um3","numBlobs",
                                    "cell_count","cell_count_resolved"), showProgress = FALSE)
d <- d[numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 & !class_name %in% cfg$junk]
d[, chain_elig := !is.na(cell_count) & cell_count != -1]                  # YOLO ran
d[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1,
                      round(cell_count_resolved), 1L)]                    # cells in this ROI
d[(chain_elig) == FALSE, nchain := 1L]                                   # non-chain taxa = 1 cell

# ---- per-sample chain metrics -----------------------------------------------
chn <- d[chain_elig == TRUE, .(
  n_chain_roi = .N,
  mean_chain  = mean(nchain),                    # cells per chain (ROI-weighted)
  med_chain   = as.numeric(median(nchain)),
  p90_chain   = as.numeric(quantile(nchain, 0.90)),
  max_chain   = max(nchain),
  frac_multi  = mean(nchain >= 2),               # share of chain-taxon ROIs that are multicellular
  cells_chain = sum(nchain)                       # cells in chain-forming diatoms
), by = sample]
allc <- d[, .(cells_tot = sum(nchain)), by = sample]                     # all classified cells
s <- merge(allc, chn, by = "sample", all.x = TRUE)
for (col in c("n_chain_roi","mean_chain","med_chain","p90_chain","max_chain","frac_multi","cells_chain"))
  s[is.na(get(col)), (col) := 0]
s[, frac_cells_chain := cells_chain / cells_tot]                         # diatom-chain share of the community
s[, reliable_chain   := n_chain_roi >= cfg$min_chain_roi]

dat <- merge(meta, s, by = "sample", all.x = TRUE)
fwrite(dat, file.path(cfg$out_dir, "03_chain_metrics.csv"))
ana <- dat[qc_pass == TRUE & reliable_chain == TRUE]
ana[, depth_f := factor(depth, levels = cfg$depths)]
cat(sprintf("chain-reliable samples: %d / %d QC-pass (>=%d chain ROIs)\n",
            nrow(ana), sum(dat$qc_pass), cfg$min_chain_roi))
cat(sprintf("overall: mean chain length %.2f cells | frac cells in chain diatoms %.2f | frac multicell %.2f\n",
            mean(ana$mean_chain), mean(ana$frac_cells_chain), mean(ana$frac_multi)))

# ---- sensitivity: non-diatom colonial classes removed -----------------------
# Dinobryon, Uroglenopsis and Peridiniella are chain-counted but not diatoms;
# samples are re-thresholded on the remaining chain-counted ROIs.
nd_pat <- "^Dinobryon|^Uroglenopsis|^Peridiniella"
nd_roi <- d[chain_elig == TRUE & sample %in% meta[qc_pass == TRUE, sample]]
nd     <- nd_roi[!grepl(nd_pat, class_name), .(mean_chain = mean(nchain), n_chain_roi = .N), by = sample]
nd     <- merge(nd[n_chain_roi >= cfg$min_chain_roi], meta[, .(sample, depth)], by = "sample")
nd_prof <- nd[, .(mean_chain = mean(mean_chain)), by = depth][order(depth)]
cat(sprintf("non-diatom colonial classes (%.2f%% of chain-counted ROIs) removed: 3 m mean %.2f (all classes %.2f) | 16 m %.2f (%.2f)\n",
            100 * nd_roi[, mean(grepl(nd_pat, class_name))],
            nd_prof[depth == 3, mean_chain], ana[depth == 3, mean(mean_chain)],
            nd_prof[depth == 16, mean_chain], ana[depth == 16, mean(mean_chain)]))
fwrite(merge(nd_prof[, .(depth, mean_chain_diatoms_only = mean_chain, n_samples_diatoms_only = nd[, .N, by = depth][order(depth), N])],
             ana[, .(mean_chain_all_classes = mean(mean_chain), n_samples = .N), by = depth], by = "depth")[
       , pct_non_diatom_colonial_roi := 100 * nd_roi[, mean(grepl(nd_pat, class_name))]],
       file.path(cfg$out_dir, "03_colonial_sensitivity.csv"))

# ---- GAMM: chain length ~ depth x season, day random effect ------------------
ana[, day := factor(date)]
m_chain <- tryCatch(
  bam(mean_chain ~ s(depth, k = 5) + s(yday, k = 8) + ti(depth, yday, k = c(5, 6)) +
        s(day, bs = "re"), data = ana, method = "fREML"),
  error = function(e) { message("chain GAMM failed: ", conditionMessage(e)); NULL })

# ---- (1) depth profile of chain length ---------------------------------------
# Raw per-sample mean with a day-cluster bootstrap SE of that same mean.
boot_se_day <- function(value, day, B = 2000, seed = 1) {
  days <- unique(day); nd <- length(days)
  idx_by_day <- split(seq_along(value), day)
  set.seed(seed)
  boot_means <- vapply(seq_len(B), function(b) {
    d_samp <- sample(days, nd, replace = TRUE)
    mean(value[unlist(idx_by_day[as.character(d_samp)], use.names = FALSE)])
  }, numeric(1))
  sd(boot_means)
}
prof <- ana[, .(mean = mean(mean_chain), n = .N, n_day = uniqueN(day),
                frac_multi = mean(frac_multi),
                se = boot_se_day(mean_chain, day)), by = depth][order(depth)]
cat("\n-- chain length & multicell fraction by depth (SE: day-cluster bootstrap) --\n"); print(prof)
fwrite(prof, file.path(cfg$out_dir, "03_chain_depth_profile.csv"))
p1 <- ggplot(prof, aes(mean, depth)) +
  geom_path(colour = "steelblue") +
  geom_errorbar(aes(xmin = mean - se, xmax = mean + se), orientation = "y", width = 0.4, colour = "steelblue") +
  geom_point(size = 2.5, colour = "steelblue") +
  geom_text(aes(x = mean + se, label = paste0("n = ", n)), hjust = -0.25, size = 3, colour = "grey30") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  scale_y_reverse(breaks = cfg$depths) +
  theme_paper() +
  labs(x = "Mean chain length (cells / chain; ± SE clustered by day)", y = "Depth (m)",
       title = "Vertical structure of diatom chain length")

# ---- (2) phenology: chain length over the season, by depth -------------------
# geom_smooth's GAM, predicted on a date grid and masked in sampling gaps > 7 days
sm <- ana[, {
  xx <- as.numeric(date)
  g <- mgcv::gam(mean_chain ~ s(x, k = 6), data = data.frame(mean_chain = mean_chain, x = xx), method = "REML")
  grid <- seq(min(xx), max(xx), length.out = 160)
  pr <- predict(g, data.frame(x = grid), se.fit = TRUE)
  dd <- sort(unique(xx)); gaps <- which(diff(dd) > 7)
  in_gap <- rep(FALSE, length(grid))
  for (i in gaps) in_gap <- in_gap | (grid > dd[i] & grid < dd[i + 1])
  .(date = as.Date(grid, origin = "1970-01-01"), fit = ifelse(in_gap, NA_real_, pr$fit),
    se = ifelse(in_gap, NA_real_, pr$se.fit))
}, by = depth_f]
p2 <- ggplot(ana, aes(date, mean_chain, colour = depth_f, fill = depth_f, linetype = depth_f)) +
  geom_ribbon(data = sm, aes(x = date, ymin = fit - se, ymax = fit + se, fill = depth_f),
              inherit.aes = FALSE, alpha = 0.12) +
  geom_line(data = sm, aes(date, fit, colour = depth_f, linetype = depth_f), inherit.aes = FALSE,
            linewidth = 0.8, na.rm = TRUE) +
  geom_rug(sides = "b", alpha = 0.3, length = unit(0.025, "npc"), linetype = 1) +
  scale_colour_viridis_d(name = "Depth (m)") + scale_fill_viridis_d(name = "Depth (m)") +
  scale_linetype_manual(values = c(1, 2, 3, 4, 5, 6), name = "Depth (m)") +
  theme_paper() +
  labs(x = NULL, y = "Mean chain length (cells / chain)",
       title = "Autumn phenology of chain length by depth")

# ---- (3) depth x week heatmap of mean chain length ---------------------------
hw <- ana[, .(mean_chain = mean(mean_chain), n = .N), by = .(week, depth_f)]
p3 <- ggplot(hw, aes(factor(week), depth_f, fill = mean_chain)) +
  geom_tile(colour = "white") +
  scale_fill_viridis_c(option = "C", name = "Mean chain\nlength (cells)") +
  scale_y_discrete(limits = rev(levels(ana$depth_f))) +
  theme_paper() +
  labs(x = "Study week (from 10 Aug 2016)", y = "Depth (m)",
       title = "Chain length across depth × time")

# ---- (4) per-taxon chain-length phenology (top chain-formers) ----------------
# legend display names from the functional-group map
disp <- fread(cfg$fg_map)[, setNames(display_name, class)]
top_taxa <- d[chain_elig == TRUE, .(n = .N), by = class_name][order(-n)][1:cfg$n_taxa_focus, class_name]
tx <- d[chain_elig == TRUE & class_name %in% top_taxa,
        .(mean_chain = mean(nchain), n = .N), by = .(sample, class_name)]
tx <- merge(tx, meta[, .(sample, date, depth, week)], by = "sample")
tx <- tx[n >= 5]                                                         # taxon must have >=5 chains in the sample
tx[, taxon := factor(disp[class_name], levels = disp[top_taxa])]
p4 <- ggplot(tx, aes(date, mean_chain, colour = taxon, fill = taxon, linetype = taxon)) +
  # ±SE ribbons as in panel (c); fitted on the log scale
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE, alpha = 0.12, linewidth = 0.8) +
  scale_colour_manual(values = oi_taxa, name = NULL, labels = ital_taxon) +
  scale_fill_manual(values = oi_taxa, name = NULL, labels = ital_taxon) +
  scale_linetype_manual(values = c(1, 2, 3, 4, 5, 6), name = NULL, labels = ital_taxon) +
  # log y-axis keeps the early Asterionellopsis peak on scale
  scale_y_log10(breaks = c(1, 2, 3, 5, 10, 20)) +
  theme_paper() +
  theme(legend.text = ggtext::element_markdown()) +
  labs(x = NULL, y = "Mean chain length (cells / chain)",
       title = "Per-taxon chain-length phenology")

# ---- (5) chain-length distribution by depth (pooled ROIs) --------------------
dist <- merge(d[chain_elig == TRUE, .(sample, nchain)],
              meta[, .(sample, depth_f)], by = "sample")
p5 <- ggplot(dist[nchain <= 30], aes(nchain, after_stat(density), fill = depth_f)) +
  geom_histogram(binwidth = 1, colour = "white", linewidth = 0.05) +
  facet_wrap(~ depth_f, ncol = 3) + scale_fill_viridis_d(guide = "none") +
  theme_paper() +
  labs(x = "Chain length (cells; truncated at 30)", y = "Density",
       title = "Chain-length distribution by depth")

ggsave(file.path(cfg$out_dir, "03_chain_depth_profile.png"),  p1, width = 5.5, height = 5,   dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "03_chain_phenology_depth.png"), p2, width = 9,   height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "03_chain_depth_week_heat.png"), p3, width = 10,  height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "03_chain_taxon_phenology.png"), p4, width = 9,   height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "03_chain_length_dist_depth.png"), p5, width = 10, height = 6,  dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "03_chain_panel.png"),
       ((p1 + labs(tag = "a") | p3 + labs(tag = "b")) /
        (p2 + labs(tag = "c") | p4 + labs(tag = "d"))) &
         theme(plot.tag = element_text(face = "bold", size = 15)),
       width = 16, height = 9, dpi = FIG_DPI)

# ---- summary ----------------------------------------------------------------
cat("\n========================= CHAIN-LENGTH SUMMARY =========================\n")
if (!is.null(m_chain)) {
  ss <- summary(m_chain)
  cat(sprintf("GAMM mean_chain ~ s(depth)+s(yday)+ti(depth,yday)+day-RE: dev.expl %.1f%%, n=%d days=%d\n",
              100 * ss$dev.expl, nrow(ana), uniqueN(ana$day)))
  st <- ss$s.table
  for (tm in rownames(st))
    cat(sprintf("  %-16s edf=%.2f  F=%.1f  p=%.3g\n", tm, st[tm,"edf"], st[tm,"F"], st[tm,"p-value"]))
}
cat(sprintf("chain length surface(3m)=%.2f vs deep(16m)=%.2f cells (raw depth means)\n",
            prof[depth == 3, mean], prof[depth == 16, mean]))
cat(sprintf("table -> 03_chain_metrics.csv | panel -> %s/03_chain_panel.png (plus the individual 03_*.png panels)\n", cfg$out_dir))
cat("===========================================================================\n")
