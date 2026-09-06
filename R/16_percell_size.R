# =============================================================================
# Per-cell size (biovolume / resolved cell count) by depth, its correlation with
# chain length, and the within-taxon size response (Fig. 4a-c; panel d is 15's gallery).
#
# Reads:  data/roi_level/roi_2016.txt.gz, outputs/tangesund/02_sample_meta.csv,
#         data/functional_group_map.csv; optional 15_chaeto_gallery_inset.png (panel d,
#         left empty if absent) and 16b_within_depth_correlation.csv (panel b r,
#         computed inline if absent).
# Writes: outputs/tangesund/16_percell_metrics.csv, 16_percell_depth_profile.csv,
#         16_within_taxon_size.csv, 16_roi_correlation.csv, 16_percell_depth_profile.png,
#         16_chain_vs_percell.png, 16_within_taxon_size.png, 16_percell_panel.png.
# Deps:   data.table, ggplot2, mgcv, patchwork, ggtext, png.
# Run:    Rscript R/16_percell_size.R   (after 14 and 15; run 16b, then this script again)
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(mgcv); library(patchwork); library(ggtext) })
source("R/lib_plotting_theme.R")

# markdown-italicize a Latin taxon name, leaving spp./sp./cf. upright (rendered by ggtext)
ital_taxon <- function(x) vapply(as.character(x), function(s) {
  toks <- strsplit(s, " ", fixed = TRUE)[[1]]
  paste(ifelse(toks %in% c("spp.", "sp.", "cf."), toks, paste0("*", toks, "*")), collapse = " ")
}, character(1))

# Okabe-Ito colors keyed by taxon, shared with 03_chain_length.R (yellow and black omitted)
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
  min_chain_roi = 20,
  n_taxa_focus  = 6
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
meta <- fread(cfg$meta_csv)
meta[, date := as.Date(date)]

# ---- chain-eligible ROIs with a per-cell biovolume --------------------------
d <- fread(cfg$roi_file, select = c("sample","class_name","biovolume_um3","numBlobs",
                                    "cell_count","cell_count_resolved"), showProgress = FALSE)
d <- d[numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 & !class_name %in% cfg$junk]
d[, chain_elig := !is.na(cell_count) & cell_count != -1]
d[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1,
                      round(cell_count_resolved), 1L)]
ce <- d[chain_elig == TRUE & nchain >= 1]
ce[, percell_bv := biovolume_um3 / nchain]                 # per-cell biovolume within the chain
esd <- function(bv) 2 * (3 / (4 * pi) * bv)^(1/3)          # equiv. spherical diameter (um)
ce[, percell_esd := esd(percell_bv)]
ce <- merge(ce, meta[, .(sample, depth, yday, week, date, qc_pass)], by = "sample")
ce <- ce[qc_pass == TRUE & depth %in% cfg$depths & is.finite(percell_bv) & percell_bv > 0]

# ---- per-sample: chain length and per-cell size on the same samples ----------
s <- ce[, .(n_chain_roi   = .N,
            mean_chain    = mean(nchain),
            log_percell   = mean(log10(percell_bv)),
            percell_esd   = esd(mean(percell_bv))), by = .(sample, depth, yday, week, date)]
s <- s[n_chain_roi >= cfg$min_chain_roi]
s[, day := factor(date)]
fwrite(s, file.path(cfg$out_dir, "16_percell_metrics.csv"))
rs_lab <- cor(s$mean_chain, s$percell_esd)                 # sample-level r, annotated on Fig. 4b
# within-depth r for Fig. 4b: from 16b when present, else computed inline
within_csv <- file.path(cfg$out_dir, "16b_within_depth_correlation.csv")
r_within <- if (file.exists(within_csv)) {
  fread(within_csv)[grepl("^within depth \\(depth", correlation), r]
} else {
  cor(resid(lm(mean_chain ~ factor(depth), s)), resid(lm(percell_esd ~ factor(depth), s)))
}
cat(sprintf("per-cell-size samples: %d  | ROI-level cor(chain, log per-cell bv) = %.3f  | sample-level cor(mean_chain, per-cell ESD) = %.3f\n",
            nrow(s), cor(ce$nchain, log10(ce$percell_bv)), rs_lab))
# ROI-level correlation with a day-cluster bootstrap interval
roi_cor <- function(method) {
  ce[, lbv := log10(percell_bv)]
  r <- cor(ce$nchain, ce$lbv, method = method)
  days <- unique(ce$date); idx <- split(seq_len(nrow(ce)), ce$date); set.seed(1)
  bt <- vapply(seq_len(2000), function(b) {
    ii <- unlist(idx[as.character(sample(days, length(days), replace = TRUE))], use.names = FALSE)
    cor(ce$nchain[ii], ce$lbv[ii], method = method)
  }, numeric(1))
  data.table(method = c(pearson = "Pearson", spearman = "Spearman")[method], r = r,
             ci_lo = quantile(bt, .025), ci_hi = quantile(bt, .975), n = nrow(ce), n_days = length(days),
             note = "cells per chain vs log10 per-cell biovolume, chain-counted ROIs in QC-pass samples at the six depths; day-cluster bootstrap, 2000 draws")
}
rc <- rbind(roi_cor("pearson"), roi_cor("spearman"))
fwrite(rc, file.path(cfg$out_dir, "16_roi_correlation.csv"))
cat(sprintf("  ROI-level r written: Pearson %.3f [%.3f, %.3f], Spearman %.3f [%.3f, %.3f], n = %d ROIs on %d days\n",
            rc$r[1], rc$ci_lo[1], rc$ci_hi[1], rc$r[2], rc$ci_lo[2], rc$ci_hi[2], rc$n[1], rc$n_days[1]))

# ---- (A) depth profile of chain length and per-cell size --------------------
# raw per-sample means with a day-cluster bootstrap SE, as in 03_chain_length.R
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
prof <- s[, .(mean_chain = mean(mean_chain), percell_esd = mean(percell_esd), n = .N,
              n_day = uniqueN(day), se_chain = boot_se_day(mean_chain, day),
              se_esd = boot_se_day(percell_esd, day)), by = depth][order(depth)]
cat("\n-- chain length & per-cell ESD by depth (SE: day-cluster bootstrap) --\n"); print(prof)
fwrite(prof, file.path(cfg$out_dir, "16_percell_depth_profile.csv"))

# ---- GAMM: per-cell size ~ depth x season, day random effect -----------------
m_size <- tryCatch(
  bam(log_percell ~ s(depth, k = 5) + s(yday, k = 8) + ti(depth, yday, k = c(5, 6)) +
        s(day, bs = "re"), data = s, method = "fREML"),
  error = function(e) { message("per-cell GAMM failed: ", conditionMessage(e)); NULL })

# ---- (B) within-taxon per-cell size by depth ---------------------------------
top <- ce[, .N, by = class_name][order(-N)][1:cfg$n_taxa_focus, class_name]
disp <- fread(cfg$fg_map)[, setNames(display_name, class)]
cet <- ce[class_name %in% top]
cet[, taxon := factor(disp[class_name], levels = disp[top])]
cet[, `:=`(day = factor(date), samp = factor(sample))]
# ROI-level model with day and sample REs (many ROIs share a sample)
m_within <- tryCatch(
  bam(log10(percell_bv) ~ taxon + s(depth, k = 5) + s(yday, k = 6) + s(day, bs = "re") +
        s(samp, bs = "re"),
      data = cet, method = "fREML", discrete = TRUE, nthreads = 4),
  error = function(e) { message("within-taxon GAMM failed: ", conditionMessage(e)); NULL })
# day-cluster bootstrap SE of any statistic; wt applies esd(mean(.)) to resampled ROIs
boot_se_stat <- function(value, day, stat_fn = mean, B = 1000, seed = 1) {
  days <- unique(day); nd <- length(days)
  idx_by_day <- split(seq_along(value), day)
  set.seed(seed)
  boot_vals <- vapply(seq_len(B), function(b) {
    d_samp <- sample(days, nd, replace = TRUE)
    stat_fn(value[unlist(idx_by_day[as.character(d_samp)], use.names = FALSE)])
  }, numeric(1))
  sd(boot_vals)
}
wt <- cet[, .(percell_esd = esd(mean(percell_bv)), n = .N, n_day = uniqueN(day),
              se_esd = boot_se_stat(percell_bv, day, function(x) esd(mean(x)))),
          by = .(taxon, depth)][order(taxon, depth)]
fwrite(wt, file.path(cfg$out_dir, "16_within_taxon_size.csv"))
cat("\n-- within-taxon per-cell ESD (um) surface(3m) vs deep(16m) --\n")
print(dcast(wt[depth %in% c(3, 16)], taxon ~ depth, value.var = "percell_esd"))

# ---- figures ----------------------------------------------------------------
# panel (a): two profiles on separate x axes sharing the depth axis
p1a <- ggplot(prof, aes(mean_chain, depth)) +
  geom_errorbar(aes(xmin = mean_chain - se_chain, xmax = mean_chain + se_chain),
                orientation = "y", colour = "steelblue", width = 0.3) +
  geom_path(colour = "steelblue") + geom_point(size = 2.5, colour = "steelblue") +
  scale_y_reverse(breaks = cfg$depths) +
  theme_paper() +
  labs(x = "Mean chain length\n(cells / chain; ± SE clustered by day)", y = "Depth (m)",
       title = "Chain length and per-cell size decline together with depth")
p1b <- ggplot(prof, aes(percell_esd, depth)) +
  geom_errorbar(aes(xmin = percell_esd - se_esd, xmax = percell_esd + se_esd),
                orientation = "y", colour = "firebrick", width = 0.3) +
  geom_path(colour = "firebrick") + geom_point(size = 2.5, colour = "firebrick") +
  scale_y_reverse(breaks = cfg$depths) +
  theme_paper() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(x = "Per-cell ESD\n(µm; ± SE clustered by day)", y = NULL)
p1 <- (p1a + labs(tag = "a")) | p1b

p2 <- ggplot(s, aes(mean_chain, percell_esd, colour = factor(depth))) +
  geom_point(alpha = 0.35, size = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.7) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3, size = 3.6,
           label = sprintf("atop(italic(r) == \"%.2f\"~(all~samples), italic(r) == \"%.2f\"~(within~depth))",
                           rs_lab, r_within), parse = TRUE) +
  scale_colour_viridis_d(name = "Depth (m)") +
  theme_paper() +
  labs(x = "Mean chain length (cells / chain)", y = "Per-cell ESD (µm)",
       title = "Longer chains are built of larger cells (sample level)")

# shape doubles the color so the six taxa stay separable in grayscale
p3 <- ggplot(wt, aes(percell_esd, depth, colour = taxon, shape = taxon)) +
  geom_errorbar(aes(xmin = percell_esd - se_esd, xmax = percell_esd + se_esd),
                orientation = "y", width = 0.3, alpha = 0.6) +
  geom_path() + geom_point(size = 2.2) +
  scale_y_reverse(breaks = cfg$depths) +
  scale_colour_manual(values = oi_taxa, name = NULL, labels = ital_taxon) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8, 3), name = NULL, labels = ital_taxon) +
  theme_paper() +
  theme(legend.text = ggtext::element_markdown()) +
  labs(x = "Per-cell ESD (µm; ± SE clustered by day)", y = "Depth (m)",
       title = "Within-taxon: cells shrink with depth in five of the six chain-formers")

ggsave(file.path(cfg$out_dir, "16_percell_depth_profile.png"),  p1, width = 6,  height = 5,   dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "16_chain_vs_percell.png"),        p2, width = 6,  height = 5,   dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "16_within_taxon_size.png"),       p3, width = 6,  height = 5,   dpi = FIG_DPI)
# panel (d): the Chaetoceros gallery inset from 14 and 15; left empty if missing
inset_png <- file.path(cfg$out_dir, "15_chaeto_gallery_inset.png")
# wrapped in a ggplot so patchwork aligns its title and tag with panel (c)
pd <- if (file.exists(inset_png)) {
        ggplot() +
          annotation_custom(grid::rasterGrob(png::readPNG(inset_png), interpolate = TRUE)) +
          ggtitle("Representative *Chaetoceros* chains, surface vs deep") +
          labs(tag = "d") +
          theme_void(base_size = 11) +
          theme(plot.title = ggtext::element_markdown(size = 13.2, hjust = 0, margin = margin(b = 4)))
      } else {
        message("NOTE: panel (d) inset missing, run R/15_chaeto_gallery.py first; ",
                "building the figure without panel d")
        plot_spacer()
      }
# Fig. 4: a depth profile, b coupling, c within-taxon, d gallery (pd has its own tag)
p3m <- p3 + theme(legend.position = "bottom",
                  legend.key.size = grid::unit(0.8, "lines"),
                  plot.margin = margin(5.5, 5.5, 5.5, 12)) +
       guides(colour = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE))
design <- "AAAABBBB\nCCCCDDDD"
# wrap_plots keeps the nested pair of panel (a) as one cell of the design
panel <- patchwork::wrap_plots(A = p1, B = p2 + labs(tag = "b"), C = p3m + labs(tag = "c"), D = pd,
                               design = design) &
  theme(plot.tag = element_text(face = "bold", size = 15))
ggsave(file.path(cfg$out_dir, "16_percell_panel.png"), panel, width = 13, height = 10, dpi = FIG_DPI)

# ---- summary ----------------------------------------------------------------
cat("\n========================= PER-CELL SIZE SUMMARY =========================\n")
if (!is.null(m_size)) {
  ss <- summary(m_size); st <- ss$s.table
  cat(sprintf("GAMM log_percell ~ s(depth)+s(yday)+ti+day-RE: dev.expl %.1f%%, n=%d days=%d\n",
              100 * ss$dev.expl, nrow(s), uniqueN(s$day)))
  for (tm in rownames(st))
    cat(sprintf("  %-16s edf=%.2f  F=%.1f  p=%.3g\n", tm, st[tm,"edf"], st[tm,"F"], st[tm,"p-value"]))
}
if (!is.null(m_within)) {
  st2 <- summary(m_within)$s.table
  cat(sprintf("Within-taxon ROI GAMM (per-cell bv ~ taxon + s(depth)+s(yday)+day-RE+sample-RE):\n"))
  cat(sprintf("  s(depth) after taxon: edf=%.2f F=%.1f p=%.3g\n",
              st2["s(depth)","edf"], st2["s(depth)","F"], st2["s(depth)","p-value"]))
}
cat(sprintf("per-cell ESD surface(3m)=%.1f vs deep(16m)=%.1f um; chain %.2f vs %.2f cells\n",
            prof[depth == 3, percell_esd], prof[depth == 16, percell_esd],
            prof[depth == 3, mean_chain], prof[depth == 16, mean_chain]))
cat("tables -> 16_percell_metrics.csv, 16_percell_depth_profile.csv, 16_within_taxon_size.csv\n")
cat("=============================================================================\n")
