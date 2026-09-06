# =============================================================================
# Functional-group shares and diversity (richness, rarefied richness, Shannon) by
# depth and season, weighted by cell abundance.
#
# Reads:  outputs/tangesund/02_sample_meta.csv, data/roi_level/roi_2016.txt.gz,
#         data/functional_group_map.csv (from 01).
# Writes: outputs/tangesund/05_composition.csv, 05_richness_rarefied_profile.csv,
#         05_surface_weekly_composition.csv, the 05_*.png panels and 05_composition_panel.png.
# Deps:   data.table, ggplot2, mgcv, patchwork; vegan if it is installed.
# Run:    Rscript R/05_composition.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")  # English month names on date axes
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(mgcv); library(patchwork) })
source("R/lib_plotting_theme.R")

cfg <- list(
  roi_file  = "data/roi_level/roi_2016.txt.gz",
  meta_csv  = "outputs/tangesund/02_sample_meta.csv",
  fg_map    = "data/functional_group_map.csv",
  out_dir   = "outputs/tangesund",
  junk      = c("Debris", "Beads", "Air_bubbles", "Zooplankton", "unclassified", "mixed", "Unicells"),
  depths    = c(3, 6, 8, 11, 13, 16)
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
meta <- fread(cfg$meta_csv); meta[, date := as.Date(date)]

# ---- functional-group map (from 01; do not infer groups from taxon names) ----
fg <- fread(cfg$fg_map)[, .(class_name = class, grp = group)]
grp_levels <- c("chain diatom","diatom (other)","dinoflagellate","ciliate","cryptophyte",
                "other flagellate/nano","cyanobacteria","heterotroph/other","other")
# display labels (Title Case, plural) for legends; internal keys stay lowercase
grp_disp <- c("chain diatom"="Chain-forming diatoms", "diatom (other)"="Diatoms (other)",
              "dinoflagellate"="Dinoflagellates", "ciliate"="Ciliates", "cryptophyte"="Cryptophytes",
              "other flagellate/nano"="Other flagellates/nano", "cyanobacteria"="Cyanobacteria",
              "heterotroph/other"="Heterotrophs/other", "other"="Other")
# Okabe-Ito colors; catch-all groups in gray
grp_pal <- setNames(c("#009E73","#56B4E9","#D55E00","#0072B2","#CC79A7",
                      "#E69F00","#000000","#999999","#cccccc"), unname(grp_disp))

# ---- ROIs -> per (sample, taxon) cell abundance and group --------------------
d <- fread(cfg$roi_file, select = c("sample","class_name","biovolume_um3","numBlobs",
                                    "cell_count","cell_count_resolved"), showProgress = FALSE)
d <- d[numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 & !class_name %in% cfg$junk]
d[, ncell := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1, round(cell_count_resolved), 1L)]
d <- merge(d, fg, by = "class_name", all.x = TRUE)
d[is.na(grp) | !grp %in% grp_levels, grp := "other"]
d[, grp := factor(grp, levels = grp_levels)]

cat("-- functional-group map (cells, share) --\n")
gm <- d[, .(cells = sum(ncell)), by = grp][order(-cells)][, share := round(cells/sum(cells), 3)][]
print(gm)

# ---- per-sample: group shares + diversity ------------------------------------
gs <- dcast(d[, .(cells = sum(ncell)), by = .(sample, grp)], sample ~ grp,
            value.var = "cells", fill = 0, drop = FALSE)
gcols <- grp_levels[grp_levels %in% names(gs)]
gs[, tot := rowSums(.SD), .SDcols = gcols]
for (g in gcols) gs[, (paste0("f_", g)) := get(g) / tot]

# Hurlbert rarefied richness at `size` cells; vegan::rarefy if installed, else the
# identical closed form
have_vegan <- requireNamespace("vegan", quietly = TRUE)
rarefy_richness <- function(n_vec, size) {
  if (sum(n_vec) < size) return(NA_real_)
  if (have_vegan) return(as.numeric(suppressWarnings(vegan::rarefy(round(n_vec), sample = size))))
  sum(1 - exp(lchoose(sum(n_vec) - n_vec, size) - lchoose(sum(n_vec), size)))
}
div <- d[, {
  a <- .SD[, .(cells = sum(ncell)), by = class_name][cells > 0]
  p <- a$cells / sum(a$cells)
  .(richness = nrow(a), rich_rare100 = rarefy_richness(a$cells, 100),
    rich_rare200 = rarefy_richness(a$cells, 200),
    shannon = -sum(p * log(p)), cells_tot = sum(a$cells))
}, by = sample]

dat <- Reduce(function(a,b) merge(a,b,by="sample",all.x=TRUE), list(meta, gs, div))
fwrite(dat, file.path(cfg$out_dir, "05_composition.csv"))
ana <- dat[qc_pass == TRUE]
ana[, `:=`(depth_f = factor(depth, levels = cfg$depths), day = factor(date),
           cd_share = `f_chain diatom`)]   # syntactic-name copy for the GAMM; not named f_* so melt skips it
# weekly surface (3 m) composition and diversity
wk_surf <- ana[depth == 3, .(n = .N, n_days = uniqueN(date), f_chain_diatom = mean(`f_chain diatom`),
                             f_dinoflagellate = mean(f_dinoflagellate), shannon = mean(shannon)), by = week][order(week)]
fwrite(wk_surf, file.path(cfg$out_dir, "05_surface_weekly_composition.csv"))
cat(sprintf("\ncomposition samples (QC-pass): %d | median richness %d, Shannon %.2f\n",
            nrow(ana), as.integer(median(ana$richness)), median(ana$shannon)))

# ---- GAMMs: richness, shannon, chain-diatom share vs depth x season ----------
fit_gamm <- function(yv) tryCatch(
  bam(reformulate(c("s(depth,k=5)","s(yday,k=8)","ti(depth,yday,k=c(5,6))","s(day,bs=\"re\")"), yv),
      data = ana, method = "fREML"),
  error = function(e) { message(yv, " GAMM failed: ", conditionMessage(e)); NULL })
m_rich <- fit_gamm("richness"); m_shan <- fit_gamm("shannon"); m_cd <- fit_gamm("cd_share")
m_rrare <- fit_gamm("rich_rare100")

# ---- (1) depth profiles of group shares (season-pooled raw means) ------------
fshare <- melt(ana, id.vars = c("depth", "day"), measure.vars = patterns("^f_"),
               variable.name = "grp", value.name = "frac")
fshare[, grp := sub("^f_", "", grp)]
plotted <- fshare[, mean(frac, na.rm = TRUE), by = grp][V1 >= 0.003, grp]  # drop negligible groups
fshare <- fshare[grp %in% plotted]
fshare[, grp := factor(grp_disp[grp], levels = grp_disp[grp_levels[grp_levels %in% plotted]])]
# raw per-sample mean with a day-cluster bootstrap SE, as in 03_chain_length.R
boot_se_day <- function(value, day, B = 2000, seed = 1) {
  keep <- !is.na(value); value <- value[keep]; day <- day[keep]
  days <- unique(day); nd <- length(days)
  idx_by_day <- split(seq_along(value), day)
  set.seed(seed)
  boot_means <- vapply(seq_len(B), function(b) {
    d_samp <- sample(days, nd, replace = TRUE)
    mean(value[unlist(idx_by_day[as.character(d_samp)], use.names = FALSE)])
  }, numeric(1))
  sd(boot_means)
}
prof <- fshare[, .(frac = mean(frac, na.rm = TRUE), se = boot_se_day(frac, day)), by = .(depth, grp)]
setorder(prof, grp, depth)  # geom_path draws in row order: sort by depth, not by first appearance
p1 <- ggplot(prof, aes(frac, depth, colour = grp, shape = grp)) +
  geom_errorbar(aes(xmin = frac - se, xmax = frac + se), orientation = "y", width = 0.3, alpha = 0.6) +
  geom_path(linewidth = 0.8) + geom_point(size = 2) +
  scale_y_reverse(breaks = cfg$depths) + scale_colour_manual(values = grp_pal, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8, 3, 4, 5, 6)[seq_along(levels(prof$grp))], name = NULL) +
  theme_paper() +
  labs(x = "Mean fraction of classified cells (± SE clustered by day)", y = "Depth (m)",
       title = "Functional-group composition by depth")

# ---- (2) phenology: group shares over the season (surface 3 m) ---------------
surf <- melt(ana[depth == 3], id.vars = "date", measure.vars = patterns("^f_"),
             variable.name = "grp", value.name = "frac")
surf[, grp := sub("^f_", "", grp)]
surf <- surf[grp %in% plotted]
surf[, grp := factor(grp_disp[grp], levels = grp_disp[grp_levels[grp_levels %in% plotted]])]
# quasibinomial GAM on a date grid, masked in sampling gaps > 7 days (as in 03)
sm_c <- surf[, {
  xx <- as.numeric(date)
  g <- mgcv::gam(frac ~ s(x, k = 6), family = quasibinomial(),
                 data = data.frame(frac = frac, x = xx), method = "REML")
  grid <- seq(min(xx), max(xx), length.out = 160)
  pr <- predict(g, data.frame(x = grid), se.fit = TRUE)          # logit scale
  dd <- sort(unique(xx)); gaps <- which(diff(dd) > 7)
  in_gap <- rep(FALSE, length(grid))
  for (i in gaps) in_gap <- in_gap | (grid > dd[i] & grid < dd[i + 1])
  .(date = as.Date(grid, origin = "1970-01-01"),
    fit = ifelse(in_gap, NA_real_, plogis(pr$fit)),
    lo  = ifelse(in_gap, NA_real_, plogis(pr$fit - pr$se.fit)),
    hi  = ifelse(in_gap, NA_real_, plogis(pr$fit + pr$se.fit)))
}, by = grp]
p2 <- ggplot(sm_c, aes(date, fit, colour = grp, fill = grp)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA, na.rm = TRUE) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = grp_pal, name = NULL) + scale_fill_manual(values = grp_pal, name = NULL) +
  theme_paper() +
  labs(x = NULL, y = "Fraction of cells (3 m)", title = "Compositional phenology at the surface")

# ---- (3) diversity depth x week --------------------------------------------
hw <- ana[, .(shannon = mean(shannon), richness = mean(richness), n = .N), by = .(week, depth_f)]
p3 <- ggplot(hw, aes(factor(week), depth_f, fill = shannon)) +
  geom_tile(colour = "white") + scale_fill_viridis_c(option = "D", name = "Shannon\n(ln)") +
  scale_y_discrete(limits = rev(levels(ana$depth_f))) + theme_paper() +
  labs(x = "Study week (from 10 Aug 2016)", y = "Depth (m)",
       title = "Diversity across depth × time")

# ---- (4) chain-diatom share depth profile -----------------------------------
p4 <- ggplot(ana, aes(`f_chain diatom`, depth)) +
  geom_point(size = 0.4, alpha = 0.25, colour = "grey50") +
  stat_summary(aes(group = depth), fun = mean, geom = "point", colour = "firebrick", size = 2.5, orientation = "y") +
  scale_y_reverse(breaks = cfg$depths) + theme_paper() +
  labs(x = "Chain-forming diatom share of cells", y = "Depth (m)",
       title = "Chain-diatom dominance by depth")

ggsave(file.path(cfg$out_dir, "05_group_composition_depth.png"), p1, width = 8, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "05_composition_phenology_surface.png"), p2, width = 9, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "05_diversity_depth_week.png"), p3, width = 10, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "05_chain_diatom_share_depth.png"), p4, width = 6, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "05_composition_panel.png"),
       ((p1 + labs(tag = "a") | p3 + labs(tag = "b")) /
        (p2 + labs(tag = "c"))) &
         theme(plot.tag = element_text(face = "bold", size = 15)),
       width = 16, height = 9, dpi = FIG_DPI)

# ---- (5) raw vs rarefied richness and Shannon by depth ------------------------
rprof <- ana[, .(cells_tot = median(cells_tot, na.rm = TRUE),
                 raw = mean(richness, na.rm = TRUE),
                 rare100 = mean(rich_rare100, na.rm = TRUE),
                 shannon = mean(shannon, na.rm = TRUE)),
             by = depth][order(depth)]
fwrite(rprof, file.path(cfg$out_dir, "05_richness_rarefied_profile.csv"))

# ---- summary ----------------------------------------------------------------
cat("\n========================= COMPOSITION SUMMARY =========================\n")
for (nm in c("richness","rich_rare100","shannon")) {
  m <- get(switch(nm, richness="m_rich", rich_rare100="m_rrare", shannon="m_shan")); if (is.null(m)) next
  ss <- summary(m); st <- ss$s.table
  cat(sprintf("%-13s dev.expl %.1f%% | s(depth) F=%.1f p=%.2g | s(yday) F=%.1f p=%.2g | ti F=%.1f p=%.2g\n",
              nm, 100*ss$dev.expl, st["s(depth)","F"], st["s(depth)","p-value"],
              st["s(yday)","F"], st["s(yday)","p-value"], st["ti(depth,yday)","F"], st["ti(depth,yday)","p-value"]))
}
cds <- ana[, .(cd = mean(`f_chain diatom`), rich = mean(richness)), by = depth][order(depth)]
cat(sprintf("chain-diatom cell share surface(3m)=%.2f vs deep(16m)=%.2f | raw richness %.0f vs %.0f\n",
            cds[depth==3,cd], cds[depth==16,cd], cds[depth==3,rich], cds[depth==16,rich]))
cat(sprintf("Effort-controlled: rarefied@100 %.1f vs %.1f | Shannon %.2f vs %.2f (surface vs deep)\n",
            rprof[depth==3,rare100], rprof[depth==16,rare100], rprof[depth==3,shannon], rprof[depth==16,shannon]))
cat(sprintf("table -> 05_composition.csv | panel -> %s/05_composition_panel.png (plus the individual 05_*.png panels)\n", cfg$out_dir))
cat("=========================================================================\n")
