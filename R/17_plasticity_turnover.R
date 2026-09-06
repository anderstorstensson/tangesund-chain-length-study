# =============================================================================
# Within- vs between-taxon partition of the 3 m - 16 m community-weighted mean
# chain length, per month and overall, with day-bootstrap intervals (Fig. 3).
#
# Reads:  data/roi_level/roi_2016.txt.gz, outputs/tangesund/02_sample_meta.csv
# Writes: outputs/tangesund/ 17_plasticity_turnover.csv, 17_partition_variants.csv, and the
#         17_*.png figures.
# Deps:   data.table, ggplot2, patchwork
# Run:    Rscript R/17_plasticity_turnover.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("R/lib_plotting_theme.R")

cfg <- list(
  roi_file  = "data/roi_level/roi_2016.txt.gz",
  meta_csv  = "outputs/tangesund/02_sample_meta.csv",
  out_dir   = "outputs/tangesund",
  junk      = c("Debris", "Beads", "Air_bubbles", "Zooplankton"),
  surf = 3, deep = 16,
  min_taxon_roi = 5,     # a taxon must have >=5 chains in a (depth,day) cell to count
  nboot = 2000,
  seed  = 1
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
meta <- fread(cfg$meta_csv)
meta[, date := as.Date(date)]
meta[, period := factor(month.abb[as.integer(format(date, "%m"))], levels = c("Aug","Sep","Oct"))]

# ---- chain-eligible ROIs, joined to depth/day/period ------------------------
d <- fread(cfg$roi_file, select = c("sample","class_name","biovolume_um3","numBlobs",
                                    "cell_count","cell_count_resolved"), showProgress = FALSE)
d <- d[numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 & !class_name %in% cfg$junk]
d <- d[!is.na(cell_count) & cell_count != -1]                          # chain-eligible only
d[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1,
                      round(cell_count_resolved), 1L)]
d <- merge(d, meta[qc_pass == TRUE, .(sample, depth, date, period)], by = "sample")
d <- d[depth %in% c(cfg$surf, cfg$deep)]
d[, day := as.integer(date)]

# ---- decomposition on a (depth x taxon) table of {L_i, f_i} ------------------
# CWM = sum_i f_i * L_i. dCWM = surf - deep =
#   sum_i fbar_i (L_i,surf - L_i,deep)   [within / plasticity]
# + sum_i Lbar_i (f_i,surf - f_i,deep)   [between / turnover]
decompose <- function(dat) {
  agg <- dat[, .(L = mean(nchain), nn = .N), by = .(depth, class_name)]
  agg <- agg[nn >= cfg$min_taxon_roi]
  agg[, f := nn / sum(nn), by = depth]
  w <- merge(agg[depth == cfg$surf, .(class_name, L_s = L, f_s = f)],
             agg[depth == cfg$deep, .(class_name, L_d = L, f_d = f)],
             by = "class_name", all = TRUE)
  for (col in c("L_s","L_d","f_s","f_d")) w[is.na(get(col)), (col) := 0]
  w[, `:=`(fbar = (f_s + f_d)/2, Lbar = (L_s + L_d)/2)]
  within  <- sum(w$fbar * (w$L_s - w$L_d))
  between <- sum(w$Lbar * (w$f_s - w$f_d))
  list(cwm_s = sum(w$f_s * w$L_s), cwm_d = sum(w$f_d * w$L_d),
       within = within, between = between, total = within + between)
}

# ---- day-level bootstrap ----------------------------------------------------
boot_ci <- function(dat, nboot = cfg$nboot) {
  days <- unique(dat$day)
  out <- data.table(within = numeric(nboot), between = numeric(nboot), frac_within = numeric(nboot))
  for (b in seq_len(nboot)) {
    bd <- sample(days, length(days), replace = TRUE)
    sub <- dat[.(bd), on = "day", allow.cartesian = TRUE, nomatch = NULL]
    r <- decompose(sub)
    out[b, `:=`(within = r$within, between = r$between,
                frac_within = r$within / (r$within + r$between))]
  }
  out
}

set.seed(cfg$seed)
rows <- list()
for (per in c("All", "Aug", "Sep", "Oct")) {
  dat <- if (per == "All") d else d[period == per]
  if (uniqueN(dat$day) < 3) next
  r  <- decompose(dat)
  bb <- boot_ci(dat)
  qs <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
  rows[[per]] <- data.table(
    period = per, cwm_surf = r$cwm_s, cwm_deep = r$cwm_d, dCWM = r$total,
    within = r$within, between = r$between,
    frac_within = 100 * r$within / r$total,
    within_lo = qs(bb$within)[1], within_hi = qs(bb$within)[2],
    between_lo = qs(bb$between)[1], between_hi = qs(bb$between)[2],
    frac_within_lo = 100 * qs(bb$frac_within)[1], frac_within_hi = 100 * qs(bb$frac_within)[2])
}
res <- rbindlist(rows)
fwrite(res, file.path(cfg$out_dir, "17_plasticity_turnover.csv"))

# ---- partition variants: weak detectors dropped; cell- instead of colony-weighted
weak_taxa <- c("Chaetoceros", "Asterionellopsis", "Pauliella")  # Table S6 "weak" (Melosira is moderate)
d[, is_weak := grepl(paste(weak_taxa, collapse = "|"), class_name)]

decompose_cells <- function(dat) {                 # f_i weighted by cells, not colonies
  agg <- dat[, .(L = mean(nchain), nn = .N, ncell = sum(nchain)), by = .(depth, class_name)]
  agg <- agg[nn >= cfg$min_taxon_roi]
  agg[, f := ncell / sum(ncell), by = depth]
  w <- merge(agg[depth == cfg$surf, .(class_name, L_s = L, f_s = f)],
             agg[depth == cfg$deep, .(class_name, L_d = L, f_d = f)],
             by = "class_name", all = TRUE)
  for (col in c("L_s","L_d","f_s","f_d")) w[is.na(get(col)), (col) := 0]
  w[, `:=`(fbar = (f_s + f_d)/2, Lbar = (L_s + L_d)/2)]
  within <- sum(w$fbar * (w$L_s - w$L_d)); between <- sum(w$Lbar * (w$f_s - w$f_d))
  list(within = within, between = between, total = within + between)
}

v_full <- decompose(d)
v_acc  <- decompose(d[is_weak == FALSE])
v_cell <- decompose_cells(d)
variants <- data.table(
  variant = c("all detectors, colony-weighted (primary)",
              "accurate/moderate detectors only, colony-weighted",
              "all detectors, cell-weighted"),
  dCWM    = c(v_full$total,  v_acc$total,  v_cell$total),
  within  = c(v_full$within, v_acc$within, v_cell$within),
  between = c(v_full$between,v_acc$between,v_cell$between))
variants[, frac_within := 100 * within / dCWM]
fwrite(variants, file.path(cfg$out_dir, "17_partition_variants.csv"))
cat("\n-- within-taxon share under two re-weightings --\n")
print(variants[, .(variant, dCWM = round(dCWM,2), within = round(within,2),
                   between = round(between,2), frac_within = round(frac_within))])
cat("\n-- surface(3m) - deep(16m) chain-length decomposition --\n")
print(res[, .(period, dCWM = round(dCWM,2), within = round(within,2), between = round(between,2),
              frac_within = round(frac_within), CI = sprintf("[%.0f,%.0f]", frac_within_lo, frac_within_hi))])

# ---- figures ----------------------------------------------------------------
long <- melt(res[, .(period, within, between)], id.vars = "period",
             variable.name = "component", value.name = "cells")
long[, component := factor(component, levels = c("within","between"),
                           labels = c("Within-taxon","Between-taxon"))]
long[, period := factor(period, levels = c("All","Aug","Sep","Oct"))]
p1 <- ggplot(long, aes(period, cells, fill = component)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = c("Within-taxon" = "#2c7fb8", "Between-taxon" = "#d95f0e"),
                    name = NULL) +
  theme_paper() +
  labs(x = NULL, y = "Surface − deep chain length (cells)",
       title = "Vertical chain-length gradient is mostly within-taxon")
p2 <- ggplot(res, aes(period, frac_within)) +
  geom_col(fill = "#2c7fb8", width = 0.6) +
  geom_errorbar(aes(ymin = frac_within_lo, ymax = frac_within_hi), width = 0.2) +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey40") +
  scale_x_discrete(limits = c("All","Aug","Sep","Oct")) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_paper() +
  labs(x = NULL, y = "% of gradient that is within-taxon",
       title = "Plasticity share (day-bootstrap 95% CI)")
ggsave(file.path(cfg$out_dir, "17_plasticity_turnover.png"), p1, width = 6.5, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "17_plasticity_share.png"),    p2, width = 5.5, height = 4.5, dpi = FIG_DPI)
# 9.5 x 4.5 in fits a 180 mm double-column figure without rescaling
ggsave(file.path(cfg$out_dir, "17_plasticity_panel.png"),
       ((p1 + labs(tag = "a", title = "Gradient: within- vs between-taxon") +
           theme(legend.position = "bottom", legend.text = element_text(size = 8))) |
        (p2 + labs(tag = "b", title = "Within-taxon share (95% CI)"))) &
         theme(plot.tag = element_text(face = "bold", size = 13)),
       width = 9.5, height = 4.5, dpi = FIG_DPI)

cat("\n===================== PLASTICITY vs TURNOVER SUMMARY =====================\n")
cat(sprintf("Overall: surface−deep Δchain = %.2f cells | within-taxon %.0f%% (95%% CI %.0f–%.0f%%)\n",
            res[period=="All", dCWM], res[period=="All", frac_within],
            res[period=="All", frac_within_lo], res[period=="All", frac_within_hi]))
cat("table -> 17_plasticity_turnover.csv | panel -> 17_plasticity_panel.png (plus the individual 17_*.png panels)\n")
cat("=============================================================================\n")
