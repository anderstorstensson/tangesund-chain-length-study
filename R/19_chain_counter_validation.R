# =============================================================================
# Chain-counter validation: per-detector count accuracy and detection metrics
# (Table S6), and the central contrast with the weak detectors dropped.
#
# Reads:  data/chains_counter_accuracy.csv, data/chains/chains_<detector>_yolo11s/results.csv,
#         data/roi_level/roi_2016.txt.gz, outputs/tangesund/02_sample_meta.csv
# Writes: outputs/tangesund/19_counter_validation.csv, 19_counter_robustness.csv,
#         19_counter_accuracy.png; the Methods numbers are printed
# Deps:   data.table, ggplot2
# Run:    Rscript R/19_chain_counter_validation.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

cfg <- list(
  acc_csv   = "data/chains_counter_accuracy.csv",
  chains_dir= "data/chains",
  roi_file  = "data/roi_level/roi_2016.txt.gz",
  meta_csv  = "outputs/tangesund/02_sample_meta.csv",
  out_dir   = "outputs/tangesund",
  junk      = c("Debris","Beads","Air_bubbles","Zooplankton"),
  depths    = c(3,6,8,11,13,16),
  min_chain_roi = 20,
  # weak detectors: held-out MAE > 0.5 cells (flag in chains_counter_accuracy.csv:
  # ok below 0.16, moderate 0.16-0.5, weak above 0.5)
  weak_pat  = "Chaetoceros|Asterionellopsis|Pauliella"
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# classifier label -> per-genus detector
map_det <- function(cn) data.table::fcase(
  grepl("^Chaetoceros", cn), "chaetoceros",
  grepl("^Skeletonema", cn), "skeletonema",
  grepl("^Pseudo-nitzschia", cn), "pseudo-nitzschia",
  grepl("^Leptocylindrus", cn), "leptocylindrus",
  grepl("^Asterionellopsis", cn), "asterionellopsis",
  grepl("^Cerataulina", cn), "cerataulina_pelagica",
  grepl("^Thalassiosira|^Cyclotella", cn), "thalassiosira",
  grepl("^Guinardia", cn), "guinardia_delicatula",
  grepl("^Thalassionema|^Diatoma", cn), "thalassionema_nitzschioides",
  grepl("^Melosira", cn), "melosira",
  grepl("^Dactyliosolen", cn), "dactyliosolen_fragilissimus",
  grepl("^Peridiniella", cn), "peridiniella_catenata",
  grepl("^Pauliella", cn), "pauliella",
  grepl("^Dinobryon", cn), "dinobryon",
  grepl("^Uroglenopsis", cn), "uroglenopsis",
  grepl("^Eucampia", cn), "eucampia",
  grepl("^Lauderia|^Detonula", cn), "lauderia_detonula",
  default = NA_character_)

# ---- (1)+(2) per-detector validation table ----------------------------------
acc <- fread(cfg$acc_csv)
det_metrics <- function(detector) {
  f <- file.path(cfg$chains_dir, sprintf("chains_%s_yolo11s", detector), "results.csv")
  if (!file.exists(f)) return(data.table(detector = detector))
  r <- fread(f); setnames(r, trimws(names(r)))
  # best.pt = highest fitness (0.1 mAP50 + 0.9 mAP50-95), 20 epochs (patience) before the last row
  best <- r[which.max(0.1 * `metrics/mAP50(B)` + 0.9 * `metrics/mAP50-95(B)`)]
  stopifnot(best$epoch == r[.N]$epoch - 20L)
  pick <- function(nm) { nm <- trimws(nm); if (nm %in% names(best)) as.numeric(best[[nm]]) else NA_real_ }
  data.table(detector = detector,
             box_P   = pick("metrics/precision(B)"), box_R = pick("metrics/recall(B)"),
             mAP50   = pick("metrics/mAP50(B)"),     mAP50_95 = pick("metrics/mAP50-95(B)"))
}
det <- rbindlist(lapply(acc$detector, det_metrics), fill = TRUE)
tab <- merge(acc, det, by = "detector", all.x = TRUE)[order(mae)]
cat("== per-detector validation (Table S6) ==\n")
print(tab[, .(display, eval_n, iou, exact_pct, within1_pct, mae,
              mAP50 = round(mAP50,3), flag)])

w <- tab[flag != "weak"]
cat(sprintf("\nAccurate/moderate detectors (n=%d): exact %.0f-%.0f%%, MAE %.02f-%.02f cells\n",
            nrow(w), min(w$exact_pct), max(w$exact_pct), min(w$mae), max(w$mae)))
cat(sprintf("Weak detectors: %s (exact %.0f-%.0f%%, MAE %.02f-%.02f; over-count, counts approximate)\n",
            paste(tab[flag=="weak", display], collapse=", "),
            min(tab[flag=="weak",exact_pct]), max(tab[flag=="weak",exact_pct]),
            min(tab[flag=="weak",mae]), max(tab[flag=="weak",mae])))

# ---- (3) central contrast without the weak detectors ------------------------
meta <- fread(cfg$meta_csv)
d <- fread(cfg$roi_file, select = c("sample","class_name","numBlobs","biovolume_um3",
                                    "cell_count","cell_count_resolved"), showProgress = FALSE)
d <- d[numBlobs>0 & is.finite(biovolume_um3) & biovolume_um3>0 & !class_name %in% cfg$junk]
d <- d[!is.na(cell_count) & cell_count != -1]
d[, nchain := fifelse(!is.na(cell_count_resolved)&cell_count_resolved>=1, round(cell_count_resolved),1L)]
d <- merge(d, meta[qc_pass==TRUE, .(sample,depth)], by="sample")[depth %in% cfg$depths]
d[, weak := grepl(cfg$weak_pat, class_name)]
frac_weak <- 100*mean(d$weak)

# ---- study abundance per detector + abundance-weighted accuracy --------------
d[, detector := map_det(class_name)]
ab <- d[, .(n_roi_study = .N), by = detector]
ab[, pct_study := 100 * n_roi_study / sum(n_roi_study)]
tab <- merge(tab, ab, by = "detector", all.x = TRUE)
tab[is.na(n_roi_study), `:=`(n_roi_study = 0L, pct_study = 0)]
setorder(tab, -pct_study)
fwrite(tab, file.path(cfg$out_dir, "19_counter_validation.csv"))   # with abundance
pr <- tab[pct_study > 0]                                            # detectors actually present
wt <- pr$pct_study / sum(pr$pct_study)
w_mae   <- sum(wt * pr$mae); w_exact <- sum(wt * pr$exact_pct); w_within1 <- sum(wt * pr$within1_pct)
cat(sprintf("\n== study abundance (chain ROIs) & abundance-weighted accuracy ==\n"))
print(tab[, .(display, taxon_group, form, pct_study = round(pct_study,2), exact_pct, mae, flag)])
cat(sprintf("Detectors exceeding 0.1%% of chain ROIs: %d of 17; the poorest (Pauliella 36%%, Eucampia) are <0.03%%\n",
            nrow(tab[pct_study >= 0.1])))
nd_pct <- tab[taxon_group != "diatom", sum(pct_study)]
cat(sprintf("Non-diatom colonial taxa (Dinobryon, Uroglenopsis, Peridiniella) = %.2f%% of chain ROIs (colonial chrysophytes + one chain-forming dinoflagellate)\n", nd_pct))
cat(sprintf("Abundance-weighted count accuracy across the study: exact %.0f%%, within-1 %.0f%%, MAE %.02f cells\n",
            w_exact, w_within1, w_mae))
cat(sprintf("  (dominated by Chaetoceros = %.0f%% of chain ROIs, itself the weak/over-counting detector)\n",
            tab[detector=="chaetoceros", pct_study]))

prof <- function(dt, lab){
  s <- dt[, .(mc=mean(nchain), .N), by=.(sample,depth)][N>=cfg$min_chain_roi]
  p <- s[, .(mean_chain=mean(mc)), by=depth][order(depth)]
  data.table(set=lab, surf3m=p[depth==3,mean_chain], deep16m=p[depth==16,mean_chain],
             contrast=p[depth==3,mean_chain]-p[depth==16,mean_chain],
             monotonic=all(diff(p$mean_chain) < 0))
}
rob <- rbind(prof(d, "all detectors"), prof(d[weak==FALSE], "accurate/moderate only"))
fwrite(rob, file.path(cfg$out_dir, "19_counter_robustness.csv"))
cat(sprintf("\n== detector robustness (%.0f%% of chain ROIs are from weak detectors) ==\n", frac_weak))
print(rob)

# ---- figure: accuracy vs abundance ------------------------------------------
fig <- tab[pct_study > 0]                                    # detectors present at Tangesund
fig[, display := factor(display, levels = fig[order(mae), display])]
pS1 <- ggplot(fig, aes(mae, display, colour = flag)) +
  geom_segment(aes(x = 0, xend = mae, yend = display), linewidth = 0.4) +
  geom_point(aes(size = pct_study)) +
  geom_text(aes(label = sprintf("%.0f%% exact", exact_pct)), hjust = 0, nudge_x = 0.17, size = 2.8) +
  scale_colour_manual(values = c(ok="#2c7fb8", moderate="#d9a441", weak="#d7301f"), name = "accuracy") +
  scale_size_area(max_size = 9, name = "% of study\nchain ROIs") +
  expand_limits(x = max(fig$mae)*1.35) +
  theme_bw(base_size = 11) +
  labs(x = "Count MAE (cells / chain)", y = NULL,
       title = "Colony-cell-counter accuracy vs study abundance")
ggsave(file.path(cfg$out_dir, "19_counter_accuracy.png"), pS1, width = 8.5, height = 5, dpi = 150)

# ---- Methods numbers ---------------------------------------------------------
cat("\n===================== NUMBERS FOR THE METHODS (DETECTOR VALIDATION) =====================\n")
cat(sprintf("N (validated chains, summed): %d across %d per-taxon detectors\n", sum(tab$eval_n), nrow(tab)))
cat(sprintf("Best detectors (Skeletonema, Pseudo-nitzschia): exact ~97%%, within-1 ~99-100%%, MAE 0.035 cells\n"))
cat(sprintf("Range across accurate/moderate detectors: MAE %.02f-%.02f, exact %.0f-%.0f%%\n",
            min(w$mae), max(w$mae), min(w$exact_pct), max(w$exact_pct)))
cat(sprintf("Weak (over-counting): Chaetoceros exact 63.5%%/MAE 0.57, Asterionellopsis 53.3%%/0.67\n"))
cat(sprintf("Robustness: surface-deep contrast %+.2f (all) -> %+.2f (weak detectors dropped, %.0f%% of ROIs); both monotonic\n",
            rob[set=="all detectors",contrast], rob[set=="accurate/moderate only",contrast], frac_weak))
cat("tables -> 19_counter_validation.csv, 19_counter_robustness.csv | fig -> 19_counter_accuracy.png\n")
cat("============================================================\n")
