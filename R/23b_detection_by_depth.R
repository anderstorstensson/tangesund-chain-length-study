# =============================================================================
# Unclassified and zero-detection rates by depth (A), an upper limit from the
# unclassified ROIs (C), and the Chaetoceros contrast with the size-tertile count
# biases of 21 applied under three assumptions for chains below the validated range (B).
#
# Reads:  data/roi_level/roi_2016.txt.gz, outputs/tangesund/02_sample_meta.csv,
#         outputs/tangesund/21_chaeto_size_bias.csv
# Writes: outputs/tangesund/ 23b_detection_by_depth.csv, 23b_unclassified_by_depth.csv,
#         23b_chaeto_bias_bound.csv
# Deps:   data.table
# Run:    Rscript R/23b_detection_by_depth.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))
cfg <- list(roi_file = "data/roi_level/roi_2016.txt.gz",
            meta_csv = "outputs/tangesund/02_sample_meta.csv",
            bias_csv = "outputs/tangesund/21_chaeto_size_bias.csv",
            out_dir = "outputs/tangesund", micron = 1/2.77,          # um per pixel
            junk = c("Debris", "Beads", "Air_bubbles", "Zooplankton"))
meta <- fread(cfg$meta_csv)[qc_pass == TRUE, .(sample, depth, date)]
d <- fread(cfg$roi_file, select = c("sample","class_name","class_name_auto","numBlobs","biovolume_um3","cell_count",
                                    "cell_count_resolved","BoundingBox_xwidth","BoundingBox_ywidth"),
           showProgress = FALSE)
d <- merge(d, meta, by = "sample")
d[, `:=`(unclassified = class_name == "unclassified",
         junk = class_name %in% cfg$junk,
         elig = !is.na(cell_count) & cell_count != -1)]
# zero detections are delivered as cell_count == 0 with the resolved count already set to 1
d[, zero_det := elig & cell_count == 0]
d[, chaeto := grepl("^Chaetoceros", class_name)]
cat(sprintf("ROIs in QC samples: %s | unclassified %.1f%% of non-junk | chain-eligible %s | zero-detection %.2f%% of eligible\n",
            format(nrow(d), big.mark = ","), 100 * d[junk == FALSE, mean(unclassified)],
            format(sum(d$elig), big.mark = ","), 100 * d[elig == TRUE, mean(zero_det)]))
cat(sprintf("  zero-detection encodings among eligible: cell_count == 0: %d | resolved NA: %d | resolved < 1: %d\n",
            d[elig == TRUE, sum(cell_count == 0, na.rm = TRUE)], d[elig == TRUE, sum(is.na(cell_count_resolved))],
            d[elig == TRUE, sum(cell_count_resolved < 1, na.rm = TRUE)]))

# ---- A: censoring by depth -----------------------------------------------------
byd <- d[junk == FALSE, .(rois = .N, unclassified_pct = 100 * mean(unclassified),
                  chain_eligible = sum(elig), zero_detection_pct = 100 * mean(zero_det[elig]),
                  chaeto_rois = sum(chaeto & elig),
                  chaeto_zero_detection_pct = 100 * mean(zero_det[chaeto & elig])), by = depth][order(depth)]
tot <- d[junk == FALSE, .(depth = NA_integer_, rois = .N, unclassified_pct = 100 * mean(unclassified),
                  chain_eligible = sum(elig), zero_detection_pct = 100 * mean(zero_det[elig]),
                  chaeto_rois = sum(chaeto & elig),
                  chaeto_zero_detection_pct = 100 * mean(zero_det[chaeto & elig]))]
A <- rbind(byd, tot)
fwrite(A, file.path(cfg$out_dir, "23b_detection_by_depth.csv"))
cat("\n-- censoring by depth (23b_detection_by_depth.csv) --\n")
print(A[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

# ---- C: unclassified ROIs by highest-probability class -----------------------------
# upper limit: each unclassified ROI whose top class is a chain-former gets that
# class's 3 m mean chain length
d[, nchain := fifelse(elig & !is.na(cell_count_resolved) & cell_count_resolved >= 1, round(cell_count_resolved), 1)]
elig_classes <- unique(d[elig == TRUE, class_name])
u <- d[unclassified == TRUE & junk == FALSE & !class_name_auto %in% cfg$junk & depth %in% c(3, 16)]
u[, chain_former_top := class_name_auto %in% elig_classes]
L_surf <- d[elig == TRUE & depth == 3, .(L_surf = mean(nchain)), by = .(class_name_auto = class_name)]
add <- merge(u[chain_former_top == TRUE, .N, by = .(depth, class_name_auto)], L_surf, by = "class_name_auto")
add <- add[, .(n_added = sum(N), cells_added = sum(N * L_surf)), by = depth]
pooled <- d[elig == TRUE & depth %in% c(3, 16), .(chain_rois = .N, cells = sum(nchain)), by = depth]
C <- u[, .(unclassified = .N, chain_former_top = sum(chain_former_top),
           chaeto_chain_top = sum(class_name_auto == "Chaetoceros_spp_chain")), by = depth]
C <- Reduce(function(a, b) merge(a, b, by = "depth"), list(C, pooled, add))[order(depth)]
C[, `:=`(chain_former_top_pct = 100 * chain_former_top / unclassified,
         pooled_mean = cells / chain_rois,
         bounded_mean = (cells + cells_added) / (chain_rois + n_added))]
fwrite(C, file.path(cfg$out_dir, "23b_unclassified_by_depth.csv"))
cat("\n-- unclassified ROIs by top class, and the bound (23b_unclassified_by_depth.csv) --\n")
print(C[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat(sprintf("  pooled contrast %.2f -> %.2f cells if every unclassified chain-former-like ROI had held its class's 3 m mean\n",
            C[depth == 3, pooled_mean] - C[depth == 16, pooled_mean],
            C[depth == 3, bounded_mean] - C[depth == 16, bounded_mean]))

# ---- B: tertile-bias bound on the Chaetoceros contrast ---------------------------
bias <- fread(cfg$bias_csv)
edges <- as.numeric(unlist(strsplit(bias$esd_um_range, "-")))
lo_support <- min(edges); b_small <- bias[bin == "percell-small", mean_bias]
b_mid <- bias[bin == "percell-mid", mean_bias]; b_large <- bias[bin == "percell-large", mean_bias]
cut_mid <- edges[2]; cut_large <- edges[4]                  # small | mid | large tertile boundaries
ch <- d[chaeto & elig & numBlobs > 0 & is.finite(biovolume_um3) & biovolume_um3 > 0 &
        class_name == "Chaetoceros_spp_chain"]
ch[, nchain := fifelse(!is.na(cell_count_resolved) & cell_count_resolved >= 1, round(cell_count_resolved), 1L)]
# footprint per-cell ESD (bounding-box area / count), as in 21
ch[, fp_esd := 2 * sqrt((BoundingBox_xwidth * BoundingBox_ywidth) / nchain / pi) * cfg$micron]
ch[, tertile := fcase(fp_esd < lo_support, "below support", fp_esd <= cut_mid, "small",
                      fp_esd <= cut_large, "mid", default = "large")]
assign_bias <- function(t, out_of_support) fcase(t == "below support", out_of_support, t == "small", b_small,
                                                 t == "mid", b_mid, t == "large", b_large)
scen <- list(`smallest tertile (−0.05)` = b_small, `no bias (0)` = 0, `mid tertile (+0.17)` = b_mid)
B <- rbindlist(lapply(names(scen), function(s) {
  ch[, bias := assign_bias(tertile, scen[[s]])]
  r <- ch[, .(n = .N, pct_below_support = 100 * mean(tertile == "below support"),
              mean_bias = mean(bias), mean_chain = mean(nchain), mean_chain_corrected = mean(nchain - bias)),
          by = depth][order(depth)]
  r[, scenario_out_of_support := s][]
}))
fwrite(B, file.path(cfg$out_dir, "23b_chaeto_bias_bound.csv"))
cat("\n-- tertile-bias bound (23b_chaeto_bias_bound.csv) --\n")
print(B[, .(scenario_out_of_support, depth, n, pct_below_support = round(pct_below_support, 1),
            mean_bias = round(mean_bias, 3), mean_chain = round(mean_chain, 2),
            corrected = round(mean_chain_corrected, 2))])
for (s in names(scen)) {
  r <- B[scenario_out_of_support == s]
  cat(sprintf("  %-26s surface-minus-deep: raw %+.2f | corrected %+.2f | differential bias %+.3f cells\n", s,
              r[depth == 3, mean_chain] - r[depth == 16, mean_chain],
              r[depth == 3, mean_chain_corrected] - r[depth == 16, mean_chain_corrected],
              r[depth == 3, mean_bias] - r[depth == 16, mean_bias]))
}
