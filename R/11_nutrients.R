# =============================================================================
# Bottle nutrients and size-fractionated chlorophyll matched to the IFCB samples
# (Tangesund station casts, nearest date within 4 days, interpolated to sample depth).
#
# Reads:  data/tångesund_water_samples/Results ... .xlsx (sheet "Data"), or its public copy
#         skagerrak_chemical_observations.csv (doi:10.5281/zenodo.18378122) in the same folder;
#         outputs/tangesund/ 02_sample_meta.csv, 03_chain_metrics.csv, 05_composition.csv,
#         06_env_matched.csv (optional, buoy cross-check).
# Writes: outputs/tangesund/ 11_nutrients_matched.csv, 11_depth_descriptives.csv, and the
#         11_*.png figures.
# Deps:   readxl, data.table, mgcv, ggplot2, patchwork
# Run:    Rscript R/11_nutrients.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")  # English month names on date axes
suppressPackageStartupMessages({ library(readxl); library(data.table); library(mgcv); library(ggplot2); library(patchwork) })
source("R/lib_plotting_theme.R")

cfg <- list(
  xlsx    = "data/tångesund_water_samples/Results Tångesund + Skagerak version 1 Oct 2019 BK.xlsx",
  csv     = "data/tångesund_water_samples/skagerrak_chemical_observations.csv",  # Zenodo doi:10.5281/zenodo.18378122
  meta    = "outputs/tangesund/02_sample_meta.csv",
  out_dir = "outputs/tangesund",
  match_tol_days = 4, depths = c(3,6,8,11,13,16)
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
num <- function(v) suppressWarnings(as.numeric(v))

# ---- read bottle data (row 2 header), Tangesund only -------------------------
# The Zenodo CSV is the same table with nutrients and chl a rounded to 2 dp.
b <- if (file.exists(cfg$xlsx)) {
  as.data.table(suppressMessages(read_excel(cfg$xlsx, sheet = "Data", skip = 1)))
} else if (file.exists(cfg$csv)) {
  fread(cfg$csv, sep = ";", skip = 1, encoding = "UTF-8", na.strings = c("", "NA"), colClasses = "character")
} else stop("bottle data not found: download skagerrak_chemical_observations.csv from ",
            "https://doi.org/10.5281/zenodo.18378122 into data/tångesund_water_samples/")
setnames(b, make.unique(names(b)))
b <- b[grepl("Tång|Tang", station_name, ignore.case = TRUE)]
b[, date := as.Date(date)]
for (v in c("depth","salinityDown","tempDown","oxyDown","chla","chla20um","chla3um",
            "PO4","NO23","SiO3","NH4","NPratio")) b[[v]] <- num(b[[v]])
b <- b[is.finite(depth) & is.finite(date)]
b[, din := NO23 + fifelse(is.finite(NH4), NH4, 0)]          # DIN = NOx + NH4 (NH4 missing -> 0)
b[, si_din := SiO3 / pmax(din, 1e-3)]
cat(sprintf("Tangesund bottle casts: %d rows | %d dates %s..%s | depths %s\n",
            nrow(b), uniqueN(b$date), min(b$date), max(b$date),
            paste(sort(unique(b$depth)), collapse=",")))

# ---- match IFCB samples to nearest bottle date, interpolate to sample depth ---
meta <- fread(cfg$meta); meta[, date := as.Date(date)]
bdates <- data.table(bdate = sort(unique(b$date)))[, dk := bdate]
m <- meta[, .(sample, depth, date)][, dk := date]
setkey(bdates, dk); setkey(m, dk)
m <- bdates[m, roll = "nearest"]
m[, dday := abs(as.numeric(date - bdate))]
m <- m[dday <= cfg$match_tol_days]
cat(sprintf("IFCB samples with a bottle match (<=%d d): %d / %d\n",
            cfg$match_tol_days, nrow(m), nrow(meta)))

bvars <- c("tempDown","salinityDown","oxyDown","chla","chla20um","chla3um","PO4","NO23","SiO3","din","NPratio","si_din")
interp_var <- function(v) {
  bp <- b[is.finite(get(v)), .(val = mean(get(v))), by = .(bdate = date, depth)]
  x <- merge(m[, .(sample, depth, bdate)], bp, by = "bdate", allow.cartesian = TRUE)
  x[, .(value = if (uniqueN(depth.y) >= 2) approx(depth.y, val, xout = depth.x[1], rule = 2)$y
                else val[which.min(abs(depth.y - depth.x[1]))]),
    by = .(sample, depth.x)][, .(sample, variable = v, value)]
}
envlong <- rbindlist(lapply(bvars, interp_var))
nut <- dcast(envlong, sample ~ variable, value.var = "value")
nut <- merge(meta[, .(sample, depth, date, week, yday)], nut, by = "sample")
nut <- merge(nut, m[, .(sample, bottle_date = bdate, dday)], by = "sample")
fwrite(nut, file.path(cfg$out_dir, "11_nutrients_matched.csv"))

# ---- per-depth medians over the sample-matched values -------------------------
# f20: >20 um chlorophyll fraction per sample
nut[, f20 := chla20um / chla]
desc_depth <- nut[depth %in% cfg$depths,
  .(n_samples = .N, chla = median(chla, na.rm = TRUE), f20 = median(f20, na.rm = TRUE),
    din = median(din, na.rm = TRUE), np = median(NPratio, na.rm = TRUE),
    sio3 = median(SiO3, na.rm = TRUE), no23 = median(NO23, na.rm = TRUE),
    po4 = median(PO4, na.rm = TRUE)), by = depth][order(depth)]
desc_depth[, depth := as.character(depth)]
desc_all <- nut[, .(depth = "all", n_samples = .N, chla = median(chla, na.rm = TRUE),
  f20 = median(f20, na.rm = TRUE), din = median(din, na.rm = TRUE), np = median(NPratio, na.rm = TRUE),
  sio3 = median(SiO3, na.rm = TRUE), no23 = median(NO23, na.rm = TRUE), po4 = median(PO4, na.rm = TRUE))]
desc <- rbind(desc_depth, desc_all)
fwrite(desc, file.path(cfg$out_dir, "11_depth_descriptives.csv"))
cat("\n-- sample-matched medians by depth (11_depth_descriptives.csv) --\n")
print(desc[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

# ---- join chain / composition / buoy -----------------------------------------
chn <- fread("outputs/tangesund/03_chain_metrics.csv")[, .(sample, qc_pass, reliable_chain, mean_chain, frac_cells_chain)]
cmp <- fread("outputs/tangesund/05_composition.csv")[, .(sample, richness, f_chain_diatom = `f_chain diatom`)]
buoy <- if (file.exists("outputs/tangesund/06_env_matched.csv"))
  fread("outputs/tangesund/06_env_matched.csv")[, .(sample, insitu_temp, insitu_sal)] else NULL
dat <- Reduce(function(a,b) merge(a,b,by="sample",all.x=TRUE), c(list(nut, chn, cmp), if(!is.null(buoy)) list(buoy)))
ana <- dat[qc_pass == TRUE & reliable_chain == TRUE]
ana[, day := factor(date)]

# ---- buoy vs bottle cross-check (independent hydrography) ---------------------
if (!is.null(buoy)) {
  cc <- ana[is.finite(insitu_temp) & is.finite(tempDown)]
  cat(sprintf("\nbuoy vs bottle CTD cross-check: temp r=%.2f, salinity r=%.2f (n=%d)\n",
              cor(cc$insitu_temp, cc$tempDown), cor(cc$insitu_sal, cc$salinityDown), nrow(cc)))
}

# ---- resource GAMMs ----------------------------------------------------------
fit <- function(yv, terms) tryCatch(
  bam(reformulate(c(terms, "s(day,bs=\"re\")"), yv), data = ana[is.finite(get(yv))], method = "fREML", select = TRUE),
  error = function(e) { message(yv, ": ", conditionMessage(e)); NULL })
m_si <- fit("mean_chain", c("s(SiO3,k=5)","s(NPratio,k=5)","s(yday,k=6)"))
if (!is.null(m_si)) {
  ss <- summary(m_si); st <- ss$s.table
  cat(sprintf("\nchain ~ SiO3 + N:P + season: dev.expl %.1f%%\n", 100*ss$dev.expl))
  for (tm in c("s(SiO3)","s(NPratio)","s(yday)")) if (tm %in% rownames(st))
    cat(sprintf("  %-12s edf=%.2f F=%.1f p=%.3g\n", tm, st[tm,"edf"], st[tm,"F"], st[tm,"p-value"]))
}

# ---- figures -----------------------------------------------------------------
p1 <- ggplot(ana, aes(SiO3, mean_chain, colour = depth)) +
  geom_point(size = 0.5, alpha = 0.35) + scale_colour_viridis_c(name = "Depth (m)") +
  geom_smooth(method = "gam", formula = y ~ s(x, k=5), colour = "firebrick", se = FALSE) +
  theme_paper() + labs(x = "Silicate SiO3 (umol/l)", y = "Mean chain length",
    title = "Chain length vs silicate")
p2 <- ggplot(ana, aes(NPratio, mean_chain, colour = depth)) +
  geom_point(size = 0.5, alpha = 0.35) + scale_colour_viridis_c(name = "Depth (m)") +
  geom_smooth(method = "gam", formula = y ~ s(x, k=5), colour = "firebrick", se = FALSE) +
  theme_paper() + labs(x = "N:P ratio", y = "Mean chain length", title = "Chain length vs N:P")
p3 <- ggplot(ana[is.finite(chla20um)], aes(chla20um, frac_cells_chain, colour = depth)) +
  geom_point(size = 0.6, alpha = 0.4) + scale_colour_viridis_c(name = "Depth (m)") +
  geom_smooth(method = "gam", formula = y ~ s(x, k=5), colour = "firebrick", fill = "firebrick", se = TRUE, level = 0.68, alpha = 0.15) +
  theme_paper() + labs(x = "Chlorophyll >20 µm (µg L⁻¹)", y = "Chain-diatom share of cells",
    title = "Large-fraction chl vs chain-diatom dominance")
# nutrient depth-time context (bottle grid)
bl <- melt(b[, .(date, depth, SiO3, din, NPratio, chla)], id.vars = c("date","depth"),
           variable.name = "var", value.name = "val")
bl[, var := factor(c(SiO3="Silicate (µM)", din="DIN (µM)", NPratio="Molar N:P", chla="Chl a (µg L⁻¹)")[var],
                   levels = c("Silicate (µM)","DIN (µM)","Molar N:P","Chl a (µg L⁻¹)"))]
# one plot per variable so each keeps its own absolute color scale
bl4 <- bl[is.finite(val)]
# facet strip names the variable (28 removes panel titles); discrete depth axis keeps
# the 1 m tile from extending to 0 m
tile_one <- function(v) ggplot(bl4[var == v], aes(date, factor(depth), fill = val)) +
  geom_tile() + scale_y_discrete(limits = rev) + facet_wrap(~ var) +
  scale_fill_viridis_c(option = "C", name = NULL) +
  theme_paper() + theme(legend.key.height = unit(0.8, "lines"), legend.position = "right") +
  labs(x = NULL, y = if (v == levels(bl4$var)[1]) "Depth (m)" else NULL,
       tag = if (v == levels(bl4$var)[1]) "b" else NULL)
p4 <- wrap_plots(lapply(levels(bl4$var), tile_one), nrow = 1)

ggsave(file.path(cfg$out_dir, "11_chain_vs_silicate.png"), p1, width = 7, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "11_chain_vs_np.png"),        p2, width = 7, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "11_largechl_vs_chaindiatom.png"), p3, width = 7, height = 4.5, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "11_nutrient_depth_time.png"), p4, width = 13, height = 3.6, dpi = FIG_DPI)
ggsave(file.path(cfg$out_dir, "11_nutrient_panel.png"),
       # p1 and p2 are written as standalone figures only
       (((plot_spacer() | (p3 + labs(tag = "a")) | plot_spacer()) + plot_layout(widths = c(1, 2, 1))) /
        p4) &
         theme(plot.tag = element_text(face = "bold", size = 15)),
       width = 13, height = 8, dpi = FIG_DPI)

cat(sprintf("table -> 11_nutrients_matched.csv | panel -> %s/11_nutrient_panel.png (plus the individual 11_*.png panels)\n", cfg$out_dir))
