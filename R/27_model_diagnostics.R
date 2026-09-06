# =============================================================================
# Sampling-design diagnostics: depth x hour coverage, and concurvity of the
# hydrography GAMM (model A of 10) and the joint depth + salinity + temperature model.
#
# Reads:  outputs/tangesund/03_chain_metrics.csv, 06_env_matched.csv
# Writes: outputs/tangesund/27_depth_hour_coverage.csv, 27_concurvity_full.csv,
#         27_concurvity_pairwise.csv, 27_concurvity_joint.csv
# Deps:   data.table, mgcv
# Run:    Rscript R/27_model_diagnostics.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table); library(mgcv) })
cfg <- list(out_dir = "outputs/tangesund", depths = c(3,6,8,11,13,16),
            lat = 58.075, lon = 11.493)

solar_elev <- function(dt, lat=cfg$lat, lon=cfg$lon) {
  doy  <- as.integer(format(dt, "%j")); hutc <- as.integer(format(dt,"%H")) + as.integer(format(dt,"%M"))/60
  decl <- 23.45*sin(2*pi*(284+doy)/365)*pi/180; H <- (15*((hutc+lon/15)-12))*pi/180
  asin(sin(lat*pi/180)*sin(decl)+cos(lat*pi/180)*cos(decl)*cos(H))*180/pi
}

## ---------- (1) Depth x hour coverage (diel design) --------------------------
dat <- fread(file.path(cfg$out_dir, "03_chain_metrics.csv"))
ts  <- regmatches(dat$sample, regexpr("D[0-9]{8}T[0-9]{6}", dat$sample))
dat[, datetime := as.POSIXct(ts, format = "D%Y%m%dT%H%M%S", tz = "UTC")]
dat[, hour := as.integer(format(datetime, "%H"))]
ana <- dat[qc_pass == TRUE & reliable_chain == TRUE]
ana[, `:=`(depth_f = factor(depth, levels = cfg$depths),
           daynight = fifelse(solar_elev(datetime) > 0, "day", "night"))]

cat("=========== (1) DEPTH x HOUR COVERAGE (diel design) ===========\n")
cat(sprintf("diel analysis samples: %d\n", nrow(ana)))
tab <- table(depth = ana$depth_f, hour = ana$hour)
cat(sprintf("cor(depth, hour) = %.4f\n", cor(ana$depth, ana$hour)))
sup <- suppressWarnings(chisq.test(tab))
cramV <- sqrt(as.numeric(sup$statistic) / (sum(tab) * (min(dim(tab)) - 1)))
cat(sprintf("chi-square(depth x hour): X2=%.0f df=%d p=%.3g ; Cramer's V=%.3f\n",
            sup$statistic, sup$parameter, sup$p.value, cramV))

covd <- ana[, .(n = .N, hours_covered = uniqueN(hour),
                n_day = sum(daynight=="day"), n_night = sum(daynight=="night")),
            by = depth_f][order(depth_f)]
covd[, cor_depth_hour := round(cor(ana$depth, ana$hour), 4)]
covd[, cramers_v := round(cramV, 3)]
print(covd)
fwrite(covd, file.path(cfg$out_dir, "27_depth_hour_coverage.csv"))

## ---------- (2) Concurvity of the mechanism GAMM ----------------------------
cat("\n=========== (2) CONCURVITY: chain ~ hydrography GAMM (10 model A) ===========\n")
env <- fread(file.path(cfg$out_dir, "06_env_matched.csv"))
chn <- fread(file.path(cfg$out_dir, "03_chain_metrics.csv"))[, .(sample, qc_pass, reliable_chain, mean_chain)]
A   <- merge(env, chn, by="sample", all.x=TRUE)[qc_pass==TRUE & reliable_chain==TRUE &
             is.finite(insitu_temp) & is.finite(mld)]
A[, day := factor(date)]
cat(sprintf("mechanism samples: %d\n", nrow(A)))
vars <- c("depth","insitu_sal","dsigma","mld","insitu_temp","yday")
cat("predictor correlations:\n"); print(round(cor(A[, ..vars]), 2))

mA <- bam(mean_chain ~ s(insitu_temp, k=6) + s(insitu_sal, k=6) + s(dsigma, k=6) +
            s(mld, k=6) + s(yday, k=8) + ti(insitu_sal, yday, k=c(5,6)) + s(day, bs="re"), data = A,
          method = "fREML")   # same as mA in 10
cc  <- concurvity(mA, full = TRUE)
ccw <- concurvity(mA, full = FALSE)$worst
cat("\nconcurvity(full=TRUE) [0=none .. 1=total]:\n"); print(round(cc, 3))
cat("\nconcurvity(full=FALSE)$worst (pairwise):\n"); print(round(ccw, 3))
fwrite(as.data.table(round(cc,3), keep.rownames="measure"),
       file.path(cfg$out_dir, "27_concurvity_full.csv"))
fwrite(as.data.table(round(ccw,3), keep.rownames="term"),
       file.path(cfg$out_dir, "27_concurvity_pairwise.csv"))
obs <- cc["observed", grep("^s\\(|^ti\\(", colnames(cc))]
obs <- obs[names(obs) != "s(day)"]
cat(sprintf("\nobserved concurvity of the fixed-effect smooths: %.2f-%.2f (%s highest)\n",
            min(obs), max(obs), names(obs)[which.max(obs)]))

# joint model of 10 (depth + salinity + temperature)
mJ <- bam(mean_chain ~ s(depth, k=5) + s(insitu_sal, k=6) + s(insitu_temp, k=6) + s(yday, k=8) +
            ti(insitu_sal, yday, k=c(5,6)) + s(day, bs="re"), data = A, method = "fREML")
ccJ <- concurvity(mJ, full = TRUE); ccJw <- concurvity(mJ, full = FALSE)$worst
cat("\njoint model concurvity(full=TRUE):\n"); print(round(ccJ, 3))
cat(sprintf("pairwise worst, s(depth) vs s(insitu_sal): %.3f\n", ccJw["s(depth)", "s(insitu_sal)"]))
fwrite(as.data.table(round(ccJ,3), keep.rownames="measure"),
       file.path(cfg$out_dir, "27_concurvity_joint.csv"))
cat("\nDONE. Tables -> outputs/tangesund/27_*.csv\n")
