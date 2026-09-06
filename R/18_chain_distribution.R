# =============================================================================
# Chain-length distribution by depth: fraction of chains with > 1 cell, median, mean
# and 90th percentile of the per-sample distribution, with day-cluster bootstrap SE,
# and their 3 m to 16 m change.
#
# Reads:  outputs/tangesund/03_chain_metrics.csv
# Writes: outputs/tangesund/18_chain_distribution_profile.csv, 18_chain_distribution_change.csv
# Deps:   data.table
# Run:    Rscript R/18_chain_distribution.R
# =============================================================================
Sys.setlocale("LC_TIME", "C")
suppressPackageStartupMessages({ library(data.table) })

cfg <- list(
  chain_csv = "outputs/tangesund/03_chain_metrics.csv",
  out_dir   = "outputs/tangesund",
  depths    = c(3, 6, 8, 11, 13, 16),
  nboot     = 2000,
  seed      = 1
)

d <- fread(cfg$chain_csv)
d[, day := factor(as.Date(date))]
ana <- d[qc_pass == TRUE & reliable_chain == TRUE]
cat(sprintf("chain-trait-thresholded samples: %d over %d days\n", nrow(ana), uniqueN(ana$day)))

# day-cluster bootstrap SE of the sample-weighted mean (as in 03 and 16)
boot_se_day <- function(value, day, B = cfg$nboot, seed = cfg$seed) {
  days <- unique(day); nd <- length(days)
  idx_by_day <- split(seq_along(value), day)
  set.seed(seed)
  sd(vapply(seq_len(B), function(b) {
    ds <- sample(days, nd, replace = TRUE)
    mean(value[unlist(idx_by_day[as.character(ds)], use.names = FALSE)])
  }, numeric(1)))
}

prof <- ana[, .(n = .N, n_day = uniqueN(day),
                frac_multi = mean(frac_multi), se_frac  = boot_se_day(frac_multi, day),
                med        = mean(med_chain),  se_med   = boot_se_day(med_chain, day),
                mean_chain = mean(mean_chain), se_mean  = boot_se_day(mean_chain, day),
                p90        = mean(p90_chain),  se_p90   = boot_se_day(p90_chain, day)),
            by = depth][order(depth)]
fwrite(prof, file.path(cfg$out_dir, "18_chain_distribution_profile.csv"))

cat("\n== chain-length distribution by depth (per-sample statistics, day-bootstrap SE) ==\n")
print(prof[, .(depth, n,
               multi = sprintf("%.3f ± %.3f", frac_multi, se_frac),
               median = sprintf("%.2f ± %.2f", med, se_med),
               mean = sprintf("%.2f ± %.2f", mean_chain, se_mean),
               p90 = sprintf("%.2f ± %.2f", p90, se_p90))])

pc <- function(a, b) 100 * (b - a) / a
s <- prof[depth == min(cfg$depths)]; dp <- prof[depth == max(cfg$depths)]
cat("\n-- surface(3 m) to deep(16 m) change --\n")
chg <- data.table(
  statistic = c("Fraction of chains with >1 cell", "Median chain length (cells)",
                "Mean chain length (cells)", "90th percentile (cells)"),
  surface = c(s$frac_multi, s$med, s$mean_chain, s$p90),
  deep    = c(dp$frac_multi, dp$med, dp$mean_chain, dp$p90))
chg[, pct_change := pc(surface, deep)]
print(chg[, .(statistic, surface = round(surface, 3), deep = round(deep, 3),
              pct_change = sprintf("%+.0f%%", pct_change))])
fwrite(chg, file.path(cfg$out_dir, "18_chain_distribution_change.csv"))

cat(sprintf("tables -> %s/18_chain_distribution_{profile,change}.csv\n", cfg$out_dir))
