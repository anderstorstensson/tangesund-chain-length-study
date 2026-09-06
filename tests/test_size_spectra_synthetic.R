# =============================================================================
# Synthetic recovery test for the PLB size-spectrum estimator in R/lib_size_spectra.R.
#
# Needs no data. Run: Rscript tests/test_size_spectra_synthetic.R
# =============================================================================
source("R/lib_size_spectra.R")
set.seed(42)

ok <- TRUE
report <- function(name, pass, detail = "") {
  cat(sprintf("[%s] %s%s\n", if (pass) "PASS" else "FAIL", name,
              if (nzchar(detail)) paste0("  (", detail, ")") else ""))
  ok <<- ok && pass
}

# IFCB range: ESD 5 to 150 um (biovolume ~65 to 1.8e6 um^3)
xmin <- (pi/6) * 5^3
xmax <- (pi/6) * 150^3

# 1) point recovery
for (b_true in c(-1.5, -2.0, -2.5)) {
  x <- rplb(2e5, b = b_true, xmin = xmin, xmax = xmax)
  bhat <- fit_plb(x, xmin = xmin, xmax = xmax)$b
  report(sprintf("point recovery b=%.1f", b_true), abs(bhat - b_true) < 0.03,
         sprintf("bhat=%.3f", bhat))
}

# 2) bootstrap CI coverage at n = 300
b_true <- -2.0
covered <- 0; reps <- 100
for (i in seq_len(reps)) {
  x <- rplb(300, b = b_true, xmin = xmin, xmax = xmax)
  ci <- plb_bootstrap(x, xmin = xmin, xmax = xmax, B = 200)
  if (!is.na(ci["lo"]) && ci["lo"] <= b_true && b_true <= ci["hi"]) covered <- covered + 1
}
cov <- covered / reps
report("bootstrap 95% CI coverage (n=300)", cov >= 0.90, sprintf("coverage=%.2f", cov))

# 3) small-n estimate
x <- rplb(300, b = -2.0, xmin, xmax)
report("small-n estimate near truth", abs(fit_plb(x, xmin = xmin, xmax = xmax)$b + 2) < 0.25,
       sprintf("bhat=%.3f", fit_plb(x, xmin = xmin, xmax = xmax)$b))

# 3b) weighted fit equals the fit on the expanded sample
set.seed(7)
xu <- rplb(500, b = -2.0, xmin, xmax)
wu <- sample(1:5, length(xu), replace = TRUE)           # integer weights = cell counts
b_w <- fit_plb(xu, w = wu, xmin = xmin, xmax = xmax)$b
b_e <- fit_plb(rep(xu, wu), xmin = xmin, xmax = xmax)$b  # explicit expansion
report("weighted fit == expanded fit", abs(b_w - b_e) < 1e-6, sprintf("w=%.5f e=%.5f", b_w, b_e))
sd_w <- size_diversity(xu, w = wu, xmin = xmin, xmax = xmax)
sd_e <- size_diversity(rep(xu, wu), xmin = xmin, xmax = xmax)
report("weighted size_diversity == expanded", abs(sd_w - sd_e) < 1e-9)

# 4) size_diversity and nbss helpers
sd_val <- size_diversity(rplb(1e4, -2, xmin, xmax), xmin = xmin, xmax = xmax)
report("size_diversity finite > 0", is.finite(sd_val) && sd_val > 0, sprintf("var(log10)=%.3f", sd_val))
nb <- nbss_binned(rplb(1e4, -2, xmin, xmax))
report("nbss linearity R2 high for clean PLB", !is.na(nb$r2) && nb$r2 > 0.95, sprintf("R2=%.3f", nb$r2))

cat(if (ok) "\nALL TESTS PASSED\n" else "\nSOME TESTS FAILED\n")
quit(status = if (ok) 0 else 1)
