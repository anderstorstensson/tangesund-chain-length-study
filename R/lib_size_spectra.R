# =============================================================================
# Size-spectrum functions: bounded power-law (PLB) MLE of the size-distribution
# exponent b on unbinned sizes (Edwards et al. 2017; same estimator as
# sizeSpectra::MLE.PLB), plus size diversity and binned NBSS helpers.
#
# PLB pdf on [xmin, xmax]:
#   f(x) = C x^b ,  C = (b+1) / (xmax^(b+1) - xmin^(b+1))   (b != -1)
#                   C = 1 / (log xmax - log xmin)           (b  = -1)
#
# Reads nothing, writes nothing. Deps: base R (stats).
# Use: source("R/lib_size_spectra.R"). Tested by tests/test_size_spectra_synthetic.R.
# =============================================================================

# negative log-likelihood of the PLB exponent b given data x in [xmin, xmax]
negLL_PLB <- function(b, x, xmin, xmax, sumlogx = sum(log(x)), n = length(x)) {
  if (abs(b + 1) < 1e-8) {
    logC <- -log(log(xmax) - log(xmin))
  } else {
    C <- (b + 1) / (xmax^(b + 1) - xmin^(b + 1))   # > 0 even when b < -1
    logC <- log(C)
  }
  -(n * logC + b * sumlogx)
}

# MLE of b. x: sizes; w: weights (e.g. cells per ROI for a per-cell spectrum).
# Returns b, fitted range, effective n (= sum w), negLL and at_bound.
fit_plb <- function(x, w = rep(1, length(x)), xmin = min(x), xmax = max(x), interval = c(-10, 5)) {
  keep <- is.finite(x) & is.finite(w) & x >= xmin & x <= xmax & x > 0
  keep[is.na(keep)] <- FALSE
  x <- x[keep]; w <- w[keep]
  n <- sum(w)
  if (length(x) < 2 || !isTRUE(xmin < xmax) || !isTRUE(n >= 2))
    return(list(b = NA_real_, xmin = xmin, xmax = xmax, n = n, negLL = NA_real_))
  slx <- sum(w * log(x))                                   # weighted sufficient stat
  opt <- optimize(negLL_PLB, interval = interval, x = x, xmin = xmin, xmax = xmax,
                  sumlogx = slx, n = n)
  # at_bound: the optimizer stopped at the edge of the search interval (not converged)
  tol <- 1e-3 * diff(interval)
  at_bound <- isTRUE(opt$minimum <= interval[1] + tol) ||
              isTRUE(opt$minimum >= interval[2] - tol)
  list(b = opt$minimum, xmin = xmin, xmax = xmax, n = n, negLL = opt$objective,
       at_bound = at_bound)
}

# Bootstrap percentile CI for b. Resamples ROIs with their weights (cluster bootstrap
# when w = cells per ROI).
plb_bootstrap <- function(x, w = rep(1, length(x)), xmin = min(x), xmax = max(x),
                          B = 1000, conf = 0.95, interval = c(-10, 5)) {
  keep <- is.finite(x) & is.finite(w) & x >= xmin & x <= xmax & x > 0
  keep[is.na(keep)] <- FALSE
  x <- x[keep]; w <- w[keep]; m <- length(x)
  if (m < 2 || !isTRUE(sum(w) >= 2)) return(c(b = NA_real_, lo = NA_real_, hi = NA_real_))
  bb <- vapply(seq_len(B), function(i) {
    s <- sample.int(m, m, replace = TRUE)
    fit_plb(x[s], w[s], xmin = xmin, xmax = xmax, interval = interval)$b
  }, numeric(1))
  a <- (1 - conf) / 2
  c(b = fit_plb(x, w, xmin, xmax, interval)$b,
    lo = unname(quantile(bb, a, na.rm = TRUE)),
    hi = unname(quantile(bb, 1 - a, na.rm = TRUE)))
}

# Equivalent spherical diameter (um) from biovolume (um^3)
esd_from_biovol <- function(v) (6 * v / pi)^(1 / 3)

# Size diversity = (weighted) variance of log10(size) over the fitted range
# (Quintana et al.). w = cells/ROI for a per-cell spectrum.
size_diversity <- function(x, w = rep(1, length(x)), xmin = min(x), xmax = max(x)) {
  keep <- is.finite(x) & is.finite(w) & x >= xmin & x <= xmax & x > 0
  keep[is.na(keep)] <- FALSE
  x <- x[keep]; w <- w[keep]
  if (length(x) < 2 || !isTRUE(sum(w) >= 2)) return(NA_real_)
  lx <- log10(x); mu <- sum(w * lx) / sum(w)
  sum(w * (lx - mu)^2) / sum(w)
}

# Binned normalized biomass size spectrum (log2 bins), for plotting and linearity
# checks. Returns bin midpoints, NBSS and the R^2 of a log-log OLS line.
nbss_binned <- function(v, bin_base = 2) {
  v <- v[v > 0 & is.finite(v)]
  if (length(v) < 8) return(list(mids = numeric(0), nbss = numeric(0), r2 = NA_real_))
  lo <- floor(log(min(v), bin_base)); hi <- ceiling(log(max(v), bin_base))
  breaks <- bin_base^(lo:hi)
  bin <- cut(v, breaks, include.lowest = TRUE)
  bw  <- diff(breaks)
  tot <- tapply(v, bin, sum)                       # total biovolume per bin
  tot[is.na(tot)] <- 0
  mids <- sqrt(head(breaks, -1) * tail(breaks, -1))
  nbss <- as.numeric(tot) / bw                     # normalize by bin width
  keep <- nbss > 0
  r2 <- if (sum(keep) >= 3)
          summary(lm(log10(nbss[keep]) ~ log10(mids[keep])))$r.squared else NA_real_
  list(mids = mids[keep], nbss = nbss[keep], r2 = r2)
}

# Random draw from a PLB distribution (for the synthetic recovery test)
rplb <- function(n, b, xmin, xmax) {
  u <- runif(n)
  if (abs(b + 1) < 1e-8) {
    xmin * (xmax / xmin)^u
  } else {
    (u * (xmax^(b + 1) - xmin^(b + 1)) + xmin^(b + 1))^(1 / (b + 1))
  }
}
