skip_on_cran()

wls_beta <- function(x, p, w, tol = 1e-8) {
  keep <- is.finite(x) & is.finite(p) & is.finite(w)
  if (!any(keep)) return(NA_real_)
  x <- x[keep]
  p <- p[keep]
  w <- w[keep]
  if (!any(w > 0)) return(NA_real_)
  wm <- sum(w)
  if (!is.finite(wm) || wm <= 0) return(NA_real_)
  pm <- sum(w * p) / wm
  xm <- sum(w * x) / wm
  num <- sum(w * (p - pm) * (x - xm))
  den <- sum(w * (p - pm)^2)
  if (!is.finite(den) || den <= tol) return(NA_real_)
  num / den
}

simulate_modulated <- function(pk, base_hz, beta_hz, tr) {
  pk <- as.numeric(pk)
  Tn <- length(pk)
  V <- length(base_hz)
  base_omega <- 2 * pi * base_hz
  beta <- 2 * pi * beta_hz
  phi0 <- runif(V, -pi, pi)
  omega_mat <- matrix(base_omega, Tn, V, byrow = TRUE) + pk %o% beta
  phi <- sweep(tr * apply(omega_mat, 2, cumsum), 2, phi0, `+`)
  list(U = sin(phi), phi = phi, omega = omega_mat, beta = beta)
}

test_that("ipppi_pwi recovers modulation slopes and anti-symmetry", {
  set.seed(1001)
  tr <- 0.6
  fs <- 1 / tr
  Tn <- 480
  t <- seq(0, by = tr, length.out = Tn)
  pk <- sin(2 * pi * 0.01 * t) + 0.2 * rnorm(Tn)
  pk <- scale(pk, center = TRUE, scale = TRUE)[, 1]
  base_hz <- seq(0.045, 0.065, length.out = 5)
  beta_hz <- c(0.004, -0.002, 0.000, 0.003, -0.0015)
  sim <- simulate_modulated(pk, base_hz, beta_hz, tr)
  band <- c(0.03, 0.09)
  fi <- ipppi_instfreq(sim$U, fs = fs, band = band)
  out <- ipppi_pwi(sim$U, pk, fs = fs, band = band, weight = "unit", robust = FALSE)
  manual <- vapply(seq_len(ncol(sim$U)), function(j) wls_beta(fi$omega[, j], pk, rep(1, length(pk))), numeric(1))
  expect_equal(out$beta, manual, tolerance = 1e-7)
  expect_lt(max(abs(manual - sim$beta)), 0.01)
  expect_equal(out$mat, outer(out$beta, out$beta, `-`), tolerance = 1e-7)
  expect_equal(diag(out$mat), rep(0, length(base_hz)))
  expect_equal(out$mat, -t(out$mat), tolerance = 1e-7)
})

test_that("ipppi_pwi is invariant to constant delays", {
  set.seed(203)
  tr <- 0.7
  fs <- 1 / tr
  Tn <- 420
  t <- seq(0, by = tr, length.out = Tn)
  pk <- sin(2 * pi * 0.008 * t)
  pk <- scale(pk, center = TRUE, scale = TRUE)[, 1]
  base_hz <- rep(0.055, 3)
  beta_hz <- rep(0.0025, 3)
  sim <- simulate_modulated(pk, base_hz, beta_hz, tr)
  delay <- pi / 3
  U <- sim$U
  U[, 2] <- sin(sim$phi[, 2] + delay)
  U[, 3] <- sin(sim$phi[, 3] - delay / 2)
  out <- ipppi_pwi(U, pk, fs = fs, band = c(0.03, 0.08), weight = "unit", robust = FALSE)
  expect_lt(max(abs(out$mat)), 1e-3)
})

test_that("ipppi_pwi weighted slopes match explicit WLS", {
  set.seed(88)
  tr <- 0.65
  fs <- 1 / tr
  Tn <- 360
  t <- seq(0, by = tr, length.out = Tn)
  pk <- cos(2 * pi * 0.015 * t) + 0.1 * rnorm(Tn)
  pk <- scale(pk, center = TRUE, scale = TRUE)[, 1]
  base_hz <- c(0.05, 0.06, 0.07, 0.055)
  beta_hz <- c(0.003, -0.002, 0.0015, 0)
  phase_mod <- runif(length(base_hz), 0, 2 * pi)
  amp <- 1 + 0.2 * outer(t, phase_mod, function(tt, ph) sin(2 * pi * 0.01 * tt + ph))
  amp[amp < 0.1] <- 0.1
  sim <- simulate_modulated(pk, base_hz, beta_hz, tr)
  U <- amp * sin(sim$phi)
  res <- ipppi_instfreq(U, fs = fs, band = c(0.03, 0.09))
  out <- ipppi_pwi(U, pk, fs = fs, band = c(0.03, 0.09), weight = "amplitude", robust = FALSE)
  manual <- vapply(seq_len(ncol(U)), function(j) {
    wls_beta(res$omega[, j], pk, res$amp[, j])
  }, numeric(1))
  expect_equal(out$beta, manual, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("Robust PWI reduces bias under phase outliers", {
  set.seed(512)
  tr <- 0.6
  fs <- 1 / tr
  Tn <- 500
  t <- seq(0, by = tr, length.out = Tn)
  pk <- sin(2 * pi * 0.01 * t)
  pk <- scale(pk, center = TRUE, scale = TRUE)[, 1]
  base_hz <- c(0.05, 0.058)
  beta_hz <- c(0.004, -0.003)
  sim <- simulate_modulated(pk, base_hz, beta_hz, tr)
  U <- sim$U
  spike_idx <- sample(60:(Tn - 60), 15)
  U[spike_idx, 1] <- sin(sim$phi[spike_idx, 1] + runif(length(spike_idx), -pi, pi))
  ols <- ipppi_pwi(U, pk, fs = fs, band = c(0.03, 0.09), weight = "unit", robust = FALSE)
  rob <- ipppi_pwi(U, pk, fs = fs, band = c(0.03, 0.09), weight = "unit", robust = TRUE)
  err_ols <- abs(ols$beta - sim$beta)
  err_rob <- abs(rob$beta - sim$beta)
  expect_lt(sum(err_rob, na.rm = TRUE), sum(err_ols, na.rm = TRUE))
})

if (exists("ipppi_iphs_phase", mode = "function")) {
  test_that("Hilbert end-to-end PWI recovers known modulation", {
    set.seed(77)
    fs <- 1.2
    tr <- 1 / fs
    Tn <- 420
    t <- seq(0, by = tr, length.out = Tn)
    pk <- sin(2 * pi * 0.012 * t)
    pk <- scale(pk, center = TRUE, scale = TRUE)[, 1]
    base_hz <- c(0.045, 0.055, 0.062)
    beta_hz <- c(0.0025, -0.0015, 0.0)
    sim <- simulate_modulated(pk, base_hz, beta_hz, tr)
    band <- c(0.03, 0.08)
    out <- ipppi_pwi(sim$U, pk, fs = fs, band = band, weight = "unit", robust = FALSE)
    expect_gt(cor(out$beta, sim$beta), 0.98)
  })
}
