test_that("PWI recovers directed modulation of instantaneous frequency", {
  skip_on_cran()
  set.seed(81)
  tr <- 0.5
  fs <- 1 / tr
  Tn <- 480
  t <- seq(0, by = tr, length.out = Tn)
  pk_raw <- sin(2 * pi * 0.01 * t)
  pk <- as.numeric(scale(pk_raw, center = TRUE, scale = TRUE))
  base_hz <- c(0.07, 0.07, 0.07)
  beta_hz <- c(0.02, 0.00, -0.01)
  phase0 <- runif(length(base_hz), 0, 2 * pi)
  signal <- vapply(seq_along(base_hz), function(i) {
    freq <- base_hz[i] + beta_hz[i] * pk
    sin(2 * pi * freq * t + phase0[i]) + 0.05 * rnorm(Tn)
  }, numeric(Tn))
  out <- ipppi_pwi(signal, pk, fs = fs, band = c(0.03, 0.1), weight = "amplitude")
  expect_equal(dim(out$mat), c(3, 3))
  expect_equal(out$mat, -t(out$mat))
  beta_true <- 2 * pi * beta_hz
  expect_gt(stats::cor(out$beta, beta_true), 0.8)
})

test_that("Instantaneous frequency resets at run boundaries", {
  skip_on_cran()
  tr <- 0.7
  fs <- 1 / tr
  n_run <- 60
  t1 <- seq(0, by = tr, length.out = n_run)
  t2 <- seq(0, by = tr, length.out = n_run)
  x1 <- sin(2 * pi * 0.05 * t1)
  x2 <- sin(2 * pi * 0.05 * t2 + pi / 3)
  U <- rbind(cbind(x1, x1^2), cbind(x2, x2^2))
  runs <- rep(1:2, each = n_run)
  fi <- ipppi_instfreq(U, fs = fs, band = c(0.02, 0.09), runs = runs)
  expect_true(all(is.na(fi$omega[c(1, n_run + 1), ])))
  expect_false(any(is.infinite(fi$omega), na.rm = TRUE))
})

test_that("PWI handles disjoint event islands without drift", {
  skip_on_cran()
  set.seed(94)
  tr <- 0.8
  fs <- 1 / tr
  blocks <- list(1:40, 141:180)
  Tn <- 220
  t <- seq(0, by = tr, length.out = Tn)
  pk <- rep(0, Tn)
  pk[blocks[[1]]] <- 1
  pk[blocks[[2]]] <- -1
  pk <- scale(pk, center = TRUE, scale = TRUE)[, 1]
  base_hz <- c(0.06, 0.06, 0.06)
  beta_hz <- c(0.015, -0.005, 0.00)
  phase0 <- runif(length(base_hz), 0, 2 * pi)
  signal <- vapply(seq_along(base_hz), function(i) {
    freq <- base_hz[i] + beta_hz[i] * pk
    sin(2 * pi * freq * t + phase0[i]) + 0.03 * rnorm(Tn)
  }, numeric(Tn))
  out <- ipppi_pwi(signal, pk, fs = fs, band = c(0.02, 0.09), weight = "amplitude")
  beta_true <- 2 * pi * beta_hz
  expect_gt(stats::cor(out$beta, beta_true), 0.7)
  expect_equal(out$mat, -t(out$mat))
})

test_that("PWI skips NA psych frames", {
  skip_on_cran()
  set.seed(123)
  tr <- 0.6
  fs <- 1 / tr
  Tn <- 360
  t <- seq(0, by = tr, length.out = Tn)
  pk <- sin(2 * pi * 0.02 * t)
  pk <- scale(pk, center = TRUE, scale = TRUE)[, 1]
  na_idx <- sample(seq_len(Tn), 10)
  pk[na_idx] <- NA_real_
  base_hz <- rep(0.05, 4)
  beta_hz <- c(0.03, -0.02, 0.0, 0.015)
  phase0 <- runif(length(base_hz), 0, 2 * pi)
  signal <- vapply(seq_along(base_hz), function(i) {
    freq <- base_hz[i] + beta_hz[i] * ifelse(is.na(pk), 0, pk)
    sin(2 * pi * freq * t + phase0[i]) + 0.02 * rnorm(Tn)
  }, numeric(Tn))
  out <- ipppi_pwi(signal, pk, fs = fs, band = c(0.02, 0.08), weight = "unit")
  beta_true <- 2 * pi * beta_hz
  expect_gt(stats::cor(out$beta, beta_true, use = "complete.obs"), 0.6)
  expect_false(any(is.nan(out$beta)))
})
