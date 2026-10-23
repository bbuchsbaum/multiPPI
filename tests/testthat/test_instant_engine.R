test_that("IPHS recovers known phase coupling", {
  fs <- 1.4
  Tn <- 400
  t <- seq(0, length.out = Tn, by = 1 / fs)
  phase <- pi / 3
  x <- sin(2 * pi * 0.05 * t)
  y <- sin(2 * pi * 0.05 * t + phase)
  U <- cbind(x, y)
  pk <- matrix(1, nrow = Tn, ncol = 1)
  corr_mat <- ipppi_iphs(U, pk, fs = fs, band = c(0.04, 0.06), tau_half = 12)[[1]]
  expect_true(abs(corr_mat[1, 2] - cos(phase)) < 0.1)
})

test_that("mppi_instant_fit produces coherent structures", {
  set.seed(2)
  Tn <- 60
  V <- 6
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(1, rnorm(Tn))
  fit <- mppi_instant_fit(Y, X, psych_idx = 2L, runs = rep(1L, Tn),
                          instant = list(method = "ewm", tau_half = 5, offsets = 0:2))
  expect_s3_class(fit, "mppi_fit")
  expect_equal(fit$engine, "instant")
  expect_length(fit$Delta, 1L)
  expect_true(all(c("0", "1", "2") %in% names(fit$lagged[[1L]])))
  M0 <- mppi_get_M(fit, 1L)
  M1 <- mppi_get_M_lag(fit, 1L, lag = 1L)
  expect_equal(dim(M0), c(V, V))
  expect_equal(dim(M1), c(V, V))
  expect_equal(fit$lag_blocklens, rle(rep(1L, Tn))$lengths)
  slopes_mean <- mppi_seed_slopes(fit, seeds = 1:2, k = 1L, collapse = "mean", lag = 1L)
  expect_equal(nrow(slopes_mean), 1L)
})
