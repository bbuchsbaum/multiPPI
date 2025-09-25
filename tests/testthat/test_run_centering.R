test_that("run centering removes run-level means", {
  set.seed(303)
  Tn <- 40L
  V <- 4L
  runs <- rep(1:2, each = Tn / 2L)
  psych <- sin(seq_len(Tn) / 6) + rep(c(-0.5, 0.5), each = Tn / 2L)
  X <- cbind(1, psych)
  Y <- matrix(rnorm(Tn * V), Tn, V)

  fit_none <- mppi_fit(Y, X, psych_idx = 2L, runs = runs,
                       center_by = "none", scale = "normalized")
  fit_run <- mppi_fit(Y, X, psych_idx = 2L, runs = runs,
                      center_by = "run", scale = "normalized")

  run_means <- tapply(fit_run$pk[[1]], runs, mean)
  expect_true(all(abs(run_means) < 1e-10))
  run_means_none <- tapply(fit_none$pk[[1]], runs, mean)
  expect_true(any(abs(run_means_none) > 1e-6))

  blk <- rep.int(as.integer(Tn / 2L), 2L)
  fit_lag <- mppi_fit(Y, X, psych_idx = 2L, runs = NULL,
                      center_by = "none", scale = "normalized",
                      lags = c(0L, 1L), lag_blocklens = blk)
  expect_equal(fit_lag$lag_blocklens, blk)
})

test_that("mppi_fit_from_fmrireg requires fmrireg package", {
  if (requireNamespace("fmrireg", quietly = TRUE)) {
    skip("fmrireg available; dependency error path not triggered")
  }
  expect_error(mppi_fit_from_fmrireg(NULL, NULL, psych = 1L), "fmrireg not available")
})
