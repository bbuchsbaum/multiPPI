test_that("fmriAR prewhitening paths align", {
  skip_if_not_installed("fmriAR")

  set.seed(2024)
  Tn <- 48L
  V <- 5L
  runs <- rep(1:3, each = 16L)
  t <- seq_len(Tn)
  X <- cbind(Intercept = 1, task = sin(t / 3), drift = cos(t / 7))
  Y <- matrix(rnorm(Tn * V), Tn, V)

  fit_helper <- mppi_fit_whitened(Y, X, runs = runs, psych_idx = 2L,
                                  scale = "normalized")
  expect_equal(fit_helper$runs, runs)
  expect_s3_class(fit_helper, "mppi_fit")
  expect_equal(dim(fit_helper$R), c(Tn, V))
  expect_false(all(is.na(fit_helper$R)))
})
