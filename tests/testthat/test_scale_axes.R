test_that("scaled interaction matrices are scale-invariant", {
  set.seed(101)
  Tn <- 60L
  V <- 4L
  X <- cbind(1, rnorm(Tn))
  colnames(X) <- c("(Intercept)", "task")
  Y <- matrix(rnorm(Tn * V), Tn, V)

  fit1 <- mppi_fit(Y, X, psych_idx = 2L, lags = c(-1L, 0L, 1L))
  fit2 <- mppi_fit(2 * Y, X, psych_idx = 2L, lags = c(-1L, 0L, 1L))

  M_norm1 <- mppi_get_M_scaled(fit1, 1L, mode = "normalized")
  M_norm2 <- mppi_get_M_scaled(fit2, 1L, mode = "normalized")
  expect_equal(M_norm1, M_norm2, tolerance = 1e-6)

  M_amp1 <- mppi_get_M_scaled(fit1, 1L, mode = "amplitude")
  M_amp2 <- mppi_get_M_scaled(fit2, 1L, mode = "amplitude")
  expect_equal(norm(M_amp2, type = "F") / norm(M_amp1, type = "F"), 2, tolerance = 1e-6)
})

test_that("gain and routing axes return bounded summaries", {
  set.seed(202)
  Tn <- 80L
  V <- 5L
  X <- cbind(1, rnorm(Tn), rnorm(Tn))
  colnames(X) <- c("(Intercept)", "cond1", "cond2")
  Y <- matrix(rnorm(Tn * V), Tn, V)

  fit <- mppi_fit(Y, X, psych_idx = 2:3, lags = -1:1)
  axes <- mppi_axes(fit, lags = -1:1)

  expect_equal(nrow(axes), length(fit$names))
  expect_true(all(!is.na(axes$gain)))
  expect_true(all(abs(axes$gain) <= 1 + 1e-8))
  expect_true(all(abs(axes$routing) <= 1 + 1e-8 | is.na(axes$routing)))
})

test_that("communication modes return projected matrices for fits", {
  set.seed(303)
  Tn <- 50L
  V <- 6L
  X <- cbind(1, rnorm(Tn))
  colnames(X) <- c("(Intercept)", "cond")
  Y <- matrix(rnorm(Tn * V), Tn, V)

  fit <- mppi_fit(Y, X, psych_idx = 2L)
  modes <- mppi_modes(fit, r = 3L)

  expect_equal(ncol(modes$vectors), 3L)
  expect_equal(length(modes$projected), length(fit$names))
  expect_true(all(vapply(modes$projected, function(M) nrow(M) == 3, logical(1))))
})
