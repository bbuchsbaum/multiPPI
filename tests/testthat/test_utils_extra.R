library(testthat)

test_that(".mppi_handle_na drops incomplete rows and still errors when requested", {
  Y <- matrix(seq_len(12), nrow = 4)
  X <- matrix(seq(2, 24, by = 2), nrow = 4)
  runs <- c(1L, 1L, 2L, 2L)
  Y[2, 1] <- NA_real_
  X[4, 2] <- NA_real_

  res <- multiPPI:::`.mppi_handle_na`(Y, X, runs = runs, na_action = "omit_tr")
  expect_equal(res$dropped, 2L)
  expect_equal(res$dropped_idx, c(2L, 4L))
  expect_equal(res$Y, Y[c(1L, 3L), ], ignore_attr = TRUE)
  expect_equal(res$X, X[c(1L, 3L), ], ignore_attr = TRUE)
  expect_equal(res$runs, runs[c(1L, 3L)])

  expect_error(
    multiPPI:::`.mppi_handle_na`(Y, X, runs = runs, na_action = "error"),
    "mPPI inputs contain NA values"
  )
})


test_that(".mppi_prepare_design_matrix adds intercepts only when needed", {
  X_no_int <- matrix(c(0, 1, 0, 1, 0, 1), nrow = 3, byrow = FALSE)
  colnames(X_no_int) <- c("taskA", "taskB")
  prep_no_int <- multiPPI:::`.mppi_prepare_design_matrix`(X_no_int, psych_idx = c(1L, 2L))
  expect_true(prep_no_int$intercept_added)
  expect_equal(colnames(prep_no_int$X)[1], "(Intercept)")
  expect_equal(prep_no_int$psych_idx, c(2L, 3L))

  X_with_int <- cbind(`(Intercept)` = 1, X_no_int)
  prep_with_int <- multiPPI:::`.mppi_prepare_design_matrix`(X_with_int, psych_idx = c(2L, 3L))
  expect_false(prep_with_int$intercept_added)
  expect_identical(prep_with_int$X, X_with_int)
  expect_equal(prep_with_int$psych_idx, c(2L, 3L))
})


test_that(".mppi_project_basis chunked path matches BLAS product", {
  set.seed(123)
  R <- matrix(rnorm(60), nrow = 10)
  V <- matrix(rnorm(12), nrow = 6)
  direct <- R %*% V
  chunked <- multiPPI:::`.mppi_project_basis`(R, V, backend = "chunked", chunk_cols = 3L)
  expect_equal(chunked, direct, tolerance = 1e-10)
})


test_that(".mppi_crossprod_lagged handles run-aware lags", {
  set.seed(99)
  R <- matrix(rnorm(40), nrow = 8)
  pk <- runif(8)
  blocklens <- c(5L, 3L)

  pos_res <- multiPPI:::`.mppi_crossprod_lagged`(R, pk, lag = 1L, blocklens = blocklens)
  start <- 1L
  accum <- matrix(0, ncol(R), ncol(R))
  denom <- 0
  for (L in blocklens) {
    if (L > 1L) {
      rng <- seq.int(start, start + L - 1L - 1L)
      R0 <- R[rng, , drop = FALSE]
      R1 <- R[rng + 1L, , drop = FALSE]
      wk <- pk[rng]
      accum <- accum + crossprod(R0, wk * R1)
      denom <- denom + sum(wk^2)
    }
    start <- start + L
  }
  expect_equal(pos_res, accum / denom, tolerance = 1e-10)

  neg_res <- multiPPI:::`.mppi_crossprod_lagged`(R, pk, lag = -1L, blocklens = blocklens)
  start <- 1L
  accum <- matrix(0, ncol(R), ncol(R))
  denom <- 0
  for (L in blocklens) {
    shift <- 1L
    if (L > shift) {
      rng <- seq.int(start + shift, start + L - 1L)
      R0 <- R[rng, , drop = FALSE]
      R1 <- R[rng - shift, , drop = FALSE]
      wk <- pk[rng - shift]
      accum <- accum + crossprod(R0, wk * R1)
      denom <- denom + sum(wk^2)
    }
    start <- start + L
  }
  expect_equal(neg_res, accum / denom, tolerance = 1e-10)
})


test_that("mppi_fuse_time_beta returns precision-weighted average", {
  D_time <- matrix(c(2, 1, 1, 3), nrow = 2)
  D_beta <- matrix(c(4, 0, 0, 2), nrow = 2)
  fused <- mppi_fuse_time_beta(D_time, D_beta, var_time = 0.5, var_beta = 2)
  w_time <- 1 / 0.5
  w_beta <- 1 / 2
  expected <- (w_time * D_time + w_beta * D_beta) / (w_time + w_beta)
  expect_equal(fused, expected)
})
