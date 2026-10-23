test_that("instant neural deconvolution is run local", {
  set.seed(17); n <- 120; runs <- rep(1:2, each = 60)
  X <- cbind(run1 = runs == 1, run2 = runs == 2, task = rnorm(n)); Y <- matrix(rnorm(n * 3), n)
  i2 <- runs == 2
  Y2 <- Y; Y2[i2, ] <- Y[i2, ] + qr.resid(qr(X[i2, ]), matrix(rnorm(sum(i2) * 3, sd = 3), ncol = 3))
  for (type in c("dct", "map")) {
    cfg <- list(hrf = c(0, 1, .4), type = type, method = "fixed", lambda = 0.1, K = 10)
    a <- mppi(Y, X, 3, runs = runs, domain = "neural", deconv = cfg, engine = "instant")
    b <- mppi(Y2, X, 3, runs = runs, domain = "neural", deconv = cfg, engine = "instant")
    expect_equal(a$U[!i2, ], b$U[!i2, ], tolerance = 1e-10)
    expect_equal(unname(a$U), unname(mppi(Y, X, 3, runs = runs, domain = "neural", deconv = cfg)$U))
  }
})

test_that("instant neural domain rejects censored scans", {
  set.seed(13); n <- 50; X <- cbind(1, task = rnorm(n)); Y <- matrix(rnorm(n * 3), n)
  Y[15, 1] <- NA
  expect_error(mppi(Y, X, 2, domain = "neural", engine = "instant"), "masked deconvolution")
})

test_that("instant evidence is finite with global and run intercepts", {
  set.seed(12); n <- 40; run <- rep(1:2, each = 20)
  X <- cbind(int = 1, run1 = run == 1, run2 = run == 2, task = sin(seq_len(n)))
  Y <- matrix(rnorm(n * 3), n)
  fit <- mppi(Y, X, 4, runs = run, engine = "instant")
  expect_equal(fit$evidence$rss_base, colSums(qr.resid(qr(X[, -4]), Y)^2))
  expect_equal(fit$evidence$rss_full, colSums(qr.resid(qr(X), Y)^2))
  expect_true(is.finite(fit$evidence$delta_AIC_total))
})

test_that("instant evidence counts the design rank", {
  set.seed(12); n <- 40; run <- rep(1:2, each = 20)
  X <- cbind(int = 1, run1 = run == 1, run2 = run == 2, task = sin(seq_len(n)))
  Y <- matrix(rnorm(n * 3), n)
  fit <- mppi(Y, X, 4, runs = run, engine = "instant")
  expect_equal(fit$evidence$df_full, qr(X)$rank)
  expect_equal(fit$evidence$df_base, qr(X[, -4])$rank)
})

test_that("instant prewhitening follows the closed whitened path", {
  set.seed(15); n <- 80; runs <- rep(1:2, each = 40)
  X <- cbind(int = 1, run1 = runs == 1, run2 = runs == 2, task = rnorm(n)); Y <- matrix(rnorm(n * 3), n)
  fc <- mppi(Y, X, psych_idx = 4, runs = runs, prewhiten = TRUE, p = 1, scale = "cov", zero_diag = FALSE)
  fi <- mppi(Y, X, psych_idx = 4, runs = runs, prewhiten = TRUE, ar_order = 1, scale = "cov",
             engine = "instant", zero_diag = FALSE, instant = list(tau_half = 1e12))
  expect_equal(c(mppi_get_M(fi, 1)), c(mppi_get_M(fc, 1)), tolerance = 1e-8)
  expect_equal(fi$domain, "whitened_bold")
  expect_error(mppi(Y, X, psych_idx = 4, runs = runs, prewhiten = TRUE, ar_order = 1, domain = "neural",
                    engine = "instant"), "unsupported")
})

test_that("instant lag pairs use acquisition distances across censored scans", {
  set.seed(13); n <- 60; runs <- rep(1:2, each = 30)
  X <- cbind(1, task = rnorm(n)); Y <- matrix(rnorm(n * 3), n)
  Y[15, 1] <- NA; Y[31, 2] <- NA
  fc <- mppi(Y, X, 2, runs = runs, lags = c(-1, 2), zero_diag = FALSE, scale = "cov")
  fi <- mppi(Y, X, 2, runs = runs, engine = "instant", scale = "cov", zero_diag = FALSE,
             instant = list(tau_half = 1e12, offsets = c(-1L, 2L)))
  for (lg in c(-1L, 2L)) expect_equal(c(mppi_get_M_lag(fi, 1, lg)), c(mppi_get_M_lag(fc, 1, lg)), tolerance = 1e-8)
  expect_equal(fi$observation_index, fc$observation_index)
})

test_that("EWM slopes converge to the closed engine as tau_half grows, lag orientation included", {
  set.seed(10); n <- 240; runs <- rep(1:2, each = n / 2)
  task <- rep(rep(c(1, 0), each = 10), n / 20)
  X <- cbind(run1 = runs == 1, run2 = runs == 2, task = task)
  Y <- matrix(rnorm(n * 3), n)
  for (t in 2:n) Y[t, 2] <- Y[t, 2] + 2 * task[t - 1] * Y[t - 1, 1]
  fc <- mppi(Y, X, 3, runs = runs, lags = c(-2, 1), zero_diag = FALSE, scale = "cov")
  fi <- mppi(Y, X, 3, runs = runs, engine = "instant", scale = "cov", zero_diag = FALSE,
             instant = list(tau_half = 1e12, offsets = c(-2L, 1L)))
  for (lg in c(-2L, 0L, 1L)) expect_equal(c(mppi_get_M_lag(fi, 1, lg)), c(mppi_get_M_lag(fc, 1, lg)), tolerance = 1e-8)
  expect_equal(c(fi$Delta_normalized[[1]]), c(fc$Delta_normalized[[1]]), tolerance = 1e-8)
  M1 <- mppi_get_M_lag(fi, 1, 1)
  expect_gt(M1[1, 2], 2 * abs(M1[2, 1]))
})

test_that("EWM slopes halve their weight every tau_half scans back from the end", {
  set.seed(21); n <- 80; runs <- rep(1:2, each = 40)
  X <- cbind(1, task = rnorm(n)); Y <- matrix(rnorm(n * 3), n); Y[30, 2] <- NA
  tau <- 7
  fi <- mppi(Y, X, 2, runs = runs, engine = "instant", scale = "cov", zero_diag = FALSE,
             instant = list(tau_half = tau, offsets = c(0L, 2L)))
  keep <- setdiff(seq_len(n), 30); R <- fi$R; p <- fi$pk[[1]]; w <- 2^((keep - n) / tau)
  expect_equal(c(mppi_get_M(fi, 1)), c(crossprod(R, w * p * R) / sum(w * p^2)))
  a <- which((keep + 2) %in% keep & runs[keep] == runs[pmin(keep + 2, n)]); b <- match(keep[a] + 2, keep)
  expect_equal(c(mppi_get_M_lag(fi, 1, 2)), c(crossprod(R[a, ], w[a] * p[a] * R[b, ]) / sum(w[a] * p[a]^2)))
})
