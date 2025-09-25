test_that("ewm_gated_cross matches brute force", {
  set.seed(1)
  Tn <- 40
  r <- 4
  K <- 2
  offsets <- c(0L, 1L, -2L)
  U <- matrix(rnorm(Tn * r), Tn, r)
  pk <- matrix(rnorm(Tn * K), Tn, K)
  tau <- 6
  eta <- log(2) / tau
  cpp_out <- ewm_gated_cross_rcpp(U, pk, eta, offsets, normalized = FALSE)
  a <- exp(-eta)
  one_ma <- 1 - a
  brute <- array(0, dim = c(K, length(offsets), r, r))
  v_ref <- rep(0, r)
  for (t in seq_len(Tn) - 1L) {
    ut <- U[t + 1L, ]
    for (k in seq_len(K)) {
      g <- pk[t + 1L, k]
      for (li in seq_along(offsets)) {
        brute[k, li, , ] <- a * brute[k, li, , ]
        if (is.na(g)) next
        off <- offsets[li]
        idx <- t - off
        if (idx < 0) idx <- 0
        if (idx >= Tn) idx <- Tn - 1L
        yt <- U[idx + 1L, ]
        brute[k, li, , ] <- brute[k, li, , ] + one_ma * g * tcrossprod(ut, yt)
      }
    }
    v_ref <- a * v_ref + one_ma * (ut^2)
  }
  for (k in seq_len(K)) {
    for (li in seq_along(offsets)) {
      idx <- (k - 1L) * length(offsets) + li
      expect_equal(cpp_out[[idx]], brute[k, li, , ], tolerance = 1e-8)
    }
  }

  cpp_norm <- ewm_gated_cross_rcpp(U, pk, eta, offsets, normalized = TRUE)
  denom <- sqrt(outer(pmax(v_ref, 1e-12), pmax(v_ref, 1e-12)))
  for (k in seq_len(K)) {
    for (li in seq_along(offsets)) {
      idx <- (k - 1L) * length(offsets) + li
      expect_equal(cpp_norm[[idx]], brute[k, li, , ] / denom, tolerance = 1e-8)
    }
  }
})

test_that("ewm_gated_cross respects blocklens boundaries", {
  set.seed(11)
  Tn <- 60
  r <- 3
  K <- 2
  offsets <- c(-2L, 0L, 2L)
  U <- matrix(rnorm(Tn * r), Tn, r)
  pk <- matrix(rnorm(Tn * K), Tn, K)
  tau <- 5
  eta <- log(2) / tau
  blocklens <- c(30L, 30L)
  cpp_out <- ewm_gated_cross_rcpp(U, pk, eta, offsets, normalized = FALSE, blocklens = blocklens)
  a <- exp(-eta)
  one_ma <- 1 - a
  brute <- array(0, dim = c(K, length(offsets), r, r))
  run_id <- rep.int(seq_along(blocklens), blocklens)
  for (t in seq_len(Tn) - 1L) {
    ut <- U[t + 1L, ]
    for (k in seq_len(K)) {
      g <- pk[t + 1L, k]
      for (li in seq_along(offsets)) {
        brute[k, li, , ] <- a * brute[k, li, , ]
        if (is.na(g)) next
        off <- offsets[li]
        idx <- t - off
        if (idx < 0) idx <- 0
        if (idx >= Tn) idx <- Tn - 1L
        if (run_id[t + 1L] != run_id[idx + 1L]) next
        yt <- U[idx + 1L, ]
        brute[k, li, , ] <- brute[k, li, , ] + one_ma * g * tcrossprod(ut, yt)
      }
    }
  }
  for (k in seq_len(K)) {
    for (li in seq_along(offsets)) {
      idx <- (k - 1L) * length(offsets) + li
      expect_equal(cpp_out[[idx]], brute[k, li, , ], tolerance = 1e-8)
    }
  }
})

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
