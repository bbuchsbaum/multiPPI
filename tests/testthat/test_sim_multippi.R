context("Simulation diagnostics for multiPPI")

testthat::skip_if_not_installed("fmrihrf")

`%||%` <- function(a, b) if (!is.null(a)) a else b
zscore <- function(x) { sx <- sd(x); if (sx == 0) x*0 else (x - mean(x))/sx }

# Canonical HRF via fmrihrf (SPM-like double gamma)
h_spm <- function(TR = 2, L = 32) {
  t <- seq(0, L, by = TR)
  fmrihrf::hrf_spmg1(t)
}

conv_hrf <- function(x, h) {
  y <- stats::filter(x, filter = h, method = "convolution", sides = 1)
  y <- y[seq_along(x)]
  y[is.na(y)] <- 0
  y
}

lag_vec <- function(x, tau) {
  T <- length(x)
  if (tau > 0) c(rep(0, tau), x[1:(T - tau)])
  else if (tau < 0) c(x[(1 - tau):T], rep(0, -tau))
  else x
}

make_events <- function(T, K, n_events = rep(20L, K), seed = 1L) {
  set.seed(seed)
  if (length(n_events) == 1L) n_events <- rep(n_events, K)
  times <- seq(4, T - 4, length.out = sum(n_events))
  times <- sample(times)
  S <- matrix(0, T, K)
  idx <- 1
  for (k in seq_len(K)) {
    tk <- sort(round(times[idx:(idx + n_events[k] - 1)]))
    tk <- pmax(1, pmin(T, tk))
    S[tk, k] <- 1
    idx <- idx + n_events[k]
  }
  colnames(S) <- paste0("cond", seq_len(K))
  S
}

simulate_pair <- function(T = 300, TR = 2, S, gamma, baseline = 0, tau = 0,
                          seed_sd = 1, targ_sd = 1, drift_amp = 0, jitter_sd = 0) {
  K <- ncol(S)
  x1 <- zscore(as.numeric(arima.sim(n = T, list(ar = 0.3), sd = seed_sd)))
  if (jitter_sd > 0) {
    j <- stats::filter(rnorm(T, 0, jitter_sd), rep(1/3, 3), sides = 2)
    j[is.na(j)] <- 0
    x1j <- zscore(x1 + j)
  } else x1j <- x1
  seedd <- lag_vec(x1j, tau)
  g_t <- as.numeric(S %*% gamma)
  eps <- as.numeric(arima.sim(n = T, list(ar = 0.2), sd = targ_sd))
  slow <- if (drift_amp > 0) drift_amp * sin(2*pi*(1:T) / (T/3)) else 0
  x2 <- eps + slow + seedd * (baseline + g_t)
  x1 <- zscore(x1); x2 <- zscore(x2)
  h  <- h_spm(TR)
  y1 <- zscore(conv_hrf(x1, h) + rnorm(T, 0, 0.2))
  y2 <- zscore(conv_hrf(x2, h) + rnorm(T, 0, 0.2))
  list(U = cbind(x1, x2), Y = cbind(y1, y2), S = S, TR = TR)
}

Mk <- function(U, s) {
  den <- sum(s^2) + 1e-12
  crossprod(U, s * U) / den
}

assert_close <- function(x, target, tol, msg) {
  testthat::expect_true(all(is.finite(x)), info = paste(msg, "non-finite"))
  diff <- max(abs(x - target))
  testthat::expect_true(diff < tol, info = sprintf("%s | diff=%.3f target=%.3f tol=%.3f", msg, diff, target, tol))
}

context("Simulation diagnostics for neural vs BOLD PPI, gPPI, lagged routing, iPPI")

test_that("Neural vs BOLD PPI recovers stronger effect", {
  T  <- 240; TR <- 2
  S  <- make_events(T, K = 2, n_events = c(18, 18), seed = 13)
  gamma <- c(0.8, 0.0)
  sim <- simulate_pair(T, TR, S, gamma, baseline = 0, tau = 0, jitter_sd = 0.5)
  sA <- S[,1]; U <- sim$U; Y <- sim$Y
  M_neu <- Mk(U, sA); M_bld <- Mk(Y, sA)
  neu <- abs(M_neu[1,2]); bld <- abs(M_bld[1,2])
  testthat::expect_gt(neu, bld + 0.05)
})

test_that("No PPI vs PPI (two conditions)", {
  T <- 260; TR <- 2; S <- make_events(T, 2, n_events = c(20, 20), seed = 7)
  sim0 <- simulate_pair(T, TR, S, gamma = c(0,0), baseline = 0, tau = 0)
  M0A  <- Mk(sim0$U, S[,1])[1,2]; M0B <- Mk(sim0$U, S[,2])[1,2]
  sim1 <- simulate_pair(T, TR, S, gamma = c(1.2,0), baseline = 0, tau = 0, seed_sd = 0.5, targ_sd = 0.4)
  M1A <- Mk(sim1$U, S[,1])[1,2]; M1B <- Mk(sim1$U, S[,2])[1,2]
  testthat::expect_gt(abs(M1A), abs(M0A) + 0.15)
  testthat::expect_lt(abs(M1B), abs(M0B) + 0.05)
})

test_that("gPPI vs sPPI (three conditions)", {
  T <- 200; TR <- 2; S <- make_events(T, 3, n_events = c(15, 15, 15), seed = 11)
  gamma <- c(0.6, 0.3, 0.5)
  sim <- simulate_pair(T, TR, S, gamma = gamma, baseline = 0, tau = 0)
  U <- sim$U
  MA <- Mk(U, S[,1]); MB <- Mk(U, S[,2])
  gppi_est <- (MA - MB)[1,2]
  target   <- gamma[1] - gamma[2]
  assert_close(gppi_est, target, tol = 0.08, msg = "gPPI estimate of A-B")
  s_coll <- S[,1] - S[,2]
  sppi_est <- Mk(U, s_coll)[1,2]
  err_gppi <- abs(gppi_est - target)
  err_sppi <- abs(sppi_est - target)
  testthat::expect_gt(err_sppi, err_gppi)
})

test_that("Lagged routing energy distinguishes direction", {
  T <- 200
  set.seed(99)
  seed <- zscore(arima.sim(n = T, list(ar = 0.4), sd = 1))
  target_ff <- zscore(lag_vec(seed, 1) + rnorm(T, 0, 0.1))
  target_fb <- zscore(lag_vec(seed, -1) + rnorm(T, 0, 0.1))
  U_ff <- cbind(seed, target_ff)
  U_fb <- cbind(seed, target_fb)
  lag_energy <- function(U, lags = -2:2) {
    pos <- neg <- 0
    for (tau in lags) {
      if (tau == 0) next
      if (tau > 0) {
        idx <- seq_len(T - tau)
        M <- crossprod(U[idx, , drop = FALSE], U[idx + tau, , drop = FALSE]) / length(idx)
        pos <- pos + sqrt(sum(M^2))
      } else {
        idx <- seq_len(T + tau)
        M <- crossprod(U[idx - tau, , drop = FALSE], U[idx, , drop = FALSE]) / length(idx)
        neg <- neg + sqrt(sum(M^2))
      }
    }
    list(pos = pos, neg = neg)
  }
  e_ff <- lag_energy(U_ff)
  e_fb <- lag_energy(U_fb)
  testthat::expect_gt(e_ff$pos, e_fb$pos)
  testthat::expect_true(e_fb$neg + 0.05 >= e_ff$neg)
})

test_that("Innovations reduce slow-drift bias", {
  skip_if_not(exists("mppi_innovations"))
  T <- 220; TR <- 2; S <- make_events(T, 1, n_events = 20, seed = 33)
  gamma <- 0.5
  sim <- simulate_pair(T, TR, S, gamma = gamma, baseline = 0, tau = 0, drift_amp = 1.0)
  U <- sim$U; s <- S[,1]
  M_raw <- Mk(U, s)[1,2]
  U_i <- mppi_innovations(U, TR = TR, protect = "lowHz", f_cut = 0.03,
                          ar_order = "auto", p = 4)
  M_i <- Mk(U_i, s)[1,2]
  err_raw <- abs(M_raw - gamma)
  err_ipi <- abs(M_i   - gamma)
  testthat::expect_true(err_ipi <= err_raw + 0.01)
})
