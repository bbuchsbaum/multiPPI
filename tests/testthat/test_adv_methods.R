context("Advanced mechanistic methods")

make_fit <- function(U, pk) {
  structure(list(U = U, pk = list(pk), names = "task"), class = "mppi_fit")
}

test_that("mppi_innovations respects AR order and low-frequency protection", {
  set.seed(1)
  U <- matrix(rnorm(80), nrow = 40, ncol = 2)
  E <- mppi_innovations(U, protect = "none", ar_order = "fixed", p = 0L, center = FALSE)
  expect_equal(E, U, check.attributes = FALSE)

  slow <- rep(1, 40)
  fast <- sin(seq(0, 4 * pi, length.out = 40))
  U2 <- cbind(slow, fast)
  E2 <- mppi_innovations(U2, protect = "lowK", K_keep = 1L,
                         ar_order = "fixed", p = 0L, center = FALSE)
  expect_equal(E2[, 1], slow, tolerance = 1e-10)

  expect_error(mppi_innovations(U2, protect = "lowHz", f_cut = 0.05,
                                ar_order = "auto"),
               "TR required", fixed = TRUE)

  Eauto <- mppi_innovations(U2, TR = 0.8, protect = "lowHz", f_cut = 0.05,
                            ar_order = "auto", center = TRUE)
  expect_equal(attr(Eauto, "orders") >= 0L, rep(TRUE, ncol(U2)))
})

test_that("spectral and precision summaries recover analytical values", {
  U <- matrix(c(1, 0,
                0, 1), nrow = 2, byrow = TRUE)
  pk <- c(1, 0)
  fit <- make_fit(U, pk)
  spec <- mppi_spectral_summary(fit, 1L)
  expect_equal(spec$lambda1_share, 1, tolerance = 1e-10)
  expect_equal(spec$effrank, 1, tolerance = 1e-8)
  expect_equal(spec$sign_balance, 1)

  prec <- mppi_precision_delta(fit, 1L, lambda = 1e-3)
  Mk <- crossprod(U, pk * U) / sum(pk^2)
  S0 <- crossprod(U) / nrow(U)
  Theta0 <- solve(S0 + diag(1e-3, ncol(S0)))
  dTheta <- -Theta0 %*% Mk %*% Theta0
  expect_equal(prec$dTheta, dTheta)
  expect_equal(prec$offdiag_energy, 0, tolerance = 1e-12)
})

test_that("lagged summaries match manual calculations", {
  U <- matrix(c(1, 0,
                0, 1), nrow = 2, byrow = TRUE)
  pk <- c(1, 0)
  fit <- make_fit(U, pk)
  attr(fit$U, "TR") <- 1

  lag_neural <- mppi_lagged(fit, 1L, lags = -1:1, mode = "neural")
  expect_equal(lag_neural$M_lag[["0"]], crossprod(U, pk * U) / sum(pk^2))
  expect_equal(lag_neural$M_lag[["1"]][1, 2], 1)
  expect_equal(lag_neural$M_lag[["-1"]][2, 1], 1)
  expect_equal(lag_neural$R, 0)

  lag_innov <- mppi_lagged(fit, 1L, lags = -1:1, mode = "innov",
                           protect = "none", ar_order = "fixed", p = 0L)
  expect_equal(lag_innov$M_lag, lag_neural$M_lag)
})

test_that("mppi_classify captures gain-dominated patterns", {
  set.seed(2)
  Tn <- 40L
  pk <- rep(c(1, 0), length.out = Tn)
  U <- cbind(pk, rnorm(Tn, sd = 1e-3))
  fit <- make_fit(U, pk)
  lab <- mppi_classify(fit, 1L, lambda = 1e-3, lags = -1:1,
                       mode = "neural", B = 30L, seed = 11)
  expect_true(lab$G > 0.8)
  expect_lt(lab$C, 1e-2)
  expect_lt(abs(lab$R), 1e-6)
  expect_equal(lab$lagged$R, lab$R)

  lab_innov <- mppi_classify(fit, 1L, lambda = 1e-3, lags = -1:1,
                             mode = "innov", B = 30L, seed = 12)
  expect_equal(lab_innov$label, lab$label)
})

test_that("state-dependent PPI returns trial summaries with nuisances", {
  set.seed(3)
  Tn <- 60L
  U <- matrix(rnorm(Tn * 3), nrow = Tn, ncol = 3)
  pk <- rep(c(1, -1), length.out = Tn)
  fit <- make_fit(U, pk)
  runs <- rep(1:3, each = 20)
  onsets <- c(10L, 20L, 30L, 40L, 50L)
  gs <- sin(seq(0, 2 * pi, length.out = Tn))
  sd <- mppi_sdppi_fit(fit, 1L, onsets_idx = onsets, runs = runs,
                       Wpre = 4L, Wpost = 4L, mode = "neural",
                       nuisance_ts = list(GS = gs))
  expect_equal(nrow(sd$data), length(onsets))
  expect_true(all(c("pre_lambda1_share", "ppi_energy", "run") %in% colnames(sd$data)))
  expect_true(ncol(sd$data) >= 10)
  expect_s3_class(sd$model_ppi, "lm")
  expect_true(is.numeric(sd$cv_r))

  attr(fit$U, "TR") <- 0.8
  sd_innov <- mppi_sdppi_fit(fit, 1L, onsets_idx = onsets, runs = runs,
                             Wpre = 4L, Wpost = 4L, mode = "innov",
                             protect = "lowHz", f_cut = 0.05)
  expect_equal(nrow(sd_innov$data), length(onsets))
})

test_that("null onsets avoid guarded windows", {
  set.seed(4)
  obs <- c(10L, 30L)
  guard <- 4L
  null_idx <- mppi_make_null_onsets(60L, obs, guard = guard)
  expect_equal(length(null_idx), length(obs))
  expect_true(all(abs(outer(null_idx, obs, "-")) > guard))

  expect_error(mppi_make_null_onsets(12L, c(4L, 8L), guard = 5L),
               "Not enough candidate TRs", fixed = TRUE)
})

test_that("state predictors provide cross-validated performance", {
  set.seed(5)
  Tn <- 60L
  U <- matrix(rnorm(Tn * 2), nrow = Tn, ncol = 2)
  pk <- rep(c(1, 0), length.out = Tn)
  fit <- make_fit(U, pk)
  runs <- rep(1:3, each = 20)
  onsets <- c(12L, 24L, 36L, 48L)
  sd <- mppi_sdppi_fit(fit, 1L, onsets_idx = onsets, runs = runs,
                       Wpre = 4L, Wpost = 4L, mode = "neural")
  res_gauss <- mppi_state_predict(sd$data, sd$data$ppi_energy, runs = sd$data$run,
                                  family = "gaussian")
  expect_true(is.numeric(res_gauss$cv_r))
  y_bin <- as.integer(sd$data$ppi_energy > median(sd$data$ppi_energy))
  res_bin <- mppi_state_predict(sd$data, y_bin, runs = sd$data$run,
                                family = "binomial")
  expect_true(is.numeric(res_bin$cv_auc) || is.na(res_bin$cv_auc))

  res_bin_flat <- mppi_state_predict(sd$data, rep(1L, nrow(sd$data)),
                                     runs = sd$data$run, family = "binomial")
  expect_true(is.na(res_bin_flat$cv_auc))
})

test_that("theta metrics and gating provide interpretable summaries", {
  Theta0 <- matrix(c(2, -1,
                     -1, 2), nrow = 2, byrow = TRUE)
  mts <- mppi_theta_metrics(Theta0)
  expect_true(mts$alg_conn >= 0)
  expect_true(mts$effrankA >= 1)
  expect_true(is.na(mts$Qmax) || mts$Qmax >= 0)
  expect_true(mts$comm_avg >= 0)

  set.seed(6)
  fits <- replicate(3, {
    base <- rnorm(30)
    U <- cbind(base, base + rnorm(30, sd = 0.2))
    pk <- rnorm(30)
    make_fit(U, pk)
  }, simplify = FALSE)
  gate <- mppi_gate_by_theta(fits, 1L, controls = data.frame(tSNR = c(90, 100, 110)))
  expect_equal(nrow(gate$data), length(fits))
  expect_true(all(c("G", "C", "alg_conn", "effrankA") %in% names(gate$data)))
  expect_s3_class(gate$model_G, "lm")
  expect_s3_class(gate$model_C, "lm")
})
