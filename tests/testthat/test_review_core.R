context("Review fixes: grouping, basis back-projection, dataset path, precision gate, GLS")

.review_basis_data <- function() {
  set.seed(3)
  Tn <- 300L; V <- 6L
  X <- cbind(`(Intercept)` = 1, A = rep(c(0, 1), each = 15, length.out = Tn))
  Y <- matrix(rnorm(Tn * V), Tn, V)
  colnames(Y) <- paste0("r", seq_len(V))
  f <- rnorm(Tn) * (1 + 1.5 * X[, 2])
  Y <- Y + outer(f, rep(1, V))
  Vb <- qr.Q(qr(matrix(rnorm(V * V), V, V)))
  rownames(Vb) <- colnames(Y)
  list(X = X, Y = Y, Vb = Vb, basis = as_mppi_basis(Vb, name = "rot", orthonormalize = FALSE))
}

test_that("neural-domain grouping keeps fmridesign conditions distinct", {
  skip_if_not_installed("fmridesign")
  suppressPackageStartupMessages(library(fmridesign))
  set.seed(11)
  sf <- fmridesign::sampling_frame(blocklens = 120, TR = 2)
  ev <- data.frame(onset = seq(4, 220, by = 12),
                   cond = factor(rep(c("task_easy", "task_hard"), length.out = 19)), run = 1)
  em1 <- fmridesign::event_model(onset ~ hrf(cond), data = ev, block = ~ run, sampling_frame = sf)
  X1 <- as.matrix(fmridesign::design_matrix(em1))
  em3 <- fmridesign::event_model(onset ~ hrf(cond, basis = "spmg3"), data = ev, block = ~ run, sampling_frame = sf)
  X3 <- as.matrix(fmridesign::design_matrix(em3))
  Y <- matrix(rnorm(120 * 4), 120, 4)
  dc <- list(tr = 2, K = 40, method = "fixed", lambda = 0.1)
  fit1 <- mppi_fit(Y, X1, psych_idx = 1:2, domain = "neural", deconv = dc)
  expect_equal(fit1$names, paste0("s:", colnames(X1)))
  g3 <- mppi_group_hrf_columns(X3)
  expect_equal(unname(g3), list(1:3, 4:6))
  fit3 <- mppi_fit(Y, X3, psych_idx = 1:6, domain = "neural", deconv = dc)
  expect_equal(fit3$names, c(paste0("s:", paste(colnames(X3)[1:3], collapse = "+")),
                             paste0("s:", paste(colnames(X3)[4:6], collapse = "+"))))
  fitm <- mppi_fit(Y, X1, psych_idx = 1:2, domain = "neural", deconv = c(dc, type = "map"))
  expect_equal(length(fitm$Delta), 2L)
})

test_that("unnamed, empty and duplicated names never merge distinct columns", {
  set.seed(12)
  Tn <- 160L
  X <- cbind(rbinom(Tn, 1, 0.3), rbinom(Tn, 1, 0.3))
  Y <- matrix(rnorm(Tn * 4), Tn, 4)
  dc <- list(tr = 2, K = 40, method = "fixed", lambda = 0.1)
  fit <- mppi_fit(Y, X, psych_idx = 1:2, domain = "neural", deconv = dc)
  expect_equal(length(fit$Delta), 2L)
  expect_equal(length(mppi_group_hrf_columns(cbind(a = 1, 2, 3))), 3L)
  expect_equal(length(mppi_group_hrf_columns(cbind(a_b01 = 1, a_b01 = 2, a_b02 = 3))), 3L)
  expect_equal(length(mppi_group_hrf_columns(cbind(task_easy = 1, task_hard = 2, task.x = 3))), 3L)
  legacy <- mppi_group_hrf_columns(cbind(hrf_task_easy = 1, hrfTD_task_easy = 2, hrf_task_hard = 3, hrfTD_task_hard = 4))
  expect_equal(legacy, list(task_easy = 1:2, task_hard = 3:4))
})

test_that("full-rank rotated basis reproduces the ROI fit by default", {
  d <- .review_basis_data()
  for (zd in c(TRUE, FALSE)) for (sc in c("normalized", "cov")) {
    roi <- mppi_fit(d$Y, d$X, 2L, scale = sc, zero_diag = zd)
    bas <- mppi_fit(d$Y, d$X, 2L, scale = sc, zero_diag = zd, basis = d$basis)
    ref <- mppi_get_M(roi, 1L)
    expect_equal(mppi_reconstruct_delta(bas, 1L), ref, tolerance = 1e-10, check.attributes = FALSE)
    expect_equal(mppi_edge(bas, 1L, 2L, 5L), ref[2, 5], tolerance = 1e-10)
    expect_equal(mppi_reconstruct_delta(bas, 1L, rows = 1:3, cols = c(2, 5)),
                 ref[1:3, c(2, 5)], tolerance = 1e-10)
  }
})

test_that("rotated basis reproduces every ROI scale for edges, blocks and seed slopes", {
  d <- .review_basis_data()
  for (zd in c(TRUE, FALSE)) {
    roi <- mppi_fit(d$Y, d$X, 2L, zero_diag = zd, lags = c(0L, 1L))
    bas <- mppi_fit(d$Y, d$X, 2L, zero_diag = zd, lags = c(0L, 1L), basis = d$basis)
    for (m in c("raw", "amplitude", "normalized")) {
      ref <- mppi_get_M_scaled(roi, 1L, m)
      expect_equal(mppi_reconstruct_delta(bas, 1L, mode = m), ref, tolerance = 1e-10)
      expect_equal(mppi_edge(bas, 1L, 2L, 5L, mode = m), ref[2, 5], tolerance = 1e-10)
      expect_equal(mppi_edge(bas, 1L, 3L, 3L, mode = m), ref[3, 3], tolerance = 1e-10)
      for (lg in 0:1) {
        expect_equal(mppi_seedset_slopes(bas, seeds = 1:2, k = 1L, mode = m, lag = lg),
                     mppi_seedset_slopes(roi, seeds = 1:2, k = 1L, mode = m, lag = lg), tolerance = 1e-10)
      }
      expect_equal(mppi_seedset_slopes(bas, seeds = 1:3, k = 1L, target_idx = 1:3, mode = m),
                   mppi_seedset_slopes(roi, seeds = 1:3, k = 1L, target_idx = 1:3, mode = m), tolerance = 1e-10)
    }
  }
})

test_that("rank-reduced basis back-projects the in-subspace slope with voxel SDs", {
  d <- .review_basis_data()
  P <- d$Vb[, 1:3] %*% t(d$Vb[, 1:3])
  bas <- mppi_fit(d$Y, d$X, 2L, basis = as_mppi_basis(d$Vb[, 1:3], orthonormalize = FALSE))
  proj <- mppi_fit(d$Y %*% P, d$X, 2L)
  raw <- mppi_reconstruct_delta(bas, 1L, mode = "raw")
  expect_equal(raw, mppi_get_M_scaled(proj, 1L, "raw"), tolerance = 1e-10, check.attributes = FALSE)
  s <- matrixStats::colSds(qr.resid(qr(d$X), d$Y))
  expect_equal(mppi_reconstruct_delta(bas, 1L), raw / tcrossprod(s), tolerance = 1e-10)
})

test_that("dataset path matches the direct design path and propagates its TR", {
  set.seed(4)
  Tn <- 200L; V <- 6L
  X <- cbind(`(Intercept)` = 1, A = rbinom(Tn, 1, 0.4))
  Y <- matrix(rnorm(Tn * V), Tn, V)
  # fmridataset is not installed here: a minimal object with the fields mppi.fmri_dataset reads
  ds <- structure(list(Y = Y, TR = 2, runs = rep(1L, Tn)), class = "fmri_dataset")
  des <- structure(list(X = X, psych_idx = 2L, runs = rep(1L, Tn), basis = NULL), class = "mppi_design")
  fit <- mppi(ds, des)
  expect_null(fit$basis)
  expect_equal(fit$Delta, mppi(des, Y = Y)$Delta)
  dc <- list(K = 40, method = "fixed", lambda = 0.1)
  fn <- mppi(ds, des, domain = "neural", deconv = dc)
  expect_equal(fn$deconv$tr, 2)
  expect_equal(fn$deconv$hrf, mppi_default_hrf(tr = 2))
  expect_equal(fn$Delta, mppi(des, Y = Y, domain = "neural", deconv = c(dc, tr = 2))$Delta)
  expect_error(mppi(ds, des, domain = "neural", deconv = c(dc, tr = 1)), "TR mismatch")
  expect_equal(mppi(ds, des, basis = list(type = "pca", r = 3L))$basis$r, 3L)
})

test_that("precision gate matches fit$partial and is unit-free when normalized", {
  set.seed(5)
  Tn <- 300L; V <- 5L
  X <- cbind(`(Intercept)` = 1, A = rbinom(Tn, 1, 0.4), B = rnorm(Tn))
  Y <- matrix(rnorm(Tn * V), Tn, V)
  Y[, 2] <- Y[, 2] + X[, 2] * Y[, 1]
  fit <- mppi_fit(Y, X, 2:3, scale = "cov")
  gate <- mppi_precision_gate(fit, mode = "raw")
  for (k in 1:2) expect_equal(gate$Delta[[k]], fit$partial[[k]]$DeltaTheta, tolerance = 1e-10)
  gate2 <- mppi_precision_gate(fit, mode = "raw", ridge = fit$partial_lambda)
  expect_equal(gate2$Delta[[1]], fit$partial[[1]]$DeltaTheta, tolerance = 1e-10)
  fit10 <- mppi_fit(10 * Y, X, 2:3, scale = "cov")
  expect_equal(mppi_precision_gate(fit10)$Delta, mppi_precision_gate(fit)$Delta, tolerance = 1e-8)
  D <- diag(fit$sigma)
  expect_equal(mppi_precision_gate(fit, ridge = 0)$Delta[[1]],
               D %*% mppi_precision_gate(fit, mode = "raw", ridge = 0)$Delta[[1]] %*% D, tolerance = 1e-8)
  expect_error(mppi_precision_gate(fit, mode = "amplitude"), "no precision-domain analogue")
})

test_that("whitened fit is GLS on the fmriAR-whitened design", {
  skip_if_not_installed("fmriAR")
  set.seed(7)
  Tn <- 200L; V <- 3L
  runs <- rep(1:2, each = 100L)
  A <- rbinom(Tn, 1, 0.4)
  Y <- matrix(rnorm(Tn * V), Tn, V)
  for (j in seq_len(V)) for (r in 1:2) {
    i <- which(runs == r)
    Y[i, j] <- stats::filter(Y[i, j], 0.6, method = "recursive")
  }
  X <- cbind(`(Intercept)` = 1, A = A)
  fit <- mppi_fit_whitened(Y, X, runs = runs, psych_idx = 2L)
  w <- fmriAR::whiten_apply(fit$noise_plan, X = X, Y = Y, runs = runs)
  gls <- lm.fit(w$X, w$Y)$residuals
  expect_equal(unname(fit$R), unname(gls), tolerance = 1e-8)
  expect_equal(ncol(fit$fitted_X), ncol(X))
  fit_noint <- mppi_fit_whitened(Y, cbind(A = A), runs = runs, psych_idx = 1L)
  w2 <- fmriAR::whiten_apply(fit_noint$noise_plan, X = cbind(1, A), Y = Y, runs = runs)
  expect_equal(unname(fit_noint$R), unname(lm.fit(w2$X, w2$Y)$residuals), tolerance = 1e-8)
  expect_equal(fit_noint$names, "A")
})
