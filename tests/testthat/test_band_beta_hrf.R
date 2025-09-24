context("Bandpass, beta-series, HRF adaptive helpers")

set.seed(123)

# Bandpass -----------------------------------------------------------------

test_that("mppi_bandpass preserves dimensions and attenuates highs", {
  Tn <- 120
  t <- seq(0, 2*pi, length.out = Tn)
  slow <- sin(t)
  fast <- sin(10 * t)
  X <- cbind(slow = slow + fast, fast = fast)
  X_lp <- mppi_bandpass(X, low = NULL, high = 0.2)
  expect_equal(dim(X_lp), dim(X))
  # Low-pass should retain more slow energy than fast column
  slow_var <- var(X_lp[,1])
  fast_var <- var(X_lp[,2])
  expect_gt(slow_var, fast_var)
})

test_that("mppi_fit_multi_band returns fits per band", {
  Tn <- 60; V <- 4
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(1, sin(seq(0, pi, length.out = Tn)))
  pidx <- 2L
  bands <- list(list(low = NULL, high = 0.1, name = "low"),
                list(low = 0.15, high = NULL, name = "high"))
  fits <- mppi_fit_multi_band(Y, X, psych_idx = pidx, bands = bands)
  expect_named(fits, c("low", "high"))
  expect_equal(length(fits$low$Delta), 1L)
  expect_equal(dim(fits$low$Delta[[1]]), c(V, V))
})

# Beta-series ---------------------------------------------------------------

test_that("mppi_fit_beta and permutations run on toy data", {
  n_trials <- 12; V <- 3
  B <- matrix(rnorm(n_trials * V), n_trials, V)
  runs <- rep(1:3, length.out = n_trials)
  Pbeta <- matrix(rnorm(n_trials), ncol = 1)
  colnames(Pbeta) <- "psych"
  fit_beta <- mppi_fit_beta(B, Pbeta, scale = "cov")
  expect_equal(length(fit_beta$Delta), 1L)
  perm_beta <- mppi_beta_permute(B, Pbeta, run_trial = runs, Bperm = 19)
  expect_equal(length(perm_beta), 1L)
  expect_true(is.finite(perm_beta[[1]]$p_global))
})

# HRF adaptive --------------------------------------------------------------

test_that("mppi_hrf_adapt outputs combined matrix and weights", {
  V <- 3
  Delta1 <- diag(V)
  Delta2 <- diag(V); Delta2[1,2] <- Delta2[2,1] <- 0.5
  combo <- mppi_hrf_adapt(list(Delta1, Delta2), pk_norms = c(2, 1))
  expect_equal(length(combo$weights), 2L)
  expect_true(any(combo$weights != 0))
  expect_equal(dim(combo$Delta), c(V, V))
})
