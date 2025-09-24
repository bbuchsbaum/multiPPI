test_that("mppi_bandpass falls back to moving-average smoothing when signal missing", {
  set.seed(1)
  Tn <- 200L
  t <- seq_len(Tn)
  raw <- sin(2 * pi * t / 10) + 0.3 * sin(2 * pi * t / 2) # fast + slow
  X <- cbind(raw = raw)

  smoothed <- with_mocked_bindings(
    mppi_bandpass(X, low = NULL, high = 0.1, order = 2L),
    requireNamespace = function(pkg, quietly = TRUE) FALSE,
    .package = "base"
  )

  expect_equal(dim(smoothed), dim(X))
  # High-frequency energy should shrink after smoothing
  fft_raw <- Mod(stats::fft(X[,1]))
  fft_smooth <- Mod(stats::fft(smoothed[,1]))
  hf_idx <- 20:40
  expect_true(mean(fft_smooth[hf_idx]) < mean(fft_raw[hf_idx]))
})

test_that("mppi_fit_multi_band returns per-band fits", {
  set.seed(2)
  Tn <- 80L; V <- 3L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  bands <- list(list(low = NULL, high = 0.1, name = "low"),
                list(low = 0.1, high = 0.3, name = "mid"))

  fits <- mppi_fit_multi_band(Y, X, psych_idx = 2L, bands = bands, scale = "cov")
  expect_named(fits, c("low", "mid"))
  expect_true(all(vapply(fits, function(f) is.list(f) && length(f$Delta) == 1, logical(1))))
})

test_that("mppi_delta_partial rescales precision deltas to partial correlations", {
  set.seed(3)
  Sigma0 <- matrix(c(1, 0.2, 0.2, 1), 2, 2)
  DeltaSigma <- matrix(c(0.05, -0.04, -0.04, 0.03), 2, 2)
  lambda <- 1e-3
  DeltaTheta <- mppi_to_partial(Sigma0, DeltaSigma, lambda = lambda)
  Theta0 <- solve(Sigma0 + lambda * diag(2))
  DeltaRho <- mppi_delta_partial(Theta0, DeltaTheta)

  # approximate finite difference check
  eps <- 1e-4
  Theta_eps <- solve(Sigma0 + eps * DeltaSigma + lambda * diag(2))
  rho0 <- -Theta0[1,2] / sqrt(Theta0[1,1] * Theta0[2,2])
  rho_eps <- -Theta_eps[1,2] / sqrt(Theta_eps[1,1] * Theta_eps[2,2])
  expect_equal(rho_eps - rho0, eps * DeltaRho[1,2], tolerance = 1e-5)
})

test_that("chunked crossprod equals blas crossprod", {
  set.seed(6)
  Tn <- 70L; V <- 6L
  R <- matrix(rnorm(Tn * V), Tn, V)
  w <- rnorm(Tn)
  A <- .mppi_wcp(R, w)
  B <- .mppi_crossprod(R, w, backend = "chunked", chunk_size = 15L)
  expect_equal(A, B, tolerance = 1e-12)
})
