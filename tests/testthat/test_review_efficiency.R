test_that("kernel slope energy equals the direct cross-product energy", {
  set.seed(51); R <- matrix(rnorm(40 * 12), 40)
  for (zd in c(TRUE, FALSE)) {
    fast <- multiPPI:::.mppi_slope_energy(R, zd, kernel = TRUE)
    slow <- multiPPI:::.mppi_slope_energy(R, zd, kernel = FALSE)
    for (i in 1:5) { z <- rnorm(40); expect_equal(fast(z), slow(z), tolerance = 1e-12) }
  }
})

test_that("omnibus nulls on the kernel path match a direct replay of the same draws", {
  set.seed(52); n <- 40; X <- cbind(1, rnorm(n)); Y <- matrix(rnorm(n * 12), n)
  fit <- mppi_fit(Y, X, 2, runs = rep(1:2, each = 20))
  om <- mppi_omnibus(fit, blksize = 5, B = 25, seed = 4)
  setup <- multiPPI:::.mppi_randomization_setup(fit, 5, 25, "block_flip", "rademacher")
  set.seed(4); pk <- fit$pk[[1]]
  ref <- replicate(25, {
    z <- multiPPI:::.mppi_random_psych(pk, setup, "block_flip", "rademacher")
    M <- crossprod(fit$R, z * fit$R) / sum(z^2); diag(M) <- 0; sum(M^2)
  })
  expect_equal(om[[1]]$Qnull, ref, tolerance = 1e-12)
})

test_that("lag selection residuals do not depend on the fit scale", {
  set.seed(53); n <- 60; X <- cbind(1, rnorm(n)); Y <- matrix(rnorm(n * 5), n)
  expect_equal(mppi_fit(Y, X, 2, scale = "cov")$R, mppi_fit(Y, X, 2, scale = "normalized")$R)
})

test_that("classify filters innovations within runs, once for all draws", {
  set.seed(4); n <- 240; runs <- rep(1:2, each = 120); task <- rep(rep(c(1, 0), each = 10), n / 20)
  fit <- mppi(matrix(rnorm(n * 4), n), cbind(1, task = task), 2, runs = runs, lags = -1:1)
  cl <- mppi_classify(fit, 1, mode = "innov", B = 5, seed = 1)
  X <- fit$R; pk <- fit$pk[[1]]; Xi <- X
  for (r in 1:2) {
    i <- runs == r
    Xi[i, ] <- mppi_innovations(X[i, ], protect = "lowK", K_keep = 8L, ar_order = "auto", p = 4L, center = FALSE)
  }
  expect_equal(cl$I, norm(crossprod(Xi, pk * Xi), "F") / norm(crossprod(X, pk * X), "F"), tolerance = 1e-8)
})

test_that("sdPPI fit window covariances use only the window rows", {
  set.seed(5); n <- 200; task <- rep(rep(c(1, 0), each = 10), n / 20)
  fit <- mppi(matrix(rnorm(n * 5), n), cbind(1, task = task), 2)
  sf <- mppi_sdppi_fit(fit, 1, onsets_idx = c(41, 101), Wpre = 8, Wpost = 12)
  expect_equal(sf$data$ppi_energy, sapply(c(41, 101), function(t0) norm(crossprod(fit$R[t0:(t0 + 11), ]) / 12, "F")))
})
