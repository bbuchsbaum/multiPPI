context("Basis handling and chunked ops")

test_that("as_mppi_basis orthonormalizes columns", {
  set.seed(1)
  V <- matrix(rnorm(30), nrow = 6, ncol = 5)
  B <- as_mppi_basis(V, name = "rand")
  XtX <- crossprod(B$V)
  expect_equal(diag(XtX), rep(1, ncol(XtX)), tolerance = 1e-8)
  expect_true(max(abs(XtX[upper.tri(XtX)])) < 1e-8)
  expect_equal(B$r, ncol(B$V))
})

test_that("mppi_edge matches reconstructed delta", {
  set.seed(2)
  Tn <- 40L; V <- 6L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  basis <- as_mppi_basis(diag(V), name = "id")
  fit <- mppi_fit(Y, X, psych_idx = 2L, basis = basis, backend = "blas")
  recon <- mppi_reconstruct_delta(fit, 1L)
  expect_equal(mppi_edge(fit, 1L, 2L, 5L), recon[2,5])
})

test_that("chunked projection matches blas projection", {
  set.seed(3)
  Tn <- 60L; V <- 7L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  basis <- as_mppi_basis(qr.Q(qr(matrix(rnorm(V * 4L), V, 4L))), name = "randPCA")
  fit_blas <- mppi_fit(Y, X, psych_idx = 2L, basis = basis,
                       backend = "blas", project_backend = "blas")
  fit_chunk <- mppi_fit(Y, X, psych_idx = 2L, basis = basis,
                        backend = "chunked", chunk_size = 16L,
                        project_backend = "chunked", project_chunk_cols = 3L)
  expect_equal(mppi_get_M_scaled(fit_chunk, 1L, mode = "normalized"),
               mppi_get_M_scaled(fit_blas, 1L, mode = "normalized"), tolerance = 1e-8)
})

test_that("neural domain matches direct conversion", {
  set.seed(7)
  Tn <- 60L; V <- 8L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  runs <- rep(1:3, each = 20)
  basis <- as_mppi_basis(diag(V), name = "id")
  h <- mppi_default_hrf(tr = 1)
  fit_neural <- mppi_fit(Y, X, psych_idx = 2L, runs = runs,
                         basis = basis, domain = "neural",
                         deconv = list(K = 32))
  fit_bold <- mppi_fit(Y, X, psych_idx = 2L, runs = runs,
                        basis = basis, domain = "bold")
  conv <- mppi_neural_from_fit(fit_bold, X, psych_idx = 2L, h = h, K = 32)
  expect_equal(mppi_get_M_scaled(fit_neural, 1L, mode = "raw"),
               mppi_get_M_scaled(conv, 1L, mode = "raw"), tolerance = 1e-6)
  expect_equal(fit_neural$deconv$hrf, h)
})

test_that("permutation seed is reproducible", {
  set.seed(4)
  Tn <- 45L; V <- 5L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  fit <- mppi_fit(Y, X, psych_idx = 2L)
  res1 <- mppi_omnibus(fit, blksize = 5L, B = 49L, seed = 123)
  res2 <- mppi_omnibus(fit, blksize = 5L, B = 49L, seed = 123)
  expect_equal(res1[[1]]$Qnull, res2[[1]]$Qnull)
  pm1 <- mppi_permute(fit, k = 1L, blksize = 5L, B = 49L, seed = 321)
  pm2 <- mppi_permute(fit, k = 1L, blksize = 5L, B = 49L, seed = 321)
  expect_equal(pm1, pm2)
})

test_that("freedman-lane permutation returns valid matrix", {
  set.seed(5)
  Tn <- 30L; V <- 4L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  fit <- mppi_fit(Y, X, psych_idx = 2L)
  pmat <- mppi_permute(fit, k = 1L, blksize = 5L, B = 29L,
                       method = "freedman_lane", seed = 77,
                       studentize = FALSE)
  expect_equal(dim(pmat), c(V, V))
  expect_true(all(is.na(diag(pmat))))
})
