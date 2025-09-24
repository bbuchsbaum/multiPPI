test_that("mppi_fit reproduces closed-form covariance regression", {
  set.seed(42)
  Tn <- 60L; V <- 4L
  conf <- as.numeric(scale(rnorm(Tn), scale = FALSE))
  psych <- as.numeric(scale(rnorm(Tn), scale = FALSE))
  X <- cbind(Intercept = 1, conf = conf, psych = psych)
  Y <- matrix(rnorm(Tn * V), Tn, V)

  fit <- mppi_fit(Y, X, psych_idx = 3L, zero_diag = FALSE, scale = "cov")

  R_manual <- multiPPI:::`.mppi_residualize`(Y, X)
  Q <- X[, 1:2, drop = FALSE]
  pk_manual <- multiPPI:::`.mppi_residualize_vec`(X[, 3], Q)
  Delta_manual <- multiPPI:::`.mppi_wcp`(R_manual, pk_manual) / sum(pk_manual^2)

  expect_equal(mppi_get_M_scaled(fit, 1L, mode = "raw"), Delta_manual, tolerance = 1e-10)
  expect_equal(fit$R, R_manual, tolerance = 1e-10)
  expect_equal(fit$pk[[1]], pk_manual, tolerance = 1e-10)
})

test_that("HRF adaptive combination follows dominant energy direction", {
  Delta1 <- diag(c(1, 0))
  Delta2 <- diag(c(0, 4))
  combo <- mppi_hrf_adapt(list(Delta1, Delta2), pk_norms = c(1, 1), normalize = TRUE)
  expect_equal(combo$weights[1], 0, tolerance = 1e-10)
  scale_factor <- sum(combo$Delta * Delta2) / sum(Delta2^2)
  expect_equal(abs(scale_factor), 1, tolerance = 1e-8)
})

test_that("mppi_modes handles requested ranks", {
  Delta <- matrix(c(2, 1, 1, 2), 2, 2)
  res1 <- mppi_modes(Delta, r = 5L)
  expect_length(res1$values, 2L)
  expect_equal(ncol(res1$vectors), 2L)

  res2 <- mppi_modes(Delta, r = 1L)
  expect_length(res2$values, 1L)
  expect_equal(ncol(res2$vectors), 1L)
})

test_that("mppi_degree returns L1 and L2 summaries", {
  Delta <- matrix(c(0, 2, -3, 4), 2, 2)
  d1 <- mppi_degree(Delta, type = "L1")
  d2 <- mppi_degree(Delta, type = "L2")
  expect_equal(d1, rowSums(abs(Delta)))
  expect_equal(d2, sqrt(rowSums(Delta^2)))
})

test_that("accumulate backend matches blas backend", {
  set.seed(555)
  Tn <- 50L; V <- 4L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  fit_blas <- mppi_fit(Y, X, psych_idx = 2L, backend = "blas")
  fit_acc  <- mppi_fit(Y, X, psych_idx = 2L, backend = "accumulate")
  expect_equal(mppi_get_M_scaled(fit_blas, 1L, mode = "normalized"),
               mppi_get_M_scaled(fit_acc, 1L, mode = "normalized"), tolerance = 1e-10)
})

test_that("chunked backend matches blas backend", {
  set.seed(556)
  Tn <- 80L; V <- 5L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  fit_blas <- mppi_fit(Y, X, psych_idx = 2L, backend = "blas")
  fit_chunk <- mppi_fit(Y, X, psych_idx = 2L, backend = "chunked", chunk_size = 16L)
  expect_equal(mppi_get_M_scaled(fit_blas, 1L, mode = "normalized"),
               mppi_get_M_scaled(fit_chunk, 1L, mode = "normalized"), tolerance = 1e-8)
})

test_that("chunked basis projection matches blas projection", {
  set.seed(557)
  Tn <- 50L; V <- 6L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  basis <- as_mppi_basis(qr.Q(qr(matrix(rnorm(V * 4L), V, 4L))), name = "randPCA")
  fit_blas <- mppi_fit(Y, X, psych_idx = 2L, basis = basis,
                       backend = "blas", project_backend = "blas")
  fit_chunk <- mppi_fit(Y, X, psych_idx = 2L, basis = basis,
                        backend = "chunked", chunk_size = 16L,
                        project_backend = "chunked", project_chunk_cols = 2L)
  expect_equal(mppi_get_M_scaled(fit_blas, 1L, mode = "normalized"),
               mppi_get_M_scaled(fit_chunk, 1L, mode = "normalized"), tolerance = 1e-8)
})

test_that("variance decomposition and partial outputs present for cov scale", {
  set.seed(42)
  Tn <- 40L; V <- 3L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  fit <- mppi_fit(Y, X, psych_idx = 2L, scale = "cov")
  expect_true(!is.null(fit$Sigma0))
  expect_true(!is.null(fit$variance))
  expect_true(!is.null(fit$partial))
  expect_true("DeltaRho" %in% names(fit$partial[[1]]))
})

test_that("packed option stores upper triangle", {
  set.seed(77)
  Tn <- 30L; V <- 3L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  fit_packed <- mppi_fit(Y, X, psych_idx = 2L, packed = TRUE)
  expect_true(fit_packed$packed)
  packed_entry <- fit_packed$Delta[[1]]
  expect_type(packed_entry, "list")
  expect_true(all(c("values","dim") %in% names(packed_entry)))
  fit_full <- mppi_fit(Y, X, psych_idx = 2L, packed = FALSE)
  unpacked <- multiPPI:::`.mppi_unpack_upper`(packed_entry)
  expect_equal(unpacked, mppi_get_M_scaled(fit_full, 1L, mode = "normalized"), tolerance = 1e-10)
  # round trip pack -> unpack -> pack
  repacked <- multiPPI:::`.mppi_pack_upper`(unpacked)
  expect_equal(repacked$values, packed_entry$values)
  expect_equal(repacked$dim, packed_entry$dim)
})

test_that("mppi matrix method returns mppi_fit", {
  set.seed(88)
  Tn <- 40L; V <- 3L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  pidx <- 2L
  fit <- mppi(Y, X, psych_idx = pidx)
  expect_s3_class(fit, "mppi_fit")
})

test_that("basis projection matches full-space Δ", {
  set.seed(101)
  Tn <- 60L; V <- 8L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  pidx <- 2L
  full <- mppi_fit(Y, X, psych_idx = pidx, backend = "blas")
  Vbasis <- as_mppi_basis(diag(V), name = "identity")
  fit_basis <- mppi_fit(Y, X, psych_idx = pidx, basis = Vbasis,
                        backend = "chunked", chunk_size = 16L,
                        project_backend = "chunked", project_chunk_cols = 4L)
  M <- mppi_get_M_scaled(fit_basis, 1L, mode = "normalized")
  Delta_full <- mppi_get_M_scaled(full, 1L, mode = "normalized")
  recon <- Vbasis$V %*% M %*% t(Vbasis$V)
  diag(recon) <- 0
  diag(Delta_full) <- 0
  expect_equal(recon, Delta_full, tolerance = 1e-6)
  submat <- mppi_reconstruct_delta(fit_basis, 1L, rows = 1:3, cols = 5:7)
  expect_equal(submat, recon[1:3, 5:7], tolerance = 1e-6)
})

test_that("Precision update matches first-order matrix inverse", {
  Sigma0 <- matrix(c(2, 0.3, 0.3, 1.5), 2, 2)
  DeltaSigma <- matrix(c(0.2, -0.1, -0.1, 0.3), 2, 2)
  lambda <- 1e-3
  Theta0 <- solve(Sigma0 + lambda * diag(2))

  DeltaTheta <- mppi_to_partial(Sigma0, DeltaSigma, lambda = lambda)

  eps <- 1e-4
  Theta_exact <- solve(Sigma0 + eps * DeltaSigma + lambda * diag(2))
  Theta_linear <- Theta0 + eps * DeltaTheta
  expect_equal(Theta_linear, Theta_exact, tolerance = 1e-6)
})

test_that("Variance decomposition reconstructs DeltaSigma", {
  set.seed(99)
  Tn <- 80L; V <- 3L
  pk <- as.numeric(scale(rnorm(Tn), scale = FALSE))
  base <- matrix(rnorm(Tn * V), Tn, V)
  R <- sweep(base, 1, 1 + 0.3 * pk, `*`)
  DeltaSigma <- multiPPI:::`.mppi_wcp`(R, pk) / sum(pk^2)

  parts <- mppi_decompose_variance(R, pk, DeltaSigma)
  recon <- (apply(R, 2, sd) %o% apply(R, 2, sd)) * parts$Delta_rho + parts$variance_terms
  expect_equal(recon, DeltaSigma, tolerance = 1e-8)
})

test_that("Rank CV prefers low-rank structure", {
  set.seed(2024)
  Tn <- 70L; V <- 5L
  tvec <- seq_len(Tn) / Tn
  coeff <- runif(V, 0.5, 1.5)
  R <- as.matrix(tvec %o% coeff)
  pk <- rep(1, Tn)
  Delta <- multiPPI:::`.mppi_wcp`(R, pk) / sum(pk^2)

  cv <- mppi_rank_cv(R, pk, Delta, holdout_frac = 0.25, R_reps = 3)
  expect_equal(cv$r, 1L)
})
