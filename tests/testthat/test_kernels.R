context("Kernel helpers")

set.seed(2025)

test_that("ppi_batch matches manual zero-lag crossproducts", {
  Tn <- 50L; r <- 4L; K <- 3L
  U <- matrix(rnorm(Tn * r), Tn, r)
  pk_mat <- matrix(rnorm(Tn * K), Tn, K)

  res <- ppi_batch(U, pk_mat, lags = 0L)
  mats <- res$matrices[[1]]
  denom <- res$denominator[, 1]

  for (k in seq_len(K)) {
    pk <- pk_mat[, k]
    d <- sum(pk^2)
    if (d < .Machine$double.eps) next
    manual <- crossprod(U, U * pk) / d
    expect_equal(mats[,, k], manual, tolerance = 1e-10)
    expect_equal(denom[k], d, tolerance = 1e-10)
  }
})

test_that("ppi_batch handles lags and blocklens", {
  Tn <- 40L; r <- 3L; K <- 2L
  U <- matrix(rnorm(Tn * r), Tn, r)
  pk_mat <- matrix(rnorm(Tn * K), Tn, K)
  lags <- c(-2L, 0L, 1L)
  runs <- rep(1:2, each = Tn/2)
  blocklens <- rep(Tn/2, 2)

  res <- ppi_batch(U, pk_mat, lags = lags, blocklens = blocklens)

  for (li in seq_along(lags)) {
    lg <- lags[li]
    cube <- res$matrices[[li]]
    for (k in seq_len(K)) {
      pk <- pk_mat[, k]
      if (lg == 0L) {
        denom <- sum(pk^2)
        Mk_manual <- multiPPI:::`.mppi_crossprod`(U, pk) / denom
      } else {
        Mk_manual <- multiPPI:::`.mppi_crossprod_lagged`(U, pk, lg, blocklens = blocklens)
      }
      expect_equal(cube[,, k], Mk_manual, tolerance = 1e-8)
    }
  }
})

test_that("precision_gate_cpp matches manual", {
  r <- 3L; K <- 2L
  Theta0 <- diag(runif(r, 0.5, 1.5))
  Ms <- array(0, dim = c(r, r, K))
  Ms[,,1] <- matrix(c(0.2, -0.1, 0.05,
                      -0.1, 0.3, 0.07,
                      0.05, 0.07, -0.2), r, r)
  Ms[,,2] <- matrix(c(-0.1, 0.04, 0,
                       0.04, 0.05, -0.02,
                       0, -0.02, 0.03), r, r)
  gate <- precision_gate_cpp(Theta0, Ms)
  for (k in seq_len(K)) {
    expect_equal(gate[,, k], -Theta0 %*% Ms[,, k] %*% Theta0, tolerance = 1e-10)
  }
})
