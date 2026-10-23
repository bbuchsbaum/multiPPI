test_that("partial correlation derivative matches independent finite differences", {
  S <- matrix(c(2, .4, .4, 1.2), 2)
  D <- matrix(c(.5, -.2, -.2, .3), 2)
  lambda <- .01; theta <- solve(S + lambda * diag(2))
  dtheta <- mppi_to_partial(S, D, lambda)
  partial <- function(x) { q <- solve(x + lambda * diag(2)); z <- -cov2cor(q); diag(z) <- 1; z }
  eps <- 1e-5
  oracle <- (partial(S + eps * D) - partial(S - eps * D)) / (2 * eps)
  expect_equal(mppi_delta_partial(theta, dtheta), oracle, tolerance = 1e-8)
  expect_equal(diag(mppi_delta_partial(theta, dtheta)), c(0, 0))
})

test_that("variance and correlation derivatives distinguish pure changes", {
  R <- rbind(c(2, 1), c(0, sqrt(3)), c(-2, -1), c(0, -sqrt(3)))
  S <- crossprod(R) / nrow(R); v <- diag(S)
  variance_change <- c(.4, -.1)
  scaling <- diag(variance_change / (2 * v))
  Dvar <- scaling %*% S + S %*% scaling
  a <- mppi_decompose_variance(R, rep(1, 4), Dvar)
  expect_equal(a$Delta_rho, matrix(0, 2, 2), tolerance = 1e-12)
  expect_equal(a$variance_terms, Dvar, tolerance = 1e-12)
  Dcor <- matrix(c(0, .2, .2, 0), 2)
  b <- mppi_decompose_variance(R, rep(1, 4), Dcor)
  expect_equal(b$variance_terms, matrix(0, 2, 2))
  eps <- 1e-5
  oracle <- (cov2cor(S + eps * Dcor) - cov2cor(S - eps * Dcor)) / (2 * eps)
  expect_equal(b$Delta_rho, oracle, tolerance = 1e-9)
  expect_error(mppi_decompose_variance(matrix(0, 4, 2), rep(1, 4), Dcor), "positive")
})

test_that("rank CV retains negative modes and respects run-local gaps", {
  set.seed(123)
  R <- matrix(rep(c(1, 2), each = 24), 24, 2)
  cv <- mppi_rank_cv(R, rep(-1, 24), runs = rep(1:2, each = 12), gap = 1, R_reps = 2)
  expect_equal(cv$r, 1)
  expect_lt(cv$losses[2], 1e-20)
  for (s in cv$splits) for (run in 1:2) {
    ids <- which(rep(1:2, each = 12) == run)
    held <- intersect(s$test, ids); train <- intersect(s$train, ids)
    expect_true(all(diff(held) == 1))
    expect_true(all(abs(outer(train, held, "-")) > 1))
  }
})

test_that("group operations reject mismatched node and condition identities", {
  A <- matrix(c(1, .2, .2, 2), 2, dimnames = list(c("a", "b"), c("a", "b")))
  B <- A[2:1, 2:1]
  expect_error(mppi_group_ebayes(list(A, B), use_limma = FALSE), "node labels")
  expect_error(pwf_stack(list(A, B)), "node labels")
  expect_error(pwf_stack(list(A, A), labels = c("b", "a")), "match")
})

test_that("one-sample limma uses one observation per subject", {
  skip_if_not_installed("limma")
  set.seed(1); mats <- replicate(8, { a <- matrix(rnorm(9), 3); a + t(a) }, simplify = FALSE)
  actual <- mppi_group_ebayes(mats)
  X <- do.call(rbind, lapply(mats, function(x) x[upper.tri(x)]))
  ref <- limma::eBayes(limma::lmFit(t(X), matrix(1, 8, 1)))
  expect_equal(unname(actual$t), unname(ref$t[, 1]))
  expect_equal(unname(actual$p), unname(ref$p.value[, 1]))
})
