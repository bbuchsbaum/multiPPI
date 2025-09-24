context("Core mPPI shapes")

test_that("mppi_fit returns matrices of correct size", {
  set.seed(1)
  Tn <- 100; V <- 10
  Y <- matrix(rnorm(Tn*V), Tn, V)
  p <- rnorm(Tn)
  X <- cbind(1, p)
  fit <- mppi_fit(Y, X, psych_idx = 2L)
  expect_equal(length(fit$Delta), 1L)
  expect_equal(dim(fit$Delta[[1]]), c(V, V))
})
