test_that("rank CV scores the off-diagonal structure that zero-diagonal fits analyse", {
  set.seed(61); Tn <- 200; V <- 10
  pk <- rep(c(-1, 1), each = Tn / 2)
  R <- matrix(rnorm(Tn * V), Tn) * ifelse(pk > 0, 2, 1)
  expect_equal(mppi_rank_cv(R, pk, R_reps = 3)$r, 0)
  expect_gt(mppi_rank_cv(R, pk, R_reps = 3, zero_diag = FALSE)$r, 0)
})

test_that("DCT deconvolution rejects prior weights of the wrong length", {
  set.seed(62); Y <- matrix(rnorm(100), 50)
  expect_error(multiPPI:::mppi_deconv_dct(Y, c(0, 1, 0.5), K = 10, method = "fixed", q = rep(1, 5)), "length K")
  expect_silent(multiPPI:::mppi_deconv_dct(Y, c(0, 1, 0.5), K = 10, method = "fixed", q = rep(1, 10)))
})
