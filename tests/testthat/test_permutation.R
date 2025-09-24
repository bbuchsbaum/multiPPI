test_that("mppi_omnibus exposes RNG metadata", {
  set.seed(999)
  Tn <- 40L; V <- 3L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  fit <- mppi_fit(Y, X, psych_idx = 2L)
  res <- mppi_omnibus(fit, blksize = 5L, B = 19L, seed = 123)
  meta <- attr(res, "rng")
  expect_type(meta, "list")
  expect_equal(meta$seed_arg, 123)
  expect_true(length(meta$kind) >= 2)
})

test_that("mppi_permute phase method returns metadata", {
  set.seed(2025)
  Tn <- 30L; V <- 3L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = sin(seq_len(Tn)/4))
  fit <- mppi_fit(Y, X, psych_idx = 2L, runs = rep(1:3, each = 10))
  pm <- mppi_permute(fit, k = 1L, blksize = 5L, B = 11L, method = "phase", seed = 321)
  meta <- attr(pm, "rng")
  expect_type(meta, "list")
  expect_equal(meta$method, "phase")
})

test_that("freedman-lane omnibus runs and returns metadata", {
  set.seed(1234)
  Tn <- 24L; V <- 3L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  fit <- mppi_fit(Y, X, psych_idx = 2L)
  res <- mppi_omnibus(fit, blksize = 4L, B = 19L, method = "freedman_lane", seed = 42)
  meta <- attr(res, "rng")
  expect_equal(meta$method, "freedman_lane")
  expect_equal(length(res[[1]]$Qnull), 19L)
})

test_that("freedman-lane edgewise permutation yields matrix", {
  set.seed(567)
  Tn <- 20L; V <- 3L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(Intercept = 1, task = rnorm(Tn))
  fit <- mppi_fit(Y, X, psych_idx = 2L)
  pm <- mppi_permute(fit, k = 1L, blksize = 4L, B = 9L, method = "freedman_lane", seed = 99,
                     studentize = FALSE)
  expect_equal(dim(pm), c(V, V))
  expect_true(all(is.na(diag(pm))))
  meta <- attr(pm, "rng")
  expect_equal(meta$method, "freedman_lane")
})
