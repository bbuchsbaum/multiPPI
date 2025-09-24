context("Neural-domain defaults and grouping")

test_that("neural fit works with default HRF", {
  set.seed(101)
  Tn <- 48L; V <- 6L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind("(Intercept)" = 1, conf = rnorm(Tn), task = rnorm(Tn))
  fit <- mppi_fit(Y, X, psych_idx = 3L, domain = "neural")
  expect_equal(fit$domain, "neural")
  expect_equal(length(fit$names), 1L)
  expect_type(fit$deconv$sticks, "double")
  expect_equal(ncol(fit$deconv$sticks), 1L)
  expect_equal(fit$deconv$hrf, mppi_default_hrf())
  expect_equal(length(fit$pk), 1L)
  expect_false(any(is.na(fit$Delta[[1]])))
})

test_that("neural fit groups HRF basis columns", {
  set.seed(102)
  Tn <- 64L; V <- 5L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind("(Intercept)" = 1,
             nuisance = rnorm(Tn),
             hrf_A = rnorm(Tn),
             hrfTD_A = rnorm(Tn),
             hrf_B = rnorm(Tn),
             hrfTD_B = rnorm(Tn))
  fit <- mppi_fit(Y, X, psych_idx = 3:6, domain = "neural")
  expect_equal(length(fit$names), 2L)
  expect_true(all(grepl("^s:", fit$names)))
  expect_equal(ncol(fit$deconv$sticks), 2L)
  expect_equal(length(fit$pk), 2L)
  expect_false(anyNA(unlist(fit$Delta)))
})

test_that("neural fit supports MAP deconvolution backend", {
  set.seed(103)
  Tn <- 40L; V <- 4L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind("(Intercept)" = 1,
             conf = rnorm(Tn),
             taskA = rnorm(Tn),
             taskB = rnorm(Tn))
  fit_map <- mppi_fit(Y, X, psych_idx = 3:4, domain = "neural",
                      deconv = list(type = "map", lambda = 8, TR = 1))
  expect_equal(fit_map$deconv$type, "map")
  expect_equal(ncol(fit_map$deconv$sticks), length(fit_map$names))
  expect_true(all(vapply(fit_map$pk, function(v) length(v) == Tn, logical(1))))
  expect_false(any(is.na(fit_map$Delta[[1]])))
})
