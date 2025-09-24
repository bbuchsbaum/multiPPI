test_that("mppi_rest_modulators produces aligned sticks", {
  set.seed(1)
  U <- matrix(rnorm(300), 100, 3)
  mods <- mppi_rest_modulators(U,
                               include = c("states", "bursts", "envelopes"),
                               states = list(K = 3L, method = "kmeans", soft = FALSE),
                               bursts = list(ar_order = 2L, thresholds = 0.9),
                               envelopes = list(transform = "power"))
  expect_gt(length(mods), 0)
  expect_true(all(vapply(mods, length, 1L) == nrow(U)))
  expect_true(all(vapply(mods, function(v) all(is.finite(v)), TRUE)))
})

test_that("mppi_rest returns an mppi_fit with attached modulators", {
  set.seed(2)
  Y <- matrix(rnorm(2000), 100, 20)
  fit <- mppi_rest(Y,
                   basis = list(type = "pca", r = 6L),
                   include = c("states", "bursts", "envelopes"),
                   prewhiten = FALSE,
                   modulator_domain = "bold",
                   domain = "bold")
  expect_s3_class(fit, "mppi_fit")
  pk <- attr(fit, "pk")
  expect_true(is.list(pk))
  expect_gt(length(pk), 0)
  expect_true(all(vapply(pk, length, 1L) == nrow(Y)))
})

test_that("mppi_rest_features aggregates by subject", {
  skip_if_not_installed("fmriAR")
  set.seed(3)
  Y <- matrix(rnorm(600), 100, 6)
  runs <- rep(1:2, each = 50)
  fits <- list(
    mppi_rest(Y, runs = runs, basis = list(type = "pca", r = 4),
              prewhiten = TRUE, include = c("states", "bursts"),
              modulator_domain = "bold", domain = "bold"),
    mppi_rest(Y, runs = runs, basis = list(type = "pca", r = 4),
              prewhiten = TRUE, include = c("states", "envelopes"),
              modulator_domain = "bold", domain = "bold")
  )
  feats <- mppi_rest_features(fits, aggregate = "mean")
  expect_equal(nrow(feats), length(fits))
  expect_true(all(c("G", "R", "mag") %in% colnames(feats)))
})
