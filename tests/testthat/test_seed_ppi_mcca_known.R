test_that("seed_ppi_mcca returns 0 on zero interaction block", {
  set.seed(42)
  Tn <- 150
  V  <- 40
  # Design with one centered column
  D  <- scale(matrix(rnorm(Tn), Tn, 1), center = TRUE, scale = FALSE)
  # Data matrix: first 12 columns are zero (neighborhood), others random
  Y  <- matrix(rnorm(Tn * V), Tn, V)
  Y[, 1:12] <- 0
  # Seed as exogenous series
  Xseed <- scale(matrix(rnorm(Tn), Tn, 1), center = TRUE, scale = FALSE)
  # Neighborhood that covers the zero block
  neighborhoods <- list(as.integer(1:12))

  res <- seed_ppi_mcca(
    data = Y,
    seed = Xseed,
    design = D,
    neighborhoods = neighborhoods,
    tr = 1,
    domain = "hemo",
    prewhiten = FALSE,
    pca = list(seed_k = 1, light_k = 5, center = TRUE, scale = TRUE),
    reg = list(l2_x = 0.05, l2_y = 0.05),
    stat = "r1",
    inference = list(method = "none"),
    return = "map"
  )
  expect_s3_class(res, "mppi_sl")
  expect_equal(length(res$stat), 1L)
  expect_equal(res$stat[1], 0, tolerance = 1e-10)
})

test_that("seed_ppi_mcca ~ null shows small effect for independent noise", {
  set.seed(7)
  Tn <- 800
  V  <- 80
  # Independent Gaussian seed, design, and data
  Y  <- matrix(rnorm(Tn * V), Tn, V)
  D  <- scale(cbind(rnorm(Tn), rnorm(Tn)), center = TRUE, scale = FALSE)
  Xseed <- scale(matrix(rnorm(Tn), Tn, 1), center = TRUE, scale = FALSE)
  neighborhoods <- list(as.integer(1:16))

  res <- seed_ppi_mcca(
    data = Y,
    seed = Xseed,
    design = D,
    neighborhoods = neighborhoods,
    tr = 1,
    domain = "hemo",
    prewhiten = FALSE,
    pca = list(seed_k = 1, light_k = 8, center = TRUE, scale = TRUE),
    reg = list(l2_x = 0.05, l2_y = 0.05),
    stat = "r1",
    inference = list(method = "none"),
    return = "map"
  )
  expect_true(is.finite(res$stat[1]))
  expect_lt(abs(res$stat[1]), 0.15)
})

test_that("Wilks with k=1 <= full Wilks and matches -log(1-r1^2)", {
  set.seed(101)
  Tn <- 300
  V  <- 60
  Y  <- matrix(rnorm(Tn * V), Tn, V)
  D  <- scale(cbind(rnorm(Tn), rnorm(Tn)), center = TRUE, scale = FALSE)
  # Seed with q=4 components
  Xseed <- scale(matrix(rnorm(Tn * 4), Tn, 4), center = TRUE, scale = FALSE)
  neighborhoods <- list(as.integer(1:24))

  # Full Wilks
  res_full <- seed_ppi_mcca(
    data = Y,
    seed = Xseed,
    design = D,
    neighborhoods = neighborhoods,
    tr = 1,
    domain = "hemo",
    prewhiten = FALSE,
    pca = list(seed_k = 4, light_k = 8, center = TRUE, scale = TRUE),
    reg = list(l2_x = 0.05, l2_y = 0.05),
    stat = "wilks",
    wilks_k = NULL,
    inference = list(method = "none"),
    return = "map"
  )
  # Wilks with k=1 (top root only)
  res_k1 <- seed_ppi_mcca(
    data = Y,
    seed = Xseed,
    design = D,
    neighborhoods = neighborhoods,
    tr = 1,
    domain = "hemo",
    prewhiten = FALSE,
    pca = list(seed_k = 4, light_k = 8, center = TRUE, scale = TRUE),
    reg = list(l2_x = 0.05, l2_y = 0.05),
    stat = "wilks",
    wilks_k = 1L,
    inference = list(method = "none"),
    return = "map"
  )
  # r1-only expected Wilks equals -log(1 - r1^2)
  res_r1 <- seed_ppi_mcca(
    data = Y,
    seed = Xseed,
    design = D,
    neighborhoods = neighborhoods,
    tr = 1,
    domain = "hemo",
    prewhiten = FALSE,
    pca = list(seed_k = 4, light_k = 8, center = TRUE, scale = TRUE),
    reg = list(l2_x = 0.05, l2_y = 0.05),
    stat = "r1",
    inference = list(method = "none"),
    return = "map"
  )
  expect_lte(res_k1$stat[1], res_full$stat[1] + 1e-10)
  expect_equal(res_k1$stat[1], -log(1 - pmin(res_r1$stat[1], 0.999)^2), tolerance = 1e-8)
})

