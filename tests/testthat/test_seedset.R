test_that("seed-set utilities support multiple selectors", {
  set.seed(7)
  Tn <- 36L
  V <- 5L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  colnames(Y) <- paste0("v", seq_len(V))
  X <- cbind(1, sin(seq_len(Tn) / 4))
  fit <- mppi_fit(Y, X, psych_idx = 2L, scale = "normalized", zero_diag = FALSE)

  slopes_num <- mppi_seedset_slopes(fit, seeds = 1:2, k = 1L)
  slopes_char <- mppi_seedset_slopes(fit, seeds = c("v1", "v2"), k = 1L)
  expect_equal(slopes_num, slopes_char)

  slopes_logical <- mppi_seedset_slopes(fit, seeds = c(TRUE, TRUE, FALSE, FALSE, FALSE), k = 1L)
  expect_equal(slopes_num, slopes_logical)

  expect_error(mppi_seedset_slopes(fit, seeds = "unknown", k = 1L), "Unknown")

  map_mean <- mppi_seedset_map(fit, seeds = 1:2, k = 1L, aggregate = "mean")
  expect_s3_class(map_mean, "mppi_seedset_map")
  expect_equal(map_mean$aggregate, "mean")
  expect_equal(length(map_mean$map), length(map_mean$targets))
})

test_that("seed-set map handles bases and precision aggregation", {
  set.seed(17)
  Tn <- 32L
  V <- 6L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  colnames(Y) <- paste0("roi", seq_len(V))
  X <- cbind(1, cos(seq_len(Tn) / 5))
  Vbasis <- qr.Q(qr(matrix(rnorm(V * 3), V, 3)))
  rownames(Vbasis) <- paste0("roi", seq_len(V))
  basis <- list(V = Vbasis, name = "rand", r = ncol(Vbasis))
  fit <- mppi_fit(Y, X, psych_idx = 2L, basis = basis, scale = "cov", zero_diag = FALSE)

  Theta0 <- diag(ncol(fit$R)) * 0.5
  map_prec <- mppi_seedset_map(fit, seeds = 1:2, k = 1L,
                               aggregate = "precision", Theta0 = Theta0,
                               mode = "amplitude")
  expect_equal(map_prec$aggregate, "precision")
  expect_equal(map_prec$mode, "amplitude")
  expect_true(all(names(map_prec$seeds) %in% c("roi1", "roi2")))
})

test_that("seed-set permutation produces reproducible metadata", {
  set.seed(101)
  Tn <- 30L
  V <- 4L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind(1, seq_len(Tn) / Tn)
  fit <- mppi_fit(Y, X, psych_idx = 2L, scale = "normalized", zero_diag = FALSE)

  res <- mppi_seedset_permute(fit, seeds = 1:2, k = 1L, aggregate = "mean",
                              mode = "normalized", blksize = 5L, B = 29L,
                              seed = 42, keep_draws = TRUE)
  expect_true(all(c("map", "weights", "statistic", "settings") %in% names(res)))
  expect_equal(res$aggregate, "mean")
  expect_equal(length(res$map), length(res$targets))
  expect_equal(length(res$statistic$p), length(res$map))
  expect_equal(nrow(res$draws), 29L)
  rng_meta <- attr(res, "rng")
  expect_equal(rng_meta$seed_arg, 42)
})
