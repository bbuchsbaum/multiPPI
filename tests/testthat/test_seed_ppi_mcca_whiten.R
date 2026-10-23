test_that("seed_ppi_mcca runs with prewhitening and run-centering", {
  skip_if_not_installed("fmriAR")
  set.seed(2025)
  Tn <- 120
  V  <- 50
  runs <- c(rep(1L, Tn/2), rep(2L, Tn/2))
  # Design: intercept + one task regressor (demeaned by run inside function)
  D  <- cbind(1, rnorm(Tn))
  colnames(D) <- c("int","task")
  # Data matrix: AR(1) structure to stabilize fmriAR estimation
  Y  <- matrix(0, Tn, V)
  phi <- 0.6
  noise <- matrix(rnorm(Tn * V), Tn, V)
  for (j in seq_len(V)) {
    for (t in 2:Tn) Y[t, j] <- phi * Y[t - 1, j] + noise[t, j]
  }
  # Seed as indices into columns of Y
  seed_idx <- 1:6
  # Simple neighborhoods (~12 voxels each)
  neighborhoods <- lapply(seq(1, V, by = 5), function(i) as.integer(intersect(i:(i+11), 1:V)))

  res <- seed_ppi_mcca(
    data = Y,
    seed = seed_idx,
    design = D,
    nuisance = NULL,
    neighborhoods = neighborhoods,
    tr = 1,
    runs = runs,
    center_by = "run",
    domain = "hemo",
    prewhiten = TRUE,
    ar_method = "ar",
    ar_order = 4L,
    pca = list(seed_k = 4, light_k = 8, center = TRUE, scale = TRUE),
    reg = list(l2_x = 0.1, l2_y = 0.1),
    stat = "r1",
    inference = list(method = "none"),
    return = "map"
  )

  expect_s3_class(res, "mppi_sl")
  expect_equal(length(res$stat), length(neighborhoods))
  expect_true(any(is.finite(res$stat)))
})
