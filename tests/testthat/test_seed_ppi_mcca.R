test_that("seed_ppi_mcca runs on matrix with provided neighborhoods", {
  skip_on_os("windows")
  set.seed(123)
  Tn <- 80
  V <- 60
  # synthetic data: two task regressors
  X <- cbind(1, rnorm(Tn), as.numeric(gl(4, Tn/4)))
  colnames(X) <- c("int","cond","block")
  Y <- matrix(rnorm(Tn * V), Tn, V)
  # seed as indices
  seed_idx <- 1:5
  # make simple neighborhoods of size 10
  neighborhoods <- lapply(seq(1, V, by = 3), function(i) as.integer(intersect(i:(i+9), 1:V)))
  # run with no prewhitening and hemo domain to avoid suggests
  res <- multiPPI:::seed_ppi_mcca(
    data = Y,
    seed = seed_idx,
    design = as.matrix(X),
    nuisance = NULL,
    neighborhoods = neighborhoods,
    tr = 1,
    domain = "hemo",
    prewhiten = FALSE,
    pca = list(seed_k = 3, light_k = 5, center = TRUE, scale = TRUE),
    reg = list(l2_x = 0.1, l2_y = 0.1),
    stat = "r1",
    inference = list(method = "none", perms = 0L),
    return = c("map","weights")
  )
  expect_s3_class(res, "mppi_sl")
  expect_equal(length(res$stat), length(neighborhoods))
  expect_true(any(is.finite(res$stat)))
  expect_length(res$weights_x, length(neighborhoods))
  expect_length(res$weights_y, length(neighborhoods))
})

