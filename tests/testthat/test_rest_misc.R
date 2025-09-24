test_that("mppi_rest_design returns normalized columns", {
  set.seed(11)
  Tn <- 200
  base <- sin(seq_len(Tn) / 10)
  pk <- list(
    ctx1 = base,
    ctx2 = base + rnorm(Tn, sd = 0.1)
  )
  design <- mppi_rest_design(pk, nuisance = matrix(rnorm(Tn), Tn, 1))
  Xc <- design$X[, design$psych_idx, drop = FALSE]
  expect_equal(ncol(Xc), length(pk))
  norms <- sqrt(colSums(Xc^2))
  expect_true(all(abs(norms - 1) < 1e-6))
  expect_true(all(design$X[, 1] == 1))
})

test_that("mppi_rest_trait_cv recovers linear trait relationship", {
  set.seed(12)
  subj <- 1:12
  trait <- scale(rnorm(length(subj)))[, 1]
  feats <- do.call(rbind, lapply(subj, function(s) {
    data.frame(
      condition = c("ctx1", "ctx2"),
      G = c(trait[s] + rnorm(1, sd = 0.01), rnorm(1, sd = 0.3)),
      R = rnorm(2, sd = 0.2),
      mag = rnorm(2, sd = 0.2),
      subj = s
    )
  }))
  cv <- mppi_rest_trait_cv(feats, y = setNames(trait, subj), k = 4)
  expect_true(is.list(cv))
  expect_gt(cv$cv_r2, 0.6)
  expect_equal(length(cv$pred), length(trait))
})

test_that("mppi_sim_rest_dataset respects custom templates", {
  set.seed(13)
  r <- 4
  template <- matrix(rnorm(r * r), r, r)
  template <- 0.5 * (template + t(template))
  sims <- mppi_sim_rest_dataset(T = 200, V = 20, r = r,
                                templates = list(template),
                                contexts = "ctx", amp = 0.5,
                                obs_noise = 0)
  truth_mat <- sims$truth$coupling[["ctx"]]
  fnorm <- sqrt(sum((0.5 * (template + t(template)))^2))
  expect_equal(truth_mat, (0.5 * (template + t(template))) / fnorm)
})

test_that("mppi_rest_modulators handles envelope groups", {
  set.seed(14)
  U <- matrix(rnorm(300), 100, 3)
  mods <- mppi_rest_modulators(U,
                               include = "envelopes",
                               envelopes = list(groups = list(a = 1:2, b = 3)))
  expect_named(mods, c("env_a", "env_b"))
  expect_equal(length(mods$env_a), nrow(U))
})

test_that("mppi_rest_extract handles matrices with attributes", {
  Y <- matrix(rnorm(100), 20, 5)
  attr(Y, "TR") <- 0.8
  attr(Y, "runs") <- rep(1:2, each = 10)
  info <- mppi_rest_extract(Y)
  expect_equal(info$TR, 0.8)
  expect_equal(nrow(info$Y), 20)
  expect_equal(info$runs, as.integer(rep(1:2, each = 10)))
})
