test_that("mppi_select_psych handles names, patterns, and errors", {
  set.seed(1L)
  X <- cbind(`(Intercept)` = 1,
             taskA = rnorm(12),
             taskB = rnorm(12),
             confound = rnorm(12))
  expect_equal(mppi_select_psych(X, names_or_idx = "taskA"), 2L)
  expect_equal(mppi_select_psych(X, names_or_idx = c(3L, 2L)), c(3L, 2L))
  expect_equal(sort(mppi_select_psych(X, patterns = "^task")), 2:3)
  expect_error(mppi_select_psych(X), "Provide either names_or_idx or patterns")
  expect_error(mppi_select_psych(X, patterns = "^missing"), "No design columns matched")
})

test_that("mppi_coerce_beta coerces frame and list inputs", {
  mat <- matrix(1:6, nrow = 3, byrow = TRUE)
  expect_identical(mppi_coerce_beta(mat), mat)

  df <- data.frame(trial = c(2, 1, 2, 1, 2, 1),
                   roi = rep(letters[1:3], each = 2),
                   beta = 1:6)
  coerced_df <- mppi_coerce_beta(df)
  expect_equal(dim(coerced_df), c(2, 3))
  expected_trial1 <- stats::setNames(df$beta[df$trial == 1], df$roi[df$trial == 1])
  expect_equal(coerced_df[1, ], expected_trial1[sort(names(expected_trial1))])

  lst <- list(c(1, 2, 3), c(4, 5, 6))
  coerced_list <- mppi_coerce_beta(lst)
  expect_equal(coerced_list, rbind(lst[[1]], lst[[2]]))

  bad_list <- list(1:3, 1:2)
  expect_error(mppi_coerce_beta(bad_list), "same length")
})

test_that("mppi_compare_models returns delta gaps and AIC differences", {
  fit_bold <- list(names = c("c1", "c2"),
                   Delta = list(matrix(c(1, 0, 0, 1), 2),
                                matrix(c(0, 1, 1, 2), 2)))
  fit_neural <- list(names = c("c1", "c2"),
                     Delta = list(matrix(c(0.5, 0, 0, 0.5), 2),
                                  matrix(c(0.5, 0.5, 0.5, 0.5), 2)))
  class(fit_bold) <- class(fit_neural) <- c("mppi_fit", "list")

  resid_bold <- matrix(c(1, 2, 3,
                         4, 5, 6), nrow = 3, byrow = TRUE)
  resid_neural <- resid_bold * 0.5
  pk <- list(c(1, 0, 1), c(0, 1, 1))

  res <- mppi_compare_models(fit_bold, fit_neural,
                             resid_bold = resid_bold,
                             resid_neural = resid_neural,
                             pk = pk)
  expect_equal(res$names, c("c1", "c2"))
  expect_equal(res$delta_gap, c(sqrt(0.5), 3), tolerance = 1e-7)
  expect_true(all(res$AIC_diff > 0))
})

test_that("mppi_hrf_ensemble blends fits with data-driven or manual weights", {
  fit1 <- list(names = "c1",
               Delta = list(matrix(c(1, 0, 0, 1), 2)),
               AIC_total = c(10),
               packed = FALSE)
  fit2 <- list(names = "c1",
               Delta = list(matrix(c(3, 0, 0, 3), 2)),
               AIC_total = c(12),
               packed = FALSE)
  class(fit1) <- class(fit2) <- c("mppi_fit", "list")

  ens_auto <- mppi_hrf_ensemble(list(fit1 = fit1, fit2 = fit2))
  weights_auto <- ens_auto$ensemble$weights
  expect_equal(sum(weights_auto), 1)
  expected_auto <- weights_auto[["fit1"]] * fit1$Delta[[1]] +
    weights_auto[["fit2"]] * fit2$Delta[[1]]
  expect_equal(ens_auto$Delta[[1]], expected_auto, tolerance = 1e-7)
  expect_equal(ens_auto$ensemble$method, "aic")

  ens_manual <- mppi_hrf_ensemble(list(fit1 = fit1, fit2 = fit2), weights = c(0.2, 0.8))
  expected_manual <- 0.2 * fit1$Delta[[1]] + 0.8 * fit2$Delta[[1]]
  expect_equal(ens_manual$Delta[[1]], expected_manual)
  expect_equal(ens_manual$ensemble$weights, c(fit1 = 0.2, fit2 = 0.8))
})

test_that("mppi_sim_rest_neural generates standardized modulators and templates", {
  sim <- mppi_sim_rest_neural(T = 60, r = 4,
                              contexts = c("ctxA", "ctxB"),
                              amp = c(0.2, 0.3),
                              freq = c(0.05, 0.08),
                              seed = 123)
  expect_equal(dim(sim$U), c(60, 4))
  expect_named(sim$modulators, c("ctxA", "ctxB"))
  mod_mat <- do.call(cbind, sim$modulators)
  expect_true(all(abs(colMeans(mod_mat)) < 1e-12))
  expect_true(all(abs(apply(mod_mat, 2, sd) - 1) < 1e-6))
  coupling_norms <- vapply(sim$coupling, function(M) sqrt(sum(M^2)), numeric(1))
  expect_true(all(abs(coupling_norms - 1) < 1e-8))
})
