context("Group-level wrappers")

set.seed(123)

build_fit <- function(seed_shift = 0) {
  set.seed(100 + seed_shift)
  Tn <- 16L; V <- 4L
  X <- cbind(1, seq_len(Tn))
  Y <- matrix(rnorm(Tn * V), Tn, V)
  mppi_fit(Y = Y, X = X, psych_idx = 2, lags = -1:1)
}

fits <- lapply(0:2, build_fit)
names(fits) <- paste0("subj", seq_along(fits))

stack_obj <- mppi_stack(fits, k = 1, lag = 0)

test_that("mppi_stack produces a complete edge matrix", {
  expect_equal(nrow(stack_obj$matrix), length(fits))
  expect_equal(ncol(stack_obj$matrix), choose(4, 2))
  expect_true(all(is.finite(stack_obj$matrix)))
})

glm_design <- cbind(Intercept = 1, group = c(-1, 0, 1))

res_glm <- mppi_group_glm(stack_obj, design = glm_design, method = "lm")

test_that("mppi_group_glm runs OLS without limma", {
  expect_equal(nrow(res_glm), ncol(stack_obj$matrix))
  expect_true(all(is.finite(res_glm$estimate)))
  expect_true(all(!is.na(res_glm$q)))
})

test_that("lagged accessors work", {
  M0 <- mppi_get_M(fits[[1]], 1)
  M1 <- mppi_get_M_lag(fits[[1]], 1, lag = 1)
  Mneg <- mppi_get_M_lag(fits[[1]], 1, lag = -1)
  expect_equal(dim(M0), dim(M1))
  expect_equal(dim(M0), dim(Mneg))
  expect_true(is.matrix(M1))
  expect_true(is.matrix(Mneg))
})
