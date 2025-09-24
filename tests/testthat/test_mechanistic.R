context("Mechanistic helpers")

mock_fit <- function(pk_vals) {
  R <- diag(length(pk_vals))
  denom <- sum(pk_vals^2)
  Delta <- if (denom < .Machine$double.eps) matrix(0, length(pk_vals), length(pk_vals))
           else crossprod(R, pk_vals * R) / denom
  structure(list(Delta = list(Delta), names = "task", pk = list(pk_vals),
                 R = R, basis = NULL), class = c("mppi_fit", "list"))
}

test_that("mppi_gain matches manual precision update", {
  set.seed(1)
  Tn <- 32L; V <- 3L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  X <- cbind("(Intercept)" = 1, task = rnorm(Tn))
  fit <- mppi_fit(Y, X, psych_idx = 2L, zero_diag = FALSE, scale = "cov")
  pk <- fit$pk[[1]]
  R <- fit$R
  S0 <- crossprod(R) / nrow(R)
  Mk <- crossprod(R, pk * R) / sum(pk^2)
  Theta0 <- solve(S0 + diag(1e-3, ncol(S0)))
  dTheta <- - Theta0 %*% Mk %*% Theta0
  gain <- mppi_gain(fit, 1L, lambda = 1e-3)
  expect_equal(gain$gain, diag(dTheta), tolerance = 1e-8)
})

test_that("mppi_routing_index detects feedforward bias", {
  fit_pos <- mock_fit(c(2, 0))
  fit_neg <- mock_fit(c(0, 2))
  routing <- mppi_routing_index(list("1" = fit_pos, "-1" = fit_neg),
                                hierarchy = c(1, 2))
  expect_equal(routing$routing, (0.5 - 2) / (0.5 + 2), tolerance = 1e-8)
  expect_true(routing$energy_neg > routing$energy_pos)
})

test_that("mppi_project_templates returns normalized weights", {
  fit <- mock_fit(c(2, 0))
  templates <- list(diag = diag(2), off = matrix(c(0, 1, 1, 0), 2, 2))
  w <- mppi_project_templates(fit, templates = templates)
  expect_equal(names(w), names(templates))
  expect_equal(w[["diag"]], (0.5) / sqrt(2), tolerance = 1e-8)
  expect_equal(w[["off"]], 0, tolerance = 1e-12)
})

test_that("template helpers build orthogonal patterns", {
  labs <- c("A", "A", "B")
  tmp <- mppi_templates_within_between(labs)
  expect_equal(dim(tmp$within), c(3L, 3L))
  expect_true(all(diag(tmp$within) == 0))
  grad <- mppi_templates_gradient(c(1, 2, 3))
  expect_equal(dim(grad$dir), c(3L, 3L))
  expect_true(abs(sum(grad$dir * grad$ortho)) < 1e-8)
})

test_that("mechanism report bundles gain, routing, templates", {
  fit <- mock_fit(c(2, 0))
  lag_fits <- list("1" = fit, "-1" = mock_fit(c(0, 2)))
  templates <- list(diag = diag(2))
  rep <- mppi_mechanism_report(fit, hierarchy = c(1, 2), lag_fits = lag_fits,
                               templates = templates)
  expect_true(is.list(rep$gain))
  expect_true(is.list(rep$routing))
  expect_true(is.numeric(rep$templates))
})
