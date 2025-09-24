library(testthat)

set.seed(42)

Tn <- 40
V <- 4
psych <- as.numeric(scale(rnorm(Tn)))
conf <- as.numeric(scale(rnorm(Tn)))
base_noise <- matrix(rnorm(Tn * V, sd = 0.3), Tn, V)
effect <- matrix(c(0.8, -0.6, 0.4, 0), nrow = Tn, ncol = V, byrow = TRUE)
Y <- base_noise + psych * effect + 0.5 * conf * matrix(c(0.2, 0, -0.1, 0.3), nrow = Tn, ncol = V, byrow = TRUE)
X <- cbind(`(Intercept)` = 1, conf, psych)

fit_bold <- mppi_fit(Y, X, psych_idx = 3L, zero_diag = FALSE, scale = "cov")

fit_innov <- mppi_fit(Y, X, psych_idx = 3L,
                      zero_diag = FALSE,
                      scale = "cov",
                      domain = "innovations",
                      deconv = list(lambda = 5,
                                    protect = "none",
                                    ar_order = "fixed",
                                    p = 2L,
                                    tr = 1))


test_that("innovations domain stores neural and innovation residuals", {
  expect_identical(fit_innov$domain, "innovations")
  expect_false(is.null(fit_innov$U))
  expect_true(is.matrix(fit_innov$R_raw))
  expect_true(is.matrix(fit_innov$R))
  expect_false(isTRUE(all.equal(fit_innov$R, fit_innov$R_raw)))
  expect_true("innov" %in% names(fit_innov$deconv))
  expect_equal(nrow(fit_innov$R), nrow(fit_innov$R_raw))
})

test_that("mppi_evidence returns per-column and aggregate summaries", {
  ev <- mppi_evidence(fit_bold)
  expect_true(is.list(ev))
  expect_true(all(c("per_column", "aggregate") %in% names(ev)))
  per_col <- ev$per_column
  agg <- ev$aggregate
  expect_s3_class(per_col, "data.frame")
  expect_true(all(c("component", "delta_AIC", "delta_BIC") %in% names(per_col)))
  expect_equal(nrow(per_col), ncol(Y))
  expect_s3_class(agg, "data.frame")
  expect_equal(nrow(agg), 1L)
  expect_gt(agg$delta_AIC, 0)
  expect_gt(agg$delta_BIC, 0)
  ev_noagg <- mppi_evidence(fit_bold, aggregate = FALSE)
  expect_s3_class(ev_noagg, "data.frame")
  expect_equal(nrow(ev_noagg), ncol(Y))
})

test_that("evidence totals align with stored metadata", {
  ev <- mppi_evidence(fit_bold)
  stored <- fit_bold$evidence
  expect_equal(sum(ev$per_column$delta_AIC), stored$delta_AIC_total, tolerance = 1e-8)
  expect_equal(sum(ev$per_column$delta_BIC), stored$delta_BIC_total, tolerance = 1e-8)
})
