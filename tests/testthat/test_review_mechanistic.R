context("Review fixes: mechanistic readouts on fitted objects")

# node 1 drives node 2 at lag 1 during task, or nodes 1 and 2 share a slow source
sim_lag_fit <- function(directed = TRUE, Tn = 300, V = 4, seed = 1, a = 2) {
  set.seed(seed)
  task <- rep(c(rep(1, 15), rep(0, 15)), length.out = Tn)
  z <- matrix(rnorm(Tn * V), Tn, V)
  if (directed) {
    for (t in 2:Tn) z[t, 2] <- z[t, 2] + a * task[t - 1] * z[t - 1, 1]
  } else {
    s <- as.numeric(arima.sim(list(ar = 0.8), n = Tn))
    z[, 1:2] <- z[, 1:2] + a * task * s
  }
  mppi_fit(z, cbind(int = 1, task = task), psych_idx = 2, lags = -2:2)
}

# three task-gated lag-1 drives (1->2, 3->4, 5->6), or a task-gated common gain
sim_route_fit <- function(seed, gain = FALSE, Tn = 400, a = 0.8, phi = 0.6) {
  set.seed(seed)
  task <- head(rep(rep(c(1, 0), 30), sample(8:20, 60, replace = TRUE)), Tn)
  z <- matrix(rnorm(Tn * 6), Tn, 6)
  if (gain) {
    z <- z + 1.5 * task * rnorm(Tn)
  } else {
    for (src in c(1, 3, 5)) {
      z[, src] <- as.numeric(arima.sim(list(ar = phi), n = Tn, sd = sqrt(1 - phi^2)))
      for (t in 2:Tn) if (task[t - 1] == 1) z[t, src + 1] <- sqrt(1 - a^2) * z[t, src + 1] + a * z[t - 1, src]
    }
  }
  mppi_fit(z, cbind(int = 1, task = task), psych_idx = 2, lags = -2:2)
}

test_that("routing indices recover a planted lagged drive and vanish for symmetric coupling", {
  fit <- sim_lag_fit(TRUE)
  expect_gt(mppi_get_M_lag(fit, 1, 1)[1, 2], 5 * abs(mppi_get_M_lag(fit, 1, 1)[2, 1]))
  up <- mppi_routing_index(fit, 1, hierarchy = 1:4)
  down <- mppi_routing_index(fit, 1, hierarchy = 4:1)
  expect_gt(up$routing, 0.5)
  expect_equal(down$routing, -up$routing)
  expect_gt(up$lead_energy[1, 2], 10 * up$lead_energy[2, 1])
  expect_equal(mppi_routing_index(list("1" = fit, "-1" = fit), 1, hierarchy = 1:4)$routing,
               mppi_routing_index(fit, 1, hierarchy = 1:4, pos = 1, neg = -1)$routing)
  expect_gt(mppi_lagged(fit, 1, lags = -2:2)$R, 0.8)
  expect_gt(mppi_axes(fit)$routing, 0.8)

  sym <- sim_lag_fit(FALSE)
  expect_lt(abs(mppi_routing_index(sym, 1, hierarchy = 1:4)$routing), 0.15)
  expect_lt(mppi_lagged(sym, 1, lags = -2:2)$R, 0.15)
  expect_lt(mppi_axes(sym)$routing, 0.15)

  feats <- mppi_rest_features(list(fit, sym), lags = -2:2, modes = 2L)
  expect_gt(feats$R[1], 0.8)
  expect_lt(feats$R[2], 0.15)
})

test_that("mppi_classify can label routing-dominant lagged coupling", {
  route <- mppi_classify(sim_route_fit(1), 1, lags = -2:2, B = 100, seed = 1)
  expect_lte(route$G, 0.5)
  expect_lte(route$pR, 0.05)
  expect_true(is.na(route$pI))
  expect_equal(route$label, "Routing-dominant")
  gain <- mppi_classify(sim_route_fit(1, gain = TRUE), 1, lags = -2:2, B = 100, seed = 1)
  expect_equal(gain$label, "Gain-dominant")
})

test_that("seed-set svd map keeps the sign of a common positive effect", {
  set.seed(11)
  Tn <- 400; V <- 12
  task <- rep(c(1, 0), each = 20, length.out = Tn)
  X <- cbind(int = 1, task = task)
  fits <- replicate(12, {
    z <- matrix(rnorm(Tn * V), Tn, V); drv <- rnorm(Tn); on <- task == 1
    z[on, 1:3] <- z[on, 1:3] + drv[on]
    z[on, 4:8] <- z[on, 4:8] + 0.8 * drv[on]
    mppi_fit(z, X, psych_idx = 2)
  }, simplify = FALSE)
  for (f in fits) {
    s <- mppi_seedset_map(f, seeds = 1:3, k = 1)
    m <- mppi_seedset_map(f, seeds = 1:3, k = 1, aggregate = "mean")
    expect_gt(sum(s$weights), 1)
    expect_gt(cor(s$map[4:8], m$map[4:8]), 0.9)
    expect_true(all(s$map[4:8] > 0))
  }
  two <- mppi_seedset_map(fits[[1]], seeds = 1:2, k = 1)
  expect_true(all(two$weights > 0))
  g <- mppi_group_seedset(fits, seeds = 1:3, k = 1, design = data.frame(i = rep(1, 12)))
  expect_true(all(g$t[4:8] > 10))
})

test_that("mechanistic helpers use the innovations series behind fit$Delta", {
  set.seed(15); Tn <- 200; V <- 4
  task <- rep(rep(c(1, 0), each = 10), Tn / 20)
  X <- cbind(int = 1, task = task)
  Y <- matrix(rnorm(Tn * V), Tn, V)
  f <- mppi_fit(Y, X, 2, domain = "innovations", scale = "cov", zero_diag = FALSE,
                deconv = list(hrf = mppi_default_hrf(1), method = "fixed", lambda = 1))
  D <- mppi_get_M_scaled(f, 1, "raw")
  expect_false(isTRUE(all.equal(f$U, f$R, check.attributes = FALSE)))
  expect_equal(mppi_lagged(f, 1, lags = 0L)$M_lag[["0"]], D, tolerance = 1e-8, check.attributes = FALSE)
  expect_equal(mppi_spectral_summary(f, 1)$lambda, eigen(D, symmetric = TRUE)$values, tolerance = 1e-8)
  Theta0 <- solve(crossprod(f$R) / nrow(f$R) + diag(1e-3, V))
  expect_equal(mppi_precision_delta(f, 1)$dTheta, -Theta0 %*% D %*% Theta0, tolerance = 1e-8, check.attributes = FALSE)
})

test_that("mppi_hrf_adapt combined matrix aligns with its inputs", {
  set.seed(8); Tn <- 300; V <- 5
  task <- rep(c(rep(1, 12), rep(0, 12)), length.out = Tn)
  h <- mppi_default_hrf(1)
  x1 <- as.numeric(stats::filter(c(rep(0, length(h) - 1), task), rev(h), sides = 1))[-seq_len(length(h) - 1)]
  x2 <- c(0, diff(x1))
  z <- matrix(rnorm(Tn * V), Tn, V)
  z[, 1:3] <- z[, 1:3] + 2 * pmax(x1, 0) * rnorm(Tn)
  for (sgn in c(1, -1)) {
    X <- cbind(int = 1, c1 = sgn * x1, c2 = sgn * x2)
    f <- mppi_fit(z, X, psych_idx = 2:3, scale = "cov")
    D1 <- mppi_get_M_scaled(f, 1, "raw"); D2 <- mppi_get_M_scaled(f, 2, "raw")
    expect_equal(sign(D1[1, 2]), sgn)
    combo <- mppi_hrf_adapt(list(D1, D2))
    expect_gt(cor(c(combo$Delta), c(D1)), 0.9)
    expect_gt(sum(combo$Delta * (D1 + D2)), 0)
    normed <- mppi_hrf_adapt(list(D1, D2), pk_norms = vapply(f$pk, function(p) sum(p^2), 0))
    expect_gt(sum(normed$Delta * (D1 + D2)), 0)
  }
})
