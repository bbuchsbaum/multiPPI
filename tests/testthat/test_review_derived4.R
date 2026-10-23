library(testthat)

# Phase 4 review items 4.1 (classify RNG), 4.5-4.9, 4.21

block_design <- function(Tn, len = 20) {
  task <- rep(rep(c(1, 0), each = len), length.out = Tn)
  X <- cbind(`(Intercept)` = 1, task = task)
  list(X = X, fd = structure(list(X = X, psych_idx = 2L, runs = rep(1L, Tn)), class = "mppi_design"))
}

test_that("mppi_classify restores the caller's RNG stream (4.1)", {
  set.seed(2); Tn <- 200
  d <- block_design(Tn, 10)
  fit <- mppi_fit(matrix(rnorm(Tn * 4), Tn, 4), d$X, psych_idx = 2, lags = -1:1)
  set.seed(42); s0 <- .Random.seed
  a <- mppi_classify(fit, 1, lags = -1:1, B = 20, seed = 1L)
  expect_true(identical(.Random.seed, s0))
  b <- mppi_classify(fit, 1, lags = -1:1, B = 20, seed = 1L)
  expect_identical(a[c("pG", "pC", "pR")], b[c("pG", "pC", "pR")])
  set.seed(7); s7 <- .Random.seed
  u1 <- mppi_classify(fit, 1, lags = -1:1, B = 20, seed = NULL)
  s_after <- .Random.seed
  expect_false(identical(s_after, s7))
  mppi_classify(fit, 1, lags = -1:1, B = 20, seed = NULL)
  expect_false(identical(.Random.seed, s_after))
  set.seed(7)
  u2 <- mppi_classify(fit, 1, lags = -1:1, B = 20, seed = NULL)
  expect_identical(u1[c("pG", "pC", "pR")], u2[c("pG", "pC", "pR")])
})

test_that("mppi_sdppi counts each event of a convolved regressor once (4.5)", {
  skip_if_not_installed("fmridesign")
  suppressPackageStartupMessages(library(fmridesign))
  sf <- fmridesign::sampling_frame(blocklens = c(100, 100), TR = 2)
  ev <- data.frame(onset = c(10, 80, 150, 20, 90, 160), cond = factor(rep(c("a", "b"), each = 3)),
                   run = rep(1:2, each = 3))
  des <- mppi_model(fmridesign::event_model(onset ~ hrf(cond), data = ev, block = ~ run, sampling_frame = sf))
  set.seed(1)
  fit <- mppi_fit(matrix(rnorm(200 * 4), 200, 4), des$X, psych_idx = des$psych_idx, runs = des$runs)
  sd <- mppi_sdppi(fit, des, lags = -1:1)
  # scans sample at TR/2 + (i - 1) * TR; trials start at the first scan of each response
  expect_equal(sd$onset, c(6, 41, 76, 111, 146, 181))
  expect_equal(sd$trial, c(1:3, 1:3))
  expect_equal(sd$run, rep(1:2, each = 3))

  ev2 <- data.frame(onset = seq(10, 180, by = 8), cond = factor("a"), run = 1)
  des2 <- mppi_model(fmridesign::event_model(onset ~ hrf(cond), data = ev2, block = ~ run,
                                             sampling_frame = fmridesign::sampling_frame(blocklens = 100, TR = 2)))
  fit2 <- mppi_fit(matrix(rnorm(300), 100, 3), des2$X, psych_idx = des2$psych_idx)
  expect_error(mppi_sdppi(fit2, des2), "overlap")
  on2 <- ev2$onset / 2 + 1
  sd2 <- mppi_sdppi(fit2, des2, onsets = on2, lags = integer(0))
  expect_equal(sd2$onset, on2)
})

test_that("mppi_sdppi trial magnitude is a modulation, zero without one (4.5)", {
  Tn <- 200; d <- block_design(Tn)
  task <- d$X[, 2]
  # residual outer products identical at every scan: activation, no modulation
  eps <- numeric(Tn); eps[task == 1] <- rep(c(1, -1), 50); eps[task == 0] <- rep(c(-1, 1), 50)
  base <- eps %o% c(3, 2, 1, 2)
  fit <- mppi_fit(3 * task + base, d$X, psych_idx = 2)
  expect_equal(fit$R, base, tolerance = 1e-10, check.attributes = FALSE)
  for (m in c("raw", "normalized")) {
    sd <- mppi_sdppi(fit, d$fd, lags = integer(0), mode = m)
    expect_equal(sd$onset, seq(1, Tn, by = 40))
    expect_true(all(sd$magnitude < 1e-8))
  }
  set.seed(3); g <- rnorm(Tn)
  fit_mod <- mppi_fit(3 * task + base + (task * g) %o% c(1, -1, 1, 0), d$X, psych_idx = 2)
  expect_true(all(mppi_sdppi(fit_mod, d$fd, lags = integer(0), mode = "raw")$magnitude > 0.5))
})

test_that("mppi_sdppi uses the series behind fit$Delta and lead-lag directedness (4.21)", {
  # node 1 drives node 2 one scan later during trials
  set.seed(4); Tn <- 400
  d <- block_design(Tn)
  Z <- matrix(rnorm(Tn * 3), Tn, 3)
  Z[-1, 2] <- Z[-1, 2] + 2 * d$X[-1, 2] * Z[-Tn, 1]
  fl <- mppi_fit(Z, d$X, psych_idx = 2)
  r <- mppi_sdppi(fl, d$fd, lags = -2:2, mode = "raw")$routing
  expect_equal(length(r), 10L)
  expect_gt(stats::median(r), 0.5)

  set.seed(15); Tn <- 200
  d <- block_design(Tn, 10)
  f <- mppi_fit(matrix(rnorm(Tn * 4), Tn, 4), d$X, 2, domain = "innovations", scale = "cov",
                zero_diag = FALSE, deconv = list(hrf = mppi_default_hrf(1), method = "fixed", lambda = 1))
  expect_false(isTRUE(all.equal(f$U, f$R, check.attributes = FALSE)))
  sd <- mppi_sdppi(f, d$fd, pre_window = 10L, lags = integer(0), mode = "raw")
  v1 <- eigen(f$Sigma0, symmetric = TRUE)$vectors[, 1]
  Sp <- stats::cov(f$R[11:20, ])
  expect_equal(sd$onset[2], 21L)
  expect_equal(sd$v1_dominance[2], drop(v1 %*% Sp %*% v1) / sum(diag(Sp)), tolerance = 1e-10)
  # one trial covering every scan reproduces the fit's own slope and mppi_lagged()
  sd1 <- mppi_sdppi(f, d$fd, onsets = 1L, mode = "raw", lags = -2:2)
  expect_equal(sd1$magnitude, sqrt(sum(mppi_get_M_scaled(f, 1, "raw")^2)), tolerance = 1e-10)
  expect_equal(sd1$routing, mppi_lagged(f, 1, lags = -2:2)$R, tolerance = 1e-10)
})

test_that("mppi_sdppi_fit windows stay within the onset's run (4.6)", {
  set.seed(6); Tn <- 200; runs <- rep(1:2, each = 100)
  d <- block_design(Tn, 10)
  fit <- mppi_fit(matrix(rnorm(Tn * 3), Tn, 3), d$X, psych_idx = 2, runs = runs)
  feats <- function(Xs, pre, post) {
    lam2 <- eigen(crossprod(Xs[pre, ]) / length(pre), symmetric = TRUE)$values^2
    c(pre_lambda1_share = lam2[1] / (sum(lam2) + 1e-12),
      ppi_energy = sqrt(sum((crossprod(Xs[post, ]) / (length(post) + 1e-12))^2)))
  }
  sd <- mppi_sdppi_fit(fit, 1, onsets_idx = c(50, 97, 101, 105), runs = runs, Wpre = 8, Wpost = 8)
  # 101 starts run 2: no pre-window from run 1, so it is skipped
  expect_equal(sd$data$run, c(1L, 1L, 2L))
  expect_equal(unlist(sd$data[2, c("pre_lambda1_share", "ppi_energy")]), feats(fit$R, 89:96, 97:100), tolerance = 1e-10)
  expect_equal(unlist(sd$data[3, c("pre_lambda1_share", "ppi_energy")]), feats(fit$R, 101:104, 105:112), tolerance = 1e-10)
  expect_equal(mppi_sdppi_fit(fit, 1, onsets_idx = c(50, 97, 101, 105), Wpre = 8, Wpost = 8)$data, sd$data)

  si <- mppi_sdppi_fit(fit, 1, onsets_idx = c(50, 105), runs = runs, Wpre = 8, Wpost = 8, mode = "innov",
                       protect = "none", ar_order = "fixed", p = 2L)
  Xi <- rbind(mppi_innovations(fit$R[1:100, ], protect = "none", ar_order = "fixed", p = 2L, center = FALSE),
              mppi_innovations(fit$R[101:200, ], protect = "none", ar_order = "fixed", p = 2L, center = FALSE))
  expect_equal(unlist(si$data[2, c("pre_lambda1_share", "ppi_energy")]), feats(Xi, 101:104, 105:112), tolerance = 1e-10)

  # censored scan 103: windows live on the acquisition grid
  Y <- matrix(rnorm(Tn * 3), Tn, 3); Y[103, ] <- NA
  fc <- mppi_fit(Y, d$X, psych_idx = 2, runs = runs)
  sc <- mppi_sdppi_fit(fc, 1, onsets_idx = 108, runs = fc$runs, Wpre = 6, Wpost = 4)
  rows <- function(scans) match(scans, fc$observation_index)
  expect_equal(unlist(sc$data[1, c("pre_lambda1_share", "ppi_energy")]),
               feats(fc$R, rows(c(102, 104:107)), rows(108:111)), tolerance = 1e-10)
})

test_that("mppi_state_predict fits real folds for every response and run type (4.7)", {
  set.seed(8); Tn <- 600; runs <- rep(1:3, each = 200)
  d <- block_design(Tn, 10)
  fit <- mppi_fit(matrix(rnorm(Tn * 4), Tn, 4), d$X, psych_idx = 2, runs = runs)
  onsets <- unlist(lapply(c(0, 200, 400), function(o) o + seq(15, 185, by = 10)))
  sd <- mppi_sdppi_fit(fit, 1, onsets_idx = onsets, runs = runs, Wpre = 6, Wpost = 6)
  expect_equal(nrow(sd$data), length(onsets))
  lin <- 2 * sd$data$pre_effrank - 3 * sd$data$pre_lambda1_share
  y <- lin + rnorm(length(lin), sd = 0.01 * sd(lin))
  g <- mppi_state_predict(sd$data, y, family = "gaussian")
  expect_gt(g$cv_r, 0.95)
  expect_equal(mppi_state_predict(sd$data, y, runs = paste0("run", sd$data$run))$cv_r, g$cv_r)
  hi <- lin + rnorm(length(lin), sd = 0.3 * sd(lin)) > stats::median(lin)
  yf <- factor(ifelse(hi, "hi", "lo"), levels = c("lo", "hi"))
  b <- suppressWarnings(mppi_state_predict(sd$data, yf, family = "binomial"))
  expect_gt(b$cv_auc, 0.8)
  expect_equal(b, suppressWarnings(mppi_state_predict(sd$data, as.integer(hi), family = "binomial")))
  expect_error(mppi_state_predict(sd$data, factor(rep(c("a", "b", "c"), length.out = nrow(sd$data))), family = "binomial"), "two classes")

  # character run labels in mppi_sdppi_fit's own cross-validation
  sc <- mppi_sdppi_fit(fit, 1, onsets_idx = onsets, runs = paste0("r", runs), Wpre = 6, Wpost = 6)
  expect_equal(sc$cv_r, sd$cv_r)
  expect_true(is.finite(sc$cv_r))
})

test_that("mppi_seedset_slopes reads lagged slopes of cov and corr fits (4.8)", {
  set.seed(1); Tn <- 200
  d <- block_design(Tn, 10)
  z <- matrix(rnorm(Tn * 5), Tn, 5)
  fc <- mppi_fit(z, d$X, psych_idx = 2, lags = 0:1, scale = "cov")
  fr <- mppi_fit(z, d$X, psych_idx = 2, lags = 0:1, scale = "corr")
  Lc <- mppi_get_M_lag(fc, 1, 1); Lr <- mppi_get_M_lag(fr, 1, 1)
  expect_equal(mppi_seedset_slopes(fc, seeds = 1:2, k = 1, mode = "raw", lag = 1), Lc[1:2, ], check.attributes = FALSE)
  expect_equal(mppi_seedset_slopes(fr, seeds = 1:2, k = 1, mode = "normalized", lag = 1), Lr[1:2, ], check.attributes = FALSE)
  for (m in c("raw", "amplitude", "normalized")) {
    expect_equal(mppi_seedset_slopes(fr, seeds = 1, k = 1, mode = m, lag = 1),
                 mppi_seedset_slopes(fc, seeds = 1, k = 1, mode = m, lag = 1))
  }
  expect_equal(mppi_seedset_slopes(fr, seeds = 1, k = 1, mode = "raw", lag = 1), Lc[1, , drop = FALSE], check.attributes = FALSE)
})

test_that("mppi_group_seedset tests the requested design column (4.9)", {
  set.seed(2); Tn <- 100; V <- 4
  d <- block_design(Tn, 10)
  fits <- replicate(6, mppi_fit(matrix(rnorm(Tn * V), Tn, V), d$X, psych_idx = 2), simplify = FALSE)
  df <- data.frame(group = rep(c("a", "b"), 3))
  M <- t(sapply(fits, function(f) mppi_seedset_map(f, 1, 1, aggregate = "mean")$map))
  ref <- lapply(seq_len(V), function(j) summary(stats::lm(M[, j] ~ group, df))$coefficients)
  g <- mppi_group_seedset(fits, seeds = 1, k = 1, aggregate = "mean", design = df, formula = ~ group)
  expect_equal(unname(g$coef), sapply(ref, function(x) x["groupb", "Estimate"]))
  expect_equal(unname(g$t), sapply(ref, function(x) x["groupb", "t value"]))
  expect_equal(unname(g$contrast), c(0, 1))
  gi <- mppi_group_seedset(fits, seeds = 1, k = 1, aggregate = "mean", design = df, formula = ~ group, coef = "(Intercept)")
  expect_equal(unname(gi$coef), sapply(ref, function(x) x["(Intercept)", "Estimate"]))
  gb <- mppi_group_seedset(fits, seeds = 1, k = 1, aggregate = "mean", design = df, formula = ~ group, contrast = c(1, 1))
  expect_equal(unname(gb$coef), unname(colMeans(M[df$group == "b", ])))
})

test_that("normalized sdPPI gain does not depend on ROI units", {
  Tn <- 200; d <- block_design(Tn); task <- d$X[, 2]
  set.seed(8); Y <- matrix(rnorm(Tn * 4), Tn) %*% chol(matrix(0.5, 4, 4) + diag(0.5, 4))
  Y <- Y + (task * rnorm(Tn)) %o% c(1, -1, 1, 0)
  a <- mppi_sdppi(mppi_fit(Y, d$X, psych_idx = 2), d$fd, lags = integer(0), mode = "normalized")
  Y2 <- Y; Y2[, 1] <- 100 * Y2[, 1]
  b <- mppi_sdppi(mppi_fit(Y2, d$X, psych_idx = 2), d$fd, lags = integer(0), mode = "normalized")
  expect_equal(b$magnitude, a$magnitude, tolerance = 1e-8)
  expect_equal(b$gain, a$gain, tolerance = 1e-8)
})
