test_that("rest seeds leave the caller's RNG stream untouched", {
  set.seed(5)
  U <- matrix(rnorm(300), 100, 3)
  feats <- data.frame(subj = rep(1:10, each = 2), condition = rep(c("a", "b"), 10),
                      G = rnorm(20), R = rnorm(20), mag = rnorm(20))
  y <- rnorm(10)
  before <- get(".Random.seed", envir = globalenv())
  st <- mppi_rest_states(U, K = 3, soft = FALSE)
  expect_true(identical(get(".Random.seed", envir = globalenv()), before))
  cv <- mppi_rest_trait_cv(feats, y, k = 5)
  expect_true(identical(get(".Random.seed", envir = globalenv()), before))
  expect_identical(mppi_rest_states(U, K = 3, soft = FALSE), st)
  expect_identical(mppi_rest_trait_cv(feats, y, k = 5), cv)
})

test_that("inst_corr skips non-finite samples and is scale invariant", {
  set.seed(1); x <- rnorm(50); y <- x + rnorm(50)
  xn <- x; xn[10] <- NaN; xa <- x; xa[10] <- NA
  r <- inst_corr(xn, y, tau_half = 5)
  expect_identical(r, inst_corr(xa, y, tau_half = 5))
  expect_true(is.na(r[10]))
  expect_true(all(is.finite(r[-10])))
  expect_equal(inst_corr(x * 1e-4, y * 1e-4, tau_half = 5), inst_corr(x, y, tau_half = 5))
  expect_true(all(is.na(inst_corr(c(0, 0, 1), c(1, 1, 1), tau_half = 2)[1:2])))
})

test_that("neural rest modulators deconvolve each run separately", {
  set.seed(4); n <- 120; runs <- rep(1:2, each = 60); i2 <- runs == 2
  Y <- matrix(rnorm(n * 4), n); attr(Y, "TR") <- 2
  Y2 <- Y; Y2[i2, ] <- Y[i2, ] + matrix(rnorm(sum(i2) * 4, sd = 3), ncol = 4)
  tc <- function(Y) attr(mppi_rest(Y, runs = runs, basis = "none", include = "envelopes",
                                   modulator_domain = "neural", domain = "bold", prewhiten = FALSE),
                         "basis_timecourses")
  a <- tc(Y); b <- tc(Y2)
  expect_equal(a[!i2, ], b[!i2, ])
  expect_false(isTRUE(all.equal(a[i2, ], b[i2, ])))
})

test_that("iphs filters each run separately", {
  set.seed(8); n <- 160; runs <- rep(1:2, each = 80); i2 <- runs == 2
  X <- cbind(run1 = !i2, run2 = i2, task = ifelse(i2, 0, sin(seq_len(n) / 4)))
  Y <- matrix(rnorm(n * 3), n)
  Y2 <- Y; Y2[i2, ] <- Y[i2, ] + matrix(rnorm(sum(i2) * 3, sd = 2), ncol = 3)
  opt <- list(method = "iphs", tau_half = 40)
  a <- mppi(Y, X, 3, runs = runs, engine = "instant", instant = opt)
  b <- mppi(Y2, X, 3, runs = runs, engine = "instant", instant = opt)
  expect_equal(a$Delta, b$Delta)
  expect_true(all(is.finite(a$Delta[[1]])))
})

test_that("a censored scan ends an iphs/pwi filtering segment", {
  set.seed(9); n <- 161; c0 <- 81; tt <- seq_len(n); post <- tt > c0; keep <- tt != c0
  X <- cbind(a = tt < c0, b = tt >= c0, task = ifelse(tt < c0, sin(tt / 4), 0))
  Y <- matrix(rnorm(n * 3), n); Y[c0, 1] <- NA
  Y2 <- Y; Y2[post, ] <- Y[post, ] + matrix(rnorm(sum(post) * 3, sd = 2), ncol = 3)
  opt <- list(method = "iphs", tau_half = 40)
  expect_equal(mppi(Y, X, 3, engine = "instant", instant = opt)$Delta,
               mppi(Y2, X, 3, engine = "instant", instant = opt)$Delta)
  # the gap acts exactly like a run boundary
  runs2 <- rep(1:2, c(c0 - 1, n - c0))
  for (opt in list(list(method = "pwi"), list(method = "iphs", tau_half = 1e9))) {
    a <- mppi(Y, X, 3, engine = "instant", instant = opt)
    b <- mppi(Y[keep, ], X[keep, ], 3, runs = runs2, engine = "instant", instant = opt)
    expect_equal(a$Delta, b$Delta)
  }
})

test_that("iphs is an end-weighted slope and pwi is amplitude-free", {
  fs <- 1.4; Tn <- 400; tt <- seq(0, length.out = Tn, by = 1 / fs)
  U <- cbind(sin(2 * pi * 0.05 * tt), sin(2 * pi * 0.05 * tt + pi / 3))
  set.seed(3); p <- matrix(rnorm(Tn), Tn)
  s1 <- ipppi_iphs(U, p, fs = fs, band = c(0.04, 0.06), tau_half = 12)[[1]]
  expect_equal(ipppi_iphs(U, 3 * p, fs = fs, band = c(0.04, 0.06), tau_half = 12)[[1]], s1 / 3)
  phi <- ipppi_iphs_phase(U, fs, c(0.04, 0.06)); w <- 2^((seq_len(Tn) - Tn) / 12)
  expect_equal(s1[1, 2], sum(w * p * cos(phi[, 1] - phi[, 2])) / sum(w * p^2))
  set.seed(6); n <- 120; X <- cbind(1, task = sin(seq_len(n) / 5)); Y <- matrix(rnorm(n * 3), n)
  a <- mppi(Y, X, 2, engine = "instant", instant = list(method = "pwi"))
  b <- mppi(10 * Y, X, 2, engine = "instant", scale = "cov", instant = list(method = "pwi"))
  expect_equal(mppi_get_M(a, 1), mppi_get_M(b, 1))
  expect_equal(a$Delta_normalized, a$Delta_raw)
})

test_that("instant engine maps lags and lag_blocklens to offsets", {
  set.seed(13); n <- 60; runs <- rep(1:2, each = 30)
  X <- cbind(1, task = rnorm(n)); Y <- matrix(rnorm(n * 3), n)
  big <- list(tau_half = 1e12)
  fc <- mppi(Y, X, 2, runs = runs, lags = c(-1, 2), zero_diag = FALSE, scale = "cov")
  fi <- mppi(Y, X, 2, runs = runs, lags = c(-1, 2), zero_diag = FALSE, scale = "cov",
             engine = "instant", instant = big)
  expect_equal(fi$lags, c(-1L, 0L, 2L))
  for (lg in c(-1L, 2L)) expect_equal(c(mppi_get_M_lag(fi, 1, lg)), c(mppi_get_M_lag(fc, 1, lg)), tolerance = 1e-8)
  expect_equal(mppi(Y, X, 2, lags = 2, engine = "instant", instant = list(offsets = 2))$lags, c(0L, 2L))
  expect_error(mppi(Y, X, 2, lags = 1, engine = "instant", instant = list(offsets = 2)), "disagree")
  fcb <- mppi(Y, X, 2, lags = 2, lag_blocklens = c(30, 30), zero_diag = FALSE, scale = "cov")
  fib <- mppi(Y, X, 2, lags = 2, lag_blocklens = c(30, 30), zero_diag = FALSE, scale = "cov",
              engine = "instant", instant = big)
  expect_equal(c(mppi_get_M_lag(fib, 1, 2)), c(mppi_get_M_lag(fcb, 1, 2)), tolerance = 1e-8)
  expect_equal(fib$lag_blocklens, c(30L, 30L))
  skip_if_not_installed("fmridesign")
  sf <- fmridesign::sampling_frame(blocklens = c(30, 30), TR = 2)
  ev <- data.frame(onset = rep(c(4, 20, 36), 2), run = rep(1:2, each = 3), cond = factor(rep(c("A", "B", "A"), 2)))
  hrf <- fmridesign::hrf
  em <- fmridesign::event_model(onset ~ hrf(cond), data = ev, block = ~run, sampling_frame = sf)
  fd <- mppi(em, Y = Y, lags = 1, engine = "instant")
  expect_equal(fd$lags, c(0L, 1L))
})

test_that("instant engine honours zero_diag and rejects options it cannot honour", {
  set.seed(14); n <- 60; X <- cbind(1, task = rnorm(n)); Y <- matrix(rnorm(n * 3), n)
  fz <- mppi(Y, X, 2, lags = 1, engine = "instant")
  fk <- mppi(Y, X, 2, lags = 1, zero_diag = FALSE, engine = "instant")
  expect_true(fz$zero_diag)
  for (nm in c("Delta", "Delta_raw", "Delta_amplitude", "Delta_normalized")) expect_equal(diag(fz[[nm]][[1]]), rep(0, 3))
  expect_equal(diag(mppi_get_M_lag(fz, 1, 1)), rep(0, 3))
  expect_true(all(diag(fk$Delta_raw[[1]]) != 0))
  M <- fk$Delta[[1]]; diag(M) <- 0
  expect_equal(fz$Delta[[1]], M)
  fc <- mppi(Y, X, 2, scale = "cov")
  fi <- mppi(Y, X, 2, scale = "cov", engine = "instant", instant = list(tau_half = 1e12))
  expect_equal(c(mppi_get_M(fi, 1)), c(mppi_get_M(fc, 1)), tolerance = 1e-8)
  expect_error(mppi(Y, X, 2, engine = "instant", backend = "chunked"), "not supported")
  expect_error(mppi(Y, X, 2, engine = "instant", instant = list(method = "pwi", offsets = 1)), "only supported")
})

test_that("mppi_rest_features labels subjects by list names", {
  set.seed(3); Y <- matrix(rnorm(600), 100, 6)
  fit <- mppi_rest(Y, basis = list(type = "pca", r = 4), prewhiten = FALSE,
                   include = c("states", "bursts"), modulator_domain = "bold", domain = "bold")
  fits <- list(subB = fit, subA = fit)
  long <- mppi_rest_features(fits, modes = 2L, aggregate = "none")
  expect_setequal(unique(long$subj), c("subB", "subA"))
  expect_equal(mppi_rest_features(fits, modes = 2L)$subj, c("subA", "subB"))
  expect_equal(unique(mppi_rest_features(unname(fits), modes = 2L, aggregate = "none")$subj), 1:2)
  expect_error(mppi_rest_features(list(a = fit, fit), modes = 2L), "unique and non-empty")
})
