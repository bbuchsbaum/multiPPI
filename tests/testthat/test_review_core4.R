context("Review fixes (phase 4): HRF convention, design names, run labels, censored blocks, datasets, basis labels")

.core4_data <- function(Tn = 160L, V = 4L, seed = 41) {
  set.seed(seed)
  X <- cbind(rbinom(Tn, 1, 0.3), rbinom(Tn, 1, 0.3))
  Y <- matrix(rnorm(Tn * V), Tn, V)
  Y[, 2] <- Y[, 2] + X[, 1] * Y[, 1]
  list(X = X, Y = Y)
}

test_that("function-of-time HRFs work in the neural domain and match their sampled kernel", {
  d <- .core4_data()
  hf <- function(t) stats::dgamma(t, 6) - stats::dgamma(t, 16) / 6
  spm <- fmrihrf::HRF_SPMG1
  for (type in c("dct", "map")) {
    dc <- list(tr = 2, type = type, K = 40, method = "fixed", lambda = 0.1)
    f_fun <- mppi_fit(d$Y, d$X, 1:2, domain = "neural", deconv = c(dc, list(hrf = hf)))
    f_vec <- mppi_fit(d$Y, d$X, 1:2, domain = "neural", deconv = c(dc, list(hrf = hf(seq(0, 32, by = 2)))))
    expect_equal(f_fun$Delta, f_vec$Delta)
    expect_equal(f_fun$deconv$hrf, hf(seq(0, 32, by = 2)))
    # fmrihrf HRF objects are functions of time evaluated over their own span
    f_obj <- mppi_fit(d$Y, d$X, 1:2, domain = "neural", deconv = c(dc, list(hrf = spm)))
    f_ker <- mppi_fit(d$Y, d$X, 1:2, domain = "neural", deconv = c(dc, list(hrf = spm(seq(0, attr(spm, "span"), by = 2)))))
    expect_equal(f_obj$Delta, f_ker$Delta)
  }
  expect_error(mppi_fit(d$Y, d$X, 1:2, domain = "neural", deconv = list(tr = 2, hrf = fmrihrf::HRF_SPMG3)), "HRF must")
})

test_that("neural-domain Delta does not depend on HRF amplitude", {
  d <- .core4_data()
  h <- mppi_default_hrf(tr = 2)
  for (cfg in list(list(type = "map"), list(type = "map", lambda = 2),
                   list(type = "dct", method = "fixed"), list(type = "dct", method = "fixed", lambda = 1),
                   list(type = "dct", method = "gcv"))) {
    dc <- c(cfg, list(tr = 2, K = 40))
    a <- mppi_fit(d$Y, d$X, 1:2, domain = "neural", deconv = c(dc, list(hrf = h)))
    b <- mppi_fit(d$Y, d$X, 1:2, domain = "neural", deconv = c(dc, list(hrf = 10 * h)))
    for (k in 1:2) expect_equal(mppi_get_M_scaled(b, k, "normalized"), mppi_get_M_scaled(a, k, "normalized"), tolerance = 1e-8)
    expect_equal(10 * b$U, a$U, tolerance = 1e-8)
  }
  R <- qr.resid(qr(cbind(1, d$X)), d$Y)
  expect_equal(mppi_deconv(R, TR = 2, hrf = 10 * h) * 10, mppi_deconv(R, TR = 2, hrf = h), tolerance = 1e-8)
})

test_that("mppi_deconv and mppi_fit share the default HRF and relative lambda", {
  d <- .core4_data()
  R <- qr.resid(qr(cbind(1, d$X)), d$Y)
  fit <- mppi_fit(d$Y, d$X, 1:2, domain = "neural", deconv = list(tr = 2, type = "map"))
  expect_equal(fit$deconv$hrf, mppi_default_hrf(tr = 2))
  expect_equal(unname(fit$U), unname(mppi_deconv(R, TR = 2)), tolerance = 1e-10)
  expect_equal(mppi_deconv(R, TR = 2), mppi_deconv(R, TR = 2, hrf = mppi_default_hrf(tr = 2)))
  fit5 <- mppi_fit(d$Y, d$X, 1:2, domain = "neural", deconv = list(tr = 2, type = "map", lambda = 5))
  expect_equal(unname(fit5$U), unname(mppi_deconv(R, TR = 2, lambda = 5)), tolerance = 1e-10)
  dct <- mppi_fit(d$Y, d$X, 1:2, domain = "neural", deconv = list(tr = 2, K = 40, method = "fixed", lambda = 0.3))
  expect_equal(unlist(dct$deconv$lambda), rep(0.3, ncol(d$Y)))
  # DCT 'fixed' ridge is lambda times the mean diagonal of the convolved-basis Gram matrix
  n <- 50L; K <- 10L; h <- mppi_default_hrf(tr = 2); y <- R[seq_len(n), 1:2]
  B <- cos(outer(seq_len(n) - 0.5, seq_len(K) - 1) * pi / n) %*% diag(c(sqrt(1 / n), rep(sqrt(2 / n), K - 1L)))
  HB <- apply(B, 2, function(b) stats::filter(c(rep(0, length(h) - 1L), b), h, sides = 1)[-seq_len(length(h) - 1L)])
  G <- crossprod(HB)
  out <- multiPPI:::mppi_deconv_dct(y, h, K = K, method = "fixed", lambda = 0.3)
  expect_equal(out$U, B %*% solve(G + 0.3 * mean(diag(G)) * diag(K), crossprod(HB, y)), tolerance = 1e-8)
})

test_that("unnamed designs name regressors psych1..K and keep every lagged matrix", {
  d <- .core4_data()
  fb <- mppi_fit(d$Y, d$X, 1:2)
  expect_equal(fb$names, c("psych1", "psych2"))
  expect_equal(colnames(fb$fitted_X), c("(Intercept)", "psych1", "psych2"))
  expect_equal(mppi_get_M(fb, "psych2"), mppi_get_M(fb, 2L))
  for (type in c("dct", "map")) {
    fn <- mppi_fit(d$Y, d$X, 1:2, domain = "neural", deconv = list(tr = 2, type = type, K = 40))
    expect_equal(fn$names, c("s:psych1", "s:psych2"))
    expect_equal(mppi_get_M(fn, "s:psych1"), mppi_get_M(fn, 1L))
  }
  named <- d$X; colnames(named) <- c("A", "B")
  lg <- mppi_lagged(d$Y, d$X, 1:2, lags = -1:1)
  expect_length(lg, 6L)
  expect_equal(names(lg), paste0(rep(c("psych1", "psych2"), each = 3), "_lag", -1:1))
  expect_equal(unname(lg), unname(mppi_lagged(d$Y, named, 1:2, lags = -1:1)))
  # confounds and partly named or duplicated designs
  Xc <- cbind(1, rnorm(nrow(d$X)), d$X)
  expect_equal(mppi_fit(d$Y, Xc, 3:4)$names, c("psych1", "psych2"))
  expect_equal(colnames(mppi_fit(d$Y, Xc, 3:4)$fitted_X), c("x1", "x2", "psych1", "psych2"))
  colnames(Xc) <- c("", "drift", "A", "")
  expect_equal(mppi_fit(d$Y, Xc, 3:4)$names, c("A", "psych2"))
  colnames(Xc) <- c("int", "drift", "A", "A")
  expect_equal(mppi_fit(d$Y, Xc, 3:4)$names, c("psych1", "psych2"))
  colnames(Xc) <- c("psych1", "drift", "", "B")
  expect_false(anyDuplicated(colnames(mppi_fit(d$Y, Xc, 3:4)$fitted_X)) > 0)
})

test_that("factor run labels fit like integer labels", {
  d <- .core4_data()
  runs <- rep(1:2, each = 80L)
  fruns <- factor(rep(c("second", "first"), each = 80L), levels = c("first", "second"))
  for (cb in c("none", "run")) {
    fi <- mppi_fit(d$Y, d$X, 1:2, runs = runs, lags = 0:2, center_by = cb)
    ff <- mppi_fit(d$Y, d$X, 1:2, runs = fruns, lags = 0:2, center_by = cb)
    expect_equal(ff$Delta, fi$Delta)
    expect_equal(ff$lagged, fi$lagged)
    expect_equal(ff$lag_blocklens, fi$lag_blocklens)
    expect_equal(ff$runs, as.character(fruns))
  }
  dc <- list(tr = 2, K = 30)
  expect_equal(mppi_fit(d$Y, d$X, 1:2, runs = fruns, domain = "neural", deconv = dc)$Delta,
               mppi_fit(d$Y, d$X, 1:2, runs = runs, domain = "neural", deconv = dc)$Delta)
  Yc <- d$Y; Yc[10, ] <- NA
  expect_equal(mppi_fit(Yc, d$X, 1:2, runs = fruns, lags = 0:1)$lagged,
               mppi_fit(Yc, d$X, 1:2, runs = runs, lags = 0:1)$lagged)
  wi <- mppi_fit_whitened(d$Y, d$X, runs = runs, psych_idx = 1:2, p = 1L)
  wf <- mppi_fit_whitened(d$Y, d$X, runs = fruns, psych_idx = 1:2, p = 1L)
  expect_equal(wf$Delta, wi$Delta)
  des <- structure(list(X = cbind(1, d$X), psych_idx = 2:3, runs = as.character(fruns), basis = NULL), class = "mppi_design")
  expect_equal(mppi_runs(des), runs)
})

test_that("censored fits store contiguous blocks of retained rows", {
  d <- .core4_data(Tn = 120L)
  runs <- rep(1:2, each = 60L)
  Y <- d$Y; Y[c(30, 90), ] <- NA
  fit <- mppi_fit(Y, d$X, 1:2, runs = runs, lags = 0:2, zero_diag = FALSE, scale = "cov")
  expect_equal(fit$lag_blocklens, c(29L, 30L, 29L, 30L))
  expect_equal(sum(fit$lag_blocklens), nrow(fit$R))
  block <- rep(seq_along(fit$lag_blocklens), fit$lag_blocklens)
  for (b in unique(block)) {
    idx <- fit$observation_index[block == b]
    expect_true(all(diff(idx) == 1L))
    expect_length(unique(runs[idx]), 1L)
  }
  # a consumer pairing adjacent rows within blocks forms the fit's own lag-1 pairs
  for (k in 1:2) {
    pairs <- multiPPI:::.mppi_crossprod_lagged(fit$R, fit$pk[[k]], 1L, blocklens = fit$lag_blocklens)
    expect_equal(pairs, attr(mppi_get_M_lag(fit, k, 1L), "raw"), tolerance = 1e-10)
  }
  expect_equal(mppi_fit(Y, d$X, 1:2)$lag_blocklens, c(29L, 59L, 30L))
  expect_equal(mppi_fit(Y, d$X, 1:2, lag_blocklens = c(60, 60))$lag_blocklens, c(29L, 30L, 29L, 30L))
  expect_equal(mppi_fit(d$Y, d$X, 1:2, runs = runs)$lag_blocklens, c(60L, 60L))
  expect_null(mppi_fit(d$Y, d$X, 1:2)$lag_blocklens)
})

test_that("dataset path runs on a real legacy fmridataset matrix_dataset", {
  md <- if (requireNamespace("fmridataset", quietly = TRUE)) tryCatch(getExportedValue("fmridataset", "matrix_dataset"), error = function(e) NULL)
  skip_if(is.null(md), "needs fmridataset <= 0.10 (fmri_dataset API)")
  set.seed(4)
  Tn <- 200L
  X <- cbind(`(Intercept)` = 1, A = rbinom(Tn, 1, 0.4))
  Y <- matrix(rnorm(Tn * 5), Tn, 5)
  des <- structure(list(X = X, psych_idx = 2L, runs = NULL, basis = NULL), class = "mppi_design")
  fit <- mppi(md(Y, TR = 2, run_length = c(120, 80)), des, domain = "neural", deconv = list(K = 30))
  expect_equal(fit$deconv$tr, 2)
  expect_equal(fit$Delta, mppi(des, Y = Y, runs = rep(1:2, c(120, 80)), domain = "neural", deconv = list(K = 30, tr = 2))$Delta)
})

test_that("basis row names survive orthonormalization and label seed and edge outputs", {
  d <- .core4_data(V = 5L)
  labs <- paste0("roi", 1:5)
  colnames(d$Y) <- labs
  Vb <- matrix(rnorm(25), 5, 5, dimnames = list(labs, NULL))
  expect_equal(rownames(as_mppi_basis(Vb)$V), labs)
  expect_equal(rownames(as_mppi_basis(Vb[, 1:3])$V), labs)
  for (b in list(as_mppi_basis(Vb), Vb, list(V = Vb, name = "raw"))) {
    fit <- mppi_fit(d$Y, d$X, 1:2, basis = b)
    expect_equal(rownames(fit$basis$V), labs)
    s <- mppi_seedset_slopes(fit, seeds = c("roi1", "roi3"), k = 1L)
    expect_equal(dimnames(s), list(c("roi1", "roi3"), labs))
    expect_equal(unname(s), unname(mppi_seedset_slopes(fit, seeds = c(1L, 3L), k = 1L)))
    blk <- mppi_reconstruct_delta(fit, 1L, rows = 1:2, cols = 4:5)
    expect_equal(dimnames(blk), list(labs[1:2], labs[4:5]))
    expect_equal(mppi_edge(fit, 1L, 2L, 5L), blk["roi2", "roi5"])
  }
  des <- structure(list(X = cbind(1, d$X), psych_idx = 2:3, runs = NULL, basis = NULL), class = "mppi_design")
  expect_equal(rownames(mppi(des, Y = d$Y, basis = list(type = "pca", r = 3L))$basis$V), labs)
})
