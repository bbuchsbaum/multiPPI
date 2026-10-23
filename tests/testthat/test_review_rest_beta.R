# Regression tests for review items 3.10, 3.11, 3.14, 3.15, 3.16, 3.17

sim_state_data <- function() {
  set.seed(310)
  s <- rep(rep(c(1L, 2L, 1L, 3L), 8), times = rep(c(13L, 15L, 12L, 10L), 8))
  Tn <- length(s)
  U <- rbind(c(0, 0, 0), c(3, 0, 0), c(0, 3, 0))[s, ] + matrix(rnorm(Tn * 3, sd = 0.5), Tn)
  Y <- matrix(rnorm(Tn * 4), Tn, 4)
  Y[, 1:2] <- Y[, 1:2] + rnorm(Tn) * (s == 2)
  Y[, 3:4] <- Y[, 3:4] + rnorm(Tn) * (s == 3)
  list(s = s, occ = outer(s, 1:3, "==") * 1, U = U, Y = Y)
}

test_that("rest state modulators track their state's occupancy", {
  d <- sim_state_data()
  Tn <- length(d$s)
  for (soft in c(FALSE, TRUE)) {
    pk <- mppi_rest_states(d$U, K = 3, soft = soft)
    expect_named(pk, c("state2", "state3"))
    design <- mppi_rest_design(pk)
    fit <- mppi_rest(d$Y, modulators = pk, basis = "none", prewhiten = FALSE,
                     domain = "bold", scale = "cov")
    for (k in 2:3) {
      nm <- paste0("state", k)
      v <- design$pk[[nm]]
      expect_gt(cor(v, d$occ[, k]), 0.95)
      expect_lt(sum(v[c(1, Tn)]^2) / sum(v^2), 0.05)
      expect_gt(cor(fit$pk[[nm]], d$occ[, k]), 0.9)
    }
    D2 <- fit$Delta$state2; D3 <- fit$Delta$state3
    expect_gt(D2[1, 2], 3 * max(abs(D2[3, 4]), abs(D3[1, 2])))
    expect_gt(D3[3, 4], 3 * max(abs(D3[1, 2]), abs(D2[3, 4])))
  }
})

test_that("hard state slopes contrast each state with the reference state", {
  d <- sim_state_data()
  pk <- mppi_rest_states(d$U, K = 3, soft = FALSE)
  design <- mppi_rest_design(pk, scale_mod = FALSE)
  fit <- mppi_fit(d$Y, design$X, design$psych_idx, scale = "cov", zero_diag = FALSE)
  S <- cbind(1 - pk$state2 - pk$state3, pk$state2, pk$state3)
  R <- qr.resid(qr(cbind(1, S[, 2:3])), d$Y)
  C <- lapply(1:3, function(k) crossprod(R * S[, k], R) / sum(S[, k]))
  expect_equal(c(fit$Delta$state2), c(C[[2]] - C[[1]]), tolerance = 1e-8)
  expect_equal(c(fit$Delta$state3), c(C[[3]] - C[[1]]), tolerance = 1e-8)
  full <- list(a = S[, 1], b = S[, 2], c = S[, 3])
  expect_error(mppi_rest_design(full), "collinear")
})

test_that("multi-band fits remove task activation before filtering", {
  set.seed(4); Tn <- 400; V <- 6
  task <- rep(c(rep(1, 20), rep(0, 20)), Tn / 40)
  h <- mppi_default_hrf(tr = 1)
  xt <- as.numeric(stats::filter(c(rep(0, length(h) - 1), task), rev(h), sides = 1)[-(1:(length(h) - 1))])
  X <- cbind(int = 1, task = xt)
  E <- matrix(rnorm(Tn * V), Tn, V)
  bands <- list(list(low = NULL, high = 0.1, name = "lo"), list(low = 0.05, high = 0.2, name = "mid"))
  D <- lapply(c(0, 3, 10), function(amp) {
    fb <- mppi_fit_multi_band(outer(xt, rep(amp, V)) + E, X, 2L, bands = bands, scale = "cov")
    lapply(fb, function(f) c(mppi_get_M_scaled(f, 1, "raw")))
  })
  for (b in c("lo", "mid")) {
    expect_equal(D[[2]][[b]], D[[1]][[b]], tolerance = 1e-8)
    expect_equal(D[[3]][[b]], D[[1]][[b]], tolerance = 1e-8)
  }
})

test_that("base-R bandpass fallback removes out-of-band frequencies within runs", {
  skip_if(requireNamespace("signal", quietly = TRUE), "tests the base-R fallback")
  Tn <- 200; t <- seq_len(Tn)
  slow <- 3 * cos(2 * pi * 0.01 * t); mid <- sin(2 * pi * 0.12 * t); fast <- sin(2 * pi * 0.4 * t)
  x <- cbind(5 + slow + mid + fast)
  hp <- mppi_bandpass(x, low = 0.06, high = NULL)
  expect_gt(cor(hp[, 1], mid + fast), 0.99)
  expect_lt(abs(mean(hp)), 1e-8)
  bp <- mppi_bandpass(x, low = 0.06, high = 0.25)
  expect_gt(cor(bp[, 1], mid), 0.99)
  runs <- rep(c("a", "b"), each = Tn / 2)
  y <- cbind(x, rnorm(Tn))
  pr <- mppi_bandpass(y, low = 0.06, high = 0.25, runs = runs)
  expect_equal(pr, rbind(mppi_bandpass(y[1:100, ], low = 0.06, high = 0.25),
                         mppi_bandpass(y[101:200, ], low = 0.06, high = 0.25)))
  y2 <- y; y2[101:200, ] <- 10 * y2[101:200, ] + 4
  expect_equal(mppi_bandpass(y2, low = 0.06, high = 0.25, runs = runs)[1:100, ], pr[1:100, ])
})

test_that("Butterworth cutoffs are read in cycles/TR", {
  skip_if_not_installed("signal")
  t <- 1:400; x <- cbind(sin(2 * pi * 0.07 * t))
  y <- mppi_bandpass(x, low = NULL, high = 0.1)
  mid <- 51:350
  expect_gt(sd(y[mid, 1]) / sd(x[mid, 1]), 0.6)
})

test_that("beta permutation null varies for a continuous regressor by default", {
  set.seed(9); n <- 80; V <- 5
  P <- matrix(runif(n), n, 1, dimnames = list(NULL, "rt"))
  B <- matrix(rnorm(n * V), n, V) + outer(sqrt(3 * P[, 1]) * rnorm(n), rep(1, V))
  r <- mppi_beta_permute(B, P, Bperm = 99, seed = 1)
  expect_gt(sd(r$rt$Qnull), 0)
  expect_lt(r$rt$p_global, 0.05)
  expect_error(mppi_beta_permute(B, P, Bperm = 9, seed = 1, method = "block_flip"), "two blocks")
  Pb <- matrix(rep(0:1, each = n / 2), n, 1, dimnames = list(NULL, "cond"))
  expect_error(mppi_beta_permute(B, Pb, Bperm = 9, seed = 1, run_trial = rep(1:2, each = n / 2)),
               "constant within every block")
  rf <- mppi_beta_permute(B, P, Bperm = 49, seed = 1, blksize = 10)
  expect_equal(rf$rt$method, "block_flip")
  expect_gt(sd(rf$rt$Qnull), 0)
})

test_that("trait CV fits an in-fold intercept and aligns subjects by name", {
  set.seed(14); n <- 40
  feat <- do.call(rbind, lapply(1:n, function(i) data.frame(condition = c("s1", "s2"), G = rnorm(2), R = 0, mag = rnorm(2, 5), subj = i)))
  y <- 30 + 3 * feat$G[feat$condition == "s1"] + rnorm(n, sd = 0.5)
  r <- mppi_rest_trait_cv(feat, y, k = 5)
  expect_gt(r$cv_r2, 0.9)
  lab <- sprintf("sub-%02d", 1:n)
  featl <- feat; featl$subj <- lab[feat$subj]
  perm <- sample(n)
  rl <- mppi_rest_trait_cv(featl, setNames(y, lab)[perm], k = 5)
  expect_gt(rl$cv_r2, 0.9)
  expect_equal(unname(rl$observed[lab]), y)
  expect_error(mppi_rest_trait_cv(feat, setNames(y, lab), k = 5), "names must match")
})

test_that("PWI phase and frequency are run-local and accept character runs", {
  set.seed(17); Tn <- 400; runs <- rep(1:2, each = 200)
  U <- matrix(rnorm(Tn * 3), Tn, 3); U2 <- U; U2[runs == 2, ] <- U2[runs == 2, ] * 5 + 3
  a <- multiPPI:::ipppi_instfreq(U, fs = 0.5, band = c(0.01, 0.1), runs = runs)
  b <- multiPPI:::ipppi_instfreq(U2, fs = 0.5, band = c(0.01, 0.1), runs = runs)
  expect_equal(b$phi[runs == 1, ], a$phi[runs == 1, ])
  expect_equal(b$omega[runs == 1, ], a$omega[runs == 1, ])
  r1 <- multiPPI:::ipppi_instfreq(U[runs == 1, ], fs = 0.5, band = c(0.01, 0.1))
  expect_equal(a$phi[runs == 1, ], r1$phi)
  cr <- multiPPI:::ipppi_instfreq(U, fs = 0.5, band = c(0.01, 0.1), runs = c("a", "b")[runs])
  expect_equal(cr$omega, a$omega)
  expect_equal(sum(is.na(cr$omega)), 2 * ncol(U))
})

test_that("mppi_coerce_beta keeps missing cells as NA and rejects duplicates", {
  df <- data.frame(trial = c(1, 1, 2, 2, 3), roi = c("A", "B", "A", "B", "A"), beta = c(1, 2, 3, NA, 5))
  B <- mppi_coerce_beta(df)
  expect_equal(unname(B), matrix(c(1, 3, 5, 2, NA, NA), 3, 2))
  dup <- rbind(df, data.frame(trial = 3, roi = "A", beta = 7))
  expect_error(mppi_coerce_beta(dup), "duplicate")
})
