test_that("joint conventional fit and contrasts agree with independent scalar GLMs", {
  set.seed(234)
  n <- 100
  seed <- rnorm(n)
  p <- cbind(A = rep(c(0, 1), each = n/2), B = rnorm(n))
  d <- ppi_design(seed, psych = p, hrf = 1, microtime = 1,
                  nuisance = cbind(motion = rnorm(n)), deconv = list(n_basis = n))
  Y <- cbind(y1 = 0.7 * seed + 0.4 * d$X[, "interaction:A:1"] + rnorm(n), y2 = rnorm(n))
  f <- ppi_fit(Y, d, noise = "iid", chunk_size = 1)
  ref <- lapply(seq_len(ncol(Y)), function(j) lm(Y[, j] ~ d$X - 1))
  for (j in seq_len(ncol(Y))) {
    expect_equal(unname(f$coefficients[, j]), unname(coef(ref[[j]])), tolerance = 1e-9)
    expect_equal(unname(f$standard_errors[, j]), unname(coef(summary(ref[[j]]))[, 2]), tolerance = 1e-9)
    expect_equal(unname(f$sigma2[j]), summary(ref[[j]])$sigma^2, tolerance = 1e-10)
  }
  z <- ppi_contrast(f, c(A = 1, B = -1))
  cc <- setNames(rep(0, ncol(d$X)), colnames(d$X)); cc[c("interaction:A:1", "interaction:B:1")] <- c(1, -1)
  expect_equal(unname(drop(z$estimate)), unname(vapply(ref, function(x) sum(cc * coef(x)), numeric(1))))
  expect_equal(unname(drop(z$standard_error)), unname(vapply(ref, function(x) sqrt(drop(t(cc) %*% vcov(x) %*% cc)), numeric(1))))
  joint <- ppi_contrast(f, structure(diag(2), dimnames = list(NULL, c("A", "B"))), test = "F")
  reduced <- d$X[, d$terms$type != "interaction", drop = FALSE]
  reference_F <- vapply(seq_len(ncol(Y)), function(j) anova(lm(Y[, j] ~ reduced - 1), ref[[j]])$F[2], numeric(1))
  expect_equal(unname(joint$statistic), unname(reference_F), tolerance = 1e-9)
  expect_equal(ppi_fit(Y, d, noise = "iid", chunk_size = 20)$coefficients, f$coefficients)
})

test_that("rank deficiency preserves estimable contrasts and marks aliases", {
  set.seed(432)
  s <- rnorm(80); p <- matrix(rnorm(80), ncol = 1, dimnames = list(NULL, "A"))
  nuisance <- cbind(one = rep(1, 80), twice = rep(2, 80))
  d <- ppi_design(s, psych = p, nuisance = nuisance, hrf = 1, microtime = 1)
  f <- ppi_fit(matrix(rnorm(160), 80), d, noise = "iid")
  expect_lt(f$rank, ncol(d$X))
  expect_true(all(is.na(f$coefficients[d$terms$type %in% c("intercept", "nuisance"), ])))
  expect_true(ppi_contrast(f, c(A = 1))$estimable)
  expect_false(ppi_contrast(f, c("intercept:1" = 1))$estimable)
  expect_true(all(is.na(ppi_contrast(f, c("intercept:1" = 1))$statistic)))
  expect_equal(f$df, nrow(d$X) - qr(d$X)$rank)
})

test_that("missing observations stay on their original grid under IID fitting", {
  set.seed(23); n <- 60
  s <- rnorm(n); p <- cbind(A = rep(c(0, 1), 30)); mask <- seq_len(n) %in% c(4, 9, 10)
  d <- ppi_design(s, psych = p, hrf = c(0, 1), microtime = 1, censor = mask)
  y <- matrix(rnorm(n * 2), n); y[mask, ] <- NA
  f <- ppi_fit(y, d, noise = "iid")
  ref <- lm(y[!mask, 1] ~ d$X[!mask, ] - 1)
  expect_equal(unname(f$coefficients[, 1]), unname(coef(ref)), tolerance = 1e-9)
  expect_error(ppi_fit(y, d), "censor gaps")
})

test_that("fixed fmriAR operator agrees with independent dense runwise GLS", {
  set.seed(122)
  n <- 120; rr <- rep(c(1, 2), c(55, 65))
  d <- ppi_design(rnorm(n), psych = cbind(A = rnorm(n), B = rnorm(n)),
                  runs = rr, hrf = 1, microtime = 1)
  y <- matrix(rnorm(n * 4), n)
  plan <- fmriAR::fit_noise(qr.resid(qr(d$X), y), runs = d$runs, p = 1L)
  f <- ppi_fit(y, d, noise = plan, chunk_size = 2L)
  phi <- as.numeric(plan$phi[[1]])
  covariance <- matrix(0, n, n)
  for (run in unique(rr)) {
    idx <- which(rr == run)
    covariance[idx, idx] <- phi^abs(outer(seq_along(idx), seq_along(idx), "-")) / (1 - phi^2)
  }
  inverse <- solve(covariance)
  vb <- solve(crossprod(d$X, inverse %*% d$X))
  beta <- vb %*% crossprod(d$X, inverse %*% y)
  e <- y - d$X %*% beta
  variance <- colSums(e * (inverse %*% e)) / (n - ncol(d$X))
  expect_equal(unname(f$coefficients), unname(beta), tolerance = 1e-9)
  expect_equal(unname(f$sigma2), unname(variance), tolerance = 1e-9)
  expect_equal(unname(f$standard_errors), unname(sqrt(outer(diag(vb), variance))), tolerance = 1e-9)
  expect_equal(ppi_fit(y, d, noise = "fmriAR", noise_control = list(p = 1L))$noise, "fmriAR")
  z <- ppi_contrast(f, c(A = 1, B = -1))
  weights <- z$contrast
  expect_equal(unname(z$standard_error), unname(sqrt(outer(diag(weights %*% vb %*% t(weights)), variance))))
})

test_that("a changing seed variance alone does not create conventional interaction", {
  set.seed(43)
  p <- rep(c(-1, 1), each = 60)
  seed <- c(scale(rnorm(60)), 2 * scale(rnorm(60)))
  d <- ppi_design(seed, psych = cbind(A = p), microtime = 1, hrf = 1, deconv = list(n_basis = 120))
  e <- qr.resid(qr(d$X), rnorm(120))
  f <- ppi_fit(matrix(0.7 * seed + e, ncol = 1), d, noise = "iid")
  expect_equal(drop(ppi_contrast(f, c(A = 1))$estimate), 0, tolerance = 1e-10)
})

test_that("README conventional example runs with the default HRF and noise model", {
  set.seed(1)
  n <- 100
  seed <- rnorm(n)
  psych <- cbind(A = rep(c(0, 1), 50), B = rep(c(0, 0, 1, 1), 25))
  design <- ppi_design(seed, psych = psych, TR = 2, run_effects = "shared")
  expect_equal(design$provenance$hrf, fmrihrf::hrf_spmg1(seq(0, 32, by = 2 / 16)))
  Y <- matrix(rnorm(n * 2), nrow = n)
  fit <- ppi_fit(Y, design, noise = "iid")
  est <- ppi_contrast(fit, c(A = 1, B = -1))$estimate
  expect_true(all(is.finite(est)))
  expect_s3_class(ppi_fit(Y, design)$noise_plan, "fmriAR_plan")
})
