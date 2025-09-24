test_that("mppi_fit_beta exposes residual regressors and permutation works", {
  set.seed(123)
  ntr <- 40L; V <- 6L
  B <- matrix(rnorm(ntr * V), ntr, V)
  Pbeta <- cbind(cond = sample(c(0, 1), ntr, replace = TRUE),
                 mod  = as.numeric(scale(rnorm(ntr), center = TRUE, scale = FALSE)))
  colnames(Pbeta) <- c("cond", "mod")

  fit <- mppi_fit_beta(B, Pbeta, zero_diag = TRUE, scale = "cov")
  expect_equal(sort(names(fit$pk)), sort(colnames(Pbeta)))
  expect_true(all(vapply(fit$pk, length, integer(1)) == ntr))
  expect_true(all(vapply(fit$meta, function(x) all(c("Q", "qr", "raw", "is_binary") %in% names(x)), logical(1))))

  perms <- mppi_beta_permute(B, Pbeta, Bperm = 20L, blksize = 10L, zero_diag = TRUE)
  expect_named(perms, colnames(Pbeta))
  expect_true(all(vapply(perms, function(x) length(x$Qnull), integer(1)) == 20L))
  expect_true(all(vapply(perms, function(x) is.matrix(x$D) && nrow(x$D) == V, logical(1))))
  expect_true(all(vapply(perms, function(x) is.na(x$p_global) || (x$p_global >= 0 && x$p_global <= 1), logical(1))))
})

test_that("packed beta output can be unpacked", {
  set.seed(456)
  ntr <- 25L; V <- 4L
  B <- matrix(rnorm(ntr * V), ntr, V)
  Pbeta <- cbind(cond = sample(c(0, 1), ntr, replace = TRUE))
  fit <- mppi_fit_beta(B, Pbeta, packed = TRUE)
  expect_true(fit$packed)
  packed <- fit$Delta[[1]]
  full <- multiPPI:::`.mppi_unpack_upper`(packed)
  fit_full <- mppi_fit_beta(B, Pbeta, packed = FALSE)
  expect_equal(full, fit_full$Delta[[1]], tolerance = 1e-10)
})

test_that("mppi_lag_select handles multiple targets and per-target output", {
  set.seed(321)
  Tn <- 60L; V <- 5L
  Y <- matrix(rnorm(Tn * V), Tn, V)
  p1 <- sin(seq_len(Tn) / 5)
  p2 <- cos(seq_len(Tn) / 7)
  X <- cbind(Intercept = 1, p1 = p1, p2 = p2)

  res_multi <- mppi_lag_select(Y, X, psych_idx = 2:3, lags = -1:1,
                                blocklens = c(30, 30), B = 20L, blksize = NULL,
                                targets = c("p1", "p2"))
  expect_named(res_multi, c("p1", "p2"))
  expect_true(all(vapply(res_multi, function(x) length(x$Qobs), integer(1)) == length(-1:1)))
  expect_true(all(vapply(res_multi, function(x) is.na(x$best_lag) || x$best_lag %in% -1:1, logical(1))))

  res_single <- mppi_lag_select(Y, X, psych_idx = 2:3, lags = -1:1,
                                 blocklens = c(30, 30), B = 10L, blksize = NULL,
                                 targets = "p1")
  expect_true(is.list(res_single))
  expect_true(all(c("best_lag", "Q", "Qnull", "p", "Qobs", "Dbest") %in% names(res_single)))
  expect_equal(length(res_single$Qnull), 10L)
})
