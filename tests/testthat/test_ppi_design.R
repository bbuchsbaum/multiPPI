test_that("masked DCT inverse and interactions match dense causal operators", {
  n <- 6L; mt <- 2L; dt <- 1; m <- n * mt
  h <- c(0, 1, 0.5, 0); seed <- c(1, 3, 9, 2, 0, 4)
  p <- cbind(A = rep(c(0, 1, 1, 0, 0, 1), each = mt), B = rep(c(0, 1, 1, 0, 1, 0), each = mt))
  keep <- c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE); k <- 4L; lambda <- 0.3
  d <- ppi_design(seed, psych = list(values = p, grid = "neural"), TR = 2,
                  microtime = mt, hrf = h, censor = !keep,
                  deconv = list(n_basis = k, lambda = lambda), center = TRUE, seed_adjust="none")
  H <- outer(seq_len(m), seq_len(m), function(i, j) ifelse(i >= j & i-j < length(h), h[pmin(pmax(i-j+1, 1), length(h))], 0)) * dt
  B <- outer(seq_len(m) - 0.5, 0:(k-1L), function(t, f) cos(pi*t*f/m))
  B[, 1] <- B[, 1] / sqrt(m); B[, -1] <- B[, -1] * sqrt(2/m)
  scan <- seq(1, m, by = mt); A <- H[scan, ] %*% B
  beta <- solve(crossprod(A[keep, ]) + diag(lambda, k), crossprod(A[keep, ], seed[keep]))
  u <- drop(B %*% beta)
  expect_equal(d$neural_seed, u, tolerance = 1e-10)
  expect_equal(unname(d$X[, d$terms$type == "task"]), unname(H[scan, ] %*% p), tolerance = 1e-10)
  expected <- H[scan, ] %*% (u * sweep(p, 2, colMeans(p)))
  expect_equal(unname(d$X[, d$terms$type == "interaction"]), unname(expected), tolerance = 1e-10)
  seed[!keep] <- c(-999, NA)
  d2 <- ppi_design(seed, psych = list(values = p, grid = "neural"), TR = 2, microtime = mt,
                   hrf = h, censor = !keep, deconv = list(n_basis = k, lambda = lambda), seed_adjust="none")
  expect_equal(d$X, d2$X)
  expect_equal(d$keep, keep)
})

test_that("runs are isolated, events respect timing, and rank is retained", {
  runs <- rep(c("a", "b"), each = 4)
  events <- data.frame(onset = c(1, 2), duration = c(2, 0), condition = c("A", "B"), run = c("a", "b"))
  d <- ppi_design(1:8, events = events, runs = runs, TR = 2, microtime = 2,
                  hrf = c(0, 1, 0), deconv = list(n_basis = 3))
  d2 <- ppi_design(c(1:4, 99:102), events = events, runs = runs, TR = 2, microtime = 2,
                   hrf = c(0, 1, 0), deconv = list(n_basis = 3))
  expect_equal(d$X[1:4, ], d2$X[1:4, ])
  expect_equal(d$psych[1:8, "A"], c(0, 1, 1, 0, 0, 0, 0, 0))
  expect_equal(d$psych[9:16, "B"], c(0, 0, 1, 0, 0, 0, 0, 0))
  expect_true(d$provenance$rank < ncol(d$X))
  expect_equal(nrow(d$terms), ncol(d$X))
  expect_true(any(d$terms$condition == "B" & d$terms$run == "a", na.rm = TRUE))
  expect_true(all(d$X[, d$terms$condition == "B" & d$terms$run == "a" & !is.na(d$terms$condition), drop = FALSE] == 0))
  expect_error(ppi_design(1:4, psych = list(values = matrix(1, 4, 1), grid = "neural")), "row count")
  expect_error(ppi_design(1:4, psych = list(values = matrix(1, 4, 1), grid = "scan"), runs = c(1, 2, 1, 2)), "contiguous")
})

test_that("scan matrix is explicitly equivalent to declared scan input", {
  p <- cbind(A = c(0, 1, 2, 0, 1), B = c(1, 0, 1, 0, 1))
  args <- list(seed = c(2, 1, 3, 0, 4), TR = 1.5, microtime = 3,
               hrf = c(0, 0.5, 1, 0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               deconv = list(n_basis = 3), run_effects = "shared")
  x <- do.call(ppi_design, c(args, list(psych = p)))
  y <- do.call(ppi_design, c(args, list(psych = list(values = p, grid = "scan"))))
  expect_equal(x$X, y$X)
  expect_identical(x$provenance$psych_grid, "scan")
  h <- args$hrf; dt <- args$TR / args$microtime; m <- nrow(p) * args$microtime
  H <- outer(seq_len(m), seq_len(m), function(i, j) ifelse(i >= j & i-j < length(h), h[pmin(pmax(i-j+1, 1), length(h))], 0)) * dt
  neural <- p[rep(seq_len(nrow(p)), each = args$microtime), ]
  expect_equal(unname(x$X[, x$terms$type == "task"]), unname(H[seq(1, m, by = args$microtime), ] %*% neural), tolerance = 1e-10)
})
