test_that("unseeded randomization advances the caller's stream; seeded calls restore it", {
  set.seed(82); n <- 40; X <- cbind(1, rnorm(n)); Y <- matrix(rnorm(n * 3), n)
  fit <- mppi_fit(Y, X, 2, runs = rep(1:2, each = 20))
  a <- mppi_omnibus(fit, blksize = 5, B = 20); b <- mppi_omnibus(fit, blksize = 5, B = 20)
  expect_false(identical(a[[1]]$Qnull, b[[1]]$Qnull))
  p1 <- mppi_permute(fit, blksize = 5, B = 20); p2 <- mppi_permute(fit, blksize = 5, B = 20)
  expect_false(identical(attr(p1, "rng")$seed, attr(p2, "rng")$seed))
  set.seed(1); before <- .Random.seed
  s1 <- mppi_omnibus(fit, blksize = 5, B = 20, seed = 9)
  expect_identical(.Random.seed, before)
  expect_identical(s1[[1]]$Qnull, mppi_omnibus(fit, blksize = 5, B = 20, seed = 9)[[1]]$Qnull)
})

test_that("mppi_group_glm accepts model.matrix designs and matches named rows", {
  fits <- lapply(0:3, function(s) {
    set.seed(100 + s)
    mppi_fit(Y = matrix(rnorm(64), 16, 4), X = cbind(1, seq_len(16)), psych_idx = 2)
  })
  names(fits) <- paste0("subj", 1:4)
  grp <- data.frame(g = factor(c("a", "a", "b", "b")))
  ref <- mppi_group_glm(fits, design = unname(model.matrix(~ g, grp)), method = "lm")
  expect_equal(mppi_group_glm(fits, design = model.matrix(~ g, grp), method = "lm")$estimate, ref$estimate)
  named <- model.matrix(~ g, grp); rownames(named) <- paste0("subj", 4:1)
  expect_equal(mppi_group_glm(fits, design = named, method = "lm")$estimate, -ref$estimate)
})

test_that("simulated rest data are defined from the first scan", {
  h <- c(0, 1, 0.5)
  sim <- mppi_sim_rest_dataset(T = 60L, V = 4L, r = 2L, hrf = h, obs_noise = 0, seed = 3)
  expect_true(all(is.finite(sim$Y)))
  expect_equal(nrow(sim$Y), 60L)
  conv <- apply(sim$U, 2, function(x) as.numeric(stats::filter(c(0, 0, x), h, sides = 1))[-(1:2)])
  expect_equal(qr.resid(qr(conv), sim$Y), matrix(0, 60, 4), tolerance = 1e-10)
})

test_that("seed CCA rejects neighborhoods outside the data instead of dropping indices", {
  set.seed(83); n <- 80; Y <- matrix(rnorm(n * 8), n); D <- cbind(1, rnorm(n), rnorm(n))
  args <- list(data = Y, design = D, seed = 1:2, runs = rep(1:2, each = 40), domain = "hemo",
               prewhiten = FALSE, inference = list(method = "none"),
               pca = list(seed_k = 2, light_k = 3, center = TRUE, scale = TRUE))
  expect_error(do.call(seed_ppi_mcca, c(args, list(neighborhoods = list(4:9)))), "outside")
  expect_error(do.call(seed_ppi_mcca, c(args, list(neighborhoods = list(c(4, 5.5, 6))))), "integer")
  expect_silent(do.call(seed_ppi_mcca, c(args, list(neighborhoods = list(4:8)))))
})

test_that("beta permutation advances an unseeded stream and restores a seeded one", {
  set.seed(91); n <- 40; B <- matrix(rnorm(n * 4), n); P <- cbind(cont = rnorm(n))
  a <- mppi_beta_permute(B, P, Bperm = 19); b <- mppi_beta_permute(B, P, Bperm = 19)
  expect_false(identical(a[[1]]$Qnull, b[[1]]$Qnull))
  set.seed(2); before <- .Random.seed
  s1 <- mppi_beta_permute(B, P, Bperm = 19, seed = 5)
  expect_identical(.Random.seed, before)
  expect_identical(s1[[1]]$Qnull, mppi_beta_permute(B, P, Bperm = 19, seed = 5)[[1]]$Qnull)
})

test_that("neural fits deconvolve with the event model's own HRF", {
  skip_if_not_installed("fmridesign")
  hrf <- fmridesign::hrf
  sf <- fmridesign::sampling_frame(c(60, 60), TR = 2)
  ev <- data.frame(onset = rep(seq(4, 100, by = 16), 2), cond = factor(rep(c("A", "B"), 7)), run = rep(1:2, each = 7))
  m <- fmridesign::event_model(onset ~ hrf(cond, basis = "gaussian"), data = ev, block = ~run, sampling_frame = sf)
  h <- attr(stats::terms(m)[[1]], "hrfspec")$hrf
  expect_equal(mppi_hrf(m)(seq(0, 20, by = 2)), h(seq(0, 20, by = 2)))
  set.seed(4); Y <- matrix(rnorm(120 * 3), 120)
  d <- as_mppi_design(event = m)
  f <- mppi(d, Y = Y, domain = "neural", deconv = list(type = "map"))
  ref <- mppi(d, Y = Y, domain = "neural", deconv = list(type = "map", hrf = h))
  expect_equal(f$Delta, ref$Delta)
  m3 <- fmridesign::event_model(onset ~ hrf(cond, basis = "spmg3"), data = ev, block = ~run, sampling_frame = sf)
  expect_silent(mppi(as_mppi_design(event = m3), Y = Y, domain = "neural", deconv = list(type = "map")))
})

test_that("models with different term HRFs need deconv$hrf only for deconvolution", {
  skip_if_not_installed("fmridesign")
  hrf <- fmridesign::hrf
  sf <- fmridesign::sampling_frame(c(60, 60), TR = 2)
  ev <- data.frame(onset = rep(seq(4, 100, by = 16), 2), a = factor(rep(c("x", "y"), 7)),
                   b = factor(rep(c("u", "v"), each = 7)), run = rep(1:2, each = 7))
  m <- fmridesign::event_model(onset ~ hrf(a) + hrf(b, basis = "gaussian"), data = ev, block = ~run, sampling_frame = sf)
  expect_null(mppi_hrf(m))
  d <- as_mppi_design(event = m)
  set.seed(4); Y <- matrix(rnorm(120 * 3), 120)
  expect_silent(mppi(d, Y = Y))
  expect_error(mppi(d, Y = Y, domain = "neural", deconv = list(type = "map")), "different HRFs")
  h <- attr(stats::terms(m)[[1]], "hrfspec")$hrf
  expect_silent(mppi(d, Y = Y, domain = "neural", deconv = list(type = "map", hrf = h)))
})
