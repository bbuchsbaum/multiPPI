test_that("seed nuisance adjustment protects task effects and removes arbitrary nuisance additions", {
  set.seed(702)
  n <- 90; rr <- rep(1:2,c(40,50)); p <- cbind(A=rnorm(n))
  nuisance <- cbind(motion=.5*p[,1]+rnorm(n))
  s <- 2*p[,1]+rnorm(n); mask <- seq_len(n) %in% c(8,9,70)
  build <- function(seed, ...) ppi_design(seed,psych=p,runs=rr,nuisance=nuisance,
    censor=mask,hrf=c(0,1,.5),microtime=1,deconv=list(n_basis=12),...)
  a <- build(s)
  b <- build(s+ifelse(rr==1,1000,-400)+ifelse(rr==1,3,-2)*nuisance[,1])
  expect_equal(a$seed_adjusted,b$seed_adjusted,tolerance=1e-9)
  expect_equal(a$neural_seed,b$neural_seed,tolerance=1e-9)
  expect_equal(a$X,b$X,tolerance=1e-9)
  for (r in 1:2) {
    idx <- which(rr==r & !mask)
    task <- a$X[idx,a$terms$type=='task' & a$terms$run==r,drop=FALSE]
    ref <- lm(s[idx] ~ scale(task,scale=FALSE) + nuisance[idx,1])
    removed <- coef(ref)[1]+coef(ref)[3]*nuisance[idx,1]
    expect_equal(a$seed_adjusted[idx],s[idx]-removed,tolerance=1e-9)
  }
  expect_gt(max(abs(build(s,seed_adjust='none')$X-build(s+1000,seed_adjust='none')$X)),1)
})

test_that("baseline overlap is removed but remaining nuisance/task overlap errors", {
  set.seed(90)
  p <- cbind(A=rep(c(0,1),20),B=rep(c(1,0),20)); s <- rnorm(40)
  d <- ppi_design(s,psych=p,hrf=1,microtime=1)
  expect_equal(mean(d$seed_adjusted),0,tolerance=1e-10)
  expect_error(ppi_design(s,psych=p,nuisance=p[,1,drop=FALSE],hrf=1,microtime=1),"unidentified")
})

test_that("contrast estimability and joint tests are invariant to row scaling", {
  set.seed(882)
  d <- ppi_design(rnorm(80),psych=cbind(A=rnorm(80),B=rnorm(80)),
    nuisance=cbind(duplicate=rep(1,80)),hrf=1,microtime=1)
  f <- ppi_fit(matrix(rnorm(160),80),d,noise='iid')
  for (s in c(1,1e-9,1e9)) {
    z <- ppi_contrast(f,c('intercept:1'=s))
    expect_false(z$estimable)
    expect_true(all(is.na(z$statistic)))
  }
  L <- structure(diag(2),dimnames=list(NULL,c('A','B')))
  a <- ppi_contrast(f,L,test='F'); L[2,] <- L[2,]*1e-9
  b <- ppi_contrast(f,L,test='F')
  expect_equal(a$numerator_df,2)
  expect_equal(a$numerator_df,b$numerator_df)
  expect_equal(a$statistic,b$statistic,tolerance=1e-10)
})

test_that("predictor units do not alter rank or transformed inference", {
  set.seed(701)
  d <- ppi_design(rnorm(80),psych=cbind(A=rnorm(80),B=rnorm(80)),hrf=1,microtime=1)
  y <- matrix(rnorm(160),80); a <- ppi_fit(y,d,noise='iid')
  idx <- which(d$terms$type=='interaction' & d$terms$condition=='B')
  d$X[,idx] <- d$X[,idx]*1e8
  b <- ppi_fit(y,d,noise='iid')
  expect_equal(a$rank,b$rank)
  expect_equal(b$coefficients[idx,]*1e8,a$coefficients[idx,],tolerance=1e-9)
  expect_equal(b$sigma2,a$sigma2,tolerance=1e-10)
  expect_equal(ppi_contrast(a,c(B=1))$statistic,ppi_contrast(b,c(B=1e8))$statistic,tolerance=1e-9)
})

test_that("deconvolution does not depend on HRF amplitude or microtime", {
  set.seed(7); n <- 120
  spm <- function(t) ifelse(t < 0, 0, dgamma(t, 6, 1) - dgamma(t, 16, 1) / 6)
  seed <- as.numeric(stats::filter(rnorm(n), 0.7, method = "recursive"))
  ps <- cbind(A = rep(rep(0:1, each = 4), length.out = n))
  d <- ppi_design(seed, psych = ps, TR = 2, hrf = spm, run_effects = "shared")
  d10 <- ppi_design(seed, psych = ps, TR = 2, hrf = function(t) 10 * spm(t), run_effects = "shared")
  expect_equal(d10$X[, "interaction:A"], d$X[, "interaction:A"], tolerance = 1e-10)
  expect_equal(d10$neural_seed, d$neural_seed / 10, tolerance = 1e-10)
  d64 <- ppi_design(seed, psych = ps, TR = 2, hrf = spm, microtime = 64, run_effects = "shared")
  expect_gt(cor(d64$X[, "interaction:A"], d$X[, "interaction:A"]), 0.999)
  expect_equal(sd(d64$X[, "interaction:A"]), sd(d$X[, "interaction:A"]), tolerance = 0.01)
})
