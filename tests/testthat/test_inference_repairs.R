test_that("experimental weights have correct moments and boundaries", {
  set.seed(81); w <- multiPPI:::.mppi_random_weights(200000,"mammen")
  expect_lt(abs(mean(w)), .01); expect_lt(abs(mean(w^2)-1), .01)
  fit <- list(R=matrix(rnorm(24),12,2), runs=rep(1:2,each=6), observation_index=c(1:3,5:7,8:13), dropped=1L)
  setup <- multiPPI:::.mppi_randomization_setup(fit, 4, 10, "block_flip", "rademacher")
  expect_equal(unname(lengths(setup$segments)),c(3L,3L,6L))
  expect_true(all(vapply(setup$blocks, function(i) length(unique(fit$runs[i]))==1 && all(diff(fit$observation_index[i])==1), logical(1))))
})

test_that("randomization streams counts and uses matching observed and null statistics", {
  set.seed(82); n <- 40; X <- cbind(1,rnorm(n)); Y <- matrix(rnorm(n*3),n)
  fit <- mppi_fit(Y,X,2,runs=rep(1:2,each=20))
  om <- mppi_omnibus(fit,blksize=5,B=30,seed=10)
  expect_gt(length(unique(om[[1]]$Qnull)),1)
  expect_equal(attr(om,"inference")$status,"experimental")
  setup <- multiPPI:::.mppi_randomization_setup(fit,5,30,"block_flip","rademacher")
  set.seed(10); pk <- fit$pk[[1]]; D <- crossprod(fit$R,pk*fit$R)/sum(pk^2); diag(D)<-0
  counts <- matrix(0,3,3); qs <- numeric(30)
  for (b in seq_len(30)) {
    z <- multiPPI:::.mppi_random_psych(pk,setup,"block_flip","rademacher")
    M <- crossprod(fit$R,z*fit$R)/sum(z^2); diag(M)<-0
    counts <- counts+(abs(M)>=abs(D)); qs[b] <- sum(M^2)
  }
  expect_equal(om[[1]]$Qnull,qs)
  p <- mppi_permute(fit,blksize=5,B=30,seed=10)
  expected <- (1+counts)/31; diag(expected)<-NA_real_
  expect_equal(unname(as.vector(p)),unname(as.vector(expected)))
  expect_error(mppi_permute(fit,studentize=TRUE),"Studentization")
  expect_error(mppi_omnibus(fit,method="freedman_lane"),"reduced product")
  expect_error(mppi_permute(fit,wild="none"),"nonconstant")
})

test_that("seed CCA supports descriptive raw seed equivalence and rejects inference", {
  set.seed(83); n <- 80; Y <- matrix(rnorm(n*8),n); D <- cbind(1,rnorm(n),rnorm(n)); D <- cbind(D,D[,2]); rr <- rep(1:2,each=40)
  args <- list(data=Y,design=D, neighborhoods=list(4:8),runs=rr,domain="hemo",prewhiten=FALSE,inference=list(method="none"),pca=list(seed_k=2,light_k=3,center=TRUE,scale=TRUE))
  a <- do.call(seed_ppi_mcca,c(args,list(seed=1:2)))
  b <- do.call(seed_ppi_mcca,c(args,list(seed=Y[,1:2,drop=FALSE])))
  expect_equal(a$stat,b$stat,tolerance=1e-10)
  args$prewhiten <- TRUE; args$ar_order <- 1L
  a <- do.call(seed_ppi_mcca,c(args,list(seed=1:2)))
  b <- do.call(seed_ppi_mcca,c(args,list(seed=Y[,1:2,drop=FALSE])))
  expect_equal(a$stat,b$stat,tolerance=1e-10)
  args$inference <- list(method="param")
  expect_error(do.call(seed_ppi_mcca,c(args,list(seed=1:2))),"inference is unavailable")
  args$inference <- list(method="perm")
  expect_error(do.call(seed_ppi_mcca,c(args,list(seed=1:2))),"inference is unavailable")
  args$inference <- list(method="none"); args$domain <- "neural"
  expect_error(do.call(seed_ppi_mcca,c(args,list(seed=1:2))),"neural deconvolution is unsupported")
})

test_that("classical CCA helper uses the all-roots Bartlett test", {
  r <- c(.4,.2); chi <- -(99-(2+2+1)/2)*sum(log1p(-r^2))
  expect_equal(multiPPI:::.cca_pvalue_param(r,100,2,2),pchisq(chi,df=4,lower.tail=FALSE))
  expect_error(multiPPI:::.cca_pvalue_param(r,100,2,2,k=1),"all canonical")
})

test_that("phase draws preserve each contiguous segment spectrum", {
  set.seed(84); fit <- list(R=matrix(rnorm(36),18,2), runs=rep(1:2,each=9), observation_index=c(1:4,6:10,11:19), dropped=1L)
  setup <- multiPPI:::.mppi_randomization_setup(fit,4,10,"phase","rademacher")
  x <- rnorm(18); z <- multiPPI:::.mppi_random_psych(x,setup,"phase","rademacher")
  for (i in setup$segments) expect_equal(Mod(fft(z[i])),Mod(fft(x[i])),tolerance=1e-10)
})
