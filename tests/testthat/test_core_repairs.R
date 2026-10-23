test_that("rank deficient designs retain intercept span and finite RSS", {
  set.seed(12); n <- 40; run <- rep(1:2, each=20)
  X <- cbind(run1=run==1, run2=run==2, task=sin(seq_len(n)), duplicate=run==1)
  Y <- matrix(rnorm(n*3), n)
  fit <- mppi_fit(Y, X, 3, runs=run)
  expect_equal(fit$evidence$df_full, qr(X)$rank)
  expect_equal(fit$evidence$rss_base, colSums(qr.resid(qr(X[, -3]), Y)^2))
  expect_equal(fit$evidence$rss_full, colSums(qr.resid(qr(X), Y)^2))
})

test_that("censored lag pairs use original acquisition distances and retain direction", {
  set.seed(13); n <- 50; X <- cbind(1, task=rnorm(n)); Y <- matrix(rnorm(n*3), n)
  Y[15, 1] <- NA
  fit <- mppi_fit(Y, X, 2, lags=1, packed=TRUE, zero_diag=FALSE, scale="normalized")
  keep <- setdiff(seq_len(n),15); a <- which(seq_len(n-1) %in% keep & (seq_len(n-1)+1) %in% keep)
  i <- match(a,keep); j <- match(a+1,keep)
  expected <- crossprod(fit$R[i,,drop=FALSE], fit$pk[[1]][i]*fit$R[j,,drop=FALSE])/sum(fit$pk[[1]][i]^2)
  expect_equal(unname(attr(fit$lagged[[1]][["1"]],"raw")), unname(expected))
  expect_true(is.matrix(fit$lagged[[1]][["1"]]))
  expect_error(mppi_fit(Y,X,2,domain="neural"), "masked deconvolution")
})

test_that("deconvolution preserves HRF timing", {
  set.seed(14); y <- matrix(rnorm(30),30)
  h <- c(0,0,1,.5,0)
  H <- outer(seq_len(30),seq_len(30),function(i,j) ifelse(i >= j & i-j < length(h), h[pmax(1,pmin(length(h),i-j+1))],0))
  D <- diff(diag(30))
  expected <- solve(crossprod(H)+2*crossprod(D),crossprod(H,y))
  expect_equal(mppi_deconv(y,TR=1,hrf=h,lambda=2), expected, tolerance=1e-9)
})

test_that("whitening dispatch supports bold and rejects unsupported order", {
  set.seed(15); n <- 80; X <- cbind(1,task=rnorm(n)); Y <- matrix(rnorm(n*3),n); runs <- rep(1:2,each=40)
  fit <- mppi(Y,X,psych_idx=2,runs=runs,prewhiten=TRUE,p=1,domain="bold")
  expect_equal(fit$domain,"whitened_bold")
  expect_error(mppi(Y,X,psych_idx=2,runs=runs,prewhiten=TRUE,p=1,domain="neural"),"unsupported")
})

test_that("current frame and event metadata propagate TR", {
  skip_if_not_installed("fmridesign"); skip_if_not_installed("fmrireg"); skip_if_not_installed("fmridataset")
  n <- 80; set.seed(16); Y <- matrix(rnorm(n*3), n)
  sf <- fmridesign::sampling_frame(blocklens=c(40,40),TR=2)
  events <- data.frame(onset=rep(c(4,16,28),2),run=rep(1:2,each=3),cond=factor(rep(c("A","B","A"),2)))
  hrf <- fmridesign::hrf
  em <- fmridesign::event_model(onset ~ hrf(cond),data=events,block=~run,sampling_frame=sf)
  bm <- fmridesign::baseline_model(basis="constant",sframe=sf)
  frame <- fmrireg::matrix_frame(Y,TR=2,run_length=c(40,40),event_table=events)
  design <- as_mppi_design(em,bm)
  expect_equal(qr(design$X)$rank,ncol(design$X))
  model <- fmrireg::fmri_model(em,bm,frame)
  expect_s3_class(mppi_fit_from_fmrireg(model,frame,psych=1:2,domain="bold"),"mppi_fit")
  expect_s3_class(mppi.fmri_frame(frame,design,domain="bold"),"mppi_fit")
  fit <- mppi(em,Y=Y,baseline_model=bm,domain="neural",deconv=list(method="fixed",K=10))
  expect_equal(fit$deconv$tr,2)
})

test_that("DCT deconvolution is run local", {
  set.seed(17); n <- 60; runs <- rep(1:2,each=30)
  X <- cbind(run1=runs==1,run2=runs==2,task=rnorm(n)); Y <- matrix(rnorm(n*2),n)
  h <- c(0,1,.4); cfg <- list(hrf=h, method="fixed", lambda=0.1, K=10, psych=matrix(X[,3],n))
  fit <- mppi_fit(Y,X,3,runs=runs,domain="neural",deconv=cfg)
  residual <- qr.resid(qr(X),Y)
  expected <- do.call(rbind,lapply(split(seq_len(n),runs),function(i) mppi_deconv_dct(residual[i,,drop=FALSE],h=h,K=10,method="fixed",lambda=0.1)$U))
  expect_equal(unname(fit$U),unname(expected))
})
