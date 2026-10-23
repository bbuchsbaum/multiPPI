test_that("neural psych handles run-local and aliased columns", {
  set.seed(91); n <- 40; runs <- rep(1:2,each=20)
  X <- cbind(a=rnorm(n),b=rnorm(n)); h <- c(0,1,.4)
  whole <- mppi_psych_neural_from_X(X,1:2,h,K=8,method="fixed",runs=runs)
  pieces <- lapply(split(seq_len(n),runs),function(i) mppi_psych_neural_from_X(X[i,,drop=FALSE],1:2,h,K=8,method="fixed")$S)
  expect_equal(whole$S,do.call(rbind,pieces))
  X[,2] <- X[,1]; aliased <- mppi_psych_neural_from_X(X,1:2,h,K=8,method="fixed")
  expect_true(all(is.finite(aliased$S)))
})

test_that("HRF ensemble weights and caches agree", {
  set.seed(92); n <- 40; Y <- matrix(rnorm(n*3),n); X <- cbind(1,task=rnorm(n))
  f1 <- mppi_fit(Y,X,2,domain="neural",deconv=list(hrf=c(0,1,.2),K=10,method="fixed"))
  f2 <- mppi_fit(Y,X,2,domain="neural",deconv=list(hrf=c(0,1,.7),K=10,method="fixed"))
  e <- mppi_hrf_ensemble(list(f1,f2)); expect_equal(unname(e$ensemble$weights),c(.5,.5))
  for (mode in c("raw","amplitude","normalized"))
    expect_equal(mppi_get_M_scaled(e,1,mode),.5*mppi_get_M_scaled(f1,1,mode)+.5*mppi_get_M_scaled(f2,1,mode))
  expect_equal(mppi_get_M_lag(e,1,0),mppi_get_M(e,1))
  expect_error(mppi_hrf_ensemble(list(f1,f2),weights=1),"one finite")
  expect_error(mppi_hrf_ensemble(list(f1,f2),weights=c(1,Inf)),"one finite")
  f2$names <- "different"; expect_error(mppi_hrf_ensemble(list(f1,f2)),"matching domain")
})

test_that("evidence is task main effects and cross-domain AIC is rejected", {
  set.seed(93); Y <- matrix(rnorm(120),40); X <- cbind(1,task=rnorm(40))
  b <- mppi_fit(Y,X,2); n <- mppi_fit(Y,X,2,domain="neural",deconv=list(K=10,method="fixed"))
  expect_equal(mppi_evidence(b)$estimand,"task_main_effect_baseline_glm")
  expect_equal(attr(mppi_evidence(b,FALSE),"estimand"),"task_main_effect_baseline_glm")
  expect_true(is.numeric(mppi_compare_models(b,n)$delta_gap))
  expect_error(mppi_compare_models(b,n,resid_bold=b$R,resid_neural=n$R,pk=b$pk),"AIC comparisons")
})

test_that("lag-selection null repeats the same maximum search", {
  set.seed(94); n <- 30; X <- cbind(1,task=rnorm(n)); Y <- matrix(rnorm(n*3),n); runs <- rep(1:2,each=15); lags <- -1:1
  set.seed(99); got <- mppi_lag_select(Y,X,2,lags=lags,blocklens=c(15,15),blksize=5,B=7,scale="corr",zero_diag=FALSE)
  fit <- mppi_fit(Y,X,2,runs=runs,scale="corr",zero_diag=FALSE); R <- sweep(fit$R,2,apply(fit$R,2,sd),"/"); pk <- fit$pk[[1]]
  blocks <- split(seq_len(n),rep(seq_len(6),each=5))
  stat <- function(x) vapply(lags,function(lg) {z<-multiPPI:::.mppi_shift_by_run(x,lg,c(15,15)); sum((crossprod(R,z*R)/sum(z^2))^2)},numeric(1))
  set.seed(99); expected <- replicate(7,{w<-sample(c(-1,1),6,replace=TRUE); z<-pk; for(j in seq_along(blocks)) z[blocks[[j]]]<-z[blocks[[j]]]*w[j]; max(stat(z))})
  expect_equal(unname(got$Qobs),stat(pk)); expect_equal(got$Qnull,expected)
  expect_equal(got$inference$status,"experimental")
  shifts <- mppi_lagged(Y,X,2,lags=-1:1,blocklens=c(15,15))
  expect_equal(attr(shifts,"lag_estimand"),"psychological_shift")
})

test_that("neural conversion preserves diagonal policy and validates its inputs", {
  set.seed(95); n <- 40; Y <- matrix(rnorm(n*3),n); X <- cbind(1,task=rnorm(n)); runs <- rep(1:2,each=20); h <- c(0,1,.3)
  bold <- mppi_fit(Y,X,2,runs=runs,zero_diag=FALSE)
  converted <- mppi_neural_from_fit(bold,X,2,h,K=8,method="fixed")
  direct <- mppi_fit(Y,X,2,runs=runs,zero_diag=FALSE,domain="neural",deconv=list(hrf=h,K=8,method="fixed"))
  expect_equal(converted$U,direct$U,tolerance=1e-8)
  expect_equal(mppi_get_M_scaled(converted,1,"raw"),mppi_get_M_scaled(direct,1,"raw"),tolerance=1e-8)
  altered <- X; altered[1,2] <- altered[1,2]+1
  expect_error(mppi_neural_from_fit(bold,altered,2,h,K=8),"do not match")
  Y[4,1] <- NA; censored <- mppi_fit(Y,X,2,runs=runs)
  expect_error(mppi_neural_from_fit(censored,X,2,h,K=8),"Censored")
})
