test_that("NeuroVec extraction preserves voxel and time identities", {
  a <- array(seq_len(2 * 2 * 2 * 7), c(2, 2, 2, 7))
  space <- neuroim2::NeuroSpace(dim(a))
  x <- neuroim2::NeuroVec(a, space)
  out <- multiPPI:::.sl_as_matrix(x)
  expect_equal(dim(out$X), c(7L, 8L))
  expect_equal(out$X[3, 5], a[1, 1, 2, 3])
  expect_equal(out$X, t(matrix(as.numeric(a), nrow = 8)))
  expect_null(out$tr)
})

test_that("rest frame extraction uses current sampling metadata", {
  skip_if_not_installed("fmrireg")
  Y <- matrix(seq_len(80), 20, 4)
  frame <- fmrireg::matrix_frame(Y, TR = 2, run_length = c(8, 12))
  extracted <- multiPPI:::mppi_rest_extract(frame)
  expect_equal(unname(extracted$Y), unname(Y))
  expect_equal(as.character(extracted$runs), as.character(rep(1:2, c(8, 12))))
  expect_equal(extracted$TR, 2)
  expect_equal(multiPPI:::mppi_rest_extract(list(Y=Y,TR=3),runs=rep('runA',20))$runs, rep('runA',20))
})

test_that("derived covariance quantities do not depend on display diagonals", {
  set.seed(33)
  Y <- matrix(rnorm(300),100,3); X <- cbind(1,task=rnorm(100))
  a <- mppi_fit(Y,X,2,scale="cov",zero_diag=TRUE)
  b <- mppi_fit(Y,X,2,scale="cov",zero_diag=FALSE)
  expect_equal(a$partial,b$partial)
  expect_equal(a$variance,b$variance)
  expect_equal(diag(a$Delta_raw[[1]]),rep(0,3))
  expect_equal(a$covariance_diagonal[[1]],diag(b$Delta_raw[[1]]))
})

test_that("callable HRFs are evaluated at their actual time coordinates", {
  h <- multiPPI:::.mppi_hrf_to_vec(function(t) as.numeric(t==2),2)
  expect_equal(which(h!=0),2L)
  expect_error(multiPPI:::.mppi_hrf_to_vec(list(bad=TRUE),2),"HRF must")
})

test_that("advanced directed lags preserve original acquisition and run boundaries", {
  X <- matrix(1:10,5,2); p <- c(1,2,3,4,5)
  index <- c(1,2,4,5,6); runs <- c(1,1,1,2,2)
  fit <- structure(list(R=X,pk=list(p),names='task',runs=runs,
                        observation_index=index),class='mppi_fit')
  out <- mppi_lagged(fit,1,lags=c(-1,1),mode='neural')
  expected <- (p[1]*outer(X[1,],X[2,])+p[4]*outer(X[4,],X[5,]))/sum(p[c(1,4)]^2)
  expect_equal(out$M_lag[['1']],expected)
  expect_equal(out$M_lag[['-1']],t(expected))
  expect_equal(out$estimand,'directed_signal_lag_covariance')
  expect_error(mppi_lagged(fit,1,mode='innov'),"censor gaps")
})
