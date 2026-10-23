test_that("banded MAP matches an independent dense causal system", {
  dense <- function(y,h,lambda) {
    n <- nrow(y); H <- matrix(0,n,n)
    for (i in seq_len(n)) for (j in seq_len(i)) {
      k <- i-j+1L
      if (k <= length(h)) H[i,j] <- h[k]
    }
    D <- if (n == 1L) matrix(0,0,1) else diff(diag(n))
    solve(crossprod(H)+lambda*crossprod(D),crossprod(H,y))
  }
  set.seed(718)
  for (n in c(1L,2L,7L,19L,60L)) for (h in list(c(1), c(1,.2,0,-.1), c(0,1,0,.2), c(.8,0,.3,0,.1,.05,0,.02))) {
    y <- matrix(rnorm(n*3),n,3)
    # Leading-zero kernels identify n=1 only if H is nonzero; the
    # first-difference prior alone cannot identify a constant.
    if (n == 1L && h[1] == 0) next
    ref <- dense(y,h,2)
    got <- multiPPI:::mppi_deconv_map(y,h,2)
    expect_equal(unname(got),unname(ref),tolerance=1e-9,info=paste('n=',n,'h=',paste(h,collapse=',')))
  }
  y <- matrix(rnorm(36),12,3); h <- c(1,.1,0,.05)
  expect_equal(multiPPI:::mppi_deconv_map(y,h,0),dense(y,h,0),tolerance=1e-10)
  expect_equal(multiPPI:::mppi_deconv_map(matrix(6,1,1),c(2,0),4),matrix(3,1,1))
})

test_that("MAP rejects invalid and singular systems without jitter", {
  map <- multiPPI:::mppi_deconv_map; y <- matrix(1,5,2)
  expect_error(map(y,c(1,.3),-1),"finite non-negative lambda")
  expect_error(map(y,c(1,.3),Inf),"finite non-negative lambda")
  expect_error(map(y,c(1,NA),1),"finite data and HRF")
  y[2,1] <- NA; expect_error(map(y,c(1,.3),1),"finite data and HRF")
  expect_error(map(matrix(1,1,1),c(0,1),2),"singular")
  expect_error(map(matrix(1,5,1),c(0,1),0),"singular")
  expect_error(map(matrix(1,5,1),c(0,0),0.1),"singular")
  expect_error(map(matrix(numeric(),0,1),1,1),"at least one observation")
})
