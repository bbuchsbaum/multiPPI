args <- commandArgs(TRUE); pkgload::load_all(args[1],quiet=TRUE)
set.seed(948); h <- c(.8,.22,.1,-.04,0,.03,.01,0,-.01,.005,0,.002); lambda <- .7
for (n in c(256L,512L,1024L,2048L,4096L)) {
  Y <- matrix(rnorm(n*2L),n,2L)
  band_time <- system.time(band <- multiPPI:::mppi_deconv_map(Y,h,lambda))['elapsed']
  dense_time <- system.time({
    H <- matrix(0,n,n)
    for (k in seq_along(h)) {
      i <- seq.int(k,n); H[cbind(i,i-k+1L)] <- h[k]
    }
    A <- crossprod(H); diag(A) <- diag(A) + lambda*c(1,rep(2,n-2L),1)
    i <- seq_len(n-1L); A[cbind(i+1L,i)] <- A[cbind(i+1L,i)]-lambda; A[cbind(i,i+1L)] <- A[cbind(i,i+1L)]-lambda
    factor <- chol(A); dense <- backsolve(factor,forwardsolve(t(factor),crossprod(H,Y)))
  })['elapsed']
  err <- max(abs(band-dense)); if (err > 1e-8) stop('dense mismatch')
  cat(sprintf('n=%d band_sec=%.6f dense_sec=%.6f max_abs_error=%.4g\n',n,band_time,dense_time,err))
  rm(H,A,factor,dense,band,Y); gc(FALSE)
}
