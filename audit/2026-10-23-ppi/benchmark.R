args <- commandArgs(TRUE)
pkgload::load_all(args[1], quiet=TRUE)
out <- args[2]
set.seed(20261023)
Tn <- 600L; V <- 400L; K <- 3L
Y <- matrix(rnorm(Tn*V),Tn,V)
X <- cbind(Intercept=1,matrix(rnorm(Tn*K),Tn,K))
colnames(X)[-1] <- paste0("task",seq_len(K))
basis <- as_mppi_basis(qr.Q(qr(matrix(rnorm(V*25),V,25))))
cases <- list(
  normalized_blas=list(),
  normalized_chunked=list(backend="chunked",chunk_size=150L),
  normalized_packed=list(packed=TRUE),
  covariance=list(scale="cov"),
  basis_rank25=list(basis=basis),
  neural_DCT=list(domain="neural",deconv=list(K=48L,method="fixed")))
summaries <- list()
for (name in names(cases)) {
  pars <- c(list(Y=Y,X=X,psych_idx=2:4),cases[[name]])
  warm <- do.call(mppi_fit,pars)
  times <- numeric(3)
  for(i in seq_along(times)) {
    invisible(gc())
    times[i] <- system.time(ans <- do.call(mppi_fit,pars))[["elapsed"]]
  }
  memfile <- file.path(out,paste0("alloc-",name,".txt"))
  Rprofmem(memfile)
  ans <- do.call(mppi_fit,pars)
  Rprofmem(NULL)
  lines <- readLines(memfile)
  bytes <- suppressWarnings(as.numeric(sub(" .*","",lines)))
  sum_alloc <- sum(bytes,na.rm=TRUE)
  unlink(memfile)
  serial <- file.path(out,paste0("size-",name,".rds"))
  saveRDS(ans,serial,compress=FALSE)
  serialized_bytes <- file.info(serial)$size
  unlink(serial)
  retained <- if(requireNamespace("lobstr",quietly=TRUE)) as.numeric(lobstr::obj_size(ans)) else NA_real_
  summaries[[name]] <- data.frame(case=name,T=Tn,V=V,K=K,min_s=min(times),median_s=median(times),max_s=max(times),
    R_allocated_MiB=sum_alloc/2^20,retained_shared_MiB=retained/2^20,uncompressed_RDS_MiB=serialized_bytes/2^20)
  print(summaries[[name]],row.names=FALSE)
}
result <- do.call(rbind,summaries)
write.csv(result,file.path(out,"benchmark.csv"),row.names=FALSE)

# Equivalent-output kernel comparisons; no reduction/basis approximation.
set.seed(20261024)
R <- matrix(rnorm(300*80),300,80); pp <- matrix(rnorm(300*3),300,3)
fn <- function(s) get(s,asNamespace("multiPPI"))
backends <- c("blas","chunked","accumulate")
ref <- fn(".mppi_crossprod")(R,pp[,1],backend="blas")
for(b in backends) {
  m <- fn(".mppi_crossprod")(R,pp[,1],backend=b,chunk_size=75)
  tm <- replicate(3,system.time(fn(".mppi_crossprod")(R,pp[,1],backend=b,chunk_size=75))[["elapsed"]])
  cat("kernel",b,"max_error",max(abs(m-ref)),"seconds",paste(tm,collapse=","),"\n")
}
batch <- fn("ppi_batch")(R,pp,lags=c(-1,0,1),blocklens=c(140,160))
oracle_errors <- c()
for(li in seq_along(c(-1,0,1))) for(k in 1:3) {
  lag <- c(-1,0,1)[li]; accum <- matrix(0,80,80); den <- 0
  for (ids in list(1:140,141:300)) {
    times <- if (lag==0) ids else head(ids,-1)
    for(t in times) {
      left <- if(lag < 0) t+1 else t
      right <- if(lag > 0) t+1 else t
      accum <- accum + pp[t,k]*outer(R[left,],R[right,])
      den <- den + pp[t,k]^2
    }
  }
  oracle_errors <- c(oracle_errors,max(abs(batch$matrices[[li]][,,k]-accum/den)))
}
cat("native_batch_independent_outer_product_max_error",max(oracle_errors),"\n")
capture.output(sessionInfo(),file=file.path(out,"benchmark-session.txt"))
