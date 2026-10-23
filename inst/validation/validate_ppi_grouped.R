# Supplemental targetwise AR screen, predeclared before execution.
# 2000 null replicates; alpha .05; one-sided 95% exact upper bound <= .075.
paths <- c('R/ppi_design.R','R/ppi_fit.R','inst/validation/validate_ppi_grouped.R')
hash <- tools::md5sum(paths)
writeLines(capture.output(print(hash)),'audit/2026-10-23-ppi/implementation-validation-grouped-start-hashes.txt')
source('R/ppi_design.R'); source('R/ppi_fit.R')
set.seed(122031); n <- 120L; runs <- rep(1:2,c(55,65))
seed <- as.numeric(arima.sim(list(ar=.4),n=n)); p <- cbind(A=rep(c(rep(0,5),rep(1,5)),12),B=rep(c(rep(1,7),rep(0,5)),10))
d <- ppi_design(seed,psych=p,runs=runs,microtime=2,hrf=c(0,1,.5,.2),run_effects='shared',nuisance=matrix(sin(seq_len(n)/9),ncol=1),deconv=list(n_basis=20))
idx <- which(d$terms$type=='interaction' & d$terms$condition=='A'); stopifnot(length(idx)==1, qr(d$X)$rank==ncol(d$X))
b <- rep(0,ncol(d$X)); b[which(d$terms$type=='physiology')] <- .3; b[which(d$terms$type=='task')] <- c(.3,-.2)
gen <- function(phi) { e <- numeric(n); for(ids in split(seq_len(n),runs)) {e[ids[1]] <- rnorm(1,sd=1/sqrt(1-phi^2)); for(j in ids[-1]) e[j] <- phi*e[j-1]+rnorm(1)}; e }
N <- 2000L; set.seed(122036); out <- matrix(NA_real_,N,3)
for (i in seq_len(N)) {
  Y <- sweep(sapply(c(.65,.1,.1,.1),gen),1,drop(d$X %*% b),'+')
  fit <- ppi_fit(Y,d,noise_control=list(pooling='parcel',parcels=1:4,p=1L))
  e <- fit$coefficients[idx,1]; se <- fit$standard_errors[idx,1]; df <- fit$df[1]
  out[i,] <- c(2*pt(abs(e/se),df,lower.tail=FALSE),e,abs(e)<=qt(.975,df)*se)
  if(i%%250L==0L) {cat(i,'complete\n');flush.console()}
}
stopifnot(identical(hash,tools::md5sum(paths)))
hits <- sum(out[,1]<.05); ci <- binom.test(hits,N)$conf.int; upper <- qbeta(.95,hits+1,N-hits)
write.csv(data.frame(scenario='heterogeneous_targetwise_ar1',n=N,rejections=hits,rate=hits/N,ci_low=ci[1],ci_high=ci[2],upper_95=upper,screen=upper<=.075,bias=mean(out[,2]),coverage=mean(out[,3])), 'audit/2026-10-23-ppi/implementation-validation-grouped.csv',row.names=FALSE)
writeLines(c('seed 122031 design; seed 122036 simulation; one-sided 95% upper <= .075; supplemental screen.',capture.output(print(hash)),R.version.string,paste('fmriAR',packageVersion('fmriAR'))),'audit/2026-10-23-ppi/implementation-validation-grouped-provenance.txt')
