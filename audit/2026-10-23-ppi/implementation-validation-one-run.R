# Supplemental one-run fixed-design IID check; not part of the four predeclared primary screens.
paths <- c('R/ppi_design.R','R/ppi_fit.R','audit/2026-10-23-ppi/implementation-validation-one-run.R')
hashes <- tools::md5sum(paths)
source('R/ppi_design.R'); source('R/ppi_fit.R')
set.seed(122035); n <- 120L
seed <- as.numeric(arima.sim(list(ar=.4), n=n))
p <- cbind(A=rep(c(rep(0,5),rep(1,5)),12), B=rep(c(rep(1,7),rep(0,5)),10))
d <- ppi_design(seed,psych=p,microtime=2,hrf=c(0,1,.5,.2),run_effects='shared',nuisance=matrix(sin(seq_len(n)/9),ncol=1),deconv=list(n_basis=20))
idx <- which(d$terms$type=='interaction' & d$terms$condition=='A')
stopifnot(length(idx)==1, qr(d$X)$rank==ncol(d$X))
b <- rep(0,ncol(d$X)); b[which(d$terms$type=='physiology')] <- .3; b[which(d$terms$type=='task')] <- c(.3,-.2)
N0 <- 2000L; N1 <- 500L; truth <- .7
Y <- matrix(rnorm(n*(N0+N1)),n) + drop(d$X %*% b)
Y[,N0+seq_len(N1)] <- Y[,N0+seq_len(N1)] + truth*d$X[,idx]
f <- ppi_fit(Y,d,noise='iid'); est <- f$coefficients[idx,]; se <- f$standard_errors[idx,]; pv <- 2*pt(abs(est/se),f$df,lower.tail=FALSE)
hits <- sum(pv[seq_len(N0)]<.05); ci <- binom.test(hits,N0)$conf.int
z <- data.frame(scenario='one_run_iid_supplemental',n_null=N0,rejections=hits,rate=hits/N0,ci_low=ci[1],ci_high=ci[2],bias_null=mean(est[seq_len(N0)]),coverage_null=mean(abs(est[seq_len(N0)])<=qt(.975,f$df)*se[seq_len(N0)]),n_effect=N1,effect=truth,power=mean(pv[N0+seq_len(N1)]<.05),bias_effect=mean(est[N0+seq_len(N1)])-truth,coverage_effect=mean(abs(est[N0+seq_len(N1)]-truth)<=qt(.975,f$df)*se[N0+seq_len(N1)]))
stopifnot(identical(hashes,tools::md5sum(paths)))
write.csv(z,'audit/2026-10-23-ppi/implementation-validation-one-run.csv',row.names=FALSE)
writeLines(c('Seed 122035; supplemental post-primary fixture; conditional fixed design.',capture.output(print(hashes)),R.version.string), 'audit/2026-10-23-ppi/implementation-validation-one-run-provenance.txt')
