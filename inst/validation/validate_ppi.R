# Fixed-design conditional regression screen, not neural/deconvolution recovery.
# Predeclared: 2000 null and 500 effect replicates per scenario, alpha=.05;
# four primary null screens, Bonferroni one-sided 95% upper bound <= .075.
paths <- c('R/ppi_design.R','R/ppi_fit.R','inst/validation/validate_ppi.R')
start_hash <- tools::md5sum(paths)
writeLines(capture.output(print(start_hash)), 'audit/2026-10-23-ppi/implementation-validation-start-hashes.txt')
source('R/ppi_design.R'); source('R/ppi_fit.R')
set.seed(122031)
n <- 120L; runs <- rep(1:2, c(55,65)); seed <- as.numeric(arima.sim(list(ar=.4), n=n))
p <- cbind(A=rep(c(rep(0,5),rep(1,5)),12), B=rep(c(rep(1,7),rep(0,5)),10))
d <- ppi_design(seed, psych=p, runs=runs, microtime=2, hrf=c(0,1,.5,.2), run_effects='shared', nuisance=matrix(sin(seq_len(n)/9),ncol=1), deconv=list(n_basis=20))
stopifnot(qr(d$X)$rank == ncol(d$X))
censored <- d; censored$keep[c(12,13,45,80)] <- FALSE
# Calibration series independent of all evaluation draws; preserve full fitted plan.
set.seed(122032)
gen_ar <- function(phi, runs) {
  e <- numeric(length(runs))
  for (ids in split(seq_along(runs),runs)) {
    e[ids[1]] <- rnorm(1, sd=1/sqrt(1-phi^2))
    for (j in ids[-1]) e[j] <- phi*e[j-1]+rnorm(1)
  }
  e
}
cal <- replicate(200, gen_ar(.5,runs))
plan <- fmriAR::fit_noise(cal,runs=d$runs,p=1L,pooling='global')
stopifnot(inherits(plan,'fmriAR_plan'), identical(plan$method,'ar'), identical(plan$pooling,'global'), length(plan$phi)==1, length(plan$phi[[1]])==1, length(plan$theta[[1]])==0, abs(plan$phi[[1]])<1)
phi <- plan$phi[[1]]
idx <- which(d$terms$type=='interaction' & d$terms$condition=='A')
stopifnot(length(idx)==1)
beta <- rep(0,ncol(d$X)); beta[which(d$terms$type=='physiology')] <- .3; beta[which(d$terms$type=='task')] <- c(.3,-.2)
truth <- .7
scenarios <- c('iid','iid_censor','known_ar1','estimated_ar1','heterogeneous_ar1')
N0 <- 2000L; N1 <- 500L; alpha <- .05; simultaneous <- 1-.05/4
results <- list(); set.seed(122033)
for (scenario in scenarios) {
  result <- matrix(NA_real_,N0+N1,4,dimnames=list(NULL,c('p','estimate','se','covered')))
  for (i in seq_len(N0+N1)) {
    effect <- if(i<=N0) 0 else truth
    b <- beta; b[idx] <- effect
    if (scenario %in% c('iid','iid_censor')) err <- matrix(rnorm(n),ncol=1)
    else if (scenario == 'known_ar1') err <- matrix(gen_ar(phi,runs),ncol=1)
    else if (scenario == 'estimated_ar1') err <- matrix(gen_ar(.5,runs),ncol=1)
    else err <- sapply(c(.65,.1,.1,.1),gen_ar,runs=runs)
    Y <- sweep(err,1,drop(d$X %*% b),'+')
    design <- if(scenario=='iid_censor') censored else d
    noise <- if(scenario %in% c('iid','iid_censor')) 'iid' else if(scenario=='known_ar1') plan else 'fmriAR'
    fit <- if (is.character(noise) && noise=='fmriAR') ppi_fit(Y,design,noise=noise,noise_control=list(p=1L,pooling='global')) else ppi_fit(Y,design,noise=noise)
    est <- fit$coefficients[idx,1]; se <- fit$standard_errors[idx,1]
    pval <- 2*pt(abs(est/se),fit$df,lower.tail=FALSE)
    result[i,] <- c(pval,est,se,abs(est-effect)<=qt(.975,fit$df)*se)
  }
  null <- result[seq_len(N0),,drop=FALSE]; alternative <- result[N0+seq_len(N1),,drop=FALSE]
  hits <- sum(null[,'p']<alpha); upper <- qbeta(simultaneous,hits+1,N0-hits)
  interval <- binom.test(hits,N0,conf.level=.95)$conf.int
  results[[scenario]] <- data.frame(scenario=scenario,n_null=N0,rejections=hits,rate=hits/N0,ci_low=interval[1],ci_high=interval[2],simultaneous_upper=upper,screen=if(scenario=='heterogeneous_ar1') NA else upper<=.075,bias_null=mean(null[,'estimate']),coverage_null=mean(null[,'covered']),n_effect=N1,effect=truth,power=mean(alternative[,'p']<alpha),bias_effect=mean(alternative[,'estimate'])-truth,coverage_effect=mean(alternative[,'covered']))
  cat(scenario, 'complete\n'); flush.console()
}
stopifnot(identical(start_hash, tools::md5sum(paths)))
out <- 'audit/2026-10-23-ppi/implementation-validation-fixed-design.csv'
write.csv(do.call(rbind,results),out,row.names=FALSE)
writeLines(c(paste('seed design/calibration/simulation:',122031,122032,122033),paste('calibrated phi:',phi),paste('simultaneous primary upper confidence:',simultaneous),paste('R:',R.version.string),paste('fmriAR:',as.character(packageVersion('fmriAR'))),paste('fmrihrf:',as.character(packageVersion('fmrihrf'))),capture.output(tools::md5sum(c('R/ppi_design.R','R/ppi_fit.R','inst/validation/validate_ppi.R'))),'Primary: iid, iid_censor, known_ar1, estimated_ar1. Heterogeneous pooled AR is an adverse misspecification probe, not a primary acceptance scenario. Independent runs; one fixed two-run design; fixed-design conditional inference only. No neural recovery, HRF mismatch, selection or multiplicity qualification.'), 'audit/2026-10-23-ppi/implementation-validation-provenance.txt')
