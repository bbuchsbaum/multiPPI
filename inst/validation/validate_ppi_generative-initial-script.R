# Frozen construction screen; run from repository root. Optional first argument is
# number of independent neural seeds (default 200); second is output prefix.
args <- commandArgs(TRUE)
ns <- if (length(args)) as.integer(args[1]) else 200L
out <- if (length(args) > 1) args[2] else 'audit/2026-10-23-ppi/implementation-generative'
stopifnot(ns > 0)
source('R/ppi_design.R'); source('R/ppi_fit.R')
set.seed(20261023)
dt <- .5; micro <- 4L; sizes <- c(90L,110L); nr <- 10L
h <- dgamma(seq(0,32,by=dt),6)-dgamma(seq(0,32,by=dt),16)/6
h <- h/(sum(h)*dt)
hdelay <- c(0,0,h)[seq_along(h)]
# Independent dense lower-triangular causal operators, sampled at first bin.
operator <- function(n,h) {
  m <- n*micro; ix <- outer(seq(1,m,by=micro),seq_len(m),'-')
  ans <- matrix(0,n,m); ok <- ix >= 0 & ix < length(h)
  ans[ok] <- h[ix[ok]+1L]*dt; ans
}
H <- lapply(sizes,operator,h=h); Hd <- lapply(sizes,operator,h=hdelay)
psych <- lapply(seq_along(sizes), function(r) {
  t <- (seq_len(sizes[r]*micro)-1)*dt
  phase <- (t-if(r==2) 10 else 0) %% 40
  cbind(A=as.numeric(phase<10), B=as.numeric(phase>=20 & phase<30))
})
runs <- rep(1:2,sizes); Z <- unlist(lapply(sizes,function(n) sin(2*pi*(seq_len(n)-1)/37)))
P <- do.call(rbind,psych); rows <- split(seq_along(runs),runs)
phi <- exp(-dt/2)
records <- vector('list',ns*4L); at <- 0L
for (s in seq_len(ns)) {
  u <- lapply(sizes,function(n) as.numeric(arima.sim(list(ar=phi),n=n*micro,sd=sqrt(1-phi^2),n.start=500)))
  sn <- lapply(seq_along(sizes),function(r) u[[r]]+1.2*psych[[r]][,1]-.6*psych[[r]][,2])
  hs <- unlist(Map(function(h,x) drop(h%*%x),H,sn))
  seednoise <- rnorm(sum(sizes)); targetnoise <- matrix(rnorm(sum(sizes)*nr,sd=.5),sum(sizes),nr)
  trueint <- lapply(seq_along(sizes),function(r) H[[r]] %*% (sweep(psych[[r]],2,colMeans(psych[[r]]),'-') * sn[[r]]))
  for (scenario in c('matched_null','matched_effect','delayed_effect','seednoise_effect')) {
    isdelay <- scenario=='delayed_effect'; effect <- scenario!='matched_null'
    seed <- hs+rep(c(1000,1200),sizes)+.3*Z
    if (scenario=='seednoise_effect') seed <- seed+.25*sd(hs)*seednoise
    d <- ppi_design(seed,psych=list(values=P,grid='neural'),TR=2,runs=runs,nuisance=matrix(Z,ncol=1,dimnames=list(NULL,'z')),hrf=h,deconv=list(lambda=.1,n_basis=24),run_effects='separate',microtime=4)
    # Oracle target mean built directly from sampled dense operators, not design columns.
    meanY <- numeric(sum(sizes));
    for (r in 1:2) {
      hh <- if(isdelay) Hd[[r]] else H[[r]]
      v <- 1.5*psych[[r]][,1]-.8*psych[[r]][,2]+.6*sn[[r]]
      if(effect) v <- v+.5*(psych[[r]][,1]-mean(psych[[r]][,1]))*sn[[r]]
      meanY[rows[[r]]] <- c(800,900)[r]+.2*Z[rows[[r]]] + drop(hh%*%v)
    }
    f0 <- ppi_fit(matrix(meanY,ncol=1),d,noise='iid')
    c0 <- ppi_contrast(f0,c(A=1,B=-1)); stopifnot(c0$estimable)
    bstar <- as.numeric(c0$estimate)
    y <- meanY+targetnoise
    fit <- ppi_fit(y,d,noise='iid'); ct <- ppi_contrast(fit,c(A=1,B=-1))
    est <- as.numeric(ct$estimate); se <- as.numeric(ct$standard_error)
    crit <- qt(.975,ct$df)
    # Residualize both constructed and true interaction columns against all
    # noninteraction columns, then compare corresponding run/condition columns.
    idx <- which(d$terms$type=='interaction'); base <- d$X[,-idx,drop=FALSE]
    Q <- qr(base); constructed <- qr.resid(Q,d$X[,idx,drop=FALSE])
    oracle <- matrix(0,nrow(d$X),length(idx))
    for (j in seq_along(idx)) {
      k <- idx[j]; r <- as.integer(d$terms$run[k]); condition <- d$terms$condition[k]
      oracle[rows[[r]],j] <- trueint[[r]][,condition]
    }
    oracle <- qr.resid(Q,oracle)
    correlations <- vapply(seq_along(idx),function(j) cor(oracle[,j],constructed[,j]),numeric(1))
    relerr <- sqrt(sum((constructed-oracle)^2)/sum(oracle^2))
    at <- at+1L
    records[[at]] <- data.frame(seed=s,scenario=scenario,bstar=bstar,mean_scale=max(abs(meanY)),rank=fit$rank,ncol=ncol(d$X),oracle_rank=qr(oracle)$rank,cor_min=min(correlations),cor_mean=mean(correlations),relerr=relerr,reject=sum(ct$p_value<.05),cover=sum(abs(est-bstar)<=crit*se),n=nr)
  }
  if(s%%10==0 || s==ns) cat('completed seed',s,'of',ns,'\n')
}
dat <- do.call(rbind,records); write.csv(dat,paste0(out,'-seeds.csv'),row.names=FALSE)
set.seed(20261024); B <- 5000L
scenarios <- unique(dat$scenario); boot <- matrix(NA_real_,B,length(scenarios),dimnames=list(NULL,scenarios))
boot_reject <- boot
for (b in seq_len(B)) {
  draw <- sample.int(ns,ns,replace=TRUE)
  for (j in seq_along(scenarios)) {
    x <- dat[dat$scenario==scenarios[j],][draw,]
    boot[b,j] <- mean(x$bstar); boot_reject[b,j] <- sum(x$reject)/sum(x$n)
  }
}
sumtab <- do.call(rbind,lapply(scenarios,function(sc) {
  x <- dat[dat$scenario==sc,]; b <- boot[,sc]; br <- boot_reject[,sc]
  data.frame(scenario=sc,seeds=ns,targets=sum(x$n),mean_bstar=mean(x$bstar),bstar_boot_low=unname(quantile(b,.025)),bstar_boot_high=unname(quantile(b,.975)),reject=sum(x$reject),rejection_rate=sum(x$reject)/sum(x$n),rejection_cluster_low=unname(quantile(br,.025)),rejection_cluster_high=unname(quantile(br,.975)),coverage=sum(x$cover)/sum(x$n),max_abs_bstar_over_scale=max(abs(x$bstar)/x$mean_scale),min_cor=min(x$cor_min),mean_cor=mean(x$cor_mean),mean_relerr=mean(x$relerr),max_relerr=max(x$relerr),min_rank=min(x$rank),min_oracle_rank=min(x$oracle_rank))
}))
write.csv(sumtab,paste0(out,'-summary.csv'),row.names=FALSE)
saveRDS(list(bstar=boot,rejection=boot_reject,seed=20261024,draws=B),paste0(out,'-bootstrap.rds'))
print(sumtab,row.names=FALSE)
