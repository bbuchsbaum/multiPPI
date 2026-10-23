args <- commandArgs(TRUE)
pkgload::load_all(args[1],quiet=TRUE)
out <- args[2]
suppressPackageStartupMessages(library(fmridesign))
set.seed(20261025)
n <- 120L; runs <- rep(1:2,each=60L)
Y <- matrix(rnorm(n*4),n,4)
sf <- sampling_frame(blocklens=c(60,60),TR=2)
events <- expand.grid(onset=c(4,24,44,64),run=1:2)
events$cond <- factor(rep(c("A","B"),4))
em <- event_model(onset ~ hrf(cond),data=events,block=~run,sampling_frame=sf)
bm <- baseline_model(basis="constant",sframe=sf)
frame <- fmrireg::matrix_frame(Y,TR=2,run_length=c(60,60),event_table=events)
model <- fmrireg::fmri_model(em,bm,frame)
cat("REAL CURRENT DATASET AND MODEL\n")
print(class(frame)); print(class(fmrireg::design_matrix(model)))
print(all.equal(unname(fmridataset::collect_assay(frame)),Y))
tryCatch(mppi_fit_from_fmrireg(model,frame,psych=1:2,domain="bold"),
         error=function(e) cat("mppi_fit_from_fmrireg:",conditionMessage(e),"\n"))
tryCatch(mppi(frame,as_mppi_design(em,bm),domain="bold"),
         error=function(e) cat("mppi(frame):",conditionMessage(e),"\n"))

cat("\nRUNWISE BASELINE RANK\n")
des <- as_mppi_design(em,bm)
cat("columns",ncol(des$X),"rank",qr(des$X)$rank,"\n")
print(colnames(des$X))
fit <- mppi(des,Y=Y)
cat("finite baseline RSS",all(is.finite(fit$evidence$rss_base)),"\n")
tryCatch(mppi_fit_whitened(Y,des$X,runs,des$psych_idx,domain="bold",p=1L),
         error=function(e) cat("direct whitening:",conditionMessage(e),"\n"))

cat("\nMAMMEN MOMENTS AND CLASSICAL CCA REFERENCE\n")
p <- (sqrt(5)+1)/(2*sqrt(5)); w <- c((1-sqrt(5))/2,(1+sqrt(5))/2)
print(c(mean=sum(w*c(1-p,p)),second_moment=sum(w^2*c(1-p,p))))
rhos <- c(.4,.2); nn <- 100; pp <- 2; qq <- 2
chi <- -((nn-1)-(pp+qq+1)/2)*sum(log(1-rhos^2))
print(c(package_p=multiPPI:::.cca_pvalue_param(rhos,nn,pp,qq),
        classical_all_roots_p=pchisq(chi,df=pp*qq,lower.tail=FALSE)))

cat("\nREFERENCE CONVENTIONAL PPI WITH FmriAR (NOT A MULTIPPI ENTRY POINT)\n")
nn <- 240L; rr <- rep(1:2,each=120L); tr <- 2
hh <- fmrihrf::hrf_spmg1(seq(0,32,by=tr)); hh <- hh/sum(hh)
causal <- function(x,h) {
  yy <- numeric(length(x))
  for(t in seq_along(x)) for(k in seq_len(min(t,length(h)))) yy[t] <- yy[t]+h[k]*x[t-k+1]
  yy
}
byrun <- function(x) {
  ans <- x
  for(ii in split(seq_len(nn),rr)) ans[ii] <- causal(x[ii],hh)
  ans
}
psych <- as.numeric(rep(rep(c(0,1),each=12),10))
neural <- as.numeric(arima.sim(list(ar=.3),n=nn))
seed_bold <- byrun(neural)
task_bold <- byrun(psych)
interaction <- byrun((neural-ave(neural,rr))* (psych-ave(psych,rr)))
D <- cbind(run1=as.numeric(rr==1),run2=as.numeric(rr==2),task=task_bold,
           physiological=seed_bold,PPI=interaction)
beta <- matrix(rnorm(ncol(D)*8),ncol(D),8)
err <- matrix(0,nn,8)
for(j in 1:8) for(ii in split(seq_len(nn),rr)) err[ii,j] <- arima.sim(list(ar=.5),n=length(ii),sd=.1)
response <- D %*% beta + err
resid <- qr.resid(qr(D),response)
plan <- fmriAR::fit_noise(resid,runs=rr,p=1L,method="ar")
xy <- fmriAR::whiten_apply(plan,X=D,Y=response,runs=rr,parallel=FALSE)
phi <- as.numeric(plan$phi[[1]])
cat("estimated common AR1",phi,"\n")
Sigma <- matrix(0,nn,nn)
for(ii in split(seq_len(nn),rr)) Sigma[ii,ii] <- outer(seq_along(ii),seq_along(ii),function(i,j) phi^abs(i-j))/(1-phi^2)
Sinv <- solve(Sigma)
expected_beta <- solve(crossprod(D,Sinv %*% D),crossprod(D,Sinv %*% response))
actual_beta <- qr.coef(qr(xy$X),xy$Y)
W <- matrix(0,nn,nn)
for(ii in split(seq_len(nn),rr)) {
  W[ii[1],ii[1]] <- sqrt(1-phi^2)
  for(j in 2:length(ii)) {W[ii[j],ii[j]] <- 1; W[ii[j],ii[j-1]] <- -phi}
}
print(c(design_whitening_error=max(abs(xy$X-W%*%D)),
        target_whitening_error=max(abs(xy$Y-W%*%response)),
        coefficient_dense_GLS_error=max(abs(actual_beta-expected_beta))))
cat("This verifies the reusable fmriAR operator, not an implemented conventional PPI API in multiPPI.\n")
capture.output(sessionInfo(),file=file.path(out,"integration-session.txt"))
