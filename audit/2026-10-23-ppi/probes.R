args <- commandArgs(TRUE)
pkg <- if (length(args)) args[1] else "/Users/bbuchsbaum/code/multiPPI"
out <- if (length(args) > 1) args[2] else tempdir()
pkgload::load_all(pkg, quiet = TRUE)
ns <- asNamespace("multiPPI")
fn <- function(name) get(name, ns)
options(width = 110)
set.seed(20261023)
selected <- Sys.getenv("AUDIT_PROBE", "")
results <- if (nzchar(selected) && file.exists(file.path(out,"probe-results.rds"))) readRDS(file.path(out,"probe-results.rds")) else list()
probe <- function(name, code) {
  if (nzchar(selected) && !grepl(selected,name)) return(invisible(NULL))
  cat("\nPROBE:", name, "\n")
  ans <- tryCatch(force(code), error = function(e) list(error = conditionMessage(e)))
  print(ans)
  results[[name]] <<- ans
  invisible(ans)
}
convmat <- function(h, n) {
  H <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(i))
    if (i-j+1 <= length(h)) H[i,j] <- h[i-j+1]
  H
}
n <- 120L
p <- rep(c(-1, 1), each = n/2)
X <- cbind(Intercept = 1, task = p)
Y <- matrix(rnorm(n*4), n, 4)
runs <- rep(1:2, each = n/2)

probe("conventional_PPI_counterexample", {
  z <- seq(-1, 1, length.out = n/2)
  z <- z / sqrt(mean(z^2))
  s <- c(z, 2*z)
  e <- qr.resid(qr(cbind(1, p, s, s*p)), rnorm(n)) * .05
  y <- .7*s + e
  conventional <- coef(lm(y ~ s*p))[["s:p"]]
  fit <- mppi_fit(cbind(seed = s, target = y), X, 2L)
  list(conventional_interaction = conventional,
       covariance_slope = fit$Delta_raw[[1]][1,2],
       normalized_covariance_slope = fit$Delta_normalized[[1]][1,2],
       condition_seed_variances = c(mean(s[p==-1]^2), mean(s[p==1]^2)))
})
fit <- mppi_fit(Y, X, 2L, runs = runs)
probe("default_permutation_is_identity", {
  omni <- mppi_omnibus(fit, B = 49L, seed = 7)[[1]]
  edge <- mppi_permute(fit, B = 49L, seed = 7)
  list(global_p = omni$p_global, unique_null_statistics = length(unique(omni$Qnull)),
       observed_Q = omni$Q, null_Q = omni$Qnull[1],
       range_edge_p = range(edge, na.rm=TRUE))
})
probe("omnibus_observed_null_diagonal_mismatch", {
  omni <- mppi_omnibus(fit, B=9L, seed=7)[[1]]
  M <- crossprod(fit$R, fit$pk[[1]]*fit$R)/sum(fit$pk[[1]]^2)
  list(null_minus_observed = omni$Qnull[1]-omni$Q,
       diagonal_energy = sum(diag(M)^2))
})
probe("CCA_permutations_preserve_test_statistic", {
  A <- scale(matrix(rnorm(n*2),n,2),scale=FALSE)
  B <- scale(A %*% matrix(c(1,.2,.1,1),2) + matrix(rnorm(n*2,sd=.1),n,2),scale=FALSE)
  cca <- fn(".partial_cca")
  r0 <- cca(A,B)$r1
  differences <- replicate(30, {
    sg <- sample(c(-1,1),n,TRUE)
    cca(A*sg,B*sg)$r1-r0
  })
  list(observed_r1=r0, max_null_deviation=max(abs(differences)),
       reported_p=fn(".cca_pvalue_perm")(A,B,perms=49L),
       unregularized_oracle_max_error=max(abs(cca(A,B,0,0)$rhos-cancor(A,B)$cor)))
})
probe("public_prewhitening_route", {
  mppi(Y, X=X, psych_idx=2L, runs=runs, prewhiten=TRUE, domain="bold")
  "completed"
})
probe("direct_prewhitening_route", {
  fw <- mppi_fit_whitened(Y,X,runs,2L,domain="bold",p=1L)
  list(dim=dim(fw$R), finite=all(is.finite(fw$R)))
})
probe("seed_neural_dependency", {
  seed_ppi_mcca(Y,seed=1L,design=matrix(p,n,1),neighborhoods=list(2:4),
                tr=2,domain="neural",prewhiten=FALSE)
  "completed"
})
probe("fmrireg_bridge", {
  fmodel <- structure(list(),class="audit_model")
  design_matrix.audit_model <- function(x,...) X
  registerS3method("design_matrix","audit_model",design_matrix.audit_model,
                   envir=asNamespace("fmridesign"))
  mppi_fit_from_fmrireg(fmodel,list(),psych=2L,domain="bold")
  "completed"
})
probe("actual_dependency_exports", {
  desc <- read.dcf(file.path(pkg,"DESCRIPTION"))
  declared <- trimws(unlist(strsplit(paste(desc[,c("Imports","Suggests")],collapse=","),",")))
  list(fmrihrf_deconvolve_series="deconvolve_series" %in% getNamespaceExports("fmrihrf"),
       fmrireg_get_data_matrix="get_data_matrix" %in% getNamespaceExports("fmrireg"),
       fmridataset_get_data_matrix="get_data_matrix" %in% getNamespaceExports("fmridataset"),
       undeclared_packages=setdiff(c("fmridesign","fmridataset","testthat"),declared))
})
probe("real_event_model_TR_propagation", {
  hrf <- fmridesign::hrf
  em <- fmridesign::event_model(onset ~ hrf(cond),
    data=data.frame(onset=c(4,24,44,64),cond=factor(c("A","B","A","B")),run=1),
    block=~run,sampling_frame=fmridesign::sampling_frame(blocklens=n,TR=2))
  ff <- mppi(em,Y=Y,domain="neural",deconv=list(K=24,method="fixed"))
  list(design_TR=mppi_tr(em),fit_TR=ff$deconv$tr,
       design_dim=dim(fmridesign::design_matrix(em)),
       task_and_psych_identical=isTRUE(all.equal(mppi_psych(as_mppi_design(em),T=n),
                                              mppi_task(as_mppi_design(em)))))
})
probe("MAP_preserves_HRF_time_origin", {
  h <- c(0,1)
  u <- sin(seq_len(30)*.6)
  yy <- matrix(convmat(h,30) %*% u,30,1)
  got <- mppi_deconv(yy,TR=1,hrf=h,lambda=.1)
  identity_h <- mppi_deconv(yy,TR=1,hrf=1,lambda=.1)
  H <- convmat(h,30)
  D <- diff(diag(30))
  expected <- solve(crossprod(H)+.1*crossprod(D),crossprod(H,yy))
  list(error_vs_causal_oracle=max(abs(got-expected)),
       difference_from_wrong_zero_delay_HRF=max(abs(got-identity_h)))
})
probe("DCT_matches_independent_dense_ridge", {
  nn <- 40L; K <- 12L; h <- c(0,.3,.8,.5,.1); lam <- .17
  yy <- matrix(rnorm(nn*3),nn,3); q <- seq(.4,2,length.out=K)
  B <- outer(seq_len(nn)-.5,0:(K-1),function(t,k) cos(pi*t*k/nn))*sqrt(2/nn)
  B[,1] <- 1/sqrt(nn)
  HB <- convmat(h,nn) %*% B
  expected <- B %*% solve(crossprod(HB)+lam*diag(q),crossprod(HB,yy))
  got <- fn("mppi_deconv_dct")(yy,h,K=K,method="fixed",lambda=lam,q=q)$U
  list(max_absolute_error=max(abs(got-expected)))
})
probe("neural_run_isolation", {
  pp <- rnorm(n)
  XX <- cbind(Intercept=1,run2=as.numeric(runs==2),task=pp)
  YY <- Y
  addition <- qr.resid(qr(cbind(1,pp[runs==2])),matrix(rnorm((n/2)*4),n/2,4))*8
  YY[runs==2,] <- YY[runs==2,]+addition
  b1 <- mppi_fit(Y,XX,3L,runs=runs)
  b2 <- mppi_fit(YY,XX,3L,runs=runs)
  opts <- list(K=n,method="fixed",lambda=.1,hrf=c(0,.5,1,.5))
  n1 <- mppi_fit(Y,XX,3L,runs=runs,domain="neural",deconv=opts)
  n2 <- mppi_fit(YY,XX,3L,runs=runs,domain="neural",deconv=opts)
  list(bold_run1_max_change=max(abs(b1$R[runs==1,]-b2$R[runs==1,])),
       neural_run1_max_change=max(abs(n1$R[runs==1,]-n2$R[runs==1,])))
})
probe("lagged_packed_preserves_direction", {
  full <- mppi_fit(Y,X,2L,lags=1L)
  pack <- mppi_fit(Y,X,2L,lags=1L,packed=TRUE)
  a <- mppi_get_M_lag(full,1L,1L)
  b <- mppi_get_M_lag(pack,1L,1L)
  list(full_asymmetry=max(abs(a-t(a))),packed_asymmetry=max(abs(b-t(b))),
       max_roundtrip_error=max(abs(a-b)))
})
probe("censor_gap_creates_false_lag_neighbor", {
  YY <- Y; YY[60,] <- NA_real_
  f <- mppi_fit(YY,X,2L,lags=1L,zero_diag=FALSE)
  original <- setdiff(seq_len(n),60L)
  ii <- which(diff(original)==1L)
  expected <- crossprod(f$R[ii,,drop=FALSE],f$pk[[1]][ii]*f$R[ii+1,,drop=FALSE])/sum(f$pk[[1]][ii]^2)
  list(max_error=max(abs(attr(f$lagged[[1]][["1"]],"raw")-expected)),
       actual_pairs=nrow(f$R)-1L,valid_pairs=length(ii))
})
probe("rank_selection_retains_negative_modes", {
  nn <- 800L; pp <- rep(c(-1,1),nn/2)
  RR <- (sqrt(2-pp)*rnorm(nn)) %o% c(1,2,3)
  delta <- crossprod(RR,pp*RR)/sum(pp^2)
  set.seed(98)
  cv <- mppi_rank_cv(RR,pp,delta,R_reps=5)
  list(eigenvalues=eigen(delta,symmetric=TRUE)$values,selected_rank=cv$r,losses=cv$losses)
})
probe("evidence_depends_on_task_activation_not_PPI", {
  YY <- matrix(rnorm(n*4,sd=.1),n,4) + 5*p
  ff <- mppi_fit(YY,X,2L)
  f1 <- mppi_fit(YY,X,2L,domain="neural",deconv=list(K=20,hrf=c(0,1,.2),method="fixed"))
  f2 <- mppi_fit(YY,X,2L,domain="neural",deconv=list(K=20,hrf=c(0,.2,1),method="fixed"))
  list(task_delta_AIC=mppi_evidence(ff)$aggregate$delta_AIC,
       HRF_AIC_difference=max(abs(f1$AIC_total-f2$AIC_total)),
       HRF_slope_difference=max(abs(f1$Delta_raw[[1]]-f2$Delta_raw[[1]])))
})
probe("seed_matrix_and_index_use_same_preprocessing", {
  DD <- cbind(task=p)
  YY <- Y; YY[,1] <- YY[,1]+3*p
  args0 <- list(data=YY,design=DD,neighborhoods=list(2:4),domain="hemo",prewhiten=FALSE,
                pca=list(seed_k=1,light_k=2,center=TRUE,scale=TRUE))
  a <- do.call(seed_ppi_mcca,c(args0,list(seed=1L)))
  b <- do.call(seed_ppi_mcca,c(args0,list(seed=YY[,1,drop=FALSE])))
  list(index_stat=a$stat,matrix_stat=b$stat,difference=a$stat-b$stat)
})
saveRDS(results,file.path(out,"probe-results.rds"))
capture.output(sessionInfo(),file=file.path(out,"probe-session.txt"))
cat("\nCompleted",length(results),"audit probes.\n")
