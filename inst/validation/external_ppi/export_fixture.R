args <- commandArgs(TRUE)
if (length(args) != 2L) stop('usage: Rscript export_fixture.R SOURCE_COPY OUTPUT_DIR')
pkgload::load_all(args[1], quiet=TRUE)
source_sha <- function(f) strsplit(system2("shasum",c("-a","256",shQuote(file.path(args[1],"R",f))),stdout=TRUE)," +")[[1]][1]
source_hashes <- list(ppi_design=source_sha("ppi_design.R"),ppi_fit=source_sha("ppi_fit.R"))
out <- args[2]; dir.create(out, recursive=TRUE, showWarnings=FALSE)
set.seed(8252026); n <- 100L; runs <- rep(c('r1','r2'),each=50L)
psych <- cbind(A=rnorm(n),B=rnorm(n)); seed <- rnorm(n)
d <- ppi_design(seed,psych=psych,runs=runs,hrf=c(0,1,.3),microtime=1L,deconv=list(lambda=.3,n_basis=25L),run_effects='separate')
Y <- d$X %*% matrix(rnorm(ncol(d$X)*3),ncol(d$X),3) + matrix(rnorm(n*3),n,3)
colnames(Y) <- paste0('target',1:3)
plan <- fmriAR::fit_noise(qr.resid(qr(d$X),Y),runs=d$runs,p=1L,method='ar',parallel=FALSE)
plan$phi[[1]] <- .42
write_case <- function(name, des, noise, yy=Y, noise_control=list()) {
  fit <- ppi_fit(yy,des,noise=noise,noise_control=noise_control)
  ct <- ppi_contrast(fit,c(A=1,B=-1))
  cf <- ppi_contrast(fit,rbind(c(A=1,B=0),c(A=0,B=1)),test='F')
  ct_small <- ppi_contrast(fit,c(A=1e-12,B=-1e-12))
  cf_small <- ppi_contrast(fit,1e-12*rbind(c(A=1,B=0),c(A=0,B=1)),test='F')
  item <- list(source_hashes=source_hashes, X=unname(des$X),Y=unname(yy),keep=unname(des$keep),runs=unname(des$runs),
    design_names=colnames(des$X),rank=fit$rank,df=fit$df,coefficients=unname(fit$coefficients),
    sigma2=unname(fit$sigma2),standard_errors=unname(fit$standard_errors),
    t_contrast=unname(ct$contrast),t_estimate=unname(ct$estimate),t_se=unname(ct$standard_error),t_stat=unname(ct$statistic),
    f_contrast=unname(cf$contrast),f_stat=unname(cf$statistic),f_df=cf$numerator_df,
    small_t_stat=unname(ct_small$statistic),small_f_stat=unname(cf_small$statistic),
    phi=if(inherits(noise,'fmriAR_plan')) as.numeric(noise$phi[[1]]) else NULL,
    groups=if (!is.null(fit$groups)) fit$groups else NULL,
    group_phi=if (!is.null(fit$groups)) vapply(fit$groups, function(g) as.numeric(fit$noise_plan[[g]]$phi[[1]]), numeric(1)) else NULL)
  jsonlite::write_json(item,file.path(out,paste0(name,'.json')),digits=17,auto_unbox=TRUE,na='null')
}
write_case('iid',d,'iid')
censored <- d; censored$keep[c(5,18,57,58)] <- FALSE; write_case('censored_iid',censored,'iid')
aliased <- d; aliased$X <- cbind(aliased$X, duplicate=aliased$X[,1]); aliased$terms <- rbind(aliased$terms,data.frame(name='duplicate',type='nuisance',condition=NA,run=NA)); write_case('rank_deficient',aliased,'iid')
write_case('fixed_ar',d,plan)
scaled <- d; scaled$X[,1] <- scaled$X[,1]*1e-8; scaled$X[,2] <- scaled$X[,2]*1e8
write_case('scaled_columns',scaled,'iid')
# Heterogeneous AR series: grouping is estimated by multiPPI, then frozen for
# the external conditional-operator comparison; this is not calibration.
groups <- c('positive','negative','positive'); phi <- c(.65,-.45,.65)
err <- matrix(0,n,3)
for (j in seq_len(3)) for (r in unique(runs)) {
  ii <- which(runs==r); z <- rnorm(length(ii)); err[ii[1],j] <- z[1]
  for (k in 2:length(ii)) err[ii[k],j] <- phi[j]*err[ii[k-1],j]+z[k]
}
Yg <- d$X %*% matrix(rnorm(ncol(d$X)*3),ncol(d$X),3) + err
colnames(Yg) <- colnames(Y)
write_case('grouped_ar',d,'fmriAR',yy=Yg,noise_control=list(pooling='parcel',parcels=groups,p=1L))
cat('fixture exported', as.character(packageVersion('multiPPI')), 'fmriAR', as.character(packageVersion('fmriAR')), '\n')
