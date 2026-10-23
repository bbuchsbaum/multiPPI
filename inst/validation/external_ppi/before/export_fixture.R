args <- commandArgs(TRUE)
if (length(args) != 2L) stop('usage: Rscript export_fixture.R SOURCE_COPY OUTPUT_DIR')
pkgload::load_all(args[1], quiet=TRUE)
out <- args[2]; dir.create(out, recursive=TRUE, showWarnings=FALSE)
set.seed(8252026); n <- 100L; runs <- rep(c('r1','r2'),each=50L)
psych <- cbind(A=rnorm(n),B=rnorm(n)); seed <- rnorm(n)
d <- ppi_design(seed,psych=psych,runs=runs,hrf=c(0,1,.3),microtime=1L,deconv=list(lambda=.3,n_basis=25L),run_effects='separate')
Y <- d$X %*% matrix(rnorm(ncol(d$X)*3),ncol(d$X),3) + matrix(rnorm(n*3),n,3)
colnames(Y) <- paste0('target',1:3)
plan <- fmriAR::fit_noise(qr.resid(qr(d$X),Y),runs=d$runs,p=1L,method='ar',parallel=FALSE)
plan$phi[[1]] <- .42
write_case <- function(name, des, noise) {
  fit <- ppi_fit(Y,des,noise=noise)
  ct <- ppi_contrast(fit,c(A=1,B=-1))
  cf <- ppi_contrast(fit,rbind(c(A=1,B=0),c(A=0,B=1)),test='F')
  item <- list(X=unname(des$X),Y=unname(Y),keep=unname(des$keep),runs=unname(des$runs),
    design_names=colnames(des$X),rank=fit$rank,df=fit$df,coefficients=unname(fit$coefficients),
    sigma2=unname(fit$sigma2),standard_errors=unname(fit$standard_errors),
    t_contrast=unname(ct$contrast),t_estimate=unname(ct$estimate),t_se=unname(ct$standard_error),t_stat=unname(ct$statistic),
    f_contrast=unname(cf$contrast),f_stat=unname(cf$statistic),f_df=cf$numerator_df,
    phi=if(inherits(noise,'fmriAR_plan')) as.numeric(noise$phi[[1]]) else NULL)
  jsonlite::write_json(item,file.path(out,paste0(name,'.json')),digits=17,auto_unbox=TRUE,na='null')
}
write_case('iid',d,'iid')
censored <- d; censored$keep[c(5,18,57,58)] <- FALSE; write_case('censored_iid',censored,'iid')
aliased <- d; aliased$X <- cbind(aliased$X, duplicate=aliased$X[,1]); aliased$terms <- rbind(aliased$terms,data.frame(name='duplicate',type='nuisance',condition=NA,run=NA)); write_case('rank_deficient',aliased,'iid')
write_case('fixed_ar',d,plan)
cat('fixture exported', as.character(packageVersion('multiPPI')), 'fmriAR', as.character(packageVersion('fmriAR')), '\n')
