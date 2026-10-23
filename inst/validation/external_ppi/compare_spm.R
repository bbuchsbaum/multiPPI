# Compare actual SPM microtime=1 output with the frozen multiPPI construction.
args <- commandArgs(TRUE)
if (length(args) != 3L) stop('usage: Rscript compare_spm.R SOURCE_COPY FIXTURE_DIR OUTPUT_JSON')
pkgload::load_all(args[1],quiet=TRUE)
source_sha <- function(f) strsplit(system2("shasum",c("-a","256",shQuote(file.path(args[1],"R",f))),stdout=TRUE)," +")[[1]][1]
source_hashes <- list(ppi_design=source_sha("ppi_design.R"),ppi_fit=source_sha("ppi_fit.R"))
a <- as.matrix(read.csv(file.path(args[2],'actual_ppi.csv'),header=FALSE)); h <- scan(file.path(args[2],'actual_hrf.csv'),quiet=TRUE)
seed <- a[,1]; psych <- matrix(a[,2],ncol=1,dimnames=list(NULL,'A'))
build <- function(seed_adjust=NULL) {
  args <- list(seed=seed,psych=psych,TR=1,microtime=1,sampling_offset=0,hrf=h,
               center=TRUE,deconv=list(n_basis=40L,lambda=.1))
  if (!is.null(seed_adjust)) args$seed_adjust <- seed_adjust
  do.call(ppi_design,args)
}
measure <- function(x,y) list(max_abs=max(abs(x-y)),rmse=sqrt(mean((x-y)^2)),cor=unname(cor(x,y)))
compare <- function(d) {
  get <- function(type) d$X[,which(d$terms$type==type)]
  list(seed_adjust=d$provenance$seed_adjust,task=measure(get('task'),a[,4]),
       physiology_raw=measure(get('physiology'),a[,3]),
       physiology_centered=measure(get('physiology')-mean(get('physiology')),a[,3]),
       neural_seed=measure(d$neural_seed,a[,5]),interaction_raw=measure(get('interaction'),a[,6]),
       interaction_centered=measure(get('interaction')-mean(get('interaction')),a[,6]))
}
raw <- build('none'); default <- build()
result <- list(source_hashes=source_hashes, spm_commit='54dc957191112dcad116e3eaa457e881fa8260de', n=nrow(a),
 conventions=list(TR=1,microtime=1,sampling_offset=0,shared_kernel='actual SPM spm_hrf(1)',
  psych='SPM detrended psychological input',spm='PEB deconvolution; padded DCT, seed and PPI detrended',
  multipppi='DCT ridge lambda=0.1, 40 basis, unpadded; neural seed not centered'),
 raw_seed=compare(raw),default_seed_adjustment=compare(default),
 qualification='Descriptive output differences only; different seed adjustment, deconvolution priors, padding, and centering preclude a numerical equivalence claim.')

jsonlite::write_json(result,args[3],digits=17,auto_unbox=TRUE,pretty=TRUE)
cat(jsonlite::toJSON(result,digits=8,auto_unbox=TRUE,pretty=TRUE),'\n')
