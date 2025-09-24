# Prewhitened estimator (fmriAR) ------------------------------------------

#' mPPI with fmriAR prewhitening
#' @param Y T x V data
#' @param X T x q design
#' @param runs vector of run labels (length T)
#' @param psych_idx indices of psychological columns
#' @param ar_method "ar" or "arma"; @param p AR order or "auto"
mppi_fit_whitened <- function(Y, X, runs, psych_idx,
                              ar_method = "ar", p = "auto",
                              zero_diag = TRUE, scale = c("normalized","cov","corr"),
                              center_by = c("none","run"), na_action = c("omit_tr","error"),
                              backend = c("blas","accumulate","chunked"), chunk_size = NULL,
                              packed = FALSE, basis = NULL,
                              domain = c("neural","bold"),
                              project_backend = c("blas","chunked"), project_chunk_cols = NULL,
                              lags = 0L, lag_blocklens = NULL) {
  if (!requireNamespace("fmriAR", quietly = TRUE)) {
    stop("fmriAR not available. Install via remotes::install_github('bbuchsbaum/fmriAR').")
  }
  stopifnot(length(runs) == nrow(Y), nrow(X) == nrow(Y))
  scale_choice <- match.arg(scale)
  domain <- match.arg(domain)
  center_by <- match.arg(center_by)
  na_action <- match.arg(na_action)
  backend <- match.arg(backend)
  project_backend <- match.arg(project_backend)
  # Initial residuals for noise fit
  qrX <- qr(X)
  res  <- Y - X %*% qr.coef(qrX, Y)
  plan <- fmriAR::fit_noise(res, runs = runs, method = ar_method, p = p)
  xyw  <- fmriAR::whiten_apply(plan, X = X, Y = Y, runs = runs)
  mppi_fit(Y = xyw$Y, X = xyw$X, psych_idx = psych_idx, runs = runs,
           zero_diag = zero_diag, scale = scale_choice,
           center_by = center_by, na_action = na_action,
           backend = backend, chunk_size = chunk_size, packed = packed,
           basis = basis, project_backend = project_backend,
           project_chunk_cols = project_chunk_cols,
           domain = domain,
           lags = lags, lag_blocklens = lag_blocklens)
}
