# Prewhitened estimator (fmriAR) ------------------------------------------

#' mPPI with fmriAR prewhitening
#'
#' @param Y Numeric matrix of observed time series (T x V).
#' @param X Numeric design matrix (T x q). An intercept is added before
#'   whitening when `X` spans none; the whitened design is then used as-is
#'   (no raw intercept is appended), so residuals and slopes are GLS.
#' @param runs Integer or factor vector of run labels (length T).
#' @param psych_idx Integer vector of psychological column indices within `X`.
#' @param ar_method String selecting AR noise model (`"ar"` or `"arma"`).
#' @param p AR order or `"auto"` when delegating to `fmriAR::fit_noise()`.
#' @param zero_diag Logical; zero the diagonal of each interaction matrix.
#' @param scale Output scale for interaction matrices (`"normalized"`,
#'   `"cov"`, or `"corr"`).
#' @param center_by Strategy for centering psychological regressors
#'   (`"none"` or `"run"`).
#' @param na_action Behaviour when encountering missing data (`"omit_tr"`
#'   drops time points, `"error"` aborts).
#' @param backend Computational backend for crossproducts (`"blas"`,
#'   `"accumulate"`, or `"chunked"`).
#' @param chunk_size Optional integer chunk size when using chunked backend.
#' @param packed Logical; request packed upper-triangular storage when
#'   supported.
#' @param basis Optional `mppi_basis` object for projecting into a subspace.
#' @param domain Domain for the fit (`"neural"` or `"bold"`).
#' @param deconv Reserved deconvolution specification; prewhitened neural
#'   construction is unsupported and rejected.
#' @param project_backend Backend for projecting residuals into the basis.
#' @param project_chunk_cols Chunk width for column-wise projection when
#'   `project_backend = "chunked"`.
#' @param lags Integer vector of lags to include (defaults to `0L`).
#' @param lag_blocklens Optional integer vector of block lengths for lagged
#'   fits.
mppi_fit_whitened <- function(Y, X, runs, psych_idx,
                              ar_method = "ar", p = "auto",
                              zero_diag = TRUE, scale = c("normalized","cov","corr"),
                              center_by = c("none","run"), na_action = c("omit_tr","error"),
                              backend = c("blas","accumulate","chunked"), chunk_size = NULL,
                              packed = FALSE, basis = NULL,
                              domain = c("bold","neural","innovations"), deconv = NULL,
                              project_backend = c("blas","chunked"), project_chunk_cols = NULL,
                              lags = 0L, lag_blocklens = NULL) {
  stopifnot(length(runs) == nrow(Y), nrow(X) == nrow(Y))
  scale_choice <- match.arg(scale)
  domain <- match.arg(domain)
  if (domain != "bold") stop("Prewhitening before neural deconvolution is unsupported; use domain=\"bold\" or disable prewhitening.", call. = FALSE)
  if (anyNA(Y) || anyNA(X)) stop("Prewhitening with missing observations is unsupported.", call. = FALSE)
  center_by <- match.arg(center_by)
  na_action <- match.arg(na_action)
  backend <- match.arg(backend)
  project_backend <- match.arg(project_backend)
  # GLS: add any intercept before whitening; the whitened design is then used as-is
  design <- .mppi_prepare_design_matrix(X, psych_idx)
  X <- design$X
  qrX <- qr(X)
  res  <- qr.resid(qrX, Y)
  plan <- fmriAR::fit_noise(res, runs = runs, method = ar_method, p = p)
  xyw  <- fmriAR::whiten_apply(plan, X = X, Y = Y, runs = runs)
  Xw <- xyw$X
  attr(Xw, "mppi_whitened") <- TRUE
  fit <- mppi_fit(Y = xyw$Y, X = Xw, psych_idx = design$psych_idx, runs = runs,
           zero_diag = zero_diag, scale = scale_choice,
           center_by = center_by, na_action = na_action,
           backend = backend, chunk_size = chunk_size, packed = packed,
           basis = basis, project_backend = project_backend,
           project_chunk_cols = project_chunk_cols,
           domain = domain,
           lags = lags, lag_blocklens = lag_blocklens)
  fit$domain <- "whitened_bold"
  fit$deconv$domain <- "whitened_bold"
  fit$noise_plan <- plan
  fit
}
