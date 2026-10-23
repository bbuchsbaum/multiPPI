# Engine state construction -------------------------------------------------

#' Internal helper: construct shared state for mPPI engines
#'
#' Runs the closed-form engine (`mppi_fit()`, or `mppi_fit_whitened()` when
#' `prewhiten = TRUE`) and returns its prepared quantities, so alternative
#' engines (e.g. instantaneous) share its run-local deconvolution, censoring
#' rules, prewhitening and rank-safe evidence. `index` and `segments` give
#' each retained row's acquisition position and run segment (from
#' `lag_blocklens` when given, as in `mppi_fit()`, else from `runs`), for
#' lag pairs and filtering on the acquisition grid.
#'
#' @keywords internal
.mppi_engine_state <- function(Y, X, psych_idx, runs = NULL,
                              basis = NULL,
                              domain = c("bold", "neural", "innovations"),
                              deconv = NULL,
                              center_by = c("none", "run"),
                              na_action = c("omit_tr", "error"),
                              project_backend = c("blas", "chunked"),
                              project_chunk_cols = NULL,
                              prewhiten = FALSE, ar_method = "ar", ar_order = "auto",
                              lag_blocklens = NULL) {
  stopifnot(is.matrix(Y), is.matrix(X))
  domain <- match.arg(domain)
  center_by <- match.arg(center_by)
  na_action <- match.arg(na_action)
  project_backend <- match.arg(project_backend)
  fit <- if (prewhiten) {
    mppi_fit_whitened(Y, X, runs, psych_idx, ar_method = ar_method, p = ar_order,
                      center_by = center_by, na_action = na_action, basis = basis,
                      domain = domain, deconv = deconv, project_backend = project_backend,
                      project_chunk_cols = project_chunk_cols, lag_blocklens = lag_blocklens)
  } else {
    mppi_fit(Y, X, psych_idx, runs = runs, center_by = center_by, na_action = na_action,
             basis = basis, project_backend = project_backend,
             project_chunk_cols = project_chunk_cols, domain = domain, deconv = deconv,
             lag_blocklens = lag_blocklens)
  }
  n <- nrow(Y)
  segments <- if (!is.null(lag_blocklens)) rep(seq_along(lag_blocklens), lag_blocklens) else
    if (is.null(runs)) rep(1L, n) else cumsum(c(TRUE, runs[-1] != runs[-n]))
  if (length(segments) != n) stop("lag_blocklens must describe the acquisition grid.", call. = FALSE)
  list(
    R_work = fit$R,
    R_unscaled = fit$R_raw,
    sigma = fit$sigma,
    sigma_roi = fit$sigma_roi,
    pk_list = fit$pk,
    pk_names = fit$names,
    denoms = fit$denom,
    Sigma0 = fit$Sigma0,
    evidence = fit$evidence,
    runs = fit$runs,
    index = fit$observation_index,
    segments = segments[fit$observation_index],
    lag_blocklens = fit$lag_blocklens,
    dropped = fit$dropped,
    center_by = fit$center_by,
    domain = fit$domain,
    deconv = fit$deconv,
    basis = fit$basis,
    project_backend = fit$project_backend,
    project_chunk_cols = fit$project_chunk_cols
  )
}
