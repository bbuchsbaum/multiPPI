# Evidence summaries -------------------------------------------------------

#' Model evidence for an mPPI fit
#'
#' Computes information criteria improvements for the GLM used inside
#' `mppi_fit`, contrasting the full model (confounds + psychological
#' regressors) against the task-only baseline (confounds only).
#'
#' @param fit Object returned by `mppi_fit()` or wrappers.
#' @param aggregate Logical; return aggregate totals alongside per-column
#'   metrics (default `TRUE`).
#' @return A list with components `per_column` (data frame) and, when
#'   `aggregate = TRUE`, `aggregate` summarising totals. Columns include
#'   per-component AIC/BIC for the full and baseline models and their
#'   differences (baseline minus full). Positive deltas indicate improved
#'   evidence for the PPI model.
#' @export
mppi_evidence <- function(fit, aggregate = TRUE) {
  if (!inherits(fit, "mppi_fit")) {
    stop("fit must be an mppi_fit object.", call. = FALSE)
  }
  ev <- fit$evidence
  if (is.null(ev) || is.null(ev$rss_full) || is.null(ev$rss_base)) {
    stop("Evidence summaries unavailable; refit with the current multiPPI version.", call. = FALSE)
  }
  n <- ev$n
  df_full <- ev$df_full
  df_base <- ev$df_base
  rss_full <- as.numeric(ev$rss_full)
  rss_base <- as.numeric(ev$rss_base)
  if (length(rss_full) != length(rss_base)) {
    stop("Evidence vectors have mismatched lengths.", call. = FALSE)
  }
  eps <- .Machine$double.eps
  rss_full_safe <- pmax(rss_full, eps)
  rss_base_safe <- pmax(rss_base, eps)
  aic_full <- n * log(rss_full_safe / n) + 2 * df_full
  aic_base <- n * log(rss_base_safe / n) + 2 * df_base
  bic_full <- n * log(rss_full_safe / n) + df_full * log(n)
  bic_base <- n * log(rss_base_safe / n) + df_base * log(n)
  delta_aic <- aic_base - aic_full
  delta_bic <- bic_base - bic_full
  cols <- ev$column
  if (is.null(cols)) {
    cols <- paste0("component", seq_along(rss_full))
  }
  per_column <- data.frame(
    component = cols,
    rss_full = rss_full,
    rss_base = rss_base,
    AIC_full = aic_full,
    AIC_base = aic_base,
    delta_AIC = delta_aic,
    BIC_full = bic_full,
    BIC_base = bic_base,
    delta_BIC = delta_bic,
    stringsAsFactors = FALSE
  )
  if (!aggregate) {
    return(per_column)
  }
  agg <- data.frame(
    n = n,
    df_full = df_full,
    df_base = df_base,
    components = length(rss_full),
    AIC_full = sum(aic_full),
    AIC_base = sum(aic_base),
    delta_AIC = sum(delta_aic),
    BIC_full = sum(bic_full),
    BIC_base = sum(bic_base),
    delta_BIC = sum(delta_bic),
    stringsAsFactors = FALSE
  )
  list(per_column = per_column, aggregate = agg)
}
