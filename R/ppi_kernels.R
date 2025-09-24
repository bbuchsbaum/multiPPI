# Internal helpers for fast kernels ---------------------------------------

#' Batch computation of interaction matrices across contexts and lags
#'
#' @param U T x r matrix of residual time-series
#' @param pk_mat T x K matrix of modulators
#' @param lags integer vector of lags (default 0)
#' @param blocklens optional integer vector of run lengths
#' @return list with elements `matrices` (list of r x r x K arrays) and
#'   `denominator` (K x length(lags) matrix of sums pk^2)
#' @keywords internal
ppi_batch <- function(U, pk_mat, lags = 0L, blocklens = NULL) {
  stopifnot(is.matrix(U), is.matrix(pk_mat), nrow(U) == nrow(pk_mat))
  res <- ppi_batch_cpp(U, pk_mat, as.integer(lags),
                       if (is.null(blocklens)) integer() else as.integer(blocklens))
  mats <- lapply(res$matrices, function(cube) {
    arr <- as.array(cube)
    dimnames(arr) <- NULL
    arr
  })
  list(matrices = mats, denominator = res$denominator)
}

#' Ridge-stabilised precision matrix
#' @keywords internal
ridge_precision <- function(S0, lambda = 0.05) {
  ridge_precision_cpp(S0, lambda)
}
