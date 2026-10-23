# One HRF convention for every deconvolution path: a numeric vector is a kernel
# sampled at TR from lag zero; a function (e.g. an fmrihrf HRF) is a function of
# time in seconds, evaluated on seq(0, span, by = TR) with its "span" attribute
# when present, else `span`. NULL is mppi_default_hrf().
#' @keywords internal
.mppi_hrf_to_vec <- function(hrf_obj, TR, span = 32) {
  if (is.null(hrf_obj)) return(mppi_default_hrf(tr = TR, duration = span))
  if (is.function(hrf_obj)) {
    s <- attr(hrf_obj, "span")
    if (is.numeric(s) && length(s) == 1L && is.finite(s) && s > 0) span <- s
    hrf_obj <- hrf_obj(seq(0, span, by = TR))
  }
  if (is.numeric(hrf_obj) && NCOL(hrf_obj) == 1L) return(as.numeric(hrf_obj))
  stop("HRF must be a numeric kernel sampled at TR or a function of time in seconds returning one kernel (multi-basis HRFs are unsupported).", call. = FALSE)
}

#' MAP/Tikhonov deconvolution in neural domain
#'
#' Solves \eqn{\min_u \|y - Hu\|^2 + \lambda s \|Du\|^2}{min ||y - Hu||^2 + lambda s ||Du||^2}
#' per column, with \eqn{H} causal convolution with the HRF, \eqn{D} first
#' differences and \eqn{s} the mean diagonal of \eqn{H'H}. The penalty is thus
#' relative to the operator scale: scaling the HRF by \eqn{c} scales the
#' estimate by \eqn{1/c} and leaves its shape unchanged.
#'
#' @param Y Numeric matrix (time x series)
#' @param TR Repetition time in seconds
#' @param hrf Numeric kernel sampled at `TR` from lag zero, or a function of time
#'   in seconds (e.g. an fmrihrf HRF) evaluated on `seq(0, span, by = TR)`, where
#'   `span` is the function's `span` attribute, else 32 s. `NULL` uses
#'   [mppi_default_hrf()].
#' @param lambda Relative Tikhonov weight (default 10): the first-difference
#'   penalty is `lambda` times the mean diagonal of \eqn{H'H}.
#' @return Time x series matrix of neural-domain estimates
#' @export
mppi_deconv <- function(Y, TR, hrf = NULL, lambda = 10) {
  if (!is.matrix(Y)) stop("mppi_deconv: Y must be a matrix.", call. = FALSE)
  if (!is.numeric(TR) || length(TR) != 1L || !is.finite(TR) || TR <= 0) {
    stop("mppi_deconv: TR must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("mppi_deconv: lambda must be a finite non-negative scalar.", call. = FALSE)
  }
  h <- .mppi_hrf_to_vec(hrf, TR)
  if (!length(h) || any(!is.finite(h)) || !any(h != 0)) stop("mppi_deconv: HRF must be finite and nonzero.", call. = FALSE)
  if (length(h) > nrow(Y)) h <- h[seq_len(nrow(Y))]
  # Column j of H holds h[1:min(L, n - j + 1)], so diag(H'H) are partial HRF energies.
  energy <- cumsum(h^2)
  mppi_deconv_map(Y, h = h, lambda = lambda * mean(energy[pmin(length(h), rev(seq_len(nrow(Y))))]))
}
