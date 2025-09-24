#' @keywords internal
.mppi_default_spm <- function(TR, L = 32) {
  t <- seq(0, L, by = TR)
  a1 <- 6; a2 <- 16; b1 <- 1; b2 <- 1; c <- 1 / 6
  h <- (t^(a1 - 1) * exp(-t / b1) / (b1^a1 * gamma(a1))) -
       c * (t^(a2 - 1) * exp(-t / b2) / (b2^a2 * gamma(a2)))
  as.numeric(h)
}

#' @keywords internal
.mppi_hrf_to_vec <- function(hrf_obj, TR) {
  if (is.null(hrf_obj)) return(.mppi_default_spm(TR))
  if (is.numeric(hrf_obj)) return(as.numeric(hrf_obj))
  if (is.function(hrf_obj)) {
    h <- try(hrf_obj(TR = TR), silent = TRUE)
    if (!inherits(h, "try-error") && is.numeric(h)) return(as.numeric(h))
  }
  if (requireNamespace("fmrihrf", quietly = TRUE)) {
    out <- try(fmrihrf::spm_hrf(TR), silent = TRUE)
    if (!inherits(out, "try-error")) return(as.numeric(out))
  }
  .mppi_default_spm(TR)
}

#' MAP/Tikhonov deconvolution in neural domain
#'
#' @param Y Numeric matrix (time x series)
#' @param TR Repetition time in seconds
#' @param hrf Optional haemodynamic response specification (numeric, function, or NULL)
#' @param lambda Tikhonov weight (default 10)
#' @return Time x series matrix of neural-domain estimates
#' @export
mppi_deconv <- function(Y, TR, hrf = NULL, lambda = 10) {
  if (!is.matrix(Y)) stop("mppi_deconv: Y must be a matrix.", call. = FALSE)
  if (!is.numeric(TR) || length(TR) != 1L || !is.finite(TR) || TR <= 0) {
    stop("mppi_deconv: TR must be a positive scalar.", call. = FALSE)
  }
  h <- .mppi_hrf_to_vec(hrf, TR)
  h <- h[h != 0]
  if (!length(h)) stop("mppi_deconv: HRF vector is empty after removing zeros.", call. = FALSE)
  if (length(h) > nrow(Y)) h <- h[seq_len(nrow(Y))]
  mppi_deconv_map(Y, h = h, lambda = lambda)
}
