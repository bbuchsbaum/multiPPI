#' DCT-regularized hemodynamic deconvolution (fast, vectorized via Rcpp)
#'
#' Implements a DCT-basis deconvolution with Tikhonov/frequency prior as advocated by
#' Gitelman et al. (2003) for forming interaction terms in neural space. Operates column-wise
#' on a T x R matrix and chooses the penalty per column by generalized cross-validation.
#' Penalties are relative to the mean diagonal \eqn{s} of the Gram matrix of the
#' HRF-convolved DCT basis (the ridge is \eqn{\lambda s Q}), and the GCV grid is
#' scaled by \eqn{s}, so the estimate's shape does not depend on HRF amplitude.
#'
#' @param Y numeric matrix T x R (prewhitened residuals, or basis-projected components)
#' @param h numeric vector HRF kernel (in TR units)
#' @param K integer, number of DCT components (default 64 or T/2)
#' @param method "gcv" (default) for per-column GCV lambda, or "fixed"
#' @param lambda relative ridge if method="fixed" (multiplied by \eqn{s})
#' @param q optional length-K vector of frequency prior weights (diagonal of Q); default uniform
#' @return list with U (T x R neural estimates), beta (K x R coefficients), lambda
#'   (R vector of relative penalties, on the scale of the `lambda` argument)
#' @examples
#' \dontrun{
#'   library(multiPPI)
#'   TR <- 2; T <- 300; t <- seq(0, (T-1))*TR
#'   h <- fmrihrf::hrf_spmg1(seq(0, 32, by = TR))
#'   Y <- matrix(rnorm(T*10), T, 10)  # 10 components or voxels
#'   out <- mppi_deconv_dct(Y, h, K = 64, method = "gcv")
#'   U   <- out$U
#' }
#' @export
mppi_deconv_dct <- function(Y, h, K = NULL, method = c("gcv","fixed"),
                            lambda = 1e-1, q = NULL) {
  stopifnot(is.matrix(Y), is.numeric(h))
  if (is.null(K)) K <- min(64L, floor(nrow(Y)/2L))
  method <- match.arg(method)
  if (length(K) != 1L || !is.finite(K) || K < 1 || K > nrow(Y)) stop("K must be in 1..T", call. = FALSE)
  q_in <- if (is.null(q)) NULL else as.numeric(q)
  s <- .mppi_dct_gram_scale(nrow(Y), as.numeric(h), as.integer(K))
  if (!is.finite(s) || s <= 0) stop("HRF has no energy on the DCT basis.", call. = FALSE)
  out <- deconv_dct_multi(Y, h = as.numeric(h), K = as.integer(K), q_in = q_in,
                          method = method, lambda_fixed = lambda * s,
                          lam_lo = 1e-6 * s, lam_hi = 1e6 * s)
  out$lambda <- out$lambda / s
  out
}

# Mean diagonal of the Gram matrix of the causally HRF-convolved orthonormal
# DCT-II basis, as built by deconv_dct_multi().
.mppi_dct_gram_scale <- function(n, h, K) {
  B <- cos(outer(seq_len(n) - 0.5, seq_len(K) - 1) * pi / n) *
    rep(c(sqrt(1 / n), rep(sqrt(2 / n), K - 1L)), each = n)
  HB <- h[1L] * B
  for (k in seq_len(min(length(h), n) - 1L))
    HB[-seq_len(k), ] <- HB[-seq_len(k), , drop = FALSE] + h[k + 1L] * B[seq_len(n - k), , drop = FALSE]
  mean(colSums(HB^2))
}
