#' DCT-regularized hemodynamic deconvolution (fast, vectorized via Rcpp)
#'
#' Implements a DCT-basis deconvolution with Tikhonov/frequency prior as advocated by
#' Gitelman et al. (2003) for forming interaction terms in neural space. Operates column-wise
#' on a T x R matrix and chooses the penalty per column by generalized cross-validation.
#'
#' @param Y numeric matrix T x R (prewhitened residuals, or basis-projected components)
#' @param h numeric vector HRF kernel (in TR units)
#' @param K integer, number of DCT components (default 64 or T/2)
#' @param method "gcv" (default) for per-column GCV lambda, or "fixed"
#' @param lambda fixed ridge if method="fixed"
#' @param q optional length-K vector of frequency prior weights (diagonal of Q); default uniform
#' @return list with U (T x R neural estimates), beta (K x R coefficients), lambda (R vector)
#' @examples
#' \dontrun{
#'   library(mppiDeconv)
#'   TR <- 2; T <- 300; t <- seq(0, (T-1))*TR
#'   h <- fmrihrf::spm_hrf(TR)  # or any HRF kernel vector
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
  q_in <- if (is.null(q)) NULL else as.numeric(q)
  deconv_dct_multi(Y, h = as.numeric(h), K = as.integer(K), q_in = q_in,
                   method = method, lambda_fixed = lambda)
}
