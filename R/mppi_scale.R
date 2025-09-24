# Scale-aware helpers -------------------------------------------------------

#' Internal: compute scale factors for component time courses
#' @keywords internal
.mppi_sigma_cols <- function(X) {
  if (is.null(X)) stop("Scale factors unavailable: residual matrix missing.", call. = FALSE)
  sigma <- matrixStats::colSds(X)
  sigma[is.na(sigma) | sigma < 1e-12] <- 1e-12
  sigma
}

#' Internal: apply Edin-style scaling transforms
#' @keywords internal
.mppi_scale_matrix <- function(M, sigma, mode) {
  mode <- match.arg(mode, c("raw", "amplitude", "normalized"))
  if (mode == "raw") return(M)
  if (is.null(sigma)) stop("Scale factors unavailable for rescaling interaction matrix.", call. = FALSE)
  Sinv <- diag(1 / sigma, length(sigma), length(sigma))
  if (mode == "amplitude") {
    M %*% Sinv
  } else {
    Sinv %*% M %*% Sinv
  }
}

#' Retrieve an interaction matrix in a selected scale
#'
#' @param fit `mppi_fit` object
#' @param k regressor index or name
#' @param mode scale view: `"raw"`, `"amplitude"`, or `"normalized"`
#' @return matrix in the requested scale
#' @export
mppi_get_M_scaled <- function(fit, k, mode = c("raw", "amplitude", "normalized")) {
  mode <- match.arg(mode)
  Mk_raw <- mppi_get_M(fit, k)
  if (mode == "raw") return(Mk_raw)
  sigma <- fit$sigma
  if (is.null(sigma)) {
    base_matrix <- if (!is.null(fit$R_raw)) fit$R_raw else fit$R
    sigma <- .mppi_sigma_cols(base_matrix)
  }
  .mppi_scale_matrix(Mk_raw, sigma, mode)
}
