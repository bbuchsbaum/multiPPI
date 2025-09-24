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
  idx <- if (is.character(k)) match(k, fit$names) else as.integer(k)
  if (is.na(idx) || idx < 1 || idx > length(fit$Delta)) {
    stop("Invalid regressor index", call. = FALSE)
  }

  get_from_store <- function(store, attr_name) {
    if (!is.null(store) && length(store) >= idx) {
      val <- store[[idx]]
      if (!is.null(val)) return(val)
    }
    entry <- fit$Delta[[idx]]
    attr(entry, attr_name)
  }

  Mk_raw <- get_from_store(fit$Delta_raw, "raw")
  if (is.null(Mk_raw)) Mk_raw <- mppi_get_M(fit, k)
  if (is.list(Mk_raw) && !is.null(Mk_raw$values)) Mk_raw <- .mppi_unpack_upper(Mk_raw)
  if (mode == "raw") return(Mk_raw)

  if (mode == "amplitude") {
    Mk_amp <- get_from_store(fit$Delta_amplitude, "amplitude")
    if (!is.null(Mk_amp)) {
      if (is.list(Mk_amp) && !is.null(Mk_amp$values)) Mk_amp <- .mppi_unpack_upper(Mk_amp)
      return(Mk_amp)
    }
  }
  if (mode == "normalized") {
    Mk_norm <- get_from_store(fit$Delta_normalized, "normalized")
    if (!is.null(Mk_norm)) {
      if (is.list(Mk_norm) && !is.null(Mk_norm$values)) Mk_norm <- .mppi_unpack_upper(Mk_norm)
      return(Mk_norm)
    }
  }

  sigma <- fit$sigma
  if (is.null(sigma)) {
    base_matrix <- if (!is.null(fit$R_raw)) fit$R_raw else fit$R
    sigma <- .mppi_sigma_cols(base_matrix)
  }
  .mppi_scale_matrix(Mk_raw, sigma, mode)
}
