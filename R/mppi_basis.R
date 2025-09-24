# Basis utilities ----------------------------------------------------------

#' Create an mPPI basis object
#' @param V matrix (V x r) whose columns span the desired subspace
#' @param name optional label for reporting
#' @param orthonormalize logical; orthonormalize columns via QR (default TRUE)
#' @return object of class 'mppi_basis'
as_mppi_basis <- function(V, name = "basis", orthonormalize = TRUE) {
  V <- as.matrix(V)
  if (orthonormalize) {
    qrV <- qr(V)
    Q <- qr.Q(qrV)
    rank <- qrV$rank
    V <- Q[, seq_len(rank), drop = FALSE]
  }
  structure(list(V = V, name = name, r = ncol(V)), class = "mppi_basis")
}

# Internal helper: PCA basis (columns span low-rank subspace in ROI space)
.mppi_pca_basis <- function(Y, r = NULL, center = TRUE, scale = FALSE, name = "pca") {
  stopifnot(is.matrix(Y))
  n_time <- nrow(Y)
  n_space <- ncol(Y)
  if (n_space == 0L) stop("Cannot compute a PCA basis on a matrix with zero columns.")
  r_max <- min(n_time, n_space)
  if (is.null(r)) r <- min(50L, r_max)
  r <- as.integer(r)
  if (!is.finite(r) || r <= 0L) r <- min(50L, r_max)
  r <- min(r, r_max)
  if (center || scale) {
    Yc <- scale(Y, center = center, scale = scale)
    if (!is.matrix(Yc)) Yc <- as.matrix(Yc)
  } else {
    Yc <- Y
  }
  sv <- La.svd(Yc, nu = 0L, nv = r)
  V <- t(sv$vt)
  if (ncol(V) > r) V <- V[, seq_len(r), drop = FALSE]
  basis_name <- sprintf("%s(r=%d)", name, ncol(V))
  as_mppi_basis(V, name = basis_name)
}

.mppi_resolve_basis <- function(Y, basis, r = NULL, dataset = NULL, design = NULL) {
  stopifnot(is.matrix(Y))
  coalesce <- function(...) {
    for (val in list(...)) {
      if (!is.null(val)) return(val)
    }
    NULL
  }

  # Already a concrete basis specification
  if (inherits(basis, "mppi_basis")) {
    if (is.null(basis$source)) basis$source <- "provided"
    return(basis)
  }
  if (is.matrix(basis) || (is.list(basis) && !is.null(basis$V))) {
    obj <- .mppi_check_basis(basis, ncol(Y))
    obj$source <- "provided"
    return(obj)
  }

  basis_spec <- basis
  basis_type <- "auto"
  basis_rank <- r

  if (is.list(basis_spec) && is.null(basis_spec$V)) {
    basis_type <- coalesce(basis_spec$type, basis_type)
    basis_rank <- coalesce(basis_spec$r, basis_rank)
  } else if (is.character(basis_spec)) {
    basis_type <- basis_spec[[1L]]
  } else if (!is.null(basis_spec)) {
    stop("Unsupported basis specification; provide a matrix, mppi_basis, character key, or list with type/r.")
  }

  basis_type <- tolower(basis_type)
  allowed_types <- c("auto", "pca", "roi", "none")
  basis_type <- match.arg(basis_type, allowed_types)

  if (is.null(basis_rank)) {
    basis_rank <- min(50L, min(nrow(Y), ncol(Y)))
  }
  basis_rank <- as.integer(basis_rank)
  if (!is.finite(basis_rank) || basis_rank <= 0L) {
    basis_rank <- min(50L, min(nrow(Y), ncol(Y)))
  }

  take_existing <- function(obj) {
    if (is.null(obj)) return(NULL)
    cand <- coalesce(obj$basis, obj$V, attr(obj, "basis"))
    if (is.null(cand) && is.list(obj) && !is.null(obj$basis) && is.list(obj$basis)) {
      cand <- obj$basis
    }
    cand
  }

  if (identical(basis_type, "auto")) {
    candidate <- coalesce(take_existing(dataset), take_existing(design))
    if (!is.null(candidate)) {
      obj <- .mppi_check_basis(candidate, ncol(Y))
      obj$source <- "auto"
      return(obj)
    }
    basis_type <- "pca"
  }

  if (basis_type %in% c("roi", "none")) {
    return(NULL)
  }

  if (identical(basis_type, "pca")) {
    obj <- .mppi_pca_basis(Y, r = basis_rank)
    obj$source <- "pca"
    return(obj)
  }

  stop(sprintf("Unsupported basis type '%s'.", basis_type))
}

#' Extract a basis-space interaction matrix
mppi_get_M <- function(fit, k) {
  if (is.character(k)) idx <- match(k, fit$names) else idx <- k
  if (is.na(idx) || idx < 1 || idx > length(fit$Delta)) stop("Invalid regressor index")
  D <- fit$Delta[[idx]]
  if (is.list(D) && !is.null(D$values)) {
    .mppi_unpack_upper(D)
  } else {
    D
  }
}

#' Extract a lagged interaction matrix
#' @param fit An `mppi_fit` object.
#' @param k Regressor index or name.
#' @param lag Integer lag τ (τ = 0 returns the contemporaneous matrix).
#' @return Matrix representing the interaction at lag τ in the working space.
#' @export
mppi_get_M_lag <- function(fit, k, lag = 0L) {
  lag <- as.integer(lag)
  if (lag == 0L) return(mppi_get_M(fit, k))
  if (is.null(fit$lagged)) {
    stop("Fit does not store lagged interactions; refit with lags != 0.", call. = FALSE)
  }
  if (is.character(k)) idx <- match(k, fit$names) else idx <- as.integer(k)
  if (is.na(idx) || idx < 1 || idx > length(fit$lagged)) {
    stop("Invalid regressor index for lag extraction.", call. = FALSE)
  }
  per_reg <- fit$lagged[[idx]]
  if (is.null(per_reg)) per_reg <- fit$lagged[[fit$names[idx]]]
  if (is.null(per_reg)) {
    stop("Requested regressor does not have lagged matrices stored.", call. = FALSE)
  }
  Mk <- per_reg[[as.character(lag)]]
  if (is.null(Mk)) {
    stop(sprintf("No lagged matrix stored for lag %s.", lag), call. = FALSE)
  }
  if (is.list(Mk) && !is.null(Mk$values)) {
    Mk <- .mppi_unpack_upper(Mk)
  }
  Mk
}

#' Reconstruct a submatrix of Δ from a basis fit
mppi_reconstruct_delta <- function(fit, k, rows = NULL, cols = NULL) {
  if (is.null(fit$basis)) {
    D <- mppi_get_M(fit, k)
    if (!is.null(rows) || !is.null(cols)) {
      if (is.null(rows)) rows <- seq_len(nrow(D))
      if (is.null(cols)) cols <- seq_len(ncol(D))
      D <- D[rows, cols, drop = FALSE]
    }
    return(D)
  }
  V <- fit$basis$V
  M <- mppi_get_M(fit, k)
  if (is.null(rows)) rows <- seq_len(nrow(V))
  if (is.null(cols)) cols <- seq_len(nrow(V))
  V_rows <- V[rows, , drop = FALSE]
  V_cols <- V[cols, , drop = FALSE]
  full <- V_rows %*% M %*% t(V_cols)
  if (length(rows) == length(cols) && all(rows == cols)) diag(full) <- 0
  full
}

#' Compute a single edge from a basis fit
mppi_edge <- function(fit, k, i, j) {
  if (is.null(fit$basis)) stop("mppi_edge() requires a basis-aware fit.")
  V <- fit$basis$V
  M <- mppi_get_M(fit, k)
  as.numeric(V[i, , drop = FALSE] %*% M %*% t(V[j, , drop = FALSE]))
}
