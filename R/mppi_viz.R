# Visualization helpers ----------------------------------------------------

#' Leading eigenmodes of interaction structure
#'
#' For matrix inputs the symmetric eigenmodes are returned. When supplied with an
#' `mppi_fit`, the communication subspace across conditions is estimated using the
#' requested scale view and optional ridge penalty.
#'
#' @param Delta interaction matrix or `mppi_fit`
#' @param r number of modes to retain
#' @param view scale view when `Delta` is an `mppi_fit`
#' @param penalty aggregate penalty for the communication subspace
#' @param lambda ridge strength when `penalty = "ridge"`
#' @return list containing eigenvalues/vectors and, for fits, per-condition projections
mppi_modes <- function(Delta, r = 3L, view = c("normalized", "raw", "amplitude"),
                       penalty = c("none", "ridge"), lambda = 1e-4) {
  if (inherits(Delta, "mppi_fit")) {
    view <- match.arg(view)
    penalty <- match.arg(penalty)
    mats <- lapply(seq_along(Delta$names), function(idx) mppi_get_M_scaled(Delta, idx, mode = view))
    if (!length(mats)) stop("Fit does not contain interaction matrices.", call. = FALSE)
    A <- Reduce(`+`, lapply(mats, function(M) M %*% t(M)))
    A <- 0.5 * (A + t(A))
    if (penalty == "ridge") {
      A <- A + lambda * diag(nrow(A))
    }
    ev <- eigen(A, symmetric = TRUE)
    keep <- seq_len(min(r, length(ev$values)))
    W <- ev$vectors[, keep, drop = FALSE]
    base_names <- rownames(mats[[1]])
    if (!is.null(base_names)) rownames(W) <- base_names
    colnames(W) <- paste0("mode", seq_along(keep))
    k <- length(keep)
    lab <- colnames(W)
    projected <- lapply(mats, function(M) {
      P <- crossprod(W, M %*% W)
      if (!is.matrix(P) || nrow(P) != k || ncol(P) != k) {
        P <- matrix(P, nrow = k, ncol = k)
      }
      if (!is.null(lab) && length(lab) == k) {
        dimnames(P) <- list(lab, lab)
      }
      P
    })
    names(projected) <- Delta$names
    list(values = ev$values[keep], vectors = W,
         projected = projected, view = view, penalty = penalty)
  } else {
    r <- max(1L, r)
    ev <- eigen((Delta + t(Delta))/2, symmetric = TRUE)
    keep <- seq_len(min(r, length(ev$values)))
    list(values = ev$values[keep], vectors = ev$vectors[, keep, drop = FALSE])
  }
}

#' Nodewise degree of modulation (L1 or L2)
mppi_degree <- function(Delta, type = c("L1","L2")) {
  type <- match.arg(type)
  A <- abs(Delta)
  if (type == "L1") rowSums(A, na.rm = TRUE) else sqrt(rowSums(A^2, na.rm = TRUE))
}

#' Project ΔΣ onto a hypothesis mask W (VxV)
mppi_project <- function(Delta, W) {
  stopifnot(all(dim(Delta) == dim(W)))
  sum(Delta * W, na.rm = TRUE)
}
