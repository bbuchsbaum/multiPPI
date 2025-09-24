# Low-rank selection / shrinkage ------------------------------------------

#' Cross-validated rank selection for ΔΣ (component space recommended)
#'
#' This helper recomputes task-modulated slopes on training/validation splits and
#' scores low-rank approximations via held-out Frobenius loss. It is agnostic to
#' whether you operate in voxel/ROI or basis space, but for high-dimensional
#' data you should supply the reduced residual matrix (e.g., `fit$Z`) together
#' with its corresponding slope matrix (`mppi_get_M(fit, k)`) so the SVD is
#' performed on an r × r matrix rather than a 50k × 50k object.
#'
#' @param R Residualized time series in the working space (T × V or T × r).
#' @param pk Residualized psychological regressor (length T).
#' @param Delta Optional observed slope matrix (V × V); retained for backward
#'   compatibility but recomputed internally from `R`/`pk` in each split.
#' @param holdout_frac Fraction of time points for hold-out blocks.
#' @param R_reps Number of random splits.
#' @return List with the selected rank `r`, per-rank losses, and the rank grid.
mppi_rank_cv <- function(R, pk, Delta, holdout_frac = 0.2, R_reps = 3) {
  Tn <- nrow(R); V <- ncol(R)
  rmax <- min(50L, V)
  ranks <- 0:rmax
  losses <- matrix(0, nrow = R_reps, ncol = length(ranks))
  for (rep in seq_len(R_reps)) {
    H <- sort(sample.int(Tn, floor(holdout_frac*Tn)))
    K <- setdiff(seq_len(Tn), H)
    denomK <- sum(pk[K]^2); denomH <- sum(pk[H]^2)
    DeltaK <- if (denomK < .Machine$double.eps) {
      matrix(0, V, V)
    } else {
      crossprod(R[K, , drop = FALSE], pk[K] * R[K, , drop = FALSE]) / denomK
    }
    targetH <- if (denomH < .Machine$double.eps) {
      matrix(0, V, V)
    } else {
      crossprod(R[H, , drop = FALSE], pk[H] * R[H, , drop = FALSE]) / denomH
    }
    sK <- svd(DeltaK, nu = min(rmax, V), nv = min(rmax, V))
    for (j in seq_along(ranks)) {
      r <- ranks[j]
      if (r == 0) {
        approx <- matrix(0, V, V)
      } else {
        approx <- sK$u[, 1:r, drop = FALSE] %*%
          diag(sK$d[1:r], nrow = r, ncol = r) %*%
          t(sK$u[, 1:r, drop = FALSE])
      }
      losses[rep, j] <- mean((targetH - approx)^2, na.rm = TRUE)
    }
  }
  rstar <- ranks[ which.min(colMeans(losses)) ]
  list(r = rstar, losses = as.numeric(colMeans(losses)), ranks = ranks)
}
