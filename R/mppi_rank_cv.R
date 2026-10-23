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
#' @param runs Optional contiguous run labels; each run contributes one held-out block.
#' @param gap Number of scans excluded from training on each side of a held-out block.
#' @details R and pk are externally preprocessed; this is not nested preprocessing CV.
#' @return List with the selected rank `r`, per-rank losses, and the rank grid.
mppi_rank_cv <- function(R, pk, Delta = NULL, holdout_frac = 0.2, R_reps = 3, runs = NULL, gap = 0L) {
  if (!is.matrix(R) || !is.numeric(R) || !nrow(R) || !ncol(R) || any(!is.finite(R)) || !is.numeric(pk) || length(pk) != nrow(R) || any(!is.finite(pk))) stop("R and pk must be finite numeric data with matching rows")
  if (length(holdout_frac) != 1L || !is.finite(holdout_frac) || holdout_frac <= 0 || holdout_frac >= 1) stop("holdout_frac must be between zero and one")
  if (length(R_reps) != 1L || !is.finite(R_reps) || R_reps < 1 || R_reps != floor(R_reps)) stop("R_reps must be a positive integer")
  if (length(gap) != 1L || !is.finite(gap) || gap < 0 || gap != floor(gap)) stop("gap must be a nonnegative integer")
  Tn <- nrow(R); V <- ncol(R)
  if (is.null(runs)) runs <- rep(1L, Tn)
  if (length(runs) != Tn || anyNA(runs) || anyDuplicated(rle(as.character(runs))$values)) stop("runs must be contiguous and identify every row")
  blocks <- split(seq_len(Tn), factor(runs, levels = unique(runs)))
  if (any(lengths(blocks) < 2L)) stop("each run must contain at least two scans")
  splits <- vector("list", R_reps)
  rmax <- min(50L, V)
  ranks <- 0:rmax
  losses <- matrix(0, nrow = R_reps, ncol = length(ranks))
  for (rep in seq_len(R_reps)) {
    H <- K <- integer()
    for (block in blocks) {
      h <- min(length(block) - 1L, max(1L, floor(holdout_frac * length(block))))
      start <- sample.int(length(block) - h + 1L, 1L)
      hold <- seq.int(start, length.out = h)
      excluded <- seq.int(max(1L, start - gap), min(length(block), start + h - 1L + gap))
      H <- c(H, block[hold]); K <- c(K, block[-excluded])
    }
    if (!length(K)) stop("gap and holdout leave no training observations")
    splits[[rep]] <- list(train = K, test = H)
    if (sum(pk[K]^2) < .Machine$double.eps || sum(pk[H]^2) < .Machine$double.eps) stop("psychological regressor has no energy in a split")
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
    effective_rank <- sum(sK$d > max(sK$d) * V * .Machine$double.eps)
    for (j in seq_along(ranks)) {
      r <- min(ranks[j], effective_rank)
      if (r == 0) {
        approx <- matrix(0, V, V)
      } else {
        approx <- sK$u[, 1:r, drop = FALSE] %*%
          diag(sK$d[1:r], nrow = r, ncol = r) %*%
          t(sK$v[, 1:r, drop = FALSE])
      }
      losses[rep, j] <- mean((targetH - approx)^2, na.rm = TRUE)
    }
  }
  rstar <- ranks[ which.min(colMeans(losses)) ]
  list(r = rstar, losses = as.numeric(colMeans(losses)), ranks = ranks, splits = splits)
}
