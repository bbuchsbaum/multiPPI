# Partial connectivity / interpretation -----------------------------------

#' First-order precision update from covariance slope
#' @param Sigma0 baseline covariance (R'R divided by T) or correlation
#' @param DeltaSigma slope matrix for regressor k
#' @param lambda ridge parameter for invertibility
#' @return matrix ΔΘ
#' @export
mppi_to_partial <- function(Sigma0, DeltaSigma, lambda = 1e-3) {
  V <- nrow(Sigma0); I <- diag(V)
  Theta0 <- solve(Sigma0 + lambda * I)
  - Theta0 %*% DeltaSigma %*% Theta0
}

#' Partial-correlation delta approximation
#' @param Theta0 baseline precision matrix
#' @param DeltaTheta precision delta matrix
#' @export
mppi_delta_partial <- function(Theta0, DeltaTheta) {
  if (!is.matrix(Theta0) || !is.matrix(DeltaTheta) || nrow(Theta0) != ncol(Theta0) || !identical(dim(Theta0), dim(DeltaTheta)) || any(!is.finite(Theta0)) || any(!is.finite(DeltaTheta)) || any(diag(Theta0) <= 0)) stop("precision matrices must be finite, square, conformable, and have positive baseline diagonal")
  d <- diag(Theta0)
  result <- (-DeltaTheta + 0.5 * Theta0 * (outer(diag(DeltaTheta) / d, rep(1, length(d))) + outer(rep(1, length(d)), diag(DeltaTheta) / d))) / sqrt(outer(d, d))
  diag(result) <- 0
  result
}
