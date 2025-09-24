# Partial connectivity / interpretation -----------------------------------

#' First-order precision update from covariance slope
#' @param Sigma0 baseline covariance (T^{-1} R'R) or correlation
#' @param DeltaSigma slope matrix for regressor k
#' @param lambda ridge parameter for invertibility
#' @return matrix ΔΘ
mppi_to_partial <- function(Sigma0, DeltaSigma, lambda = 1e-3) {
  V <- nrow(Sigma0); I <- diag(V)
  Theta0 <- solve(Sigma0 + lambda * I)
  - Theta0 %*% DeltaSigma %*% Theta0
}

#' Partial-correlation delta approximation
#' @param Theta0 baseline precision; @param DeltaTheta precision delta
mppi_delta_partial <- function(Theta0, DeltaTheta) {
  V <- nrow(Theta0)
  d <- 1/sqrt(diag(Theta0))
  S <- d %o% d
  - DeltaTheta * S
}
