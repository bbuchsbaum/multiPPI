# Variance-correlation decomposition --------------------------------------

#' Decompose ΔCov into correlation change and variance terms
#' @param R residualized time series (optionally per-condition bins)
#' @param pk residualized psychological regressor (validated for compatibility; DeltaSigma is supplied)
#' @param DeltaSigma raw covariance slope matrix including its diagonal
#' @return list with Delta_rho (approx), var_terms (matrix), and note
mppi_decompose_variance <- function(R, pk, DeltaSigma) {
  if (!is.matrix(R) || !is.numeric(R) || !nrow(R) || any(!is.finite(R)) || !is.numeric(pk) || length(pk) != nrow(R) || any(!is.finite(pk)) || !is.matrix(DeltaSigma) || !identical(dim(DeltaSigma), c(ncol(R), ncol(R))) || any(!is.finite(DeltaSigma))) stop("R, pk and DeltaSigma must be finite and conformable")
  Sigma0 <- crossprod(R) / nrow(R)
  variance <- diag(Sigma0)
  if (any(variance <= 0)) stop("baseline variances must be positive")
  scale <- sqrt(outer(variance, variance))
  rel <- diag(DeltaSigma) / variance
  variance_terms <- 0.5 * Sigma0 * (outer(rel, rep(1, length(rel))) + outer(rep(1, length(rel)), rel))
  Delta_rho <- (DeltaSigma - variance_terms) / scale
  diag(Delta_rho) <- 0
  list(Delta_rho = Delta_rho, variance_terms = variance_terms,
       Sigma0 = Sigma0, note = "DeltaSigma = (sqrt(diag(Sigma0)) outer sqrt(diag(Sigma0))) * Delta_rho + variance_terms")
}
