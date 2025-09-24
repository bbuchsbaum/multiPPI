# Variance-correlation decomposition --------------------------------------

#' Decompose ΔCov into correlation change and variance terms
#' @param R residualized time series (optionally per-condition bins)
#' @param pk residualized psychological regressor
#' @param DeltaSigma covariance slope matrix
#' @return list with Delta_rho (approx), var_terms (matrix), and note
mppi_decompose_variance <- function(R, pk, DeltaSigma) {
  # Approximate Δrho by standardizing R first
  s <- apply(R, 2, sd); s[s == 0] <- 1
  Rc <- sweep(R, 2, s, "/")
  pk_vec <- as.numeric(pk)
  denom <- sum(pk_vec^2)
  Dcorr <- if (denom < .Machine$double.eps) matrix(0, ncol(R), ncol(R)) else crossprod(Rc, pk_vec * Rc) / denom
  # Variance contributions are residual term
  var_terms <- DeltaSigma - (Dcorr * (s %o% s))
  list(Delta_rho = Dcorr, variance_terms = var_terms,
       note = "DeltaSigma ≈ (s ⊗ s) * Delta_rho + variance_terms")
}
