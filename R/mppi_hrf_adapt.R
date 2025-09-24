# HRF-adaptive combination -------------------------------------------------

#' Combine multiple ΔΣ from an HRF basis using energy-maximizing weights
#' @param deltas_list list of VxV slope matrices for the same psychological effect but different HRF basis columns
#' @param pk_norms optional vector of ||pk||^2 to normalize columns
#' @param normalize logical; if TRUE, normalize Gram by pk_norms
#' @return list with Delta (combined), weights, Gram matrix M
mppi_hrf_adapt <- function(deltas_list, pk_norms = NULL, normalize = TRUE) {
  stopifnot(is.list(deltas_list), length(deltas_list) >= 1)
  B <- length(deltas_list)
  M <- matrix(0, B, B)
  for (i in seq_len(B)) for (j in i:B) {
    v <- sum(deltas_list[[i]] * deltas_list[[j]], na.rm = TRUE)
    M[i,j] <- M[j,i] <- v
  }
  if (normalize && !is.null(pk_norms)) {
    D <- diag(1/sqrt(pk_norms), nrow = B, ncol = B)
    M <- D %*% M %*% D
  }
  ev <- eigen(M, symmetric = TRUE)
  a  <- ev$vectors[,1,drop = TRUE]
  if (normalize && !is.null(pk_norms)) a <- a / sqrt(pk_norms)
  Delta_star <- Reduce(`+`, Map(function(w, D) w*D, a, deltas_list))
  list(Delta = Delta_star, weights = a, M = M)
}
