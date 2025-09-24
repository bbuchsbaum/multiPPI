# Fusion of time-domain and beta-domain -----------------------------------

#' Precision-weighted fusion of ΔΣ maps from time and trial domains
#' @param D_time VxV matrix
#' @param D_beta VxV matrix
#' @param var_time scalar variance estimate for D_time (use omnibus null var)
#' @param var_beta scalar variance estimate for D_beta
mppi_fuse_time_beta <- function(D_time, D_beta, var_time, var_beta) {
  w_t <- 1/var_time; w_b <- 1/var_beta
  (w_t*D_time + w_b*D_beta) / (w_t + w_b)
}
