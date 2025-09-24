# Gain-routing axes ---------------------------------------------------------

#' Summarise conditions along gain and routing axes
#'
#' @param fit `mppi_fit` object
#' @param lags integer vector of lags (excluding zero) used for routing; defaults to stored lags
#' @param mode scale view for the interaction matrices
#' @param use_cached reuse stored lagged matrices when available
#' @return data.frame with condition name, gain, and routing indices
#' @export
mppi_axes <- function(fit, lags = NULL, mode = c("normalized", "amplitude", "raw"), use_cached = TRUE) {
  stopifnot(inherits(fit, "mppi_fit"))
  mode <- match.arg(mode)
  cond_names <- fit$names
  if (is.null(cond_names)) cond_names <- paste0("k", seq_along(fit$Delta))
  Sigma0 <- fit$Sigma0
  if (is.null(Sigma0)) {
    stop("Baseline covariance (Sigma0) not stored in fit; refit with scale='cov'.", call. = FALSE)
  }
  Sigma_sym <- 0.5 * (Sigma0 + t(Sigma0))
  eig_base <- eigen(Sigma_sym, symmetric = TRUE)
  v1 <- eig_base$vectors[, 1]
  if (is.null(lags)) {
    if (!is.null(fit$lags)) {
      lags <- setdiff(sort(unique(as.integer(fit$lags))), 0L)
    } else {
      lags <- integer(0)
    }
  } else {
    lags <- setdiff(sort(unique(as.integer(lags))), 0L)
  }
  eps <- 1e-12
  rows <- vector("list", length(cond_names))

  get_lag_mats <- function(k_idx) {
    if (!length(lags)) return(list())
    if (use_cached && !is.null(fit$lagged) && !is.null(fit$lags) && all(lags %in% as.integer(fit$lags))) {
      lag_store <- fit$lagged[[k_idx]]
      if (is.null(lag_store)) lag_store <- fit$lagged[[cond_names[k_idx]]]
      if (is.null(lag_store)) return(list())
      lapply(as.character(lags), function(lg) lag_store[[lg]])
    } else {
      lagged <- mppi_lagged(fit, k = k_idx, lags = sort(unique(c(0L, lags))))
      lapply(as.character(lags), function(lg) lagged$M_lag[[lg]])
    }
  }

  for (k_idx in seq_along(cond_names)) {
    Mk <- mppi_get_M_scaled(fit, k_idx, mode = mode)
    Sk <- 0.5 * (Mk + t(Mk))
    fnorm <- sqrt(sum(Sk^2, na.rm = TRUE))
    gain <- if (fnorm < eps) NA_real_ else as.numeric((t(v1) %*% Sk %*% v1) / fnorm)
    if (!is.na(gain)) gain <- max(min(gain, 1), -1)

    lag_mats <- get_lag_mats(k_idx)
    if (length(lag_mats)) {
      pos_sum <- 0
      neg_sum <- 0
      for (i in seq_along(lags)) {
        lg <- lags[i]
        Mk_lag <- lag_mats[[i]]
        if (is.null(Mk_lag)) next
        e <- sqrt(sum(Mk_lag^2, na.rm = TRUE))
        if (lg > 0) pos_sum <- pos_sum + e
        if (lg < 0) neg_sum <- neg_sum + e
      }
      routing <- (pos_sum - neg_sum) / (pos_sum + neg_sum + eps)
    } else {
      routing <- NA_real_
    }
    rows[[k_idx]] <- data.frame(condition = cond_names[k_idx],
                                gain = gain,
                                routing = routing,
                                scale = mode,
                                stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
