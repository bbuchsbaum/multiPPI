# Lagged / directed variants ----------------------------------------------

#' Lagged Matrix-PPI (directed)
#' @param Y T x V data
#' @param X T x q design
#' @param psych_idx indices of psychological main effects
#' @param lags integer vector of lags (e.g., -2:2). Positive = delay p_k.
#' @param blocklens integer vector of run lengths (TR units) to prevent cross-run spillover
mppi_lagged <- function(x, ...) UseMethod("mppi_lagged")

#' @rdname mppi_lagged
#' @export
mppi_lagged.default <- function(Y, X, psych_idx, lags = -2:2, blocklens = NULL,
                                zero_diag = TRUE, scale = c("cov","corr")) {
  scale <- match.arg(scale)
  design_prep <- .mppi_prepare_design_matrix(X, psych_idx)
  X <- design_prep$X
  p_idx <- sort(unique(design_prep$psych_idx))
  all_idx <- seq_len(ncol(X))
  base_idx <- setdiff(all_idx, p_idx)
  # Residualize once
  R <- .mppi_residualize(Y, X)
  if (scale == "corr") {
    s <- apply(R, 2, sd); s[s == 0] <- 1
    R <- sweep(R, 2, s, "/")
  }
  out <- list()
  for (ii in seq_along(p_idx)) {
    k  <- p_idx[ii]
    Q  <- X[, c(base_idx, setdiff(p_idx, k)), drop = FALSE]
    pk0 <- .mppi_residualize_vec(X[, k], Q)
    for (lg in lags) {
      pk <- if (is.null(blocklens)) {
        # simple shift without wrap
        if (lg == 0) pk0 else {
          z <- rep(0, length(pk0))
          if (lg > 0) { z[(lg+1):length(z)] <- pk0[1:(length(z)-lg)] }
          else        { z[1:(length(z)+lg)] <- pk0[(1-lg):length(pk0)] }
          z
        }
      } else .mppi_shift_by_run(pk0, lg, blocklens)
      denom <- sum(pk^2)
      Dk <- if (denom < .Machine$double.eps) matrix(NA_real_, ncol(Y), ncol(Y))
            else .mppi_wcp(R, pk) / denom
      if (zero_diag) diag(Dk) <- 0
      out[[paste0(colnames(X)[k], "_lag", lg)]] <- Dk
    }
  }
  out
}

#' Select lag by maximizing omnibus statistic with permutations
#' @param lags vector of candidate lags
#' @param B number of permutations
#' @param targets optional subset of psychological regressors (names or indices)
#' @param zero_diag logical; zero the diagonal of returned matrices
mppi_lag_select <- function(Y, X, psych_idx, lags = -2:2, blocklens = NULL, B = 499L,
                            blksize = 10L, scale = c("cov","corr"), targets = NULL,
                            zero_diag = TRUE) {
  scale <- match.arg(scale)
  fit <- mppi_fit(Y, X, psych_idx, zero_diag = FALSE, scale = scale)
  Tn <- nrow(fit$R)
  V <- ncol(fit$R)
  p_names <- fit$names
  if (length(p_names) == 0) stop("No psychological regressors supplied via 'psych_idx'.")
  pick_targets <- function(x) {
    if (is.null(x)) return(seq_along(p_names))
    if (is.numeric(x)) {
      if (!all(x %in% seq_along(p_names))) stop("'targets' indices out of range.")
      return(x)
    }
    match_idx <- match(x, p_names)
    if (any(is.na(match_idx))) stop("Unable to match targets: ", paste(x[is.na(match_idx)], collapse = ", "))
    match_idx
  }
  tgt_idx <- pick_targets(targets)
  if (!is.null(blocklens)) {
    if (sum(blocklens) != Tn) stop("Sum(blocklens) must equal number of time points.")
  }
  shift_no_block <- function(x, lag) {
    if (lag == 0) return(x)
    n <- length(x)
    out <- numeric(n)
    if (lag > 0) {
      if (lag < n) out[(lag+1):n] <- x[1:(n-lag)]
    } else {
      lag <- abs(lag)
      if (lag < n) out[1:(n-lag)] <- x[(lag+1):n]
    }
    out
  }
  make_perm_blocks <- function() {
    if (!is.null(blksize)) {
      split(seq_len(Tn), ceiling(seq_len(Tn)/blksize))
    } else if (!is.null(blocklens)) {
      ends <- cumsum(blocklens)
      starts <- c(1, head(ends, -1) + 1)
      Map(seq, starts, ends)
    } else {
      list(seq_len(Tn))
    }
  }
  idx_blocks <- make_perm_blocks()
  signflip_block <- function(vec) {
    out <- vec
    sgn <- sample(c(-1, 1), length(idx_blocks), replace = TRUE)
    jj <- 1L
    for (g in idx_blocks) {
      out[g] <- sgn[jj] * out[g]
      jj <- jj + 1L
    }
    out
  }
  results <- vector("list", length(tgt_idx))
  names(results) <- p_names[tgt_idx]
  for (ii in seq_along(tgt_idx)) {
    pk0 <- fit$pk[[tgt_idx[ii]]]
    lag_vals <- setNames(rep(NA_real_, length(lags)), paste0("lag", lags))
    D_per_lag <- vector("list", length(lags))
    for (jj in seq_along(lags)) {
      lg <- lags[jj]
      pk_lag <- if (is.null(blocklens)) shift_no_block(pk0, lg) else .mppi_shift_by_run(pk0, lg, blocklens)
      denom <- sum(pk_lag^2)
      if (denom < .Machine$double.eps) {
        D_per_lag[[jj]] <- matrix(NA_real_, V, V)
        next
      }
      Dk <- .mppi_wcp(fit$R, pk_lag) / denom
      if (zero_diag) diag(Dk) <- 0
      D_per_lag[[jj]] <- Dk
      lag_vals[jj] <- sum(Dk^2, na.rm = TRUE)
    }
    if (all(is.na(lag_vals))) {
      results[[ii]] <- list(best_lag = NA_integer_, Q = NA_real_, Qnull = rep(NA_real_, B),
                            p = NA_real_, Qobs = lag_vals, Dbest = matrix(NA_real_, V, V))
      next
    }
    best_idx <- which.max(lag_vals)
    best_lag <- lags[best_idx]
    Dbest <- D_per_lag[[best_idx]]
    Qbest <- lag_vals[best_idx]
    pk_best <- if (is.null(blocklens)) shift_no_block(pk0, best_lag) else .mppi_shift_by_run(pk0, best_lag, blocklens)
    denom_best <- sum(pk_best^2)
    if (denom_best < .Machine$double.eps) {
      results[[ii]] <- list(best_lag = best_lag, Q = NA_real_, Qnull = rep(NA_real_, B),
                            p = NA_real_, Qobs = lag_vals, Dbest = Dbest)
      next
    }
    Qnull <- rep(NA_real_, B)
    for (b in seq_len(B)) {
      pk_perm <- signflip_block(pk_best)
      denomb <- sum(pk_perm^2)
      if (denomb < .Machine$double.eps) next
      Db <- .mppi_wcp(fit$R, pk_perm) / denomb
      if (zero_diag) diag(Db) <- 0
      Qnull[b] <- sum(Db^2, na.rm = TRUE)
    }
    valid <- !is.na(Qnull)
    pval <- if (!any(valid)) NA_real_ else (1 + sum(Qnull[valid] >= Qbest)) / (sum(valid) + 1)
    results[[ii]] <- list(best_lag = best_lag, Q = Qbest, Qnull = Qnull,
                          p = pval, Qobs = lag_vals, Dbest = Dbest)
  }
  if (length(results) == 1) results[[1]] else results
}
