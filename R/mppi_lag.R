# Lagged / directed variants ----------------------------------------------

#' Psychological-regressor shifts for covariance modulation
#' @param x object containing a design or fit to dispatch on
#' @param ... additional arguments passed to methods
#' @export
mppi_lagged <- function(x, ...) UseMethod("mppi_lagged")

#' @rdname mppi_lagged
#' @param x T x V data matrix
#' @param X T x q design matrix
#' @param psych_idx indices of psychological main effects
#' @param lags integer vector of lags (e.g., -2:2). Positive values delay `pk`
#' @param blocklens integer vector of run lengths (TR units) to prevent cross-run spillover
#' @param zero_diag logical; zero the diagonal of returned matrices
#' @param scale compute covariance (`"cov"`) or correlation (`"corr"`) slopes
#' @export
mppi_lagged.default <- function(x, X, psych_idx, lags = -2:2, blocklens = NULL,
                                zero_diag = TRUE, scale = c("cov","corr"), ...) {
  scale <- match.arg(scale)
  Y <- x
  if (any(!is.finite(Y)) || any(!is.finite(X))) stop("Psychological shifts require a complete acquisition grid.", call. = FALSE)
  if (any(!is.finite(lags)) || any(lags != as.integer(lags))) stop("lags must contain integers.", call. = FALSE)
  if (!is.null(blocklens) && (any(blocklens <= 0) || sum(blocklens) != nrow(Y))) stop("blocklens must describe the acquisition grid.", call. = FALSE)
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
      pk <- .mppi_shift_by_run(pk0, lg, blocklens %||% length(pk0))
      denom <- sum(pk^2)
      Dk <- if (denom < .Machine$double.eps) matrix(NA_real_, ncol(Y), ncol(Y))
            else .mppi_wcp(R, pk) / denom
      if (zero_diag) diag(Dk) <- 0
      out[[paste0(colnames(X)[k], "_lag", lg)]] <- Dk
    }
  }
  attr(out, "lag_estimand") <- "psychological_shift"
  out
}

#' Select lag by maximizing omnibus statistic with permutations
#' @param Y T x V data matrix
#' @param X T x q design matrix
#' @param psych_idx indices of psychological main effects to test
#' @param lags vector of candidate lags
#' @param blocklens optional run lengths to respect when shifting regressors
#' @param B number of permutations
#' @param blksize optional block size used for permutation sign flips
#' @param scale compute covariance (`"cov"`) or correlation (`"corr"`) slopes
#' @param targets optional subset of psychological regressors (names or indices)
#' @param zero_diag logical; zero the diagonal of returned matrices
#' @export
mppi_lag_select <- function(Y, X, psych_idx, lags = -2:2, blocklens = NULL, B = 499L,
                            blksize = 10L, scale = c("cov","corr"), targets = NULL,
                            zero_diag = TRUE) {
  scale <- match.arg(scale)
  if (any(!is.finite(Y)) || any(!is.finite(X))) stop("Lag selection requires a complete acquisition grid.", call. = FALSE)
  if (!length(lags) || any(!is.finite(lags)) || any(lags != as.integer(lags))) stop("lags must be nonempty integers.", call. = FALSE)
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B)) stop("B must be a positive integer.", call. = FALSE)
  if (!is.null(blksize) && (length(blksize) != 1L || !is.finite(blksize) || blksize < 1 || blksize != as.integer(blksize))) stop("blksize must be a positive integer.", call. = FALSE)
  n <- nrow(Y)
  if (is.null(blocklens)) blocklens <- n
  if (any(!is.finite(blocklens)) || any(blocklens <= 0) || any(blocklens != as.integer(blocklens)) || sum(blocklens) != n) stop("blocklens must describe the acquisition grid.", call. = FALSE)
  runs <- rep(seq_along(blocklens), blocklens)
  # Only residuals and regressors are used; the cheapest fit scale avoids V^3 diagnostics.
  fit <- mppi_fit(Y, X, psych_idx, runs=runs, zero_diag=zero_diag, scale="normalized", na_action="error")
  p_names <- fit$names
  targets_idx <- if (is.null(targets)) seq_along(p_names) else if (is.numeric(targets)) as.integer(targets) else match(targets,p_names)
  if (!length(targets_idx) || anyNA(targets_idx) || any(targets_idx < 1 | targets_idx > length(p_names)) || (is.numeric(targets) && any(targets != targets_idx))) stop("Invalid targets.", call. = FALSE)
  R <- fit$R
  if (scale == "corr") {
    sd <- apply(R, 2, stats::sd); sd[sd == 0] <- 1
    R <- sweep(R, 2, sd, "/")
  }
  segments <- split(seq_len(n), runs)
  blocks <- if (is.null(blksize)) segments else unlist(lapply(segments, function(i) split(i, ceiling(seq_along(i)/blksize))), recursive=FALSE)
  if (length(blocks) < 2L) stop("Lag selection requires at least two randomization blocks.", call. = FALSE)
  evaluate <- function(pk) {
    mats <- lapply(lags, function(lg) {
      shifted <- .mppi_shift_by_run(pk, lg, blocklens)
      den <- sum(shifted^2)
      if (den <= .Machine$double.eps) return(matrix(NA_real_, ncol(R), ncol(R)))
      m <- .mppi_wcp(R,shifted)/den
      if (zero_diag) diag(m) <- 0
      m
    })
    q <- vapply(mats,function(m) if (all(is.na(m))) NA_real_ else sum(m^2),numeric(1))
    list(mats=mats, q=q)
  }
  energy <- .mppi_slope_energy(R, zero_diag)
  results <- setNames(vector("list",length(targets_idx)),p_names[targets_idx])
  for (i in seq_along(targets_idx)) {
    pk <- fit$pk[[targets_idx[i]]]; obs <- evaluate(pk)
    names(obs$q) <- paste0("lag",lags)
    if (all(is.na(obs$q))) stop("No candidate lag has nonzero regressor energy.", call. = FALSE)
    best <- which.max(obs$q); null <- numeric(B)
    for (b in seq_len(B)) {
      signs <- sample(c(-1,1),length(blocks),replace=TRUE); perm <- pk
      for (j in seq_along(blocks)) perm[blocks[[j]]] <- signs[j] * perm[blocks[[j]]]
      q <- vapply(lags, function(lg) {
        shifted <- .mppi_shift_by_run(perm, lg, blocklens)
        if (sum(shifted^2) <= .Machine$double.eps) NA_real_ else energy(shifted)
      }, numeric(1))
      null[b] <- if (all(is.na(q))) NA_real_ else max(q,na.rm=TRUE)
    }
    valid <- is.finite(null)
    results[[i]] <- list(best_lag=lags[best], Q=obs$q[best], Qnull=null,
                         p=if (any(valid)) (1+sum(null[valid]>=obs$q[best]))/(sum(valid)+1) else NA_real_,
                         Qobs=obs$q, Dbest=obs$mats[[best]],
                         estimand="psychological_shift", inference=list(status="experimental", assumptions="Blockwise sign symmetry; calibration is unqualified.", selection="maximum over all candidate lags in every draw"))
  }
  if (length(results)==1L) results[[1L]] else results
}
