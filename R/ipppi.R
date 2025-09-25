# Instantaneous estimators --------------------------------------------------

ipppi_ewm <- function(U, pk, tau_half, offsets,
                      normalized = FALSE, blocklens = NULL) {
  stopifnot(is.matrix(U), is.matrix(pk))
  if (nrow(U) != nrow(pk)) stop("U and pk must share the same number of rows.", call. = FALSE)
  eta <- log(2) / tau_half
  if (!is.null(blocklens)) {
    blocklens <- as.integer(blocklens)
    if (any(!is.finite(blocklens)) || any(blocklens <= 0L)) {
      stop("blocklens must contain positive run lengths.", call. = FALSE)
    }
  }
  res <- ewm_gated_cross_rcpp(U, pk, eta, as.integer(offsets), normalized,
                              if (is.null(blocklens)) integer() else blocklens)
  K <- attr(res, "K") %||% ncol(pk)
  L <- attr(res, "L") %||% length(offsets)
  offsets <- attr(res, "offsets") %||% as.integer(offsets)
  out <- vector("list", K)
  for (k in seq_len(K)) {
    offset_list <- vector("list", L)
    for (li in seq_len(L)) {
      idx <- (k - 1L) * L + li
      offset_list[[li]] <- res[[idx]]
    }
    names(offset_list) <- as.character(offsets)
    out[[k]] <- offset_list
  }
  names(out) <- colnames(pk)
  list(mats = out, offsets = as.integer(offsets))
}

ipppi_bandpass_fft <- function(U, fs, band) {
  U <- as.matrix(U)
  T <- nrow(U)
  r <- ncol(U)
  freq <- (seq_len(T) - 1) * fs / T
  keep <- (freq >= band[1] & freq <= band[2]) |
    (freq >= fs - band[2] & freq <= fs - band[1])
  out <- matrix(0, T, r)
  for (j in seq_len(r)) {
    spec <- fft(U[, j])
    spec[!keep] <- 0
    out[, j] <- Re(fft(spec, inverse = TRUE) / T)
  }
  scale(out, center = TRUE, scale = TRUE)
}

ipppi_analytic_signal <- function(x) {
  N <- length(x)
  X <- fft(x)
  h <- rep(0, N)
  if (N %% 2 == 0) {
    h[1] <- 1
    h[N / 2 + 1] <- 1
    h[2:(N / 2)] <- 2
  } else {
    h[1] <- 1
    h[2:((N + 1) / 2)] <- 2
  }
  fft(X * h, inverse = TRUE) / N
}

ipppi_iphs_phase <- function(U, fs, band) {
  Ubp <- ipppi_bandpass_fft(U, fs, band)
  apply(Ubp, 2, function(col) Arg(ipppi_analytic_signal(col)))
}

#' Instantaneous phase synchrony aggregator
#'
#' @keywords internal
ipppi_iphs <- function(U, pk, fs, band, tau_half = 10) {
  stopifnot(is.matrix(U), is.matrix(pk))
  if (nrow(U) != nrow(pk)) stop("U and pk must share the same number of rows.", call. = FALSE)
  phi <- ipppi_iphs_phase(U, fs, band)
  T <- nrow(phi)
  r <- ncol(phi)
  K <- ncol(pk)
  eta <- log(2) / tau_half
  const_a <- exp(-eta)
  one_minus_a <- 1 - const_a
  out <- vector("list", K)
  for (k in seq_len(K)) {
    M <- matrix(0, r, r)
    for (t in seq_len(T)) {
      g <- pk[t, k]
      if (is.na(g)) {
        M <- const_a * M
        next
      }
      ct <- cos(outer(phi[t, ], phi[t, ], `-`))
      M <- const_a * M + one_minus_a * g * ct
    }
    out[[k]] <- M
  }
  if (!is.null(colnames(pk))) names(out) <- colnames(pk)
  out
}

#' Exponentially weighted instantaneous correlation
#'
#' @param x Numeric vector (time series).
#' @param y Numeric vector (time series) of the same length as `x`.
#' @param tau_half Half-life (in samples) of the exponential kernel.
#' @param offset Optional lag applied to `y` (positive values shift `y` forward).
#' @param warmup Optional number of initial samples to treat as warm-up (returned as `NA`).
#' @param fill Fill strategy for out-of-range samples when applying `offset`:
#'   `"zero"` (default), `"edge"`, or `"na"`.
#' @return Numeric vector of length `length(x)` containing the streaming correlation.
#' @export
inst_corr <- function(x, y, tau_half, offset = 0L, warmup = NULL, fill = c("zero", "edge", "na")) {
  fill <- match.arg(fill)
  warmup_idx <- if (is.null(warmup)) -1L else as.integer(warmup)
  inst_corr_rcpp(as.numeric(x), as.numeric(y), tau_half = tau_half,
                 offset = as.integer(offset), warmup = warmup_idx,
                 fill = fill)
}
