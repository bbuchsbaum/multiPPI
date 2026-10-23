# Instantaneous estimators --------------------------------------------------

# Final state of the exponentially weighted slope filter, in closed form:
# sum_t w_t p_t U[t, i] U[t + s, j] / sum_t w_t p_t^2 with w_t = 2^((t - t_end) / tau_half),
# transposed for s < 0 (the closed engine's lag convention). Rows sit at acquisition
# positions `index`; pairs must be s scans apart there and share a segment.
ipppi_ewm <- function(U, pk, tau_half, offsets,
                      index = seq_len(nrow(U)), segments = rep(1L, nrow(U))) {
  stopifnot(is.matrix(U), is.matrix(pk), length(index) == nrow(U), length(segments) == nrow(U))
  if (nrow(U) != nrow(pk)) stop("U and pk must share the same number of rows.", call. = FALSE)
  offsets <- as.integer(offsets)
  w <- 2^((index - max(index)) / tau_half)
  pos <- integer(max(index)); pos[index] <- seq_along(index)
  pairs <- lapply(offsets, function(s) {
    b <- pos[index + abs(s)]
    a <- which(!is.na(b) & b > 0L)
    b <- b[a]
    ok <- segments[a] == segments[b]
    list(a = a[ok], b = b[ok])
  })
  V <- ncol(U)
  out <- lapply(seq_len(ncol(pk)), function(k) {
    mats <- Map(function(s, pr) {
      a <- pr$a; b <- pr$b
      wp <- w[a] * pk[a, k]
      denom <- sum(wp * pk[a, k])
      if (!is.finite(denom) || denom < .Machine$double.eps) return(matrix(NA_real_, V, V))
      if (s >= 0L) crossprod(U[a, , drop = FALSE], wp * U[b, , drop = FALSE]) / denom
      else crossprod(U[b, , drop = FALSE], wp * U[a, , drop = FALSE]) / denom
    }, offsets, pairs)
    names(mats) <- as.character(offsets)
    mats
  })
  names(out) <- colnames(pk)
  list(mats = out, offsets = offsets)
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

# Labels of maximal runs of rows that are consecutive scans of one segment.
.ipppi_contiguous <- function(index, segments) {
  cumsum(c(TRUE, diff(index) != 1L | diff(segments) != 0L))
}

# Band-limited Hilbert phase, computed within each contiguous segment (NA in
# single-row segments).
ipppi_iphs_phase <- function(U, fs, band, segments = rep(1L, nrow(U))) {
  phi <- matrix(NA_real_, nrow(U), ncol(U))
  for (idx in split(seq_len(nrow(U)), segments)) {
    if (length(idx) < 2L) next
    Ubp <- ipppi_bandpass_fft(U[idx, , drop = FALSE], fs, band)
    phi[idx, ] <- apply(Ubp, 2, function(col) Arg(ipppi_analytic_signal(col)))
  }
  phi
}

#' Instantaneous phase synchrony aggregator
#'
#' End-weighted slope of instantaneous phase locking on each psychological
#' regressor, \eqn{\sum_t w_t p_t \cos(\phi_{t,i} - \phi_{t,j}) / \sum_t w_t p_t^2}
#' with \eqn{w_t = 2^{(t - t_{end})/h}} on the acquisition grid `index`, as for
#' the EWM estimator. Phases are band-limited Hilbert phases computed within
#' each contiguous segment (consecutive `index` values sharing a `segments`
#' label); rows without a phase are skipped pairwise. With `pk = 1` this is the
#' end-weighted mean phase locking.
#'
#' @keywords internal
ipppi_iphs <- function(U, pk, fs, band, tau_half = 10,
                       index = seq_len(nrow(U)), segments = rep(1L, nrow(U))) {
  stopifnot(is.matrix(U), is.matrix(pk), length(index) == nrow(U), length(segments) == nrow(U))
  if (nrow(U) != nrow(pk)) stop("U and pk must share the same number of rows.", call. = FALSE)
  phi <- ipppi_iphs_phase(U, fs, band, .ipppi_contiguous(index, segments))
  ok <- is.finite(phi) + 0
  C <- ifelse(ok > 0, cos(phi), 0)
  S <- ifelse(ok > 0, sin(phi), 0)
  w <- 2^((index - max(index)) / tau_half)
  out <- lapply(seq_len(ncol(pk)), function(k) {
    p <- pk[, k]
    p[!is.finite(p)] <- 0
    wp <- w * p
    denom <- crossprod(ok, wp * p * ok)
    M <- (crossprod(C, wp * C) + crossprod(S, wp * S)) / denom
    M[!is.finite(denom) | denom < .Machine$double.eps] <- NA_real_
    M
  })
  if (!is.null(colnames(pk))) names(out) <- colnames(pk)
  out
}

.ipppi_unwrap <- function(theta) {
  y <- as.numeric(theta)
  if (!length(y)) return(y)
  dy <- diff(y)
  dy <- (dy + pi) %% (2 * pi) - pi
  c(y[1], y[1] + cumsum(dy))
}

ipppi_instfreq <- function(U, fs, band, runs = NULL) {
  stopifnot(is.matrix(U), length(fs) == 1L, is.finite(fs), fs > 0)
  if (length(band) != 2L || any(!is.finite(band))) {
    stop("band must be a finite length-2 numeric vector.", call. = FALSE)
  }
  band <- sort(as.numeric(band))
  if (diff(band) <= 0) stop("band must satisfy band[1] < band[2].", call. = FALSE)
  U <- as.matrix(U)
  Tn <- nrow(U)
  V <- ncol(U)
  if (is.null(runs)) runs <- rep.int(1L, Tn)
  if (length(runs) != Tn) stop("runs must have length nrow(U).", call. = FALSE)
  if (!is.atomic(runs) || anyNA(runs)) stop("runs must be non-missing run labels.", call. = FALSE)
  runs <- as.character(runs)
  # filtering, Hilbert transform and scaling are run-local
  seg <- cumsum(c(TRUE, runs[-1L] != runs[-Tn]))
  phi <- amp <- omega <- matrix(NA_real_, Tn, V)
  for (idx in split(seq_len(Tn), seg)) {
    if (length(idx) < 2L) {
      next
    }
    Ubp <- ipppi_bandpass_fft(U[idx, , drop = FALSE], fs, band)
    for (j in seq_len(V)) {
      analytic <- ipppi_analytic_signal(Ubp[, j])
      amp[idx, j] <- Mod(analytic)
      phi[idx, j] <- Arg(analytic)
      omega[idx, j] <- fs * c(NA_real_, diff(.ipppi_unwrap(phi[idx, j])))
    }
  }
  list(phi = phi, amp = amp, omega = omega)
}

ipppi_pwi <- function(U, pk, fs, band, runs = NULL,
                      weight = c("amplitude", "unit"), robust = FALSE, tol = 1e-8) {
  stopifnot(is.matrix(U))
  pk <- as.numeric(pk)
  if (nrow(U) != length(pk)) stop("pk must have length nrow(U).", call. = FALSE)
  weight <- match.arg(weight)
  fi <- ipppi_instfreq(U, fs = fs, band = band, runs = runs)
  Tn <- nrow(U)
  V <- ncol(U)
  pk_vec <- pk
  weighted_beta <- function(x, w) {
    keep <- is.finite(x) & is.finite(pk_vec) & is.finite(w)
    if (!any(keep)) return(NA_real_)
    x <- x[keep]
    p <- pk_vec[keep]
    w <- w[keep]
    if (!any(w > 0)) return(NA_real_)
    wm <- sum(w)
    if (!is.finite(wm) || wm <= 0) return(NA_real_)
    pm <- sum(w * p) / wm
    xm <- sum(w * x) / wm
    num <- sum(w * (p - pm) * (x - xm))
    den <- sum(w * (p - pm)^2)
    if (!is.finite(den) || den <= tol) return(NA_real_)
    beta <- num / den
    if (!robust) return(beta)
    kappa <- 1.345
    for (it in seq_len(10L)) {
      r <- x - (xm + beta * (p - pm))
      s0 <- stats::mad(r, center = 0, constant = 1.4826, na.rm = TRUE)
      if (!is.finite(s0) || s0 <= tol) break
      wh <- pmin(1, (kappa * s0) / pmax(abs(r), .Machine$double.eps))
      w2 <- w * wh
      if (!any(w2 > 0)) break
      wm2 <- sum(w2)
      if (!is.finite(wm2) || wm2 <= 0) break
      pm <- sum(w2 * p) / wm2
      xm <- sum(w2 * x) / wm2
      num <- sum(w2 * (p - pm) * (x - xm))
      den <- sum(w2 * (p - pm)^2)
      if (!is.finite(den) || den <= tol) break
      beta_new <- num / den
      if (!is.finite(beta_new)) break
      if (abs(beta_new - beta) < 1e-8) {
        beta <- beta_new
        break
      }
      beta <- beta_new
    }
    beta
  }
  w_mat <- if (weight == "amplitude") fi$amp else matrix(1, Tn, V)
  betas <- vapply(seq_len(V), function(j) weighted_beta(fi$omega[, j], w_mat[, j]), numeric(1))
  mat <- outer(betas, betas, `-`)
  diag(mat) <- 0
  list(mat = mat, beta = betas, divergence = rowSums(mat, na.rm = FALSE), omega = fi$omega)
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
#' @return Numeric vector of length `length(x)` containing the streaming
#'   (uncentered) correlation. Samples where `x` or the shifted `y` is not finite
#'   (`NA`, `NaN`, `Inf`) are skipped: they return `NA` and only decay the
#'   running moments. Values are `NA` while either running moment is zero.
#' @export
inst_corr <- function(x, y, tau_half, offset = 0L, warmup = NULL, fill = c("zero", "edge", "na")) {
  fill <- match.arg(fill)
  warmup_idx <- if (is.null(warmup)) -1L else as.integer(warmup)
  inst_corr_rcpp(as.numeric(x), as.numeric(y), tau_half = tau_half,
                 offset = as.integer(offset), warmup = warmup_idx,
                 fill = fill)
}
