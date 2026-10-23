# Frequency-specific variants ---------------------------------------------

#' Simple bandpass using signal::butter if available, else a DCT-domain filter
#'
#' With the signal package, a zero-phase Butterworth filter (`filtfilt`);
#' otherwise an ideal filter that keeps the DCT components whose frequencies
#' lie in `[low, high]`. Each run is filtered separately.
#' @param x T x V matrix
#' @param low highpass cutoff in (0, 0.5) cycles/TR; NULL for none
#' @param high lowpass cutoff in (0, 0.5) cycles/TR; NULL for none
#' @param order filter order for Butterworth if 'signal' present
#' @param runs optional run labels; filtering does not cross run boundaries
mppi_bandpass <- function(x, low = NULL, high = 0.1, order = 2L, runs = NULL) {
  stopifnot(is.matrix(x))
  if (is.null(low) && is.null(high)) return(x)
  ok <- function(f) is.null(f) || (length(f) == 1L && is.finite(f) && f > 0 && f < 0.5)
  if (!ok(low) || !ok(high)) stop("Band cutoffs must lie in (0, 0.5) cycles/TR.", call. = FALSE)
  if (!is.null(low) && !is.null(high) && low >= high) stop("Band requires low < high.", call. = FALSE)
  if (anyNA(x)) stop("Bandpass filtering requires complete data.", call. = FALSE)
  Tn <- nrow(x)
  if (is.null(runs)) runs <- rep.int(1L, Tn)
  if (length(runs) != Tn || anyNA(runs)) stop("runs must label every row of x.", call. = FALSE)
  use_signal <- requireNamespace("signal", quietly = TRUE)
  if (use_signal) {
    # butter() takes cutoffs relative to Nyquist (0.5 cycles/TR)
    bf <- if (is.null(low)) signal::butter(order, 2 * high, type = "low")
          else if (is.null(high)) signal::butter(order, 2 * low, type = "high")
          else signal::butter(order, 2 * c(low, high), type = "pass")
  }
  seg <- cumsum(c(TRUE, as.character(runs[-1L]) != as.character(runs[-Tn])))
  y <- x
  for (idx in split(seq_len(Tn), seg)) {
    xi <- x[idx, , drop = FALSE]
    y[idx, ] <- if (use_signal) {
      apply(xi, 2, function(col) as.numeric(signal::filtfilt(bf, col)))
    } else {
      .mppi_dct_filter(xi, low %||% 0, high %||% 0.5)
    }
  }
  y
}

# Ideal band filter in the DCT domain: the DFT of the half-sample symmetric
# extension, whose bin j carries DCT component j at j / (2n) cycles/sample.
.mppi_dct_filter <- function(x, low, high) {
  n <- nrow(x)
  f <- c(0:n, rev(seq_len(n - 1L))) / (2 * n)
  keep <- f >= low & f <= high
  apply(x, 2, function(col) {
    Re(stats::fft(stats::fft(c(col, rev(col))) * keep, inverse = TRUE))[seq_len(n)] / (2 * n)
  })
}

#' Fit ΔΣ per frequency band and optionally combine
#'
#' As in `mppi_fit()`, `Y` is residualized on the full design `X` before the
#' residuals are band-limited (within runs), so task main effects cannot leak
#' into ΔΣ. The band-limited residuals are then fitted with the unfiltered `X`,
#' whose psychological regressors weight them as in the broadband fit.
#' @param Y T x V data matrix
#' @param X T x q design matrix
#' @param psych_idx psychological regressor indices
#' @param bands list of band definitions, each with `low`, `high`, and `name`
#' @param zero_diag logical; zero the diagonal of ΔΣ for each band
#' @param scale summary scale for returned fits (`"normalized"`, `"cov"`, or `"corr"`)
#' @param runs optional run labels, used for filtering and passed to `mppi_fit()`
mppi_fit_multi_band <- function(Y, X, psych_idx, bands,
                                zero_diag = TRUE, scale = c("normalized","cov","corr"),
                                runs = NULL) {
  scale_choice <- match.arg(scale)
  stopifnot(is.matrix(Y), is.matrix(X), nrow(X) == nrow(Y))
  if (anyNA(Y) || anyNA(X)) stop("Band-limited fits require complete data.", call. = FALSE)
  R <- .mppi_residualize(Y, .mppi_prepare_design_matrix(X)$X)
  res <- list()
  for (bd in bands) {
    Rb <- mppi_bandpass(R, low = bd$low, high = bd$high, runs = runs)
    fit <- mppi_fit(Rb, X, psych_idx, runs = runs, zero_diag = zero_diag, scale = scale_choice)
    res[[bd$name]] <- fit
  }
  res
}
