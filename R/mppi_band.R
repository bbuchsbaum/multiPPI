# Frequency-specific variants ---------------------------------------------

#' Simple bandpass using stats::filter or signal::butter if available
#' @param x T x V matrix
#' @param low highpass cutoff in (0, 0.5) cycles/TR; NULL for none
#' @param high lowpass cutoff in (0, 0.5); NULL for none
#' @param order filter order for Butterworth if 'signal' present
mppi_bandpass <- function(x, low = NULL, high = 0.1, order = 2L) {
  stopifnot(is.matrix(x))
  if (!is.null(low) || !is.null(high)) {
    if (requireNamespace("signal", quietly = TRUE)) {
      fs <- 1  # sample freq (per TR unit)
      Wc <- c(ifelse(is.null(low), 0, low), ifelse(is.null(high), 0.499, high))
      if (is.null(low)) { # lowpass
        bf <- signal::butter(order, Wc[2], type = "low")
      } else if (is.null(high)) { # highpass
        bf <- signal::butter(order, Wc[1], type = "high")
      } else { # bandpass
        bf <- signal::butter(order, Wc, type = "pass")
      }
      y <- apply(x, 2, function(col) as.numeric(signal::filtfilt(bf, col)))
      return(matrix(y, nrow = nrow(x), ncol = ncol(x)))
    } else {
      # Fallback: simple moving-average detrend / low-pass via running mean
      if (!is.null(high)) {
        k <- max(3L, round(1/(2*high)))
        w <- rep(1/k, k)
        y <- apply(x, 2, function(col) as.numeric(stats::filter(col, w, sides = 2, circular = FALSE)))
        y[is.na(y)] <- 0
        return(matrix(y, nrow = nrow(x), ncol = ncol(x)))
      } else {
        return(x)
      }
    }
  }
  x
}

#' Fit ΔΣ per frequency band and optionally combine
#' @param bands list of lists: each with low, high, name
mppi_fit_multi_band <- function(Y, X, psych_idx, bands,
                                zero_diag = TRUE, scale = c("cov","corr")) {
  res <- list()
  for (bd in bands) {
    Yb <- mppi_bandpass(Y, low = bd$low, high = bd$high)
    fit <- mppi_fit(Yb, X, psych_idx, zero_diag = zero_diag, scale = scale)
    res[[bd$name]] <- fit
  }
  res
}
