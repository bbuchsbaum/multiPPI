# Internal helpers ---------------------------------------------------------

# Residualize a matrix Y on design Z using QR
.mppi_residualize <- function(Y, Z) {
  stopifnot(is.matrix(Y), is.matrix(Z), nrow(Y) == nrow(Z))
  qr_obj <- qr(Z)
  Y - qr.fitted(qr_obj, Y)
}

# Residualize a vector pk on Q using QR
.mppi_residualize_vec <- function(pk, Q) {
  stopifnot(is.numeric(pk), is.matrix(Q), length(pk) == nrow(Q))
  qr_obj <- qr(Q)
  as.numeric(pk - qr.fitted(qr_obj, pk))
}

# Center a vector within runs (if provided)
.mppi_center_by_run <- function(pk, runs) {
  if (is.null(runs)) return(pk)
  stopifnot(length(pk) == length(runs))
  split_idx <- split(seq_along(pk), runs)
  for (idx in split_idx) {
    pk[idx] <- pk[idx] - mean(pk[idx])
  }
  pk
}

.mppi_check_basis <- function(basis, V_dim) {
  if (is.null(basis)) return(NULL)
  if (inherits(basis, "mppi_basis")) {
    V <- basis$V
    if (nrow(V) != V_dim) stop(sprintf("basis$V must have %d rows to match data columns (got %d).", V_dim, nrow(V)))
    return(basis)
  }
  V <- if (is.list(basis) && !is.null(basis$V)) as.matrix(basis$V) else as.matrix(basis)
  if (nrow(V) != V_dim) stop(sprintf("basis$V must have %d rows to match data columns (got %d).", V_dim, nrow(V)))
  name <- if (is.list(basis) && !is.null(basis$name)) basis$name else "basis"
  # Orthonormalize columns (QR) for stability
  qrV <- qr(V)
  Q <- qr.Q(qrV)
  rank <- qrV$rank
  Q <- Q[, seq_len(rank), drop = FALSE]
  structure(list(V = Q, name = name, r = ncol(Q)), class = "mppi_basis")
}

.mppi_project_basis <- function(R, V, backend = c("blas","chunked"), chunk_cols = NULL) {
  backend <- match.arg(backend)
  if (backend == "blas" || is.null(chunk_cols)) {
    return(R %*% V)
  }
  if (chunk_cols <= 0) chunk_cols <- min(2048L, ncol(R))
  Z <- matrix(0, nrow(R), ncol(V))
  starts <- seq.int(1L, ncol(R), by = chunk_cols)
  for (st in starts) {
    ed <- min(ncol(R), st + chunk_cols - 1L)
    idx <- st:ed
    Z <- Z + R[, idx, drop = FALSE] %*% V[idx, , drop = FALSE]
  }
  Z
}

# Handle missing values across Y and X according to na_action
.mppi_handle_na <- function(Y, X, runs = NULL, na_action = c("omit_tr", "error")) {
  na_action <- match.arg(na_action)
  na_y <- which(rowSums(!is.na(Y)) != ncol(Y))
  na_x <- which(rowSums(!is.na(X)) != ncol(X))
  bad <- sort(unique(c(na_y, na_x)))
  n_drop <- length(bad)
  if (n_drop == 0) {
    return(list(Y = Y, X = X, runs = runs, dropped = 0L, dropped_idx = integer(0)))
  }
  if (na_action == "error") {
    stop(sprintf("mPPI inputs contain NA values in %d time points.", n_drop))
  }
  keep <- setdiff(seq_len(nrow(Y)), bad)
  list(
    Y = Y[keep, , drop = FALSE],
    X = X[keep, , drop = FALSE],
    runs = if (is.null(runs)) NULL else runs[keep],
    dropped = n_drop,
    dropped_idx = bad
  )
}

# Weighted crossproduct: R' diag(w) R
.mppi_wcp <- function(R, w) {
  stopifnot(is.matrix(R))
  w_vec <- as.numeric(w)
  stopifnot(length(w_vec) == nrow(R))
  crossprod(R, w_vec * R)
}

.mppi_crossprod <- function(R, pk, backend = c("blas", "accumulate", "chunked"), chunk_size = NULL) {
  backend <- match.arg(backend)
  if (backend == "blas") {
    return(.mppi_wcp(R, pk))
  }
  V <- ncol(R)
  D <- matrix(0, V, V)
  if (backend == "accumulate") {
    for (t in seq_len(nrow(R))) {
      rt <- R[t, , drop = FALSE]
      D <- D + pk[t] * crossprod(rt)
    }
  } else {
    if (is.null(chunk_size) || chunk_size <= 0) chunk_size <- min(2048L, nrow(R))
    starts <- seq.int(1L, nrow(R), by = chunk_size)
    for (st in starts) {
      ed <- min(st + chunk_size - 1L, nrow(R))
      idx <- st:ed
      Rblock <- R[idx, , drop = FALSE]
      D <- D + crossprod(Rblock, pk[idx] * Rblock)
    }
  }
  D <- (D + t(D)) / 2
  D
}

.mppi_crossprod_lagged <- function(R, pk, lag, blocklens = NULL,
                                    backend = c("blas", "accumulate", "chunked"),
                                    chunk_size = NULL) {
  lag <- as.integer(lag)
  if (lag == 0L) return(.mppi_crossprod(R, pk, match.arg(backend), chunk_size))
  backend <- match.arg(backend)
  Tn <- nrow(R)
  V <- ncol(R)
  if (abs(lag) >= Tn) {
    return(matrix(NA_real_, V, V))
  }
  if (is.null(blocklens)) {
    if (lag > 0L) {
      idx <- seq_len(Tn - lag)
      R0 <- R[idx, , drop = FALSE]
      R1 <- R[idx + lag, , drop = FALSE]
      wk <- pk[idx]
    } else {
      shift <- abs(lag)
      idx <- seq_len(Tn - shift)
      R0 <- R[idx + shift, , drop = FALSE]
      R1 <- R[idx, , drop = FALSE]
      wk <- pk[idx]
    }
    denom <- sum(wk^2)
    if (denom < .Machine$double.eps) {
      return(matrix(NA_real_, V, V))
    }
    return(crossprod(R0, wk * R1) / denom)
  }
  blocklens <- as.integer(blocklens)
  accum <- matrix(0, V, V)
  denom_total <- 0
  start <- 1L
  for (L in blocklens) {
    if (lag > 0L) {
      if (L > lag) {
        rng <- seq.int(start, start + L - lag - 1L)
        R0 <- R[rng, , drop = FALSE]
        R1 <- R[rng + lag, , drop = FALSE]
        wk <- pk[rng]
        accum <- accum + crossprod(R0, wk * R1)
        denom_total <- denom_total + sum(wk^2)
      }
    } else {
      shift <- abs(lag)
      if (L > shift) {
        rng <- seq.int(start + shift, start + L - 1L)
        R0 <- R[rng, , drop = FALSE]
        R1 <- R[rng - shift, , drop = FALSE]
        wk <- pk[rng - shift]
        accum <- accum + crossprod(R0, wk * R1)
        denom_total <- denom_total + sum(wk^2)
      }
    }
    start <- start + L
  }
  if (denom_total < .Machine$double.eps) {
    return(matrix(NA_real_, V, V))
  }
  accum / denom_total
}

.mppi_pack_upper <- function(M) {
  stopifnot(is.matrix(M), nrow(M) == ncol(M))
  V <- ncol(M)
  idx <- upper.tri(M, diag = TRUE)
  list(values = M[idx], dim = V)
}

.mppi_unpack_upper <- function(packed) {
  V <- packed$dim
  M <- matrix(0, V, V)
  idx <- upper.tri(M, diag = TRUE)
  M[idx] <- packed$values
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  M
}

.mppi_row_outer_packed <- function(R) {
  V <- ncol(R)
  idx <- upper.tri(matrix(0, V, V), diag = TRUE)
  P <- sum(idx)
  out <- matrix(0, nrow(R), P)
  for (t in seq_len(nrow(R))) {
    Rt <- R[t, , drop = FALSE]
    out[t, ] <- (crossprod(Rt))[idx]
  }
  out
}

.mppi_frob2 <- function(M) {
  if (is.list(M) && all(c("values","dim") %in% names(M))) {
    sum(M$values^2, na.rm = TRUE)
  } else {
    sum(M^2, na.rm = TRUE)
  }
}

.mppi_phase_randomize_seg <- function(x) {
  n <- length(x)
  if (n <= 2L) return(x)
  Fx <- fft(x)
  half <- floor(n / 2)
  for (k in seq.int(2L, half)) {
    angle <- runif(1, 0, 2 * pi)
    mag <- Mod(Fx[k])
    Fx[k] <- mag * exp(1i * angle)
    conj_idx <- n - k + 2L
    if (conj_idx <= length(Fx)) {
      Fx[conj_idx] <- Conj(Fx[k])
    }
  }
  if (n %% 2 == 0) {
    idx <- half + 1L
    Fx[idx] <- Mod(Fx[idx]) * exp(1i * runif(1, 0, 2 * pi))
  }
  Re(fft(Fx, inverse = TRUE) / n)
}

.mppi_phase_randomize <- function(pk, runs = NULL) {
  if (is.null(runs)) runs <- rep(1L, length(pk))
  stopifnot(length(pk) == length(runs))
  out <- pk
  idx <- split(seq_along(pk), runs)
  for (g in idx) {
    out[g] <- .mppi_phase_randomize_seg(pk[g])
  }
  out
}

# Internal: detect if a design already contains an intercept-like column
.mppi_has_intercept <- function(X) {
  stopifnot(is.matrix(X))
  if ("(Intercept)" %in% colnames(X)) return(TRUE)
  any(apply(X, 2, function(z) {
    sd_z <- stats::sd(z)
    (is.na(sd_z) || sd_z < .Machine$double.eps) && abs(mean(z) - 1) < 1e-8
  }))
}

# Prepare a design matrix, ensuring an intercept column and adjusting indices
.mppi_prepare_design_matrix <- function(X, psych_idx = NULL) {
  stopifnot(is.matrix(X))
  if (.mppi_has_intercept(X)) {
    return(list(X = X,
                psych_idx = if (is.null(psych_idx)) psych_idx else sort(unique(as.integer(psych_idx))),
                intercept_added = FALSE))
  }
  X_new <- cbind(`(Intercept)` = 1, X)
  idx_new <- if (is.null(psych_idx)) psych_idx else sort(unique(as.integer(psych_idx))) + 1L
  list(X = X_new, psych_idx = idx_new, intercept_added = TRUE)
}

# Backwards-compatible helper used by legacy code paths
.mppi_with_intercept <- function(X) {
  .mppi_prepare_design_matrix(X)$X
}

# Public helper: select psychological regressors by name or regex
#' Select psychological regressors from a design matrix
#' @param X design matrix (T x q)
#' @param patterns character vector of regex patterns (applied to colnames(X))
#' @param names_or_idx specific names or integer indices
#' @return integer vector of column indices
mppi_select_psych <- function(X, patterns = NULL, names_or_idx = NULL) {
  stopifnot(is.matrix(X))
  if (!is.null(names_or_idx)) {
    if (is.character(names_or_idx)) return(match(names_or_idx, colnames(X)))
    if (is.numeric(names_or_idx))   return(names_or_idx)
  }
  if (is.null(patterns)) stop("Provide either names_or_idx or patterns to select psychological regressors.")
  keep <- unique(unlist(lapply(patterns, function(p) grep(p, colnames(X), perl = TRUE))))
  if (length(keep) == 0) stop("No design columns matched your patterns.")
  keep
}

# Shift a vector within run blocks by an integer lag
.mppi_shift_by_run <- function(x, lag, blocklens) {
  if (lag == 0 || length(blocklens) == 0) return(x)
  out <- x
  start <- 1L
  for (L in blocklens) {
    rng <- start:(start + L - 1L)
    seg <- x[rng]
    if (abs(lag) >= L) {
      out[rng] <- 0
    } else if (lag > 0) {
      out[rng] <- c(rep(0, lag), seg[1:(L - lag)])
    } else {
      out[rng] <- c(seg[(1 - lag):L], rep(0, -lag))
    }
    start <- start + L
  }
  out
}
