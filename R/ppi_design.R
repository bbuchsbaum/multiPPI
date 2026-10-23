#' Construct a conventional run-local PPI design
#'
#' Adjust and deconvolve the BOLD seed independently in each run using a fixed-ridge
#' DCT basis, multiply by neural psychological inputs, and reconvolve. Task main
#' effects use the uncentered psychological inputs. By default seed nuisance
#' effects are estimated conditional on centered task main effects, preserving
#' task effects while removing run baselines and supplied nuisance contributions.
#' Seed units are retained. Rank is recorded without removing aliased columns.
#'
#' @param seed Numeric raw BOLD vector on the complete acquisition grid.
#' @param psych A numeric matrix on the scan grid, or a list with `values`, a
#'   numeric matrix, and `grid`, either "scan" or "neural". Scan values are
#'   held constant across the microtime bins in each scan. Neural rows are
#'   concatenated by run, with `microtime` rows per scan. Column names identify
#'   conditions. Supply exactly one of `psych` and `events`.
#' @param events Data frame with onset and duration in seconds from the run
#'   start and condition columns; optional run and amplitude columns. Run is
#'   required for multiple runs. Duration events are overlap-weighted boxcars;
#'   zero-duration events have area equal to amplitude in their containing bin.
#' @param TR Positive sampling interval in seconds, common to all runs.
#' @param runs Vector of contiguous run labels, one per scan; default one run.
#' @param nuisance Optional finite numeric matrix with one row per scan.
#' @param hrf Numeric kernel sampled at `TR / microtime` from lag zero, or a
#'   function evaluated on that grid from zero to 32 seconds. The default is
#'   `fmrihrf::hrf_spmg1` evaluated on the same grid. Zeros are preserved.
#' @param deconv List with optional positive `lambda` (default 1) and integer
#'   `n_basis` (default: the number of scans in each run, one DCT component per
#'   scan). The ridge penalty is `lambda` times the mean diagonal of the
#'   observed basis Gram matrix, so the deconvolved seed and interaction terms
#'   do not depend on HRF amplitude or the microtime grid.
#' @param censor Logical vector marking excluded scans. The full time grid is
#'   retained and deconvolution uses only observed seed values.
#' @param run_effects Whether task, physiology and interaction terms are
#'   separate by run or shared; nuisance and intercepts are always run-specific.
#' @param center Center each psychological input on its run's neural grid before
#'   forming the interaction. Task main effects remain uncentered.
#' @param microtime Positive number of equal bins per TR (default 16).
#' @param sampling_offset Scan sampling time in seconds from each TR's start;
#'   must be a multiple of the microtime bin width and less than TR.
#' @param seed_adjust `"nuisance"` (default) removes run intercept and nuisance
#'   contributions estimated conditional on sampled task effects centered over
#'   observed scans. Remaining task/nuisance overlap is rejected as unidentified.
#'   `"none"` explicitly deconvolves the supplied seed without adjustment.
#' @return A `ppi_design` list with full-grid `X`, `keep`, `runs`, `TR`, term
#'   metadata, inputs, neural seed, deconvolution and provenance. Convolution is
#'   an unnormalized causal Riemann sum; scans sample at `sampling_offset` within each TR.
#'   Physiological values on excluded rows are zero placeholders, and must not
#'   be fitted. Rank metadata describes the retained rows.
#' @export
ppi_design <- function(seed, psych = NULL, events = NULL, TR = 1, runs = NULL,
                       nuisance = NULL, hrf = NULL, deconv = list(), censor = NULL,
                       run_effects = c("separate", "shared"), center = TRUE,
                       microtime = 16L, sampling_offset = 0,
                       seed_adjust = c("nuisance", "none")) {
  run_effects <- match.arg(run_effects)
  seed_adjust <- match.arg(seed_adjust)
  if (!is.numeric(seed) || !is.null(dim(seed)) || !length(seed)) stop("seed must be a numeric vector")
  n <- length(seed)
  if (length(TR) != 1L || !is.numeric(TR) || !is.finite(TR) || TR <= 0) stop("TR must be a positive scalar")
  if (length(microtime) != 1L || !is.numeric(microtime) || !is.finite(microtime) || microtime < 1 || microtime != floor(microtime)) stop("microtime must be a positive integer")
  microtime <- as.integer(microtime); dt <- TR / microtime
  if (!is.numeric(sampling_offset) || length(sampling_offset) != 1L || !is.finite(sampling_offset) || sampling_offset < 0 || sampling_offset >= TR || abs(sampling_offset / dt - round(sampling_offset / dt)) > 1e-8) stop("sampling_offset must be in [0, TR) and a multiple of TR/microtime")
  offset_bins <- as.integer(round(sampling_offset / dt))
  if (length(center) != 1L || is.na(center) || !is.logical(center)) stop("center must be TRUE or FALSE")
  if (is.null(runs)) runs <- rep(1L, n)
  if (length(runs) != n || anyNA(runs)) stop("runs must identify every scan")
  runs <- as.character(runs)
  if (any(!nzchar(runs)) || anyDuplicated(rle(runs)$values)) stop("runs must be nonempty contiguous blocks")
  labels <- unique(runs); blocks <- lapply(labels, function(x) which(runs == x))
  if (is.null(censor)) censor <- rep(FALSE, n)
  if (!is.logical(censor) || length(censor) != n || anyNA(censor)) stop("censor must be a logical vector marking excluded scans")
  keep <- !censor
  if (any(!is.finite(seed[keep]))) stop("observed seed values must be finite")
  if (any(vapply(blocks, function(i) !any(keep[i]), logical(1)))) stop("each run needs observed seed values")
  if (is.null(psych) == is.null(events)) stop("supply exactly one of psych or events")
  if (!is.null(psych)) {
    if (is.matrix(psych)) psych <- list(values = psych, grid = "scan")
    if (!is.list(psych) || !all(c("values", "grid") %in% names(psych)) || length(psych$grid) != 1L || !psych$grid %in% c("scan", "neural")) stop("psych must be list(values = matrix, grid = 'scan' or 'neural')")
    p <- psych$values
    if (!is.matrix(p) || !is.numeric(p) || !ncol(p) || any(!is.finite(p))) stop("psych values must be a finite numeric matrix")
    expected <- if (psych$grid == "scan") n else n * microtime
    if (nrow(p) != expected) stop("psych row count does not match its declared grid")
    if (is.null(colnames(p))) colnames(p) <- paste0("condition", seq_len(ncol(p)))
    if (anyNA(colnames(p)) || any(!nzchar(colnames(p))) || anyDuplicated(colnames(p))) stop("psych condition names must be unique and nonempty")
    neural_psych <- if (psych$grid == "scan") p[rep(seq_len(n), each = microtime), , drop = FALSE] else p
  } else {
    if (!is.data.frame(events) || !all(c("onset", "duration", "condition") %in% names(events))) stop("events must have onset, duration, and condition columns")
    if (!nrow(events)) stop("events must contain at least one condition")
    if (!is.numeric(events$onset) || !is.numeric(events$duration) || any(!is.finite(events$onset)) || any(!is.finite(events$duration)) || any(events$onset < 0) || any(events$duration < 0)) stop("event onset and duration must be finite and nonnegative")
    condition <- as.character(events$condition)
    if (anyNA(condition) || any(!nzchar(condition))) stop("event conditions must be nonempty")
    erun <- if ("run" %in% names(events)) as.character(events$run) else if (length(labels) == 1L) rep(labels, nrow(events)) else stop("events require run for multiple runs")
    if (anyNA(erun) || any(!erun %in% labels)) stop("event run must identify an existing run")
    amp <- if ("amplitude" %in% names(events)) events$amplitude else rep(1, nrow(events))
    if (!is.numeric(amp) || any(!is.finite(amp))) stop("event amplitude must be finite numeric")
    conditions <- unique(condition)
    neural_psych <- matrix(0, n * microtime, length(conditions), dimnames = list(NULL, conditions))
    for (j in seq_len(nrow(events))) {
      ri <- match(erun[j], labels); m <- length(blocks[[ri]]) * microtime
      start <- events$onset[j]; end <- start + events$duration[j]
      if (start >= m * dt || end > m * dt + 1e-10) stop("event extends outside its run")
      offset <- (blocks[[ri]][1L] - 1L) * microtime
      if (events$duration[j] == 0) {
        bin <- floor(start / dt) + 1L
        neural_psych[offset + bin, condition[j]] <- neural_psych[offset + bin, condition[j]] + amp[j] / dt
      } else {
        left <- (seq_len(m) - 1L) * dt
        overlap <- pmax(0, pmin(left + dt, end) - pmax(left, start))
        neural_psych[offset + seq_len(m), condition[j]] <- neural_psych[offset + seq_len(m), condition[j]] + amp[j] * overlap / dt
      }
    }
  }
  if (is.null(hrf)) hrf <- fmrihrf::hrf_spmg1(seq(0, 32, by = dt))
  else if (is.function(hrf)) hrf <- hrf(seq(0, 32, by = dt))
  if (!is.numeric(hrf) || !length(hrf) || any(!is.finite(hrf)) || !any(hrf != 0)) stop("hrf must be a finite nonzero numeric kernel or function")
  if (!is.list(deconv) || is.null(names(deconv)) && length(deconv) || any(!names(deconv) %in% c("lambda", "n_basis"))) stop("deconv accepts only lambda and n_basis")
  lambda <- if (is.null(deconv$lambda)) 1 else deconv$lambda
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) stop("deconv lambda must be positive")
  if (!is.null(deconv$n_basis) && (!is.numeric(deconv$n_basis) || length(deconv$n_basis) != 1L || !is.finite(deconv$n_basis) || deconv$n_basis < 1 || deconv$n_basis != floor(deconv$n_basis))) stop("n_basis must be a positive integer")
  if (is.null(nuisance)) nuisance <- matrix(numeric(), n, 0L)
  if (!is.matrix(nuisance) || !is.numeric(nuisance) || nrow(nuisance) != n || any(!is.finite(nuisance))) stop("nuisance must be a finite numeric matrix with one row per scan")
  if (ncol(nuisance) && is.null(colnames(nuisance))) colnames(nuisance) <- paste0("nuisance", seq_len(ncol(nuisance)))
  if (ncol(nuisance) && (anyNA(colnames(nuisance)) || any(!nzchar(colnames(nuisance))) || anyDuplicated(colnames(nuisance)))) stop("nuisance names must be unique and nonempty")
  # Convolve every column of x; cache the padded kernel transform at each run length.
  kernel_fft <- new.env(parent = emptyenv())
  conv <- function(x) {
    x <- as.matrix(x); len <- nrow(x)
    size <- stats::nextn(len + length(hrf) - 1L, factors = 2L)
    key <- as.character(size)
    if (!exists(key, envir = kernel_fft, inherits = FALSE))
      assign(key, stats::fft(c(hrf, rep(0, size - length(hrf)))), envir = kernel_fft)
    padded <- rbind(x, matrix(0, size - len, ncol(x)))
    transformed <- stats::mvfft(padded) * get(key, envir = kernel_fft)
    Re(stats::mvfft(transformed, inverse = TRUE))[seq_len(len), , drop = FALSE] * dt / size
  }
  # Orthonormal DCT-II columns for frequencies f on an m-bin grid.
  dct <- function(f, m) {
    b <- outer(seq_len(m) - 0.5, f, function(t, f) cos(pi * t * f / m))
    sweep(b, 2L, ifelse(f == 0, sqrt(1 / m), sqrt(2 / m)), "*")
  }
  task <- interaction <- matrix(0, n, ncol(neural_psych), dimnames = list(NULL, colnames(neural_psych)))
  neural_seed <- numeric(n * microtime); adjusted_seed <- seed
  info <- vector("list", length(blocks)); names(info) <- labels
  for (r in seq_along(blocks)) {
    i <- blocks[[r]]; mi <- seq.int((i[1L] - 1L) * microtime + 1L, i[length(i)] * microtime)
    m <- length(mi); sample <- seq.int(1L + offset_bins, m, by = microtime)
    task[i, ] <- conv(neural_psych[mi, , drop = FALSE])[sample, , drop = FALSE]
    obs <- keep[i]
    if (seed_adjust == "nuisance") {
      # Protect variation in task effects, with baseline overlap removed first.
      protected <- task[i[obs], , drop = FALSE]
      original_scale <- apply(abs(protected), 2L, max)
      protected <- sweep(protected, 2L, colMeans(protected), "-")
      negligible <- apply(abs(protected), 2L, max) <= 1e-10 * original_scale
      protected[, negligible] <- 0
      qt <- .ppi_design_basis(protected)
      qn <- .ppi_design_basis(cbind(1, nuisance[i[obs], , drop = FALSE]))
      rn <- qn - qt %*% crossprod(qt, qn)
      singular <- svd(rn, nu = 0L, nv = 0L)$d
      if (length(singular) < ncol(qn) || min(singular) <= 1e-7)
        stop("Seed nuisance adjustment is unidentified because nuisance overlaps protected task effects; revise the design or explicitly choose seed_adjust='none'.")
      adjusted_seed[i[obs]] <- seed[i[obs]] - drop(qn %*% qr.solve(rn, seed[i[obs]], tol = 1e-7))
    }
    k <- if (is.null(deconv$n_basis)) length(i) else deconv$n_basis
    if (k > m) stop("n_basis cannot exceed the microtime grid length of any run")
    # Build the basis in column blocks so memory stays O(m) per block, not O(m * k).
    freq_blocks <- split(seq_len(k) - 1L, (seq_len(k) - 1L) %/% 64L)
    A <- matrix(0, length(i), k)
    for (f in freq_blocks) A[, f + 1L] <- conv(dct(f, m))[sample, , drop = FALSE]
    ao <- A[obs, , drop = FALSE]
    gram <- crossprod(ao)
    # The ridge is relative to the mean basis energy, so the solution does not
    # depend on HRF amplitude or the microtime grid.
    ridge <- lambda * mean(diag(gram))
    if (!is.finite(ridge) || ridge <= 0) stop("deconvolution basis has no energy at the observed scans")
    coef <- solve(gram + diag(ridge, k), crossprod(ao, adjusted_seed[i][obs]))
    u <- numeric(m)
    for (f in freq_blocks) u <- u + drop(dct(f, m) %*% coef[f + 1L])
    neural_seed[mi] <- u
    pp <- neural_psych[mi, , drop = FALSE]
    if (center) pp <- sweep(pp, 2L, colMeans(pp), "-")
    interaction[i, ] <- conv(u * pp)[sample, , drop = FALSE]
    info[[r]] <- list(n_basis = k, observed = sum(obs), ridge = ridge,
                      condition = kappa(gram + diag(ridge, k)))
  }
  cols <- list(); term <- list()
  add <- function(x, name, type, condition = NA_character_, run = NA_character_) {
    cols[[length(cols) + 1L]] <<- x
    term[[length(term) + 1L]] <<- data.frame(name = name, type = type, condition = condition, run = run, stringsAsFactors = FALSE)
  }
  for (r in seq_along(blocks)) add(as.numeric(runs == labels[r]), paste0("intercept:", labels[r]), "intercept", run = labels[r])
  phys <- adjusted_seed; phys[!keep] <- 0
  for (type in c("nuisance", "task", "physiology", "interaction")) {
    mat <- switch(type, nuisance = nuisance, task = task, physiology = matrix(phys, ncol = 1L, dimnames = list(NULL, "seed")), interaction = interaction)
    if (!ncol(mat)) next
    for (j in seq_len(ncol(mat))) {
      cond <- if (type %in% c("task", "interaction")) colnames(mat)[j] else NA_character_
      if (run_effects == "shared" && type != "nuisance") add(mat[, j], paste0(type, ":", colnames(mat)[j]), type, cond)
      else for (r in seq_along(blocks)) {
        x <- numeric(n); x[blocks[[r]]] <- mat[blocks[[r]], j]
        add(x, paste0(type, ":", colnames(mat)[j], ":", labels[r]), type, cond, labels[r])
      }
    }
  }
  terms <- do.call(rbind, term); terms$name <- make.unique(terms$name)
  X <- do.call(cbind, cols); colnames(X) <- terms$name
  rank <- qr(X[keep, , drop = FALSE])$rank
  structure(list(X = X, runs = runs, TR = TR, keep = keep, terms = terms,
                 seed = seed, seed_adjusted = phys, psych = neural_psych, neural_seed = neural_seed,
                 deconv = list(method = "dct_ridge", lambda = lambda, runs = info),
                 provenance = list(microtime = microtime, dt = dt, hrf = hrf,
                                   sampling_offset = sampling_offset, convolution = "causal Riemann sum",
                                   psych_grid = if (is.null(events)) psych$grid else "events",
                                   center = center, run_effects = run_effects, seed_adjust = seed_adjust,
                                   rank = rank, n_columns = ncol(X))), class = "ppi_design")
}

.ppi_design_basis <- function(x) {
  if (!ncol(x)) return(matrix(0, nrow(x), 0L))
  scale <- apply(abs(x), 2L, max); scale[scale == 0] <- 1
  z <- sweep(x, 2L, scale, "/")
  s <- svd(z, nu = min(dim(z)), nv = 0L)
  rank <- if (length(s$d) && s$d[1L] > 0) sum(s$d > 1e-7 * s$d[1L]) else 0L
  s$u[, seq_len(rank), drop = FALSE]
}
