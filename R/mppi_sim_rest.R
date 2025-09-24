# Synthetic rest-state generators ----------------------------------------------

#' Simulate low-rank neural time-series with known modulators
#'
#' Generates latent component activity that depends linearly on supplied
#' oscillator-like modulators. The simulator constructs modulators with zero mean
#' and unit variance and injects symmetric coupling patterns so that
#' `mppi_rest()` can recover them.
#'
#' @param T Number of time points.
#' @param r Latent dimensionality (components).
#' @param contexts Character vector of context names.
#' @param amp Vector of coupling amplitudes per context.
#' @param freq Vector of temporal frequencies (in cycles / TR) for modulators.
#' @param noise_sd Standard deviation of additive white noise on latent series.
#' @param seed Optional seed for reproducibility.
#' @param templates Optional list of r x r matrices (one per context) specifying
#'   the coupling templates. Each matrix is symmetrised and scaled to unit
#'   Frobenius norm. When omitted, random symmetric templates are generated.
#' @return List with latent series `U`, modulators, and true coupling matrices.
#' @export
mppi_sim_rest_neural <- function(T = 600L, r = 6L,
                                 contexts = c("ctx1", "ctx2"),
                                 amp = rep(0.4, length(contexts)),
                                 freq = seq(0.01, length.out = length(contexts), by = 0.005),
                                 noise_sd = 0.02,
                                 seed = NULL,
                                 templates = NULL) {
  stopifnot(T > 0, r > 0)
  if (length(contexts) != length(amp)) {
    stop("Length of 'amp' must equal number of contexts.", call. = FALSE)
  }
  if (length(freq) != length(contexts)) {
    stop("Length of 'freq' must equal number of contexts.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  K <- length(contexts)
  t_idx <- seq_len(T)
  # Build zero-mean unit-variance modulators (sine waves with random phases)
  phases <- runif(K, 0, 2 * pi)
  g_raw <- vapply(seq_len(K), function(k) {
    g <- sin(2 * pi * freq[k] * t_idx + phases[k])
    as.numeric(scale(g))
  }, numeric(T))
  colnames(g_raw) <- contexts
  modulators <- as.list(as.data.frame(g_raw))

  # Symmetric coupling templates (unit Frobenius)
  B_list <- vector("list", K)
  if (!is.null(templates)) {
    if (length(templates) != K) {
      stop("templates list must match number of contexts.", call. = FALSE)
    }
  }
  for (k in seq_len(K)) {
    if (!is.null(templates)) {
      Tk <- templates[[k]]
      if (!is.matrix(Tk) || any(dim(Tk) != c(r, r))) {
        stop("Each template must be an r x r matrix.", call. = FALSE)
      }
      S <- 0.5 * (Tk + t(Tk))
    } else {
      A <- matrix(rnorm(r * r), r, r)
      S <- 0.5 * (A + t(A))
    }
    fnorm <- sqrt(sum(S * S))
    if (fnorm < .Machine$double.eps) {
      S <- diag(r)
      fnorm <- sqrt(sum(S * S))
    }
    B_list[[k]] <- S / fnorm
  }
  names(B_list) <- contexts

  # Generate latent activity
  Z <- matrix(rnorm(T * r), T, r)
  U <- Z
  for (k in seq_len(K)) {
    U <- U + amp[k] * g_raw[, k] * (Z %*% B_list[[k]])
  }
  if (noise_sd > 0) {
    U <- U + matrix(rnorm(T * r, sd = noise_sd), T, r)
  }
  list(U = U,
       modulators = modulators,
       coupling = B_list,
       amp = amp,
       freq = freq)
}

#' Simulate observed resting-state dataset from latent components
#'
#' @param T Number of time points.
#' @param V Observed dimensionality (voxels/ROIs).
#' @param r Latent dimensionality.
#' @param hrf Optional HRF vector to convolve each latent component.
#' @param obs_noise Standard deviation of observation noise.
#' @param templates Optional list of coupling templates forwarded to
#'   [mppi_sim_rest_neural()].
#' @param ... Passed to [mppi_sim_rest_neural()].
#' @return List with observed matrix `Y`, latent `U`, modulators, loadings, and truth.
#' @export
mppi_sim_rest_dataset <- function(T = 600L, V = 40L, r = 6L,
                                  hrf = NULL, obs_noise = 0.05,
                                  templates = NULL, ...) {
  sim <- mppi_sim_rest_neural(T = T, r = r, templates = templates, ...)
  U <- sim$U
  if (!is.null(hrf)) {
    h <- as.numeric(hrf)
    if (!length(h)) stop("HRF must be a numeric vector.", call. = FALSE)
    conv_col <- function(x) stats::filter(c(x, rep(0, length(h))), h, sides = 1)[seq_len(length(x))]
    U <- apply(U, 2, conv_col)
  }
  L <- qr.Q(qr(matrix(rnorm(V * r), V, r)))
  Y <- U %*% t(L)
  if (obs_noise > 0) {
    Y <- Y + matrix(rnorm(T * V, sd = obs_noise), T, V)
  }
  mods <- sim$modulators
  truth <- list(coupling = sim$coupling,
                amp = sim$amp,
                freq = sim$freq,
                U = sim$U)
  list(Y = Y,
       U = sim$U,
       loadings = L,
       modulators = mods,
       truth = truth)
}

#' Simulate a cohort of resting datasets with trait-linked coupling
#'
#' @param n_subj Number of subjects.
#' @param trait Optional numeric vector of length `n_subj`.
#' @param trait_effect Numeric vector of per-context slopes mapping trait to
#'        coupling amplitude.
#' @param amp_base Baseline coupling amplitudes per context.
#' @param amp_noise SD of subject-level amplitude noise.
#' @param seed Optional seed for reproducibility.
#' @param templates Optional list of coupling templates applied to every subject.
#' @param ... Additional arguments passed to [mppi_sim_rest_dataset()] (shared across subjects).
#' @return List with `datasets`, `trait`, and subject-level amplitude table.
#' @export
mppi_sim_rest_cohort <- function(n_subj = 10L,
                                 trait = NULL,
                                 contexts = NULL,
                                 trait_effect = c(0.5, 0),
                                 amp_base = rep(0.3, length(trait_effect)),
                                 amp_noise = 0.05,
                                 seed = NULL,
                                 templates = NULL,
                                 ...) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(trait)) trait <- scale(rnorm(n_subj))[, 1]
  if (length(trait) != n_subj) stop("Trait length must equal n_subj.", call. = FALSE)
  if (is.null(contexts)) {
    contexts <- paste0("ctx", seq_along(trait_effect))
  }
  if (length(contexts) != length(trait_effect)) {
    stop("contexts and trait_effect must have the same length.", call. = FALSE)
  }
  if (length(amp_base) != length(trait_effect)) {
    stop("amp_base and trait_effect must have same length.", call. = FALSE)
  }
  amps <- matrix(NA_real_, n_subj, length(trait_effect), dimnames = list(NULL, contexts))
  datasets <- vector("list", n_subj)
  for (i in seq_len(n_subj)) {
    subj_amp <- amp_base + trait_effect * trait[i] + rnorm(length(trait_effect), sd = amp_noise)
    amps[i, ] <- subj_amp
    datasets[[i]] <- mppi_sim_rest_dataset(amp = subj_amp,
                                          contexts = contexts,
                                          templates = templates,
                                          ...)
  }
  list(datasets = datasets,
       trait = as.numeric(trait),
       amplitudes = amps)
}
