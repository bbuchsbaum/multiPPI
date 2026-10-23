# Neural-domain helpers ----------------------------------------------------

#' Default canonical HRF used for neural-domain fits
#'
#' Provides a double-gamma haemodynamic response similar to SPM's canonical
#' kernel so `mppi_fit(..., domain = "neural")` can operate without an
#' explicit HRF from the caller. Oversampling lets callers request smoother
#' kernels but the default returns one value per TR.
#'
#' @param tr Repetition time (seconds).
#' @param duration Total length of the HRF support (seconds).
#' @param oversample Optional oversampling factor relative to `tr`.
#' @return Numeric vector containing the canonical HRF sampled at `tr`.
#' @export
mppi_default_hrf <- function(tr = 1, duration = 32, oversample = 1) {
  stopifnot(tr > 0, duration > 0)
  oversample <- as.integer(oversample)
  if (oversample < 1) stop("oversample must be >= 1")
  dt <- tr / oversample
  t <- seq(0, duration, by = dt)
  a1 <- 6; a2 <- 16
  b1 <- 1; b2 <- 1
  c <- 1 / 6
  h <- stats::dgamma(t, shape = a1, scale = b1) - c * stats::dgamma(t, shape = a2, scale = b2)
  if (sum(h) != 0) h <- h / sum(h)
  idx <- seq(1, length(h), by = oversample)
  h[idx]
}

# Internal helpers ------------------------------------------------------

.mppi_match_regressor <- function(fit, k) {
  idx <- if (is.character(k)) match(k, fit$names) else as.integer(k)
  if (is.na(idx) || idx < 1 || idx > length(fit$pk)) {
    stop("Invalid regressor index", call. = FALSE)
  }
  idx
}

.mppi_residual_matrix <- function(fit) {
  if (!is.null(fit$basis) && !is.null(fit$Z)) return(fit$Z)
  if (!is.null(fit$R)) return(fit$R)
  if (!is.null(fit$U)) return(fit$U)
  stop("Fit must expose residuals in $R, $Z, or $U for mechanistic helpers.", call. = FALSE)
}

.mppi_full_slope <- function(fit, idx) {
  pk <- fit$pk[[idx]]
  if (is.null(pk)) stop("Fit does not store pk for requested regressor.", call. = FALSE)
  X <- .mppi_residual_matrix(fit)
  if (length(pk) != nrow(X)) stop("pk length does not match residual rows.", call. = FALSE)
  denom <- sum(pk^2)
  if (denom < .Machine$double.eps) {
    return(matrix(NA_real_, ncol(X), ncol(X)))
  }
  crossprod(X, pk * X) / denom
}

#' Deconvolve psychological regressors into neural-domain "sticks"
#' 
#' Converts the psychological portions of the design matrix to neural-space
#' predictors using the same DCT deconvolution used inside `mppi_fit`. When
#' groups are supplied the resulting sticks are grouped/averaged by condition
#' (useful for HRF basis expansions).
#'
#' @param X Full design matrix used in the MPPI fit.
#' @param psych_idx Integer indices of the psychological columns in `X`.
#' @param h Haemodynamic response function vector.
#' @param K Number of DCT components to retain (defaults to 64 or `T/2`).
#' @param method Deconvolution method (`"gcv"` or `"fixed"`).
#' @param lambda Ridge parameter when `method = "fixed"`.
#' @param groups Optional list grouping columns (relative to `psych_idx`).
#' @param runs Optional run labels; deconvolution is independent within each run.
#' @return List with `S` (orthogonalised neural sticks), `names`, and `lambda`.
#' @export
mppi_psych_neural_from_X <- function(X, psych_idx, h, K = NULL,
                                     method = c("gcv", "fixed"), lambda = 1e-1,
                                     groups = NULL, runs = NULL) {
  method <- match.arg(method)
  if (any(!is.finite(X))) stop("Censored or nonfinite design is unsupported for neural deconvolution.", call. = FALSE)
  if (!is.null(runs)) {
    if (length(runs) != nrow(X) || anyNA(runs)) stop("runs must have one nonmissing label per row.", call. = FALSE)
    segments <- split(seq_len(nrow(X)), cumsum(c(TRUE, runs[-1L] != runs[-length(runs)])))
    if (length(segments) > 1L) {
      if (is.null(K)) K <- min(64L, floor(min(lengths(segments))/2L))
      parts <- lapply(segments, function(i) mppi_psych_neural_from_X(X[i,,drop=FALSE], psych_idx, h, K, method, lambda, groups))
      return(list(S=do.call(rbind,lapply(parts, `[[`, "S")), names=parts[[1]]$names, lambda=lapply(parts, `[[`, "lambda")))
    }
  }
  Xp <- X[, psych_idx, drop = FALSE]
  if (is.null(K)) K <- min(64L, floor(nrow(Xp) / 2L))
  dec <- mppi_deconv_dct(Xp, h = as.numeric(h), K = K, method = method, lambda = lambda)
  U <- dec$U
  if (is.null(groups)) {
    S <- U
    colnames(S) <- if (!is.null(colnames(Xp))) paste0("s:", colnames(Xp)) else paste0("s", seq_len(ncol(U)))
  } else {
    Slist <- vector("list", length(groups))
    nms   <- character(length(groups))
    for (g in seq_along(groups)) {
      idx <- groups[[g]]
      idx <- idx[idx >= 1 & idx <= ncol(U)]
      if (!length(idx)) stop("Empty group in 'groups'.")
      Ug <- U[, idx, drop = FALSE]
      norms <- sqrt(colSums(Ug^2)); norms[norms == 0] <- 1
      Ug <- sweep(Ug, 2, norms, "/")
      Slist[[g]] <- rowMeans(Ug)
      if (!is.null(colnames(Xp))) {
        nms[g] <- paste0("s:", paste(colnames(Xp)[idx], collapse = "+"))
      } else {
        nms[g] <- paste0("s:grp", g)
      }
    }
    S <- do.call(cbind, Slist)
    colnames(S) <- nms
  }
  S_res <- S
  for (ii in seq_len(ncol(S))) {
    Qs <- cbind(1, S[, setdiff(seq_len(ncol(S)), ii), drop = FALSE])
    S_res[, ii] <- .mppi_residualize_vec(S[, ii], Qs)
  }
  list(S = S_res, names = colnames(S_res), lambda = dec$lambda)
}

#' Convert an existing bold-domain fit to neural domain
#' @param fit result of `mppi_fit()` fit in the BOLD domain
#' @param X design matrix used for the original fit (time x predictors)
#' @param psych_idx indices of psychological regressors within `X`
#' @param h canonical HRF or vector used for deconvolution
#' @param K optional number of DCT components to retain during deconvolution
#' @param method smoothing strategy for HRF estimation (`"gcv"` or `"fixed"`)
#' @param lambda penalty strength when `method = "fixed"`
#' @param groups optional grouping of psychological columns for shared HRFs
#' @export
mppi_neural_from_fit <- function(fit, X, psych_idx, h, K = NULL,
                                 method = c("gcv", "fixed"), lambda = 1e-1, groups = NULL) {
  method <- match.arg(method)
  if (!inherits(fit, "mppi_fit") || !identical(fit$domain, "bold")) stop("Conversion requires a BOLD-domain mppi_fit.", call. = FALSE)
  if ((fit$dropped %||% 0L) > 0L || any(!is.finite(X))) stop("Censored observations are unsupported for neural conversion.", call. = FALSE)
  if (is.null(fit$R_raw) || nrow(X) != nrow(fit$R_raw)) stop("Original-grid BOLD residuals and matching design are required.", call. = FALSE)
  if (is.null(fit$zero_diag) || is.null(fit$fitted_X) || is.null(fit$psych_idx)) stop("Original design and zero_diag metadata are required for neural conversion; refit with current multiPPI.", call. = FALSE)
  supplied <- .mppi_prepare_design_matrix(X, psych_idx)
  if (!isTRUE(all.equal(supplied$X, fit$fitted_X, tolerance=1e-12)) || !identical(supplied$psych_idx, fit$psych_idx)) stop("Supplied design or psychological indices do not match the original fit.", call. = FALSE)
  if (!is.null(fit$ensemble)) stop("Ensemble slopes have no single residual series to convert.", call. = FALSE)
  cfg <- list(hrf=h, K=K, method=method, lambda=lambda, groups=groups)
  converted <- mppi_fit(fit$R_raw, X, psych_idx, runs=fit$runs,
                        zero_diag=fit$zero_diag, scale=fit$scale,
                        center_by=fit$center_by, backend=fit$backend,
                        packed=fit$packed, chunk_size=fit$chunk_size,
                        domain="neural", deconv=cfg, lags=fit$lags %||% 0L,
                        lag_blocklens=fit$lag_blocklens)
  converted$bold_residuals <- fit$R_raw
  converted$evidence <- fit$evidence
  converted$AIC_total <- if (length(fit$AIC_total)) setNames(rep(fit$AIC_total[[1L]], length(converted$names)), converted$names) else NULL
  converted$BIC_total <- if (length(fit$BIC_total)) setNames(rep(fit$BIC_total[[1L]], length(converted$names)), converted$names) else NULL
  converted$basis <- fit$basis
  if (!is.null(fit$basis)) converted$Z <- converted$R
  converted
}

#' Compare bold vs neural-domain fits by descriptive slope gap
#' @param fit_bold output of `mppi_fit()` in the BOLD domain
#' @param fit_neural output of `mppi_fit()` or `mppi_neural_from_fit()` in the neural domain
#' @param resid_bold optional residual matrix from the BOLD fit (time x V)
#' @param resid_neural optional residual matrix from the neural fit (time x V)
#' @param pk optional list of psychological regressors (aligned with fits) used for information criteria
#' @export
mppi_compare_models <- function(fit_bold, fit_neural, resid_bold = NULL, resid_neural = NULL, pk = NULL) {
  stopifnot(length(fit_bold$names) == length(fit_neural$names))
  gap <- numeric(length(fit_bold$names))
  for (i in seq_along(fit_bold$names)) {
    B <- mppi_get_M(fit_bold, i)
    N <- mppi_get_M(fit_neural, i)
    gap[i] <- .mppi_frob2(B - N) / (sqrt(.mppi_frob2(N)) + 1e-12)
  }
  out <- list(names = fit_bold$names, delta_gap = gap, estimand = "descriptive_cross_domain_slope_gap")
  if (!is.null(resid_bold) || !is.null(resid_neural) || !is.null(pk))
    stop("AIC comparisons across BOLD and neural response domains are unsupported.", call. = FALSE)
  out
}

#' HRF ensemble average with equal or caller-specified weights
#' @param fits list of neural-domain fits to combine
#' @param weights Optional non-negative vector with one weight per fit.
#' @export
mppi_hrf_ensemble <- function(fits, weights = NULL) {
  if (!is.list(fits) || length(fits) < 2L || !all(vapply(fits, inherits, logical(1), "mppi_fit"))) stop("Provide at least two mppi_fit objects.", call. = FALSE)
  ref <- fits[[1]]
  if (!identical(ref$domain, "neural")) stop("Ensemble requires neural-domain fits.", call. = FALSE)
  if (any(vapply(fits, function(f) any((f$lags %||% 0L) != 0L), logical(1)))) stop("Multi-lag HRF ensembles are unsupported.", call. = FALSE)
  compatible <- vapply(fits, function(f) identical(f$domain,ref$domain) && identical(f$names,ref$names) && identical(f$scale,ref$scale) && identical(f$zero_diag,ref$zero_diag) && identical(f$basis, ref$basis) && all(vapply(seq_along(ref$names), function(i) identical(dim(mppi_get_M(f,i)),dim(mppi_get_M(ref,i))), logical(1))), logical(1))
  if (!all(compatible)) stop("Ensemble fits must have matching domain, names, dimensions, scale and basis.", call. = FALSE)
  method_tag <- if (is.null(weights)) "uniform" else "user"
  if (is.null(weights)) weights <- rep(1, length(fits))
  if (!is.numeric(weights) || length(weights) != length(fits) || any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) stop("weights must contain one finite non-negative value per fit, with a positive sum.", call. = FALSE)
  if (!is.null(names(weights))) {
    if (is.null(names(fits)) || anyDuplicated(names(fits)) || anyDuplicated(names(weights)) || !setequal(names(weights), names(fits))) stop("Named weights must match unique fit names.", call. = FALSE)
    weights <- weights[names(fits)]
  }
  weights <- weights / sum(weights)
  names(weights) <- names(fits) %||% paste0("fit", seq_along(fits))
  out <- ref
  blend <- function(i, mode) Reduce(`+`, Map(function(w,f) w * mppi_get_M_scaled(f,i,mode=mode), weights, fits))
  out$Delta_raw <- lapply(seq_along(ref$names), blend, mode="raw")
  out$Delta_amplitude <- lapply(seq_along(ref$names), blend, mode="amplitude")
  out$Delta_normalized <- lapply(seq_along(ref$names), blend, mode="normalized")
  out$Delta <- lapply(seq_along(ref$names), function(i) {
    m <- switch(ref$scale, cov=out$Delta_raw[[i]], normalized=out$Delta_normalized[[i]], corr=out$Delta_normalized[[i]])
    attr(m,"raw") <- out$Delta_raw[[i]]; attr(m,"amplitude") <- out$Delta_amplitude[[i]]; attr(m,"normalized") <- out$Delta_normalized[[i]]
    if (isTRUE(ref$packed)) .mppi_pack_upper(m) else m
  })
  names(out$Delta) <- names(out$Delta_raw) <- names(out$Delta_amplitude) <- names(out$Delta_normalized) <- ref$names
  out$lagged <- setNames(lapply(out$Delta, function(m) setNames(list(m), "0")), ref$names)
  unavailable <- c("R","R_raw","U","Z","pk","sigma","Sigma0","Theta0","variance","partial","denom","evidence","AIC_total","BIC_total","covariance_diagonal","fitted_X","psych_idx","deconv")
  out[unavailable] <- rep(list(NULL), length(unavailable))
  out$ensemble <- list(weights=weights, method=method_tag, estimand="weighted_covariance_modulation_slopes")
  out
}

#' Gain/precision index from ΔTheta
#' @param fit result of `mppi_fit()` (basis aware fits supported)
#' @param k condition index or name within `fit`
#' @param lambda ridge penalty applied to the baseline precision
#' @param roi_idx optional subset of ROI indices when mapping a basis fit back to voxel/ROI space
#' @export
mppi_gain <- function(fit, k = 1L, lambda = 1e-3, roi_idx = NULL) {
  idx <- .mppi_match_regressor(fit, k)
  X <- .mppi_residual_matrix(fit)
  S0 <- crossprod(X) / nrow(X)
  Mk <- .mppi_full_slope(fit, idx)
  if (all(is.na(Mk))) {
    comp_gain <- rep(NA_real_, ncol(S0))
    dTheta <- matrix(NA_real_, ncol(S0), ncol(S0))
  } else {
    Theta0 <- solve(S0 + diag(lambda, ncol(S0)))
    dTheta <- - Theta0 %*% Mk %*% Theta0
    comp_gain <- diag(dTheta)
  }
  if (is.null(fit$basis)) {
    if (!is.null(roi_idx)) comp_gain <- comp_gain[roi_idx]
    list(gain = comp_gain, space = "roi")
  } else {
    if (!is.null(roi_idx)) {
      V <- fit$basis$V[roi_idx, , drop = FALSE]
      roi_gain <- rowSums((V %*% dTheta) * V)
      list(gain_component = comp_gain, gain_roi = roi_gain, idx = roi_idx, space = "roi+basis")
    } else {
      list(gain_component = comp_gain, space = "basis")
    }
  }
}

# Lead-lag energy E[i, j]: squared lagged slopes with i leading j, pooled over
# nonzero lags. Signal-lag matrices have rows leading at tau > 0 and satisfy
# M(-tau) = t(M(tau)), so negative lags are transposed before pooling.
.mppi_lead_energy <- function(mats, lags) {
  E <- NULL
  for (ii in seq_along(lags)) {
    lg <- as.integer(lags[ii]); M <- mats[[ii]]
    if (lg == 0L || is.null(M)) next
    M <- matrix(M, nrow(M), ncol(M))
    if (lg < 0L) M <- t(M)
    M[!is.finite(M)] <- 0
    E <- if (is.null(E)) M^2 else E + M^2
  }
  if (!is.null(E)) diag(E) <- 0
  E
}

# Directedness of lagged coupling in [0, 1]: sum over pairs of |E_ij - E_ji|
# divided by sum over pairs of (E_ij + E_ji). 0 = every pair symmetric, 1 = one-way.
.mppi_directedness <- function(E) {
  if (is.null(E)) return(NA_real_)
  sum(abs(E - t(E))) / (2 * sum(E) + 1e-12)
}

#' Routing asymmetry index using lagged fits and a hierarchy vector
#'
#' Measures whether task-modulated lagged coupling flows up or down a
#' hierarchy. For signal lags \eqn{\tau > 0}, \eqn{M(\tau)_{ij}} couples node
#' \eqn{i} at time \eqn{t} with node \eqn{j} at \eqn{t + \tau} (i leads j), and
#' \eqn{M(-\tau) = M(\tau)^T}; negative-lag matrices are therefore transposed
#' and pooled with the positive ones into the lead-lag energy
#' \eqn{E_{ij} = \sum_\tau M(\tau)_{ij}^2} (diagonal excluded). With hierarchy
#' \eqn{u}, upward energy is \eqn{\sum_{ij} E_{ij} (u_j - u_i)_+}, downward energy
#' is \eqn{\sum_{ij} E_{ij} (u_i - u_j)_+}, and the index is
#' `(up - down) / (up + down)`, in \eqn{[-1, 1]}: positive when lagged coupling
#' runs from low to high \eqn{u}, zero for symmetric (undirected) coupling.
#' @param fits_by_lag either an `mppi_fit` fitted with lagged matrices or a named
#'   list of fits (names coercible to the lag), each storing that lag in `$lagged`
#' @param k condition index or name to evaluate
#' @param hierarchy numeric hierarchy vector (length = number of components/ROIs)
#' @param space coordinate space to evaluate (`"auto"`, `"basis"`, or `"roi"`)
#' @param pos positive integer lags pooled into the lead-lag energy
#' @param neg negative integer lags pooled into the lead-lag energy (transposed)
#' @param eps numerical guard to avoid division by zero
#' @return List with `routing`, `energy_up`, `energy_down`, the pairwise
#'   `lead_energy` matrix, `pos`, `neg` and `lag_values`.
#' @export
mppi_routing_index <- function(fits_by_lag, k = 1L, hierarchy,
                               space = c("auto", "basis", "roi"),
                               pos = c(1, 2), neg = c(-1, -2), eps = 1e-8) {
  space <- match.arg(space)
  use_single_fit <- inherits(fits_by_lag, "mppi_fit")
  if (use_single_fit) {
    fit_ref <- fits_by_lag
    if (is.null(fit_ref$lags) || is.null(fit_ref$lagged)) {
      stop("Fit does not contain lagged outputs; refit with lags specified.", call. = FALSE)
    }
    lag_vals <- as.integer(fit_ref$lags)
    lag_vals <- sort(unique(lag_vals))
    if (space == "auto") space <- if (!is.null(fit_ref$basis)) "basis" else "roi"
    X0 <- .mppi_residual_matrix(fit_ref)
    fetch_matrix <- function(lag) mppi_get_M_lag(fit_ref, k, lag)
  } else {
    if (!length(fits_by_lag)) stop("Provide at least one lagged fit.")
    lag_names <- names(fits_by_lag)
    if (is.null(lag_names)) {
      stop("fits_by_lag must be a named list of lagged fits.", call. = FALSE)
    }
    lag_vals <- suppressWarnings(as.integer(lag_names))
    if (any(is.na(lag_vals))) {
      stop("Lag names must coerce to integers (e.g., '-2', '1').", call. = FALSE)
    }
    lag_vals <- sort(unique(lag_vals))
    fit_ref <- fits_by_lag[[1]]
    if (space == "auto") space <- if (!is.null(fit_ref$basis)) "basis" else "roi"
    X0 <- .mppi_residual_matrix(fit_ref)
    fit_lookup <- fits_by_lag
    fetch_matrix <- function(lag) {
      fit <- fit_lookup[[as.character(lag)]]
      if (is.null(fit)) {
        stop(sprintf("No fit supplied for lag %s", lag), call. = FALSE)
      }
      # contemporaneous slopes are symmetric and carry no direction
      mppi_get_M_lag(fit, k, lag)
    }
  }
  u <- as.numeric(hierarchy)
  if (length(u) != ncol(X0)) stop("Hierarchy vector length must match fit dimension.")
  use <- c(intersect(lag_vals, pos[pos > 0]), intersect(lag_vals, neg[neg < 0]))
  E <- .mppi_lead_energy(lapply(use, fetch_matrix), use)
  if (is.null(E)) E <- matrix(0, length(u), length(u))
  W <- outer(u, u, function(a, b) b - a)
  Eup <- sum(E * pmax(W, 0))
  Edown <- sum(E * pmax(-W, 0))
  R <- (Eup - Edown) / (Eup + Edown + eps)
  list(routing = R, energy_up = Eup, energy_down = Edown, lead_energy = E,
       pos = pos, neg = neg, lag_values = lag_vals)
}

#' Project an interaction matrix onto hypothesis templates
#'
#' @param fit An `mppi_fit` object (ROI or basis).
#' @param k Regressor index or name.
#' @param templates Named list of template matrices matching the fit dimension.
#' @param normalize Logical; divide by the Frobenius norm of each template.
#' @return Named numeric vector of template weights.
#' @export
mppi_project_templates <- function(fit, k = 1L, templates, normalize = TRUE) {
  if (!length(templates)) stop("templates must be a non-empty list.")
  idx <- .mppi_match_regressor(fit, k)
  Mk <- .mppi_full_slope(fit, idx)
  dim_fit <- ncol(Mk)
  apply_template <- function(W) {
    W <- as.matrix(W)
    if (!all(dim(W) == c(dim_fit, dim_fit))) stop("Template dimensions must match fit dimension.")
    num <- sum(Mk * W, na.rm = TRUE)
    if (!normalize) return(num)
    denom <- sqrt(sum(W^2, na.rm = TRUE))
    if (denom < .Machine$double.eps) return(NA_real_)
    num / denom
  }
  vapply(templates, apply_template, numeric(1))
}

#' Construct within/between-network templates
#'
#' @param labels Factor or character vector labelling each component.
#' @return Named list with `within` and `between` templates (unit Frobenius norm).
#' @export
mppi_templates_within_between <- function(labels) {
  labs <- as.character(labels)
  n <- length(labs)
  same <- outer(labs, labs, `==`)
  diag_mask <- diag(1, n)
  within <- same & !diag_mask
  between <- (!same) & !diag_mask
  normalize <- function(M) {
    denom <- sqrt(sum(M^2))
    if (denom < .Machine$double.eps) matrix(0, n, n) else M / denom
  }
  list(within = normalize(within), between = normalize(between))
}

#' Construct gradient-aligned and orthogonal templates
#'
#' @param u Numeric hierarchy/gradient vector.
#' @return List with `dir` (rank-1 projector) and `ortho` components.
#' @export
mppi_templates_gradient <- function(u) {
  u <- as.numeric(u)
  if (anyNA(u)) stop("Gradient vector contains NA values.")
  if (all(u == 0)) stop("Gradient vector must have non-zero length.")
  u <- u / sqrt(sum(u^2))
  W_dir <- tcrossprod(u, u)
  W_ortho <- diag(length(u)) - W_dir
  normalize <- function(M) {
    denom <- sqrt(sum(M^2))
    if (denom < .Machine$double.eps) matrix(0, nrow(M), ncol(M)) else M / denom
  }
  list(dir = normalize(W_dir), ortho = normalize(W_ortho))
}

#' Compile mechanistic gain/routing/template summaries
#'
#' @param fit `mppi_fit` object for the target regressor.
#' @param k Regressor index or name.
#' @param hierarchy Optional vector for routing index.
#' @param lag_fits Optional lagged `mppi_fit` (or named list of lagged fits) for
#'   [mppi_routing_index()].
#' @param templates Optional list of template matrices for projections.
#' @param lambda Ridge term for gain computation.
#' @param roi_idx Optional ROI indices for gain back-projection.
#' @param normalize_templates Logical; pass to `mppi_project_templates`.
#' @return List with gain, routing, and template summaries (any can be `NULL`).
#' @export
mppi_mechanism_report <- function(fit, k = 1L, hierarchy = NULL, lag_fits = NULL,
                                  templates = NULL, lambda = 1e-3, roi_idx = NULL,
                                  normalize_templates = TRUE) {
  gain <- mppi_gain(fit, k = k, lambda = lambda, roi_idx = roi_idx)
  routing <- if (!is.null(hierarchy) && !is.null(lag_fits)) {
    mppi_routing_index(lag_fits, k = k, hierarchy = hierarchy)
  } else NULL
  template_weights <- if (!is.null(templates)) {
    mppi_project_templates(fit, k = k, templates = templates, normalize = normalize_templates)
  } else NULL
  list(gain = gain, routing = routing, templates = template_weights)
}
