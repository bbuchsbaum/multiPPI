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
#' @return List with `S` (orthogonalised neural sticks), `names`, and `lambda`.
#' @export
mppi_psych_neural_from_X <- function(X, psych_idx, h, K = NULL,
                                     method = c("gcv", "fixed"), lambda = 1e-1,
                                     groups = NULL) {
  method <- match.arg(method)
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
    S_res[, ii] <- S[, ii] - Qs %*% qr.coef(qr(Qs), S[, ii])
  }
  list(S = S_res, names = colnames(S_res), lambda = dec$lambda)
}

#' Convert an existing bold-domain fit to neural domain
#' @export
mppi_neural_from_fit <- function(fit, X, psych_idx, h, K = NULL,
                                 method = c("gcv", "fixed"), lambda = 1e-1, groups = NULL) {
  method <- match.arg(method)
  if (is.null(groups)) {
    groups <- tryCatch(mppi_group_hrf_columns(X[, psych_idx, drop = FALSE]),
                       error = function(e) NULL)
  }
  ps <- mppi_psych_neural_from_X(X, psych_idx, h, K, method, lambda, groups)
  S <- ps$S
  if (!is.null(fit$basis)) {
    if (is.null(fit$Z)) stop("fit must store basis residuals in $Z for conversion.")
    U <- mppi_deconv_dct(fit$Z, h = h, K = K, method = method, lambda = lambda)$U
    out <- fit
    out$Delta <- lapply(seq_len(ncol(S)), function(ii) {
      denom <- sum(S[, ii]^2)
      if (denom < .Machine$double.eps) matrix(NA_real_, ncol(fit$Z), ncol(fit$Z)) else {
        Mk <- crossprod(U, S[, ii] * U) / denom
        diag(Mk) <- 0
        Mk
      }
    })
    out$names <- ps$names
    out$U <- U
    out$pk <- lapply(seq_len(ncol(S)), function(ii) S[, ii])
    out$domain <- "neural"
    out
  } else {
    if (is.null(fit$R)) stop("fit must store residuals in $R for conversion.")
    Ufull <- mppi_deconv_dct(fit$R, h = h, K = K, method = method, lambda = lambda)$U
    out <- fit
    out$Delta <- lapply(seq_len(ncol(S)), function(ii) {
      denom <- sum(S[, ii]^2)
      if (denom < .Machine$double.eps) matrix(NA_real_, ncol(fit$R), ncol(fit$R)) else {
        Dk <- crossprod(Ufull, S[, ii] * Ufull) / denom
        diag(Dk) <- 0
        Dk
      }
    })
    out$names <- ps$names
    out$U <- Ufull
    out$pk <- lapply(seq_len(ncol(S)), function(ii) S[, ii])
    out$domain <- "neural"
    out
  }
}

#' Compare bold vs neural-domain fits (Δ-gap, AIC difference)
#' @export
mppi_compare_models <- function(fit_bold, fit_neural, resid_bold = NULL, resid_neural = NULL, pk = NULL) {
  stopifnot(length(fit_bold$names) == length(fit_neural$names))
  gap <- numeric(length(fit_bold$names))
  for (i in seq_along(fit_bold$names)) {
    B <- mppi_get_M(fit_bold, i)
    N <- mppi_get_M(fit_neural, i)
    gap[i] <- .mppi_frob2(B - N) / (sqrt(.mppi_frob2(N)) + 1e-12)
  }
  out <- list(names = fit_bold$names, delta_gap = gap)
  if (!is.null(resid_bold) && !is.null(resid_neural) && !is.null(pk)) {
    rss_b <- rss_n <- numeric(length(fit_bold$names))
    for (i in seq_along(fit_bold$names)) {
      denom <- sum(pk[[i]]^2)
      Bb <- crossprod(resid_bold, pk[[i]] * resid_bold) / denom
      Bn <- crossprod(resid_neural, pk[[i]] * resid_neural) / denom
      rss_b[i] <- sum((Bb - mppi_get_M(fit_bold, i))^2, na.rm = TRUE)
      rss_n[i] <- sum((Bn - mppi_get_M(fit_neural, i))^2, na.rm = TRUE)
    }
    n <- nrow(resid_bold)
    AIC_b <- n * log(rss_b / n)
    AIC_n <- n * log(rss_n / n)
    out$AIC_bold <- AIC_b
    out$AIC_neural <- AIC_n
    out$AIC_diff <- AIC_b - AIC_n
  }
  out
}

#' Evidence-weighted HRF ensemble average
#' @export
mppi_hrf_ensemble <- function(fits, weights = NULL) {
  stopifnot(is.list(fits), length(fits) >= 2)
  fit_names <- names(fits)
  method_tag <- "user"
  aic_scores <- NULL
  aic_list <- lapply(fits, function(f) f$AIC_total)
  len <- vapply(aic_list, length, integer(1))
  have_aic <- length(len) && length(unique(len)) == 1L && len[1] > 0L
  aic_mat <- if (have_aic) do.call(cbind, aic_list) else NULL
  if (is.null(weights)) {
    method_tag <- "uniform"
    if (have_aic) {
      if (all(is.finite(aic_mat))) {
        aic_scores <- colSums(aic_mat)
        w <- exp(-0.5 * (aic_scores - min(aic_scores)))
        if (sum(w) > 0) {
          weights <- w / sum(w)
          method_tag <- "aic"
        }
      }
    }
    if (is.null(weights)) {
      weights <- rep(1 / length(fits), length(fits))
    }
  } else {
    weights <- as.numeric(weights)
  }
  if (anyNA(weights) || any(weights < 0)) {
    stop("weights must be non-negative and finite.", call. = FALSE)
  }
  if (sum(weights) == 0) stop("At least one weight must be positive.", call. = FALSE)
  weights <- weights / sum(weights)
  if (is.null(fit_names)) {
    fit_names <- paste0("fit", seq_along(fits))
  }
  names(weights) <- fit_names
  out <- fits[[1]]
  combine_matrix <- function(idx) {
    Mk <- Reduce(`+`, Map(function(w, f) w * mppi_get_M(f, idx), weights, fits))
    if (isTRUE(out$packed)) .mppi_pack_upper(Mk) else Mk
  }
  out$Delta <- lapply(seq_along(out$names), combine_matrix)
  if (have_aic && all(dim(aic_mat))) {
    if (is.null(aic_scores)) aic_scores <- colSums(aic_mat)
    out$AIC_total <- as.numeric(aic_mat %*% weights)
  }
  out$ensemble <- list(weights = weights, method = method_tag, scores = aic_scores)
  out
}

#' Gain/precision index from ΔTheta
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

#' Routing asymmetry index using lagged fits and a hierarchy vector
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
      mppi_get_M(fit, k)
    }
  }
  u <- as.numeric(hierarchy)
  if (length(u) != ncol(X0)) stop("Hierarchy vector length must match fit dimension.")
  Wdir <- tcrossprod(u, u)
  pos_lags <- intersect(lag_vals, pos)
  neg_lags <- intersect(lag_vals, neg)
  proj_energy <- function(lag) {
    Mk <- fetch_matrix(lag)
    sum(Mk * Wdir, na.rm = TRUE)
  }
  Epos <- if (length(pos_lags)) sum(vapply(pos_lags, proj_energy, numeric(1))) else 0
  Eneg <- if (length(neg_lags)) sum(vapply(neg_lags, proj_energy, numeric(1))) else 0
  R <- (Epos - Eneg) / (Epos + Eneg + eps)
  list(routing = R, energy_pos = Epos, energy_neg = Eneg,
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
#' @param lag_fits Optional named list of lagged fits for routing.
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
