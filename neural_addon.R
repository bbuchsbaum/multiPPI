
# multiPPI neural-domain add-on -----------------------------------------------
# Depends: mppiDeconv (Rcpp), your existing multiPPI outputs

#' Deconvolve psychological columns (design) to neural "stick" vectors
#' @param X T x q design matrix (same one used in multiPPI fit)
#' @param psych_idx integer vector of psychological columns in X
#' @param h HRF kernel vector
#' @param K DCT components
#' @param method "gcv" or "fixed"
#' @param lambda ridge parameter if method="fixed"
#' @param groups optional list of integer vectors, grouping HRF-basis columns by condition
#' @return list with S (T x Kp), names
mppi_psych_neural_from_X <- function(X, psych_idx, h, K = NULL,
                                     method = c("gcv","fixed"), lambda = 1e-1,
                                     groups = NULL) {
  method <- match.arg(method)
  Xp <- X[, psych_idx, drop = FALSE]
  if (is.null(K)) K <- min(64L, floor(nrow(Xp)/2L))
  U <- mppiDeconv::mppi_deconv_dct(Xp, h = as.numeric(h), K = K, method = method, lambda = lambda)$U
  if (is.null(groups)) {
    S <- U
    colnames(S) <- if (!is.null(colnames(Xp))) paste0("s:", colnames(Xp)) else paste0("s", seq_len(ncol(U)))
  } else {
    Slist <- vector("list", length(groups))
    nms   <- character(length(groups))
    for (g in seq_along(groups)) {
      idx <- groups[[g]]
      G   <- intersect(idx, seq_len(ncol(U)))
      if (!length(G)) stop("Empty group in 'groups'.")
      Ug <- U[, G, drop = FALSE]
      nn <- sqrt(colSums(Ug^2)); nn[nn == 0] <- 1
      Ug <- sweep(Ug, 2, nn, "/")
      Slist[[g]] <- rowMeans(Ug)
      nms[g] <- paste0("s:", if (!is.null(colnames(Xp))) paste(colnames(Xp)[G], collapse = "+") else paste0("grp", g))
    }
    S <- do.call(cbind, Slist); colnames(S) <- nms
  }
  # residualize each s_k against the others
  S_res <- S
  for (ii in seq_len(ncol(S))) {
    Qs <- cbind(1, S[, setdiff(seq_len(ncol(S)), ii), drop = FALSE])
    S_res[, ii] <- S[, ii] - Qs %*% qr.coef(qr(Qs), S[, ii])
  }
  list(S = S_res, names = colnames(S_res))
}

#' Produce neural-domain M_k from an existing multiPPI fit (basis or ROI)
#' @param fit an mppi_fit (from your multiPPI package), either ROI-space or basis-space
#' @param X the design matrix used in that fit
#' @param psych_idx integer vector of psychological columns in X
#' @param h HRF kernel (vector in TR units)
#' @param K DCT components (default: min(64, T/2))
#' @param method "gcv" or "fixed"
#' @param lambda ridge if method="fixed"
#' @param groups optional list of indices (relative to psych_idx) to group HRF-basis columns by condition
#' @return an object like fit, with $M (basis) or $Delta (ROI) computed in neural domain,
#'         plus $U (neural residuals/components), $pk (the neural s_k)
mppi_neural_from_fit <- function(fit, X, psych_idx, h, K = NULL,
                                 method = c("gcv","fixed"), lambda = 1e-1, groups = NULL) {
  stopifnot(is.matrix(X), length(psych_idx) >= 1)
  method <- match.arg(method)
  # 1) neural psych vectors from design
  ps <- mppi_psych_neural_from_X(X, psych_idx, h, K, method, lambda, groups)
  S <- ps$S; pNms <- ps$names

  # 2) deconvolve residual BOLD to neural residuals/components
  if (!is.null(fit$basis)) {
    # basis-space
    stopifnot(!is.null(fit$Z))
    U <- mppiDeconv::mppi_deconv_dct(fit$Z, h = as.numeric(h), K = K, method = method, lambda = lambda)$U
    outM <- vector("list", ncol(S))
    for (ii in seq_len(ncol(S))) {
      denom <- sum(S[, ii]^2)
      if (denom < .Machine$double.eps) outM[[ii]] <- matrix(NA_real_, ncol(fit$Z), ncol(fit$Z))
      else {
        Mk <- crossprod(U, S[, ii] * U) / denom
        diag(Mk) <- 0
        outM[[ii]] <- Mk
      }
    }
    structure(list(M = outM, names = pNms, Z = fit$Z, U = U, pk = lapply(seq_len(ncol(S)), function(i) S[, i]),
                   basis = fit$basis, scale = fit$scale, domain = "neural"),
              class = c("mppi_fit","mppi_fit_basis"))
  } else {
    # ROI-space
    stopifnot(!is.null(fit$R))
    Ufull <- mppiDeconv::mppi_deconv_dct(fit$R, h = as.numeric(h), K = K, method = method, lambda = lambda)$U
    outD <- vector("list", ncol(S))
    for (ii in seq_len(ncol(S))) {
      denom <- sum(S[, ii]^2)
      if (denom < .Machine$double.eps) outD[[ii]] <- matrix(NA_real_, ncol(fit$R), ncol(fit$R))
      else {
        Dk <- crossprod(Ufull, S[, ii] * Ufull) / denom
        diag(Dk) <- 0
        outD[[ii]] <- Dk
      }
    }
    structure(list(Delta = outD, names = pNms, R = fit$R, U = Ufull, pk = lapply(seq_len(ncol(S)), function(i) S[, i]),
                   scale = fit$scale, domain = "neural"),
              class = "mppi_fit")
  }
}

#' Compare bold vs neural fits (Δ-gap and optional AIC)
mppi_compare_models <- function(fit_bold, fit_neural, resid_bold = NULL, resid_neural = NULL, pk = NULL) {
  nms <- fit_bold$names
  stopifnot(all(nms == fit_neural$names))
  gap <- numeric(length(nms))
  for (i in seq_along(nms)) {
    B <- if (!is.null(fit_bold$basis)) fit_bold$M[[i]] else fit_bold$Delta[[i]]
    N <- if (!is.null(fit_neural$basis)) fit_neural$M[[i]] else fit_neural$Delta[[i]]
    gap[i] <- sqrt(sum((B - N)^2, na.rm = TRUE)) / (sqrt(sum(N^2, na.rm = TRUE)) + 1e-12)
  }
  out <- list(delta_gap = gap, names = nms)
  if (!is.null(resid_bold) && !is.null(resid_neural) && !is.null(pk)) {
    rss_b <- rss_n <- numeric(length(nms))
    for (i in seq_along(nms)) {
      denom <- sum(pk[[i]]^2)
      Db <- crossprod(resid_bold, pk[[i]] * resid_bold) / denom
      Dn <- crossprod(resid_neural, pk[[i]] * resid_neural) / denom
      B  <- if (!is.null(fit_bold$basis)) fit_bold$M[[i]] else fit_bold$Delta[[i]]
      N  <- if (!is.null(fit_neural$basis)) fit_neural$M[[i]] else fit_neural$Delta[[i]]
      rss_b[i] <- sum((Db - B)^2, na.rm = TRUE)
      rss_n[i] <- sum((Dn - N)^2, na.rm = TRUE)
    }
    k <- 0; n <- length(resid_bold)
    AIC_b <- 2*k + n*log(rss_b/n)
    AIC_n <- 2*k + n*log(rss_n/n)
    out$AIC_bold   <- AIC_b
    out$AIC_neural <- AIC_n
    out$AIC_diff   <- AIC_b - AIC_n
  }
  out
}

#' Evidence-weighted HRF ensemble average of fits
mppi_hrf_ensemble <- function(fits, weights = NULL) {
  stopifnot(is.list(fits), length(fits) >= 2)
  if (is.null(weights)) weights <- rep(1/length(fits), length(fits))
  weights <- weights / sum(weights)
  K <- length(fits[[1]]$names)
  if (!is.null(fits[[1]]$basis)) {
    Mbar <- vector("list", K)
    for (k in seq_len(K)) Mbar[[k]] <- Reduce(`+`, Map(function(w, f) w * f$M[[k]], weights, fits))
    out <- fits[[1]]; out$M <- Mbar; out
  } else {
    Dbar <- vector("list", K)
    for (k in seq_len(K)) Dbar[[k]] <- Reduce(`+`, Map(function(w, f) w * f$Delta[[k]], weights, fits))
    out <- fits[[1]]; out$Delta <- Dbar; out
  }
}


# -------- Mechanistic readouts (neural-aware) -------------------------------

#' Precision/gain index in bold or neural domain (basis or ROI)
mppi_gain <- function(fit, k = 1L, lambda = 1e-3, roi_idx = NULL) {
  if (!is.null(fit$basis)) {
    X <- if (!is.null(fit$U)) fit$U else fit$Z
    S0 <- crossprod(X) / nrow(X)
    Mk <- if (!is.null(fit$U)) crossprod(fit$U, fit$pk[[k]] * fit$U) / sum(fit$pk[[k]]^2)
          else                 crossprod(fit$Z, fit$pk[[k]] * fit$Z) / sum(fit$pk[[k]]^2)
  } else {
    X <- if (!is.null(fit$U)) fit$U else fit$R
    S0 <- crossprod(X) / nrow(X)
    Mk <- if (!is.null(fit$U)) crossprod(fit$U, fit$pk[[k]] * fit$U) / sum(fit$pk[[k]]^2)
          else                 crossprod(fit$R, fit$pk[[k]] * fit$R) / sum(fit$pk[[k]]^2)
  }
  diag(Mk) <- diag(Mk)  # keep diag here; we need it for precision change
  Theta0 <- solve(S0 + diag(lambda, nrow(S0)))
  dTheta <- - Theta0 %*% Mk %*% Theta0
  comp_gain <- diag(dTheta)
  if (is.null(fit$basis)) {
    if (!is.null(roi_idx)) comp_gain <- comp_gain[roi_idx]
    return(list(gain = comp_gain, space = "roi"))
  } else {
    if (!is.null(roi_idx)) {
      V <- fit$basis$V[roi_idx, , drop = FALSE]
      g_roi <- rowSums((V %*% dTheta) * V)
      return(list(gain_component = comp_gain, gain_roi = g_roi, idx = roi_idx, space = "roi+basis"))
    } else {
      return(list(gain_component = comp_gain, space = "basis"))
    }
  }
}

#' Routing asymmetry index using lagged fits and a hierarchy vector (basis or ROI)
mppi_routing_index <- function(fits_by_lag, k = 1L, hierarchy, pos = c(1,2), neg = c(-1,-2), eps = 1e-8) {
  f0 <- fits_by_lag[[1]]
  is_basis <- !is.null(f0$basis)
  u <- as.numeric(hierarchy)
  if (is_basis) {
    stopifnot(length(u) == ncol(f0$Z) || (!is.null(f0$U) && length(u) == ncol(f0$U)))
    Wdir <- tcrossprod(u, u)
    proj <- function(M) sum(M * Wdir, na.rm = TRUE)
    getM <- function(fit) {
      if (!is.null(fit$U)) crossprod(fit$U, fit$pk[[k]] * fit$U) / sum(fit$pk[[k]]^2)
      else                 crossprod(fit$Z, fit$pk[[k]] * fit$Z) / sum(fit$pk[[k]]^2)
    }
  } else {
    stopifnot(length(u) == ncol(f0$R) || (!is.null(f0$U) && length(u) == ncol(f0$U)))
    Wdir <- tcrossprod(u, u)
    proj <- function(D) sum(D * Wdir, na.rm = TRUE)
    getM <- function(fit) {
      if (!is.null(fit$U)) crossprod(fit$U, fit$pk[[k]] * fit$U) / sum(fit$pk[[k]]^2)
      else                 crossprod(fit$R, fit$pk[[k]] * fit$R) / sum(fit$pk[[k]]^2)
    }
  }
  Epos <- 0; Eneg <- 0
  for (nm in names(fits_by_lag)) {
    τ <- as.integer(nm)
    Mk <- getM(fits_by_lag[[nm]])
    if (τ %in% pos)  Epos <- Epos + proj(Mk)
    if (τ %in% neg)  Eneg <- Eneg + proj(Mk)
  }
  R <- (Epos - Eneg) / (Epos + Eneg + eps)
  list(R = R, Epos = Epos, Eneg = Eneg, pos = pos, neg = neg)
}

#' Template projection ("DCM-lite")
mppi_project_templates <- function(fit, k = 1L, templates, normalize = TRUE) {
  get_Mfull <- function(f) {
    if (!is.null(f$basis)) {
      if (!is.null(f$U)) crossprod(f$U, f$pk[[k]] * f$U) / sum(f$pk[[k]]^2)
      else               crossprod(f$Z, f$pk[[k]] * f$Z) / sum(f$pk[[k]]^2)
    } else {
      if (!is.null(f$U)) crossprod(f$U, f$pk[[k]] * f$U) / sum(f$pk[[k]]^2)
      else               crossprod(f$R, f$pk[[k]] * f$R) / sum(f$pk[[k]]^2)
    }
  }
  Mfull <- get_Mfull(fit)
  sapply(templates, function(W) {
    if (normalize) sum(Mfull * W, na.rm = TRUE) / sqrt(sum(W^2))
    else           sum(Mfull * W, na.rm = TRUE)
  })
}
