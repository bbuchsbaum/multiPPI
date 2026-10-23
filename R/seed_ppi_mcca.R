#' Multivariate seed–PPI via partial CCA (searchlight-friendly)
#'
#' Computes a searchlight map of multivariate seed–PPI using ridge-regularized
#' partial CCA between a seed set (time-series or indices) and the local
#' interaction block Y:Design within each neighborhood. Neighborhoods can be
#' provided directly (as a list of integer column indices over the data matrix)
#' or constructed via neuroim2::searchlight on a NeuroVol mask.
#'
#' Data input can be either a T x V matrix or a complete-grid neuroim2::NeuroVec.
#' For clustered data, supply the time by cluster matrix and explicit neighborhoods.
#' When using neuroim2::searchlight, provide a compatible NeuroVol mask and ensure
#' the NeuroVec shares the same space.
#'
#' @param data   T x V numeric matrix, or a complete-grid neuroim2::NeuroVec.
#' @param seed   Seed specification: integer indices, logical mask (over V), label(s),
#'               or a T x q numeric matrix.
#' @param design Numeric design matrix (T x m) or fmridesign::baseline_model.
#' @param nuisance Optional numeric nuisance matrix (T x r) or baseline_model.
#' @param neighborhoods Optional list of integer vectors (column indices into data).
#' @param mask  Optional neuroim2::NeuroVol mask used to build neighborhoods via
#'              neuroim2::searchlight when neighborhoods is NULL.
#' @param searchlight list(radius=6, eager=TRUE, nonzero=TRUE) controls neighborhood building
#'                    when using neuroim2::searchlight.
#' @param tr     TR (seconds). If NULL and `data` is a NeuroVec, attempts to read from header.
#' @param runs Run labels, one per observation.
#' @param center_by Whether to center nonconstant design columns within runs.
#' @param domain Working domain; currently only descriptive `hemo` is supported.
#' @param prewhiten logical; if TRUE, prewhiten time-series via fmriAR before modeling.
#' @param ar_method Noise model delegated to fmriAR.
#' @param ar_order AR order or `"auto"`, delegated to fmriAR.
#' @param pca    list(seed_k=5, light_k=10, center=TRUE, scale=TRUE) PCA reductions.
#' @param reg    list(l2_x=0.05, l2_y=0.05) ridge penalties for seed and interaction blocks.
#' @param stat   c("r1","wilks") – report first canonical correlation or true multi-root Wilks (−∑log(1−ρ_j²)).
#' @param wilks_k Optional number of roots included in the descriptive Wilks statistic.
#' @param inference Only `method="none"` is supported; CCA inference is unqualified.
#' @param return  character subset of c("map","weights") controlling outputs.
#' @return list with class 'mppi_sl' containing: stat, pval (optional), centers,
#'         neighborhoods (indices), seeds (if resolved), weights_x/y (optional), meta.
#' @export
seed_ppi_mcca <- function(data,
                          seed,
                          design,
                          nuisance = NULL,
                          neighborhoods = NULL,
                          mask = NULL,
                          searchlight = list(radius = 6, eager = TRUE, nonzero = TRUE),
                          tr = NULL,
                          runs = NULL,
                          center_by = c("none","run"),
                          domain = c("hemo","neural"),
                          prewhiten = TRUE,
                          ar_method = "ar",
                          ar_order = "auto",
                          pca = list(seed_k = 5, light_k = 10, center = TRUE, scale = TRUE),
                          reg = list(l2_x = 0.05, l2_y = 0.05),
                          stat = c("r1","wilks"),
                          wilks_k = NULL,
                          inference = list(method = c("none","param","perm"), perms = 0L, block = NULL),
                          return = c("map","weights")) {

  stat <- match.arg(stat)
  domain <- match.arg(domain)
  center_by <- match.arg(center_by)
  method_allowed <- c("none","param","perm")
  requested_method <- inference$method %||% "none"
  if (length(requested_method) == 1L && requested_method %in% c("parametric", "permutation")) requested_method <- if (requested_method == "parametric") "param" else "perm"
  inference$method <- match.arg(requested_method, method_allowed)
  if (inference$method != "none") stop("Seed CCA inference is unavailable: data-dependent nuisance, ridge and temporal dependence require a qualified null model.", call. = FALSE)
  if (domain == "neural") stop("Seed CCA neural deconvolution is unsupported pending a qualified run-local construction.", call. = FALSE)
  if (!is.null(wilks_k)) {
    if (!is.numeric(wilks_k) || length(wilks_k) != 1L || !is.finite(wilks_k) || wilks_k < 1) {
      stop("wilks_k must be NULL or a positive integer.", call. = FALSE)
    }
    wilks_k <- as.integer(wilks_k)
  }
  return <- intersect(return, c("map","weights"))
  if (!length(return)) return <- "map"

  D0 <- .sl_as_matrix(data)
  Xmat <- D0$X
  V <- ncol(Xmat)
  TR <- tr %||% D0$tr %||% NA_real_

  # Build neighborhoods
  neigh <- neighborhoods
  centers <- NULL
  if (is.null(neigh)) {
    if (!requireNamespace("neuroim2", quietly = TRUE)) {
      stop("Provide 'neighborhoods' or install 'neuroim2' to build via searchlight.", call. = FALSE)
    }
    if (is.null(mask)) stop("'mask' (NeuroVol) must be provided when neighborhoods=NULL.", call. = FALSE)
    rad <- as.numeric(searchlight$radius %||% 6)
    eager <- isTRUE(searchlight$eager %||% TRUE)
    nonzero <- isTRUE(searchlight$nonzero %||% TRUE)
    sl <- neuroim2::searchlight(mask, radius = rad, eager = eager, nonzero = nonzero)
    neigh <- .sl_searchlight_to_indices(sl, mask)
    # Centers correspond to mask nonzero linear indices in order
    centers <- which(mask != 0)
    if (length(neigh) != length(centers)) {
      # fallback: derive centers from neighborhood first element
      centers <- vapply(neigh, function(id) if (length(id)) as.integer(id[[1]]) else NA_integer_, 1L)
    }
  } else {
    # Validate provided neighborhoods
    if (!is.list(neigh) || !all(vapply(neigh, function(x) (is.integer(x) || is.numeric(x)) && all(x == round(x), na.rm = TRUE), logical(1)))) {
      stop("'neighborhoods' must be a list of integer vectors (column indices).", call. = FALSE)
    }
    neigh <- lapply(neigh, as.integer)
    centers <- vapply(neigh, function(id) if (length(id)) as.integer(id[[1]]) else NA_integer_, 1L)
  }

  # External seed matrices are raw data and receive the same fitted operator.
  external <- is.matrix(seed)
  if (external && (nrow(seed) != nrow(Xmat) || !is.numeric(seed))) stop("External seed must be a numeric matrix with one row per observation.", call. = FALSE)
  X_input <- if (external) cbind(Xmat, seed) else Xmat
  # Preprocess
  prep <- .sl_preprocess(X_input, design, nuisance,
                         prewhiten = prewhiten, TR = TR, domain = domain,
                         runs = runs, center_by = center_by,
                         ar_method = ar_method, ar_order = ar_order, noise_cols = V)
  Y_all <- prep$Y[, seq_len(V), drop = FALSE]
  seed_processed <- if (external) prep$Y[, V + seq_len(ncol(seed)), drop = FALSE] else seed
  D_mat <- prep$D
  N_mat <- prep$N

  # Seed matrix X (T x q)
  X_seed <- .sl_seed_matrix(seed_processed, Y_all, pca = pca)

  keep_weights <- "weights" %in% return
  out_stat <- rep(NA_real_, length(neigh))
  out_p <- if (inference$method != "none") rep(NA_real_, length(neigh)) else NULL
  WX <- if (keep_weights) vector("list", length(neigh)) else NULL
  WY <- if (keep_weights) vector("list", length(neigh)) else NULL

  # Pre-center D and N once
  D_cent <- D_mat
  N_cent <- if (!is.null(N_mat)) N_mat else NULL

  outside <- vapply(neigh, function(id) anyNA(id) || any(id < 1L | id > ncol(Y_all)), logical(1))
  if (any(outside)) stop(sprintf("%d neighborhood(s) contain column indices outside 1..%d.", sum(outside), ncol(Y_all)), call. = FALSE)

  for (i in seq_along(neigh)) {
    idx <- neigh[[i]]
    if (length(idx) < 2L) next

    # Searchlight data and PCA
    Y_l <- Y_all[, idx, drop = FALSE]
    Y_l <- .pca_reduce(Y_l, k = pca$light_k, center = pca$center, scale = pca$scale)

    # Interaction block W = Y_l : D, with informative column names
    # name searchlight PCs then broadcast
    if (is.null(colnames(Y_l))) colnames(Y_l) <- paste0("YPC", seq_len(ncol(Y_l)))
    W <- .kr_broadcast(Y_l, D_cent)

    # Residualize X and W against Z = [Y_l, D, N] using QR (avoids T×T projector)
    # Center Y_l by strategy
    Yl_cent <- if (identical(center_by, "run") && !is.null(runs)) {
      apply(Y_l, 2L, function(col) .mppi_center_by_run(col, runs))
    } else {
      .col_center(Y_l)
    }
    Z <- if (is.null(N_cent)) cbind(Yl_cent, D_cent) else cbind(Yl_cent, D_cent, N_cent)
    if (ncol(Z) > 0L) {
      qrZ <- qr(Z)
      # Guardrails: residual DOF and size ratios
      n_eff <- nrow(Z) - qrZ$rank
      p_eff <- ncol(X_seed)
      q_eff <- ncol(W)
      if (n_eff <= max(p_eff, q_eff) + 2L) {
        warning(sprintf("Low effective DOF after partialling (n_eff=%d) relative to dims (p=%d, q=%d). Consider reducing PCA dims or ridge.",
                        n_eff, p_eff, q_eff), call. = FALSE)
      }
      if ((qrZ$rank / nrow(Z)) > 0.8) {
        warning(sprintf("Projector rank (rank(Z)=%d) is high relative to T=%d; stability may suffer.", qrZ$rank, nrow(Z)), call. = FALSE)
      }
      Xr <- .col_center(X_seed) - qr.fitted(qrZ, .col_center(X_seed))
      Wr <- .col_center(W)      - qr.fitted(qrZ, .col_center(W))
    } else {
      Xr <- .col_center(X_seed)
      Wr <- .col_center(W)
    }

    cca <- .partial_cca(Xr, Wr, l2_x = reg$l2_x %||% 0.05, l2_y = reg$l2_y %||% 0.05)
    # True Wilks uses all roots: −∑log(1−ρ_j²)
    s <- if (identical(stat, "r1")) {
      cca$r1
    } else {
      rhos <- pmin(cca$rhos, 0.999)
      k_use <- if (is.null(wilks_k)) length(rhos) else min(wilks_k, length(rhos))
      rhos <- rhos[seq_len(k_use)]
      -sum(log(pmax(1 - rhos^2, 1e-12)))
    }
    out_stat[i] <- s

    if (keep_weights) {
      WX[[i]] <- drop(cca$ax)
      WY[[i]] <- drop(cca$bw)
    }

    if (identical(inference$method, "param")) {
      out_p[i] <- .cca_pvalue_param(cca$rhos, n = nrow(Xr), p = ncol(Xr), q = ncol(Wr), k = wilks_k)
    } else if (identical(inference$method, "perm") && (inference$perms %||% 0L) > 0L) {
      out_p[i] <- .cca_pvalue_perm(Xr, Wr, perms = as.integer(inference$perms), block = inference$block, stat = stat, k = wilks_k)
    }
  }

  res <- list(
    stat = out_stat,
    pval = out_p,
    centers = centers,
    neighborhoods = neigh,
    seeds = attr(X_seed, "seed_idx") %||% NA_integer_,
    weights_x = WX,
    weights_y = WY,
    meta = list(
      stat = stat,
      domain = domain,
      prewhiten = prewhiten,
      pca = pca,
      reg = reg,
      TR = TR,
      dims = D0$dims,
      design_cols = colnames(D_mat),
      nuisance_cols = if (!is.null(N_mat)) colnames(N_mat) else NULL,
      inference = "descriptive_only", seed_input = if (external) "raw_matrix" else "indexed_data"
    )
  )
  class(res) <- c("mppi_sl","list")
  res
}


# ---------- Small helpers (file-local) --------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

.sl_as_matrix <- function(data) {
  if (is.matrix(data)) {
    list(X = data, dims = NULL, tr = attr(data, "TR") %||% NULL)
  } else if (inherits(data, "NeuroVec") || inherits(data, "ClusteredNeuroVec")) {
    if (!requireNamespace("neuroim2", quietly = TRUE)) {
      stop("'neuroim2' is required to use NeuroVec inputs.", call. = FALSE)
    }
    if (inherits(data, "ClusteredNeuroVec")) stop("ClusteredNeuroVec uses cluster indices; pass its time by cluster matrix and explicit neighborhoods.", call. = FALSE)
    X <- t(neuroim2::as.matrix(data))
    if (nrow(X) != dim(data)[4L] || ncol(X) != prod(dim(data)[1:3])) stop("NeuroVec matrix does not expose a complete voxel grid.", call. = FALSE)
    spacing <- neuroim2::spacing(neuroim2::space(data))
    tr <- if (length(spacing) >= 4L) spacing[4L] else NULL
    list(X = X, dims = dim(data), tr = tr)
  } else {
    stop("Unsupported 'data' type. Provide T x V matrix or NeuroVec.", call. = FALSE)
  }
}

.sl_searchlight_to_indices <- function(sl, mask) {
  # sl: either deflist or eager list of ROI-like objects
  n <- tryCatch(length(sl), error = function(e) NA_integer_)
  if (!is.finite(n) || n <= 0L) stop("Empty searchlight result.", call. = FALSE)
  get_idx <- function(roi) {
    # most robust attempts to obtain integer indices
    idx <- tryCatch(neuroim2::indices(roi), error = function(e) NULL)
    if (is.null(idx)) idx <- tryCatch(roi$indices, error = function(e) NULL)
    if (is.null(idx)) {
      # final fallback: assume NeuroVol-like and use which(roi != 0)
      idx <- tryCatch(which(roi != 0), error = function(e) NULL)
    }
    if (is.null(idx)) stop("Could not extract indices from searchlight ROI.", call. = FALSE)
    as.integer(idx)
  }
  out <- vector("list", n)
  for (i in seq_len(n)) {
    roi <- sl[[i]]  # deflist will materialize lazily
    out[[i]] <- get_idx(roi)
  }
  out
}

.design_as_matrix <- function(obj) {
  if (inherits(obj, c("baseline_model","event_model","fmridesign","mppi_design"))) {
    if (!requireNamespace("fmrireg", quietly = TRUE)) {
      stop("Install 'fmrireg' or pass a numeric design matrix.", call. = FALSE)
    }
    as.matrix(fmrireg::design_matrix(obj))
  } else if (is.matrix(obj)) {
    obj
  } else {
    stop("Provide fmridesign::baseline_model or numeric matrix for 'design'.", call. = FALSE)
  }
}

.sl_preprocess <- function(X, design, nuisance,
                           prewhiten = TRUE, TR = NULL, domain = "neural",
                           runs = NULL, center_by = c("none","run"),
                           ar_method = "ar", ar_order = "auto", noise_cols = ncol(X)) {
  center_by <- match.arg(center_by)
  D <- .design_as_matrix(design)
  if (nrow(D) != nrow(X) || any(!is.finite(X)) || any(!is.finite(D))) stop("Data and design must have matching rows and finite values; censoring is unsupported.", call. = FALSE)
  if (!is.null(runs) && (length(runs) != nrow(X) || anyNA(runs))) stop("Runs must have one nonmissing label per observation.", call. = FALSE)
  if (identical(center_by, "run") && !is.null(runs)) {
    for (j in seq_len(ncol(D))) D[, j] <- .mppi_center_by_run(D[, j], runs)
  } else {
    D <- .col_center(D)
  }
  N <- if (!is.null(nuisance)) {
    Nm <- .design_as_matrix(nuisance)
    if (nrow(Nm) != nrow(X) || any(!is.finite(Nm))) stop("Nuisance must have matching rows and finite values.", call. = FALSE)
    if (identical(center_by, "run") && !is.null(runs)) {
      for (j in seq_len(ncol(Nm))) Nm[, j] <- .mppi_center_by_run(Nm[, j], runs)
      Nm
    } else {
      .col_center(Nm)
    }
  } else NULL

  Y <- X
  # Prewhiten Y and Z0=[D,N] together for consistency
  if (isTRUE(prewhiten)) {
    if (is.null(runs)) stop("'runs' must be provided when prewhiten = TRUE.", call. = FALSE)
    Z0 <- if (is.null(N)) D else cbind(D, N)
    # Fit AR on residuals of Y after projecting Z0 (per standard practice)
    qrZ <- if (ncol(Z0)) qr(Z0) else NULL
    res <- if (!is.null(qrZ)) qr.resid(qrZ, Y[, seq_len(noise_cols), drop = FALSE]) else Y[, seq_len(noise_cols), drop = FALSE]
    plan <- fmriAR::fit_noise(res, runs = runs, method = ar_method, p = ar_order)
    xyw <- fmriAR::whiten_apply(plan, X = if (is.null(Z0)) matrix(0, nrow(Y), 0) else Z0, Y = Y, runs = runs)
    Y <- xyw$Y
    if (!is.null(Z0) && ncol(Z0)) {
      Z0w <- xyw$X
      # Split back into D and N
      nd <- ncol(D)
      if (is.null(N)) {
        D <- Z0w
      } else {
        D <- Z0w[, seq_len(nd), drop = FALSE]
        N <- Z0w[, -(seq_len(nd)), drop = FALSE]
      }
    }
  }

  if (identical(domain, "neural")) stop("Seed CCA neural deconvolution is unsupported pending a qualified run-local construction.", call. = FALSE)

  # Globally remove [D,N] main effects from Y via QR residuals
  if (!is.null(N)) {
    qrZN <- qr(cbind(D, N))
    Y <- Y - qr.fitted(qrZN, Y)
  } else if (!is.null(D)) {
    qrD <- qr(D)
    Y <- Y - qr.fitted(qrD, Y)
  }
  list(Y = Y, D = D, N = N)
}

.sl_seed_matrix <- function(seed, Y_all, pca) {
  Tn <- nrow(Y_all)
  if (is.matrix(seed) && nrow(seed) == Tn) {
    X <- seed
    seed_idx <- NULL
  } else if (is.numeric(seed) || is.logical(seed)) {
    idx <- if (is.logical(seed)) which(seed) else as.integer(seed)
    X <- Y_all[, idx, drop = FALSE]
    seed_idx <- idx
  } else {
    stop("Unsupported 'seed'. Use indices, logical mask, or T x q matrix.", call. = FALSE)
  }
  X <- .pca_reduce(X, k = pca$seed_k, center = pca$center, scale = pca$scale)
  attr(X, "seed_idx") <- seed_idx
  X
}

.col_center <- function(M) sweep(M, 2L, colMeans(M), `-`)

.pca_reduce <- function(M, k, center = TRUE, scale = TRUE) {
  if (isTRUE(center)) M <- .col_center(M)
  if (isTRUE(scale)) {
    # Robust column scaling: avoid division by zero for constant columns
    sds <- sqrt(pmax(matrixStats::colVars(as.matrix(M)), .Machine$double.eps))
    M <- sweep(M, 2L, sds, "/")
  }
  k <- as.integer(k)
  if (!is.finite(k) || k <= 0L) return(M)
  if (ncol(M) <= k) return(M)
  k_use <- min(k, ncol(M))
  sv <- suppressWarnings(svd(M, nu = 0, nv = k_use))
  if (k_use == 0L) return(M)
  M %*% sv$v
}

.kr_broadcast <- function(Y, D) {
  blocks <- lapply(seq_len(ncol(D)), function(k) {
    Y * D[, k]
  })
  out <- do.call(cbind, blocks)
  # informative column names: YPCj:DesignName
  y_nms <- colnames(Y)
  d_nms <- colnames(D)
  if (!is.null(y_nms) && !is.null(d_nms)) {
    colnames(out) <- as.vector(outer(y_nms, d_nms, function(y, d) paste0(y, ":", d)))
  }
  out
}

.projector_from_qr <- function(Z) {
  qrZ <- qr(Z)
  Q <- qr.Q(qrZ)
  diag(nrow(Z)) - Q %*% t(Q)
}

.partial_cca <- function(X, W, l2_x = 0.05, l2_y = 0.05) {
  n <- nrow(X); q <- ncol(X); r <- ncol(W)
  Sxx <- crossprod(X) / (n - 1) + diag(l2_x, q)
  Sww <- crossprod(W) / (n - 1) + diag(l2_y, r)
  Sxw <- crossprod(X, W) / (n - 1)
  cx <- chol(Sxx); cy <- chol(Sww)
  K <- solve(t(cx), Sxw) %*% solve(cy)
  sv <- svd(K, nu = 1, nv = 1)
  list(r1 = sv$d[1], rhos = sv$d, ax = solve(cx, sv$u[, 1, drop = FALSE]), bw = solve(cy, sv$v[, 1, drop = FALSE]))
}

.cca_pvalue_param <- function(rhos, n, p, q, k = NULL) {
  if (length(n) != 1L || length(p) != 1L || length(q) != 1L || any(!is.finite(c(n,p,q))) || p < 1 || q < 1 || p != as.integer(p) || q != as.integer(q) || n <= 1 + (p+q+1)/2) stop("Invalid classical CCA dimensions or sample size.", call. = FALSE)
  roots <- min(p,q)
  if (length(rhos) != roots || any(!is.finite(rhos)) || any(rhos < 0 | rhos >= 1) || (!is.null(k) && (length(k) != 1L || k != roots))) stop("Classical CCA test requires all canonical correlations in [0, 1).", call. = FALSE)
  chi <- -((n-1)-(p+q+1)/2) * sum(log1p(-rhos^2))
  stats::pchisq(chi, df=p*q, lower.tail=FALSE)
}

.cca_pvalue_perm <- function(Xr, Wr, perms = 500L, block = NULL, stat = c("r1","wilks"), k = NULL) {
  stat <- match.arg(stat)
  cca_obs <- .partial_cca(Xr, Wr)
  obs <- if (identical(stat, "r1")) cca_obs$r1 else {
    rhos <- pmin(cca_obs$rhos, 0.999)
    k_use <- if (is.null(k)) length(rhos) else min(k, length(rhos))
    rhos <- rhos[seq_len(k_use)]
    -sum(log(pmax(1 - rhos^2, 1e-12)))
  }
  Tn <- nrow(Xr)
  rvals <- numeric(perms)
  for (b in seq_len(perms)) {
    if (is.null(block)) {
      sgn <- sample(c(-1, 1), size = Tn, replace = TRUE)
      Xb <- Xr * sgn; Wb <- Wr * sgn
    } else {
      # crude block permutation
      idx <- split(seq_len(Tn), ceiling(seq_len(Tn) / block))
      pb <- sample(seq_along(idx))
      ord <- unlist(idx[pb], use.names = FALSE)
      Xb <- Xr[ord, , drop = FALSE]
      Wb <- Wr[ord, , drop = FALSE]
    }
    cca_b <- .partial_cca(Xb, Wb)
    rvals[b] <- if (identical(stat, "r1")) cca_b$r1 else {
      rhos <- pmin(cca_b$rhos, 0.999)
      k_use <- if (is.null(k)) length(rhos) else min(k, length(rhos))
      rhos <- rhos[seq_len(k_use)]
      -sum(log(pmax(1 - rhos^2, 1e-12)))
    }
  }
  mean(rvals >= obs)
}
