# Seed-set utilities -------------------------------------------------------

.mppi_seedset_dim <- function(fit) {
  if (!is.null(fit$basis)) {
    return(nrow(fit$basis$V))
  }
  R_like <- fit$R_raw
  if (is.null(R_like)) R_like <- fit$R
  if (is.null(R_like)) stop("Unable to determine spatial dimension from fit.", call. = FALSE)
  ncol(R_like)
}

.mppi_seedset_labels <- function(fit) {
  labs <- NULL
  if (!is.null(fit$basis)) {
    labs <- rownames(fit$basis$V)
    if (is.null(labs) && !is.null(dimnames(fit$basis$V)[[1L]])) {
      labs <- dimnames(fit$basis$V)[[1L]]
    }
  }
  if (is.null(labs) && !is.null(fit$R_raw)) {
    labs <- colnames(fit$R_raw)
  }
  if (is.null(labs) && !is.null(fit$R)) {
    labs <- colnames(fit$R)
  }
  labs
}

.mppi_as_index <- function(x, n, labels = NULL, what = "index") {
  if (is.null(x)) stop(sprintf("%s must be provided", what), call. = FALSE)
  if (is.logical(x)) {
    if (length(x) != n) stop(sprintf("Logical %s must have length %d", what, n), call. = FALSE)
    idx <- which(x)
  } else if (is.numeric(x)) {
    idx <- as.integer(round(x))
    if (any(idx < 1L | idx > n)) stop(sprintf("Numeric %s outside 1:%d", what, n), call. = FALSE)
  } else if (is.character(x)) {
    if (is.null(labels)) {
      stop(sprintf("Character %s requested but no labels available", what), call. = FALSE)
    }
    idx <- match(x, labels)
    if (anyNA(idx)) {
      missing <- x[is.na(idx)]
      stop(sprintf("Unknown %s label(s): %s", what, paste(missing, collapse = ", ")), call. = FALSE)
    }
  } else {
    stop(sprintf("Unsupported %s specification", what), call. = FALSE)
  }
  if (!length(idx)) stop(sprintf("No %s selected", what), call. = FALSE)
  unique(idx)
}

.mppi_seedset_slice <- function(M, fit, seed_idx, target_idx) {
  if (is.list(M) && !is.null(M$values)) {
    M <- .mppi_unpack_upper(M)
  }
  if (is.null(fit$basis)) {
    out <- M[seed_idx, target_idx, drop = FALSE]
  } else {
    V <- fit$basis$V
    V_seed <- V[seed_idx, , drop = FALSE]
    V_target <- V[target_idx, , drop = FALSE]
    out <- V_seed %*% M %*% t(V_target)
  }
  if (length(seed_idx) == length(target_idx) && all(seed_idx == target_idx)) {
    diag(out) <- 0
  }
  out
}

.mppi_seedset_precision_diag <- function(fit, Theta, seed_idx) {
  if (!is.null(fit$basis)) {
    V <- fit$basis$V
    V_seed <- V[seed_idx, , drop = FALSE]
    rowSums((V_seed %*% Theta) * V_seed)
  } else {
    diag(Theta)[seed_idx]
  }
}

#' Seed-set → whole-brain slopes for a condition/contrast
#'
#' @param fit `mppi_fit` object.
#' @param seeds Indices, logical mask, or labels identifying the seed rows.
#' @param k Condition/contrast name or index.
#' @param target_idx Optional target indices (defaults to all locations).
#' @param mode Scale of the interaction matrix: `"normalized"`, `"amplitude"`, or `"raw"`.
#' @param lag Optional lag (offset) to evaluate; defaults to contemporaneous (0).
#' @return Matrix of dimension |seeds| × |targets|.
#' @export
mppi_seedset_slopes <- function(fit, seeds, k, target_idx = NULL,
                                mode = c("normalized", "amplitude", "raw"),
                                lag = 0L) {
  mode <- match.arg(mode)
  V <- .mppi_seedset_dim(fit)
  labels <- .mppi_seedset_labels(fit)
  seed_idx <- .mppi_as_index(seeds, V, labels, what = "seeds")
  if (is.null(target_idx)) {
    target_idx <- seq_len(V)
  }
  target_idx <- .mppi_as_index(target_idx, V, labels, what = "targets")
  lag <- as.integer(lag)
  if (lag == 0L) {
    Mk <- mppi_get_M_scaled(fit, k, mode)
  } else {
    stored_mode <- fit$scale %||% "normalized"
    Mk_lag <- mppi_get_M_lag(fit, k, lag)
    Mk <- if (mode == stored_mode) {
      Mk_lag
    } else {
      .mppi_convert_scale(Mk_lag, fit$sigma, stored_mode, mode)
    }
  }
  slopes <- .mppi_seedset_slice(Mk, fit, seed_idx, target_idx)
  seed_names <- if (!is.null(labels)) labels[seed_idx] else paste0("seed_", seed_idx)
  target_names <- if (!is.null(labels)) labels[target_idx] else paste0("v", target_idx)
  rownames(slopes) <- seed_names
  colnames(slopes) <- target_names
  attr(slopes, "seed_idx") <- seed_idx
  attr(slopes, "target_idx") <- target_idx
  k_label <- if (is.character(k)) k else fit$names[[k]]
  if (is.null(k_label)) k_label <- as.character(k)
  attr(slopes, "k") <- k_label
  attr(slopes, "mode") <- mode
  attr(slopes, "lag") <- lag
  slopes
}

#' Collapse seed-set slopes into a single map
#'
#' @param fit `mppi_fit` object.
#' @param seeds Seed specification (as in [mppi_seedset_slopes()]).
#' @param k Condition/contrast name or index.
#' @param aggregate Aggregation rule: `"svd"`, `"mean"`, or `"precision"`.
#' @param target_idx Optional target indices (defaults to all locations).
#' @param mode Scale for the interaction matrix (`"normalized"`, `"amplitude"`, `"raw"`).
#' @param Theta0 Optional baseline precision matrix (component or voxel space) for
#'   `aggregate = "precision"`.
#' @param lag Optional lag (offset) to evaluate; defaults to contemporaneous (0).
#' @return Object of class `mppi_seedset_map` containing the aggregated map and weights.
#' @export
mppi_seedset_map <- function(fit, seeds, k,
                             aggregate = c("svd", "mean", "precision"),
                             target_idx = NULL,
                             mode = c("normalized", "amplitude", "raw"),
                             Theta0 = NULL,
                             lag = 0L) {
  aggregate <- match.arg(aggregate)
  mode <- match.arg(mode)
  slopes <- mppi_seedset_slopes(fit, seeds, k, target_idx = target_idx, mode = mode, lag = lag)
  seed_idx <- attr(slopes, "seed_idx")
  target_idx_vals <- attr(slopes, "target_idx")
  seed_names <- rownames(slopes)
  k_label <- attr(slopes, "k")
  if (is.null(k_label)) k_label <- as.character(k)
  if (nrow(slopes) == 1L) {
    weights <- 1
    map <- as.numeric(slopes)
  } else if (aggregate == "mean") {
    weights <- rep(1 / nrow(slopes), nrow(slopes))
    map <- drop(colMeans(slopes))
  } else if (aggregate == "svd") {
    centered <- scale(slopes, center = TRUE, scale = FALSE)
    sv <- La.svd(centered)
    weights <- sv$u[, 1L]
    # Normalize weights to unit L2 norm for stability
    weights <- weights / sqrt(sum(weights^2))
    map <- drop(weights %*% slopes)
  } else {
    Theta_use <- Theta0
    if (is.null(Theta_use)) {
      Theta_use <- fit$Theta0
    }
    if (is.null(Theta_use)) {
      Sigma0 <- fit$Sigma0
      if (is.null(Sigma0)) stop("Baseline covariance unavailable; provide Theta0 for precision weighting.", call. = FALSE)
      lambda <- fit$partial_lambda
      if (is.null(lambda)) lambda <- 1e-3
      Theta_use <- tryCatch(solve(Sigma0 + lambda * diag(nrow(Sigma0))),
                            error = function(e) stop("Failed to invert baseline covariance for precision weights.", call. = FALSE))
    }
    diag_vals <- .mppi_seedset_precision_diag(fit, Theta_use, seed_idx)
    diag_vals[diag_vals <= 0] <- min(diag_vals[diag_vals > 0], na.rm = TRUE)
    weights <- 1 / diag_vals
    weights <- weights / sum(weights)
    map <- drop(weights %*% slopes)
  }
  target_names <- colnames(slopes)
  out <- list(
    map = setNames(as.numeric(map), target_names),
    seeds = seed_idx,
    target_idx = target_idx_vals,
    weights = setNames(as.numeric(weights), seed_names),
    targets = target_names,
    k = k_label,
    mode = mode,
    aggregate = aggregate,
    lag = attr(slopes, "lag")
  )
  class(out) <- c("mppi_seedset_map", "list")
  out
}

#' Seed slopes convenience wrapper
#'
#' @param fit `mppi_fit` object (including instantaneous engine).
#' @param seeds Seed definition passed to [mppi_seedset_slopes()].
#' @param k Condition/contrast index or name.
#' @param collapse `"none"` to return one row per seed, `"mean"` to average seeds.
#' @param mode Scale of the interaction matrix (`"normalized"`, `"amplitude"`, `"raw"`).
#' @param lag Optional lag offset (defaults to 0).
#' @return Matrix of seed slopes (potentially collapsed to a single row).
#' @export
mppi_seed_slopes <- function(fit, seeds, k = 1L,
                             collapse = c("none", "mean"),
                             mode = c("normalized", "amplitude", "raw"),
                             lag = 0L) {
  collapse <- match.arg(collapse)
  mode <- match.arg(mode)
  slopes <- mppi_seedset_slopes(fit, seeds, k, mode = mode, lag = lag)
  if (collapse == "mean" && nrow(slopes) > 1L) {
    collapsed <- matrix(colMeans(slopes), nrow = 1L)
    colnames(collapsed) <- colnames(slopes)
    rownames(collapsed) <- "mean"
    attr(collapsed, "seed_idx") <- attr(slopes, "seed_idx")
    attr(collapsed, "target_idx") <- attr(slopes, "target_idx")
    attr(collapsed, "k") <- attr(slopes, "k")
    attr(collapsed, "mode") <- attr(slopes, "mode")
    attr(collapsed, "lag") <- attr(slopes, "lag")
    return(collapsed)
  }
  slopes
}

#' @export
print.mppi_seedset_map <- function(x, ...) {
  cat(sprintf("mPPI seed-set map [%s] aggregate=%s |seeds|=%d targets=%d\n",
              x$mode, x$aggregate, length(x$seeds), length(x$map)))
  cat(sprintf("condition/contrast: %s\n", as.character(x$k)))
  if (!is.null(x$weights)) {
    cat("seed weights:\n")
    print(utils::head(sort(x$weights, decreasing = TRUE), 10L))
  }
  invisible(x)
}

#' Freedman–Lane permutations for a seed-set map
#'
#' @param fit `mppi_fit` object.
#' @param seeds Seed specification (as in [mppi_seedset_slopes()]).
#' @param k Condition/contrast index or name.
#' @param aggregate Aggregation rule passed to [mppi_seedset_map()].
#' @param mode Scale of the interaction matrices (`"normalized"`, `"amplitude"`, `"raw"`).
#' @param target_idx Optional target indices.
#' @param Theta0 Optional precision matrix for `aggregate = "precision"`.
#' @param blksize Block length for permutations (TRs).
#' @param B Number of permutations.
#' @param seed Optional RNG seed (restored afterwards).
#' @param keep_draws Logical; store the seed-set maps for each permutation.
#' @return List with observed map, p-values, z-scores, and permutation metadata.
#' @export
mppi_seedset_permute <- function(fit, seeds, k,
                                 aggregate = c("svd", "mean", "precision"),
                                 mode = c("normalized", "amplitude", "raw"),
                                 target_idx = NULL,
                                 Theta0 = NULL,
                                 blksize = 10L,
                                 B = 499L,
                                 seed = NULL,
                                 keep_draws = FALSE) {
  aggregate <- match.arg(aggregate)
  mode <- match.arg(mode)
  V <- .mppi_seedset_dim(fit)
  labels <- .mppi_seedset_labels(fit)
  seed_idx <- .mppi_as_index(seeds, V, labels, what = "seeds")
  target_idx_spec <- if (is.null(target_idx)) seq_len(V) else target_idx
  target_idx_idx <- .mppi_as_index(target_idx_spec, V, labels, what = "targets")
  map <- mppi_seedset_map(fit, seed_idx, k, aggregate = aggregate,
                          target_idx = target_idx_idx, mode = mode, Theta0 = Theta0)
  target_idx_idx <- if (!is.null(map$target_idx)) map$target_idx else target_idx_idx
  weights <- as.numeric(map$weights)
  seed_names <- names(map$weights)
  target_names <- map$targets
  obs_map <- as.numeric(map$map)
  k_idx <- if (is.character(k)) match(k, fit$names) else as.integer(k)
  if (is.na(k_idx) || k_idx < 1L || k_idx > length(fit$pk)) {
    stop("Invalid regressor index for permutation.", call. = FALSE)
  }
  Tn <- nrow(fit$R)
  if (blksize <= 0L) stop("blksize must be positive", call. = FALSE)
  idx_blocks <- split(seq_len(Tn), ceiling(seq_len(Tn) / blksize))
  outer_centered <- .mppi_row_outer_packed(fit$R)
  col_means <- colMeans(outer_centered)
  outer_centered <- sweep(outer_centered, 2, col_means, "-")
  pk0 <- fit$pk[[k_idx]]
  denom0 <- sum(pk0^2)
  if (denom0 < .Machine$double.eps) {
    stop("Psychological regressor has near-zero variance; cannot run permutations.", call. = FALSE)
  }
  pkc <- pk0 - mean(pk0)
  obs_vec <- drop(crossprod(pkc, outer_centered) / denom0)
  Vc <- ncol(fit$R)
  M_obs <- .mppi_unpack_upper(list(values = obs_vec, dim = Vc))
  diag(M_obs) <- 0
  sigma <- fit$sigma
  if (is.null(sigma)) {
    base_mat <- if (!is.null(fit$R_raw)) fit$R_raw else fit$R
    sigma <- matrixStats::colSds(base_mat)
    sigma[sigma < 1e-12] <- 1e-12
  }
  scale_matrix <- function(M) {
    if (mode == "raw") return(M)
    .mppi_scale_matrix(M, sigma, mode)
  }
  target_idx_idx <- as.integer(target_idx_idx)
  slopes_obs <- .mppi_seedset_slice(scale_matrix(M_obs), fit, seed_idx, target_idx_idx)
  obs_map_check <- drop(weights %*% slopes_obs)
  if (length(obs_map) == length(obs_map_check) &&
      max(abs(obs_map - obs_map_check)) > 1e-6) {
    obs_map <- obs_map_check
  }
  perm_maps <- matrix(NA_real_, nrow = B, ncol = length(obs_map))
  rng_pre <- RNGkind()
  seed_pre <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  if (!is.null(seed)) set.seed(seed)
  rng_used <- RNGkind()
  init_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    do.call(RNGkind, as.list(rng_pre))
    if (is.null(seed_pre)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", seed_pre, envir = .GlobalEnv)
    }
  }, add = TRUE)
  for (b in seq_len(B)) {
    perm_idx <- unlist(sample(idx_blocks))
    perm_vec <- drop(crossprod(pkc, outer_centered[perm_idx, , drop = FALSE]) / denom0)
    M_perm <- .mppi_unpack_upper(list(values = perm_vec, dim = Vc))
    diag(M_perm) <- 0
    slopes_perm <- .mppi_seedset_slice(scale_matrix(M_perm), fit, seed_idx, target_idx_idx)
    perm_maps[b, ] <- drop(weights %*% slopes_perm)
  }
  perm_mean <- colMeans(perm_maps)
  perm_sd <- apply(perm_maps, 2L, stats::sd)
  z <- (obs_map - perm_mean) / (perm_sd + 1e-8)
  abs_obs <- abs(obs_map)
  abs_perm <- abs(perm_maps)
  thresh <- matrix(abs_obs, nrow = B, ncol = length(abs_obs), byrow = TRUE)
  p <- (1 + colSums(abs_perm >= thresh)) / (B + 1)
  result <- list(
    map = setNames(obs_map, target_names),
    weights = setNames(weights, seed_names),
    seeds = seed_idx,
    target_idx = target_idx_idx,
    targets = target_names,
    k = map$k,
    mode = mode,
    aggregate = aggregate,
    statistic = list(z = setNames(z, target_names),
                     p = setNames(p, target_names),
                     mean = setNames(perm_mean, target_names),
                     sd = setNames(perm_sd, target_names)),
    settings = list(B = B, blksize = blksize, seed = seed)
  )
  if (keep_draws) {
    colnames(perm_maps) <- target_names
    result$draws <- perm_maps
  }
  attr(result, "rng") <- list(kind = rng_used, seed = init_seed,
                               final_seed = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL,
                               method = "freedman_lane", seed_arg = seed)
  class(result) <- c("mppi_seedset_perm", "list")
  result
}
#' Group-level inference on seed-set maps
#'
#' @param fits List of `mppi_fit` objects (subjects).
#' @param seeds Seed specification shared across subjects.
#' @param k Condition/contrast index or name.
#' @param aggregate Aggregation rule passed to [mppi_seedset_map()].
#' @param mode Interaction matrix scale (`"normalized"`, `"amplitude"`, `"raw"`).
#' @param design Data frame containing subject-level covariates.
#' @param formula Model formula evaluated against `design`.
#' @param Theta0_list Optional list of precision matrices (one per subject) used when
#'   `aggregate = "precision"`.
#' @param target_idx Optional shared target indices.
#' @return List with coefficient, t- and p-maps plus metadata.
#' @export
mppi_group_seedset <- function(fits, seeds, k,
                               aggregate = c("svd", "mean", "precision"),
                               mode = c("normalized", "amplitude", "raw"),
                               design,
                               formula = ~ 1,
                               Theta0_list = NULL,
                               target_idx = NULL) {
  if (!length(fits)) stop("Provide at least one fit for group analysis.", call. = FALSE)
  aggregate <- match.arg(aggregate)
  mode <- match.arg(mode)
  n_subj <- length(fits)
  if (nrow(design) != n_subj) {
    stop("design must have the same number of rows as fits.", call. = FALSE)
  }
  if (!is.null(Theta0_list) && !length(Theta0_list) %in% c(1L, n_subj)) {
    stop("Theta0_list must have length 1 or match the number of fits.", call. = FALSE)
  }
  maps <- vector("list", n_subj)
  weights <- vector("list", n_subj)
  for (ii in seq_len(n_subj)) {
    Theta_i <- if (is.null(Theta0_list)) NULL else if (length(Theta0_list) == 1L) Theta0_list[[1L]] else Theta0_list[[ii]]
    map_i <- mppi_seedset_map(fits[[ii]], seeds, k, aggregate = aggregate,
                              target_idx = target_idx, mode = mode, Theta0 = Theta_i)
    maps[[ii]] <- map_i$map
    weights[[ii]] <- map_i$weights
  }
  target_names <- names(maps[[1L]])
  if (is.null(target_names)) {
    target_names <- paste0("v", seq_along(maps[[1L]]))
  }
  M <- do.call(rbind, maps)
  colnames(M) <- target_names
  X <- stats::model.matrix(formula, design)
  if (nrow(X) != n_subj) {
    stop("Model matrix rows must match number of subjects.", call. = FALSE)
  }
  qrX <- qr(X)
  betahat <- qr.coef(qrX, M)
  fitted_vals <- X %*% betahat
  resid <- M - fitted_vals
  df <- nrow(X) - qrX$rank
  sigma2 <- colSums(resid^2) / pmax(df, 1L)
  contrast <- rep(0, ncol(X))
  contrast[1L] <- 1
  XtX_inv <- chol2inv(qr.R(qrX))
  c_var <- drop(t(contrast) %*% XtX_inv %*% contrast)
  coef <- drop(contrast %*% betahat)
  se <- sqrt(pmax(sigma2, 0) * c_var)
  tval <- coef / se
  pval <- 2 * stats::pt(-abs(tval), df = df)
  list(
    coef = setNames(coef, target_names),
    t = setNames(tval, target_names),
    p = setNames(pval, target_names),
    df = df,
    targets = target_names,
    aggregate = aggregate,
    mode = mode,
    design = design,
    formula = formula,
    weights = weights,
    maps = maps
  )
}
