# Resting-state multiPPI helpers --------------------------------------------

#' Extract time series, runs, and TR for resting-state analyses
#'
#' This utility accepts either a matrix (`T x V`) or an `fmri_dataset` and
#' returns the pieces required by the rest pipeline.
#'
#' @param data Matrix or fmri_dataset-like object.
#' @param runs Optional run labels (overrides dataset metadata when supplied).
#' @return List with `Y`, `runs`, and `TR` entries.
#' @keywords internal
mppi_rest_extract <- function(data, runs = NULL) {
  if (is.matrix(data)) {
    Y <- data
    TR <- attr(data, "TR")
    runs_ds <- runs %||% attr(data, "runs")
  } else if (inherits(data, "fmri_dataset") && requireNamespace("fmridataset", quietly = TRUE)) {
    Y <- fmridataset::get_data_matrix(data)
    TR <- tryCatch(fmridataset::get_TR(data), error = function(e) NULL)
    runs_ds <- runs %||% tryCatch(fmridataset::get_run_lengths(data), error = function(e) NULL)
    if (length(runs_ds) && length(runs_ds) != nrow(Y)) {
      runs_ds <- rep(seq_along(runs_ds), runs_ds)
    }
    if (!length(runs_ds)) {
      runs_ds <- tryCatch(fmridataset::blockids(data), error = function(e) NULL)
    }
  } else {
    Y <- data$Y %||% data$timeseries %||% stop("Provide a matrix or fmri_dataset-like object.", call. = FALSE)
    TR <- runs <- NULL
    runs_ds <- runs %||% data$runs %||% attr(data, "runs")
  }
  if (!is.null(runs_ds) && length(runs_ds) != nrow(Y)) {
    stop("Runs vector length must match number of time points.", call. = FALSE)
  }
  list(
    Y = as.matrix(Y),
    runs = if (is.null(runs_ds)) NULL else as.integer(runs_ds),
    TR = TR %||% attr(data, "TR")
  )
}

#' Endogenous microstate modulators
#'
#' Builds a list of T-length modulators describing latent states obtained from a
#' PCA basis time-series. Methods include k-means (default) and an optional
#' sticky HMM (via `depmixS4`).
#'
#' @param U Low-rank time-series matrix (`T x r`).
#' @param K Number of states.
#' @param method "kmeans" (default) or "hmm" (requires depmixS4).
#' @param soft Return soft posteriors (TRUE) or hard one-hot sticks (FALSE).
#' @param sticky Self-transition boost for the HMM (ignored for kmeans).
#' @param seed Seed for reproducibility.
#' @param zscore Whether to z-score columns of `U` before clustering.
#' @return Named list of modulators.
mppi_rest_states <- function(U, K = 6L, method = c("kmeans", "hmm"),
                             soft = TRUE, sticky = 2, seed = 1L, zscore = TRUE) {
  method <- match.arg(method)
  set.seed(seed)
  X <- if (zscore) scale(U) else U
  if (method == "kmeans") {
    km <- stats::kmeans(X, centers = K, nstart = 8L)
    sticks <- stats::model.matrix(~ 0 + factor(km$cluster))
    colnames(sticks) <- paste0("state", seq_len(ncol(sticks)))
    if (!soft) return(as.list(as.data.frame(sticks)))
    # simple temporal smoothing for pseudo-soft output
    smooth <- vapply(seq_len(ncol(sticks)), function(j) {
      sm <- stats::filter(sticks[, j], rep(1/3, 3), sides = 2)
      sm[is.na(sm)] <- 0
      as.numeric(sm)
    }, numeric(nrow(sticks)))
    colnames(smooth) <- colnames(sticks)
    return(as.list(as.data.frame(smooth)))
  }
  if (!requireNamespace("depmixS4", quietly = TRUE)) {
    stop("Install depmixS4 to enable HMM-based microstates.", call. = FALSE)
  }
  responses <- lapply(seq_len(ncol(X)), function(j) X[, j] ~ 1)
  fams <- replicate(ncol(X), depmixS4::gaussian(), simplify = FALSE)
  mod <- depmixS4::depmix(response = responses, data = as.data.frame(X),
                          nstates = K, family = fams)
  fit <- suppressWarnings(depmixS4::fit(mod, verbose = FALSE))
  post <- depmixS4::posterior(fit)
  gamma <- as.matrix(post[, paste0("S", seq_len(K)), drop = FALSE])
  colnames(gamma) <- paste0("state", seq_len(K))
  if (soft) return(as.list(as.data.frame(gamma)))
  hard <- stats::model.matrix(~ 0 + factor(max.col(gamma)))
  colnames(hard) <- paste0("state", seq_len(ncol(hard)))
  as.list(as.data.frame(hard))
}

#' Innovation burst modulators
#'
#' Builds sparse stick functions marking high-magnitude innovations derived from
#' AR residuals across basis time-series.
#'
#' @param U Low-rank time-series matrix (`T x r`).
#' @param ar_order AR order used when computing residuals.
#' @param thresholds Quantile thresholds used to define bursts.
#' @return Named list of binary modulators.
mppi_rest_bursts <- function(U, ar_order = 4L, thresholds = c(0.90, 0.95, 0.98)) {
  Tn <- nrow(U); r <- ncol(U)
  eps <- matrix(0, Tn, r)
  for (j in seq_len(r)) {
    ts <- U[, j]
    if (sd(ts) == 0) next
    fit <- try(stats::ar(ts, order.max = ar_order, aic = FALSE, method = "yw"), silent = TRUE)
    if (inherits(fit, "try-error") || !length(fit$ar)) {
      eps[, j] <- ts
    } else {
      res <- ts
      for (lag in seq_along(fit$ar)) {
        res[(lag + 1):Tn] <- res[(lag + 1):Tn] - fit$ar[lag] * ts[seq_len(Tn - lag)]
      }
      eps[, j] <- res
    }
  }
  mag <- sqrt(rowSums(eps^2))
  lapply(thresholds, function(q) as.numeric(mag >= stats::quantile(mag, q, na.rm = TRUE))) |>
    setNames(paste0("burst_q", sprintf("%02d", round(100 * thresholds))))
}

#' Envelope modulators
#'
#' Creates continuous modulators reflecting power or amplitude envelopes of
#' chosen component groups.
#'
#' @param U Low-rank time-series matrix (`T x r`).
#' @param groups Named list of column indices. `NULL` uses all columns.
#' @param transform "power" (default) or "abs" for absolute-value envelopes.
#' @param zscore Whether to z-score the resulting envelope.
#' @return Named list of continuous modulators.
mppi_rest_envelopes <- function(U, groups = list(all = NULL), transform = c("power", "abs"),
                                zscore = TRUE) {
  transform <- match.arg(transform)
  out <- list()
  cols <- seq_len(ncol(U))
  for (nm in names(groups)) {
    idx <- groups[[nm]] %||% cols
    block <- U[, idx, drop = FALSE]
    env <- switch(transform,
                  power = rowSums(block^2),
                  abs = rowSums(abs(block)))
    if (zscore) env <- as.numeric(scale(env))
    out[[paste0("env_", nm)]] <- env
  }
  out
}

#' Assemble endogenous modulators for rest multiPPI
#'
#' @param U Low-rank time-series matrix (`T x r`).
#' @param include Character vector selecting which families to include.
#' @param states,bursts,envelopes Named lists of arguments passed to the
#'   respective helper.
#' @param extra Optional list of additional modulators to append.
#' @return Named list of modulators ready for design construction.
mppi_rest_modulators <- function(U,
                                 include = c("states", "bursts", "envelopes"),
                                 states = list(),
                                 bursts = list(),
                                 envelopes = list(),
                                 extra = NULL) {
  include <- intersect(include, c("states", "bursts", "envelopes"))
  pk <- list()
  if ("states" %in% include) {
    pk <- c(pk, do.call(mppi_rest_states, c(list(U = U), states)))
  }
  if ("bursts" %in% include) {
    pk <- c(pk, do.call(mppi_rest_bursts, c(list(U = U), bursts)))
  }
  if ("envelopes" %in% include) {
    pk <- c(pk, do.call(mppi_rest_envelopes, c(list(U = U), envelopes)))
  }
  if (!is.null(extra)) pk <- c(pk, extra)
  if (!length(pk)) stop("No modulators were created.", call. = FALSE)
  stopifnot(all(vapply(pk, length, 1L) == nrow(U)))
  pk
}

#' Build a gPPI-style design for endogenous modulators
#'
#' Residualizes each modulator against all others, optional nuisance regressors,
#' and an intercept to span the full context space (gPPI logic).
#'
#' @param pk Named list of modulators.
#' @param nuisance Optional nuisance matrix (`T x q`).
#' @param runs Optional run labels used for within-run centering.
#' @param include_intercept Include an intercept column (default TRUE).
#' @param scale_mod Normalize each modulator to unit L2 norm after residualizing.
#' @return List containing design matrix `X`, psychological indices, and the
#'   residualized modulators.
mppi_rest_design <- function(pk, nuisance = NULL, runs = NULL,
                             include_intercept = TRUE, scale_mod = TRUE) {
  pk <- lapply(pk, as.numeric)
  Tn <- length(pk[[1L]])
  stopifnot(all(vapply(pk, length, 1L) == Tn))
  if (!is.null(nuisance)) {
    nuisance <- as.matrix(nuisance)
    if (nrow(nuisance) != Tn) stop("Nuisance rows must match data length.", call. = FALSE)
  }
  N <- length(pk)
  resid_pk <- vector("list", N)
  intercept <- if (include_intercept) matrix(1, Tn, 1) else NULL
  base_Q <- if (is.null(nuisance) && is.null(intercept)) NULL else cbind(nuisance, intercept)
  for (j in seq_len(N)) {
    others <- pk[setdiff(seq_len(N), j)]
    Q <- do.call(cbind, c(others, list(base_Q)))
    if (!is.null(Q)) {
      resid <- .mppi_residualize_vec(pk[[j]], as.matrix(Q))
    } else {
      resid <- pk[[j]]
    }
    if (!is.null(runs)) resid <- .mppi_center_by_run(resid, runs)
    if (scale_mod) {
      sdv <- sqrt(sum(resid^2))
      if (sdv < .Machine$double.eps) {
        warning(sprintf("Modulator '%s' has near-zero variance after residualization.", names(pk)[j]), call. = FALSE)
      }
      resid <- resid / max(sdv, .Machine$double.eps)
    }
    resid_pk[[j]] <- resid
  }
  names(resid_pk) <- names(pk)
  P <- do.call(cbind, resid_pk)
  colnames(P) <- names(resid_pk)
  X <- if (is.null(nuisance)) P else cbind(nuisance, P)
  if (include_intercept) {
    X <- cbind(`(Intercept)` = 1, X)
  }
  psych_idx <- seq(ncol(X) - length(resid_pk) + 1L, ncol(X))
  list(X = X, psych_idx = psych_idx, pk = resid_pk)
}

#' Resting-state multiPPI wrapper
#'
#' @param data Matrix or fmri_dataset providing the time series.
#' @param runs Optional run labels.
#' @param modulators Optional pre-built modulators; if NULL they are derived from
#'   the basis time-series.
#' @param include Which modulator families to construct when `modulators` is NULL.
#' @param basis Basis specification passed to `.mppi_resolve_basis()`.
#' @param modulator_domain Domain for building modulators ("auto", "bold", "neural").
#' @param domain Domain for the final mPPI fit.
#' @param prewhiten Use fmriAR prewhitening (default TRUE when available).
#' @param ar_method,ar_order Controls for fmriAR noise modelling.
#' @param nuisance Optional nuisance regressors aligned with the data.
#' @param attach_pk Attach residualized modulator list to the returned fit.
#' @param ... Passed through to `mppi_fit()` / `mppi_fit_whitened()`.
#' @return `mppi_fit` object.
mppi_rest <- function(data,
                      runs = NULL,
                      modulators = NULL,
                      include = c("states", "bursts", "envelopes"),
                      basis = list(type = "pca", r = 120L),
                      modulator_domain = c("auto", "bold", "neural"),
                      domain = c("neural", "bold", "innovations"),
                      prewhiten = TRUE,
                      ar_method = "ar",
                      ar_order = "auto",
                      nuisance = NULL,
                      attach_pk = TRUE,
                      ...) {
  dom <- match.arg(domain)
  mod_dom <- match.arg(modulator_domain)
  info <- mppi_rest_extract(data, runs = runs)
  Y <- info$Y
  runs_vec <- info$runs
  if (is.null(runs_vec)) runs_vec <- runs %||% attr(data, "runs")
  basis_obj <- .mppi_resolve_basis(Y, basis = basis, dataset = if (is.list(data)) data else NULL)
  U_bold <- if (is.null(basis_obj)) Y else Y %*% basis_obj$V
  U_for_pk <- U_bold
  if (identical(mod_dom, "auto")) {
    mod_dom <- if (!is.null(info$TR) && dom == "neural") "neural" else "bold"
  }
  if (identical(mod_dom, "neural")) {
    hrf <- mppi_default_hrf(tr = info$TR %||% 1, duration = 32)
    U_for_pk <- mppi_deconv(U_bold, TR = info$TR %||% 1, hrf = hrf, lambda = 10)
  }
  pk <- modulators %||% mppi_rest_modulators(U_for_pk, include = include)
  design <- mppi_rest_design(pk, nuisance = nuisance, runs = runs_vec, include_intercept = TRUE)
  X <- design$X
  psych_idx <- design$psych_idx
  fit <- if (prewhiten) {
    mppi_fit_whitened(Y = Y, X = X, runs = runs_vec,
                      psych_idx = psych_idx, ar_method = ar_method, p = ar_order,
                      basis = basis_obj, domain = dom, ...)
  } else {
    mppi_fit(Y = Y, X = X, psych_idx = psych_idx, runs = runs_vec,
             basis = basis_obj, domain = dom, ...)
  }
  if (attach_pk) {
    attr(fit, "pk") <- design$pk
    attr(fit, "runs") <- runs_vec
    attr(fit, "basis_timecourses") <- U_for_pk
    attr(fit, "TR") <- info$TR
  }
  fit
}

#' Subject-level features for resting-state traits
#'
#' Aggregates Gain, Routing, and mode energies across contexts for a list of
#' rest fits.
#'
#' @param fits List of `mppi_fit` objects.
#' @param lags Lags used when computing routing.
#' @param modes Number of communication modes to retain.
#' @param aggregate "mean" to average over contexts per subject, otherwise
#'   returns long-format rows.
#' @return Data frame of subject-level features.
mppi_rest_features <- function(fits, lags = -2:2, modes = 10L, aggregate = c("mean", "none")) {
  aggregate <- match.arg(aggregate)
  rows <- vector("list", length(fits))
  for (i in seq_along(fits)) {
    fit <- fits[[i]]
    ax <- mppi_axes(fit, lags = lags)
    if (!"G" %in% names(ax) && "gain" %in% names(ax)) ax$G <- ax$gain
    if (!"R" %in% names(ax) && "routing" %in% names(ax)) ax$R <- ax$routing
    md <- mppi_modes(fit, r = modes)
    proj <- md$projected
    mag <- vapply(proj, function(M) sqrt(sum(M * M)), 0)
    ax$mag <- mag[match(ax$condition, names(mag))]
    ax$subj <- i
    rows[[i]] <- ax
  }
  feats <- do.call(rbind, rows)
  if (aggregate == "mean") {
    cols <- intersect(c("G", "R", "mag"), names(feats))
    stats::aggregate(feats[, cols, drop = FALSE],
                     by = list(subj = feats$subj),
                     FUN = mean)
  } else {
    feats
  }
}

#' Simple nested CV ridge for trait prediction
#'
#' @param features Output from `mppi_rest_features(aggregate = "none")`.
#' @param y Named numeric trait vector (names correspond to subject indices).
#' @param covariates Optional covariate matrix keyed by subject name.
#' @param k Number of outer CV folds.
#' @param seed RNG seed.
#' @param lambda Ridge penalty.
#' @return List with predictions, observed values, and CV R^2.
mppi_rest_trait_cv <- function(features, y, covariates = NULL,
                               k = 5L, seed = 1L, lambda = 1e-2) {
  stopifnot(is.numeric(y))
  subj_ids <- sort(unique(features$subj))
  if (is.null(names(y))) names(y) <- subj_ids
  if (length(y) != length(subj_ids)) {
    stop("Trait vector length must equal number of subjects in features.", call. = FALSE)
  }
  y <- y[order(as.numeric(names(y)))]
  set.seed(seed)
  folds <- split(sample(subj_ids), rep_len(seq_len(k), length(subj_ids)))
  Xwide <- stats::reshape(features,
                          idvar = "subj", timevar = "condition",
                          direction = "wide")
  rownames(Xwide) <- Xwide$subj
  pred_cols <- grepl("^(G|R|mag)\\.", colnames(Xwide))
  X <- as.matrix(Xwide[, pred_cols, drop = FALSE])
  y_aligned <- y[match(rownames(X), names(y))]
  if (!is.null(covariates)) {
    covariates <- covariates[match(rownames(X), rownames(covariates)), , drop = FALSE]
    X <- cbind(X, covariates)
  }
  preds <- rep(NA_real_, length(y_aligned)); names(preds) <- rownames(X)
  for (fold in seq_len(k)) {
    test_ids <- folds[[fold]]
    train_ids <- setdiff(subj_ids, test_ids)
    Xtr <- X[as.character(train_ids), , drop = FALSE]
    ytr <- y_aligned[match(as.character(train_ids), names(y_aligned))]
    XtX <- crossprod(Xtr)
    beta <- solve(XtX + diag(lambda, ncol(Xtr)), crossprod(Xtr, ytr))
    preds[as.character(test_ids)] <- X[as.character(test_ids), , drop = FALSE] %*% beta
  }
  cv_r2 <- 1 - sum((y_aligned - preds)^2, na.rm = TRUE) /
    sum((y_aligned - mean(y_aligned, na.rm = TRUE))^2, na.rm = TRUE)
  list(pred = preds, observed = y_aligned, cv_r2 = cv_r2)
}
