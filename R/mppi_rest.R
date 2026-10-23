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
  } else if (inherits(data, "fmri_frame") && requireNamespace("fmridataset", quietly = TRUE)) {
    Y <- fmridataset::collect_assay(data)
    schema <- fmridataset::temporal_schema(data)
    TR <- unique(as.numeric(schema$TR))
    if (length(TR) != 1L || !is.finite(TR) || !isTRUE(schema$contiguous)) stop("Resting-state frame requires common regular sampling.", call. = FALSE)
    if (!is.null(schema$censor) && any(schema$censor)) stop("Censored resting-state frames require a supported temporal missing-data model.", call. = FALSE)
    runs_ds <- runs %||% schema$block_ids
  } else {
    Y <- data$Y %||% data$timeseries %||% stop("Provide a matrix or fmri_dataset-like object.", call. = FALSE)
    TR <- data$TR %||% attr(data, "TR")
    runs_ds <- runs %||% data$runs %||% attr(data, "runs")
  }
  if (!is.null(runs_ds) && length(runs_ds) != nrow(Y)) {
    stop("Runs vector length must match number of time points.", call. = FALSE)
  }
  list(
    Y = as.matrix(Y),
    runs = runs_ds,
    TR = TR %||% attr(data, "TR")
  )
}

#' Endogenous microstate modulators
#'
#' Builds a list of T-length modulators describing latent states obtained from a
#' PCA basis time-series. Methods include k-means (default) and an optional
#' sticky HMM (via `depmixS4`).
#'
#' Occupancies sum to one at every time point, so a full set spans the
#' intercept and cannot be fitted alongside it. States are labelled by
#' decreasing occupancy and treatment-coded: `state1`, the most occupied, is the
#' reference and is omitted, so the coupling slope of `statek` contrasts state
#' k with state 1.
#'
#' @param U Low-rank time-series matrix (`T x r`).
#' @param K Number of states (at least 2).
#' @param method "kmeans" (default) or "hmm" (requires depmixS4).
#' @param soft Return soft occupancies (TRUE) or hard one-hot sticks (FALSE).
#'   For k-means, soft occupancies are a 3-TR moving average of the sticks.
#' @param sticky Self-transition boost for the HMM (ignored for kmeans).
#' @param seed Seed for reproducibility; the caller's RNG state is restored
#'   afterwards. `NULL` draws from (and advances) the caller's stream.
#' @param zscore Whether to z-score columns of `U` before clustering.
#' @return Named list of `K - 1` occupancy modulators (`state2`, ..., `stateK`).
mppi_rest_states <- function(U, K = 6L, method = c("kmeans", "hmm"),
                             soft = TRUE, sticky = 2, seed = 1L, zscore = TRUE) {
  method <- match.arg(method)
  if (K < 2) stop("K must be at least 2.", call. = FALSE)
  restore_rng <- .mppi_seed_rng(seed); on.exit(restore_rng(), add = TRUE)
  X <- if (zscore) scale(U) else U
  if (method == "kmeans") {
    km <- stats::kmeans(X, centers = K, nstart = 8L)
    occ <- stats::model.matrix(~ 0 + factor(km$cluster))
    if (soft) {
      # edge-replicated moving average, so rows still sum to one
      n <- nrow(occ)
      occ <- (occ[c(1L, seq_len(n - 1L)), , drop = FALSE] + occ +
                occ[c(seq_len(n)[-1L], n), , drop = FALSE]) / 3
    }
  } else {
    if (!requireNamespace("depmixS4", quietly = TRUE)) {
      stop("Install depmixS4 to enable HMM-based microstates.", call. = FALSE)
    }
    responses <- lapply(seq_len(ncol(X)), function(j) X[, j] ~ 1)
    fams <- replicate(ncol(X), stats::gaussian(), simplify = FALSE)
    mod <- depmixS4::depmix(response = responses, data = as.data.frame(X),
                            nstates = K, family = fams)
    fit <- suppressWarnings(depmixS4::fit(mod, verbose = FALSE))
    post <- depmixS4::posterior(fit)
    gamma <- as.matrix(post[, paste0("S", seq_len(K)), drop = FALSE])
    occ <- if (soft) gamma else stats::model.matrix(~ 0 + factor(max.col(gamma)))
  }
  occ <- unname(as.matrix(occ))[, order(-colSums(occ)), drop = FALSE]
  if (ncol(occ) < 2L) stop("Fewer than two states are occupied.", call. = FALSE)
  colnames(occ) <- paste0("state", seq_len(ncol(occ)))
  as.list(as.data.frame(occ[, -1L, drop = FALSE]))
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
#' Centers each modulator on the intercept and optional nuisance regressors
#' (and within runs when `runs` is given). Modulators are not residualized on
#' each other here: the fit residualizes each one on all other design columns
#' (gPPI logic), so its slope reflects its unique variance. Modulators collinear
#' with each other or with the intercept/nuisance, such as a full set of state
#' occupancies, are an error; drop a reference modulator.
#'
#' @param pk Named list of modulators.
#' @param nuisance Optional nuisance matrix (`T x q`).
#' @param runs Optional run labels used for within-run centering.
#' @param include_intercept Include an intercept column (default TRUE).
#' @param scale_mod Normalize each modulator to unit L2 norm after centering.
#' @return List containing design matrix `X`, psychological indices, and the
#'   centered modulators as they enter `X`.
mppi_rest_design <- function(pk, nuisance = NULL, runs = NULL,
                             include_intercept = TRUE, scale_mod = TRUE) {
  pk <- lapply(pk, as.numeric)
  Tn <- length(pk[[1L]])
  stopifnot(all(vapply(pk, length, 1L) == Tn))
  if (!is.null(nuisance)) {
    nuisance <- as.matrix(nuisance)
    if (nrow(nuisance) != Tn) stop("Nuisance rows must match data length.", call. = FALSE)
  }
  if (any(vapply(pk, anyNA, TRUE))) stop("Modulators must not contain NA.", call. = FALSE)
  N <- length(pk)
  resid_pk <- vector("list", N)
  intercept <- if (include_intercept) matrix(1, Tn, 1) else NULL
  base_Q <- if (is.null(nuisance) && is.null(intercept)) NULL else cbind(nuisance, intercept)
  for (j in seq_len(N)) {
    resid <- if (is.null(base_Q)) pk[[j]] else .mppi_residualize_vec(pk[[j]], base_Q)
    if (!is.null(runs)) resid <- .mppi_center_by_run(resid, runs)
    sdv <- sqrt(sum(resid^2))
    if (sdv <= 1e-8 * sqrt(sum(pk[[j]]^2))) {
      stop(sprintf("Modulator '%s' has no variance beyond the intercept and nuisance.", names(pk)[j]), call. = FALSE)
    }
    if (scale_mod) resid <- resid / sdv
    resid_pk[[j]] <- resid
  }
  names(resid_pk) <- names(pk)
  P <- do.call(cbind, resid_pk)
  colnames(P) <- names(resid_pk)
  X <- if (is.null(nuisance)) P else cbind(nuisance, P)
  if (include_intercept) {
    X <- cbind(`(Intercept)` = 1, X)
  }
  if (qr(X)$rank < N + if (is.null(base_Q)) 0L else qr(base_Q)$rank) {
    stop("Modulators are collinear with each other or with the intercept/nuisance (e.g. a full set of state occupancies); drop a reference modulator.", call. = FALSE)
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
#'   Neural modulators deconvolve each run separately.
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
    Tn <- nrow(U_bold)
    seg <- if (is.null(runs_vec)) rep(1L, Tn) else cumsum(c(TRUE, runs_vec[-1L] != runs_vec[-Tn]))
    U_for_pk <- do.call(rbind, lapply(split(seq_len(Tn), seg), function(i)
      mppi_deconv(U_bold[i, , drop = FALSE], TR = info$TR %||% 1, hrf = hrf, lambda = 10)))
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
#' @param fits List of `mppi_fit` objects; its names, when present, label the
#'   subjects (otherwise `1..n`).
#' @param lags Lags used when computing routing.
#' @param modes Number of communication modes to retain.
#' @param aggregate "mean" to average over contexts per subject, otherwise
#'   returns long-format rows.
#' @return Data frame of subject-level features.
mppi_rest_features <- function(fits, lags = -2:2, modes = 10L, aggregate = c("mean", "none")) {
  aggregate <- match.arg(aggregate)
  ids <- names(fits)
  if (is.null(ids)) {
    ids <- seq_along(fits)
  } else if (anyNA(ids) || !all(nzchar(ids)) || anyDuplicated(ids)) {
    stop("Names of 'fits' must be unique and non-empty (or absent).", call. = FALSE)
  }
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
    ax$subj <- ids[[i]]
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

#' Cross-validated ridge regression for trait prediction
#'
#' Single-level k-fold cross-validation with a fixed ridge penalty (`lambda` is
#' not tuned; there is no inner loop). In each fold, features are centered and
#' scaled and the trait is centered (an unpenalized intercept) using
#' training-fold statistics only, which are then applied to the held-out
#' subjects.
#'
#' @param features Output from `mppi_rest_features(aggregate = "none")`.
#' @param y Numeric trait vector. When named, names are matched to
#'   `features$subj`; otherwise values follow `sort(unique(features$subj))`.
#' @param covariates Optional numeric covariate matrix, matched to subjects by
#'   row names when present, otherwise by position.
#' @param k Number of CV folds.
#' @param seed RNG seed for the fold split; the caller's RNG state is restored
#'   afterwards. `NULL` draws from (and advances) the caller's stream.
#' @param lambda Ridge penalty on the standardized features.
#' @return List with predictions, observed values (both named by subject), and
#'   CV R^2.
mppi_rest_trait_cv <- function(features, y, covariates = NULL,
                               k = 5L, seed = 1L, lambda = 1e-2) {
  stopifnot(is.numeric(y))
  subj_ids <- as.character(sort(unique(features$subj)))
  n <- length(subj_ids)
  # rows of a subject-keyed object: by names when present, else by position
  align <- function(ids, len, what) {
    if (is.null(ids)) {
      if (len != n) stop(sprintf("%s must have one entry per subject in features.", what), call. = FALSE)
      return(seq_len(n))
    }
    m <- match(subj_ids, ids)
    if (len != n || anyNA(m) || anyDuplicated(ids)) {
      stop(sprintf("%s names must match features$subj one-to-one.", what), call. = FALSE)
    }
    m
  }
  y <- setNames(unname(y)[align(names(y), length(y), "y")], subj_ids)
  restore_rng <- .mppi_seed_rng(seed); on.exit(restore_rng(), add = TRUE)
  folds <- split(subj_ids[sample.int(n)], rep_len(seq_len(k), n))
  Xwide <- stats::reshape(features,
                          idvar = "subj", timevar = "condition",
                          direction = "wide")
  rownames(Xwide) <- as.character(Xwide$subj)
  pred_cols <- grepl("^(G|R|mag)\\.", colnames(Xwide))
  X <- as.matrix(Xwide[subj_ids, pred_cols, drop = FALSE])
  if (!is.null(covariates)) {
    covariates <- as.matrix(covariates)
    stopifnot(is.numeric(covariates))
    X <- cbind(X, covariates[align(rownames(covariates), nrow(covariates), "covariates"), , drop = FALSE])
  }
  preds <- setNames(rep(NA_real_, n), subj_ids)
  for (test_ids in folds) {
    train_ids <- setdiff(subj_ids, test_ids)
    Xtr <- X[train_ids, , drop = FALSE]
    mu <- colMeans(Xtr)
    s <- apply(Xtr, 2, stats::sd)
    s[!is.finite(s) | s == 0] <- 1
    Ztr <- sweep(sweep(Xtr, 2, mu), 2, s, "/")
    ym <- mean(y[train_ids])
    beta <- solve(crossprod(Ztr) + diag(lambda, ncol(Ztr)), crossprod(Ztr, y[train_ids] - ym))
    Zte <- sweep(sweep(X[test_ids, , drop = FALSE], 2, mu), 2, s, "/")
    preds[test_ids] <- ym + Zte %*% beta
  }
  cv_r2 <- 1 - sum((y - preds)^2, na.rm = TRUE) /
    sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  list(pred = preds, observed = y, cv_r2 = cv_r2)
}
