#' Fit conventional PPI or generalized PPI
#'
#' Fits target BOLD against a completed [ppi_design()]. All psychological,
#' physiological, interaction, and nuisance columns are fitted jointly.
#' fmriAR estimates noise from full-model residuals, then applies the same
#' operator to the completed design and responses. Target chunks share a
#' factorization when they share a noise operator. No target is deconvolved.
#'
#' @param Y Numeric time by target matrix, or an `fmri_frame`.
#' @param design A `ppi_design` object.
#' @param noise `"fmriAR"` (default), `"iid"`, or a fitted `fmriAR_plan`.
#' @param noise_control Named arguments for `fmriAR::fit_noise`. Global or
#'   run pooling assumes a shared noise operator. For grouped operators use
#'   `list(pooling = 'parcel', parcels = labels)`, one nonmissing label per
#'   target. `parcels = seq_len(ncol(Y))` estimates each target separately and
#'   can be computationally expensive. Shared operators assume sufficiently
#'   homogeneous target noise; heterogeneous pooling can inflate error rates.
#' @param chunk_size Positive target chunk width.
#' @param tol Relative singular-value tolerance for design rank.
#' @return A `ppi_fit` with coefficients, standard errors, variance, degrees
#'   of freedom, estimability, design metadata and noise plan. Inference with
#'   estimated noise is approximate and conditional on the supplied seed/design.
#' @details Censored observations are supported with IID errors. Temporal
#'   whitening with censor gaps requires a marginal observed-data noise
#'   operator; unsupported gap handling is rejected rather than silently
#'   interpreting gaps as independent noise segments. Individual aliased
#'   coefficients are `NA`; estimable contrasts remain available.
#'
#'   Targets whose full-model residuals are numerically zero, such as a target
#'   identical to the seed, carry no information about temporal noise. They
#'   are flagged in `degenerate`, excluded from fmriAR noise estimation with a
#'   warning naming them, and whitened with the operator estimated from the
#'   remaining targets of their group; a group with no other targets uses IID
#'   errors. A grouped noise fit that fails reports the group and its targets.
#' @export
ppi_fit <- function(Y, design, noise = "fmriAR", noise_control = list(),
                    chunk_size = 2048L, tol = 1e-7) {
  if (!inherits(design, "ppi_design")) stop("design must be a ppi_design")
  if (inherits(Y, "fmri_frame")) {
    if (!requireNamespace("fmridataset", quietly = TRUE)) stop("fmridataset is required for fmri_frame input")
    schema <- fmridataset::temporal_schema(Y)
    if (!isTRUE(schema$contiguous) || length(schema$run_ids) != nrow(design$X)) stop("frame must preserve the contiguous full acquisition grid")
    if (length(schema$TR) == 0L || any(!is.finite(schema$TR)) || any(abs(schema$TR - design$TR) > 1e-8)) stop("frame TR does not match design")
    if (!identical(as.character(schema$run_ids), as.character(design$runs))) stop("frame run labels do not match design")
    offset <- design$provenance$sampling_offset
    if (!is.null(schema$start_time) && (is.null(offset) || any(!is.finite(schema$start_time)) || any(abs(schema$start_time - offset) > 1e-8))) stop("frame sampling offset does not match design")
    censor <- schema$censor
    if (is.null(censor)) censor <- rep(FALSE, length(schema$run_ids))
    if (!is.logical(censor) || length(censor) != length(design$keep) || anyNA(censor) || !identical(unname(!censor), unname(design$keep))) stop("frame censor mask does not match design")
    Y <- as.matrix(fmridataset::collect_assay(Y))
  }
  X <- design$X; keep <- design$keep
  if (!is.matrix(Y) || !is.numeric(Y) || nrow(Y) != nrow(X) || ncol(Y) < 1L)
    stop("Y must be a numeric time by target matrix matching design")
  if (any(!is.finite(Y[keep, , drop = FALSE]))) stop("observed target values must be finite; use a common design censor mask")
  if (!is.logical(keep) || length(keep) != nrow(X) || anyNA(keep) || any(!is.finite(X))) stop("invalid design observation mask or matrix")
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L || !is.finite(chunk_size) || chunk_size < 1 || chunk_size != floor(chunk_size)) stop("chunk_size must be a positive integer")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0 || tol >= 1) stop("tol must lie between zero and one")
  fixed <- inherits(noise, "fmriAR_plan")
  if (!fixed && (!is.character(noise) || length(noise) != 1L || !noise %in% c("fmriAR", "iid"))) stop("noise must be 'fmriAR', 'iid', or an fmriAR_plan")
  use_ar <- fixed || identical(noise, "fmriAR")
  if (!is.list(noise_control) || (length(noise_control) && (is.null(names(noise_control)) || any(!nzchar(names(noise_control))) || anyDuplicated(names(noise_control))))) stop("noise_control must be a named list")
  if (length(noise_control) && (!use_ar || fixed)) stop("noise_control applies only when estimating fmriAR noise")
  if (any(names(noise_control) %in% c("resid", "Y", "X", "design", "runs", "censor"))) stop("noise_control must not override data, design, runs, or censor")
  grouped <- !fixed && identical(noise, "fmriAR") && identical(noise_control$pooling, "parcel")
  if (grouped) {
    parcels <- noise_control$parcels
    if ((!is.numeric(parcels) && !is.character(parcels) && !is.factor(parcels)) || length(parcels) != ncol(Y) || anyNA(parcels) || (is.numeric(parcels) && any(!is.finite(parcels))) || any(!nzchar(as.character(parcels)))) stop("parcel pooling requires one valid nonmissing label per target")
    labels <- as.character(parcels); groups <- unique(labels)
    control <- noise_control; control$parcels <- NULL; control$pooling <- "global"
    target_label <- function(cols) .ppi_target_label(colnames(Y), cols)
    degenerate <- list()
    subfits <- lapply(groups, function(g) {
      cols <- which(labels == g)
      withCallingHandlers(
        tryCatch(ppi_fit(Y[, cols, drop = FALSE], design, noise = "fmriAR", noise_control = control, chunk_size = chunk_size, tol = tol),
                 error = function(e) stop(sprintf("noise estimation failed for parcel '%s' (targets %s): %s", g, target_label(cols), conditionMessage(e)), call. = FALSE)),
        ppi_degenerate_targets = function(w) {
          degenerate[[g]] <<- cols[w$targets]
          invokeRestart("muffleWarning")
        })
    })
    names(subfits) <- groups
    if (length(degenerate)) {
      white <- names(degenerate)[vapply(names(degenerate), function(g) subfits[[g]]$noise == "iid", logical(1))]
      warning(sprintf("%d target(s) have numerically zero full-model residuals (e.g. a target identical to the seed) and were excluded from noise estimation: %s.%s",
                      length(unlist(degenerate)),
                      paste(sprintf("%s (parcel '%s')", vapply(degenerate, target_label, character(1)), names(degenerate)), collapse = ", "),
                      if (length(white)) paste0(" Parcels with no other targets use IID errors: ", paste0("'", white, "'", collapse = ", "), ".") else ""),
              call. = FALSE)
    }
    result <- subfits[[1L]]
    result$coefficients <- result$standard_errors <- matrix(NA_real_, nrow = ncol(X), ncol = ncol(Y), dimnames = list(colnames(X), colnames(Y)))
    result$beta_solution <- result$coefficients
    result$sigma2 <- setNames(numeric(ncol(Y)), colnames(Y))
    result$degenerate <- setNames(logical(ncol(Y)), colnames(Y))
    for (i in seq_along(groups)) {
      cols <- which(labels == groups[i]); f <- subfits[[i]]
      result$coefficients[, cols] <- f$coefficients
      result$standard_errors[, cols] <- f$standard_errors
      result$beta_solution[, cols] <- f$beta_solution
      result$sigma2[cols] <- f$sigma2
      result$degenerate[cols] <- f$degenerate
    }
    result$groups <- labels; result$group_fits <- subfits
    result$noise_plan <- lapply(subfits, `[[`, "noise_plan")
    result$covariance_unscaled <- result$row_space <- NULL
    result$df <- vapply(labels, function(g) subfits[[g]]$df, numeric(1))
    result$rank <- vapply(labels, function(g) subfits[[g]]$rank, numeric(1))
    result$estimable <- do.call(cbind, lapply(labels, function(g) subfits[[g]]$estimable))
    colnames(result$estimable) <- colnames(Y)
    result$inference <- "conditional approximate grouped GLS"
    return(result)
  }
  if (!is.null(noise_control$parcels)) stop("parcels requires pooling='parcel'")
  if (use_ar && any(!keep)) stop("fmriAR with censor gaps requires observed-data whitening; use noise='iid' only if IID errors are appropriate")
  observed <- which(keep)
  Xobs <- X[observed, , drop = FALSE]
  initial <- .ppi_factor(Xobs, tol)
  if (length(observed) <= initial$rank) stop("no residual degrees of freedom")
  plan <- NULL
  degenerate <- logical(ncol(Y))
  if (use_ar) {
    if (fixed) plan <- noise else {
      if (!is.null(noise_control$pooling) && !noise_control$pooling %in% c("global", "run")) stop("only shared global/run noise operators are supported")
      # This is the full-model residual, including interaction columns.
      residual <- Y - initial$U %*% crossprod(initial$U, Y)
      # Targets the design explains exactly (e.g. identical to the seed) leave
      # rounding noise that carries no information about temporal noise.
      degenerate <- sqrt(colSums(residual^2)) <= sqrt(.Machine$double.eps) * sqrt(colSums(Y^2))
      if (any(degenerate)) {
        warning(structure(class = c("ppi_degenerate_targets", "warning", "condition"), list(
          message = sprintf("%d target(s) have numerically zero full-model residuals (e.g. a target identical to the seed) and were excluded from noise estimation: %s.%s",
                            sum(degenerate), .ppi_target_label(colnames(Y), which(degenerate)),
                            if (all(degenerate)) " No other targets remain; using IID errors." else ""),
          call = NULL, targets = which(degenerate))))
      }
      if (all(degenerate)) use_ar <- FALSE else {
        args <- c(list(resid = residual[, !degenerate, drop = FALSE], runs = design$runs, design = X), noise_control)
        plan <- do.call(fmriAR::fit_noise, args)
      }
      rm(residual)
    }
  }
  if (use_ar) {
    if (!plan$pooling %in% c("global", "run")) stop("only shared global/run noise operators are supported")
    if (!is.null(plan$n_time) && plan$n_time != nrow(X)) stop("noise plan has incompatible time length")
    if (!is.null(plan$runs) && !identical(as.character(plan$runs), as.character(design$runs))) stop("noise plan runs do not match design")
    if (!is.null(plan$censor) && any(plan$censor)) stop("noise plan must not contain censor gaps")
    Xw <- fmriAR::whiten_apply(plan, X = X, Y = matrix(0, nrow(X), 1L),
                               runs = design$runs, parallel = FALSE)$X
  } else Xw <- Xobs
  fact <- .ppi_factor(Xw, tol)
  df <- nrow(Xw) - fact$rank
  if (df <= 0L) stop("no residual degrees of freedom after whitening")
  q <- ncol(X); v <- ncol(Y)
  solution <- matrix(0, q, v, dimnames = list(colnames(X), colnames(Y)))
  sigma2 <- numeric(v)
  for (start in seq.int(1L, v, by = chunk_size)) {
    cols <- seq.int(start, min(v, start + chunk_size - 1L))
    yy <- Y[observed, cols, drop = FALSE]
    if (use_ar) yy <- fmriAR::whiten_apply(plan, X = X, Y = yy, runs = design$runs, parallel = FALSE)$Y
    uy <- crossprod(fact$U, yy)
    solution[, cols] <- fact$inverse %*% uy
    residual <- yy - fact$U %*% uy
    sigma2[cols] <- colSums(residual^2) / df
  }
  coefficients <- solution
  coefficients[!fact$estimable, ] <- NA_real_
  se <- sqrt(outer(diag(fact$covariance), sigma2))
  dimnames(se) <- dimnames(solution); se[!fact$estimable, ] <- NA_real_
  names(sigma2) <- colnames(Y)
  names(degenerate) <- colnames(Y)
  structure(list(coefficients = coefficients, standard_errors = se,
                 sigma2 = sigma2, degenerate = degenerate, df = df, rank = fact$rank,
                 estimable = stats::setNames(fact$estimable, colnames(X)),
                 design = design, noise_plan = plan,
                 noise = if (use_ar) "fmriAR" else "iid",
                 noise_estimated = use_ar && !fixed,
                 inference = if (use_ar && !fixed) "conditional approximate GLS" else "conditional fixed-covariance GLM",
                 estimand = "conventional PPI", schema_version = 1L,
                 beta_solution = solution, covariance_unscaled = fact$covariance,
                 row_space = fact$V, column_scale = fact$column_scale,
                 tolerance = tol), class = "ppi_fit")
}

.ppi_target_label <- function(names, cols) {
  shown <- if (is.null(names)) as.character(cols) else names[cols]
  if (length(shown) > 10L) shown <- c(shown[1:10], sprintf("... (%d more)", length(shown) - 10L))
  paste(shown, collapse = ", ")
}

.ppi_factor <- function(X, tol) {
  column_scale <- apply(abs(X), 2L, max)
  column_scale[column_scale == 0] <- 1
  scaled <- sweep(X, 2L, column_scale, "/")
  s <- svd(scaled, nu = min(dim(X)), nv = min(dim(X)))
  rank <- if (length(s$d) && s$d[1L] > 0) sum(s$d > tol * s$d[1L]) else 0L
  if (!rank) stop("design has zero rank")
  V <- s$v[, seq_len(rank), drop = FALSE]
  U <- s$u[, seq_len(rank), drop = FALSE]
  inverse <- sweep(sweep(V, 2L, s$d[seq_len(rank)], "/"), 1L, column_scale, "/")
  projection_error <- sqrt(rowSums((diag(ncol(X)) - tcrossprod(V))^2))
  list(rank = rank, V = V, U = U, inverse = inverse,
       covariance = tcrossprod(inverse), column_scale = column_scale,
       estimable = projection_error <= tol)
}

#' Conventional PPI contrasts
#'
#' @param fit A `ppi_fit`.
#' @param contrast Named numeric vector or matrix with one contrast per row.
#'   Names may be complete design-column names, or condition names. Condition
#'   weights select interaction columns and average equally across the runs in
#'   which that column is not identically zero; runs where the condition is
#'   absent, or constant so that centering removes it, are skipped.
#'   Unspecified weights are zero. Mixing the two naming conventions is an error.
#' @param test `"t"` for individual rows or `"F"` for a joint test.
#' @return List of estimates, standard errors (t tests), statistics, p-values,
#'   degrees of freedom, estimability, and the resolved contrast matrix.
#' @details Nonestimable contrasts return `NA`, never a contrast of whichever
#'   coefficients happened to survive a pivot. P-values are uncorrected and
#'   inherit the conditional assumptions and noise qualification of the fit.
#' @export
ppi_contrast <- function(fit, contrast, test = c("t", "F")) {
  if (!inherits(fit, "ppi_fit")) stop("fit must be a ppi_fit")
  test <- match.arg(test)
  if (!is.null(fit$group_fits)) {
    parts <- lapply(fit$group_fits, ppi_contrast, contrast = contrast, test = test)
    result <- parts[[1L]]
    n <- length(fit$groups)
    matrix_fields <- if (test == "t") c("estimate", "standard_error", "statistic", "p_value") else character()
    if (test == "t") {
      for (field in matrix_fields) {
        value <- matrix(NA_real_, nrow(parts[[1L]][[field]]), n, dimnames = list(rownames(parts[[1L]][[field]]), colnames(fit$coefficients)))
        for (g in names(parts)) value[, fit$groups == g] <- parts[[g]][[field]]
        result[[field]] <- value
      }
      result$df <- vapply(fit$groups, function(g) parts[[g]]$df, numeric(1))
      result$estimable <- do.call(cbind, lapply(fit$groups, function(g) parts[[g]]$estimable))
    } else {
      for (field in c("statistic", "p_value")) {
        value <- rep(NA_real_, n)
        for (g in names(parts)) value[fit$groups == g] <- parts[[g]][[field]]
        names(value) <- colnames(fit$coefficients)
        result[[field]] <- value
      }
      result$numerator_df <- vapply(fit$groups, function(g) parts[[g]]$numerator_df, numeric(1))
      result$denominator_df <- vapply(fit$groups, function(g) parts[[g]]$denominator_df, numeric(1))
      result$estimable <- vapply(fit$groups, function(g) parts[[g]]$estimable, logical(1))
    }
    result$inference <- fit$inference
    return(result)
  }
  if (is.numeric(contrast) && is.null(dim(contrast))) {
    contrast <- matrix(contrast, nrow = 1L, dimnames = list(NULL, names(contrast)))
  }
  if (!is.matrix(contrast) || !is.numeric(contrast) || !nrow(contrast) || !ncol(contrast) || any(!is.finite(contrast))) stop("contrast must be a finite named numeric vector or matrix")
  nm <- colnames(contrast)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) stop("contrast columns must have unique nonempty names")
  X <- fit$design$X; terms <- fit$design$terms
  L <- matrix(0, nrow(contrast), ncol(X), dimnames = list(rownames(contrast), colnames(X)))
  conditions <- unique(terms$condition[terms$type == "interaction"])
  if (all(nm %in% colnames(X))) L[, nm] <- contrast
  else if (all(nm %in% conditions)) {
    observed <- X[fit$design$keep, , drop = FALSE]
    for (j in seq_along(nm)) {
      idx <- which(terms$type == "interaction" & terms$condition == nm[j])
      # Skip runs where the condition has no interaction term at all.
      present <- idx[colSums(abs(observed[, idx, drop = FALSE])) > 0]
      if (length(present)) idx <- present
      L[, idx] <- contrast[, j] / length(idx)
    }
  } else stop("contrast names must all be design columns or all be interaction conditions")
  Lscaled <- sweep(L, 2L, fit$column_scale, "/")
  row_scale <- apply(abs(Lscaled), 1L, max)
  row_scale[row_scale == 0] <- 1
  Lunit <- sweep(Lscaled, 1L, row_scale, "/")
  projection <- (Lunit %*% fit$row_space) %*% t(fit$row_space)
  estimable <- sqrt(rowSums((Lunit - projection)^2)) <= fit$tolerance * sqrt(rowSums(Lunit^2))
  estimate <- L %*% fit$beta_solution
  covariance <- L %*% fit$covariance_unscaled %*% t(L)
  if (test == "t") {
    se <- sqrt(outer(pmax(0, diag(covariance)), fit$sigma2))
    statistic <- estimate / se
    valid <- estimable & diag(covariance) > 0
    estimate[!valid, ] <- se[!valid, ] <- statistic[!valid, ] <- NA_real_
    p <- 2 * stats::pt(abs(statistic), df = fit$df, lower.tail = FALSE)
    return(list(estimate = estimate, standard_error = se, statistic = statistic,
                p_value = p, df = fit$df, estimable = estimable, contrast = L,
                correction = "none", inference = fit$inference))
  }
  # Row scaling must not change the joint hypothesis or its numerator rank.
  fscale <- apply(abs(L), 1L, max); fscale[fscale == 0] <- 1
  Ltest <- sweep(L, 1L, fscale, "/")
  estimate <- Ltest %*% fit$beta_solution
  covariance <- Ltest %*% fit$covariance_unscaled %*% t(Ltest)
  sd <- sqrt(pmax(0, diag(covariance))); sd[sd == 0] <- 1
  covariance <- covariance / outer(sd, sd)
  estimate <- sweep(estimate, 1L, sd, "/")
  ee <- eigen(covariance, symmetric = TRUE)
  good <- ee$values > fit$tolerance * max(ee$values, 0)
  numerator_df <- sum(good)
  statistic <- rep(NA_real_, ncol(fit$beta_solution))
  if (all(estimable) && numerator_df > 0) {
    z <- crossprod(ee$vectors[, good, drop = FALSE], estimate)
    statistic <- colSums(sweep(z^2, 1L, ee$values[good], "/")) / numerator_df / fit$sigma2
  }
  list(statistic = statistic, p_value = if (numerator_df) stats::pf(statistic, numerator_df, fit$df, lower.tail = FALSE) else rep(NA_real_, length(statistic)),
       numerator_df = numerator_df, denominator_df = fit$df,
       estimable = all(estimable), contrast = L, correction = "none", inference = fit$inference)
}
