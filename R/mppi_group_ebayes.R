# Group-level empirical Bayes (limma) -------------------------------------

#' Moderated one-sample tests on vectorized ΔΣ across subjects
#' @param Delta_list list of VxV matrices (same V), one per subject
#' @param use_limma logical; if TRUE and limma present, use eBayes
#' @return list with t, p, and (if limma) moderated stats
mppi_group_ebayes <- function(Delta_list, use_limma = TRUE) {
  stopifnot(is.list(Delta_list), length(Delta_list) >= 2)
  V <- nrow(Delta_list[[1]])
  U <- upper.tri(matrix(0, V, V), diag = FALSE)
  X <- do.call(rbind, lapply(Delta_list, function(D) as.numeric(D[U])))
  # One-sample vs zero
  if (use_limma && requireNamespace("limma", quietly = TRUE)) {
    design <- cbind(Intercept = 1)
    fit <- limma::lmFit(t(X), design)
    fit <- limma::eBayes(fit)
    tstat <- fit$t[,1]; pval <- fit$p.value[,1]
  } else {
    m <- colMeans(X); s <- apply(X, 2, sd); n <- nrow(X)
    tstat <- m / (s / sqrt(n))
    # two-sided p-values using Student t with n-1 df
    pval <- 2*pt(-abs(tstat), df = n-1)
  }
  list(t = tstat, p = pval, idx_upper = which(U, arr.ind = FALSE), V = V)
}
#' Stack subject-level ΔΣ matrices into an edge matrix
#' @param fits List of `mppi_fit` objects (one per subject).
#' @param k Regressor index or name to extract.
#' @param lag Integer lag τ to extract (defaults to contemporaneous, τ = 0).
#' @param include_diag Logical; include diagonal elements when stacking.
#' @return Object of class `mppi_stack` with the stacked matrix and edge metadata.
#' @export
mppi_stack <- function(fits, k = 1L, lag = 0L, include_diag = FALSE) {
  stopifnot(is.list(fits), length(fits) >= 1)
  lag <- as.integer(lag)
  first <- fits[[1]]
  Mk0 <- mppi_get_M_lag(first, k, lag)
  dims <- dim(Mk0)
  mask <- upper.tri(Mk0, diag = include_diag)
  edges <- which(mask, arr.ind = TRUE)
  edge_labels <- apply(edges, 1, function(idx) {
    rn <- rownames(Mk0); cn <- colnames(Mk0)
    left <- if (!is.null(rn)) rn[idx[1]] else idx[1]
    right <- if (!is.null(cn)) cn[idx[2]] else idx[2]
    paste(left, right, sep = "|")
  })
  stack_mat <- matrix(NA_real_, nrow = length(fits), ncol = sum(mask))
  for (i in seq_along(fits)) {
    Mk <- mppi_get_M_lag(fits[[i]], k, lag)
    if (!all(dim(Mk) == dims)) {
      stop("All fits must have matching matrix dimensions.", call. = FALSE)
    }
    stack_mat[i, ] <- Mk[mask]
  }
  subj_names <- names(fits)
  if (is.null(subj_names)) subj_names <- paste0("subj", seq_along(fits))
  rownames(stack_mat) <- subj_names
  colnames(stack_mat) <- edge_labels
  res <- list(matrix = stack_mat,
              edges = edges,
              mask = mask,
              lag = lag,
              include_diag = include_diag,
              regressor = if (is.character(k)) k else first$names[as.integer(k)],
              dimension = dims,
              subjects = subj_names)
  class(res) <- c("mppi_stack", "list")
  res
}

#' Group-level GLM/limma inference on stacked ΔΣ edges
#' @param fits Either an `mppi_stack` object or list of `mppi_fit` objects.
#' @param design Design matrix (subjects × predictors).
#' @param contrast Optional contrast vector (length = ncol(design)).
#' @param coef Optional coefficient index/name when `contrast` is `NULL`.
#' @param method Fitting backend: `"limma"` (default) or plain `"lm"`.
#' @param k Regressor index/name (used when `fits` is a list of fits).
#' @param lag Lag τ to analyse (passed to `mppi_stack`).
#' @param include_diag Logical; include diagonal elements when stacking.
#' @param adjust Multiple-comparison correction method for `p.adjust`.
#' @return Data frame with edge-level estimates, test statistics, p- and q-values.
#' @export
mppi_group_glm <- function(fits, design, contrast = NULL, coef = NULL,
                           method = c("limma", "lm"), k = 1L, lag = 0L,
                           include_diag = FALSE, adjust = "fdr") {
  method <- match.arg(method)
  stack_obj <- if (inherits(fits, "mppi_stack")) {
    fits
  } else {
    mppi_stack(fits, k = k, lag = lag, include_diag = include_diag)
  }
  Y <- stack_obj$matrix
  design <- as.matrix(design)
  if (nrow(design) != nrow(Y)) {
    stop("design must have one row per subject in the stack.", call. = FALSE)
  }
  subj_names <- stack_obj$subjects
  if (!is.null(subj_names)) {
    rownames(Y) <- subj_names
    if (!is.null(rownames(design))) {
      idx <- match(subj_names, rownames(design))
      if (any(is.na(idx))) stop("Design rows must be named for all subjects in the stack.", call. = FALSE)
      design <- design[idx, , drop = FALSE]
    }
  }
  edge_labels <- colnames(Y)
  edges <- stack_obj$edges
  if (!is.null(contrast)) {
    contrast <- as.numeric(contrast)
    if (length(contrast) != ncol(design)) {
      stop("contrast must have length equal to ncol(design).", call. = FALSE)
    }
  }
  if (!is.null(coef) && !is.null(contrast)) {
    warning("Ignoring 'coef' because 'contrast' was supplied.", call. = FALSE)
    coef <- NULL
  }
  if (method == "limma") {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("limma package is required for method='limma'.", call. = FALSE)
    }
    fit <- limma::lmFit(t(Y), design)
    if (!is.null(contrast)) {
      fit <- limma::contrasts.fit(fit, contrast)
      fit <- limma::eBayes(fit)
      beta <- fit$coefficients[, 1, drop = TRUE]
      tstat <- fit$t[, 1, drop = TRUE]
      pval <- fit$p.value[, 1, drop = TRUE]
    } else {
      coef_idx <- if (is.null(coef)) ncol(design) else {
        if (is.character(coef)) match(coef, colnames(design)) else as.integer(coef)
      }
      if (is.na(coef_idx) || coef_idx < 1 || coef_idx > ncol(design)) {
        stop("Invalid 'coef' specification for limma output.", call. = FALSE)
      }
      fit <- limma::eBayes(fit)
      beta <- fit$coefficients[, coef_idx]
      tstat <- fit$t[, coef_idx]
      pval <- fit$p.value[, coef_idx]
    }
    estimate <- beta
    df_val <- if (length(fit$df.total) == 1) fit$df.total else fit$df.total[1]
  } else {
    qrX <- qr(design)
    if (qrX$rank < ncol(design)) {
      stop("Design matrix is rank deficient.", call. = FALSE)
    }
    beta_mat <- qr.coef(qrX, Y)
    resid <- Y - design %*% beta_mat
    df_val <- nrow(Y) - qrX$rank
    sigma2 <- colSums(resid^2) / df_val
    XtX_inv <- chol2inv(qr.R(qrX))
    if (!is.null(contrast)) {
      cvec <- contrast
    } else {
      coef_idx <- if (is.null(coef)) ncol(design) else {
        if (is.character(coef)) match(coef, colnames(design)) else as.integer(coef)
      }
      if (is.na(coef_idx) || coef_idx < 1 || coef_idx > ncol(design)) {
        stop("Invalid 'coef' specification for lm output.", call. = FALSE)
      }
      cvec <- rep(0, ncol(design)); cvec[coef_idx] <- 1
    }
    estimate <- as.numeric(cvec %*% beta_mat)
    var_c <- as.numeric(t(cvec) %*% XtX_inv %*% cvec)
    se <- sqrt(pmax(var_c, 0) * sigma2)
    zero_var <- se == 0
    se[zero_var] <- NA_real_
    tstat <- estimate / se
    pval <- 2 * stats::pt(-abs(tstat), df = df_val)
  }
  padj <- stats::p.adjust(pval, method = adjust)
  result <- data.frame(edge = edge_labels,
                       i = edges[, 1],
                       j = edges[, 2],
                       estimate = estimate,
                       statistic = tstat,
                       df = rep(df_val, length(estimate)),
                       p = pval,
                       q = padj,
                       lag = stack_obj$lag,
                       regressor = stack_obj$regressor,
                       stringsAsFactors = FALSE)
  attr(result, "stack") <- stack_obj
  attr(result, "method") <- method
  result
}
