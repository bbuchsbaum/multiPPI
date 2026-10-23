# Group-level empirical Bayes (limma) -------------------------------------

#' Moderated one-sample tests on vectorized Δ matrices across subjects
#' @param Delta_list list of VxV matrices (same V), one per subject
#' @param use_limma logical; if TRUE and limma present, use eBayes
#' @param directed logical; treat matrices as directed (no symmetrisation)
#' @param triangle which triangle/orientation to retain when vectorising
#' @param include_diag logical; include diagonal elements
#' @return list with t, p, and (if limma) moderated stats
mppi_group_ebayes <- function(Delta_list, use_limma = TRUE,
                              directed = FALSE,
                              triangle = c("upper", "lower", "both"),
                              include_diag = FALSE) {
  stopifnot(is.list(Delta_list), length(Delta_list) >= 2)
  V <- nrow(Delta_list[[1]])
  .mppi_check_group_matrices(Delta_list)
  triangle <- match.arg(triangle)
  idx_matrix <- .mppi_edge_mask(V, include_diag = include_diag,
                                directed = directed, triangle = triangle)
  X <- do.call(rbind, lapply(Delta_list, function(D) as.numeric(D[idx_matrix])))
  # One-sample vs zero
  if (use_limma && requireNamespace("limma", quietly = TRUE)) {
    design <- matrix(1, nrow = nrow(X), ncol = 1, dimnames = list(NULL, "Intercept"))
    fit <- limma::lmFit(t(X), design)
    fit <- limma::eBayes(fit)
    tstat <- fit$t[,1]; pval <- fit$p.value[,1]
  } else {
    m <- colMeans(X); s <- apply(X, 2, sd); n <- nrow(X)
    tstat <- m / (s / sqrt(n))
    # two-sided p-values using Student t with n-1 df
    pval <- 2*pt(-abs(tstat), df = n-1)
  }
  list(t = tstat, p = pval, idx = which(idx_matrix, arr.ind = FALSE), V = V,
       directed = directed, triangle = triangle, include_diag = include_diag)
}
#' Stack subject-level ΔΣ matrices into an edge matrix
#' @param fits List of `mppi_fit` objects (one per subject).
#' @param k Regressor index or name to extract.
#' @param lag Integer lag τ to extract (defaults to contemporaneous, τ = 0).
#' @param include_diag Logical; include diagonal elements when stacking.
#' @param directed Logical; treat matrices as directed (keep orientation).
#' @param triangle Which triangle/orientation to retain (`"upper"`, `"lower"`, `"both"`).
#' @return Object of class `mppi_stack` with the stacked matrix and edge metadata.
#' @export
mppi_stack <- function(fits, k = 1L, lag = 0L, include_diag = FALSE,
                       directed = FALSE, triangle = c("upper", "lower", "both")) {
  triangle <- match.arg(triangle)
  stopifnot(is.list(fits), length(fits) >= 1)
  lag <- as.integer(lag)
  first <- fits[[1]]
  chosen <- function(fit) {
    if (is.character(k)) {
      if (length(k) != 1L || is.null(fit$names) || !k %in% fit$names) stop("Condition name must exist in every fit.")
      k
    } else {
      if (length(k) != 1L || is.na(k) || k < 1 || k != floor(k) || is.null(fit$names) || k > length(fit$names)) stop("Condition index must identify a named condition in every fit.")
      fit$names[k]
    }
  }
  condition <- chosen(first)
  if (any(!vapply(fits, function(fit) identical(chosen(fit), condition), logical(1)))) stop("Condition identities must match across fits.")
  Mk0 <- mppi_get_M_lag(first, k, lag)
  dims <- dim(Mk0)
  mask <- .mppi_edge_mask(dims[1], include_diag = include_diag,
                          directed = directed, triangle = triangle)
  edges <- which(mask, arr.ind = TRUE)
  edge_labels <- apply(edges, 1, function(idx) {
    rn <- rownames(Mk0); cn <- colnames(Mk0)
    left <- if (!is.null(rn)) rn[idx[1]] else idx[1]
    right <- if (!is.null(cn)) cn[idx[2]] else idx[2]
    sep <- if (directed) "->" else "|"
    paste(left, right, sep = sep)
  })
  stack_mat <- matrix(NA_real_, nrow = length(fits), ncol = sum(mask))
  for (i in seq_along(fits)) {
    Mk <- mppi_get_M_lag(fits[[i]], k, lag)
    if (!is.matrix(Mk) || !all(dim(Mk) == dims) || !identical(dimnames(Mk), dimnames(Mk0))) {
      stop("All fits must have matching matrix dimensions and node labels.", call. = FALSE)
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
              directed = directed,
              triangle = triangle,
              regressor = if (is.character(k)) k else first$names[as.integer(k)],
              dimension = dims,
              subjects = subj_names)
  class(res) <- c("mppi_stack", "list")
  res
}

#' Stack a list of directed phase-winding matrices
#'
#' @param W_list List of square numeric matrices (same dimension).
#' @param labels Optional vector of node labels.
#' @param include_diag Logical; include diagonal terms.
#' @param triangle Which triangle/orientation to retain when stacking.
#' @return Object of class `mppi_stack` suitable for `mppi_group_glm()`.
#' @export
pwf_stack <- function(W_list, labels = NULL, include_diag = FALSE,
                      triangle = c("upper", "lower", "both")) {
  stopifnot(is.list(W_list), length(W_list) >= 1)
  triangle <- match.arg(triangle)
  V <- nrow(W_list[[1]])
  dims <- c(V, V)
  .mppi_check_group_matrices(W_list)
  if (!is.null(labels) && (length(labels) != V || anyNA(labels) || anyDuplicated(labels))) stop("labels must uniquely identify each node")
  if (!is.null(labels) && !is.null(rownames(W_list[[1]])) && !identical(as.character(labels), rownames(W_list[[1]]))) stop("labels must match matrix node labels")
  if (is.null(labels)) {
    labels <- rownames(W_list[[1]])
    if (is.null(labels)) labels <- paste0("v", seq_len(V))
  }
  mask <- .mppi_edge_mask(V, include_diag = include_diag,
                          directed = TRUE, triangle = triangle)
  edges <- which(mask, arr.ind = TRUE)
  edge_labels <- apply(edges, 1, function(idx) {
    paste(labels[idx[1]], labels[idx[2]], sep = "->")
  })
  stack_mat <- matrix(NA_real_, nrow = length(W_list), ncol = sum(mask))
  for (i in seq_along(W_list)) {
    stack_mat[i, ] <- as.numeric(W_list[[i]][mask])
  }
  subj_names <- names(W_list)
  if (is.null(subj_names)) subj_names <- paste0("subj", seq_along(W_list))
  rownames(stack_mat) <- subj_names
  colnames(stack_mat) <- edge_labels
  res <- list(matrix = stack_mat,
              edges = edges,
              mask = mask,
              lag = 0L,
              include_diag = include_diag,
              directed = TRUE,
              triangle = triangle,
              regressor = "PWF",
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
#' @param directed Logical; treat matrices as directed (no symmetrisation).
#' @param triangle Which triangle/orientation to retain (`"upper"`, `"lower"`, `"both"`).
#' @param adjust Multiple-comparison correction method for `p.adjust`.
#' @return Data frame with edge-level estimates, test statistics, p- and q-values.
#' @export
mppi_group_glm <- function(fits, design, contrast = NULL, coef = NULL,
                           method = c("limma", "lm"), k = 1L, lag = 0L,
                           include_diag = FALSE, adjust = "fdr",
                           directed = FALSE, triangle = c("upper", "lower", "both")) {
  triangle <- match.arg(triangle)
  method <- match.arg(method)
  stack_obj <- if (inherits(fits, "mppi_stack")) {
    fits
  } else {
    mppi_stack(fits, k = k, lag = lag, include_diag = include_diag,
               directed = directed, triangle = triangle)
  }
  Y <- stack_obj$matrix
  design <- as.matrix(design)
  if (nrow(design) != nrow(Y)) {
    stop("design must have one row per subject in the stack.", call. = FALSE)
  }
  subj_names <- stack_obj$subjects
  if (!is.null(subj_names)) {
    rownames(Y) <- subj_names
    # Default row names (as from model.matrix) carry no subject identity.
    if (!is.null(rownames(design)) && !identical(rownames(design), as.character(seq_len(nrow(design))))) {
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
                       directed = stack_obj$directed,
                       triangle = stack_obj$triangle,
                       include_diag = stack_obj$include_diag,
                       stringsAsFactors = FALSE)
  attr(result, "stack") <- stack_obj
  attr(result, "method") <- method
  result
}
.mppi_edge_mask <- function(V, include_diag = FALSE, directed = FALSE,
                           triangle = c("upper", "lower", "both")) {
  triangle <- match.arg(triangle)
  if (directed) {
    r <- row(matrix(0, V, V))
    c <- col(matrix(0, V, V))
    if (triangle == "both") {
      mask <- r != c
      if (include_diag) mask <- mask | (r == c)
      mask
    } else if (triangle == "upper") {
      mask <- r < c
      if (include_diag) mask <- mask | (r == c)
      mask
    } else {
      mask <- r > c
      if (include_diag) mask <- mask | (r == c)
      mask
    }
  } else {
    upper.tri(matrix(0, V, V), diag = include_diag)
  }
}

.mppi_check_group_matrices <- function(matrices) {
  first <- matrices[[1L]]
  if (!is.matrix(first) || !is.numeric(first) || nrow(first) != ncol(first)) stop("Matrices must be numeric and square.")
  dn <- dimnames(first)
  if (xor(is.null(dn[[1L]]), is.null(dn[[2L]])) || (!is.null(dn[[1L]]) && (!identical(dn[[1L]], dn[[2L]]) || anyNA(dn[[1L]]) || anyDuplicated(dn[[1L]])))) stop("Matrix node labels must be unique and match on both axes.")
  for (x in matrices) if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), dim(first)) || !identical(dimnames(x), dn) || any(!is.finite(x))) stop("Matrices must have matching dimensions, node labels, and finite values.")
  invisible(NULL)
}
