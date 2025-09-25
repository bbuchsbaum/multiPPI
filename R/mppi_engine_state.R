# Engine state construction -------------------------------------------------

#' Internal helper: construct shared state for mPPI engines
#'
#' Mirrors the preprocessing path inside `mppi_fit()`, returning the
#' intermediate quantities so alternative engines (e.g. instantaneous) can
#' operate on the same prepared data without duplicating logic.
#'
#' @keywords internal
.mppi_engine_state <- function(Y, X, psych_idx, runs = NULL,
                              basis = NULL,
                              domain = c("bold", "neural", "innovations"),
                              deconv = NULL,
                              center_by = c("none", "run"),
                              na_action = c("omit_tr", "error"),
                              project_backend = c("blas", "chunked"),
                              project_chunk_cols = NULL) {
  stopifnot(is.matrix(Y), is.matrix(X))
  domain <- match.arg(domain)
  center_by <- match.arg(center_by)
  na_action <- match.arg(na_action)
  project_backend_choice <- match.arg(project_backend)

  inputs <- .mppi_prepare_fit_inputs(Y, X, psych_idx, runs = runs,
                                     basis = basis, na_action = na_action)

  Y_use <- inputs$Y
  X_use <- inputs$X
  runs_use <- inputs$runs
  p_idx <- sort(unique(inputs$psych_idx))
  base_idx <- inputs$base_idx
  basis_info <- inputs$basis

  if (center_by == "run" && is.null(runs_use)) {
    warning("center_by='run' requested but 'runs' is NULL; skipping run centering.", call. = FALSE)
    center_by <- "none"
  }

  X_base <- if (length(base_idx)) X_use[, base_idx, drop = FALSE] else NULL
  if (!is.null(X_base)) {
    qr_base <- qr(X_base)
    R_base_full <- Y_use - X_base %*% qr.coef(qr_base, Y_use)
  } else {
    R_base_full <- Y_use
  }

  R_res <- .mppi_residualize(Y_use, X_use)
  rss_full_raw <- colSums(R_res^2)
  rss_base_raw <- colSums(R_base_full^2)

  R_unscaled <- if (is.null(basis_info)) {
    R_res
  } else {
    .mppi_project_basis(R_res, basis_info$V, backend = project_backend_choice,
                        chunk_cols = project_chunk_cols)
  }

  deconv_cfg <- .mppi_resolve_deconv_config(domain, deconv, default_tr = 1,
                                            default_duration = 32,
                                            default_oversample = 1)
  deconv_info <- list(domain = domain)
  R_work <- R_unscaled
  n_psych <- length(p_idx)
  pNms <- if (!is.null(colnames(X_use))) colnames(X_use)[p_idx] else paste0("psych", seq_len(n_psych))
  groups <- deconv_cfg$groups

  if (isTRUE(deconv_cfg$active)) {
    if (is.null(groups)) {
      groups <- tryCatch(mppi_group_hrf_columns(X_use[, p_idx, drop = FALSE]),
                         error = function(e) NULL)
    }
    deconv_type <- tolower(deconv_cfg$type %||% "dct")
    if (identical(deconv_type, "map")) {
      lambda_map <- deconv_cfg$lambda
      if (is.null(lambda_map)) lambda_map <- 10
      R_unscaled <- mppi_deconv(R_unscaled, TR = deconv_cfg$tr, hrf = deconv_cfg$hrf, lambda = lambda_map)
      R_work <- R_unscaled
      if (is.null(deconv_cfg$sticks)) {
        S_map <- mppi_deconv(X_use[, p_idx, drop = FALSE], TR = deconv_cfg$tr,
                             hrf = deconv_cfg$hrf, lambda = lambda_map)
        if (!is.null(groups)) {
          grp_names <- character(length(groups))
          Slist <- vector("list", length(groups))
          for (g in seq_along(groups)) {
            idx <- intersect(groups[[g]], seq_len(ncol(S_map)))
            if (!length(idx)) stop("Empty group in 'groups'.")
            Ug <- S_map[, idx, drop = FALSE]
            norms <- sqrt(colSums(Ug^2)); norms[norms == 0] <- 1
            Ug <- sweep(Ug, 2, norms, "/")
            Slist[[g]] <- rowMeans(Ug)
            if (!is.null(colnames(X_use))) {
              grp_names[g] <- paste0("s:", paste(colnames(X_use)[p_idx][idx], collapse = "+"))
            } else {
              grp_names[g] <- paste0("s:grp", g)
            }
          }
          deconv_psych <- do.call(cbind, Slist)
          colnames(deconv_psych) <- grp_names
        } else {
          deconv_psych <- S_map
          colnames(deconv_psych) <- if (!is.null(colnames(X_use))) paste0("s:", colnames(X_use)[p_idx]) else paste0("s", seq_len(ncol(S_map)))
        }
        for (ii in seq_len(ncol(deconv_psych))) {
          Qs <- cbind(1, deconv_psych[, setdiff(seq_len(ncol(deconv_psych)), ii), drop = FALSE])
          deconv_psych[, ii] <- deconv_psych[, ii] - Qs %*% qr.coef(qr(Qs), deconv_psych[, ii])
        }
        pNms <- colnames(deconv_psych)
      } else {
        deconv_psych <- as.matrix(deconv_cfg$sticks)
        if (nrow(deconv_psych) != nrow(X_use)) stop("deconv$psych must have same rows as design.")
        if (!is.null(colnames(deconv_psych))) pNms <- colnames(deconv_psych)
      }
      neural_lambda <- rep(lambda_map, ncol(deconv_psych))
      n_psych <- ncol(deconv_psych)
      deconv_info <- c(deconv_info,
                       list(type = "map", lambda = lambda_map, lambda_psych = neural_lambda,
                            hrf = deconv_cfg$hrf, groups = groups, psych_names = pNms,
                            tr = deconv_cfg$tr, duration = deconv_cfg$duration,
                            oversample = deconv_cfg$oversample,
                            sticks = deconv_psych))
    } else {
      K <- deconv_cfg$K
      if (is.null(K)) K <- min(64L, floor(nrow(R_unscaled) / 2L))
      method_deconv <- deconv_cfg$method
      if (is.null(method_deconv)) method_deconv <- "gcv"
      lambda_fixed <- deconv_cfg$lambda
      if (is.null(lambda_fixed)) lambda_fixed <- 1e-1
      q_vec <- deconv_cfg$q
      dec_out <- mppi_deconv_dct(R_unscaled, h = deconv_cfg$hrf, K = K,
                                 method = method_deconv, lambda = lambda_fixed, q = q_vec)
      R_unscaled <- dec_out$U
      R_work <- R_unscaled
      if (is.null(deconv_cfg$sticks)) {
        ps <- mppi_psych_neural_from_X(X_use, p_idx, h = deconv_cfg$hrf, K = K,
                                       method = method_deconv, lambda = lambda_fixed,
                                       groups = groups)
        deconv_psych <- ps$S
        pNms <- ps$names
        neural_lambda <- ps$lambda
      } else {
        deconv_psych <- as.matrix(deconv_cfg$sticks)
        if (nrow(deconv_psych) != nrow(X_use)) stop("deconv$psych must have same rows as design.")
        if (!is.null(colnames(deconv_psych))) pNms <- colnames(deconv_psych)
        neural_lambda <- rep(lambda_fixed, ncol(deconv_psych))
      }
      n_psych <- ncol(deconv_psych)
      deconv_info <- c(deconv_info,
                       list(type = "dct", method = method_deconv, K = K, lambda = dec_out$lambda,
                            lambda_psych = neural_lambda, hrf = deconv_cfg$hrf, q = q_vec,
                            groups = groups, psych_names = pNms, tr = deconv_cfg$tr,
                            duration = deconv_cfg$duration, oversample = deconv_cfg$oversample,
                            sticks = deconv_psych))
    }
    if (identical(deconv_cfg$domain, "innovations") && exists("deconv_psych")) {
      innov_cfg <- deconv_cfg$innov
      if (is.null(innov_cfg)) {
        innov_cfg <- list(protect = "lowK", K_keep = 8L, f_cut = 0.03,
                          ar_order = "auto", p = 4L, center = TRUE)
      }
      coln <- colnames(deconv_psych)
      deconv_psych <- mppi_innovations(deconv_psych,
                                       TR = deconv_cfg$tr,
                                       protect = innov_cfg$protect,
                                       K_keep = innov_cfg$K_keep,
                                       f_cut = innov_cfg$f_cut,
                                       ar_order = innov_cfg$ar_order,
                                       p = innov_cfg$p,
                                       center = innov_cfg$center)
      if (!is.null(coln)) colnames(deconv_psych) <- coln
      deconv_info$innov <- c(innov_cfg, list(tr = deconv_cfg$tr))
    }
  }

  if (identical(deconv_cfg$domain, "innovations")) {
    innov_cfg <- deconv_cfg$innov
    if (is.null(innov_cfg)) {
      innov_cfg <- list(protect = "lowK", K_keep = 8L, f_cut = 0.03,
                        ar_order = "auto", p = 4L, center = TRUE)
    }
    R_work <- mppi_innovations(R_work,
                               TR = deconv_cfg$tr,
                               protect = innov_cfg$protect,
                               K_keep = innov_cfg$K_keep,
                               f_cut = innov_cfg$f_cut,
                               ar_order = innov_cfg$ar_order,
                               p = innov_cfg$p,
                               center = innov_cfg$center)
    attr(R_work, "TR") <- deconv_cfg$tr
    if (is.null(deconv_info$innov)) {
      deconv_info$innov <- c(innov_cfg, list(tr = deconv_cfg$tr))
    }
  }

  sigma_unscaled <- matrixStats::colSds(R_work)

  pk_list <- vector("list", n_psych)
  names(pk_list) <- pNms
  denoms <- numeric(n_psych)
  names(denoms) <- pNms

  maybe_center <- function(vec) {
    if (center_by == "run" && !is.null(runs_use)) {
      .mppi_center_by_run(vec, runs_use)
    } else {
      vec
    }
  }

  if (isTRUE(deconv_cfg$active) && exists("deconv_psych")) {
    idx_seq <- seq_len(ncol(deconv_psych))
    for (ii in idx_seq) {
      sk <- deconv_psych[, ii]
      S_neural <- cbind(X_use[, base_idx, drop = FALSE],
                        if (length(idx_seq) > 1) deconv_psych[, setdiff(idx_seq, ii), drop = FALSE] else NULL)
      pk_vec <- if (is.null(S_neural) || ncol(S_neural) == 0) sk else .mppi_residualize_vec(sk, as.matrix(S_neural))
      pk_vec <- maybe_center(pk_vec)
      denoms[ii] <- sum(pk_vec^2)
      pk_list[[ii]] <- pk_vec
    }
  } else {
    for (ii in seq_len(n_psych)) {
      k <- p_idx[ii]
      Q <- X_use[, c(base_idx, setdiff(p_idx, k)), drop = FALSE]
      pk_vec <- .mppi_residualize_vec(X_use[, k], Q)
      pk_vec <- maybe_center(pk_vec)
      denoms[ii] <- sum(pk_vec^2)
      pk_list[[ii]] <- pk_vec
    }
  }

  cov_matrix_source <- if (identical(deconv_cfg$domain, "innovations")) R_work else R_unscaled
  Sigma0 <- crossprod(cov_matrix_source) / nrow(cov_matrix_source)

  n_obs <- nrow(Y_use)
  df_full <- ncol(X_use)
  df_base <- if (is.null(X_base)) 0L else ncol(X_base)
  rss_full_safe <- pmax(rss_full_raw, .Machine$double.eps)
  rss_base_safe <- pmax(rss_base_raw, .Machine$double.eps)
  aic_full_vec <- n_obs * log(rss_full_safe / n_obs) + 2 * df_full
  aic_base_vec <- n_obs * log(rss_base_safe / n_obs) + 2 * df_base
  bic_full_vec <- n_obs * log(rss_full_safe / n_obs) + df_full * log(n_obs)
  bic_base_vec <- n_obs * log(rss_base_safe / n_obs) + df_base * log(n_obs)

  evidence_info <- list(
    n = n_obs,
    df_full = df_full,
    df_base = df_base,
    rss_full = rss_full_raw,
    rss_base = rss_base_raw,
    column = colnames(Y_use),
    delta_AIC_total = sum(aic_base_vec - aic_full_vec),
    delta_BIC_total = sum(bic_base_vec - bic_full_vec),
    AIC_full_total = sum(aic_full_vec),
    AIC_base_total = sum(aic_base_vec),
    BIC_full_total = sum(bic_full_vec),
    BIC_base_total = sum(bic_base_vec)
  )

  list(
    Y = Y_use,
    X = X_use,
    runs = runs_use,
    p_idx = p_idx,
    base_idx = base_idx,
    basis = basis_info,
    R_work = R_work,
    R_unscaled = R_unscaled,
    sigma = sigma_unscaled,
    pk_list = pk_list,
    pk_names = pNms,
    denoms = denoms,
    Sigma0 = Sigma0,
    evidence = evidence_info,
    dropped = inputs$dropped,
    dropped_idx = inputs$dropped_idx,
    center_by = center_by,
    domain = domain,
    deconv = deconv_info,
    cov_source = cov_matrix_source,
    n_obs = n_obs,
    project_backend = project_backend_choice,
    project_chunk_cols = project_chunk_cols
  )
}

