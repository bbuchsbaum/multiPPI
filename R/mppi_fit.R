# Core estimators ----------------------------------------------------------

.mppi_prepare_fit_inputs <- function(Y, X, psych_idx, runs = NULL,
                                     basis = NULL,
                                     na_action = c("omit_tr", "error")) {
  stopifnot(is.matrix(Y), is.matrix(X))
  if (is.null(psych_idx) || !length(psych_idx)) {
    stop("psych_idx must be a non-empty vector of column indices.", call. = FALSE)
  }
  if (!is.null(runs) && length(runs) != nrow(Y)) {
    stop("Length of 'runs' must match number of rows in Y.", call. = FALSE)
  }
  na_action <- match.arg(na_action)
  design_prep <- .mppi_prepare_design_matrix(X, psych_idx)
  X_use <- design_prep$X
  p_idx <- design_prep$psych_idx
  prep <- .mppi_handle_na(Y, X_use, runs, na_action = na_action)
  basis_info <- .mppi_check_basis(basis, ncol(prep$Y))
  list(
    Y = prep$Y,
    X = prep$X,
    runs = prep$runs,
    psych_idx = p_idx,
    base_idx = setdiff(seq_len(ncol(prep$X)), p_idx),
    dropped = prep$dropped,
    dropped_idx = prep$dropped_idx,
    basis = basis_info
  )
}

.mppi_resolve_deconv_config <- function(domain, deconv, default_tr = 1,
                                        default_duration = 32,
                                        default_oversample = 1) {
  domain <- match.arg(domain, c("bold", "neural", "innovations"))
  if (identical(domain, "bold")) {
    return(list(active = FALSE, domain = "bold"))
  }
  cfg <- if (is.null(deconv)) list() else deconv
  tr_val <- cfg$tr
  if (is.null(tr_val)) tr_val <- cfg$TR
  if (is.null(tr_val)) tr_val <- default_tr
  duration_val <- cfg$duration
  if (is.null(duration_val)) duration_val <- default_duration
  oversample_val <- cfg$oversample
  if (is.null(oversample_val)) oversample_val <- default_oversample
  hrf_input <- cfg$hrf
  if (is.function(hrf_input)) {
    hrf_vec <- as.numeric(hrf_input(tr = tr_val, duration = duration_val, oversample = oversample_val))
  } else if (is.null(hrf_input)) {
    hrf_vec <- mppi_default_hrf(tr = tr_val, duration = duration_val, oversample = oversample_val)
  } else {
    hrf_vec <- as.numeric(hrf_input)
  }
  if (!length(hrf_vec) || any(!is.finite(hrf_vec))) {
    stop("Resolved HRF for neural domain must be finite and non-empty.", call. = FALSE)
  }
  out <- list(
    active = TRUE,
    domain = domain,
    type = tolower(cfg$type %||% "dct"),
    hrf = hrf_vec,
    tr = tr_val,
    duration = duration_val,
    oversample = oversample_val,
    groups = cfg$groups,
    lambda = cfg$lambda,
    method = cfg$method,
    K = cfg$K,
    q = cfg$q,
    sticks = cfg$psych
  )
  if (identical(domain, "innovations")) {
    prot <- cfg$protect %||% "lowK"
    prot <- match.arg(prot, c("none", "lowK", "lowHz"))
    ar_mode <- cfg$ar_order %||% cfg$ar %||% "auto"
    ar_mode <- match.arg(ar_mode, c("auto", "fixed"))
    out$innov <- list(
      protect = prot,
      K_keep = as.integer(cfg$K_keep %||% 8L),
      f_cut = cfg$f_cut %||% 0.03,
      ar_order = ar_mode,
      p = as.integer(cfg$p %||% 4L),
      center = isTRUE(cfg$center %||% TRUE)
    )
  }
  out
}

#' Matrix-PPI estimator from raw matrices
#' @param Y T x V matrix of time series (preprocessed)
#' @param X T x q design matrix (intercept + confounds + psychological mains)
#' @param psych_idx integer indices of psychological columns in X
#' @param zero_diag logical; zero the diagonal of ΔΣ
#' @param scale "cov" (default) or "corr" to standardize R to unit SD
#' @param lags Integer vector of lags τ to compute (τ = 0 yields contemporaneous slopes).
#' @param lag_blocklens Optional run lengths for lagged computations when `runs` is unavailable.
#' @return list with contemporaneous and lagged slope matrices plus residual metadata.
mppi_fit <- function(Y, X, psych_idx, runs = NULL,
                     zero_diag = TRUE, scale = c("normalized","cov","corr"),
                     center_by = c("none","run"), na_action = c("omit_tr","error"),
                     backend = c("blas","accumulate","chunked"), packed = FALSE,
                     chunk_size = NULL, basis = NULL,
                     project_backend = c("blas","chunked"), project_chunk_cols = NULL,
                     domain = c("bold","neural","innovations"), deconv = NULL,
                     lags = 0L, lag_blocklens = NULL) {
  stopifnot(is.matrix(Y), is.matrix(X))
  scale_choice <- match.arg(scale)
  center_by <- match.arg(center_by)
  na_action_choice <- match.arg(na_action)
  backend <- match.arg(backend)
  project_backend_choice <- match.arg(project_backend)
  domain <- match.arg(domain)

  inputs <- .mppi_prepare_fit_inputs(Y, X, psych_idx, runs = runs,
                                     basis = basis, na_action = na_action_choice)
  Y <- inputs$Y
  X <- inputs$X
  runs <- inputs$runs
  p_idx <- sort(unique(inputs$psych_idx))
  base_idx <- inputs$base_idx
  Tn <- nrow(X)
  stopifnot(nrow(Y) == Tn)
  basis_info <- inputs$basis
  if (center_by == "run" && is.null(runs)) {
    warning("center_by='run' requested but 'runs' is NULL; skipping run centering.", call. = FALSE)
    center_by <- "none"
  }
  project_backend <- if (is.null(basis_info)) "blas" else project_backend_choice
  lags <- sort(unique(as.integer(lags)))
  if (!length(lags)) lags <- 0L
  if (!0L %in% lags) lags <- c(0L, lags)
  if (!is.null(lag_blocklens)) {
    lag_blocklens <- as.integer(lag_blocklens)
    if (any(lag_blocklens <= 0)) {
      stop("lag_blocklens must contain positive lengths.", call. = FALSE)
    }
    if (sum(lag_blocklens) != Tn) {
      stop("Sum of lag_blocklens must match number of usable time points.", call. = FALSE)
    }
  } else if (!is.null(runs)) {
    lag_blocklens <- rle(runs)$lengths
  }

  if (length(base_idx)) {
    X_base <- X[, base_idx, drop = FALSE]
    qr_base <- qr(X_base)
    R_base_full <- Y - X_base %*% qr.coef(qr_base, Y)
  } else {
    X_base <- NULL
    R_base_full <- Y
  }

  R <- .mppi_residualize(Y, X)
  rss_full_raw <- colSums(R^2)
  rss_base_raw <- colSums(R_base_full^2)
  if (!is.null(basis_info)) {
    R_unscaled <- .mppi_project_basis(R, basis_info$V, backend = project_backend,
                                      chunk_cols = project_chunk_cols)
  } else {
    R_unscaled <- R
  }
  R_work <- R_unscaled
  deconv_cfg <- .mppi_resolve_deconv_config(domain, deconv)
  deconv_info <- list(domain = domain)

  denom_warn <- function(name) {
    warning(sprintf("Regressor '%s' has near-zero variance after residualization; returning NA matrix.",
                    name), call. = FALSE)
  }

  pNms <- if (is.null(colnames(X))) paste0("psych", seq_along(p_idx)) else colnames(X)[p_idx]
  deconv_psych <- NULL
  neural_lambda <- NULL
  n_psych <- length(p_idx)
  if (isTRUE(deconv_cfg$active)) {
    groups <- deconv_cfg$groups
    if (is.null(groups)) {
      groups <- tryCatch(mppi_group_hrf_columns(X[, p_idx, drop = FALSE]),
                         error = function(e) NULL)
    }
    deconv_type <- tolower(deconv_cfg$type %||% "dct")
    if (identical(deconv_type, "map")) {
      lambda_map <- deconv_cfg$lambda
      if (is.null(lambda_map)) lambda_map <- 10
      R_unscaled <- mppi_deconv(R_unscaled, TR = deconv_cfg$tr, hrf = deconv_cfg$hrf, lambda = lambda_map)
      R_work <- R_unscaled
      if (is.null(deconv_cfg$sticks)) {
        S_map <- mppi_deconv(X[, p_idx, drop = FALSE], TR = deconv_cfg$tr, hrf = deconv_cfg$hrf, lambda = lambda_map)
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
            if (!is.null(colnames(X))) {
              grp_names[g] <- paste0("s:", paste(colnames(X)[p_idx][idx], collapse = "+"))
            } else {
              grp_names[g] <- paste0("s:grp", g)
            }
          }
          deconv_psych <- do.call(cbind, Slist)
          colnames(deconv_psych) <- grp_names
        } else {
          deconv_psych <- S_map
          colnames(deconv_psych) <- if (!is.null(colnames(X))) paste0("s:", colnames(X)[p_idx]) else paste0("s", seq_len(ncol(S_map)))
        }
        for (ii in seq_len(ncol(deconv_psych))) {
          Qs <- cbind(1, deconv_psych[, setdiff(seq_len(ncol(deconv_psych)), ii), drop = FALSE])
          deconv_psych[, ii] <- deconv_psych[, ii] - Qs %*% qr.coef(qr(Qs), deconv_psych[, ii])
        }
        pNms <- colnames(deconv_psych)
      } else {
        deconv_psych <- as.matrix(deconv_cfg$sticks)
        if (nrow(deconv_psych) != nrow(X)) stop("deconv$psych must have same rows as design.")
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
        ps <- mppi_psych_neural_from_X(X, p_idx, h = deconv_cfg$hrf, K = K,
                                       method = method_deconv, lambda = lambda_fixed,
                                       groups = groups)
        deconv_psych <- ps$S
        pNms <- ps$names
        neural_lambda <- ps$lambda
      } else {
        deconv_psych <- as.matrix(deconv_cfg$sticks)
        if (nrow(deconv_psych) != nrow(X)) stop("deconv$psych must have same rows as design.")
        if (!is.null(colnames(deconv_psych))) pNms <- colnames(deconv_psych)
      }
      n_psych <- ncol(deconv_psych)
      deconv_info <- c(deconv_info,
                       list(type = "dct", method = method_deconv, K = K, lambda = dec_out$lambda,
                            lambda_psych = neural_lambda, hrf = deconv_cfg$hrf, q = q_vec,
                            groups = groups, psych_names = pNms, tr = deconv_cfg$tr,
                            duration = deconv_cfg$duration, oversample = deconv_cfg$oversample,
                            sticks = deconv_psych))
    }
    if (identical(deconv_cfg$domain, "innovations") && !is.null(deconv_psych)) {
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

  out  <- vector("list", n_psych)
  pks  <- vector("list", n_psych)
  denoms <- numeric(n_psych)
  raw_list <- vector("list", n_psych)
  amp_list <- vector("list", n_psych)
  norm_list <- vector("list", n_psych)

  compute_entry <- function(pk_vec, idx, name) {
    if (center_by == "run" && !is.null(runs)) {
      pk_vec <- .mppi_center_by_run(pk_vec, runs)
    }
    denom <- sum(pk_vec^2)
    denoms[idx] <<- denom
    if (denom < .Machine$double.eps) {
      denom_warn(name)
      raw <- matrix(NA_real_, ncol(R_work), ncol(R_work))
    } else {
      raw <- .mppi_crossprod(R_work, pk_vec, backend = backend, chunk_size = chunk_size) / denom
      if (zero_diag && !all(is.na(raw))) diag(raw) <- 0
    }
    amp <- if (all(is.na(raw))) raw else .mppi_scale_matrix(raw, sigma_unscaled, "amplitude")
    norm <- if (all(is.na(raw))) raw else .mppi_scale_matrix(raw, sigma_unscaled, "normalized")
    raw_list[[idx]] <<- raw
    amp_list[[idx]] <<- amp
    norm_list[[idx]] <<- norm
    pks[[idx]] <<- pk_vec

    chosen <- switch(scale_choice,
                     normalized = norm,
                     cov = raw,
                     corr = norm)
    attr(chosen, "raw") <- raw
    attr(chosen, "amplitude") <- amp
    attr(chosen, "normalized") <- norm

    if (packed) {
      pack_obj <- .mppi_pack_upper(chosen)
      attr(pack_obj, "raw") <- raw
      attr(pack_obj, "amplitude") <- amp
      attr(pack_obj, "normalized") <- norm
      out[[idx]] <<- pack_obj
    } else {
      out[[idx]] <<- chosen
    }
  }

  if (isTRUE(deconv_cfg$active)) {
    idx_seq <- seq_len(n_psych)
    for (ii in idx_seq) {
      sk <- deconv_psych[, ii]
      S_neural <- cbind(X[, base_idx, drop = FALSE],
                        if (n_psych > 1) deconv_psych[, setdiff(idx_seq, ii), drop = FALSE] else NULL)
      pk_vec <- if (is.null(S_neural) || ncol(S_neural) == 0) sk else .mppi_residualize_vec(sk, as.matrix(S_neural))
      compute_entry(pk_vec, ii, pNms[ii])
    }
  } else {
    for (ii in seq_len(n_psych)) {
      k  <- p_idx[ii]
      Q  <- X[, c(base_idx, setdiff(p_idx, k)), drop = FALSE]
      pk_vec <- .mppi_residualize_vec(X[, k], Q)
      compute_entry(pk_vec, ii, pNms[ii])
    }
  }
  names(out) <- pNms
  names(pks) <- pNms
  names(denoms) <- pNms
  n_obs <- nrow(Y)
  df_full <- ncol(X)
  df_base <- if (is.null(X_base)) 0L else ncol(X_base)
  rss_full_safe <- pmax(rss_full_raw, .Machine$double.eps)
  rss_base_safe <- pmax(rss_base_raw, .Machine$double.eps)
  aic_full_vec <- n_obs * log(rss_full_safe / n_obs) + 2 * df_full
  aic_base_vec <- n_obs * log(rss_base_safe / n_obs) + 2 * df_base
  bic_full_vec <- n_obs * log(rss_full_safe / n_obs) + df_full * log(n_obs)
  bic_base_vec <- n_obs * log(rss_base_safe / n_obs) + df_base * log(n_obs)
  total_AIC_full <- sum(aic_full_vec)
  total_BIC_full <- sum(bic_full_vec)
  evidence_info <- list(
    n = n_obs,
    df_full = df_full,
    df_base = df_base,
    rss_full = rss_full_raw,
    rss_base = rss_base_raw,
    column = colnames(Y),
    delta_AIC_total = sum(aic_base_vec - aic_full_vec),
    delta_BIC_total = sum(bic_base_vec - bic_full_vec),
    AIC_full_total = total_AIC_full,
    AIC_base_total = sum(aic_base_vec),
    BIC_full_total = total_BIC_full,
    BIC_base_total = sum(bic_base_vec)
  )
  lag_store <- setNames(vector("list", length(out)), pNms)
  zero_idx <- match(0L, lags)
  for (ii in seq_len(n_psych)) {
    lag_list <- vector("list", length(lags))
    names(lag_list) <- as.character(lags)
    if (!is.na(zero_idx)) lag_list[[zero_idx]] <- out[[ii]]
    pk_vec <- pks[[ii]]
    for (li in seq_along(lags)) {
      lg <- lags[li]
      if (!is.na(zero_idx) && li == zero_idx) next
      raw_lag <- .mppi_crossprod_lagged(R_work, pk_vec, lg,
                                        blocklens = lag_blocklens,
                                        backend = backend, chunk_size = chunk_size)
      if (!all(is.na(raw_lag)) && zero_diag) diag(raw_lag) <- 0
      amp_lag <- if (all(is.na(raw_lag))) raw_lag else .mppi_scale_matrix(raw_lag, sigma_unscaled, "amplitude")
      norm_lag <- if (all(is.na(raw_lag))) raw_lag else .mppi_scale_matrix(raw_lag, sigma_unscaled, "normalized")
      chosen_lag <- switch(scale_choice,
                           normalized = norm_lag,
                           cov = raw_lag,
                           corr = norm_lag)
      attr(chosen_lag, "raw") <- raw_lag
      attr(chosen_lag, "amplitude") <- amp_lag
      attr(chosen_lag, "normalized") <- norm_lag
      if (packed) {
        pack_lag <- .mppi_pack_upper(chosen_lag)
        attr(pack_lag, "raw") <- raw_lag
        attr(pack_lag, "amplitude") <- amp_lag
        attr(pack_lag, "normalized") <- norm_lag
        lag_list[[li]] <- pack_lag
      } else {
        lag_list[[li]] <- chosen_lag
      }
    }
    lag_store[[ii]] <- lag_list
  }
  cov_matrix_source <- if (identical(deconv_cfg$domain, "innovations")) R_work else R_unscaled
  Sigma0 <- crossprod(cov_matrix_source) / nrow(cov_matrix_source)
  partial_lambda <- 1e-3
  var_list <- NULL
  partial_list <- NULL
  Theta0 <- NULL
  if (scale_choice == "cov") {
    Theta0 <- solve(Sigma0 + partial_lambda * diag(ncol(Sigma0)))
    var_list <- vector("list", length(out))
    partial_list <- vector("list", length(out))
    for (ii in seq_along(out)) {
      Delta_full <- if (packed) .mppi_unpack_upper(out[[ii]]) else out[[ii]]
      if (all(is.na(Delta_full))) {
        var_list[ii] <- list(NULL)
        partial_list[ii] <- list(NULL)
        next
      }
      dec <- mppi_decompose_variance(R_unscaled, pks[[ii]], Delta_full)
      var_list[[ii]] <- dec
      DeltaTheta <- mppi_to_partial(Sigma0, Delta_full, lambda = partial_lambda)
      partial_list[[ii]] <- list(DeltaTheta = DeltaTheta,
                                 DeltaRho = mppi_delta_partial(Theta0, DeltaTheta))
    }
    names(var_list) <- pNms
    names(partial_list) <- pNms
  }
  if (!is.null(basis_info)) {
    basis_info$project_backend <- project_backend
    basis_info$project_chunk_cols <- project_chunk_cols
  }
  res <- list(Delta = out, names = pNms, R = R_work, R_raw = R_unscaled,
              Delta_raw = raw_list, Delta_amplitude = amp_list, Delta_normalized = norm_list,
              sigma = sigma_unscaled, pk = pks,
              scale = scale_choice, dropped = inputs$dropped, center_by = center_by,
              runs = runs, denom = denoms, backend = backend,
              n_used = nrow(R_work), Sigma0 = Sigma0, partial_lambda = partial_lambda,
              Theta0 = Theta0, variance = var_list, partial = partial_list,
              lags = lags, lagged = lag_store, lag_blocklens = lag_blocklens,
              AIC_total = rep(total_AIC_full, length(out)),
              BIC_total = rep(total_BIC_full, length(out)),
              packed = packed, chunk_size = chunk_size, basis = basis_info,
              project_backend = project_backend,
              project_chunk_cols = project_chunk_cols,
              domain = domain, deconv = deconv_info,
              evidence = evidence_info)
  if (!is.null(basis_info)) res$Z <- R_work
  if (identical(domain, "neural") || identical(domain, "innovations")) res$U <- R_unscaled
  if (!is.null(res$AIC_total)) names(res$AIC_total) <- pNms
  if (!is.null(res$BIC_total)) names(res$BIC_total) <- pNms
  class(res) <- c("mppi_fit", "list")
  res
}

#' Convenience wrapper when using fmrireg objects
#' @param fmodel fmri_model (from fmrireg)
#' @param dataset fmrireg dataset supplying Y
#' @param psych regex patterns or indices
#' @param scale "normalized", "cov", or "corr"
mppi_fit_from_fmrireg <- function(fmodel, dataset, psych, scale = c("normalized","cov","corr"), zero_diag = TRUE,
                                  runs = NULL, center_by = c("none","run"), na_action = c("omit_tr","error"),
                                  backend = c("blas","accumulate","chunked"), chunk_size = NULL,
                                  packed = FALSE, basis = NULL,
                                  project_backend = c("blas","chunked"), project_chunk_cols = NULL,
                                  domain = c("neural","bold","innovations"), deconv = NULL,
                                  lags = 0L, lag_blocklens = NULL) {
  if (!requireNamespace("fmrireg", quietly = TRUE)) {
    stop("fmrireg not available. Install via remotes::install_github('bbuchsbaum/fmrireg').")
  }
  X <- fmrireg::design_matrix(fmodel)
  if (is.null(colnames(X)) && is.character(psych)) stop("Design has no column names; supply indices for 'psych'.")
  pidx <- if (is.numeric(psych)) psych else mppi_select_psych(X, patterns = psych)
  Y <- fmrireg::get_data_matrix(dataset)
  scale_choice <- match.arg(scale)
  mppi_fit(Y, X, pidx, runs = runs, zero_diag = zero_diag, scale = scale_choice,
           center_by = match.arg(center_by), na_action = match.arg(na_action),
           backend = match.arg(backend), packed = packed, chunk_size = chunk_size,
           basis = basis, project_backend = match.arg(project_backend),
           project_chunk_cols = project_chunk_cols,
           domain = match.arg(domain), deconv = deconv,
           lags = lags, lag_blocklens = lag_blocklens)
}
