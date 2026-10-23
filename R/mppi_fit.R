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
  hrf_vec <- if (is.null(cfg$hrf)) mppi_default_hrf(tr = tr_val, duration = duration_val, oversample = oversample_val) else .mppi_hrf_to_vec(cfg$hrf, tr_val, span = duration_val)
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
#' @param X T x q design matrix (intercept + confounds + psychological mains).
#'   Missing, empty or duplicated column names are replaced by unique defaults:
#'   `psych1`, ..., `psychK` for the psychological columns (in index order) and
#'   `x<j>` for column `j` otherwise; these name the returned regressors.
#' @param psych_idx integer indices of psychological columns in X
#' @param runs optional run labels (integer, character or factor), one per row
#'   of `Y`; factors are stored as character labels
#' @param zero_diag logical; zero the diagonal of ΔΣ. With a `basis`, the stored
#'   r x r matrices are component-space summaries; ROI-space edges come from
#'   [mppi_reconstruct_delta()], [mppi_edge()] and [mppi_seedset_slopes()], which
#'   back-project the full component slope and apply `zero_diag` and scaling in ROI space.
#' @param scale summary scale for returned matrices; choose among `"normalized"`, `"cov"`, or `"corr"`
#' @param center_by strategy for mean-centering psychological regressors (`"none"` or per-`runs`)
#' @param na_action how to handle missing trials (`"omit_tr"` to drop them or `"error"` to abort)
#' @param backend computational backend for cross-products (`"blas"`, `"accumulate"`, or `"chunked"`)
#' @param packed logical; return upper-triangular results in packed vector form
#' @param chunk_size optional column chunk size when using the chunked backend
#' @param basis optional basis definition used to project residuals before fitting
#' @param project_backend backend for projecting onto the supplied basis
#' @param project_chunk_cols optional number of columns per projection chunk
#' @param domain working domain (`"bold"`, `"neural"`, or `"innovations"`); enables deconvolution if not bold
#' @param deconv optional deconvolution configuration list used when `domain != "bold"`:
#'   `tr`; `hrf`, a numeric kernel sampled at `tr` from lag zero or a function of
#'   time in seconds (e.g. an fmrihrf HRF) evaluated on `seq(0, span, by = tr)`,
#'   `span` being its `span` attribute, else `duration` (default 32 s); default
#'   [mppi_default_hrf()]; `type` (`"dct"` or `"map"`); `lambda`, relative to the
#'   operator scale so results do not depend on HRF amplitude: for `"map"` the
#'   first-difference penalty is `lambda` (default 10) times the mean diagonal of
#'   \eqn{H'H} (see [mppi_deconv()]), for `"dct"` with `method = "fixed"` the ridge
#'   is `lambda` (default 0.1) times the mean diagonal of the convolved-basis Gram
#'   matrix (GCV is scale-free); `method`, `K`, `q`, `groups`, `psych`.
#' @param lags Integer vector of lags τ to compute (τ = 0 yields contemporaneous slopes)
#' @param lag_blocklens Optional run lengths on the acquisition grid (summing to
#'   `nrow(Y)`) for lagged computations when `runs` is unavailable
#' @return list with contemporaneous and lagged slope matrices plus residual
#'   metadata. `lag_blocklens` holds the lengths of the contiguous blocks of the
#'   retained rows of `R` (`NULL` for one block): blocks break at run boundaries
#'   and at censored scans, so adjacent rows within a block are adjacent scans;
#'   `observation_index` gives each retained row's acquisition position.
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
  if (is.factor(runs)) runs <- as.character(runs)

  original_n <- nrow(Y)
  original_runs <- runs
  original_blocklens <- lag_blocklens
  inputs <- .mppi_prepare_fit_inputs(Y, X, psych_idx, runs = runs,
                                     basis = basis, na_action = na_action_choice)
  if (inputs$dropped > 0L && domain != "bold")
    stop("Neural and innovations domains do not support censored observations; a masked deconvolution operator is required.", call. = FALSE)
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
  obs_idx <- setdiff(seq_len(original_n), inputs$dropped_idx)
  if (inputs$dropped > 0L) {
    if (!is.null(lag_blocklens) && (any(original_blocklens <= 0) || sum(original_blocklens) != original_n)) stop("lag_blocklens must describe the original acquisition grid.", call. = FALSE)
    segments <- if (!is.null(original_blocklens)) rep(seq_along(original_blocklens), original_blocklens) else if (!is.null(original_runs)) cumsum(c(TRUE, original_runs[-1] != original_runs[-original_n])) else rep(1L, original_n)
    # Stored blocks of retained rows break at run boundaries and censored gaps.
    lag_blocklens <- rle(cumsum(c(TRUE, diff(obs_idx) != 1L | diff(segments[obs_idx]) != 0L)))$lengths
  } else if (!is.null(lag_blocklens)) {
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
    R_base_full <- .mppi_residualize(Y, X_base)
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

  pNms <- colnames(X)[p_idx]
  deconv_psych <- NULL
  neural_lambda <- NULL
  n_psych <- length(p_idx)
  if (isTRUE(deconv_cfg$active)) {
    if (domain == "innovations" && !is.null(runs) && length(unique(runs)) > 1L) stop("Multi-run innovations are unsupported: run-local filtering is required.", call. = FALSE)
    run_indices <- if (is.null(runs)) list(seq_len(Tn)) else split(seq_len(Tn), cumsum(c(TRUE, diff(as.integer(as.factor(runs))) != 0)))
    groups <- deconv_cfg$groups
    if (is.null(groups)) {
      groups <- tryCatch(mppi_group_hrf_columns(X[, p_idx, drop = FALSE]),
                         error = function(e) NULL)
    }
    deconv_type <- tolower(deconv_cfg$type %||% "dct")
    if (identical(deconv_type, "map")) {
      lambda_map <- deconv_cfg$lambda
      if (is.null(lambda_map)) lambda_map <- 10
      R_unscaled <- do.call(rbind, lapply(run_indices, function(i) mppi_deconv(R_unscaled[i, , drop=FALSE], TR = deconv_cfg$tr, hrf = deconv_cfg$hrf, lambda = lambda_map)))
      R_work <- R_unscaled
      if (is.null(deconv_cfg$sticks)) {
        S_map <- do.call(rbind, lapply(run_indices, function(i) mppi_deconv(X[i, p_idx, drop = FALSE], TR = deconv_cfg$tr, hrf = deconv_cfg$hrf, lambda = lambda_map)))
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
          deconv_psych[, ii] <- .mppi_residualize_vec(deconv_psych[, ii], Qs)
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
      if (is.null(K)) K <- min(64L, floor(min(lengths(run_indices)) / 2L))
      method_deconv <- deconv_cfg$method
      if (is.null(method_deconv)) method_deconv <- "gcv"
      lambda_fixed <- deconv_cfg$lambda
      if (is.null(lambda_fixed)) lambda_fixed <- 1e-1
      q_vec <- deconv_cfg$q
      dec_parts <- lapply(run_indices, function(i) mppi_deconv_dct(R_unscaled[i, , drop=FALSE], h = deconv_cfg$hrf, K = K, method = method_deconv, lambda = lambda_fixed, q = q_vec))
      dec_out <- list(U = do.call(rbind, lapply(dec_parts, `[[`, "U")), lambda = lapply(dec_parts, `[[`, "lambda"))
      R_unscaled <- dec_out$U
      R_work <- R_unscaled
      if (is.null(deconv_cfg$sticks)) {
        ps_parts <- lapply(run_indices, function(i) mppi_psych_neural_from_X(X[i, , drop=FALSE], p_idx, h = deconv_cfg$hrf, K = K, method = method_deconv, lambda = lambda_fixed, groups = groups))
        deconv_psych <- do.call(rbind, lapply(ps_parts, `[[`, "S"))
        pNms <- ps_parts[[1]]$names
        neural_lambda <- lapply(ps_parts, `[[`, "lambda")
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
  # ROI-space SDs for back-projected basis slopes: voxel residuals in BOLD; otherwise the
  # back-projected working series, since deconvolution/filtering runs on components.
  sigma_roi <- if (is.null(basis_info)) NULL else if (domain == "bold") matrixStats::colSds(R) else
    sqrt(rowSums((basis_info$V %*% stats::cov(R_work)) * basis_info$V))

  out  <- vector("list", n_psych)
  pks  <- vector("list", n_psych)
  denoms <- numeric(n_psych)
  raw_list <- vector("list", n_psych)
  amp_list <- vector("list", n_psych)
  norm_list <- vector("list", n_psych)
  covariance_diagonal <- vector("list", n_psych)

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
    }
    covariance_diagonal[[idx]] <<- diag(raw)
    if (zero_diag && !all(is.na(raw))) diag(raw) <- 0
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
  df_full <- qr(X)$rank
  df_base <- if (is.null(X_base)) 0L else qr(X_base)$rank
  rss_full_safe <- pmax(rss_full_raw, .Machine$double.eps)
  rss_base_safe <- pmax(rss_base_raw, .Machine$double.eps)
  aic_full_vec <- n_obs * log(rss_full_safe / n_obs) + 2 * df_full
  aic_base_vec <- n_obs * log(rss_base_safe / n_obs) + 2 * df_base
  bic_full_vec <- n_obs * log(rss_full_safe / n_obs) + df_full * log(n_obs)
  bic_base_vec <- n_obs * log(rss_base_safe / n_obs) + df_base * log(n_obs)
  total_AIC_full <- sum(aic_full_vec)
  total_BIC_full <- sum(bic_full_vec)
  evidence_info <- list(
    estimand = "task_main_effects",
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
      if (inputs$dropped > 0L) {
        # Pairs are formed on the acquisition grid, so a lag may span a censored scan.
        grid_R <- matrix(0, original_n, ncol(R_work)); grid_R[obs_idx, ] <- R_work
        grid_pk <- numeric(original_n); grid_pk[obs_idx] <- pk_vec
        valid <- rep(FALSE, original_n); valid[obs_idx] <- TRUE
        a <- if (abs(lg) >= original_n) integer(0) else seq_len(original_n - abs(lg)); b <- a + abs(lg)
        ok <- valid[a] & valid[b] & segments[a] == segments[b]
        a <- a[ok]; b <- b[ok]
        denom <- sum(grid_pk[a]^2)
        raw_lag <- if (denom < .Machine$double.eps) matrix(NA_real_, ncol(R_work), ncol(R_work)) else if (lg > 0) crossprod(grid_R[a, , drop=FALSE], grid_pk[a] * grid_R[b, , drop=FALSE]) / denom else crossprod(grid_R[b, , drop=FALSE], grid_pk[a] * grid_R[a, , drop=FALSE]) / denom
      } else raw_lag <- .mppi_crossprod_lagged(R_work, pk_vec, lg,
                                        blocklens = lag_blocklens,
                                        backend = backend, chunk_size = chunk_size)
      lag_diag <- diag(raw_lag)
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
      attr(chosen_lag, "diagonal") <- lag_diag
      lag_list[[li]] <- chosen_lag
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
      Delta_full <- raw_list[[ii]]
      diag(Delta_full) <- covariance_diagonal[[ii]]
      if (all(is.na(Delta_full))) {
        var_list[ii] <- list(NULL)
        partial_list[ii] <- list(NULL)
        next
      }
      dec <- mppi_decompose_variance(cov_matrix_source, pks[[ii]], Delta_full)
      var_list[[ii]] <- dec
      DeltaTheta <- -Theta0 %*% Delta_full %*% Theta0
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
              sigma = sigma_unscaled, sigma_roi = sigma_roi, pk = pks,
              scale = scale_choice, dropped = inputs$dropped, center_by = center_by,
              observation_index = obs_idx,
              covariance_diagonal = covariance_diagonal, zero_diag = zero_diag,
              fitted_X = X, psych_idx = p_idx,
              runs = runs, denom = denoms, backend = backend,
              n_used = nrow(R_work), Sigma0 = Sigma0, partial_lambda = partial_lambda,
              Theta0 = Theta0, variance = var_list, partial = partial_list,
              lags = lags, lagged = lag_store, lag_blocklens = lag_blocklens,
              AIC_total = rep(total_AIC_full, length(out)),
              BIC_total = rep(total_BIC_full, length(out)),
              packed = packed, chunk_size = chunk_size, basis = basis_info,
              project_backend = project_backend,
              project_chunk_cols = project_chunk_cols,
              domain = domain, estimand = "covariance_modulation", schema_version = 2L, deconv = deconv_info,
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
#' @param psych regex patterns or column indices identifying psychological regressors
#' @inheritParams mppi_fit
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
  X <- as.matrix(fmrireg::design_matrix(fmodel))
  if (is.null(colnames(X)) && is.character(psych)) stop("Design has no column names; supply indices for 'psych'.")
  pidx <- if (is.numeric(psych)) psych else mppi_select_psych(X, patterns = psych)
  Y <- as.matrix(fmridataset::collect_assay(dataset))
  schema <- fmridataset::temporal_schema(dataset)
  tr <- unique(as.numeric(schema$TR))
  if (length(tr) != 1L || !is.finite(tr) || tr <= 0 || !isTRUE(schema$contiguous)) stop("Frame requires regular, common-TR sampling.", call. = FALSE)
  if (!is.null(schema$censor) && any(schema$censor)) stop("Censored frames are unsupported by this adapter.", call. = FALSE)
  runs <- runs %||% schema$block_ids
  deconv <- deconv %||% list()
  supplied <- deconv$tr %||% deconv$TR
  if (!is.null(supplied) && (length(supplied) != 1L || abs(supplied - tr) > 1e-6)) stop("TR mismatch between frame and deconvolution.", call. = FALSE)
  deconv$tr <- tr
  scale_choice <- match.arg(scale)
  mppi_fit(Y, X, pidx, runs = runs, zero_diag = zero_diag, scale = scale_choice,
           center_by = match.arg(center_by), na_action = match.arg(na_action),
           backend = match.arg(backend), packed = packed, chunk_size = chunk_size,
           basis = basis, project_backend = match.arg(project_backend),
           project_chunk_cols = project_chunk_cols,
           domain = match.arg(domain), deconv = deconv,
           lags = lags, lag_blocklens = lag_blocklens)
}
