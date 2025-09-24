# Helper utilities -----------------------------------------------------------

# Internal lower-triangle vectoriser
.mppi_vec_lower <- function(M, diag = FALSE) {
  if (!is.matrix(M)) stop("Expected a matrix input.", call. = FALSE)
  M[lower.tri(M, diag = diag)]
}

.mppi_sym <- function(M) 0.5 * (M + t(M))

.mppi_asym <- function(M) 0.5 * (M - t(M))

.mppi_frob <- function(M) sqrt(sum(M * M, na.rm = TRUE))

.mppi_safe_pos <- function(x, eps = 1e-12) {
  x[is.na(x)] <- 0
  x[x < eps] <- eps
  x
}

.mppi_effective_rank <- function(M) {
  if (!is.matrix(M)) stop("Expected a matrix input.", call. = FALSE)
  sv <- suppressWarnings(svd(M, nu = 0L, nv = 0L)$d)
  if (!length(sv)) return(0)
  sv <- abs(sv)
  total <- sum(sv)
  if (!is.finite(total) || total <= 0) return(0)
  p <- sv / total
  p <- p[p > 0]
  if (!length(p)) return(0)
  exp(-sum(p * log(.mppi_safe_pos(p))))
}

.mppi_first_eigvec <- function(S) {
  if (!is.matrix(S)) stop("Expected a matrix input.", call. = FALSE)
  eig <- suppressWarnings(eigen(.mppi_sym(S), symmetric = TRUE))
  eig$vectors[, 1, drop = FALSE]
}

# Parametric psychological modulators ----------------------------------------

#' Build parametric psychological vectors for multi-condition gPPI
#'
#' @param fd `fmridesign` object supplying design information
#' @param T integer number of TRs (rows)
#' @param mods named list of numeric vectors (length `T`) containing trial-wise modulators
#' @param select optional regular expression to restrict conditions (columns of `psych(fd, ...)`)
#' @param center logical; demean each modulator within run
#' @param span_space logical; residualize against other conditions (gPPI default)
#' @param use_nuisance logical; include nuisance regressors when residualising
#' @param include_intercept logical; include an intercept during residualisation
#' @return named list of vectors ready to append to `pk()` output
#' @export
mppi_parametric <- function(fd, T, mods,
                            select = NULL,
                            center = TRUE,
                            span_space = TRUE,
                            use_nuisance = TRUE,
                            include_intercept = TRUE) {
  if (!inherits(fd, c("fmridesign", "event_model", "mppi_design"))) {
    stop("fd must be an fmridesign/event_model/mppi_design object.", call. = FALSE)
  }
  if (!is.list(mods) || !length(mods)) stop("mods must be a non-empty named list.", call. = FALSE)
  runs <- tryCatch(mppi_runs(fd, T = T), error = function(e) NULL)
  TR <- tryCatch(mppi_tr(fd), error = function(e) NULL)
  if (is.null(runs)) runs <- rep_len(1L, T)
  if (is.null(TR) || !is.finite(TR)) TR <- 1
  S <- mppi_psych(fd, T = T, TR = TR, runs = runs)
  if (!is.null(select)) {
    keep <- grepl(select, colnames(S))
    if (!any(keep)) stop("No psychological regressors match 'select'.", call. = FALSE)
    S <- S[, keep, drop = FALSE]
  }
  if (!ncol(S)) return(list())
  N <- if (use_nuisance) tryCatch(mppi_nuisance(fd, T = T), error = function(e) matrix(0, T, 0)) else matrix(0, T, 0)
  out <- list()
  for (mname in names(mods)) {
    v <- mods[[mname]]
    if (length(v) != T) stop(sprintf("Modulator '%s' has length %d (expected %d).", mname, length(v), T), call. = FALSE)
    v <- as.numeric(v)
    if (center) v <- .mppi_center_by_run(v, runs)
    for (j in seq_len(ncol(S))) {
      s_col <- S[, j]
      pk_vec <- s_col * v
      others <- if (span_space && ncol(S) > 1) S[, -j, drop = FALSE] else NULL
      Xres <- cbind(others, if (ncol(N)) N else NULL,
                    if (include_intercept) matrix(1, T, 1) else NULL)
      if (ncol(Xres)) {
        pk_vec <- pk_vec - qr.fitted(qr(Xres), pk_vec)
      }
      out[[paste0(colnames(S)[j], "::", mname)]] <- as.numeric(pk_vec)
    }
  }
  out
}

# Communication reinstatement -------------------------------------------------

#' Encode/recall reinstatement of task-gated connectivity patterns
#'
#' @param fit `mppi_fit` object
#' @param enc_pattern regex identifying encoding conditions
#' @param rec_pattern regex identifying recall conditions
#' @param mode interaction scale (`"normalized"`, `"amplitude"`, or `"raw"`)
#' @param id_fun function mapping condition names to identity labels
#' @return list with within-identity, between-identity similarities and their contrast
#' @export
mppi_reinstatement <- function(fit,
                               enc_pattern = "^enc_",
                               rec_pattern = "^rec_",
                               mode = c("normalized", "amplitude", "raw"),
                               id_fun = function(s) sub(".*?(\\d+)$", "\\1", s)) {
  mode <- match.arg(mode)
  if (!inherits(fit, "mppi_fit")) stop("fit must be an mppi_fit object.", call. = FALSE)
  cond_names <- fit$names %||% paste0("k", seq_along(fit$Delta))
  enc_names <- grep(enc_pattern, cond_names, value = TRUE)
  rec_names <- grep(rec_pattern, cond_names, value = TRUE)
  if (!length(enc_names) || !length(rec_names)) {
    stop("No encoding/recall conditions matched the supplied patterns.", call. = FALSE)
  }
  enc_ids <- setNames(enc_names, vapply(enc_names, id_fun, character(1)))
  rec_ids <- setNames(rec_names, vapply(rec_names, id_fun, character(1)))
  common <- intersect(names(enc_ids), names(rec_ids))
  if (!length(common)) {
    stop("No shared identities between encoding and recall conditions.", call. = FALSE)
  }
  vec_lower <- function(name) {
    Mk <- mppi_get_M_scaled(fit, name, mode = mode)
    .mppi_vec_lower(Mk, diag = FALSE)
  }
  within <- vapply(common, function(id) {
    enc_vec <- vec_lower(enc_ids[[id]])
    rec_vec <- vec_lower(rec_ids[[id]])
    stats::cor(enc_vec, rec_vec, use = "pairwise.complete.obs")
  }, numeric(1))
  between <- numeric(0)
  if (length(common) > 1) {
    comb <- utils::combn(common, 2L)
    between <- apply(comb, 2L, function(pair) {
      enc_vec <- vec_lower(enc_ids[[pair[1]]])
      rec_vec <- vec_lower(rec_ids[[pair[2]]])
      stats::cor(enc_vec, rec_vec, use = "pairwise.complete.obs")
    })
  }
  list(within = within,
       between = between,
       stat = mean(within, na.rm = TRUE) - mean(between, na.rm = TRUE))
}

# Trial-level sdPPI summaries -------------------------------------------------

#' Trial-level state-dependent PPI summaries
#'
#' @param fit `mppi_fit` object
#' @param fd `fmridesign` supplying event sticks (used to locate trials)
#' @param select optional regex restricting conditions
#' @param pre_window number of TRs before event onset for baseline features
#' @param lags integer lags for routing index (excluding zero)
#' @param mode interaction scale (`"normalized"`, `"amplitude"`, or `"raw"`)
#' @return data frame with per-trial summaries (condition, trial index, onset, magnitude, gain, routing, effective rank, v1 dominance)
#' @export
mppi_sdppi <- function(fit, fd,
                       select = NULL,
                       pre_window = 10L,
                       lags = -2:2,
                       mode = c("normalized", "amplitude", "raw")) {
  mode <- match.arg(mode)
  if (!inherits(fit, "mppi_fit")) stop("fit must be an mppi_fit object.", call. = FALSE)
  if (!inherits(fd, c("fmridesign", "event_model", "mppi_design"))) {
    stop("fd must be an fmridesign/event_model/mppi_design object.", call. = FALSE)
  }
  U <- fit$U %||% fit$R
  if (is.null(U)) stop("Trial-level summaries require the fit to store time-series (fit$U or fit$R).", call. = FALSE)
  Tn <- nrow(U)
  runs <- fit$runs %||% tryCatch(mppi_runs(fd, T = Tn), error = function(e) NULL) %||% rep_len(1L, Tn)
  if (length(runs) != Tn) stop("Run vector length mismatch.", call. = FALSE)
  TR <- fit$deconv$tr %||% fit$deconv$TR %||% tryCatch(mppi_tr(fd), error = function(e) NULL)
  if (is.null(TR) || !is.finite(TR)) TR <- 1
  S <- mppi_psych(fd, T = Tn, TR = TR, runs = runs)
  if (!is.null(select)) {
    keep <- grepl(select, colnames(S))
    if (!any(keep)) stop("No psychological regressors match 'select'.", call. = FALSE)
    S <- S[, keep, drop = FALSE]
  }
  if (!ncol(S)) stop("Design has zero matching psychological regressors.", call. = FALSE)
  sigma <- fit$sigma %||% matrixStats::colSds(fit$R)
  Sigma0 <- fit$Sigma0 %||% crossprod(U) / Tn
  v1 <- .mppi_first_eigvec(Sigma0)
  lag_vec <- setdiff(sort(unique(as.integer(lags))), 0L)
  blocklens <- fit$lag_blocklens %||% if (!is.null(fit$runs)) rle(fit$runs)$lengths else NULL
  results <- list()
  for (j in seq_len(ncol(S))) {
    stick <- S[, j]
    cond_name <- colnames(S)[j]
    starts <- which(stick > 0 & c(TRUE, diff(stick)) > 0)
    if (!length(starts)) {
      starts <- which(stick > 0 & c(1, diff(stick)) != 0)
    }
    if (!length(starts)) next
    trial_counter <- 0L
    for (s_idx in starts) {
      end_idx <- s_idx
      while (end_idx <= Tn && stick[end_idx] > 0) end_idx <- end_idx + 1L
      idx_range <- s_idx:(end_idx - 1L)
      s_event <- numeric(Tn)
      s_event[idx_range] <- stick[idx_range]
      denom <- sum(s_event^2)
      if (denom < 1e-10) next
      trial_counter <- trial_counter + 1L
      M_raw <- crossprod(U, s_event * U) / denom
      M_scaled <- if (mode == "raw") M_raw else .mppi_scale_matrix(M_raw, sigma, mode)
      Sk <- .mppi_sym(M_scaled)
      magnitude <- .mppi_frob(M_scaled)
      gain <- if (magnitude < 1e-8) NA_real_ else as.numeric((t(v1) %*% Sk %*% v1) / (.mppi_frob(Sk) + 1e-12))
      routing <- NA_real_
      if (length(lag_vec)) {
        pos <- neg <- 0
        for (lag in lag_vec) {
          Mk_lag <- .mppi_crossprod_lagged(U, s_event, lag, blocklens = blocklens)
          Mk_lag <- if (mode == "raw") Mk_lag else .mppi_scale_matrix(Mk_lag, sigma, mode)
          e <- .mppi_frob(Mk_lag)
          if (lag > 0) pos <- pos + e else neg <- neg + e
        }
        routing <- (pos - neg) / (pos + neg + 1e-12)
      }
      run_label <- runs[s_idx]
      run_idx <- which(runs == run_label)
      pre_start <- max(min(run_idx), s_idx - pre_window)
      pre_end <- max(pre_start, s_idx - 1L)
      if (pre_end - pre_start + 1L >= 2L) {
        Upre <- U[pre_start:pre_end, , drop = FALSE]
        Sp <- stats::cov(Upre)
        erank <- .mppi_effective_rank(Sp)
        v1dom <- as.numeric((t(v1) %*% Sp %*% v1) / (sum(diag(Sp)) + 1e-12))
      } else {
        erank <- NA_real_
        v1dom <- NA_real_
      }
      results[[length(results) + 1L]] <- data.frame(
        condition = cond_name,
        trial = trial_counter,
        onset = s_idx,
        magnitude = magnitude,
        gain = gain,
        routing = routing,
        effective_rank = erank,
        v1_dominance = v1dom,
        run = run_label,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(results)) {
    stop("No trials detected for the requested conditions.", call. = FALSE)
  }
  do.call(rbind, results)
}

# Precision gating -----------------------------------------------------------

#' Precision gating summary \eqn{\Delta\Theta_k = -\Theta_0 M_k \Theta_0}
#'
#' @param fit `mppi_fit` object
#' @param mode interaction scale (`"normalized"`, `"amplitude"`, or `"raw"`)
#' @param Theta0 optional baseline precision matrix; defaults to the inverse of `fit$Sigma0`
#' @param ridge small positive ridge added to the baseline covariance before inversion
#' @return list with per-condition precision deltas and a summary data frame
#' @export
mppi_precision_gate <- function(fit,
                                mode = c("normalized", "amplitude", "raw"),
                                Theta0 = NULL,
                                ridge = 1e-4) {
  mode <- match.arg(mode)
  if (!inherits(fit, "mppi_fit")) stop("fit must be an mppi_fit object.", call. = FALSE)
  Sigma0 <- fit$Sigma0
  if (is.null(Sigma0)) stop("fit does not store baseline covariance (Sigma0).", call. = FALSE)
  if (is.null(Theta0)) {
    p <- ncol(Sigma0)
    Theta0 <- solve(Sigma0 + ridge * diag(p))
  }
  cond_names <- fit$names %||% paste0("k", seq_along(fit$Delta))
  Delta_list <- vector("list", length(cond_names))
  summary_rows <- vector("list", length(cond_names))
  names(Delta_list) <- cond_names
  for (ii in seq_along(cond_names)) {
    Mk <- mppi_get_M_scaled(fit, ii, mode = mode)
    DeltaTheta <- -Theta0 %*% Mk %*% Theta0
    Delta_list[[ii]] <- DeltaTheta
    off <- DeltaTheta
    diag(off) <- 0
    summary_rows[[ii]] <- data.frame(
      condition = cond_names[ii],
      dtheta_energy = .mppi_frob(DeltaTheta),
      offdiag_ratio = sum(abs(off)) / (sum(abs(diag(DeltaTheta))) + 1e-12),
      stringsAsFactors = FALSE
    )
  }
  list(Delta = Delta_list,
       summary = do.call(rbind, summary_rows),
       Theta0 = Theta0,
       mode = mode)
}
