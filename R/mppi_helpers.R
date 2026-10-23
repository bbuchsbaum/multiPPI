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
#' Trials of each selected psychological column of `fd` start at `onsets` or,
#' by default, at rising edges of the regressor from (near) zero: the first
#' scan of a run, or after a scan at or below it, whose value exceeds `1e-3`
#' of the column's maximum. Each isolated event on a stick, boxcar or
#' convolved regressor is counted once (on a convolved regressor the onset is
#' the first scan of the response). When the responses of successive events
#' overlap, so the regressor rises again before returning to zero, events
#' cannot be told apart and the function stops; supply `onsets`, as also for
#' neural-domain fits, whose working regressors are deconvolved. A trial spans
#' its onset to the scan before the next onset of that condition in the run,
#' or the run end.
#'
#' Each trial's matrix is a covariance modulation of the series behind
#' `fit$Delta` (`fit$R`, which for `domain = "innovations"` holds the
#' innovations; `fit$Z` for basis fits): the slope of \eqn{x_t x_t^T} on the
#' fit's regressor `fit$pk[[k]]` restricted to the trial, after residualizing
#' that weight on an intercept and on the rest of the regressor (least-squares
#' separate). Data with activation but no modulation thus give matrices near
#' zero, and with an intercept in the design a trial covering every scan gives
#' `fit$Delta[[k]]`. The diagonal
#' follows the fit's `zero_diag`. `magnitude` is the Frobenius norm of the
#' matrix on the scale of `mode` and `gain` its alignment with the leading
#' eigenvector of `fit$Sigma0` (on the correlation scale unless
#' `mode = "raw"`). `routing` is the directedness of the trial's
#' lead-lag energy over the nonzero `lags`, in \eqn{[0, 1]}, as `R` in
#' [mppi_lagged()]. Pre-onset features use the `pre_window` scans before the
#' onset in the same run.
#'
#' @param fit `mppi_fit` object (uncensored)
#' @param fd `fmridesign`, `event_model` or `mppi_design` supplying the
#'   psychological regressors used to locate trials; column names must match
#'   `fit$names` (or their `"s:"`-prefixed neural-domain names)
#' @param select optional regex restricting conditions
#' @param pre_window number of TRs before event onset for baseline features
#' @param lags integer lags for routing index (excluding zero)
#' @param mode interaction scale (`"normalized"`, `"amplitude"`, or `"raw"`)
#' @param onsets optional list of onset scans (row indices of the fit) named by
#'   condition, or a vector when one condition is selected; replaces detection
#'   (conditions without an entry are skipped)
#' @return data frame with per-trial summaries (condition, trial index, onset, magnitude, gain, routing, effective rank, v1 dominance)
#' @export
mppi_sdppi <- function(fit, fd,
                       select = NULL,
                       pre_window = 10L,
                       lags = -2:2,
                       mode = c("normalized", "amplitude", "raw"),
                       onsets = NULL) {
  mode <- match.arg(mode)
  if (!inherits(fit, "mppi_fit")) stop("fit must be an mppi_fit object.", call. = FALSE)
  if (!inherits(fd, c("fmridesign", "event_model", "mppi_design"))) {
    stop("fd must be an fmridesign/event_model/mppi_design object.", call. = FALSE)
  }
  X <- .mppi_residual_matrix(fit)
  Tn <- nrow(X)
  index <- fit$observation_index %||% seq_len(Tn)
  if (length(index) != Tn || any(index != seq_len(Tn))) stop("Trial summaries need an uncensored fit.", call. = FALSE)
  runs <- fit$runs %||% tryCatch(mppi_runs(fd, T = Tn), error = function(e) NULL) %||% rep_len(1L, Tn)
  if (length(runs) != Tn || anyNA(runs)) stop("Run vector length mismatch.", call. = FALSE)
  TR <- fit$deconv$tr %||% fit$deconv$TR %||% tryCatch(mppi_tr(fd), error = function(e) NULL)
  if (is.null(TR) || !is.finite(TR)) TR <- 1
  S <- mppi_psych(fd, T = Tn, TR = TR, runs = runs)
  if (!is.null(select)) {
    keep <- grepl(select, colnames(S))
    if (!any(keep)) stop("No psychological regressors match 'select'.", call. = FALSE)
    S <- S[, keep, drop = FALSE]
  }
  if (!ncol(S)) stop("Design has zero matching psychological regressors.", call. = FALSE)
  if (!is.null(onsets) && !is.list(onsets)) {
    if (ncol(S) != 1L) stop("Supply onsets as a list named by condition when several conditions are selected.", call. = FALSE)
    onsets <- setNames(list(onsets), colnames(S))
  }
  sigma <- fit$sigma %||% .mppi_sigma_cols(X)
  Sigma0 <- fit$Sigma0 %||% crossprod(X) / Tn
  # Leading baseline mode on the scale of the trial matrices it is compared with.
  v1 <- .mppi_first_eigvec(if (mode == "raw") Sigma0 else .mppi_scale_matrix(Sigma0, sigma, "normalized"))
  lag_vec <- setdiff(sort(unique(as.integer(lags))), 0L)
  scale_m <- function(M) if (mode == "raw") M else .mppi_scale_matrix(M, sigma, mode)
  results <- list()
  for (j in seq_len(ncol(S))) {
    cond_name <- colnames(S)[j]
    k <- match(cond_name, fit$names)
    if (is.na(k)) k <- match(paste0("s:", cond_name), fit$names)
    if (is.na(k)) stop(sprintf("Design column '%s' matches no regressor in fit$names.", cond_name), call. = FALSE)
    w <- as.numeric(fit$pk[[k]])
    if (is.null(onsets)) {
      starts <- .mppi_trial_onsets(S[, j], runs, cond_name)
    } else {
      starts <- onsets[[cond_name]]
      if (is.null(starts)) next
      starts <- sort(unique(as.integer(starts)))
      if (anyNA(starts) || any(starts < 1L | starts > Tn)) stop("onsets must be scan indices within 1:T.", call. = FALSE)
    }
    trial_counter <- 0L
    for (s_idx in starts) {
      run_label <- runs[s_idx]
      run_idx <- which(runs == run_label)
      nxt <- starts[starts > s_idx & runs[starts] == run_label]
      in_trial <- seq_len(Tn) %in% seq.int(s_idx, if (length(nxt)) min(nxt) - 1L else max(run_idx))
      pw <- .mppi_residualize_vec(w * in_trial, cbind(1, w * !in_trial))
      denom <- sum(pw^2)
      if (denom < 1e-10) next
      trial_counter <- trial_counter + 1L
      M_raw <- crossprod(X, pw * X) / denom
      if (isTRUE(fit$zero_diag)) diag(M_raw) <- 0
      M_scaled <- scale_m(M_raw)
      Sk <- .mppi_sym(M_scaled)
      magnitude <- .mppi_frob(M_scaled)
      gain <- if (magnitude < 1e-8) NA_real_ else as.numeric((t(v1) %*% Sk %*% v1) / (.mppi_frob(Sk) + 1e-12))
      routing <- if (!length(lag_vec)) NA_real_ else
        .mppi_directedness(.mppi_lead_energy(lapply(.mppi_adv_lag_mats(X, pw, lag_vec, runs, index), scale_m), lag_vec))
      pre_start <- max(min(run_idx), s_idx - pre_window)
      pre_end <- s_idx - 1L
      if (pre_end - pre_start + 1L >= 2L) {
        Spre <- stats::cov(X[pre_start:pre_end, , drop = FALSE])
        erank <- .mppi_effective_rank(Spre)
        v1dom <- as.numeric((t(v1) %*% Spre %*% v1) / (sum(diag(Spre)) + 1e-12))
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

# Trial onsets on a design regressor: rising edges from (near) zero within runs;
# stops when a response rises again before returning to zero (overlapping events).
.mppi_trial_onsets <- function(s, runs, name) {
  tol <- 1e-3 * max(abs(s))
  if (!is.finite(tol) || tol <= 0) return(integer(0))
  n <- length(s)
  first <- c(TRUE, runs[-1L] != runs[-n])
  on <- s > tol
  starts <- which(on & (first | c(TRUE, s[-n] <= tol)))
  for (st in starts) {
    e <- st
    while (e < n && on[e + 1L] && !first[e + 1L]) e <- e + 1L
    v <- s[st:e]
    m <- which.max(v)
    if (any(diff(v[seq_len(m)]) < -tol) || any(diff(v[m:length(v)]) > tol)) {
      stop(sprintf("Responses of successive '%s' events overlap, so trials cannot be located on the regressor; supply `onsets`.", name), call. = FALSE)
    }
  }
  starts
}

# Precision gating -----------------------------------------------------------

#' Precision gating summary \eqn{\Delta\Theta_k = -\Theta_0 M_k \Theta_0}
#'
#' \eqn{M_k} is the full slope (diagonal kept, whatever `zero_diag` was), as in
#' `fit$partial`. In `"raw"` mode \eqn{\Sigma_0} and \eqn{M_k} are in data units
#' and, with the default ridge, \eqn{\Delta\Theta_k} equals
#' `fit$partial[[k]]$DeltaTheta` of a `scale = "cov"` fit. In `"normalized"` mode
#' both are rescaled to correlation scale by `fit$sigma`
#' (\eqn{D^{-1} \Sigma_0 D^{-1}}, \eqn{D^{-1} M_k D^{-1}}), so the result is unitless.
#' `"amplitude"` slopes are not symmetric covariances and are rejected.
#'
#' @param fit `mppi_fit` object
#' @param mode interaction scale (`"normalized"` or `"raw"`; `"amplitude"` errors)
#' @param Theta0 optional baseline precision matrix on the scale of `mode`;
#'   defaults to the ridged inverse of `fit$Sigma0` on that scale
#' @param ridge ridge added to the baseline covariance before inversion; defaults
#'   to `fit$partial_lambda`, the ridge used for `fit$partial`
#' @return list with per-condition precision deltas and a summary data frame
#' @export
mppi_precision_gate <- function(fit,
                                mode = c("normalized", "amplitude", "raw"),
                                Theta0 = NULL,
                                ridge = NULL) {
  mode <- match.arg(mode)
  if (!inherits(fit, "mppi_fit")) stop("fit must be an mppi_fit object.", call. = FALSE)
  if (mode == "amplitude") stop("mode = \"amplitude\" has no precision-domain analogue; use \"raw\" or \"normalized\".", call. = FALSE)
  Sigma0 <- fit$Sigma0
  if (is.null(Sigma0)) stop("fit does not store baseline covariance (Sigma0).", call. = FALSE)
  if (is.null(fit$covariance_diagonal)) stop("fit does not store the slope diagonal (covariance_diagonal); refit with mppi_fit().", call. = FALSE)
  if (mode == "normalized") {
    if (is.null(fit$sigma)) stop("fit does not store residual SDs (sigma) for the normalized scale.", call. = FALSE)
    Sigma0 <- .mppi_scale_matrix(Sigma0, fit$sigma, "normalized")
  }
  if (is.null(Theta0)) {
    ridge <- ridge %||% fit$partial_lambda %||% 1e-3
    p <- ncol(Sigma0)
    Theta0 <- solve(Sigma0 + ridge * diag(p))
  }
  cond_names <- fit$names %||% paste0("k", seq_along(fit$Delta))
  Delta_list <- vector("list", length(cond_names))
  summary_rows <- vector("list", length(cond_names))
  names(Delta_list) <- cond_names
  for (ii in seq_along(cond_names)) {
    Mk <- mppi_get_M_scaled(fit, ii, mode = "raw")
    diag(Mk) <- fit$covariance_diagonal[[ii]]
    if (mode == "normalized") Mk <- .mppi_scale_matrix(Mk, fit$sigma, "normalized")
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
