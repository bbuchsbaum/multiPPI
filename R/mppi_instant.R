# Instantaneous engine ------------------------------------------------------

.mppi_instant_defaults <- function() {
  list(
    method = "ewm",
    tau_half = 10,
    offsets = 0L,
    band = c(0.03, 0.07),
    normalize = "normalized",
    warmup = NULL
  )
}

.mppi_resolve_instant_options <- function(instant, scale_choice) {
  defaults <- .mppi_instant_defaults()
  if (is.null(instant)) instant <- list()
  opts <- modifyList(defaults, instant, keep.null = TRUE)
  opts$method <- match.arg(opts$method, c("ewm", "iphs"))
  opts$normalize <- match.arg(opts$normalize, c("normalized", "amplitude", "raw"))
  if (is.null(opts$offsets)) opts$offsets <- 0L
  opts$offsets <- as.integer(opts$offsets)
  if (!0L %in% opts$offsets) opts$offsets <- c(0L, opts$offsets)
  opts$offsets <- sort(unique(opts$offsets))
  if (opts$method == "ewm" && (is.null(opts$tau_half) || !is.finite(opts$tau_half) || opts$tau_half <= 0)) {
    stop("instant$tau_half must be a positive value for method='ewm'.", call. = FALSE)
  }
  if (opts$method == "iphs") {
    if (length(opts$band) != 2L || any(!is.finite(opts$band))) {
      stop("instant$band must be a numeric length-2 vector when method='iphs'.", call. = FALSE)
    }
    opts$band <- sort(as.numeric(opts$band))
    if (diff(opts$band) <= 0) stop("instant$band must have band[1] < band[2].", call. = FALSE)
  }
  opts$scale_choice <- scale_choice
  opts
}

.mppi_instant_prewhiten <- function(Y, X, runs, ar_method = "ar", ar_order = "auto") {
  if (!requireNamespace("fmriAR", quietly = TRUE)) {
    stop("fmriAR not available. Install via remotes::install_github('bbuchsbaum/fmriAR').")
  }
  stopifnot(length(runs) == nrow(Y), nrow(X) == nrow(Y))
  qrX <- qr(X)
  res <- Y - X %*% qr.coef(qrX, Y)
  plan <- fmriAR::fit_noise(res, runs = runs, method = ar_method, p = ar_order)
  fmriAR::whiten_apply(plan, X = X, Y = Y, runs = runs)
}

#' Instantaneous mPPI engine
#'
#' @keywords internal
mppi_instant_fit <- function(Y, X, psych_idx, runs = NULL,
                             prewhiten = FALSE,
                             basis = NULL,
                             domain = c("bold", "neural", "innovations"),
                             deconv = NULL,
                             center_by = c("none", "run"),
                             na_action = c("omit_tr", "error"),
                             scale = c("normalized", "cov", "corr"),
                             ar_method = "ar", ar_order = "auto",
                             instant = list(), ...) {
  stopifnot(is.matrix(Y), is.matrix(X))
  domain <- match.arg(domain)
  center_by <- match.arg(center_by)
  na_action_choice <- match.arg(na_action)
  scale_choice <- match.arg(scale)
  opts <- .mppi_resolve_instant_options(instant, scale_choice)

  if (prewhiten) {
    if (is.null(runs)) stop("'runs' must be supplied when prewhiten = TRUE.", call. = FALSE)
    xyw <- .mppi_instant_prewhiten(Y, X, runs, ar_method = ar_method, ar_order = ar_order)
    Y <- xyw$Y
    X <- xyw$X
  }

  state <- .mppi_engine_state(Y, X, psych_idx, runs = runs,
                              basis = basis,
                              domain = domain,
                              deconv = deconv,
                              center_by = center_by,
                              na_action = na_action_choice,
                              project_backend = "blas",
                              project_chunk_cols = NULL)

  core <- .mppi_instant_core(state, opts, ...)

  .mppi_instant_build_fit(state, core, opts, prewhiten = prewhiten)
}

.mppi_instant_core <- function(state, opts, ...) {
  method <- opts$method
  pk_mat <- do.call(cbind, state$pk_list)
  colnames(pk_mat) <- state$pk_names
  sigma <- state$sigma
  scale_choice <- opts$scale_choice

  if (method == "ewm") {
    offsets <- opts$offsets
    if (length(offsets) == 0L) offsets <- 0L
    blk <- if (is.null(state$runs)) NULL else rle(state$runs)$lengths
    ewm <- ipppi_ewm(state$R_work, pk_mat,
                     tau_half = opts$tau_half,
                     offsets = offsets,
                     normalized = FALSE,
                     blocklens = blk)
    offsets <- ewm$offsets
    zero_idx <- match(0L, offsets)
    if (is.na(zero_idx)) stop("EWM engine requires offsets to include 0.", call. = FALSE)
    lagged <- vector("list", length(state$pk_list))
    Delta_raw <- vector("list", length(state$pk_list))
    Delta_amp <- vector("list", length(state$pk_list))
    Delta_norm <- vector("list", length(state$pk_list))
    names(lagged) <- names(Delta_raw) <- names(Delta_amp) <- names(Delta_norm) <- state$pk_names
    for (k in seq_along(ewm$mats)) {
      offset_list <- ewm$mats[[k]]
      lag_sel <- vector("list", length(offset_list))
      names(lag_sel) <- names(offset_list)
      for (nm in names(offset_list)) {
        raw_mat <- offset_list[[nm]]
        amp_mat <- if (all(is.na(raw_mat))) raw_mat else .mppi_scale_matrix(raw_mat, sigma, "amplitude")
        norm_mat <- if (all(is.na(raw_mat))) raw_mat else .mppi_scale_matrix(raw_mat, sigma, "normalized")
        lag_sel[[nm]] <- .mppi_instant_select_scale(raw_mat, amp_mat, norm_mat, scale_choice)
        if (nm == "0") {
          Delta_raw[[k]] <- raw_mat
          Delta_amp[[k]] <- amp_mat
          Delta_norm[[k]] <- norm_mat
        }
      }
      lagged[[k]] <- lag_sel
    }
    return(list(Delta_raw = Delta_raw,
                Delta_amplitude = Delta_amp,
                Delta_normalized = Delta_norm,
                lagged = lagged,
                offsets = offsets,
                lag_blocklens = blk))
  }

  fs <- opts$fs
  if (is.null(fs)) {
    tr_attr <- attr(state$R_work, "TR") %||% state$deconv$tr
    if (!is.null(tr_attr) && is.finite(tr_attr) && tr_attr > 0) {
      fs <- 1 / tr_attr
    }
  }
  if (is.null(fs)) fs <- 1
  band <- opts$band
  tau_half <- opts$tau_half %||% 10
  iphs_list <- ipppi_iphs(state$R_work, pk_mat, fs = fs, band = band, tau_half = tau_half)
  lagged <- vector("list", length(iphs_list))
  Delta_raw <- vector("list", length(iphs_list))
  Delta_amp <- vector("list", length(iphs_list))
  Delta_norm <- vector("list", length(iphs_list))
  names(lagged) <- names(Delta_raw) <- names(Delta_amp) <- names(Delta_norm) <- state$pk_names
  for (k in seq_along(iphs_list)) {
    corr_mat <- iphs_list[[k]]
    if (!is.null(sigma)) {
      diag_sigma <- diag(sigma, length(sigma))
      raw_mat <- diag_sigma %*% corr_mat %*% diag_sigma
    } else {
      raw_mat <- corr_mat
    }
    amp_mat <- if (is.null(sigma)) raw_mat else .mppi_scale_matrix(raw_mat, sigma, "amplitude")
    norm_mat <- if (is.null(sigma)) corr_mat else .mppi_scale_matrix(raw_mat, sigma, "normalized")
    selected <- .mppi_instant_select_scale(raw_mat, amp_mat, norm_mat, scale_choice)
    lagged[[k]] <- list(`0` = selected)
    Delta_raw[[k]] <- raw_mat
    Delta_amp[[k]] <- amp_mat
    Delta_norm[[k]] <- norm_mat
  }
  list(Delta_raw = Delta_raw,
       Delta_amplitude = Delta_amp,
       Delta_normalized = Delta_norm,
       lagged = lagged,
       offsets = 0L,
       lag_blocklens = NULL)
}

.mppi_instant_select_scale <- function(raw, amp, norm, scale_choice) {
  switch(scale_choice,
         normalized = norm,
         cov = raw,
         corr = norm)
}

.mppi_instant_build_fit <- function(state, core, opts, prewhiten) {
  scale_choice <- opts$scale_choice
  Delta_raw <- core$Delta_raw
  Delta_amplitude <- core$Delta_amplitude
  Delta_normalized <- core$Delta_normalized
  Delta_main <- .mppi_instant_select_scale(Delta_raw, Delta_amplitude, Delta_normalized, scale_choice)
  names(Delta_main) <- state$pk_names
  names(Delta_raw) <- state$pk_names
  names(Delta_amplitude) <- state$pk_names
  names(Delta_normalized) <- state$pk_names

  res <- list(
    Delta = Delta_main,
    names = state$pk_names,
    R = state$R_work,
    R_raw = state$R_unscaled,
    Delta_raw = Delta_raw,
    Delta_amplitude = Delta_amplitude,
    Delta_normalized = Delta_normalized,
    sigma = state$sigma,
    pk = state$pk_list,
    scale = scale_choice,
    dropped = state$dropped,
    center_by = state$center_by,
    runs = state$runs,
    denom = state$denoms,
    backend = "instant",
    n_used = nrow(state$R_work),
    Sigma0 = state$Sigma0,
    partial_lambda = NA_real_,
    Theta0 = NULL,
    variance = NULL,
    partial = NULL,
    lags = core$offsets,
    lagged = core$lagged,
    lag_blocklens = core$lag_blocklens %||% NULL,
    AIC_total = rep(state$evidence$AIC_full_total, length(state$pk_list)),
    BIC_total = rep(state$evidence$BIC_full_total, length(state$pk_list)),
    packed = FALSE,
    chunk_size = NULL,
    basis = state$basis,
    project_backend = state$project_backend,
    project_chunk_cols = state$project_chunk_cols,
    domain = state$domain,
    deconv = state$deconv,
    evidence = state$evidence,
    engine = "instant",
    instant = opts,
    prewhiten = prewhiten
  )

  if (!is.null(state$basis)) res$Z <- state$R_work
  if (identical(state$domain, "neural") || identical(state$domain, "innovations")) res$U <- state$R_unscaled
  names(res$AIC_total) <- state$pk_names
  names(res$BIC_total) <- state$pk_names
  class(res) <- c("mppi_fit", "mppi_instant", "list")
  res
}
