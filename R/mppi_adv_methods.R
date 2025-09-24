# Advanced mechanistic helpers ------------------------------------------------

#' Protected-band innovations (iPPI companion)
#'
#' @param U matrix of neural/basis time-series (T x r)
#' @param TR repetition time in seconds (needed for protect = "lowHz")
#' @param protect slow-band strategy: "none", "lowK", or "lowHz"
#' @param K_keep number of DCT components to retain when protect = "lowK"
#' @param f_cut low-frequency cutoff in Hz when protect = "lowHz"
#' @param ar_order "auto" (AIC) or "fixed"
#' @param p maximum AR order when ar_order = "auto" or fixed order otherwise
#' @param center logical; center columns before whitening
#' @return matrix of innovations with attributes `orders` and `phis`
#' @export
mppi_innovations <- function(U, TR = NULL,
                             protect = c("none", "lowK", "lowHz"),
                             K_keep = 8L, f_cut = 0.03,
                             ar_order = c("auto", "fixed"), p = 4L,
                             center = TRUE) {
  stopifnot(is.matrix(U))
  protect <- match.arg(protect)
  ar_order <- match.arg(ar_order)
  Tn <- nrow(U); r <- ncol(U)
  if (center) U <- scale(U, center = TRUE, scale = FALSE)
  if (protect == "lowK") {
    U_low <- .mppi_adv_project_lowK(U, K_keep)
    U_high <- U - U_low
  } else if (protect == "lowHz") {
    if (is.null(TR)) stop("TR required for protect='lowHz'.")
    U_low <- .mppi_adv_bandlimit(U, TR, f_cut, which = "low")
    U_high <- U - U_low
  } else {
    U_low <- 0 * U
    U_high <- U
  }
  E_high <- matrix(0, Tn, r)
  orders <- integer(r)
  phis <- vector("list", r)
  for (j in seq_len(r)) {
    x <- U_high[, j]
    if (ar_order == "auto") {
      fit <- try(stats::ar(x, aic = TRUE, order.max = max(1L, p), method = "yw"), silent = TRUE)
      if (inherits(fit, "try-error")) {
        orders[j] <- 0L
        phi <- numeric(0)
      } else {
        orders[j] <- fit$order
        phi <- if (orders[j] > 0) as.numeric(fit$ar) else numeric(0)
      }
    } else {
      orders[j] <- as.integer(p)
      if (p > 0) {
        fit <- try(stats::ar(x, aic = FALSE, order.max = p, method = "yw"), silent = TRUE)
        phi <- if (inherits(fit, "try-error")) rep(0, p) else as.numeric(fit$ar)
      } else {
        phi <- numeric(0)
      }
    }
    phis[[j]] <- phi
    E_high[, j] <- .mppi_adv_ar_resid(x, phi)
  }
  E <- U_low + E_high
  attr(E, "orders") <- orders
  attr(E, "phis") <- phis
  E
}

#' Spectral summary of task-modulated covariance
#'
#' @param fit object returned by `mppi_fit` (basis or ROI space)
#' @param k condition index (1-based)
#' @param templates optional named list of template vectors for eigenvector alignment
#' @return list with eigenvalues, first eigenvector, lambda1 share, effective rank, etc.
#' @export
mppi_spectral_summary <- function(fit, k, templates = NULL) {
  Mk <- .mppi_adv_Mk(fit, k)
  Mk <- 0.5 * (Mk + t(Mk))
  eig <- eigen(Mk, symmetric = TRUE)
  lam <- eig$values
  V <- eig$vectors
  lam2 <- lam^2
  total <- sum(lam2) + 1e-12
  share1 <- lam2[1] / total
  eff_rank <- exp(-sum((lam2 / total) * log((lam2 / total) + 1e-12)))
  v1 <- V[, 1]
  sign_balance <- max(mean(v1 >= 0), mean(v1 <= 0))
  net_proj <- NULL
  if (!is.null(templates)) {
    net_proj <- vapply(templates, function(tmpl) {
      u <- as.numeric(tmpl)
      u <- u / sqrt(sum(u^2) + 1e-12)
      sum(v1 * u) / sqrt(sum(v1^2) + 1e-12)
    }, numeric(1))
  }
  list(lambda = lam, v1 = v1, lambda1_share = share1,
       effrank = eff_rank, sign_balance = sign_balance,
       net_proj = net_proj)
}

#' Precision-domain perturbation via baseline gating
#'
#' @param fit multiPPI fit
#' @param k condition index
#' @param lambda ridge parameter for baseline precision
#' @return list with `dTheta`, `Theta0`, and off-diagonal energy
#' @export
mppi_precision_delta <- function(fit, k, lambda = 1e-3) {
  X <- .mppi_adv_X(fit)
  S0 <- crossprod(X) / nrow(X)
  p <- ncol(S0)
  Theta0 <- solve(S0 + diag(lambda, p))
  Mk <- .mppi_adv_Mk(fit, k, X)
  dTheta <- -Theta0 %*% Mk %*% Theta0
  off <- dTheta - diag(diag(dTheta))
  list(dTheta = dTheta, Theta0 = Theta0,
       offdiag_energy = sqrt(sum(off^2)))
}

#' Short-lag task-gated coupling matrices
#'
#' @param fit multiPPI fit
#' @param k condition index
#' @param lags integer vector of TR lags
#' @param mode "neural" uses stored series; "innov" applies `mppi_innovations`
#' @param protect strategy for innovations
#' @param K_keep DCT components when protect = "lowK"
#' @param f_cut Hz cutoff when protect = "lowHz"
#' @param TR repetition time for protect = "lowHz"
#' @param ar_order AR order mode for innovations
#' @param p maximum AR order
#' @return list with matrices per lag and routing index R
#' @rdname mppi_lagged
#' @method mppi_lagged mppi_fit
#' @export
mppi_lagged.mppi_fit <- function(fit, k, lags = -2:2, mode = c("neural", "innov"),
                                 protect = "lowK", K_keep = 8L, f_cut = 0.03,
                                 TR = NULL, ar_order = "auto", p = 4L) {
  mode <- match.arg(mode)
  X <- .mppi_adv_X(fit)
  if (!is.null(attr(X, "TR"))) TR <- attr(X, "TR")
  pk <- .mppi_adv_pk(fit, k)
  if (mode == "innov") {
    X <- mppi_innovations(X, TR = TR, protect = protect,
                          K_keep = K_keep, f_cut = f_cut,
                          ar_order = ar_order, p = p, center = FALSE)
  }
  M_by_lag <- setNames(vector("list", length(lags)), as.character(lags))
  pos_energy <- 0
  neg_energy <- 0
  for (ii in seq_along(lags)) {
    lg <- as.integer(lags[ii])
    Mk <- .mppi_adv_lagged_cross(X, pk, lg)
    M_by_lag[[ii]] <- Mk
    e <- sqrt(sum(Mk^2, na.rm = TRUE))
    if (lg > 0) pos_energy <- pos_energy + e
    if (lg < 0) neg_energy <- neg_energy + e
  }
  denom <- pos_energy + neg_energy + 1e-12
  R <- (pos_energy - neg_energy) / denom
  list(M_lag = M_by_lag, R = R,
       energy_pos = pos_energy, energy_neg = neg_energy)
}

#' Gain vs routing classifier (G/C/I/R indices)
#'
#' @param fit multiPPI fit
#' @param k condition index
#' @param templates optional list of template vectors
#' @param lambda ridge parameter for precision delta
#' @param lags lags for routing index
#' @param mode innovations mode passed to `mppi_lagged`
#' @param B number of permutations for null distribution
#' @param seed RNG seed
#' @return list with indices, p-values, and label
#' @export
mppi_classify <- function(fit, k, templates = NULL, lambda = 1e-3,
                          lags = -2:2, mode = c("neural", "innov"),
                          B = 200L, seed = 1L) {
  mode <- match.arg(mode)
  spec <- mppi_spectral_summary(fit, k, templates)
  prec <- mppi_precision_delta(fit, k, lambda)
  lag <- mppi_lagged(fit, k, lags = lags, mode = mode)
  X <- .mppi_adv_X(fit)
  pk <- .mppi_adv_pk(fit, k)
  Mk_neural <- crossprod(X, pk * X) / (sum(pk^2) + 1e-12)
  Mk_innov <- if (mode == "innov") {
    Xi <- mppi_innovations(X, TR = attr(X, "TR") %||% NULL,
                           protect = "lowK", K_keep = 8L,
                           ar_order = "auto", p = 4L,
                           center = FALSE)
    crossprod(Xi, pk * Xi) / (sum(pk^2) + 1e-12)
  } else Mk_neural
  I_ratio <- sqrt(sum(Mk_innov^2)) / (sqrt(sum(Mk_neural^2)) + 1e-12)
  set.seed(seed)
  perm <- .mppi_adv_perm_shift(pk, B)
  G_null <- C_null <- I_null <- R_null <- numeric(B)
  TR_attr <- attr(X, "TR") %||% NULL
  for (b in seq_len(B)) {
    sb <- perm[, b]
    Mb <- crossprod(X, sb * X) / (sum(sb^2) + 1e-12)
    lam_b <- eigen(0.5 * (Mb + t(Mb)), symmetric = TRUE)$values
    lam2_b <- lam_b^2
    G_null[b] <- lam2_b[1] / (sum(lam2_b) + 1e-12)
    dTheta_b <- -prec$Theta0 %*% Mb %*% prec$Theta0
    off_b <- dTheta_b - diag(diag(dTheta_b))
    C_null[b] <- sqrt(sum(off_b^2))
    if (mode == "innov") {
      Xi_b <- mppi_innovations(X, TR = TR_attr, protect = "lowK",
                               K_keep = 8L, ar_order = "auto",
                               p = 4L, center = FALSE)
      Mb_i <- crossprod(Xi_b, sb * Xi_b) / (sum(sb^2) + 1e-12)
      I_null[b] <- sqrt(sum(Mb_i^2)) / (sqrt(sum(Mb^2)) + 1e-12)
    } else {
      I_null[b] <- 1
    }
    pos <- neg <- 0
    for (lg in lags) {
      Ml <- .mppi_adv_lagged_cross(X, sb, lg)
      e <- sqrt(sum(Ml^2, na.rm = TRUE))
      if (lg > 0) pos <- pos + e
      if (lg < 0) neg <- neg + e
    }
    R_null[b] <- (pos - neg) / (pos + neg + 1e-12)
  }
  G_obs <- spec$lambda1_share
  C_obs <- prec$offdiag_energy
  R_obs <- lag$R
  pG <- mean(G_null >= G_obs)
  pC <- mean(C_null >= C_obs)
  pI <- mean(I_null >= I_ratio)
  pR <- mean(abs(R_null) >= abs(R_obs))
  label <- if (G_obs > 0.5 && pC > 0.05 && pI > 0.05 && pR > 0.05) {
    "Gain-dominant"
  } else if (G_obs <= 0.5 && pC <= 0.05 && pI <= 0.05 && pR <= 0.05) {
    "Routing-dominant"
  } else {
    "Mixed/Ambiguous"
  }
  list(G = G_obs, pG = pG,
       C = C_obs, pC = pC,
       I = I_ratio, pI = pI,
       R = R_obs, pR = pR,
       label = label,
       spectrum = spec, precision = prec, lagged = lag)
}

#' Trial-resolved state-dependent PPI analysis
#'
#' @param fit multiPPI fit (neural or BOLD domain)
#' @param k condition index
#' @param onsets_idx integer TR indices for events
#' @param runs run label per TR
#' @param Wpre number of TRs in pre-event window
#' @param Wpost number of TRs in post-event window
#' @param mode "neural" or "innov"
#' @param TR repetition time (needed when innovations protect = "lowHz")
#' @param protect innovations strategy
#' @param K_keep DCT components when protect = "lowK"
#' @param f_cut Hz cutoff when protect = "lowHz"
#' @param ar_order innovations AR mode
#' @param p maximum AR order
#' @param nuisance_ts named list of nuisance timeseries (length T)
#' @return list with trial data, models, and cross-validated correlation
#' @export
mppi_sdppi_fit <- function(fit, k, onsets_idx, runs,
                           Wpre = 8L, Wpost = 12L,
                           mode = c("neural", "innov"),
                           TR = NULL, protect = "lowK", K_keep = 8L,
                           f_cut = 0.03, ar_order = "auto", p = 4L,
                           nuisance_ts = NULL) {
  mode <- match.arg(mode)
  X <- .mppi_adv_X(fit)
  if (!is.null(attr(X, "TR"))) TR <- attr(X, "TR")
  stopifnot(length(runs) == nrow(X))
  if (mode == "innov") {
    X <- mppi_innovations(X, TR = TR, protect = protect,
                          K_keep = K_keep, f_cut = f_cut,
                          ar_order = ar_order, p = p, center = FALSE)
  }
  Tn <- nrow(X)
  r <- ncol(X)
  Theta0 <- solve(crossprod(X) / Tn + diag(1e-3, r))
  rows <- list()
  out_runs <- integer(0)
  for (idx in seq_along(onsets_idx)) {
    t0 <- as.integer(onsets_idx[idx])
    pre_idx <- .mppi_adv_window(Tn, t0, Wpre, pre = TRUE)
    post_idx <- .mppi_adv_window(Tn, t0, Wpost, pre = FALSE)
    if (length(pre_idx) < 3 || length(post_idx) < 3) next
    Sw <- crossprod(X[pre_idx, , drop = FALSE]) / length(pre_idx)
    eig <- eigen(0.5 * (Sw + t(Sw)), symmetric = TRUE)
    lam <- eig$values
    lam2 <- lam^2
    total <- sum(lam2) + 1e-12
    pre_lambda1_share <- lam2[1] / total
    pre_effrank <- exp(-sum((lam2 / total) * log((lam2 / total) + 1e-12)))
    v1 <- eig$vectors[, 1]
    pre_sign_balance <- max(mean(v1 >= 0), mean(v1 <= 0))
    Mk_post <- crossprod(X, (.mppi_adv_indicator(Tn, post_idx) * X)) /
      (length(post_idx) + 1e-12)
    dTheta_post <- -Theta0 %*% Mk_post %*% Theta0
    Mk_pre <- crossprod(X, (.mppi_adv_indicator(Tn, pre_idx) * X)) /
      (length(pre_idx) + 1e-12)
    dTheta_pre <- -Theta0 %*% Mk_pre %*% Theta0
    ppi_energy <- sqrt(sum(Mk_post^2))
    cond_energy <- sqrt(sum((dTheta_post - diag(diag(dTheta_post)))^2))
    ppi_delta <- ppi_energy - sqrt(sum(Mk_pre^2))
    cond_delta <- cond_energy - sqrt(sum((dTheta_pre - diag(diag(dTheta_pre)))^2))
    nuis <- .mppi_adv_nuisance(nuisance_ts, pre_idx)
    row <- c(pre_lambda1_share = pre_lambda1_share,
             pre_effrank = pre_effrank,
             pre_sign_balance = pre_sign_balance,
             ppi_energy = ppi_energy,
             cond_energy = cond_energy,
             ppi_energy_delta = ppi_delta,
             cond_energy_delta = cond_delta,
             nuis)
    rows[[length(rows) + 1L]] <- row
    out_runs <- c(out_runs, runs[t0])
  }
  if (!length(rows)) stop("No valid trials for sdPPI analysis.")
  mat <- do.call(rbind, rows)
  data <- as.data.frame(mat)
  data$run <- out_runs
  resp_vars <- c("ppi_energy", "cond_energy", "ppi_energy_delta", "cond_energy_delta")
  predictors <- setdiff(colnames(data), c("run", resp_vars))
  lm_ppi <- stats::lm(ppi_energy ~ ., data = data[, c("ppi_energy", predictors), drop = FALSE])
  lm_cond <- stats::lm(cond_energy ~ ., data = data[, c("cond_energy", predictors), drop = FALSE])
  lm_ppi_delta <- stats::lm(ppi_energy_delta ~ ., data = data[, c("ppi_energy_delta", predictors), drop = FALSE])
  lm_cond_delta <- stats::lm(cond_energy_delta ~ ., data = data[, c("cond_energy_delta", predictors), drop = FALSE])
  cv_r <- .mppi_adv_cv_r(data, runs = data$run, response = "ppi_energy", predictors = predictors)
  list(data = data,
       model_ppi = lm_ppi,
       model_cond = lm_cond,
       model_ppi_delta = lm_ppi_delta,
       model_cond_delta = lm_cond_delta,
       cv_r = cv_r)
}

#' Generate matched null onsets away from real events
#'
#' @param Tlen number of TRs
#' @param onsets_idx observed event indices
#' @param n number of null events
#' @param guard guard width in TRs
#' @return sorted integer vector of null onset indices
#' @export
mppi_make_null_onsets <- function(Tlen, onsets_idx, n = length(onsets_idx), guard = 8L) {
  forbid <- rep(FALSE, Tlen)
  for (t0 in onsets_idx) {
    lo <- max(1L, t0 - guard)
    hi <- min(Tlen, t0 + guard)
    forbid[lo:hi] <- TRUE
  }
  cand <- which(!forbid)
  if (length(cand) < n) stop("Not enough candidate TRs for null onsets; reduce guard.")
  sort(sample(cand, n, replace = FALSE))
}

#' Cross-validated mapping from pre-state predictors to anchors
#'
#' @param sd_df output data frame from `mppi_sdppi_fit`
#' @param y response vector (numeric or binary factor)
#' @param runs run ids for cross-validation
#' @param family "gaussian" for lm, "binomial" for logistic
#' @return list with cv_r (gaussian) or cv_auc (binomial)
#' @export
mppi_state_predict <- function(sd_df, y, runs = NULL,
                               family = c("gaussian", "binomial")) {
  family <- match.arg(family)
  X <- sd_df[, grepl("^pre_", names(sd_df)), drop = FALSE]
  if (is.null(runs)) runs <- sd_df$run %||% rep(1L, nrow(sd_df))
  runs <- as.integer(runs)
  uniq <- sort(unique(runs))
  preds <- numeric(0)
  obs <- y[0]
  for (rr in uniq) {
    idx_tr <- which(runs != rr)
    idx_te <- which(runs == rr)
    if (length(idx_tr) < 5 || length(idx_te) < 1) next
    dat_tr <- data.frame(y = y[idx_tr], X[idx_tr, , drop = FALSE])
    dat_te <- data.frame(y = y[idx_te], X[idx_te, , drop = FALSE])
    if (family == "gaussian") {
      fit <- stats::lm(y ~ ., data = dat_tr)
      mm <- model.matrix(~ ., data = dat_te)[, -1, drop = FALSE]
      preds <- c(preds, as.numeric(coef(fit)[1] + mm %*% coef(fit)[-1]))
      obs <- c(obs, y[idx_te])
    } else {
      fit <- stats::glm(y ~ ., family = stats::binomial(), data = dat_tr)
      preds <- c(preds, as.numeric(stats::predict(fit, newdata = dat_te, type = "response")))
      obs <- c(obs, y[idx_te])
    }
  }
  if (family == "gaussian") {
    cv_r <- if (length(preds) > 1) stats::cor(preds, obs) else NA_real_
    list(cv_r = cv_r)
  } else {
    if (length(unique(obs)) != 2) return(list(cv_auc = NA_real_))
    ranks <- rank(preds)
    pos <- obs == max(obs)
    n_pos <- sum(pos)
    n_neg <- sum(!pos)
    if (n_pos == 0 || n_neg == 0) return(list(cv_auc = NA_real_))
    auc <- (sum(ranks[pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    list(cv_auc = auc)
  }
}

#' Baseline precision metrics
#'
#' @param Theta0 precision matrix (positive definite)
#' @return list with algebraic connectivity, effective rank, Qmax, and mean communicability
#' @export
mppi_theta_metrics <- function(Theta0) {
  p <- nrow(Theta0)
  A <- -Theta0
  diag(A) <- 0
  A[A < 0] <- 0
  deg <- rowSums(A)
  L <- diag(deg, p) - A
  eigL <- eigen(0.5 * (L + t(L)), symmetric = TRUE)$values
  alg_conn <- if (length(eigL) >= 2) eigL[2] else NA_real_
  eigA <- eigen(0.5 * (A + t(A)), symmetric = TRUE)$values
  eigA[eigA < 0] <- 0
  lam2 <- eigA^2
  total <- sum(lam2) + 1e-12
  effrankA <- exp(-sum((lam2 / total) * log((lam2 / total) + 1e-12)))
  m <- sum(A) / 2
  if (m > 0) {
    P <- (deg %*% t(deg)) / (2 * m)
    B <- A - P
    lamB <- eigen(0.5 * (B + t(B)), symmetric = TRUE)$values
    Qmax <- sum(pmax(lamB, 0)) / (2 * m + 1e-12)
  } else {
    Qmax <- NA_real_
  }
  C <- .mppi_adv_expm_sym(A)
  comm_avg <- (sum(C) - sum(diag(C))) / (p * (p - 1) + 1e-12)
  list(alg_conn = alg_conn,
       effrankA = effrankA,
       Qmax = Qmax,
       comm_avg = comm_avg)
}

#' Examine how baseline precision gates gain vs routing
#'
#' @param fits list of `mppi_fit` objects
#' @param k condition index
#' @param lambda ridge parameter for Theta0
#' @param controls optional data frame of control variables (rows align with fits)
#' @return list with combined data and lm objects for G and C
#' @export
mppi_gate_by_theta <- function(fits, k, lambda = 1e-3, controls = NULL) {
  stopifnot(is.list(fits), length(fits) >= 1)
  rows <- vector("list", length(fits))
  for (i in seq_along(fits)) {
    fit <- fits[[i]]
    X <- .mppi_adv_X(fit)
    S0 <- crossprod(X) / nrow(X)
    Theta0 <- solve(S0 + diag(lambda, ncol(S0)))
    mts <- mppi_theta_metrics(Theta0)
    spec <- mppi_spectral_summary(fit, k)
    prec <- mppi_precision_delta(fit, k, lambda = lambda)
    rows[[i]] <- c(G = spec$lambda1_share,
                   C = prec$offdiag_energy,
                   alg_conn = mts$alg_conn,
                   effrankA = mts$effrankA,
                   Qmax = mts$Qmax,
                   comm_avg = mts$comm_avg)
  }
  df <- as.data.frame(do.call(rbind, rows))
  if (!is.null(controls)) {
    stopifnot(nrow(controls) == nrow(df))
    df <- cbind(df, controls)
  }
  fm <- reformulate(setdiff(names(df), c("G", "C")), response = "G")
  model_G <- stats::lm(fm, data = df)
  fmC <- reformulate(setdiff(names(df), c("G", "C")), response = "C")
  model_C <- stats::lm(fmC, data = df)
  list(data = df, model_G = model_G, model_C = model_C)
}

# Internal utilities ---------------------------------------------------------

.mppi_adv_X <- function(fit) {
  if (!is.null(fit$U)) return(fit$U)
  if (!is.null(fit$Z)) return(fit$Z)
  if (!is.null(fit$R)) return(fit$R)
  stop("fit must carry U, Z, or R time-series.")
}

.mppi_adv_pk <- function(fit, k) {
  stopifnot(!is.null(fit$pk), k >= 1, k <= length(fit$pk))
  as.numeric(fit$pk[[k]])
}

.mppi_adv_Mk <- function(fit, k, X = NULL) {
  pk <- .mppi_adv_pk(fit, k)
  if (is.null(X)) X <- .mppi_adv_X(fit)
  denom <- sum(pk^2) + 1e-12
  crossprod(X, pk * X) / denom
}

.mppi_adv_lagged_cross <- function(X, pk, lag) {
  lag <- as.integer(lag)
  Tn <- nrow(X)
  if (lag == 0L) return(.mppi_adv_Mk_list(X, pk))
  if (abs(lag) >= Tn) return(matrix(NA_real_, ncol(X), ncol(X)))
  if (lag > 0) {
    idx <- seq_len(Tn - lag)
    X0 <- X[idx, , drop = FALSE]
    X1 <- X[idx + lag, , drop = FALSE]
    w <- pk[idx]
  } else {
    shift <- abs(lag)
    idx <- seq_len(Tn - shift)
    X0 <- X[idx + shift, , drop = FALSE]
    X1 <- X[idx, , drop = FALSE]
    w <- pk[idx]
  }
  denom <- sum(w^2) + 1e-12
  crossprod(X0, w * X1) / denom
}

.mppi_adv_Mk_list <- function(X, pk) {
  denom <- sum(pk^2) + 1e-12
  crossprod(X, pk * X) / denom
}

.mppi_adv_perm_shift <- function(pk, B) {
  Tn <- length(pk)
  vapply(seq_len(B), function(i) {
    s <- sample.int(Tn, 1)
    c(tail(pk, -s), head(pk, s))
  }, numeric(Tn))
}

.mppi_adv_window <- function(Tn, center, width, pre = TRUE) {
  if (pre) {
    seq.int(max(1L, center - width), max(1L, center - 1L))
  } else {
    seq.int(center, min(Tn, center + width - 1L))
  }
}

.mppi_adv_indicator <- function(Tn, idx) {
  v <- numeric(Tn)
  v[idx] <- 1
  v
}

.mppi_adv_nuisance <- function(nuisance_ts, idx) {
  if (is.null(nuisance_ts) || !length(nuisance_ts)) return(numeric(0))
  out <- vapply(nuisance_ts, function(ts) {
    vals <- ts[idx]
    c(mean(vals), stats::sd(vals), if (length(vals) > 1) stats::coef(stats::lm(vals ~ seq_along(vals)))[2] else 0)
  }, numeric(3))
  as.numeric(out)
}

.mppi_adv_cv_r <- function(data, runs, response, predictors) {
  runs <- as.integer(runs)
  uniq <- sort(unique(runs))
  preds <- numeric(0)
  obs <- numeric(0)
  for (rr in uniq) {
    idx_tr <- which(runs != rr)
    idx_te <- which(runs == rr)
    if (length(idx_tr) < 5 || length(idx_te) < 1) next
    fit <- stats::lm(reformulate(predictors, response), data = data[idx_tr, , drop = FALSE])
    mm <- model.matrix(reformulate(predictors), data = data[idx_te, , drop = FALSE])
    preds <- c(preds, as.numeric(mm %*% coef(fit)))
    obs <- c(obs, data[[response]][idx_te])
  }
  if (length(preds) > 1) stats::cor(preds, obs) else NA_real_
}

.mppi_adv_dct_matrix <- function(Tn, K) {
  k <- 0:(K - 1)
  n <- 0:(Tn - 1)
  M <- outer(n + 0.5, k, function(nn, kk) cos(pi * nn * kk / Tn))
  M[, 1] <- M[, 1] / sqrt(Tn)
  if (K > 1) M[, -1] <- M[, -1] * sqrt(2 / Tn)
  M
}

.mppi_adv_project_lowK <- function(X, K) {
  Tn <- nrow(X)
  K <- min(as.integer(K), Tn)
  if (K <= 0) return(0 * X)
  D <- .mppi_adv_dct_matrix(Tn, K)
  D %*% crossprod(D, X)
}

.mppi_adv_bandlimit <- function(X, TR, f_cut, which = c("low", "high")) {
  which <- match.arg(which)
  Tn <- nrow(X)
  K <- max(1L, min(Tn, floor(2 * Tn * TR * f_cut)))
  XL <- .mppi_adv_project_lowK(X, K)
  if (which == "low") XL else X - XL
}

.mppi_adv_ar_resid <- function(x, phi) {
  p <- length(phi)
  if (p == 0) return(x)
  Tn <- length(x)
  if (p >= Tn) p <- Tn - 1L
  res <- numeric(Tn)
  res[seq_len(p)] <- x[seq_len(p)]
  for (t in (p + 1):Tn) {
    res[t] <- x[t] - sum(phi * x[(t - 1):(t - p)])
  }
  res
}

.mppi_adv_expm_sym <- function(A) {
  eig <- eigen(0.5 * (A + t(A)), symmetric = TRUE)
  V <- eig$vectors
  L <- eig$values
  V %*% (diag(exp(pmax(pmin(L, 50), -50)), length(L))) %*% t(V)
}
