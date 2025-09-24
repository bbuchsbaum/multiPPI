# Beta-series versions -----------------------------------------------------

#' Coerce fmrilss output to a trial x V beta matrix
#' @param obj matrix, data.frame with (trial, roi, beta), or list of numeric vectors
mppi_coerce_beta <- function(obj) {
  if (is.matrix(obj)) return(obj)
  if (is.data.frame(obj)) {
    req <- c("trial","roi","beta")
    if (!all(req %in% names(obj))) stop("data.frame must have columns: trial, roi, beta")
    tab <- xtabs(beta ~ trial + roi, data = obj)
    B <- as.matrix(tab)
    # Order trials numerically if rownames are numeric
    rn <- suppressWarnings(as.integer(rownames(B)))
    if (all(!is.na(rn))) B <- B[order(rn), , drop = FALSE]
    return(B)
  }
  if (is.list(obj) && all(vapply(obj, is.numeric, TRUE))) {
    lens <- vapply(obj, length, 1L)
    if (length(unique(lens)) != 1L) stop("All vectors in the list must have the same length (V).")
    return(do.call(rbind, obj))
  }
  stop("Can't coerce object to a trial x V beta matrix.")
}

#' Beta-series Matrix-PPI (β-mPPI)
#' @param B nTrials x V beta matrix
#' @param Pbeta nTrials x K psychological regressors (trial-level)
#' @param Cbeta optional nTrials x q trial-level confounds
#' @return list with Delta, names, RB, pk (residualized regressors), and meta information per regressor
mppi_fit_beta <- function(B, Pbeta, Cbeta = NULL, runs = NULL,
                          zero_diag = TRUE, scale = c("cov","corr"),
                          center_by = c("none","run"), na_action = c("omit_tr","error"),
                          backend = c("blas","accumulate","chunked"), packed = FALSE,
                          chunk_size = NULL) {
  stopifnot(is.matrix(B), is.matrix(Pbeta), nrow(B) == nrow(Pbeta))
  if (!is.null(runs)) stopifnot(length(runs) == nrow(B))
  scale <- match.arg(scale)
  center_by <- match.arg(center_by)
  na_action <- match.arg(na_action)
  backend <- match.arg(backend)
  if (center_by == "run" && is.null(runs)) {
    warning("center_by='run' requested but 'runs' is NULL; skipping run centering.", call. = FALSE)
    center_by <- "none"
  }
  Zb <- if (is.null(Cbeta)) cbind(1, Pbeta) else cbind(1, Cbeta, Pbeta)
  prep <- .mppi_handle_na(B, Zb, runs, na_action = na_action)
  B <- prep$Y; Zb <- prep$X; runs <- prep$runs
  if (prep$dropped > 0) {
    keep <- setdiff(seq_len(nrow(Pbeta)), prep$dropped_idx)
    Pbeta <- Pbeta[keep, , drop = FALSE]
    if (!is.null(Cbeta)) Cbeta <- Cbeta[keep, , drop = FALSE]
  }
  RB <- B - Zb %*% qr.coef(qr(Zb), B)
  if (scale == "corr") {
    s <- matrixStats::colSds(RB)
    s[s == 0] <- 1
    RB <- sweep(RB, 2, s, "/")
  }
  K <- ncol(Pbeta)
  out <- vector("list", K)
  pks <- vector("list", K)
  meta <- vector("list", K)
  denoms <- numeric(K)
  pNms <- colnames(Pbeta)
  if (is.null(pNms)) pNms <- paste0("psych", seq_len(K))
  denom_warn <- function(name) {
    warning(sprintf("Beta regressor '%s' has near-zero variance after residualization; returning NA matrix.",
                    name), call. = FALSE)
  }
  for (ii in seq_len(K)) {
    Qb <- if (is.null(Cbeta)) cbind(1, Pbeta[, -ii, drop = FALSE])
          else                cbind(1, Cbeta, Pbeta[, -ii, drop = FALSE])
    Qqr <- qr(Qb)
    raw_col <- Pbeta[, ii]
    pk <- as.numeric(raw_col - Qb %*% qr.coef(Qqr, raw_col))
    if (center_by == "run" && !is.null(runs)) pk <- .mppi_center_by_run(pk, runs)
    denom <- sum(pk^2)
    denoms[ii] <- denom
    Dk <- if (denom < .Machine$double.eps) {
      denom_warn(pNms[ii])
      naMat <- matrix(NA_real_, ncol(B), ncol(B))
      if (packed) .mppi_pack_upper(naMat) else naMat
    } else {
      Dfull <- .mppi_crossprod(RB, pk, backend = backend, chunk_size = chunk_size) / denom
      if (packed) .mppi_pack_upper(Dfull) else Dfull
    }
    if (!packed && zero_diag) diag(Dk) <- 0
    if (packed && zero_diag) {
      Dtemp <- .mppi_unpack_upper(Dk)
      diag(Dtemp) <- 0
      Dk <- .mppi_pack_upper(Dtemp)
    }
    out[[ii]] <- Dk
    pks[[ii]] <- pk
    meta[[ii]] <- list(Q = Qb, qr = Qqr, raw = as.numeric(raw_col),
                       is_binary = all(raw_col %in% c(0, 1)))
  }
  names(out) <- pNms
  names(pks) <- pNms
  names(meta) <- pNms
  names(denoms) <- pNms
  res <- list(Delta = out, names = pNms, RB = RB, pk = pks, meta = meta,
              dropped = prep$dropped, runs = runs, center_by = center_by,
              scale = scale, denom = denoms, backend = backend,
              n_used = nrow(RB), packed = packed, chunk_size = chunk_size)
  class(res) <- c("mppi_fit_beta", "list")
  res
}

#' Permutation for β-mPPI omnibus per regressor
mppi_beta_permute <- function(B, Pbeta, Cbeta = NULL, run_trial = NULL,
                              blksize = NULL, Bperm = 999L, zero_diag = TRUE,
                              seed = NULL, method = c("auto","block_flip","permute")) {
  method <- match.arg(method)
  fit <- mppi_fit_beta(B, Pbeta, Cbeta, zero_diag = zero_diag, scale = "cov")
  Tn <- nrow(B)
  make_blocks <- function() {
    if (!is.null(run_trial)) {
      split(seq_len(Tn), run_trial)
    } else if (!is.null(blksize)) {
      split(seq_len(Tn), ceiling(seq_len(Tn)/blksize))
    } else {
      list(seq_len(Tn))
    }
  }
  idx_blocks <- make_blocks()
  res <- vector("list", length(fit$names))
  names(res) <- fit$names
  rng_pre <- RNGkind()
  seed_pre <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  if (!is.null(seed)) set.seed(seed)
  rng_used <- RNGkind()
  init_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    do.call(RNGkind, as.list(rng_pre))
    if (is.null(seed_pre)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", seed_pre, envir = .GlobalEnv)
    }
  }, add = TRUE)
  for (nm in fit$names) {
    pk0 <- fit$pk[[nm]]
    denom0 <- sum(pk0^2)
    if (denom0 < .Machine$double.eps) {
      V <- ncol(fit$RB)
      res[[nm]] <- list(D = matrix(NA_real_, V, V), Q = NA_real_, Qnull = rep(NA_real_, Bperm),
                        p_global = NA_real_, note = "pk has near-zero variance; cannot evaluate permutations")
      next
    }
    D0 <- crossprod(fit$RB, pk0 * fit$RB) / denom0
    if (zero_diag) diag(D0) <- 0
    Q0 <- sum(D0^2, na.rm = TRUE)
    meta <- fit$meta[[nm]]
    Qb <- rep(NA_real_, Bperm)
    permute_block <- function(vec) {
      out <- vec
      for (g in idx_blocks) {
        out[g] <- sample(out[g])
      }
      out
    }
    signflip_block <- function(vec) {
      out <- vec
      sgn <- sample(c(-1,1), length(idx_blocks), replace = TRUE)
      jj <- 1L
      for (g in idx_blocks) {
        out[g] <- sgn[jj] * out[g]
        jj <- jj + 1L
      }
      out
    }
    for (b in seq_len(Bperm)) {
      if (meta$is_binary || method == "permute") {
        pk_raw <- permute_block(meta$raw)
        pkb <- as.numeric(pk_raw - meta$Q %*% qr.coef(meta$qr, pk_raw))
      } else {
        pkb <- signflip_block(pk0)
      }
      denomb <- sum(pkb^2)
      if (denomb < .Machine$double.eps) {
        next
      }
      Db <- crossprod(fit$RB, pkb * fit$RB) / denomb
      if (zero_diag) diag(Db) <- 0
      Qb[b] <- sum(Db^2, na.rm = TRUE)
    }
    valid <- !is.na(Qb)
    p_global <- if (!any(valid)) {
      NA_real_
    } else {
      (1 + sum(Qb[valid] >= Q0)) / (sum(valid) + 1)
    }
    res[[nm]] <- list(D = D0, Q = Q0, Qnull = Qb, p_global = p_global)
  }
  final_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  attr(res, "rng") <- list(kind = rng_used, seed = init_seed, final_seed = final_seed,
                            method = method, seed_arg = seed)
  res
}
