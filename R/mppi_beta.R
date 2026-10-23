# Beta-series versions -----------------------------------------------------

#' Coerce fmrilss output to a trial x V beta matrix
#'
#' For a data.frame, trial/ROI cells with no row, and NA betas, are NA in the
#' result; duplicate trial/ROI cells are an error.
#' @param obj matrix, data.frame with (trial, roi, beta), or list of numeric vectors
#' @export
mppi_coerce_beta <- function(obj) {
  if (is.matrix(obj)) return(obj)
  if (is.data.frame(obj)) {
    req <- c("trial","roi","beta")
    if (!all(req %in% names(obj))) stop("data.frame must have columns: trial, roi, beta")
    tr <- factor(obj$trial); roi <- factor(obj$roi)
    if (anyNA(tr) || anyNA(roi)) stop("trial and roi must not be NA.", call. = FALSE)
    cell <- (as.integer(roi) - 1L) * nlevels(tr) + as.integer(tr)
    if (anyDuplicated(cell)) stop("data.frame has duplicate trial/roi cells.", call. = FALSE)
    B <- matrix(NA_real_, nlevels(tr), nlevels(roi),
                dimnames = list(trial = levels(tr), roi = levels(roi)))
    B[cell] <- as.numeric(obj$beta)
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
#' @param runs optional vector tagging the run for each trial (used for centering)
#' @param zero_diag logical; zero the diagonal of the returned Δ matrices
#' @param scale compute covariance (`"cov"`) or correlation (`"corr"`) slopes
#' @param center_by strategy for mean-centering the psychological regressors
#' @param na_action how to handle rows with missing entries (`"omit_tr"` or stop)
#' @param backend backend used for cross-products (`"blas"`, `"accumulate"`, or `"chunked"`)
#' @param packed logical; return packed upper-triangle matrices
#' @param chunk_size optional chunk size when using the chunked backend
#' @return list with Delta, names, RB, pk (residualized regressors), and meta information per regressor
#' @export
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
#'
#' Blocks are the runs in `run_trial`, else consecutive chunks of `blksize`
#' trials, else one block of all trials. `"permute"` shuffles the raw regressor
#' within blocks and re-residualizes it; `"block_flip"` flips the sign of the
#' residualized regressor blockwise. `"auto"` flips continuous regressors when
#' at least two blocks carry them, and permutes otherwise. Configurations whose
#' null cannot vary (flips with fewer than two such blocks, or permutations of a
#' regressor constant within every block) are errors.
#' @inheritParams mppi_fit_beta
#' @param run_trial optional vector of run labels per trial (overrides `runs` for permutations)
#' @param blksize optional block length (trials) used when `run_trial` is absent
#' @param Bperm number of permutations to draw
#' @param seed optional RNG seed for reproducibility
#' @param method permutation strategy (`"auto"`, `"block_flip"`, or `"permute"`)
#' @export
mppi_beta_permute <- function(B, Pbeta, Cbeta = NULL, run_trial = NULL,
                              blksize = NULL, Bperm = 999L, zero_diag = TRUE,
                              seed = NULL, method = c("auto","block_flip","permute")) {
  method <- match.arg(method)
  if (length(Bperm) != 1L || !is.finite(Bperm) || Bperm < 1 || Bperm != as.integer(Bperm)) stop("Bperm must be a positive integer.", call. = FALSE)
  fit <- mppi_fit_beta(B, Pbeta, Cbeta, zero_diag = zero_diag, scale = "cov")
  if (fit$dropped > 0) stop("Permutation with incomplete trials is unsupported; remove them first.", call. = FALSE)
  Tn <- nrow(B)
  make_blocks <- function() {
    if (!is.null(run_trial)) {
      if (length(run_trial) != Tn || anyNA(run_trial)) stop("run_trial must label every trial.", call. = FALSE)
      split(seq_len(Tn), run_trial)
    } else if (!is.null(blksize)) {
      if (length(blksize) != 1L || !is.finite(blksize) || blksize < 1 || blksize != as.integer(blksize)) stop("blksize must be a positive integer.", call. = FALSE)
      split(seq_len(Tn), ceiling(seq_len(Tn)/blksize))
    } else {
      list(seq_len(Tn))
    }
  }
  idx_blocks <- make_blocks()
  res <- vector("list", length(fit$names))
  names(res) <- fit$names
  restore_rng <- .mppi_seed_rng(seed); on.exit(restore_rng(), add = TRUE)
  rng_used <- RNGkind()
  init_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
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
    n_live <- sum(vapply(idx_blocks, function(g) sum(pk0[g]^2) > 1e-12 * denom0, TRUE))
    flip <- method == "block_flip" || (method == "auto" && !meta$is_binary && n_live >= 2L)
    if (flip && n_live < 2L) {
      stop(sprintf("Block sign-flips of '%s' require at least two blocks carrying the regressor; supply run_trial or blksize, or use method = \"permute\".", nm), call. = FALSE)
    }
    if (!flip && all(vapply(idx_blocks, function(g) length(unique(meta$raw[g])) < 2L, TRUE))) {
      stop(sprintf("Within-block permutation cannot change '%s': it is constant within every block.", nm), call. = FALSE)
    }
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
      if (!flip) {
        pk_raw <- permute_block(meta$raw)
        pkb <- as.numeric(pk_raw - meta$Q %*% qr.coef(meta$qr, pk_raw))
      } else {
        pkb <- signflip_block(pk0)
      }
      denomb <- sum(pkb^2)
      if (!is.finite(denomb) || denomb < .Machine$double.eps) {
        stop("Randomized regressor has zero or invalid energy.", call. = FALSE)
      }
      Db <- crossprod(fit$RB, pkb * fit$RB) / denomb
      if (zero_diag) diag(Db) <- 0
      Qb[b] <- sum(Db^2, na.rm = TRUE)
    }
    p_global <- (1 + sum(Qb >= Q0)) / (Bperm + 1)
    res[[nm]] <- list(D = D0, Q = Q0, Qnull = Qb, p_global = p_global,
                      method = if (flip) "block_flip" else "permute")
  }
  final_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  attr(res, "rng") <- list(kind = rng_used, seed = init_seed, final_seed = final_seed,
                            method = method, seed_arg = seed)
  res
}
