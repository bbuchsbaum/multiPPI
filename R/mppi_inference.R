# Inference: omnibus + permutations ---------------------------------------

#' Omnibus test for each psychological regressor via block sign-flips
#' @param fit result of mppi_fit()
#' @param blksize integer block length (TRs)
#' @param B number of permutations
#' @param wild "none","rademacher","mammen"
#' @return list per regressor with D, Q, Qnull, p_global
mppi_omnibus <- function(fit, blksize = 10L, B = 999L, wild = c("none","rademacher","mammen"),
                         method = c("block_flip","phase","freedman_lane"), seed = NULL) {
  wild <- match.arg(wild)
  method <- match.arg(method)
  if (method %in% c("phase","freedman_lane") && wild != "none") {
    warning("wild weights ignored when method='phase' or 'freedman_lane'.", call. = FALSE)
    wild <- "none"
  }
  out <- vector("list", length(fit$names)); names(out) <- fit$names
  Tn <- nrow(fit$R)
  idx <- split(seq_len(Tn), ceiling(seq_len(Tn)/blksize))
  runs_phase <- if (!is.null(fit$runs)) fit$runs else rep(1L, Tn)
  gen_weights <- function() {
    if (wild == "none") return(rep(1, length(idx)))
    if (wild == "rademacher") return(sample(c(-1,1), length(idx), replace = TRUE))
    # Mammen two-point weights
    p <- (sqrt(5)+1)/(2*sqrt(5))
    a <- (1 - sqrt(5))/2; b <- (1 + sqrt(5))/2
    sample(c(a,b), length(idx), replace = TRUE, prob = c(1-p, p))
  }
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
  if (method == "freedman_lane") {
    V <- ncol(fit$R)
    outer_centered <- .mppi_row_outer_packed(fit$R)
    col_means <- colMeans(outer_centered)
    outer_centered <- sweep(outer_centered, 2, col_means, "-")
  }
  for (i in seq_along(fit$names)) {
    pk0 <- fit$pk[[i]]
    denom0 <- sum(pk0^2)
    if (method == "freedman_lane") {
      pkc <- pk0 - mean(pk0)
      beta_obs_vec <- drop(crossprod(pkc, outer_centered) / denom0)
      D0_full <- .mppi_unpack_upper(list(values = beta_obs_vec, dim = ncol(fit$R)))
      diag(D0_full) <- 0
      D0 <- if (isTRUE(fit$packed)) .mppi_pack_upper(D0_full) else D0_full
      Q0 <- .mppi_frob2(D0)
      perm_fun <- function() {
        perm_idx <- unlist(sample(idx))
        beta_vec <- drop(crossprod(pkc, outer_centered[perm_idx, , drop = FALSE]) / denom0)
        full <- .mppi_unpack_upper(list(values = beta_vec, dim = ncol(fit$R)))
        diag(full) <- 0
        if (isTRUE(fit$packed)) .mppi_pack_upper(full) else full
      }
      Qb <- numeric(B)
      for (b in seq_len(B)) {
        Db <- perm_fun()
        Qb[b] <- .mppi_frob2(Db)
      }
      p_global <- (1 + sum(Qb >= Q0)) / (B + 1)
      out[[i]] <- list(D = D0, Q = Q0, Qnull = Qb, p_global = p_global)
      next
    }
    D0 <- .mppi_wcp(fit$R, pk0) / denom0; diag(D0) <- 0
    Q0 <- sum(D0^2, na.rm = TRUE)
    Qb <- numeric(B)
    for (b in seq_len(B)) {
      if (method == "phase") {
        pkb <- .mppi_phase_randomize(pk0, runs_phase)
      } else {
        w <- gen_weights()
        pkb <- numeric(Tn); jj <- 1L
        for (g in idx) { pkb[g] <- w[jj]*pk0[g]; jj <- jj + 1L }
      }
      Db <- .mppi_wcp(fit$R, pkb) / sum(pkb^2)
      Qb[b] <- sum(Db^2, na.rm = TRUE)
    }
    p_global <- (1 + sum(Qb >= Q0)) / (B + 1)
    out[[i]] <- list(D = D0, Q = Q0, Qnull = Qb, p_global = p_global)
  }
  final_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  attr(out, "rng") <- list(kind = rng_used, seed = init_seed, final_seed = final_seed,
                            method = method, seed_arg = seed)
  out
}

#' Edgewise permutations with optional studentization (experimental)
#' Returns p-values matrix for a given regressor index
mppi_permute <- function(fit, k = 1L, blksize = 10L, B = 999L,
                         wild = c("none","rademacher","mammen"), studentize = TRUE,
                         method = c("block_flip","phase","freedman_lane"), seed = NULL) {
  wild <- match.arg(wild)
  method <- match.arg(method)
  if (method %in% c("phase","freedman_lane") && wild != "none") {
    warning("wild weights ignored when method='phase' or 'freedman_lane'.", call. = FALSE)
    wild <- "none"
  }
  Tn <- nrow(fit$R)
  idx <- split(seq_len(Tn), ceiling(seq_len(Tn)/blksize))
  runs_phase <- if (!is.null(fit$runs)) fit$runs else rep(1L, Tn)
  gen_weights <- function() {
    if (wild == "none") return(rep(1, length(idx)))
    if (wild == "rademacher") return(sample(c(-1,1), length(idx), replace = TRUE))
    p <- (sqrt(5)+1)/(2*sqrt(5)); a <- (1 - sqrt(5))/2; b <- (1 + sqrt(5))/2
    sample(c(a,b), length(idx), replace = TRUE, prob = c(1-p, p))
  }
  pk0 <- fit$pk[[k]]
  denom0 <- sum(pk0^2)
  V <- ncol(fit$R)
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
  if (method == "freedman_lane") {
    if (studentize) warning("studentize ignored for method='freedman_lane'.", call. = FALSE)
    outer_centered <- .mppi_row_outer_packed(fit$R)
    col_means <- colMeans(outer_centered)
    outer_centered <- sweep(outer_centered, 2, col_means, "-")
    pkc <- pk0 - mean(pk0)
    obs_vec <- drop(crossprod(pkc, outer_centered) / denom0)
    D0_full <- .mppi_unpack_upper(list(values = obs_vec, dim = V))
    diag(D0_full) <- 0
    count <- matrix(0, V, V)
    abs_obs <- abs(D0_full)
    for (b in seq_len(B)) {
      perm_idx <- unlist(sample(idx))
      perm_vec <- drop(crossprod(pkc, outer_centered[perm_idx, , drop = FALSE]) / denom0)
      Dperm <- .mppi_unpack_upper(list(values = perm_vec, dim = V))
      diag(Dperm) <- 0
      count <- count + (abs(Dperm) >= abs_obs)
    }
    pmat <- (1 + count) / (B + 1)
    diag(pmat) <- NA_real_
    attr(pmat, "rng") <- list(kind = rng_used, seed = init_seed, final_seed = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL,
                               method = method, seed_arg = seed)
    return(pmat)
  }
  D0 <- .mppi_wcp(fit$R, pk0) / denom0; diag(D0) <- 0
  if (studentize) {
    S <- .mppi_wcp(fit$R, pk0^2) / denom0
    S[S <= 0] <- min(S[S > 0], na.rm = TRUE)
    Z0 <- D0 / sqrt(S)
  } else {
    Z0 <- D0
  }
  Zb <- array(0, dim = c(V,V,B))
  for (b in seq_len(B)) {
    if (method == "phase") {
      pkb <- .mppi_phase_randomize(pk0, runs_phase)
    } else {
      w <- gen_weights()
      pkb <- numeric(Tn); jj <- 1L
      for (g in idx) { pkb[g] <- w[jj]*pk0[g]; jj <- jj + 1L }
    }
    Db <- .mppi_wcp(fit$R, pkb) / sum(pkb^2); diag(Db) <- 0
    if (studentize) {
      Sb <- .mppi_wcp(fit$R, pkb^2) / sum(pkb^2)
      Sb[Sb <= 0] <- min(Sb[Sb > 0], na.rm = TRUE)
      Zb[,,b] <- Db / sqrt(Sb)
    } else {
      Zb[,,b] <- Db
    }
  }
  pmat <- matrix(NA_real_, V, V)
  for (i in seq_len(V)) for (j in seq_len(V)) {
    if (i == j) { pmat[i,j] <- NA_real_; next }
    z0 <- Z0[i,j]
    zb <- Zb[i,j,]
    pmat[i,j] <- (1 + sum(abs(zb) >= abs(z0))) / (B + 1)
  }
  final_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  attr(pmat, "rng") <- list(kind = rng_used, seed = init_seed, final_seed = final_seed,
                             method = method, seed_arg = seed)
  pmat
}
