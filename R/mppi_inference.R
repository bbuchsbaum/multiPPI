# Experimental covariance-modulation randomization ------------------------

.mppi_randomization_setup <- function(fit, blksize, B, method, wild) {
  if (method == "freedman_lane") stop("Freedman-Lane is unavailable: a qualified reduced product-nuisance model is required.", call. = FALSE)
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B)) stop("B must be a positive integer.", call. = FALSE)
  if (length(blksize) != 1L || !is.finite(blksize) || blksize < 1 || blksize != as.integer(blksize)) stop("blksize must be a positive integer.", call. = FALSE)
  if (method == "block_flip" && wild == "none") stop("Block randomization requires nonconstant wild weights.", call. = FALSE)
  n <- nrow(fit$R)
  runs <- fit$runs %||% rep(1L, n)
  if (length(runs) != n || anyNA(runs)) stop("Invalid fit run labels.", call. = FALSE)
  index <- fit$observation_index
  if (is.null(index)) {
    if ((fit$dropped %||% 0L) > 0L) stop("Censored fit lacks original observation indices for randomization.", call. = FALSE)
    index <- seq_len(n)
  }
  if (length(index) != n || anyNA(index) || any(diff(index) <= 0)) stop("Invalid original observation indices.", call. = FALSE)
  segment <- cumsum(c(TRUE, diff(index) != 1 | runs[-1L] != runs[-n]))
  segments <- split(seq_len(n), segment)
  blocks <- unlist(lapply(segments, function(i) split(i, ceiling(seq_along(i) / blksize))), recursive = FALSE)
  if (method == "block_flip" && length(blocks) < 2L && wild == "rademacher") stop("Block sign-flips require at least two blocks for a nondegenerate statistic.", call. = FALSE)
  if (method == "phase" && all(lengths(segments) < 3L)) stop("Phase randomization has no eligible segment.", call. = FALSE)
  list(blocks = blocks, segments = segments)
}

# With a seed, draw from a reproducible stream and return a function restoring
# the caller's RNG state. Without one, the caller's stream simply advances.
.mppi_seed_rng <- function(seed) {
  if (is.null(seed)) return(function() invisible(NULL))
  kind <- RNGkind()
  saved <- if (exists(".Random.seed", .GlobalEnv, inherits=FALSE)) get(".Random.seed", .GlobalEnv) else NULL
  set.seed(seed)
  function() {
    do.call(RNGkind, as.list(kind))
    if (!is.null(saved)) assign(".Random.seed", saved, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits=FALSE)) rm(".Random.seed", envir=.GlobalEnv)
  }
}

.mppi_random_weights <- function(n, wild) {
  if (wild == "rademacher") return(sample(c(-1, 1), n, replace = TRUE))
  a <- (1 - sqrt(5))/2; b <- (1 + sqrt(5))/2
  sample(c(a, b), n, replace = TRUE, prob = c((sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5))))
}

.mppi_random_psych <- function(pk, setup, method, wild) {
  if (method == "block_flip") {
    w <- .mppi_random_weights(length(setup$blocks), wild)
    for (j in seq_along(setup$blocks)) pk[setup$blocks[[j]]] <- pk[setup$blocks[[j]]] * w[j]
  } else {
    for (i in setup$segments) {
      n <- length(i)
      if (n < 3L) next
      f <- fft(pk[i]); max_k <- floor((n - 1L)/2L)
      for (k in seq_len(max_k)) {
        angle <- runif(1L, 0, 2*pi); f[k+1L] <- f[k+1L] * exp(1i*angle)
        f[n-k+1L] <- Conj(f[k+1L])
      }
      pk[i] <- Re(fft(f, inverse=TRUE)/n)
    }
  }
  pk
}

.mppi_random_metadata <- function(method, wild, blksize, setup) {
  list(status = "experimental", method = method, weights = if (method == "block_flip") wild else NULL,
       block_size = blksize, segment_lengths = lengths(setup$segments),
       assumptions = if (method == "block_flip") "Blockwise sign or wild-weight symmetry; calibration is unqualified." else "Within-segment stationarity and Fourier phase exchangeability; calibration is unqualified.")
}

# Energy sum M_ij^2 of M = R' diag(z) R / z'z, off the diagonal if zero_diag.
# With G = RR' it equals (z'(G o G)z - [||(R o R)'z||^2]) / (z'z)^2, so
# precomputing G o G costs O(T^2) per draw instead of O(T V^2) when T <= V^2.
.mppi_slope_energy <- function(R, zero_diag = TRUE, kernel = nrow(R) <= ncol(R)^2 && nrow(R) <= 4000L) {
  if (!kernel) return(function(z) { M <- .mppi_wcp(R, z)/sum(z^2); if (zero_diag) diag(M) <- 0; sum(M^2) })
  K <- tcrossprod(R)^2; S <- R^2
  function(z) (sum(z * (K %*% z)) - if (zero_diag) sum(crossprod(S, z)^2) else 0) / sum(z^2)^2
}

#' Experimental omnibus randomization of covariance-modulation slopes
#' @param fit A covariance-modulation fit.
#' @param blksize Positive block length.
#' @param B Number of draws.
#' @param wild Block weights; the default is Rademacher.
#' @param method Randomization method; Freedman-Lane is unavailable.
#' @param seed Optional random seed.
#' @return Per-regressor statistics and experimental randomization metadata.
mppi_omnibus <- function(fit, blksize = 10L, B = 999L, wild = c("rademacher","mammen","none"),
                         method = c("block_flip","phase","freedman_lane"), seed = NULL) {
  wild <- match.arg(wild); method <- match.arg(method)
  setup <- .mppi_randomization_setup(fit, blksize, B, method, wild)
  restore_rng <- .mppi_seed_rng(seed); on.exit(restore_rng(), add=TRUE)
  initial <- if (exists(".Random.seed", .GlobalEnv, inherits=FALSE)) get(".Random.seed", .GlobalEnv) else NULL
  energy <- .mppi_slope_energy(fit$R)
  out <- setNames(vector("list", length(fit$names)), fit$names)
  for (i in seq_along(out)) {
    pk <- fit$pk[[i]]; denom <- sum(pk^2)
    if (!is.finite(denom) || denom <= .Machine$double.eps) stop("Psychological regressor has zero or invalid energy.", call.=FALSE)
    D <- .mppi_wcp(fit$R, pk)/denom; diag(D) <- 0
    Q <- sum(D^2); null <- numeric(B)
    for (b in seq_len(B)) {
      z <- .mppi_random_psych(pk, setup, method, wild); dz <- sum(z^2)
      if (!is.finite(dz) || dz <= .Machine$double.eps) stop("Randomized regressor has zero or invalid energy.", call.=FALSE)
      null[b] <- energy(z)
    }
    out[[i]] <- list(D=D, Q=Q, Qnull=null, p_global=(1+sum(null>=Q))/(B+1))
  }
  attr(out,"rng") <- list(kind=RNGkind(), seed=initial, final_seed=get(".Random.seed", .GlobalEnv), method=method, seed_arg=seed)
  attr(out,"inference") <- .mppi_random_metadata(method, wild, blksize, setup)
  out
}

#' Experimental edgewise covariance-modulation randomization
#' @inheritParams mppi_omnibus
#' @param k Regressor index or name.
#' @param studentize Must be FALSE until a valid variance estimator is available.
#' @return Matrix of experimental edgewise randomization values.
mppi_permute <- function(fit, k = 1L, blksize = 10L, B = 999L,
                         wild = c("rademacher","mammen","none"), studentize = FALSE,
                         method = c("block_flip","phase","freedman_lane"), seed = NULL) {
  wild <- match.arg(wild); method <- match.arg(method)
  if (!identical(studentize, FALSE)) stop("Studentization is unavailable pending a valid variance estimator.", call.=FALSE)
  setup <- .mppi_randomization_setup(fit, blksize, B, method, wild)
  pk <- fit$pk[[k]]; denom <- sum(pk^2)
  if (!is.finite(denom) || denom <= .Machine$double.eps) stop("Psychological regressor has zero or invalid energy.", call.=FALSE)
  restore_rng <- .mppi_seed_rng(seed); on.exit(restore_rng(), add=TRUE)
  initial <- if (exists(".Random.seed", .GlobalEnv, inherits=FALSE)) get(".Random.seed", .GlobalEnv) else NULL
  D <- .mppi_wcp(fit$R, pk)/denom; diag(D) <- 0
  count <- matrix(0, ncol(fit$R), ncol(fit$R))
  for (b in seq_len(B)) {
    z <- .mppi_random_psych(pk, setup, method, wild); dz <- sum(z^2)
    if (!is.finite(dz) || dz <= .Machine$double.eps) stop("Randomized regressor has zero or invalid energy.", call.=FALSE)
    M <- .mppi_wcp(fit$R, z)/dz; diag(M) <- 0
    count <- count + (abs(M) >= abs(D))
  }
  p <- (1+count)/(B+1); diag(p) <- NA_real_
  attr(p,"rng") <- list(kind=RNGkind(), seed=initial, final_seed=get(".Random.seed", .GlobalEnv), method=method, seed_arg=seed)
  attr(p,"inference") <- .mppi_random_metadata(method, wild, blksize, setup)
  p
}
