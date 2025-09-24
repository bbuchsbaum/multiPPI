# High-level design interface ---------------------------------------------

.mppi_null_coalesce <- function(a, b) if (!is.null(a)) a else b
.mppi_stopif <- function(cond, msg) if (!cond) stop(msg, call. = FALSE)

.mppi_as_matrix <- function(obj) {
  if (is.null(obj)) return(NULL)
  if (is.matrix(obj)) return(obj)
  if (inherits(obj, "tbl_df")) obj <- as.data.frame(obj)
  as.matrix(obj)
}

# Generics -----------------------------------------------------------------

mppi_tr <- function(x, ...) UseMethod("mppi_tr")
mppi_runs <- function(x, ...) UseMethod("mppi_runs")
mppi_psych <- function(x, ...) UseMethod("mppi_psych")
mppi_task <- function(x, ...) UseMethod("mppi_task")
mppi_nuisance <- function(x, ...) UseMethod("mppi_nuisance")
mppi_hrf <- function(x, ...) UseMethod("mppi_hrf")
mppi_pk <- function(x, ...) UseMethod("mppi_pk")
mppi_design <- function(x, ...) UseMethod("mppi_design")

mppi_tr.default <- function(x, ...) stop(sprintf("No mppi_tr() method for class %s", paste(class(x), collapse = ", ")), call. = FALSE)
mppi_runs.default <- function(x, ...) stop(sprintf("No mppi_runs() method for class %s", paste(class(x), collapse = ", ")), call. = FALSE)
mppi_psych.default <- function(x, ...) stop(sprintf("No mppi_psych() method for class %s", paste(class(x), collapse = ", ")), call. = FALSE)
mppi_task.default <- function(x, ...) stop(sprintf("No mppi_task() method for class %s", paste(class(x), collapse = ", ")), call. = FALSE)
mppi_nuisance.default <- function(x, ...) stop(sprintf("No mppi_nuisance() method for class %s", paste(class(x), collapse = ", ")), call. = FALSE)
mppi_hrf.default <- function(x, ...) NULL
mppi_pk.default <- function(x, ...) stop(sprintf("No mppi_pk() method for class %s", paste(class(x), collapse = ", ")), call. = FALSE)
mppi_design.default <- function(x, ...) stop(sprintf("No mppi_design() method for class %s", paste(class(x), collapse = ", ")), call. = FALSE)

mppi_tr.mppi_design <- function(x, ...) mppi_tr(x$event_model)

mppi_runs.mppi_design <- function(x, T = NULL, ...) {
  if (!is.null(x$runs)) return(as.integer(x$runs))
  mppi_runs(x$event_model, T = T)
}

mppi_psych.mppi_design <- function(x, T, TR = mppi_tr(x), runs = mppi_runs(x, T), ...) {
  mat <- x$X[, x$psych_idx, drop = FALSE]
  if (!is.null(T)) .mppi_stopif(nrow(mat) == T, "Design has incompatible number of rows.")
  mat
}

mppi_task.mppi_design <- function(x, ...) x$X[, x$psych_idx, drop = FALSE]

mppi_nuisance.mppi_design <- function(x, T, ...) {
  other_idx <- setdiff(seq_len(ncol(x$X)), x$psych_idx)
  if (!length(other_idx)) return(matrix(0, nrow(x$X), 0))
  mat <- x$X[, other_idx, drop = FALSE]
  if (!is.null(T)) .mppi_stopif(nrow(mat) == T, "Design has incompatible number of rows.")
  mat
}

mppi_hrf.mppi_design <- function(x, ...) mppi_hrf(x$event_model)

mppi_design.mppi_design <- function(x, T = NULL, ...) {
  TR <- mppi_tr(x)
  runs_vec <- mppi_runs(x, T)
  res <- list(TR = TR,
              runs = runs_vec,
              X_task = mppi_task(x),
              hrf = mppi_hrf(x))
  if (!is.null(T)) {
    res$S <- mppi_psych(x, T = T, TR = TR, runs = runs_vec)
    res$nuisance <- mppi_nuisance(x, T = T)
  }
  res
}

# event_model / baseline_model support ------------------------------------

mppi_tr.event_model <- function(x, ...) {
  sf <- .mppi_null_coalesce(x$sampling_frame, attr(x, "sampling_frame"))
  if (!is.null(sf) && !is.null(sf$TR)) return(as.numeric(sf$TR[1]))
  NA_real_
}

mppi_runs.event_model <- function(x, T = NULL, ...) {
  if (!requireNamespace("fmridesign", quietly = TRUE)) return(NULL)
  sf <- .mppi_null_coalesce(x$sampling_frame, attr(x, "sampling_frame"))
  if (!is.null(sf$blocklens)) {
    return(as.integer(rep(seq_along(sf$blocklens), sf$blocklens)))
  }
  if (!is.null(T)) rep_len(1L, T) else NULL
}

mppi_psych.event_model <- function(x, T, TR = mppi_tr(x), runs = mppi_runs(x, T), ...) {
  if (!requireNamespace("fmridesign", quietly = TRUE)) stop("fmridesign not available.")
  S <- .mppi_as_matrix(fmridesign::design_matrix(x))
  if (!is.null(T)) .mppi_stopif(nrow(S) == T, "event_model design has incompatible number of rows.")
  S
}

mppi_hrf.event_model <- function(x, ...) {
  if (!requireNamespace("fmridesign", quietly = TRUE)) return(NULL)
  spec <- tryCatch(fmridesign::hrf(x), error = function(e) NULL)
  if (!is.null(spec$hrf)) spec$hrf else NULL
}

mppi_tr.baseline_model <- function(x, ...) {
  sf <- .mppi_null_coalesce(x$sampling_frame, attr(x, "sampling_frame"))
  if (!is.null(sf) && !is.null(sf$TR)) return(as.numeric(sf$TR[1]))
  NA_real_
}

mppi_runs.baseline_model <- function(x, T = NULL, ...) {
  if (!requireNamespace("fmridesign", quietly = TRUE)) return(NULL)
  sf <- .mppi_null_coalesce(x$sampling_frame, attr(x, "sampling_frame"))
  if (!is.null(sf$blocklens)) return(as.integer(rep(seq_along(sf$blocklens), sf$blocklens)))
  if (!is.null(T)) rep_len(1L, T) else NULL
}

mppi_nuisance.baseline_model <- function(x, T, ...) {
  if (!requireNamespace("fmridesign", quietly = TRUE)) stop("fmridesign not available.")
  N <- .mppi_as_matrix(fmridesign::design_matrix(x))
  if (!is.null(T)) .mppi_stopif(nrow(N) == T, "baseline_model design has incompatible number of rows.")
  N
}

# fmridesign composite -----------------------------------------------------

mppi_tr.fmridesign <- function(x, ...) {
  if (!is.null(x$event_model)) return(mppi_tr(x$event_model))
  out <- .mppi_null_coalesce(attr(x, "TR"), .mppi_null_coalesce(x$TR, x$tr))
  .mppi_stopif(length(out) > 0, "fmridesign object is missing TR.")
  as.numeric(out)
}

mppi_runs.fmridesign <- function(x, T = NULL, ...) {
  if (!is.null(x$event_model)) {
    r <- mppi_runs(x$event_model, T)
    if (!is.null(r)) return(r)
  }
  r <- .mppi_null_coalesce(attr(x, "runs"), .mppi_null_coalesce(x$runs, x$run))
  if (!length(r)) {
    if (!is.null(T)) return(rep_len(1L, T))
    return(NULL)
  }
  as.integer(r)
}

mppi_psych.fmridesign <- function(x, T, TR = mppi_tr(x), runs = mppi_runs(x, T), ...) {
  if (!is.null(x$event_model)) return(mppi_psych(x$event_model, T = T, TR = TR, runs = runs))
  S <- .mppi_null_coalesce(x$S, attr(x, "S"))
  if (!is.null(S)) {
    .mppi_stopif(nrow(S) == T, "fmridesign$S has incompatible number of rows.")
    return(as.matrix(S))
  }
  stop("Unable to recover psychological design from fmridesign object.")
}

mppi_task.fmridesign <- function(x, ...) {
  if (!is.null(x$event_model)) return(mppi_psych(x$event_model, T = NULL))
  X <- .mppi_null_coalesce(x$X_task, .mppi_null_coalesce(attr(x, "X_task"), .mppi_null_coalesce(x$X, attr(x, "X"))))
  .mppi_stopif(is.matrix(X), "fmridesign object is missing task design (X_task/X).")
  as.matrix(X)
}

mppi_nuisance.fmridesign <- function(x, T, ...) {
  if (!is.null(x$baseline_model)) return(mppi_nuisance(x$baseline_model, T = T))
  N <- .mppi_null_coalesce(x$nuisance, .mppi_null_coalesce(attr(x, "nuisance"), .mppi_null_coalesce(x$confounds, attr(x, "confounds"))))
  if (is.null(N)) return(matrix(0, T, 0))
  .mppi_stopif(nrow(N) == T, "fmridesign nuisance/confounds have incompatible number of rows.")
  as.matrix(N)
}

mppi_hrf.fmridesign <- function(x, ...) {
  if (!is.null(x$event_model)) return(mppi_hrf(x$event_model))
  .mppi_null_coalesce(x$hrf, attr(x, "hrf"))
}

mppi_pk.fmridesign <- function(x, T, TR = mppi_tr(x), runs = mppi_runs(x, T),
                          center = TRUE, span_space = TRUE,
                          include_nuisance = TRUE, include_intercept = TRUE, ...) {
  if (is.null(T)) stop("pk(): 'T' (number of TRs) must be supplied.", call. = FALSE)
  S <- mppi_psych(x, T = T, TR = TR, runs = runs)
  N <- if (include_nuisance) mppi_nuisance(x, T = T) else matrix(0, T, 0)
  K <- ncol(S)
  out <- vector("list", K)
  names(out) <- colnames(S)
  for (kk in seq_len(K)) {
    s_k <- S[, kk]
    if (center) s_k <- s_k - mean(s_k)
    others <- if (span_space && K > 1) S[, setdiff(seq_len(K), kk), drop = FALSE] else NULL
    intercept_col <- if (include_intercept) matrix(1, T, 1) else NULL
    Qs <- cbind(others, N, intercept_col)
    if (ncol(Qs)) {
      qr_obj <- qr(Qs)
      s_k <- s_k - qr.fitted(qr_obj, s_k)
    }
    out[[kk]] <- as.numeric(s_k)
  }
  out
}

mppi_design.fmridesign <- function(x, T = NULL, ...) {
  TR <- mppi_tr(x)
  runs_vec <- mppi_runs(x, T)
  res <- list(TR = TR,
              runs = runs_vec,
              X_task = mppi_task(x),
              hrf = mppi_hrf(x))
  if (!is.null(T)) {
    res$S <- mppi_psych(x, T = T, TR = TR, runs = runs_vec)
    res$nuisance <- mppi_nuisance(x, T = T)
  }
  res
}

# mppi_design constructor --------------------------------------------------

as_mppi_design <- function(event, baseline = NULL, confounds = NULL,
                           include_intercept = TRUE, runs = NULL, basis = NULL) {
  if (!requireNamespace("fmridesign", quietly = TRUE)) {
    stop("fmridesign not available. Install via remotes::install_github('bbuchsbaum/fmridesign').")
  }
  event_X <- .mppi_as_matrix(fmridesign::design_matrix(event))
  event_colnames <- colnames(event_X)
  base_X <- if (!is.null(baseline)) .mppi_as_matrix(fmridesign::design_matrix(baseline)) else NULL
  conf_X <- .mppi_as_matrix(confounds)
  if (!is.null(base_X) && nrow(base_X) != nrow(event_X)) {
    stop("baseline design must have the same number of rows as the event design.")
  }
  if (!is.null(conf_X) && nrow(conf_X) != nrow(event_X)) {
    stop("confounds matrix must have the same number of rows as the event design.")
  }
  pieces <- Filter(Negate(is.null), list(base_X, conf_X))
  base_conf <- if (length(pieces)) do.call(cbind, pieces) else NULL
  X <- if (is.null(base_conf)) event_X else cbind(base_conf, event_X)
  psych_idx <- seq(ncol(X) - ncol(event_X) + 1L, ncol(X))
  if (include_intercept) {
    has_int <- any(apply(X, 2, function(col) all(abs(col - 1) < 1e-8)))
    if (!has_int) {
      X <- cbind(`(Intercept)` = 1, X)
      psych_idx <- psych_idx + 1L
    }
  }
  if (is.null(runs)) runs <- mppi_runs(event, T = nrow(event_X))
  structure(list(
    X = as.matrix(X),
    psych_idx = psych_idx,
    event_model = event,
    baseline_model = baseline,
    runs = runs,
    basis = basis
  ), class = "mppi_design")
}

mppi_model <- function(event_model, baseline_model = NULL,
                       nuisance = NULL, runs = NULL,
                       include_intercept = TRUE) {
  as_mppi_design(event = event_model, baseline = baseline_model,
                 confounds = nuisance, include_intercept = include_intercept,
                 runs = runs)
}

# Print / summary ----------------------------------------------------------

print.mppi_design <- function(x, ...) {
  cat(sprintf("mppi_design: %d time points x %d regressors (%d psychological)\n",
              nrow(x$X), ncol(x$X), length(x$psych_idx)))
  if (!is.null(colnames(x$X))) {
    preview <- head(colnames(x$X)[x$psych_idx], 5L)
    cat("Psychological columns: ", paste(preview, collapse = ", "))
    if (length(x$psych_idx) > length(preview)) cat(", ...")
    cat("\n")
  }
  invisible(x)
}
