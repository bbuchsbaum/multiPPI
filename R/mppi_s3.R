# S3 interface ---------------------------------------------------------------

#' MultiPPI generic
#' @param x object containing either time series or a prepared design
#' @param ... additional arguments passed to class-specific methods
mppi <- function(x, ...) {
  UseMethod("mppi")
}

#' @export
mppi.matrix <- function(x, X, psych_idx, runs = NULL, prewhiten = FALSE,
                        basis = NULL, domain = c("bold","neural","innovations"),
                        deconv = NULL, lags = 0L, lag_blocklens = NULL,
                        engine = c("closed", "instant"), instant = list(), ...) {
  domain <- match.arg(domain)
  engine <- match.arg(engine)
  if (identical(engine, "instant")) {
    res <- mppi_instant_fit(Y = x, X = X, psych_idx = psych_idx, runs = runs,
                            prewhiten = prewhiten, basis = basis, domain = domain,
                            deconv = deconv, instant = instant, ...)
  } else if (prewhiten) {
    if (is.null(runs)) stop("'runs' must be supplied when prewhiten = TRUE.")
    res <- mppi_fit_whitened(Y = x, X = X, runs = runs, psych_idx = psych_idx,
                             basis = basis, domain = domain, deconv = deconv,
                             lags = lags, lag_blocklens = lag_blocklens, ...)
  } else {
    res <- mppi_fit(Y = x, X = X, psych_idx = psych_idx, runs = runs,
                    basis = basis, domain = domain, deconv = deconv,
                    lags = lags, lag_blocklens = lag_blocklens, ...)
  }
  res
}

#' @export
mppi.mppi_design <- function(x, Y, runs = NULL, prewhiten = FALSE,
                             basis = NULL, domain = c("bold","neural","innovations"),
                             deconv = NULL, lags = 0L, lag_blocklens = NULL,
                             engine = c("closed", "instant"), instant = list(), ...) {
  if (!is.matrix(Y)) stop("'Y' must be a matrix (time x V).")
  domain <- match.arg(domain)
  engine <- match.arg(engine)
  basis <- basis %||% x$basis
  if (is.null(basis)) {
    basis_obj <- NULL
    basis_source <- NULL
  } else {
    basis_obj <- .mppi_resolve_basis(Y, basis = basis, design = x)
    if (is.null(basis_obj)) {
      basis_source <- if (is.character(basis)) tolower(basis[[1L]]) else NULL
    } else {
      basis_source <- basis_obj$source %||% "provided"
    }
  }
  if (domain != "bold" && !is.null(x$event_model)) {
    tr <- mppi_tr(x$event_model)
    if (length(tr) == 1L && is.finite(tr)) {
      deconv <- deconv %||% list()
      supplied <- deconv$tr %||% deconv$TR
      if (!is.null(supplied) && (length(supplied) != 1L || abs(supplied - tr) > 1e-6)) stop("TR mismatch between design and deconvolution.", call. = FALSE)
      deconv$tr <- tr
    }
  }
  runs <- runs %||% x$runs
  if (is.null(runs) && !is.null(x$event_model)) {
    runs <- mppi_runs(x$event_model, T = nrow(Y))
  }
  if (identical(engine, "instant")) {
    res <- mppi_instant_fit(Y = Y, X = x$X, psych_idx = x$psych_idx, runs = runs,
                            prewhiten = prewhiten, basis = basis_obj, domain = domain,
                            deconv = deconv, instant = instant, ...)
  } else if (prewhiten) {
    if (is.null(runs)) runs <- x$runs
    if (is.null(runs)) stop("'runs' must be provided for prewhitening.")
    res <- mppi_fit_whitened(Y = Y, X = x$X, runs = runs, psych_idx = x$psych_idx,
                             basis = basis_obj, domain = domain, deconv = deconv,
                             lags = lags, lag_blocklens = lag_blocklens, ...)
  } else {
    res <- mppi_fit(Y = Y, X = x$X, psych_idx = x$psych_idx, runs = runs,
                    basis = basis_obj, domain = domain, deconv = deconv,
                    lags = lags, lag_blocklens = lag_blocklens, ...)
  }
  if (!is.null(basis_source)) {
    res$basis_source <- basis_source
    if (!is.null(res$basis)) res$basis$source <- basis_source
  }
  res
}

#' @export
mppi.event_model <- function(x, Y = NULL, baseline_model = NULL, nuisance = NULL, dataset = NULL,
                             runs = NULL, prewhiten = FALSE, include_intercept = TRUE, basis = NULL,
                             domain = c("bold","neural","innovations"), deconv = NULL,
                             lags = 0L, lag_blocklens = NULL,
                             engine = c("closed", "instant"), instant = list(), ...) {
  domain <- match.arg(domain)
  engine <- match.arg(engine)
  design <- as_mppi_design(event = x, baseline = baseline_model, confounds = nuisance,
                           include_intercept = include_intercept, runs = runs, basis = basis)
  basis <- basis %||% design$basis
  runs <- runs %||% design$runs
  if (!is.null(dataset) && is.null(Y)) {
    return(mppi(dataset, design, runs = runs, prewhiten = prewhiten,
                basis = basis, domain = domain, deconv = deconv,
                lags = lags, lag_blocklens = lag_blocklens,
                engine = engine, instant = instant, ...))
  }
  if (is.null(Y)) stop("Provide either 'Y' or 'dataset' when calling mppi() with an event_model.")
  mppi(design, Y = Y, runs = runs, prewhiten = prewhiten,
       basis = basis, domain = domain, deconv = deconv,
       lags = lags, lag_blocklens = lag_blocklens,
       engine = engine, instant = instant, ...)
}

#' @export
mppi.fmridesign <- function(x, Y, basis = NULL, domain = c("bold","neural","innovations"),
                            prewhiten = FALSE, runs = NULL, deconv = NULL,
                            engine = c("closed", "instant"), instant = list(), ...) {
  if (!is.matrix(Y)) stop("'Y' must be a matrix (time x V).")
  domain <- match.arg(domain)
  engine <- match.arg(engine)
  T <- nrow(Y)
  info <- mppi_design(x, T = T)
  basis_input <- basis %||% x$basis
  if (!is.null(basis_input)) {
    basis_obj <- .mppi_resolve_basis(Y, basis = basis_input, design = x)
  } else {
    basis_obj <- NULL
  }
  runs_vec <- runs %||% info$runs
  if (is.null(runs_vec)) runs_vec <- rep_len(1L, T)
  if (length(runs_vec) != T) stop("Length of 'runs' must match number of time points.")

  nuisance_mat <- info$nuisance
  if (is.null(nuisance_mat)) nuisance_mat <- matrix(0, T, 0)
  X_task <- info$X_task
  if (nrow(X_task) != T) stop("Task design has incompatible number of rows.")
  X <- if (ncol(nuisance_mat)) cbind(nuisance_mat, X_task) else X_task
  psych_idx <- if (ncol(nuisance_mat)) seq.int(ncol(nuisance_mat) + 1L, ncol(X)) else seq_len(ncol(X))

  design_obj <- structure(list(X = X,
                               psych_idx = psych_idx,
                               runs = runs_vec,
                               basis = NULL),
                          class = "mppi_design")

  deconv_args <- deconv %||% list()
  if (!is.null(info$hrf) && is.null(deconv_args$hrf)) deconv_args$hrf <- info$hrf
  if (!is.null(info$TR)) {
    if (is.null(deconv_args$TR)) deconv_args$TR <- info$TR
    if (is.null(deconv_args$tr)) deconv_args$tr <- info$TR
  }

  fit <- mppi(design_obj, Y = Y, runs = runs_vec, prewhiten = prewhiten,
              basis = basis_obj, domain = domain, deconv = deconv_args,
              engine = engine, instant = instant, ...)
  fit$design <- info
  fit$design_source <- x
  if (!is.null(basis_obj)) {
    source_tag <- basis_obj$source %||% "provided"
    if (!is.null(fit$basis)) fit$basis$source <- source_tag
    fit$basis_source <- source_tag
  } else if (!is.null(basis_input) && is.character(basis_input)) {
    fit$basis_source <- tolower(basis_input[[1L]])
  }
  fit
}

#' @export
mppi.fmri_dataset <- function(x, fd, basis = NULL, domain = c("bold","neural","innovations"),
                              prewhiten = FALSE, runs = NULL, deconv = NULL,
                              r = NULL, engine = c("closed", "instant"), instant = list(), ...) {
  domain <- match.arg(domain)
  engine <- match.arg(engine)
  Y <- NULL
  TR_ds <- NULL
  runs_ds <- runs
  if (is.null(Y)) {
    Y <- x$Y %||% x$timeseries %||% x$X
  }
  if (!is.matrix(Y)) stop("fmri_dataset must supply a matrix field (or support fmridataset accessors).")
  T <- nrow(Y)
  if (!is.null(runs_ds) && length(runs_ds) != T) stop("Dataset runs have incompatible length.")
  if (is.null(runs_ds)) runs_ds <- x$runs %||% attr(x, "runs")
  if (!is.null(runs_ds) && length(runs_ds) != T) stop("Dataset runs have incompatible length.")
  if (!is.null(runs_ds)) runs_ds <- as.integer(runs_ds)
  TR_ds <- TR_ds %||% x$TR %||% attr(x, "TR")
  TR_fd <- tryCatch(mppi_tr(fd), error = function(e) NULL)
  if (!is.null(TR_ds) && !is.null(TR_fd) && abs(TR_ds - TR_fd) > 1e-6) {
    stop("TR mismatch between dataset and design.", call. = FALSE)
  }
  tr <- TR_ds %||% TR_fd
  if (domain != "bold" && !is.null(tr)) {
    deconv <- deconv %||% list()
    supplied <- deconv$tr %||% deconv$TR
    if (!is.null(supplied) && (length(supplied) != 1L || abs(supplied - tr) > 1e-6)) stop("TR mismatch between dataset and deconvolution.", call. = FALSE)
    deconv$tr <- tr
  }
  # As in mppi(design, Y): no basis unless one is requested or carried by the design.
  basis_obj <- if (is.null(basis)) NULL else .mppi_resolve_basis(Y, basis = basis, r = r, dataset = x, design = fd)
  fit <- mppi(fd, Y = Y, basis = basis_obj, domain = domain,
              prewhiten = prewhiten, runs = runs_ds, deconv = deconv,
              engine = engine, instant = instant, ...)
  fit$dataset <- x
  if (!is.null(basis)) {
    basis_source <- if (is.null(basis_obj)) {
      if (is.character(basis)) tolower(basis[[1L]]) else if (is.list(basis) && !is.null(basis$type)) tolower(basis$type) else "roi"
    } else {
      basis_obj$source %||% "provided"
    }
    fit$basis_source <- basis_source
    if (!is.null(fit$basis)) fit$basis$source <- basis_source
  }
  fit
}

#' @export
mppi.fmri_frame <- function(x, fd, runs = NULL, deconv = NULL, ...) {
  if (!requireNamespace("fmridataset", quietly = TRUE)) stop("fmridataset is required.", call. = FALSE)
  Y <- as.matrix(fmridataset::collect_assay(x))
  schema <- fmridataset::temporal_schema(x)
  tr <- unique(as.numeric(schema$TR))
  if (length(tr) != 1L || !is.finite(tr) || tr <= 0 || !isTRUE(schema$contiguous)) stop("Frame requires regular, common-TR sampling.", call. = FALSE)
  if (!is.null(schema$censor) && any(schema$censor)) stop("Censored frames are unsupported by this adapter.", call. = FALSE)
  design_tr <- tryCatch(mppi_tr(fd), error = function(e) NA_real_)
  if (length(design_tr) == 1L && is.finite(design_tr) && abs(design_tr - tr) > 1e-6) stop("TR mismatch between frame and design.", call. = FALSE)
  deconv <- deconv %||% list()
  supplied <- deconv$tr %||% deconv$TR
  if (!is.null(supplied) && (length(supplied) != 1L || abs(supplied - tr) > 1e-6)) stop("TR mismatch between frame and deconvolution.", call. = FALSE)
  deconv$tr <- tr
  mppi(fd, Y = Y, runs = runs %||% schema$block_ids, deconv = deconv, ...)
}

#' @export
mppi.default <- function(x, ...) {
  stop(sprintf("No mppi() method for objects of class %s", paste(class(x), collapse = ", ")), call. = FALSE)
}

"%||%" <- function(a, b) if (is.null(a)) b else a

# Constructors --------------------------------------------------------------

# Print / summary ----------------------------------------------------------

#' @export
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

#' @export
print.mppi_fit <- function(x, ...) {
  cat(sprintf("mppi_fit: %d time points, %d regions, %d regressors\n",
              x$n_used, ncol(x$R), length(x$names)))
  basis_txt <- if (!is.null(x$basis)) sprintf(", basis=%s (r=%d)", x$basis$name, x$basis$r) else ""
  domain_txt <- sprintf(", domain=%s", x$domain %||% "bold")
  cat(sprintf("scale = %s, backend = %s, packed = %s%s%s\n",
              x$scale, x$backend, isTRUE(x$packed), basis_txt, domain_txt))
  energies <- vapply(x$Delta, .mppi_frob2, numeric(1))
  cat("Frobenius energy per regressor:\n")
  cat(paste0("  ", names(energies), ": ", format(energies, digits = 3)), sep = "\n")
  first <- x$Delta[[1]]
  if (!is.null(first)) {
    full <- if (is.list(first)) .mppi_unpack_upper(first) else first
    if (!all(is.na(full))) {
      eig <- eigen((full + t(full))/2, only.values = TRUE)$values
      cat("Top eigenvalues (regressor 1): ", paste(format(head(eig, 3), digits = 3), collapse = ", "), "\n")
    }
  }
  invisible(x)
}

#' @export
summary.mppi_fit <- function(object, ...) {
  energies <- vapply(object$Delta, .mppi_frob2, numeric(1))
  res <- list(
    timepoints = object$n_used,
    regions = ncol(object$R),
    regressors = length(object$names),
    scale = object$scale,
    backend = object$backend,
    packed = isTRUE(object$packed),
    domain = object$domain %||% "bold",
    basis = if (!is.null(object$basis)) object$basis$name else NULL,
    basis_rank = if (!is.null(object$basis)) object$basis$r else NULL,
    energy = energies
  )
  class(res) <- "summary.mppi_fit"
  res
}

#' @export
print.summary.mppi_fit <- function(x, ...) {
  cat(sprintf("mppi_fit summary: T=%d, V=%d, K=%d\n",
              x$timepoints, x$regions, x$regressors))
  base_line <- if (!is.null(x$basis)) sprintf(", basis=%s (r=%d)", x$basis, x$basis_rank) else ""
  cat(sprintf("scale=%s, backend=%s, packed=%s%s, domain=%s\n",
              x$scale, x$backend, x$packed, base_line, x$domain))
  cat("Frobenius energy:\n")
  cat(paste0("  ", names(x$energy), ": ", format(x$energy, digits = 3)), sep = "\n")
  invisible(x)
}

# HRF grouping --------------------------------------------------------------

#' Group HRF basis columns belonging to each condition.
#'
#' Columns are grouped only when they are recognisably the basis functions of
#' one condition: fmridesign names basis function `j` of a multi-basis HRF
#' `<condition>_b<j>` (zero-padded, e.g. `cond_cond.A_b01`, `cond_cond.A_b02`),
#' so a non-empty stem whose suffixes run exactly `1..n` (`n >= 2`) forms one
#' group; so do two or more distinct legacy prefixes `hrf`, `hrfT`, `hrfD`,
#' `hrfTD` of the same full condition (`hrf_A`, `hrfTD_A`, or `hrf(A)`). Every
#' other column, including unnamed or duplicated ones, is its own group;
#' distinct conditions are never merged.
#'
#' @param design An `mppi_design` object or design matrix with column names.
#' @return A named list mapping condition keys to column indices, in design order.
#' @export
mppi_group_hrf_columns <- function(design) {
  if (inherits(design, "mppi_design")) {
    idx <- design$psych_idx
    names_vec <- colnames(design$X)[idx]
  } else if (is.matrix(design)) {
    idx <- seq_len(ncol(design))
    names_vec <- colnames(design)
    if (is.null(names_vec)) stop("Provide column names when design is a matrix.")
  } else {
    stop("design must be an mppi_design or matrix with column names.")
  }
  n <- length(idx)
  if (is.null(names_vec)) names_vec <- rep("", n)
  names_vec[is.na(names_vec)] <- ""
  stem <- tag <- rep(NA_character_, n)
  fd <- regexpr("_b[0-9]{2,}$", names_vec)
  hit <- fd > 0
  stem[hit] <- paste0("b:", substr(names_vec[hit], 1L, fd[hit] - 1L))
  tag[hit] <- as.integer(substring(names_vec[hit], fd[hit] + 2L))
  re <- "^hrf([TtDd]*)(?:[._:](.+)|\\((.+)\\))$"
  lg <- !hit & grepl(re, names_vec, perl = TRUE)
  stem[lg] <- paste0("h:", sub(re, "\\2\\3", names_vec[lg], perl = TRUE))
  tag[lg] <- toupper(sub(re, "\\1", names_vec[lg], perl = TRUE))
  grp <- seq_len(n)
  key <- ifelse(nzchar(names_vec), names_vec, paste0("col", idx))
  for (s in unique(stem[!is.na(stem)])) {
    pos <- which(stem == s)
    cond <- substring(s, 3L)
    ok <- if (startsWith(s, "b:")) identical(sort(as.integer(tag[pos])), seq_along(pos)) else !anyDuplicated(tag[pos])
    if (nzchar(cond) && length(pos) >= 2L && ok) {
      grp[pos] <- pos[1L]
      key[pos] <- cond
    }
  }
  ids <- unique(grp)
  out <- lapply(ids, function(g) idx[grp == g])
  names(out) <- make.unique(key[ids])
  out
}
