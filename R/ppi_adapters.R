#' Construct a PPI design from a categorical fmridesign event model
#'
#' Uses event-level onsets, durations and condition labels, not the convolved
#' event design as a neural input. Only unmodulated categorical, single-HRF
#' event terms are supported. Unsupported terms require explicit `events` or
#' neural `psych` input to [ppi_design()]. The builder's microtime binning and
#' causal Riemann convolution may differ from fmridesign's numerical integration.
#'
#' @param seed Numeric raw BOLD vector on the complete acquisition grid.
#' @param event A categorical `fmridesign` event_model.
#' @param baseline Optional finite numeric matrix or fmridesign baseline_model;
#'   its columns are supplied as nuisance regressors.
#' @param ... Additional arguments to [ppi_design()], except `psych`, `events`,
#'   `TR`, `runs`, `nuisance`, `hrf`, and `sampling_offset`, which are determined
#'   by the model. `microtime` must represent the model's sampling offset exactly.
#' @return A `ppi_design` object.
#' @export
ppi_design_from_model <- function(seed, event, baseline = NULL, ...) {
  if (!inherits(event, "event_model")) stop("event must be an fmridesign event_model")
  dots <- list(...)
  forbidden <- intersect(names(dots), c("psych", "events", "TR", "runs", "nuisance", "hrf", "sampling_offset"))
  if (length(forbidden)) stop("model determines: ", paste(forbidden, collapse = ", "))
  frame <- event$sampling_frame
  if (!inherits(frame, "sampling_frame")) stop("event model has no supported sampling frame")
  lens <- fmridesign::blocklens(frame)
  if (sum(lens) != length(seed)) stop("seed length does not match the sampling frame")
  samples <- fmridesign::samples(frame)
  run <- rep(seq_along(lens), lens)
  times <- split(samples, run)
  offsets <- unname(vapply(times, `[`, numeric(1), 1L))
  steps <- unlist(lapply(times, diff), use.names = FALSE)
  if (!length(steps) || any(!is.finite(steps)) || any(abs(steps - steps[1L]) > 1e-8) || any(abs(offsets - offsets[1L]) > 1e-8)) stop("sampling frame must have a common regular TR and sampling offset")
  tr <- steps[1L]
  terms <- stats::terms(event)
  if (!length(terms)) stop("event model contains no event terms")
  tables <- list(); kernels <- list()
  for (j in seq_along(terms)) {
    term <- terms[[j]]
    fail <- function() stop("unsupported continuous, modulated, subset, normalized, or multibasis model term; supply explicit events or psych to ppi_design")
    if (!inherits(term, "event_term")) fail()
    if (!length(term$events) || !all(vapply(term$events, fmridesign::is_categorical, logical(1)))) fail()
    tab <- fmridesign::event_table(term)
    if (!ncol(tab) || !all(vapply(tab, function(x) is.factor(x) || is.character(x), logical(1)))) fail()
    if (!is.null(term$subset) && !all(term$subset)) fail()
    spec <- attr(term, "hrfspec")
    if (is.null(spec) || !is.function(spec$hrf) || fmrihrf::nbasis(spec$hrf) != 1L || isTRUE(spec$normalize) || !isTRUE(spec$summate) || !is.null(spec$hrf_fun)) fail()
    ev <- fmridesign::events(term)
    if (!all(c("onset", "duration", "block", "condition") %in% names(ev)) || anyNA(ev$condition)) fail()
    tables[[j]] <- data.frame(onset = ev$onset, duration = ev$duration,
                              condition = as.character(ev$condition), run = as.character(ev$block))
    kernels[[j]] <- spec$hrf
  }
  grid <- seq(0, 32, length.out = 257L)
  reference <- kernels[[1L]](grid)
  if (!is.numeric(reference) || is.matrix(reference) || any(!is.finite(reference)) ||
      any(!vapply(kernels, function(h) isTRUE(all.equal(as.numeric(h(grid)), as.numeric(reference), tolerance = 1e-12)), logical(1)))) stop("all terms must share one scalar HRF; supply explicit events or psych to ppi_design")
  ev <- do.call(rbind, tables)
  if (anyDuplicated(unlist(lapply(tables, function(x) unique(x$condition))))) stop("model terms have overlapping condition identities; supply explicit events or psych to ppi_design")
  if (inherits(baseline, "baseline_model")) baseline <- as.matrix(fmridesign::design_matrix(baseline))
  args <- c(list(seed = seed, events = ev, TR = tr, runs = run, nuisance = baseline,
                 hrf = kernels[[1L]], sampling_offset = offsets[1L]), dots)
  do.call(ppi_design, args)
}
