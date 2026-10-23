context("Design interface generics")

skip_if_not_installed("fmridesign")

test_that("fmridesign helpers produce coherent bundles", {
  suppressPackageStartupMessages(library(fmridesign))
  set.seed(11)
  Tn <- 24L
  sf <- sampling_frame(blocklens = Tn, TR = 2)
  events <- data.frame(
    onset = c(seq(0, by = 4, length.out = Tn/4),
              seq(2, by = 4, length.out = Tn/4)),
    cond = factor(rep(c("taskA", "taskB"), each = Tn/4)),
    run = 1
  )
  events <- events[order(events$run, events$onset), ]
  emod <- event_model(onset ~ hrf(cond), data = events,
                      block = ~ run, sampling_frame = sf)
  runs_vec <- rep(1L, Tn)
  design <- mppi_model(emod, runs = runs_vec, include_intercept = FALSE)
  runs_vec <- design$runs
  S <- mppi_psych(design, T = Tn, TR = 2)
  X_task <- mppi_task(design)

  expect_equal(mppi_tr(emod), 2)
  expect_equal(mppi_runs(emod, T = Tn), runs_vec)
  expect_equal(mppi_psych(design, T = Tn, TR = 2), S)
  expect_equal(mppi_nuisance(design, T = Tn), matrix(0, Tn, 0))
  expect_true(is.function(mppi_hrf(design)))

})

test_that("mppi() handles fmridesign and fmri_frame inputs", {
  skip_if_not_installed("fmridataset")
  skip_if_not_installed("fmrireg")
  suppressPackageStartupMessages(library(fmridesign))
  set.seed(22)
  Tn <- 24L
  V <- 4L
  sf <- sampling_frame(blocklens = Tn, TR = 1.5)
  events <- data.frame(
    onset = c(seq(0, by = 4, length.out = Tn/4),
              seq(2, by = 4, length.out = Tn/4)),
    cond = factor(rep(c("taskA", "taskB"), each = Tn/4)),
    run = 1
  )
  events <- events[order(events$run, events$onset), ]
  emod <- event_model(onset ~ hrf(cond), data = events,
                      block = ~ run, sampling_frame = sf)
  runs_vec <- rep(1L, Tn)
  design <- mppi_model(emod, runs = runs_vec, include_intercept = FALSE)
  runs_vec <- design$runs
  X_task <- mppi_task(design)
  Y <- matrix(rnorm(Tn * V), Tn, V)

  manual <- mppi_fit(Y = Y, X = X_task, psych_idx = 1:2, runs = runs_vec, domain = "bold")
  via_design <- mppi(design, Y = Y, domain = "bold")

  manual_raw <- lapply(seq_along(manual$Delta), function(i) mppi_get_M_scaled(manual, i, mode = "raw"))
  design_raw <- lapply(seq_along(via_design$Delta), function(i) mppi_get_M_scaled(via_design, i, mode = "raw"))
  expect_equal(design_raw, manual_raw)
  expect_equal(via_design$pk, manual$pk)

  ds <- fmrireg::matrix_frame(datamat = Y, TR = 1.5, run_length = c(Tn))
  via_ds <- mppi.fmri_frame(ds, design, domain = "bold")

  ds_raw <- lapply(seq_along(via_ds$Delta), function(i) mppi_get_M_scaled(via_ds, i, mode = "raw"))
  expect_equal(ds_raw, manual_raw)
  expect_equal(via_ds$pk, manual$pk)

})
