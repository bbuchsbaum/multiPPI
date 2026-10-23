test_that("real categorical event models preserve sampling and event timing", {
  skip_if_not_installed("fmridesign")
  hrf <- fmridesign::hrf
  sf <- fmridesign::sampling_frame(c(12, 10), TR = 2)
  ev <- data.frame(onset = c(2, 8, 4), cond = factor(c("A", "B", "A")), run = c(1, 1, 2))
  model <- fmridesign::event_model(onset ~ hrf(cond), data = ev, block = ~run, sampling_frame = sf)
  d <- ppi_design_from_model(seq_len(22), model, run_effects = "shared")
  expect_equal(d$TR, 2)
  expect_equal(d$provenance$sampling_offset, 1)
  expect_equal(d$runs, as.character(rep(1:2, c(12, 10))))
  expect_equal(d$provenance$psych_grid, "events")
  baseline <- matrix(seq_len(22), ncol = 1, dimnames = list(NULL, "drift"))
  with_baseline <- ppi_design_from_model(seq_len(22), model, baseline = baseline, run_effects = "shared")
  nuisance <- with_baseline$X[, with_baseline$terms$type == "nuisance",drop=FALSE]
  expect_equal(ncol(nuisance),2L)
  expect_equal(unname(rowSums(nuisance)),as.numeric(baseline))
  expect_true(all(nuisance[1:12,2]==0) && all(nuisance[13:22,1]==0))
  task <- d$X[, d$terms$type == "task", drop = FALSE]
  ref <- as.matrix(fmridesign::design_matrix(model))
  expect_equal(unname(task), unname(ref), tolerance = 0.001)
  expect_error(ppi_design_from_model(seq_len(22), model, microtime = 3), "sampling_offset")
  expect_error(ppi_design_from_model(seq_len(22), model, TR = 1), "model determines")
  mod <- transform(ev, strength = c(1, 2, 3))
  modmodel <- fmridesign::event_model(onset ~ hrf(strength), data = mod, block = ~run, sampling_frame = sf)
  expect_error(ppi_design_from_model(seq_len(22), modmodel), "unsupported continuous")
})

test_that("sampling offset shifts samples without shifting the original grid", {
  p <- matrix(c(1, 0, 0, 0), ncol = 1, dimnames = list(NULL, "A"))
  d <- ppi_design(1:4, psych = p, TR = 2, microtime = 2, sampling_offset = 1,
                  hrf = c(0, 1, 0), center = FALSE)
  expect_equal(unname(d$X[, d$terms$type == "task"]), c(1, 0, 0, 0), tolerance = 1e-12)
  expect_error(ppi_design(1:4, psych = p, TR = 2, microtime = 2, sampling_offset = 0.5), "sampling_offset")
})
