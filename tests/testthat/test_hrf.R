context("HRF grouping and adaptive blending")

test_that("HRF grouping matches expected names", {
  skip_if_not_installed("fmridesign")
  library(fmridesign)
  sframe <- sampling_frame(blocklens = 20, TR = 1)
  ev <- event_model(onset ~ hrf(cond) + hrfTD(cond),
                    data = data.frame(onset = c(0,5,10,15), cond = factor(c("A","A","B","B"))),
                    block = ~1, sampling_frame = sframe)
  D <- as_mppi_design(ev)
  groups <- mppi_group_hrf_columns(D)
  expect_true("cond" %in% names(groups))
  expect_equal(length(groups[["cond"]]), 2L)
})

test_that("HRF adaptive weights recover dominant basis", {
  skip_if_not_installed("fmridesign")
  library(fmridesign)
  V <- 4L
  Delta1 <- matrix(0, V, V); Delta1[1,2] <- Delta1[2,1] <- 3
  Delta2 <- matrix(0, V, V); Delta2[3,4] <- Delta2[4,3] <- 0.5
  deltas <- list(Delta1, Delta2)
  combo <- mppi_hrf_adapt(deltas, pk_norms = c(1, 1))
  expect_equal(which.max(abs(combo$weights)), 1L)
})
