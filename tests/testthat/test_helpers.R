library(testthat)

skip_if_not_installed("fmridesign")

suppressPackageStartupMessages(library(fmridesign))

set.seed(42)

Tn <- 60L
sf <- sampling_frame(blocklens = Tn, TR = 1)
enc_onsets <- seq(5, 45, by = 10)
rec_onsets <- seq(10, 50, by = 10)
trial_tbl <- data.frame(
  onset = c(enc_onsets, enc_onsets + 2, rec_onsets, rec_onsets + 2),
  cond = factor(c(rep("enc_vid01", length(enc_onsets)),
                  rep("enc_vid02", length(enc_onsets)),
                  rep("rec_vid01", length(rec_onsets)),
                  rep("rec_vid02", length(rec_onsets)))),
  run = 1
)
trial_tbl <- trial_tbl[order(trial_tbl$onset), ]
cond_names <- levels(trial_tbl$cond)
emod <- event_model(onset ~ hrf(cond), data = trial_tbl,
                    block = ~ run, sampling_frame = sf)
design <- mppi_model(emod, runs = rep(1L, Tn), include_intercept = FALSE)
runs <- design$runs
S_mat <- mppi_psych(design, T = Tn, TR = 1)

X <- cbind(`(Intercept)` = 1, S_mat)
Y <- matrix(rnorm(Tn * 8, sd = 0.5), Tn, 8)
fit_obj <- mppi_fit(Y, X, psych_idx = 2:5, runs = runs)


test_that("mppi_parametric produces gPPI-style modulators", {
  mods <- list(vivid = seq_len(Tn) / Tn)
  pk_mod <- mppi_parametric(design, T = Tn, mods = mods, select = "^cond_cond.rec")
  expect_named(pk_mod, c("cond_cond.rec_vid01::vivid", "cond_cond.rec_vid02::vivid"))
  expect_equal(length(pk_mod), 2L)
  expect_equal(length(pk_mod[[1]]), Tn)
  expect_true(all(abs(colMeans(do.call(cbind, pk_mod))) < 1e-6))
})


test_that("mppi_reinstatement returns within/between similarities", {
  rein <- mppi_reinstatement(fit_obj,
                             enc_pattern = "^cond_cond.enc_vid",
                             rec_pattern = "^cond_cond.rec_vid[0-9]{2}$",
                             mode = "raw")
  expect_length(rein$within, 2L)
  expect_length(rein$between, 1L)
  expect_true(is.finite(rein$stat))
})


test_that("mppi_sdppi yields per-trial summaries", {
  sd_df <- mppi_sdppi(fit_obj, design, select = "^cond_cond.rec",
                      pre_window = 3L, lags = -1:1, mode = "raw")
  expect_s3_class(sd_df, "data.frame")
  expect_gt(nrow(sd_df), 0L)
  expect_true(all(c("condition", "trial", "onset", "magnitude", "gain", "routing") %in% names(sd_df)))
})


test_that("mppi_precision_gate returns per-condition deltas", {
  gate <- mppi_precision_gate(fit_obj, mode = "raw", ridge = 1e-4)
  expect_equal(length(gate$Delta), length(cond_names))
  expect_s3_class(gate$summary, "data.frame")
  expect_true(all(c("condition", "dtheta_energy", "offdiag_ratio") %in% names(gate$summary)))
})
