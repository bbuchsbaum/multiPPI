test_that("simulated rest dataset yields recoverable coupling", {
  set.seed(42)
  sim <- mppi_sim_rest_dataset(T = 800, V = 50, r = 5,
                               amp = c(0.4, -0.35),
                               contexts = c("ctx1", "ctx2"),
                               obs_noise = 0.01)
  basis_mat <- sim[["loadings"]]
  mods <- sim[["modulators"]]
  fit <- mppi_rest(sim[["Y"]],
                   modulators = mods,
                   basis = basis_mat,
                   prewhiten = FALSE,
                   domain = "bold")
  # Ground truth from latent components (noise-free reference)
  fit_true <- mppi_rest(sim[["U"]],
                        modulators = mods,
                        basis = diag(ncol(sim[["U"]])),
                        prewhiten = FALSE,
                        domain = "bold")
  ctx_names <- names(mods)
  delta_est <- fit[["Delta"]]
  truth <- sim[["truth"]]
  for (nm in ctx_names) {
    est <- delta_est[[nm]]
    amp_val <- truth$amp[match(nm, names(truth$coupling))]
    ref <- truth$coupling[[nm]] * amp_val
    diag(ref) <- 0
    expect_equal(dim(est), dim(ref))
    cc <- stats::cor(as.vector(est), as.vector(ref), use = "pairwise.complete.obs")
    expect_gt(cc, 0.8)
  }
})

test_that("trait effect correlates with recovered gain", {
  skip_if_not_installed("fmriAR")
  set.seed(7)
  cohort <- mppi_sim_rest_cohort(n_subj = 8,
                                 contexts = c("ctx1", "ctx2"),
                                 amp_base = c(0.3, 0.2),
                                 trait_effect = c(0.5, 0),
                                 amp_noise = 0.01,
                                 obs_noise = 0.005,
                                 T = 800, V = 40, r = 5)
  datasets <- cohort[["datasets"]]
  fits <- lapply(datasets, function(ds) {
    mods <- ds[["modulators"]]
    basis_mat <- ds[["loadings"]]
    mppi_rest(ds[["Y"]],
              modulators = mods,
              basis = basis_mat,
              prewhiten = FALSE,
              domain = "bold")
  })
  feats <- mppi_rest_features(fits, aggregate = "none")
  ctx1_rows <- feats[feats$condition == "ctx1", ]
  trait <- cohort[["trait"]][ctx1_rows[["subj"]]]
  expect_gt(cor(ctx1_rows[["mag"]], trait), 0.8)
})

test_that("gain metric tracks template-directed coupling", {
  skip_if_not_installed("fmriAR")
  set.seed(2027)
  r <- 6
  v <- matrix(rnorm(r), r, 1)
  v <- v / sqrt(sum(v^2))
  gain_template <- v %*% t(v)
  null_template <- diag(r)
  templates <- list(gain = gain_template, control = null_template)
  cohort <- mppi_sim_rest_cohort(n_subj = 12,
                                 contexts = c("gain", "control"),
                                 amp_base = c(0.2, 0.2),
                                 trait_effect = c(0.6, 0),
                                 amp_noise = 0.005,
                                 obs_noise = 0.003,
                                 T = 900, V = 50, r = r,
                                 templates = templates)
  fits <- lapply(cohort$datasets, function(ds) {
    mppi_rest(ds[["Y"]],
              modulators = ds[["modulators"]],
              basis = ds[["loadings"]],
              prewhiten = FALSE,
              domain = "bold")
  })
  feats <- mppi_rest_features(fits, aggregate = "none")
  gain_rows <- feats[feats$condition == "gain", ]
  control_rows <- feats[feats$condition == "control", ]
  trait <- cohort$trait
  cor_gain <- cor(gain_rows$G, trait[gain_rows$subj])
  cor_control <- cor(control_rows$G, trait[control_rows$subj])
  expect_gt(cor_gain, 0.6)
  expect_gt(abs(cor_gain), abs(cor_control) + 0.2)
})
