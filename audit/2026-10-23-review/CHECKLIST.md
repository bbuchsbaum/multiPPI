# multiPPI review checklist

Source: the review of HEAD `fd4f82d` (2026-10-23). Each phase is one priority
level, worked in order. Every fix lands with a regression test that fails on
`fd4f82d` and passes after the fix, unless the item is packaging or deletion.

Legend: `[ ]` open, `[x]` done, `[~]` resolved by failing closed (explicit error
instead of a silently wrong result), `[-]` deferred with a reason.

## Phase 1 — P0: install and headline-workflow blockers

- [x] 1.1 `src/Makevars`, `src/Makevars.win`: link `$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` (Linux load fails with `undefined symbol: dpotrf_`).
- [x] 1.2 Remove committed build objects (`src/*.o`, `src/multiPPI.so`, `..Rcheck/`) and add a `.gitignore` so fresh clones build.
- [x] 1.3 Headline defaults (`ppi_design` HRF via fmrihrf, `ppi_fit` noise via fmriAR) must work from a documented install, or fail with an actionable message; README example must run.

## Phase 2 — P1: conventional PPI correctness

- [x] 2.1 `ppi_design`: deconvolution ridge invariant to HRF scale, microtime and TR (currently absolute `lambda = 0.1`).
- [x] 2.2 `ppi_design`: default `n_basis` must not low-pass long runs (currently `min(64, scans)`).
- [x] 2.3 `ppi_contrast`: condition-named contrasts average only over runs where the condition occurs (currently `NA` if absent from any run).

Phase 2 note: `deconv$lambda` is now relative to the mean diagonal of the
observed basis Gram matrix (default 1, chosen by simulation across TR, run
length, neural timescale and seed noise), and `n_basis` defaults to one
component per scan. Receipts under `audit/2026-10-23-ppi/` and validation
scripts that pass an explicit `lambda` were produced under the old absolute
definition; rerunning them will not reproduce those numbers exactly.

## Phase 3 — P2: silently wrong scientific output (legacy and peripheral)

- [ ] 3.1 Neural-domain HRF grouping merges distinct conditions (`mppi_group_hrf_columns`, `mppi_s3.R:319`).
- [ ] 3.2 Basis fits zero the diagonal and normalize in component space before back-projection (`mppi_fit.R:342`, `mppi_basis.R:181`, `mppi_seedset.R:58`).
- [ ] 3.3 `mppi.fmri_dataset` silently picks a PCA basis and drops the dataset TR (`mppi_s3.R:191`).
- [ ] 3.4 Instant engine: run-local deconvolution, censoring rejection, rank-safe baseline, rank df, lag blocks on the acquisition grid (`mppi_engine_state.R`, `mppi_instant.R`).
- [ ] 3.5 Instant engine EWM output is a final filter state, not a slope, with transposed lag orientation (`mppi_instant.R:103`, `src/instant.cpp`).
- [ ] 3.6 `mppi_precision_gate` disagrees with `fit$partial`; normalized mode mixes units (`mppi_helpers.R:279`).
- [ ] 3.7 `mppi_fit_whitened` adds a raw intercept to the whitened design, so it is not GLS (`mppi_fit_whitened.R:59`).
- [ ] 3.8 Routing/asymmetry indices are identically zero (`mppi_neural.R:302`, `mppi_adv_methods.R:168,231`, `mppi_axes.R:59`).
- [ ] 3.9 Seed-set `aggregate = "svd"` (default) is a sign-arbitrary between-seed difference (`mppi_seedset.R:163`).
- [ ] 3.10 Rest state modulators are degenerate after residualizing one-hot states on each other plus intercept (`mppi_rest.R:60,203`).
- [ ] 3.11 `mppi_fit_multi_band` filters Y but not X; `mppi_bandpass` ignores `low` without `signal`; check `butter` cutoff units (`mppi_band.R`).
- [ ] 3.12 Mechanistic helpers use pre-innovation `fit$U` for innovations fits (`mppi_adv_methods.R:500`).
- [ ] 3.13 `mppi_hrf_adapt` combined matrix has arbitrary sign (`mppi_hrf_adapt.R:109`).
- [ ] 3.14 `mppi_beta_permute` default null is invariant, p = 1 (`mppi_beta.R:189`).
- [ ] 3.15 `mppi_rest_trait_cv`: no intercept or in-fold standardization; subject-named `y` gives NaN (`mppi_rest.R:352,372`).
- [ ] 3.16 PWI phase and frequency computed across concatenated runs; character runs give NA (`ipppi.R:118`).
- [ ] 3.17 `mppi_coerce_beta` turns NA into 0 and sums duplicate cells (`mppi_beta.R:10`).

## Phase 4 — P3: crashes or wrong results on plausible paths

- [ ] 4.1 Unseeded randomization restores `.Random.seed`, so repeated calls give identical nulls (`mppi_inference.R:69,104`).
- [ ] 4.2 Two incompatible HRF-function conventions; default HRFs differ in normalization (`mppi_fit.R:47`, `mppi_helpers.R:316`).
- [ ] 4.3 Unnamed designs: lagged result names collide; regressors named `""` (`mppi_lag.R:46`, `mppi_fit.R:190`).
- [ ] 4.4 Factor run labels crash `rle()` (`mppi_fit.R:161` and others).
- [ ] 4.5 `mppi_sdppi` counts each event several times on convolved regressors (`mppi_helpers.R:195`).
- [ ] 4.6 sdPPI pre/post windows cross run boundaries (`mppi_adv_methods.R:302`).
- [ ] 4.7 `mppi_state_predict`: gaussian path crashes, factor `y` crashes, character runs give NA (`mppi_adv_methods.R:400`).
- [ ] 4.8 `mppi_seedset_slopes(lag != 0)` errors for `scale = "cov"/"corr"` fits (`mppi_seedset.R:107`).
- [ ] 4.9 `mppi_group_seedset` reports design column 1 (intercept) (`mppi_seedset.R:439`).
- [ ] 4.10 `mppi_group_glm` rejects plain `model.matrix()` designs (`mppi_group_ebayes.R:183`).
- [ ] 4.11 `mppi_sim_rest_dataset` pads the wrong end, leaving leading NA rows (`mppi_sim_rest.R:199`).
- [ ] 4.12 `ewm_gated_cross_rcpp` clamps out-of-range lag indices instead of skipping (`src/instant.cpp:189`).
- [ ] 4.13 `inst_corr_rcpp` treats NaN as data; one NaN poisons the rest of the series (`src/instant.cpp:302`).
- [ ] 4.14 `seed_ppi_mcca` silently drops out-of-range neighbourhood indices (`seed_ppi_mcca.R:141`).
- [ ] 4.15 Frame adapters read `schema$block_ids` where `ppi_fit` uses `schema$run_ids` (`mppi_s3.R:221`, `mppi_fit.R:533`) — verify against the fmridataset contract.
- [ ] 4.16 Neural rest modulators deconvolve concatenated runs (`mppi_rest.R:275`).
- [ ] 4.17 Censored fits store compressed run lengths in `fit$lag_blocklens` (`mppi_fit.R:160`).

## Phase 5 — P4: efficiency

- [ ] 5.1 `mppi_omnibus`: exact precomputed statistic instead of O(B·T·V²) per draw.
- [ ] 5.2 `mppi_lag_select`: do not default to `scale = "cov"` (V³ inverse and unused diagnostics).
- [ ] 5.3 `ppi_fit` grouped/targetwise: share the unwhitened factorization and residuals across groups.
- [ ] 5.4 `ppi_design`: batch the basis convolutions instead of per-column FFT calls.
- [ ] 5.5 `mppi_classify` innovations mode recomputes identical innovations per draw.
- [ ] 5.6 `mppi_sdppi_fit` multiplies the full T×r matrix per trial instead of the window rows.

## Phase 6 — P5: quality and hygiene

- [ ] 6.1 NAMESPACE hand-maintained and drifted from roxygen (`mppi_deconv_dct`, `pwf_stack`).
- [ ] 6.2 Remove generated clutter (`repomix-output.*`, rendered `vignettes/*.html`) and settle the remaining loose files.
- [ ] 6.3 Dead or stale code: `src/init.cpp` (`R_init_mppiDeconv`), `.mppi_phase_randomize`, `.mppi_with_intercept`, duplicate `print.mppi_design`, unreachable branches.
- [ ] 6.4 `deconv_dct_multi` silently ignores a wrong-length `q`.
- [ ] 6.5 `mppi_rank_cv` loss includes the diagonal although fits zero it.
- [ ] 6.6 Replace mock-only tests that hid real bugs (routing, state prediction, sdPPI) with tests on fitted objects.
