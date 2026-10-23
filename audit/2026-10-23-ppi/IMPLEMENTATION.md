# multiPPI repair implementation and evidence

This report follows the [original audit](REPORT.md) and [authorized repair plan](FIX_PLAN.md). The baseline was HEAD `9d199cd53d484317c94a3cbaa6dfee972cb47f99` plus the pre-existing dirty/untracked work recorded in [the original manifest](source-manifest.json). Those original receipts and unrelated files were preserved. No commit, push, merge, site publication, or hosted CI execution is claimed.

## Supported conventional workflow

`ppi_design()` → `ppi_fit()` → `ppi_contrast()` now implements conventional seed-based PPI/gPPI. It is separate from the existing `mppi_*` covariance-modulation estimand. The [README](../../README.md) and [runnable conventional vignette](../../vignettes/vignette_09_conventional_ppi.Rmd) describe use.

The design accepts a raw BOLD seed, neural psychological inputs or events, nuisance regressors, run identities and censoring on the original acquisition grid. It uses fmrihrf by default, preserves the HRF's zeros and sample spacing, performs run-local ridge DCT deconvolution with masked observations, and reconvolves neural interaction terms before final-model whitening. Task and physiological main effects are included. The default adjusts the seed for nuisance conditional on task effects, preserving seed units. A categorical fmridesign event-model adapter preserves TR, run and sampling-offset metadata; unsupported modulated, normalized, subset or multiple-basis constructions are rejected explicitly.

Runs have separate scientific coefficients by default; shared coefficients are an explicit alternative. Nuisance terms and intercepts always remain run-specific. Condition-named contrasts average interactions equally across runs; complete design-column names allow explicit weights. Column-equilibrated SVD fitting retains alias information, uses rank for residual degrees of freedom, and reports nonestimable coefficients/contrasts as unavailable. Contrast estimability and joint F statistics are invariant to contrast-row scaling.

The default noise provider is fmriAR. Estimated noise uses full-model residuals; the same fitted operator transforms targets and the complete design. A saved fmriAR plan supports conditional fixed-noise fits. Shared global/run operators allow target chunks to reuse the factorization. Parcel labels allow separate group operators; one label per target requests targetwise estimation. IID censoring retains original time coordinates. **AR with censor gaps is rejected**, because the exact marginal observed-data whitening contract has not been implemented. A gap is never silently treated as a new independent run.

## Repair-plan disposition

| Work | Delivered evidence and remaining scope |
|---|---|
| F01: harness/contracts | Standard installed-package test harness; explicit estimands, sampling/seed contracts and rank policy; original failures retained. |
| F02: adapters/time | Current fmridataset/fmrireg APIs, TR/run propagation, rank-safe intercept handling, raw-seed equivalence, original-distance censor/lag pairs; real dependency and NeuroVec fixtures. |
| F03: deconvolution/construction | Dense forward/inverse witnesses, run isolation, masked inverse, HRF timing, nuisance-addition invariance, actual SPM comparison; no SPM-prior equivalence or latent-neural recovery claim. |
| F04: fitting/inference | IID/fixed/estimated/grouped fmriAR fitting, named t/F contrasts, independent GLS and Nilearn witnesses, fixed-design statistical screens. Censored AR remains unavailable. |
| F05: covariance/derived methods | Signed low-rank reconstruction, contiguous held-out CV blocks, directed lag storage, partial-correlation derivative, variance decomposition, task-main-effect evidence labeling, synchronized HRF-ensemble metadata and caches. |
| F06: legacy inference | Corrected randomization weights/boundaries/diagonals and streaming counts. Studentization, Freedman–Lane and seed-set permutation fail explicitly pending a qualified reduced model. Remaining covariance randomizations and lag selection are experimental, not calibrated procedures. |
| F07: broader APIs | Descriptive seed-CCA, real NeuroVec orientation, current rest/frame extraction, one-observation-per-subject limma fitting and named group alignment. Neural seed-CCA, CCA inference and clustered-volume extraction are unavailable. Existing beta/band/rest/instant/mechanistic families have executed regression checks, not universal scientific qualification. |
| F08: efficiency | Chunked conventional fitting, shared factors, cached forward convolution, streaming randomization counts, exact banded MAP solve. Measured IID seed-map and MAP workloads below; the plan's entire performance Cartesian product was not run. |
| F09: package/docs | Updated dependencies, exports, license format, portable Makevars, Rd and executable vignettes; local full-check receipts recorded below. Hosted cross-platform CI and site publication remain unexecuted. |

## Independent external comparisons

The [external comparison report](external-comparison.md) retains all scripts, synthetic fixtures, dependency versions, source hashes and adverse differences. Nilearn 0.14.1 checks six completed-design cases: IID, censored IID, rank deficiency, fixed AR(1), extreme column units, and grouped AR. Maximum discrepancies are below `2e-13` across the reported metrics. The rank-deficient case explicitly corrects Nilearn's column-count dispersion to rank-based dispersion; AR filtering explicitly matches fmriAR's stationary first-row boundary. The grouped comparison holds fitted operators fixed. These are conditional regression witnesses, not independent validation of noise estimation or interaction construction.

Actual `spm_peb_ppi` from pinned official SPM source ran under Octave. On the bounded 40-scan fixture, task convolution agrees to `7.63e-17`, the default adjusted physiological input to `4.44e-16`, and interaction correlation is `0.99916` (RMSE `0.02325`). Neural estimates differ substantially (RMSE `0.1832`). SPM uses a padded PEB inverse and detrending; multiPPI uses unpadded fixed ridge. The methods are deliberately not described as equivalent.

## Conditional statistical screens

Final conventional source SHA-256 values are `a710631460dac0ec6713b0a005d9d6d57a1e332bc1d3c67bb01de0f31d46f641` (`ppi_design.R`) and `f9255db93f93bf855b4862a5cc7c4184fa303f78afd5b7c06b7c3cd886e7f930` (`ppi_fit.R`). Start/end hashes match in the validation receipts. R 4.5.1, fmriAR 0.4.2 and fmrihrf 0.4.0.9000 were used.

The [fixed-design results](implementation-validation-final-fixed-design.csv) contain 2,000 independent null draws and 500 nonzero-effect draws per scenario. For four prespecified primary scenarios, simultaneous one-sided 95% binomial upper bounds meet the declared tolerated-inflation screen of 7.5% at nominal alpha 5%:

| Scenario | Null rejection | Simultaneous upper bound | Null interval coverage | Effect power |
|---|---:|---:|---:|---:|
| IID | 98/2000 = 4.90% | 6.09% | 95.10% | 52.4% |
| IID with censoring | 102/2000 = 5.10% | 6.31% | 94.90% | 51.6% |
| Fixed AR(1) | 102/2000 = 5.10% | 6.31% | 94.90% | 36.0% |
| Estimated AR(1) | 117/2000 = 5.85% | 7.14% | 94.15% | 41.2% |

The heterogeneous-target **adverse pooled-AR case fails**: 236/2000 = 11.8% rejection, 95% binomial CI [10.42%, 13.30%], 88.2% coverage. This evidence must remain visible. The [separate targetwise screen](implementation-validation-final-grouped.csv) gives 112/2000 = 5.6%, CI [4.63%, 6.70%], one-sided upper 6.52%, and 94.4% coverage. It supports the tested targetwise alternative; it does not establish arbitrary-noise robustness. Grouping targets does not solve unmodeled run-specific innovation variance. Estimated-noise inference remains approximate and conditional on the seed/design.

## Independent neural-to-BOLD construction screen

The [generative report](implementation-generative-report.md) and [script](../../inst/validation/validate_ppi_generative.R) add an independently constructed dense forward model, two unequal runs, a deliberately truncated inverse, 200 neural realizations and ten target-noise draws per realization. Both activation-only and constant-coupling matched nulls have negligible noiseless interaction projections and reject 4.7% of targets. The true neural interaction coefficient is 0.5, while its average projection onto the constructed design is **0.6673** (seed-cluster bootstrap 95% interval [0.6540, 0.6811]). Matched-effect rejection is 53.85%; the paired effect-minus-constant-null difference has interval [0.4595, 0.5250]. Independent QR projection agrees with the implementation within `2.39e-12`.

This screen exposes substantial construction distortion: residualized estimated-versus-oracle interaction mean correlation is 0.7542 and mean relative difference is 0.659. A one-second target-HRF delay raises the constant-coupling null rejection to **7.95%**. Seed measurement noise and delayed-HRF effect results are retained separately. Confidence intervals target the constructed-design projection, not the generating 0.5 neural coefficient; deterministic model-misspecification residuals are also recorded. The uncertainty summary resamples neural-seed clusters, retaining the ten nested targets, rather than pretending all targets independently redraw a seed. These results support conditional PPI and the stated construction sanity checks, while expressly failing to establish unbiased neural effect recovery or HRF-mismatch robustness. The first incomplete generator run is preserved and marked superseded; the corrected seven-scenario recipe was frozen before its final run.

The matched-effect mean also has a nonzero residual relative norm of 0.0534; its projection coverage is empirical under mean misspecification. Identical matched-null rejection counts are expected because both means lie in the main-effect span and share the same target-noise draws. They are not independent calibration confirmations. Independent scientific review found no critical defect in the corrected generator, projection oracle or clustered bootstrap; the broader scientific qualification remains incomplete.

## Equivalent-output performance

On this macOS arm64 machine, R 4.5.1/Accelerate, T=600, K=3, microtime=16 and one seed, the [fresh-process receipts](implementation-final-benchmark-rss.json) report:

| Targets | Design seconds | First fit seconds | Repeated fit seconds | Process peak RSS |
|---|---:|---:|---:|---:|
| 400 | 0.881 | 0.104 | 0.003–0.004 | 0.331 GiB |
| 1,000 | 0.894 | 0.105 | 0.007–0.008 | 0.369 GiB |
| 50,000 | 0.881 | 0.532 | 0.340–0.344 | 2.004 GiB |

These are individual IID workload measurements, not benchmark distributions or estimates for AR/targetwise fitting. Whole-process peak RSS includes input generation and three simultaneous returned fits, measured with `wait4`; it is not incremental solver allocation. Changing chunks between 256 and 2048 preserves coefficients, standard errors and residual variance to tolerance `1e-10`.

The [native MAP reference](implementation-native/map-benchmark.log) compares banded solving against a separately materialized dense causal-convolution/first-difference system, two targets, HRF length 12, T=256–4096. Maximum absolute error is `1.56e-15`. At T=4096 the measured banded call took 0.009 seconds and dense construction/solve 1.364 seconds. This is a single bounded reference measurement, not a general speedup claim. Replay with `Rscript inst/benchmarks/benchmark_map.R /path/to/frozen/source`.

## Package qualification and remaining limits

The fresh source snapshot passed `R CMD build --no-manual` and `R CMD check --as-cran --no-manual`, both exit zero, using `LC_ALL=C`, `_R_CHECK_FORCE_SUGGESTS_=true`, `NOT_CRAN=true`, and the complete dependency library. The [full check log](implementation-package/check-release.log) reports **zero errors, zero warnings, two notes**. The [installed test output](implementation-package/testthat-release.Rout) reports **622 passes and one skip**; the skip checks a missing-fmrireg error branch when fmrireg is installed. PWI cases executed. All package vignettes were rebuilt. The two notes concern incoming CRAN metadata/non-mainstream suggested packages and the host's inability to verify current time. These are not numerical test failures. [Source hashes](implementation-package/release-source-sha256.txt) and the [built tarball hash](implementation-package/release-tarball-sha256.txt) bind the receipt to the tested snapshot; production R and C++ sources matched the workspace after checking. Temporary build/check artifacts were removed after retaining logs.

Important limits remain explicit: no exact AR censor-gap support; no universal HRF or seed-measurement-error robustness; no causal or latent-neural coupling claim; no seed-overlap validity claim; conditional CV does not establish a fully nested prediction pipeline; legacy resampling is not generally calibrated; no multiple-comparison correction is supplied by conventional t/F contrasts; no executed hosted/platform release gate or regenerated published site. The pre-existing generated `docs/` tree was preserved and must not be mistaken for the updated source documentation.
