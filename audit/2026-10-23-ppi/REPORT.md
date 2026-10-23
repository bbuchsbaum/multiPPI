# multiPPI assessment — 23 October 2026

**The repository does not yet meet the requested goal of correct, efficient conventional PPI.** It has a useful covariance-modulation engine, several sound numerical kernels, and a substantial collection of exploratory methods. Conventional seed-based PPI/gPPI needs a distinct regression path. Existing inference, temporal preprocessing, and package adapters have reproducible defects that should be repaired before scientific use of the affected paths.

This assessment covers HEAD `9d199cd53d484317c94a3cbaa6dfee972cb47f99` **plus the existing working-tree changes**, including the untracked `seed_ppi_mcca.R`. Those changes are part of the evaluated implementation; the findings are not attributed solely to committed HEAD. No package implementation was changed. The added files are this report and its audit evidence.

Source identity is recorded in [source-manifest.json](source-manifest.json), whose SHA-256 is `69a74abe09c20315bb8f420b222747a74971d608b95ef14d900544bc5df5452f`. [receipt.json](receipt.json) records commands and the built tarball hash. All data used in the probes were synthetic. Installed dependency versions and their recorded source revisions are in [dependency-versions.txt](dependency-versions.txt); the runtime was R 4.5.1 on arm64 macOS 14.3 with Accelerate BLAS. Thread-count environment overrides were unset.

## What the main estimator actually estimates

The implementation in `R/mppi_fit.R:164` and `R/mppi_fit.R:330` residualizes the regional signals against task main effects and nuisance regressors, residualizes each task column against the other design columns, and computes

\[
R=M_XY,\qquad p_k=M_{X_{-k}}X_k,\qquad
\Delta_{ij,k}=\frac{\sum_t p_{kt}R_{ti}R_{tj}}{\sum_t p_{kt}^2}.
\]

This estimates task modulation of residual cross-products. It is not generally equal to the seed×task coefficient in a target regression. `mppi_seed_slopes()` and the seed-set helpers extract entries or projections from these matrices; they do not supply the missing conventional regression. `seed_ppi_mcca()` implements another, explicitly multivariate construction: association between a seed set and local-data×design components after partialling main effects.

Conventional PPI fits target BOLD using physiological, psychological, interaction, and nuisance terms. In a deconvolved formulation, the interaction is constructed in neural space and reconvolved for the BOLD regression. gPPI includes the condition-specific terms needed for the requested experimental contrasts. These are the relevant reference targets, rather than equality to another implementation of the covariance formula. [SPM PPI manual](https://www.fil.ion.ucl.ac.uk/spm/docs/manual/ppi/ppi/), [McLaren et al., 2012](https://pubmed.ncbi.nlm.nih.gov/22484411/)

**Decisive counterexample:** the seed has within-condition variances 1 and 4, while the target's seed coefficient is fixed at 0.7. The fixture adds target noise orthogonal to the regression design. `lm(target ~ seed * task)` gives an interaction coefficient of `3.77e-17`; multiPPI gives raw covariance modulation **1.05**, or **0.59436** after its normalization. This is a valid distinction between estimands, not a rounding error. The covariance result can be useful, but cannot be advertised as an exact acceleration of conventional PPI.

## Reproduced findings

The priorities below describe proposed repair urgency, not filed tickets. P0 identifies silent scientific-output defects; P1 identifies a required missing capability or a common-path failure. The [consolidated probe outputs](probe-results.txt), [integration output](integration.log), and accompanying scripts contain the numerical witnesses.

| Priority / claim | Evidence and impact | Source / acceptance target |
|---|---|---|
| **P1 — Conventional PPI/gPPI is missing from the inspected public paths** | The counterexample above distinguishes covariance modulation from an interaction coefficient. Neural-domain fitting still computes cross-products; it does not reconvolve a seed×task regressor and fit target BOLD. | `R/mppi_fit.R:330`, `R/mppi_seedset.R:102`, `R/seed_ppi_mcca.R:133`. Add an explicit conventional engine and compare coefficients, standard errors, and contrasts with independent scalar GLMs/GLS. Retain the covariance engine under an explicit name. |
| **P0 — Default matrix permutations do not randomize** | `method="block_flip"` with default `wild="none"` produces all-one weights. All 49 null statistics were identical; global and edgewise p-values were 1. | `R/mppi_inference.R:24`, `R/mppi_inference.R:123`. Test null variation and calibrated behavior, not just matrix shape or RNG metadata. |
| **P0 — Observed and null omnibus statistics use different diagonals** | The observed matrix has zero diagonal; the null matrix retains it. The observed/null discrepancy was **0.1029173**, exactly the diagonal energy. This also affects explicitly requested random weights and phase resampling. | `R/mppi_inference.R:77`, `R/mppi_inference.R:88`. Use exactly the same statistic and diagonal policy for observation and null. |
| **P0 — Seed-CCA permutation null is invariant** | Both inputs receive the same row signs or permutation. Every Gram matrix is preserved. A strong association with root **0.958053** had zero null deviation and p=1. The helper also ignores caller-specified ridge penalties. | `R/seed_ppi_mcca.R:415`. Establish a justified nuisance-aware permutation procedure and use the requested regularization in every draw. |
| **P1 — Public prewhitening dispatch fails** | `mppi(..., prewhiten=TRUE)` errors with `unused argument (deconv = deconv)`. The direct full-rank `mppi_fit_whitened(..., domain="bold", p=1)` call succeeds. This localizes the immediate fault to multiPPI dispatch rather than fmriAR. | `R/mppi_s3.R:24`, `R/mppi_fit_whitened.R:30`. Align argument/domain contracts and exercise the public matrix, design, and event-model routes. |
| **P1/P0 — Standard runwise baselines become rank deficient** | `as_mppi_design(event, baseline)` adds a global intercept alongside two run intercepts. The real two-run fixture had 5 columns and rank 4. Baseline RSS became nonfinite; direct whitening rejected the resulting NA residuals. | `R/mppi_design_interface.R:303`, `R/mppi_fit.R:158`, `R/mppi_fit_whitened.R:49`. Detect an intercept in the column span, use rank-safe fitted/residual operations, and count rank rather than column count for degrees of freedom. |
| **P1 — Stack adapters use obsolete or nonexistent APIs** | A real `fmrireg::matrix_frame()` and `fmri_model()` reach an error for nonexistent `fmrireg::get_data_matrix`. `mppi(frame, ...)` has no method for current `fmri_frame`. Neural seed-CCA calls nonexistent `fmrihrf::deconvolve_series`. | `R/mppi_fit.R:511`, `R/mppi_s3.R:160`, `R/seed_ppi_mcca.R:317`. Integrate the current `fmri_frame`/`collect_assay`/`temporal_schema` contract, coerce the actual tibble design to a numeric matrix, and use a real deconvolution implementation. |
| **P0 — Event-model TR is silently lost** | An actual fmridesign event model reports TR=2; `mppi(event, Y=..., domain="neural")` records and uses TR=1. The design knows the correct TR, but this dispatch path does not pass it into deconvolution. | `R/mppi_s3.R:81`, `R/mppi_fit.R:31`. Propagate sampling metadata through every entry point and check consistency. The corrected event probe is in `event-probe.log`. |
| **P0 — Neural estimation crosses run boundaries** | A perturbation confined to run 2 and orthogonal to the fitted design changed run-1 BOLD residuals by only `9.99e-16`, but changed run-1 neural estimates by **9.4397**. The concatenated deconvolution has no run separation. | `R/mppi_fit.R:196`, `R/mppi_fit.R:248`. Fit each independent temporal block with its own forward operator; test isolation by perturbing a different run. |
| **P0 — MAP deconvolution changes HRF timing** | `h <- h[h != 0]` removes meaningful zero samples. With HRF `c(0,1)`, the output is identical to using HRF `1`, and differs from an independent causal solver by **0.5710**. | `R/deconv.R:39`. Preserve HRF sample positions; compare against a dense causal convolution reference, including leading and interior zeros. |
| **P0 — Packed lagged results lose direction** | The full lag-1 result had asymmetry **0.199289**. Packing and unpacking made it symmetric, with a maximum error of the same magnitude. Triangular packing is valid only for symmetric matrices. | `R/mppi_fit.R:428`, `R/mppi_utils.R:182`. Use full/rectangular storage for directed lags or reject the incompatible option. |
| **P0 — Censoring creates false temporal neighbors** | Omitting one NA row compresses time. A lag-1 calculation used 118 pairs where only 117 were truly adjacent, changing the result by **0.0243833**. Neural deconvolution has the same time-grid problem after row omission. | `R/mppi_utils.R:63`. Preserve acquisition indices and valid contiguous segments through lagging, filtering, whitening, and deconvolution. |
| **P0 — Rank selection destroys negative modes** | The approximation uses left singular vectors on both sides with nonnegative singular values. A strong negative rank-one covariance slope with eigenvalue **−16.2102** was assigned rank 0; rank-1 loss was 126.23 versus 37.26 for zero. | `R/mppi_rank_cv.R:44`. Preserve signed eigenvalues or both SVD factors and compare held-out losses with an independent reference. Random individual-TR splitting also needs reconsideration for serially dependent data. |

Additional reproduced scientific-output problems:

- **Mammen probabilities are reversed:** the implemented two-point distribution has mean 1 and second moment 2, instead of 0 and 1. See `R/mppi_inference.R:27` and `integration.log`.
- **CCA parametric degrees of freedom are wrong for the classical full-root test:** for p=q=2 the code uses 5 rather than 4. The probe returns p=`0.0008974` versus the classical p=`0.0003526`. Default ridge roots, nuisance rank, and partial-root inference require separate justification; changing the degrees-of-freedom expression alone is insufficient.
- **Evidence measures task main effects, not PPI:** a synthetic activation-only fixture produced a task-model ΔAIC of 3709.25. Two different deconvolution HRFs received identical `AIC_total` despite different interaction matrices. The RSS is computed before deconvolution and before any interaction regression. Consequently these AIC values cannot discriminate the HRFs or establish evidence for coupling modulation (`R/mppi_fit.R:165`, `R/mppi_evidence.R:3`, `R/mppi_neural.R:256`).
- **Seed representation changes the analysis:** supplying a raw seed column as a matrix gave statistic 0.18495; specifying that same column by index gave 0.21908. Indexed seeds come from preprocessed data while matrix seeds bypass that pipeline (`R/seed_ppi_mcca.R:112`, `R/seed_ppi_mcca.R:338`). The API needs an explicit raw-versus-preprocessed seed contract.

Further concerns supported by implementation inspection, without empirical error-rate qualification in this audit:

- “Studentization” uses a weighted cross-product that can be negative, rather than a sampling variance, and replaces negative values with another entry's positive minimum. With Rademacher weights this factor is constant across draws and cancels in an individual edge's ranking, so this observation alone does not establish a changed p-value in that specific mode.
- The “Freedman–Lane” path demeans cross-products rather than fitting their complete reduced nuisance model. Residualized time series do not imply residualized products. Resampling blocks ignore run boundaries.
- `mppi_lagged()` shifts the psychological regressor and returns symmetric contemporaneous cross-products, whereas `mppi_fit(lags=...)` pairs signals at different times. These are different lag estimands. `mppi_lag_select()` chooses the observed best lag and then permutes only that lag, without repeating selection in the null.
- Seed-CCA's neural path whitens before deconvolution and builds interactions after whitening; it also leaves the psychological design in its supplied domain. This is not the conventional deconvolve–interact–reconvolve–fit sequence. Whitening and multiplication do not generally commute. An SPM development discussion identifies the same prewhitening-before-deconvolution concern. [SPM issue 68](https://github.com/spm/spm/issues/68)

## What is sound and worth retaining

- The fixed-penalty DCT solver matched an independent dense causal-convolution/DCT ridge solution with maximum error **2.12e-15**, including a nonuniform penalty vector.
- The native batch kernel matched an explicit per-time-point outer-product reference across zero, positive, and negative lags and two unequal run lengths, with maximum error **1.67e-16**.
- BLAS, chunked, and accumulating weighted cross-products agreed within **7.11e-14** on the measured fixture.
- The CCA Cholesky orientation is correct. With centered data and zero ridge, canonical roots matched `stats::cancor()` within **1.11e-16**. It should be protected by a test, not “corrected.”
- **fmriAR is the appropriate whitening component.** In an independent conventional-PPI reference harness, a completed BOLD design was passed to `fmriAR::fit_noise()` and `whiten_apply()`. Design and response whitening matched an explicit run-block AR(1) operator within `4.44e-16`; fitted coefficients matched dense GLS within **5.33e-15**. This verifies the operator for that fixture, not AR-model calibration, all noise structures, or an existing conventional multiPPI API.

## Efficiency

The core is fast for moderate ROI matrices. The final measurement used T=600, V=400, three psychological columns, one warmup and three timed repetitions per case. Allocation profiling was a separate invocation. Timings are local elapsed observations, not a general performance guarantee; the default timing had a 0.197-second outlier. [Benchmark data](benchmark.csv), [script](benchmark.R), [log](benchmark-isolated.log).

| Configuration | Median elapsed | Range | Cumulative R allocations | Uncompressed saved fit |
|---|---:|---:|---:|---:|
| Default normalized BLAS | 0.025 s | 0.024–0.197 s | 64.2 MiB | 45.2 MiB |
| Chunked | 0.031 s | 0.030–0.032 s | 88.0 MiB | 45.2 MiB |
| Packed | 0.024 s | 0.023–0.024 s | 74.2 MiB | 41.5 MiB |
| Covariance with automatic diagnostics | 0.060 s | 0.059–0.064 s | 149.0 MiB | 72.0 MiB |
| Supplied rank-25 basis | 0.007 s | 0.007–0.011 s | 19.5 MiB | 0.606 MiB |
| Neural DCT, 48 components | 0.042 s | 0.040–0.066 s | 66.7 MiB | 47.0 MiB |

Allocation totals are **not peak resident memory**. Serialization repeats some shared R objects and is also not a resident-memory measurement. The basis case estimates a reduced-space quantity, so its speed is not evidence of an exact full-space speedup. The input matrix itself is about 1.83 MiB.

The main scaling problems are explicit in the implementation:

1. Full-space fitting requires dense V×V outputs and roughly O(KTV²) cross-product work. At V=50,000, **one** double-precision matrix needs 20 GB (18.6 GiB). Raw, amplitude, and normalized matrices are retained, along with baseline covariance and optional diagnostics. `packed=TRUE` does not remove those dense backing matrices.
2. `mppi_permute()` allocates `array(0, c(V,V,B))`. At V=1,000 and B=999, that allocation alone is **7.992 GB**. Stream exceedance counts or requested edge statistics instead.
3. “Freedman–Lane” and seed-map permutation paths materialize T×V(V+1)/2 cross-products. For T=600 and V=1,000, one such matrix is **2.4024 GB**. Algebraically moving a row permutation to the predictor permits weighted cross-products without this outcome expansion, once the statistical procedure is correctly specified.
4. Covariance output unconditionally computes precision and variance diagnostics. It solves for the same baseline precision once in `mppi_fit()` and again inside `mppi_to_partial()` for each condition. Cache that factorization and make diagnostics optional.
5. The MAP solver starts sparse but converts its T×T system to dense for Cholesky. Its temporal memory is O(T²) and factorization O(T³). Preserve band structure or prefer the validated DCT formulation when appropriate; this requires an equivalent-output oracle.
6. DCT preparation is repeated for data and psychological columns. Cache a preparation object keyed by run length, TR, HRF, basis dimension, and penalty. For conventional seed PPI, only the seed needs deconvolution; deconvolving all target voxels is unnecessary.

The first optimizations should preserve the estimator: stream permutation counts; avoid redundant matrices and inversions; compute selected seed-to-target outputs directly; reuse design and deconvolution factorizations. Approximate spatial bases should remain explicit options with separate accuracy qualification.

## Recommended implementation boundary

Use the ecosystem components according to their actual responsibilities:

| Component | Responsibility |
|---|---|
| **fmrihrf** | HRF evaluation/convolution and sampling conventions. Its installed API does not supply the deconvolution function currently called by seed-CCA. |
| **fmridesign** | Task main effects, condition identity, durations/amplitudes, run-aware sampling, nuisance design, and contrasts. Obtain neural psychological inputs from event timing or explicit user input rather than calling an already convolved matrix “sticks.” |
| **multiPPI** | Seed extraction contract; run-local seed deconvolution; neural interaction construction; reconvolution; conventional PPI/gPPI orchestration; the separately named covariance and multivariate methods. |
| **fmriAR** | Estimate temporal noise from appropriate model residuals and apply the resulting operator to the **completed design and target data together**, respecting runs and censoring. Expose the relevant noise/pooling controls and retain the fitted plan as provenance. |
| **fmrireg, optional** | Reuse first-level fitting, contrasts, and result infrastructure where its current public interfaces fit the contract. Do not reimplement its dataset layer or require it merely to evaluate an HRF. |

For a conventional seed-to-target analysis, construct

\[
D=[\text{run baselines},\ \text{confounds},\ Hp_1,\ldots,Hp_K,\ s_{BOLD},\ H(p_1\odot\hat s_n),\ldots,H(p_K\odot\hat s_n)].
\]

Here H is run-local convolution with a specified sampling convention and the neural task/seed centering policy is explicit. Fit target BOLD against D. With a common whitening operator, one QR factorization per seed/design can serve all target columns, processed in chunks. With parcel-specific or voxel-specific whitening, reuse is limited to targets sharing the operator. Return condition coefficients, standard errors, degrees of freedom, estimability, and named contrasts. These maps need O(KV) output storage per seed rather than mandatory V² storage.

The first implementation slice should be a matrix-input conventional engine plus independent GLM/GLS oracle tests. Then connect real event/frame objects, repair inference, and optimize the verified computation. Concrete gates:

1. Recover known interaction coefficients and A−B contrasts; the variance-change counterexample must produce zero conventional interaction. Cover correlated conditions and nuisance regressors.
2. Agree with an independent scalar target-by-target fit, including standard errors and contrasts, under both IID noise and a fixed AR whitening operator.
3. Preserve run isolation, acquisition time, HRF zeros, durations, TR, and censored gaps. Exercise real fmridesign/fmri_frame inputs and all public whitening dispatches.
4. Demonstrate nondegenerate null statistics and empirical error control on bounded simulations with serial dependence, nuisance effects, and multiple conditions. Specify and test multiplicity/lag-selection correction separately.
5. Benchmark equivalent outputs at ROI scale and a realistic seed-to-voxel workload, recording peak memory as well as elapsed time. Accept no speedup obtained by silently changing the estimator.

## Test and packaging evidence

The explicit `testthat::test_dir()` run produced **377 passing expectations, one error**, and skip notices. The error is the obsolete `fmridataset::matrix_dataset()` call in `test_design_interface.R`. Enabling the initially skipped PWI files with `NOT_CRAN=true` added **18 passing expectations**. Two dependency-related skips remained: a missing-dependency error test skipped because fmrireg is present, and the seed-whitening test.

The seed-whitening skip is particularly misleading: the test centers its intercept to zero, computes residuals using rank-deficient `qr.coef()`, and labels the ensuing fmriAR rejection as an environment failure (`tests/testthat/test_seed_ppi_mcca_whiten.R:22`). That guard can conceal a local preprocessing defect.

The package has **no `tests/testthat.R` harness**, and `testthat` is absent from Suggests. The ordinary package check therefore does not run the testthat suite. Several simulation tests additionally compute their own local `Mk()` cross-products instead of calling the package's fitted workflow; their passing results cannot qualify that workflow.

`R CMD build --no-build-vignettes --no-manual` succeeded. The normal dependency-enforcing `R CMD check --no-manual --no-build-vignettes` stopped with **one error because limma is unavailable**. A separately labelled diagnostic check with `_R_CHECK_FORCE_SUGGESTS_=false` completed with **0 errors, 5 warnings, and 5 notes**. It confirmed undeclared dependencies, missing/unexported APIs, the invalid prewhitening call, documentation and Makevars issues, and missing vignette configuration. This reduced check is not a complete package qualification. Full logs are retained in [check.log](check.log) and [check-reduced.log](check-reduced.log).

A bounded independent mathematical review corroborated the estimand and inference findings and rejected a suspected CCA Cholesky-orientation issue. Actual probes, rather than reviewer agreement, establish the runtime findings above. This audit does not establish whole-package statistical validity, group-level calibration, neuroimaging volume/searchlight indexing correctness, hosted CI status, or cross-platform readiness. The resting-state, classification, and mechanistic interpretation surfaces received only incidental coverage through the existing suite.

## Reproduction

Run `python3 audit/2026-10-23-ppi/reproduce.py` from this checkout. It copies current source to a new temporary directory, records differences from the assessed source manifest, runs the suite and probes there, and prints the result locations. It does not modify the package implementation or install dependencies. The probe scripts deliberately capture errors and numeric outputs: a script exit of zero does not mean every assessed claim passed. Inspect the named results and acceptance targets above.

Original output is preserved, including the first event-model probe's missing formula-environment binding. That probe was corrected by binding the real `fmridesign::hrf` function and rerun; `event-probe.log` and `probe-results.txt` contain the corrected result. This harness correction did not change package source.
