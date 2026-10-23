# Comprehensive multiPPI repair plan

**Status: proposed implementation plan, 23 October 2026. No repairs have been implemented by this planning pass.**

Implementation subsequently proceeded; see the [implementation and qualification report](IMPLEMENTATION.md) for delivered changes, executed checks and remaining restrictions. The plan below is retained as the original acceptance target, not as a claim that every proposed experiment has run.

The goal is a correct, usable conventional PPI/gPPI workflow, backed by fmrihrf/fmridesign and fmriAR, alongside an explicitly identified covariance-modulation workflow. Completion also requires repairing the existing silent scientific-output defects, qualifying supported inference, and demonstrating efficient execution on equivalent outputs.

This plan follows the [source-bound audit](REPORT.md), its [reproductions](probes.R), [integration witnesses](integration.R), and [source manifest](source-manifest.json). The baseline is HEAD `9d199cd53d484317c94a3cbaa6dfee972cb47f99` **plus the audited working-tree changes**. Preserve those changes and the original evidence. A checkout of HEAD alone is not the assessed baseline.

## 1. Decisions that guide the implementation

| Area | Proposed decision |
|---|---|
| Conventional estimator | Add a separate `ppi_design()` → `ppi_fit()` → `ppi_contrast()` workflow. These are proposed names, not existing functions. Support classic contrast-coded PPI and jointly fitted condition-specific gPPI. |
| Multiple runs | Default to run-specific first-level terms and explicit within-subject contrast weights. Offer shared coefficients as a separate model choice. Define equal-run versus precision-weighted aggregation and handle absent/nonestimable conditions explicitly; runs are not independent subjects. |
| Existing estimator | Preserve `mppi()`/`mppi_fit()` as covariance modulation. Add explicit estimand/domain/storage metadata and accurate documentation; never silently change these functions to conventional regression. Correct scientifically wrong outputs with migration notes and a result-schema version. |
| HRF and design | Use fmrihrf for supported HRF/sampling operations and fmridesign for event, nuisance, condition, and contrast information. multiPPI owns the seed deconvolution and interaction orchestration absent from those public interfaces. |
| Whitening | Use **fmriAR**, with an explicit noise specification and saved fitted plan. The new conventional API defaults to fmriAR prewhitening once its qualification gates pass; `noise="iid"` is an explicit supported alternative. Missing fmriAR must produce an informative error, never an unnoticed IID fallback. |
| fmrireg | Keep it an optional adapter/fitting integration. Reuse its public infrastructure where it satisfies the contracts; do not let its availability block the matrix reference implementation or recreate its dataset layer. |
| Inference | Distinguish conventional regression, covariance slopes, and seed-CCA hypotheses. A resampling option is available for scientific use only after its null, nuisance treatment, and dependence assumptions have been specified and tested. |
| Efficiency | First preserve the mathematical result: chunk targets, reuse factorizations, stream null summaries, and avoid redundant dense matrices. Spatial low-rank approximations remain explicit and receive separate accuracy evidence. |
| Compatibility | Repair public dispatch and retain correct numerical kernels. Deprecate ambiguous arguments/labels with migration examples. Unsupported combinations fail clearly rather than return plausible p-values or silently alter the model. |

The minimal conventional model is

\[
D=[N,\ Hp_1,\ldots,Hp_K,\ s_{BOLD},\ H(p_1\odot\widehat{s}_n),\ldots,H(p_K\odot\widehat{s}_n)],\qquad Y=D\beta+\epsilon.
\]

Here `N` includes run baselines and confounds; `p` is defined on the neural event grid; `H` includes run-local convolution and sampling to scans. The physiological BOLD series, its adjustment, neural deconvolution, centering, and scaling must all be recorded. Task and physiological main effects are part of the model. The interaction-construction rationale follows the [SPM PPI manual](https://www.fil.ion.ucl.ac.uk/spm/docs/manual/ppi/ppi/); gPPI condition modelling follows [McLaren et al.](https://pubmed.ncbi.nlm.nih.gov/22484411/). This is a defined PPI implementation, not a promise of numerical identity to SPM's particular deconvolution prior.

Build the interaction columns **before final regression whitening**. Fit noise from appropriate full-model residuals, then apply the same fitted temporal operator to the complete design and each corresponding target. Noise-aware seed deconvolution may be a separate, explicitly defined inverse problem; whitening a seed and multiplying it by a task vector is not a substitute.

## 2. Work packages and dependencies

| ID | Deliverable | Depends on | Main implementation scope |
|---|---|---|---|
| F01 | Executing regression harness and explicit model contracts | Baseline audit | `DESCRIPTION`, `tests/testthat.R`, `tests/testthat/`, API contract documentation |
| F02 | Unified sampling, design, seed, and current-stack adapters | F01 | `R/mppi_design_interface.R`, `R/mppi_s3.R`, `R/mppi_utils.R`, adapter parts of `R/mppi_fit.R` |
| F03 | Correct run-local deconvolution and PPI interaction construction | F02 | `R/deconv.R`, `R/mppi_deconv_dct.R`, `R/mppi_neural.R`, `src/deconv_*.cpp`, proposed `R/ppi_design.R` |
| F04 | Conventional PPI/gPPI fitting, contrasts, and fmriAR integration | F02–F03 | Proposed `R/ppi_fit.R`, `R/ppi_contrast.R`, shared whitening adapter; `R/mppi_fit_whitened.R` |
| F05 | Correct legacy covariance, lags, selection, and derived quantities | F02; F03 for neural paths | `R/mppi_fit.R`, `R/mppi_lag.R`, `R/mppi_rank_cv.R`, `R/mppi_partial.R`, `R/mppi_decompose.R`, evidence/HRF helpers |
| F06 | Estimator-specific inference with scientific qualification | F04–F05; F07's CCA contract | `R/mppi_inference.R`, `R/mppi_seedset.R`, inference in `R/seed_ppi_mcca.R`, beta/group inference |
| F07 | Seed-CCA, imaging, group, and remaining public-API coverage | F02; F03/F04 as applicable | `R/seed_ppi_mcca.R`, group/beta/band/rest/instant/basis/classification helpers and callers |
| F08 | Exact memory/compute improvements and scale measurements | Correctness gate for each affected path | Fit/deconvolution kernels, output storage, resampling, benchmark scripts |
| F09 | Documentation, packaging, CI, and final qualification | All relevant gates above | `DESCRIPTION`, `NAMESPACE`, `man/`, `vignettes/`, `src/Makevars*`, package exclusions, CI |

Work in reviewable slices. The critical path is **F01 → F02 → F03 → F04**, giving an end-to-end conventional workflow early. F05 and the deterministic parts of F07 can follow independently once the shared contracts stabilize. F06 must not be reduced to changing random-number defaults. Optimize only computations already protected by independent references. Packaging work begins in F01 and finishes in F09.

### F01 — Make failures visible and lock down the contracts

1. Snapshot the actual dirty source and dependency versions before implementation. Use a source copy or explicitly preserved patch set if isolating work; do not discard untracked seed-CCA work or copy only committed HEAD.
2. Add the standard `tests/testthat.R` harness and declare testthat. Replace the seed-whitening test's broad error-to-skip guard with checks for genuinely unavailable dependencies. A local rank/preprocessing error must fail.
3. Convert the audit's deterministic counterexamples into regression expectations, attached to their corresponding repair slices. Keep the original audit outputs unchanged. Test package entry points, not local replicas of the implementation.
4. Specify the public model contract: time units, run IDs, original acquisition indices, event durations/amplitudes, neural versus BOLD psychological inputs, seed source/adjustment/scaling, noise model, censor policy, condition identities, contrast weights, rank policy, and result units.
5. Define explicit unsupported combinations. An already convolved task matrix alone cannot identify the original event-level neural task. Require event timing or a separately supplied neural psychological input for deconvolved PPI.

**Gate:** a built package executes the test suite; known defects are visible in a baseline failure record; no implementation failure is hidden as an environmental skip. Record the real optional-dependency matrix, including a working limma environment for group checks.

### F02 — Preserve time and integrate the actual package APIs

1. Normalize matrix, `event_model`, `mppi_design`, current `fmri_frame`, and supported fmrireg model inputs into one validated internal representation. Use current public `collect_assay`, `temporal_schema`, and sampling-frame interfaces; coerce returned design tibbles to numeric matrices deliberately. Replace nonexistent accessors and outdated class dispatch.
2. Carry TR and within-run sampling through every method. Reject conflicting metadata. Handle supported run-specific sampling explicitly; reject irregular sampling if the chosen HRF/noise implementation cannot represent it. Never substitute TR=1 for known sampling information.
3. Detect an intercept in the **column span**, rather than only a literal column of ones. Use rank-revealing QR/SVD for fitted values and residuals. Retain column identities/alias information and test contrast estimability; do not silently discard aliased scientific terms. Count rank for residual degrees of freedom.
4. Use one seed extraction/preprocessing path for indexed, externally supplied, and imaging-derived seeds. Default external series to raw BOLD; any preadjusted input must say what adjustment occurred. Preserve task effects when the documented adjustment removes nuisance. Define ROI mean/eigenvariate weighting and deterministic orientation; mask seed-overlap targets or mark them as ineligible for inferential claims.
5. Preserve the full acquisition grid and an observation mask. Censoring must not compress time. Distinguish a missing observation, an outlier handled by a spike regressor, and a true run boundary; these are different operations. Define a common temporal mask initially, with explicit rejection or separate handling of target-specific missingness.
6. Prefer drift/high-pass regressors in the complete model for the initial conventional path. If a separate temporal filter is supported, define its order, apply it consistently to design and response, and account for its rank. Prevent accidental double filtering.

**Gate:** matrix/frame/event/model representations of the same data agree; a TR=2 event model stays TR=2; two run intercepts do not acquire a redundant global intercept; equivalent raw seed inputs agree; column reordering with named contrasts preserves answers. Censored lag pairs retain original sample distances.

### F03 — Construct the right interaction on the right temporal grid

1. Preserve leading/interior/trailing HRF zeros and physical sample spacing. Specify convolution normalization, microtime/scan offset, event discretization, duration handling, downsampling, padding, and run-edge treatment.
2. Deconvolve each seed independently by run. Keep the validated fixed-penalty DCT solver and its dense reference. Repair the MAP solver against an independent causal convolution/penalty system before changing its factorization.
3. With missing seed observations, use a masked observation operator on the original grid. Do not insert observed zero-valued BOLD or restart the HRF at a censored scan. A missing scan does not erase the preceding neural response. Report unidentified or poorly conditioned cases rather than inventing information.
4. Define seed/task centering and scaling by run, including their effects on interpretation. Default to centered neural psychological inputs with means defined on the declared run grid. Build condition-specific neural interactions and reconvolve them using the declared HRF/sampling operator. Distinguish classic task-contrast PPI from a gPPI model containing the condition-specific interactions jointly.
5. Preserve condition/basis identities. Include required task main effects and an explicit reference-condition policy. Detect exhausted condition codings and other rank deficiencies; return only estimable coefficients/contrasts. Do not infer condition grouping from fragile column-name regular expressions.

Task centering is consequential when reconvolving the estimated neural seed does not recover the measured physiological regressor exactly: a coding shift can add a component the measured-seed main effect cannot absorb. Include an explicitly qualified reconvolved-seed covariate only if supporting a construction that needs it; do not silently add it to every model. Test activation-only and constant-coupling nulls with imperfect deconvolution. [Di, Reynolds and Biswal, 2017](https://pubmed.ncbi.nlm.nih.gov/28105655/)

**Gate:** dense forward/inverse references agree at declared tolerances; delayed HRFs remain delayed; changing run 2 does not affect run 1; censored seed handling agrees with a masked dense solver and is invariant to the discarded values; HRF tails persist across single/adjacent censor gaps; short runs, unequal runs, nonunit TR, overlapping/duration events, and correlated conditions are covered. Fixed-deconvolution numerical correctness and recovery of the unobserved true neural seed are separate claims.

### F04 — Fit conventional PPI/gPPI with fmriAR

1. Implement an independently checkable matrix engine first. Solve the complete design with rank-revealing QR; batch/chunk targets sharing the design. Return interaction coefficients, residual variance, residual degrees of freedom, estimability, standard errors, named t/F contrasts, and sufficient model/noise provenance. Make retention of full residuals/fitted values optional.
2. Build a narrow fmriAR adapter around its real noise-estimation and application interfaces. Estimate noise from full-model residuals with the corresponding full design; pass run/censor information explicitly. Verify design-residual orthogonality without mistaking it for proof that the correct model was supplied.
3. Expose AR order/selection, pooling, and supported design-bias correction controls. Qualify pooled autocovariances for shared estimation rather than averaging spatial time series by default. Keep an adequate covariance-tail budget separate from AR order when using design correction. Reuse one factorization only for targets sharing the same temporal operator; target/parcel-specific operators require corresponding whitened designs.
4. Resolve censoring **before** calling the implementation exact GLS. The separately inspected fmriAR source restarts filtering at gaps; confirm behavior against the exact supported installed/source revision. A filter restart does not by itself reproduce the marginal covariance of a continuous AR process across missing frames. For the default exclusion policy, require observed-data covariance/whitening with original time distances; add a bounded fmriAR capability if the public API lacks it. A restart/independent-segment model may be offered separately, explicitly labelled and calibrated. Remove censored rows from the final fit; retention of rows by an application helper is not their exclusion. Exclude ARMA-with-gaps inference until its estimation and application paths both meet this contract.
5. Define whether estimated noise parameters are frozen or refitted in each inferential procedure. Separate exact fixed-covariance GLS checks from approximate estimated-noise t/F calibration. Do not claim exact finite-sample inference solely because whitening matches an operator.
6. Repair existing `mppi(..., prewhiten=TRUE)` dispatch, `deconv`/domain argument mismatches, and rank-unsafe residual creation using the same plumbing where mathematically applicable. The covariance method still needs its own explicit definition of whether it estimates BOLD, neural, or innovation covariance; sharing the adapter does not make it conventional PPI.
7. Attach event/frame adapters and optional fmrireg fitting/contrast integration only after scalar GLM/GLS equivalence passes. Keep the reference route available for regression testing.

**Gate:** beta, standard errors, t/F contrasts, and degrees of freedom match target-by-target `lm` and an independently constructed dense GLS reference for fixed covariance. Cover two unequal runs, nuisance correlation, nonunit TR, aliased columns, absent conditions, and censored observations. The variance-change counterexample gives zero conventional interaction while retaining the distinct expected covariance slope. Exercise every supported public whitening entry point and actual installed dependency objects.

### F05 — Repair the covariance workflow and derived results

1. Preserve the validated weighted-cross-product kernels and test their public callers. Separate raw covariance quantities from normalized/display versions. Do not zero diagonals before computations that need variances or precision normalization; display/storage diagonal policy must not change the estimator silently.
2. Give psychological-regressor shifts and signal-to-signal lags distinct names and metadata. Use valid pairs at their original within-run sample distances. Keep directed lag outputs full/rectangular or reject symmetric packing; round trips must preserve asymmetry.
3. Preserve signed eigenvalues for symmetric rank approximations, or both SVD factors for general matrices. Use run/block held-out validation appropriate to dependent samples. Learn all adaptive transformations within training folds and define cross-boundary lag exclusion.
4. Rename current pre-interaction RSS/AIC as task-main-effect evidence. Conventional PPI evidence must compare models for the same observed BOLD response, with an appropriate covariance likelihood, parameter accounting, and ML/REML restrictions. Do not rank deconvolution HRFs using identical baseline task RSS. Where a valid complexity/likelihood calculation is unavailable, use properly nested predictive validation rather than an invented information criterion.
5. Check the full derivative from covariance to precision and normalized partial correlation, including diagonal-normalizer derivatives. Verify variance-versus-correlation decomposition against independently constructed pure-variance and pure-correlation changes. These two helper concerns require reproductions before being labelled confirmed bugs.
6. Make HRF adaptation, rank selection, basis reconstruction, and ensemble weights explicit about what is fitted and which data choose it. Carry selection into subsequent uncertainty calculations or restrict the result to descriptive/predictive use.

**Gate:** audit lag/rank/evidence counterexamples are resolved; negative rank-one structure is recovered; packed round trips are exact where allowed; finite-difference derivative checks agree with analytic quantities; zero-lag symmetry, signed-lag identities, and basis/full-space equivalence hold under their documented assumptions.

### F06 — Rebuild inference around stated hypotheses

Write a short method specification before coding each supported inferential path: null hypothesis, sampling unit, nuisance model, dependence model, statistic, transformations, treatment of fitted parameters, and multiplicity family. A variable null distribution is necessary but does not establish validity.

1. **Conventional PPI:** qualify fixed-design analytic GLM/GLS inference first. Treat a fixed, prespecified seed and task conditionally; do not rerun seed-only operations in every target resample unnecessarily. If seed/ROI/HRF/rank selection uses the tested response, perform selection within resampling/held-out data or disallow the corresponding inferential claim. State the conditional-mean/exogeneity assumptions: conditional inference does not correct seed measurement error or shared seed–target artifacts, or establish latent neural/causal coupling.
2. **Covariance modulation:** derive the reduced model for the cross-product outcome actually being tested. Residualizing the original time series does not automatically remove nuisance from their products. Repair or replace the nominal Freedman–Lane path accordingly. Specify heteroscedasticity/studentization; do not turn negative cross-products into variances by clipping or borrowing another edge's value.
3. **Resampling mechanics:** correct the all-one default and Mammen probabilities; use identical observation/null statistics, diagonal policy, and tuning. Respect true run/exchangeability boundaries and censor semantics. Ordinary time-point permutations, arbitrary blocks, and phase surrogates require their own assumptions; a run label alone does not make serial observations exchangeable.
4. **Seed-CCA:** stop applying the same row transformation to both inputs. Preserve caller ridge settings. Protect the correct Cholesky orientation. Restrict classical parametric formulas to justified unregularized, rank/df-valid cases; do not apply them to ridge roots by analogy. Nuisance-adjusted, higher-root, regularized, and serially dependent cases need an explicit suitable procedure. In this construction the nuisance block and interaction can both depend on local target data, so a generic partial-CCA recipe is not automatically valid. Deactivating an invalid p-value is an interim repair, not evidence of a replacement method's validity.
5. **Selection and multiplicity:** repeat the relevant lag/HRF/rank search in the null or use independent selection. For selected-lag claims use the maximum over the same lag family. Define correction over targets/edges, conditions, seeds, roots, and searchlight centers; apply identical draws across the tested family when needed for max statistics. Return the tested family and correction method.
6. **Implementation checks:** use correct Monte Carlo counting with observation inclusion, ties, tails, and finite resolution; detect a degenerate randomization space. Fix RNG streams so chunking/parallel execution does not silently change the analysis. Store null draws only on request; default to streamed counts and required maxima.

Nuisance-aware GLM resampling and CCA resampling have different requirements. The [GLM permutation reference](https://pubmed.ncbi.nlm.nih.gov/24530839/) does not by itself validate arbitrary first-level serial fMRI resampling. The [CCA permutation reference](https://arxiv.org/abs/2002.10046) specifically addresses failures after nuisance residualization and at later roots; extending it to this package's temporal/ridge setting remains a qualification task.

**Gate:** tiny enumerated reference problems agree; no advertised randomized method is invariant by construction; the bounded simulation campaign below passes for each supported inferential mode. Unsupported combinations give precise errors. Restoring every old inference option is not a prerequisite for offering a valid conventional workflow, but unresolved options remain explicitly open in comprehensive completion reporting.

### F07 — Qualify the rest of the exposed analysis surface

Create one coverage table from the actual exports: estimator/units, input contracts, independently checked behavior, inference status, and supported combinations. Group related functions rather than creating a test for every trivial wrapper.

| Surface | Required substantive checks |
|---|---|
| Seed-CCA/searchlights | Indexed/external seed equivalence; correct nuisance rank and regularization; `cancor` agreement in the classical case; no hidden preprocessing asymmetry; target/seed overlap policy; neuroim2 public API and spatial/time axis mapping. Use a tiny image with unique voxel/time values and known neighborhoods. |
| Group/limma/seed-set | Real multi-subject fits with limma installed; subject design dimensions; named condition/edge/feature alignment rather than dimensions alone; missing subjects/contrasts; rank and df; directed maps; subject as sampling unit. Distinguish fixed-effect run aggregation from between-subject inference. |
| Basis/axes/templates | Reconstruction and projection identities, sign/orientation conventions, rank truncation, and explicit approximate-output semantics. Do not average independently oriented latent axes across subjects without a defined alignment. |
| Trial/beta/band/fusion | Correct trial/condition identity, temporal filtering boundaries and edge behavior, beta-estimation uncertainty assumptions, and independent direct references for supported covariance quantities. |
| Rest/instant/state/classification/mechanistic helpers | End-to-end calls using package outputs, correct data orientation and labels, run/subject-held-out evaluation, all learned preprocessing inside folds, and an explicit boundary between descriptive indices, prediction, and inferential or causal claims. |

The audit did not establish new defects in every row. This is a bounded closure pass to avoid declaring the whole package correct from a new conventional engine alone. Reproduce suspected failures before repairing them; preserve behavior that passes meaningful independent checks. Include zero variance, insufficient observations, nonfinite inputs, unequal run lengths, empty neighborhoods, and absent conditions where relevant.

**Gate:** every exported analysis family has a documented supported contract and a meaningful executed check. Inferential features either pass their own F06 gate or have explicitly restricted availability and a recorded unresolved method requirement. Geometry and group-level qualification are reported separately from matrix correctness.

### F08 — Remove avoidable work and measure realistic scaling

1. For conventional maps, deconvolve only the seeds. Reuse one QR per completed design and shared noise operator, process target chunks, and retain requested coefficients/contrasts. Do not build V×V matrices for seed-to-target results. Separate per-seed models from an explicitly requested jointly conditioned multi-seed model.
2. Cache deconvolution preparations keyed by all relevant temporal, HRF, basis, penalty, observation-mask, and boundary settings. Cache the baseline precision factorization once where mathematically shared. Make expensive diagnostics and large intermediates opt-in.
3. Stream permutation exceedance counts, omnibus statistics, and maximum-statistic summaries. Tile selected edges when needed. Avoid both V×V×B null arrays and T×V(V+1)/2 product matrices; prove algebraic equivalence against a small materialized reference after fixing the statistical procedure.
4. Preserve sparse/banded structure in MAP solving if it is equivalent to the specified penalty system. Include factorization preparation and solve costs separately. Keep the validated DCT option; do not present it as identical to a different MAP prior.
5. Design storage around requested outputs. Symmetric packing saves memory only if dense duplicates are not retained. Directed results need directed storage. Estimate unavoidable output size and reject infeasible requests before allocation; require explicit large-output choices.
6. Benchmark serially in isolated processes, recording source/dependency identity, thread settings, cold/setup versus repeated solve time, peak resident memory, cumulative R allocations, and output size. Compare equal estimands, precision, and requested outputs.

**Benchmark matrix:** retain the audited T=600/V=400/K=3 ROI case; add conventional seed maps at T=600 and 1,200, V=1,000 and 50,000, K=3, and one/ten independently fitted seeds. Measure shared and grouped whitening separately. Exercise B=999 inference on a small correctness fixture and a feasible larger streaming case; do not make the largest Cartesian combination mandatory if it adds no scaling evidence.

**Gate:** optimized results meet fixed numerical tolerances against the reference; chunk size changes do not change results beyond expected roundoff; peak memory follows input + chunk + requested-output size for seed maps; null storage does not grow as B×V². Record measured curves and any regressions. No hardware-independent speedup claim or numeric memory promise is made before measurement. Low-rank approximations must show their own accuracy/performance tradeoff.

### F09 — Ship an honest, executable package

1. Declare dependencies actually used, including fmridesign/fmridataset and conditional packages where applicable. Choose Imports versus Suggests consistently with runtime behavior and installation sources. Test both complete and deliberately minimal dependency environments; use supported public APIs and documented version constraints.
2. Fix license metadata, non-ASCII source/check issues, undefined globals, Rd mismatches, seed argument documentation, portable Makevars, and vignette configuration. Exclude local tooling, temporary rendered files, audit working artifacts, and unrelated material from release builds without deleting the user's files.
3. Supply runnable examples for classic PPI, gPPI A−B contrasts, two runs, nuisance regressors, fmriAR whitening, censoring, a seed-to-voxel map, and the separate covariance method. Show coefficient units and supported inference assumptions. Explain migration from the old seed-slope naming, ambiguous lag names, and misleading evidence fields.
4. Run the standard test suite and dependency-enforcing `R CMD build`/`R CMD check`, with vignettes in the full qualification environment. Ensure the check log actually records test execution. Test supported operating systems/compilers with executed CI jobs; a configured or blocked job is not a pass.
5. Save numerical, integration, simulation, and benchmark receipts bound to the final source and dependency revisions. Regenerate reference documentation/site content only after source documentation and examples settle. Publishing, merging, or deploying remains a separate authorized action.

**Gate:** no unexplained errors/warnings in the full required package check, notes individually resolved or justified, runnable installed-package examples, all required dependencies tested, and explicit platform/scientific/performance qualification boundaries.

## 3. Verification tiers and scientific acceptance

| Tier | Evidence | When to run |
|---|---|---|
| Deterministic regression | Audit counterexamples; direct matrix/GLM/GLS references; finite differences; run isolation; rank/estimability; storage identities | With each relevant repair, in ordinary tests |
| Real-stack integration | Actual fmrihrf/fmridesign/fmri_frame/fmrireg/fmriAR objects; matrix-versus-adapter equivalence; real imaging and limma fixtures | With adapter changes and supported dependency revisions |
| Statistical qualification | Null rejection, confidence-interval coverage, bias/power, covariance estimation, selection/multiplicity behavior | Explicit validation scripts outside fast package tests; rerun after statistical changes |
| Performance qualification | Equal-output elapsed/setup time, peak memory and scaling | After correctness; rerun after affected optimizations |
| Packaging/platform | Built and installed package tests, complete check, executable vignettes, executed CI matrix | Before declaring the repair complete |

The statistical campaign must distinguish two data-generating levels: **fixed-design regression** with known interaction columns, and **end-to-end neural→BOLD simulation** with seed deconvolution. The former tests the solver/inference; the latter reveals deconvolution error and HRF sensitivity. Passing one does not imply the other.

Predeclare a compact scenario set: IID; known AR(1); estimated AR under correct and misspecified order; unequal independent runs; correlated task/nuisance effects and condition-dependent seed variance; censor gaps; HRF mismatch/imperfect deconvolution; and a separate lag/HRF selection and multiplicity case. Include short runs, high-leverage designs, and heterogeneous target noise in the estimated-AR tier. Re-estimate the declared noise model independently in each replicate. Add CCA/group scenarios only for methods being qualified. Cross scenarios only where an interaction is substantively necessary.

Use fixed seeds and a small pilot to establish runtime, then freeze the final simulation design and acceptance criteria before seeing qualification outcomes. A practical initial release criterion is 2,000 independent null replicates per primary scenario, nominal alpha 0.05, and a binomial upper bound no greater than 0.075, using simultaneous one-sided 95% coverage across the prespecified primary scenarios. This is an explicit tolerated-inflation screen, not proof of exact calibration. Report the full interval and raw counts, interval coverage, and matched nonzero-effect results so a conservative/broken procedure cannot pass simply by returning p=1. Use B=999 where resampling is under test, enumerate tiny cases exactly, and expand inconclusive runs under a prespecified rule rather than tuning thresholds after failure.

Document narrower supported settings when a method fails. Residual whiteness diagnostics support diagnosis but do not replace coefficient/coverage/error-rate evidence. A correct fixed AR operator, a calibrated estimated-noise procedure, and robustness to HRF/noise misspecification are distinct claims.

## 4. Audit finding closure map

| Audit finding or open coverage area | Required work/gate |
|---|---|
| Missing conventional PPI; covariance/seed-slope estimand confusion | F01, F03–F04; variance-change counterexample; F09 migration |
| Broken whitening dispatch and rank-deficient baseline residuals | F02, F04; all public paths and rank/df oracles |
| Obsolete frame/model APIs; nonexistent seed deconvolution API | F02–F03, F07, F09; actual dependency integration |
| TR loss, cross-run deconvolution, deleted HRF zeros | F02–F03; temporal perturbation and dense convolution references |
| Compressed time after censoring; AR gap handling | F02–F04; original-distance lag and observed-data GLS references |
| Same seed supplied two ways gives different results | F02, F07; one preprocessing contract |
| All-one randomization; wrong Mammen probabilities; diagonal mismatch | F06; exact moments, null mechanics and calibration |
| Unsupported studentization/reduced product model/run resampling | F06; method derivation, reference statistics and error control |
| CCA invariant null, ignored ridge, invalid parametric df/application | F06–F07; classical oracle plus qualified replacement/restriction |
| Packed directed lags; incompatible lag meanings | F05; asymmetric round trips and distinct documented estimands |
| Lag selection outside null; adaptive HRF/rank/basis selection | F05–F07; nested selection and multiplicity qualification |
| Lost negative modes; temporally leaking validation folds | F05; signed-mode and held-out-block references |
| AIC compares baseline task fits rather than interactions/HRFs | F05; same-response model comparison or nested prediction |
| Partial-correlation/variance-decomposition concerns (not yet reproduced) | F05; finite differences and pure-effect fixtures before repair |
| V² outputs, V²B nulls, product expansion, repeated solves, dense MAP | F08; equal-output memory/time evidence |
| Group, geometry, beta/band/rest/classification incidental coverage | F07; export-family coverage matrix and relevant F06 gates |
| Missing test harness, concealed skip, missing deps, package warnings | F01, F09; actual installed-package test/check execution |

## 5. First implementation slice and completion boundary

**Start with F01–F02 and the matrix subset of F03–F04:** establish the harness, normalize time/design/seed metadata, repair rank handling and run-local interaction construction, then deliver a matrix PPI/gPPI fit with named contrasts and independent IID/fixed-fmriAR GLS witnesses. Keep the first change reviewable; adapters, estimated-noise calibration, and arbitrary censoring close their own gates immediately afterward. A matrix-only milestone is useful, but is not the comprehensive completion claim.

The highest-risk dependencies are exact censor/noise behavior in fmriAR, calibrated first-level resampling, and regularized temporal CCA inference. Resolve these through bounded mathematical specifications and runnable references. If fmriAR needs a change, isolate the provider contract and cross-package fixture; do not add a second private whitening implementation to multiPPI as an unnoticed fallback.

The comprehensive fix is complete when the supported conventional workflows work through real package inputs; all reproduced silent-output defects are corrected or their invalid operations explicitly unavailable with a recorded replacement decision; the public analysis families have their stated coverage; statistical claims pass their own gates; equivalent-output performance is measured; and the built package passes its required checks. Any remaining restricted research method, unexecuted platform, or unmet inference gate must remain visible in the completion report. No implementation, scientific qualification, CI execution, or performance improvement is claimed by this plan itself.
