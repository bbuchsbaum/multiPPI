# multiPPI

Seedless, matrix-form Psychophysiological Interactions (PPI) for fMRI—fast, closed-form,
and integrated with the `bbuchsbaum` stack (`fmrireg`, `fmriAR`, `fmrilss`, `fmridesign`, `neuorim2`, `fmridataset`).

## Key features

- **Matrix-PPI**: closed-form ΔΣ per psychological regressor (no seed loops).
- **fmriAR prewhitening**: accuracy with serial autocorrelation.
- **Lagged / directed** variants with validated lag selection.
- **HRF-adaptive** weights (canonical+derivatives+lags) via small eigenproblem.
- **Neural-domain** fits via automatic DCT deconvolution and canonical HRF defaults.
- **Low-rank** shrinkage with cross-validated rank selection.
- **Partial-correlation** deltas (ΔΘ ≈ −Θ₀ ΔΣ Θ₀).
- **Variance–correlation decomposition** for interpretation.
- **β-mPPI**: trial-domain analogue using `fmrilss` betas.
- **Frequency-specific** ΔΣ with simple bandpass.
- **Group eBayes** (limma) or classical t as fallback.
- **Mechanistic readouts**: gain/precision, routing asymmetry, and template projections.

## Install

```r
# devtools::install_local("multiPPI")  # once you unzip this archive locally
# or build with R CMD build / R CMD INSTALL
```

## Neural-domain workflow

Set `domain = "neural"` to work directly with latent neural activity. The
solver now deconvolves residual BOLD and psychological regressors using a
canonical double-gamma HRF (`mppi_default_hrf()`) when you do not supply one
explicitly. Additional deconvolution controls live in a `deconv` list—override
the HRF shape, DCT dimension (`K`), generalized cross-validation vs fixed ridge,
oversampling, or even pass precomputed neural “sticks”. When the design matrix
contains HRF basis expansions, columns are grouped automatically so the fit
produces one neural regressor per condition.

```r
h <- mppi_default_hrf(tr = 0.8, duration = 32)
fit_neural <- mppi_fit(Y, X, psych_idx = psych_cols,
                       domain = "neural",
                       deconv = list(hrf = h, K = 48))

# Compare against a BOLD-domain fit
fit_bold <- mppi_fit(Y, X, psych_idx = psych_cols, domain = "bold")
neural_diff <- mppi_compare_models(fit_bold, fit_neural)

# Access neural sticks and mechanistic readouts
sticks <- fit_neural$deconv$sticks
gain   <- mppi_gain(fit_neural, k = 1L)
route  <- mppi_routing_index(list("lag0" = fit_neural), k = 1L,
                             hierarchy = rep(1, ncol(Y)))
```

For bespoke workflows, `mppi_psych_neural_from_X()` generates neural sticks
from any design matrix, `mppi_neural_from_fit()` converts existing BOLD fits,
and `mppi_hrf_ensemble()` blends fits from multiple HRFs. These helpers are now
exported with documentation alongside the core API.

## Mechanistic readouts

multiPPI now ships lightweight “DCM-lite” diagnostics that operate entirely in
the fitted subspace (ROI or basis):

- `mppi_gain()` — converts task slopes into first-order precision updates and
  reports component gain with optional ROI back-projection.
- `mppi_routing_index()` — combines lagged fits with a cortical hierarchy to
  score feedforward vs feedback routing.
- `mppi_project_templates()` — projects interaction matrices onto hypotheses
  (within/between, gradients, bespoke mechanisms) via simple inner products.
- `mppi_mechanism_report()` — bundles the above into a single, ready-to-plot
  summary list.

Template helpers (`mppi_templates_within_between()`, `mppi_templates_gradient()`)
provide quick starting points, but any r×r (or ROI×ROI) matrices work. Because
everything happens in the reduced space, the workflow scales effortlessly to
high-resolution data.
