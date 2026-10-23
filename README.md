# multiPPI

multiPPI provides conventional seed-based PPI/gPPI regression and a separate
family of task-modulated covariance estimators for fMRI. The two families have
different estimands: conventional PPI fits seed-by-task interaction regressors;
`mppi_fit()` estimates task-modulated covariance slopes.

## Installation

```r
# install.packages("remotes")
remotes::install_github("bbuchsbaum/multiPPI")
```

The default HRF (`fmrihrf`) and temporal noise model (`fmriAR`) are installed
from GitHub automatically. Optional adapters for `fmridesign`, `fmridataset`
and `fmrireg` objects need `dependencies = TRUE`.

## Conventional PPI

```r
library(multiPPI)
set.seed(1)
n <- 100
seed <- rnorm(n)                     # raw BOLD, full acquisition grid
psych <- cbind(A = rep(c(0, 1), 50), B = rep(c(0, 0, 1, 1), 25))
design <- ppi_design(seed, psych = psych, TR = 2, run_effects = "shared")
Y <- matrix(rnorm(n * 2), nrow = n) # target BOLD, same scan rows
fit <- ppi_fit(Y, design, noise = "iid")
ppi_contrast(fit, c(A = 1, B = -1))$estimate
```

A matrix supplied as `psych` represents scan-grid neural psychological inputs,
not an already-convolved task design. Event timings or explicit microtime-grid
inputs are also supported. By default, the raw seed input is adjusted for run
baselines and nuisance effects conditional on task main effects before
deconvolution. See the [conventional PPI vignette](vignettes/vignette_09_conventional_ppi.Rmd) and `?ppi_design`
for timing, centering, and convolution conventions.

The default `fmriAR` fit estimates a shared temporal noise operator. Shared
pooling assumes similar target noise; it can be inappropriate for heterogeneous
targets. Use `noise_control = list(pooling = "parcel", parcels = labels)` to
fit a separate operator per target group, or one label per target for targetwise
estimation. Estimated-noise inference is approximate and conditional on the
seed and design. Censored fits currently support IID errors only; AR with gaps is
rejected. These results do not establish causal or latent neural coupling.

## Covariance analysis and validation

The `mppi_*` functions expose the distinct covariance workflow, lagged slopes,
and derived quantities. Unsupported legacy resampling options fail closed;
legacy task-main-effect evidence should not be interpreted as interaction
model evidence. See the [validation guide](inst/validation/README.md) and
[implementation evidence](audit/2026-10-23-ppi/IMPLEMENTATION.md) for external
comparisons, simulation scope and limitations. Passing those checks does not establish
end-to-end neural recovery, full inferential calibration, or release readiness.
