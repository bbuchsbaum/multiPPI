
# multiPPI Neural-Domain Add-on

This file provides **drop-in functions** to compute **neural-domain multiPPI** from any existing
`multiPPI` first-level fit (basis or ROI). It uses the fast Rcpp/Armadillo deconvolution in
**mppiDeconv**.

## Install `mppiDeconv`

```r
install.packages("mppiDeconv_0.1.0.zip", repos = NULL, type = "source")
library(mppiDeconv)
```

## Use with a basis-space fit

```r
# 1) run your existing bold-domain basis fit
fit_bold <- mppi(ev, Y, base = base, runs = runs, prewhiten = TRUE,
                 basis = as_mppi_basis(V_shared), scale = "corr")

# 2) build neural-domain from that fit:
source("neural_addon.R")
h <- fmrihrf::spm_hrf(TR)  # or another kernel
neural <- mppi_neural_from_fit(fit_bold, X = design$X, psych_idx = design$psych_idx, h = h)

# 3) compare bold vs neural (Δ-gap and optional AIC):
cmp <- mppi_compare_models(fit_bold, neural, resid_bold = fit_bold$Z, resid_neural = neural$U, pk = neural$pk)
cmp$delta_gap
```

## Use with ROI-space fit

Same, but `fit_bold` is your ROI-space `mppi()` result. The function returns `Delta` (V×V) lists.

## HRF ensembles

```r
neural_canonical <- mppi_neural_from_fit(fit_bold, X, psych_idx, h = fmrihrf::spm_hrf(TR))
neural_delay     <- mppi_neural_from_fit(fit_bold, X, psych_idx, h = fmrihrf::spm_hrf(TR, onset = 0.5))
neural_disp      <- mppi_neural_from_fit(fit_bold, X, psych_idx, h = fmrihrf::spm_hrf(TR, p = c(6,18,1,1,6,0,32)))
neural_avg       <- mppi_hrf_ensemble(list(neural_canonical, neural_delay, neural_disp))
```

## Notes

- This implements the **Gitelman principle** (form interactions in **neural** space) but keeps the
  multiPPI estimator and your **shared basis** (scalable). See the original note for why deconvolution
  matters (esp. event-related designs). [Gitelman et al., 2003].  
- It preserves the **gPPI principle** (span the experimental space) by residualizing the **unconvolved**
  psychological vectors against each other (multi-condition). [McLaren et al., 2012].

