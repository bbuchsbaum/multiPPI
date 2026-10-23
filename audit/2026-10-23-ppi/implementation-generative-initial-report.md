> Superseded pilot receipt: this initial run omitted three required scenarios, used the wrong nuisance and seed-noise scaling, and did not use an independent projection oracle. It does not execute the full approved recipe. See implementation-generative-report.md for the corrected run.

# Bounded generative PPI construction screen

Executed 2026-10-23 with R 4.5.1, aarch64-apple-darwin20, Apple Accelerate BLAS, 26 GiB free before run. A three-seed pilot completed first; the full run exited 0. The script was saved before the pilot and its SHA-256 was unchanged at completion: `8addba67d83152eb5f4f8676133935d77ab43a68528c6486dce795b8e01cc415`. Both production source hashes were unchanged from start to finish: `ppi_design.R` `a710631460dac0ec6713b0a005d9d6d57a1e332bc1d3c67bb01de0f31d46f641`; `ppi_fit.R` `f9255db93f93bf855b4862a5cc7c4184fa303f78afd5b7c06b7c3cd886e7f930`.

The script fixes 200 independent stationary neural AR realizations, with ten independent target-noise draws per realization, two runs of 90 and 110 scans, TR 2 s and microtime 4. The block conditions, run shift, nuisance, seed and target construction, fixed ridge, nuisance adjustment, IID fit and A-minus-B interaction contrast are specified in the executable script. The causal sampled operator for the generating model is built as a dense matrix independently of the package convolution. The delayed target HRF is shifted by 1 s within the fixed 0–32 s support; seed measurement noise has SD one quarter of the pooled sampled neural seed SD. Latent seeds and target-noise draws are paired across scenarios. The seed-noise draws are also generated for each seed regardless of scenario.

For each seed and scenario, `b*` is the A-minus-B coefficient of the full constructed design projected onto the noiseless target mean, not the generating interaction amplitude 0.5. Coverage is against this scenario-specific projection. True and constructed run-by-condition interaction columns are separately residualized against all noninteraction columns; correlations and relative Frobenius error compare these residualized columns. Full design rank was 14 throughout; residualized oracle interaction rank was 4 throughout. All contrasts were estimable. Correlation and error evaluate construction fidelity, and are not statistical calibration measures.

| Scenario | Mean b* (seed bootstrap 95%) | Rejections / 2000 | Rejection proportion (seed bootstrap 95%) | Coverage of b* |
|---|---:|---:|---:|---:|
| Matched null | approximately 0 | 91 | 0.0455 (0.0365, 0.0555) | 0.9545 |
| Matched effect | 0.6686 (0.6550, 0.6825) | 1104 | 0.5520 (0.5225, 0.5825) | 0.9590 |
| Delayed target HRF | 0.6072 (0.5617, 0.6522) | 787 | 0.3935 (0.3510, 0.4360) | 0.9805 |
| Noisy seed | 0.6447 (0.6204, 0.6702) | 1020 | 0.5100 (0.4750, 0.5460) | 0.9675 |

Matched-null maximum |b*| divided by maximum absolute target mean was 1.89e-15, below the prespecified 1e-8 bound. Matched effect mean b* has positive lower bootstrap limit and its rejection proportion exceeds the null by 0.5065; the paired seed-cluster bootstrap interval for that difference is [0.4755, 0.5385]. These are bounded screen criteria, not proof of general validity.

Construction accuracy was imperfect even in the matched setup: mean residualized interaction correlation was 0.7589, the minimum over seeds and columns was 0.0559, mean relative error was 0.6571, and maximum relative error was 0.8141. Noisy-seed mean correlation was 0.7476 and mean relative error 0.6656. These adverse results are retained rather than treating detectable effects as accurate neural reconstruction.

The 5,000 fixed-seed bootstrap draws resample whole neural-seed clusters, retaining all ten targets together and pairing scenarios. Their intervals describe the finite simulated seed bank and its cluster variability; the 2,000 target outcomes are not independent draws of new seeds, and the rejection rates are not an IID binomial calibration guarantee. The two adverse scenarios contain conditional mean misspecification, so projection coverage and rejection results are descriptive for this recipe only. No conclusion is made about estimated temporal noise, censoring, population data, or alternative designs.

Companion `implementation-generative-seeds.csv` contains all 800 seed/scenario records and rejection/coverage counts. `implementation-generative-summary.csv` contains aggregate metrics; `implementation-generative-bootstrap.rds` retains all 5,000 bootstrap mean and rejection draws; `implementation-generative-run.log` records execution. Reproduce with `OPENBLAS_NUM_THREADS=1 Rscript inst/validation/validate_ppi_generative.R 200 audit/2026-10-23-ppi/implementation-generative` from the repository root.
