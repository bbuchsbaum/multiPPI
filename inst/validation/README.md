# Validation scripts

`validate_ppi.R` runs fixed-design conditional regression simulations with IID,
IID censoring, a fixed AR(1) operator, estimated AR(1), and a heterogeneous
shared-operator misspecification probe. `validate_ppi_grouped.R` runs a separate
targetwise heterogeneous-AR null screen. Both use synthetic data, fixed seeds,
source hashes, exact binomial intervals and prespecified thresholds; run them
from the source root with the declared dependencies installed. They write
receipts under `audit/2026-10-23-ppi/`, which is excluded from package builds.

These scripts test conditional regression for specific fixed designs and noise
scenarios. They do not validate latent neural recovery, HRF mismatch, arbitrary
serial noise, selection or multiplicity procedures, or causal inference.
Targetwise fitting is substantially more expensive than shared pooling.

`validate_ppi_generative.R` separately uses an independent dense neural-to-BOLD
generator to screen interaction construction, with matched activation-only and
constant-coupling nulls, a neural interaction effect, and adverse HRF-delay and
seed-noise probes. Its 200 neural realizations each have ten target-noise draws;
uncertainty is summarized by seed-cluster bootstrap. It distinguishes the true
neural coefficient from its projection onto the constructed design. A positive
effect check does not establish accurate recovery of the generating coefficient.

The [external oracle guide](external_ppi/README.md) explains replaying the
Nilearn regression comparisons and the actual SPM PPI construction probe. These
use pinned versions and synthetic inputs, with differences and qualifications
retained in the linked audit receipts.
