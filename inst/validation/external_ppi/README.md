# External PPI regression oracle

These fixtures are synthetic. The test compares Nilearn regression and contrasts conditional on the exported multiPPI design; it does not validate construction of that design. See the [audit receipt](../../../audit/2026-10-23-ppi/external-comparison.md) for interpretation, version-bound findings, source hashes, and the SPM construction comparison.

To regenerate against a frozen source copy, install the R dependencies and run:

```sh
Rscript inst/validation/external_ppi/export_fixture.R /path/to/frozen/source /tmp/multippi-oracle-fixture
python3 -m venv /tmp/multippi-nilearn-env
/tmp/multippi-nilearn-env/bin/pip install nilearn==0.14.1 numpy==2.5.4 scipy==1.18.1
OPENBLAS_NUM_THREADS=1 /tmp/multippi-nilearn-env/bin/python inst/validation/external_ppi/compare_nilearn.py /tmp/multippi-oracle-fixture /tmp/multippi-nilearn-results.json
```

For a regression-only replay of the retained fixtures, replace `/tmp/multippi-oracle-fixture` with `inst/validation/external_ppi/fixtures`. The comparison exits nonzero on a failed check. `spm_octave_probe.m` calls the actual pinned SPM routine with dependencies on the Octave path; it writes its outputs to `/tmp/multippi-spm`. It succeeded with the bounded official dependencies and no graphics toolkit in this session. Compare its saved vectors with the frozen multiPPI source using:

```sh
Rscript inst/validation/external_ppi/compare_spm.R /path/to/frozen/source inst/validation/external_ppi/fixtures /tmp/multippi-spm-comparison.json
```

Both raw-seed and default nuisance-adjusted paths are compared. The intentionally different PEB and ridge constructions are described in the audit receipt; matching the task convolution does not establish interaction equivalence.
