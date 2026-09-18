# Validation and reproducibility

## Completed release validation

- Package SHA256 verification, Python syntax checks and shell entry-point syntax checks passed.
- All ten portable analysis steps and final cross-step validation passed on the packaged synthetic dataset.
- Hand-calculated ALT edge cases passed, including mixed SNP/`*` records, exclusion of `*` carried only by a removed sample, partial called genotypes, lowercase bases and identical REF/ALT exclusion.
- Transitive PI_HAT connected components, lexicographic representative selection, out-of-range GT rejection, overwrite protection and four-portfolio loss/recovery algebra passed the packaged tests.
- Before public release, the consolidated workflow was rerun on the complete research inputs, including the full 601-sample configuration and the research VCFs required by the individual steps. All ten analysis steps and the final cross-step/manuscript-anchor validation completed successfully.
- The full research-data rerun reproduced the configured study anchors used by this release, including the wild minor-allele target universe, Current/Core/Plan A/Plan B representation, the 286-accession de-redundant core, Plan A/Plan B selection totals, SNP ALT-state and copy-weighted metrics, region-matched random-control summaries and selected-sample QC checks.

## Research-data provenance and distribution scope

The adopted research anchors were cross-checked against the reviewed September 2026 study archives and reproduced by the consolidated public workflow before release. Expected values in the configuration are used only for validation after calculation; they are never substituted for calculated results.

The public source distribution contains code, configuration templates and synthetic test data only. The full research VCFs, real sample metadata and other large or derived research inputs are not bundled and must be provided separately as described in [Input schema](INPUT_SCHEMA.md).
