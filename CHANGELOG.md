# v1.0.0

Initial public release.

- Provides portable analysis entry points for SNP-level wild minor-allele representation, inter-garden complementarity, AC-threshold and geographic-gap summaries, sample-size-standardized retention, callability/carrier quality, PI_HAT connected components and deterministic core selection, SNP ALT-state representation, quota-constrained supplementation, region-matched random controls, selected-sample QC, and direct final-portfolio measurement.
- Uses PI_HAT >=0.80 connected components with lexicographic SampleID representative selection and the study's fixed supplementation, regional-quota and random-stream conventions.
- Distinguishes complete diploid genotypes used during supplementation selection from called-allele counting used for SNP ALT and final-portfolio measurement.
- Includes synthetic end-to-end and allele edge-case tests, package-integrity checks, and protection against silent overwrite of nonempty output directories.
- Before public release, the complete research workflow was rerun on the full research inputs and reproduced the configured manuscript anchors.
- Excludes plotting code, manuscripts, large research VCFs and real sample metadata from the public source distribution.
