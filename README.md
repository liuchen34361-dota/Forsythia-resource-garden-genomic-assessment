# Forsythia resource-garden genomic assessment — v1.0.0

Custom code for SNP-level allele representation, de-redundancy and targeted supplementation.
Repository: https://github.com/liuchen34361-dota/Forsythia-resource-garden-genomic-assessment

Scientific terminology here uses resource-garden collections and field-collected wild samples. Some input keys and TSV column names containing `ex_situ` are retained only for compatibility with the study files.

## Scope

The ten analysis steps cover wild minor-allele representation and inter-garden complementarity; AC-threshold and geographic-gap summaries; sample-size-standardized retention; callability/carrier quality; PI_HAT connected components and deterministic core selection; within-collection SNP ALT-state representation; quota-constrained greedy supplementation; region-matched random controls; selected-sample QC; and direct final-portfolio representation.

This is a portable source-code release, with a synthetic example. It does not include sequencing data, research VCFs, real sample metadata, plotting code, manuscripts, software binaries, or R packages. Read [Input schema](docs/INPUT_SCHEMA.md) before running the research configuration. Required inputs must be supplied separately; accession numbers alone do not provide all derived inputs.

## Quick test

Python 3.11 or newer and Bash are sufficient for the packaged synthetic tests:

```bash
python tests/run_synthetic_smoke_test.py
python tests/test_allele_edge_cases.py
python scripts/verify_package.py
```

The tests create temporary working directories and do not modify the example data.

## Research run

Use Linux/WSL. Python scripts use the standard library. Step 06 uses bcftools to stream the Current sample subset; genotype-derived counts are calculated in Python. Its slower direct reader is available explicitly via the configuration.

```bash
conda env create -f environment.yml
conda activate forsythia_assessment
bash /absolute/path/to/code/run_pipeline.sh /absolute/path/to/code/config/analysis_config.json /absolute/path/to/project
```

First edit a copy of `config/analysis_config.json` to match your input paths. With `project_root` set to `.`, paths resolve against the project directory supplied to the shell entry point. An absolute `project_root` in JSON overrides that location for data resolution. Each output step refuses an existing nonempty output directory. Use a new `output_root` for a rerun. These portable scripts do not provide the chromosome-checkpoint resume mechanism of the study's internal execution wrappers.

Large VCF steps can take many hours; do not use the research configuration to test installation. Each script accepts `--config`; invoke individual steps only after their prerequisites have completed. Expected results validate outputs after computation; they are never substituted for calculated values. For an independent dataset, change the quotas and metadata configuration and clear the two study-specific expected-result dictionaries.

## Key definitions

- Wild baseline: polymorphic biallelic sites in the no-MAF diversity dataset; minor = min(AC, AN−AC), with ALT chosen on a tie.
- Core: connected components of the resource-garden-only PI_HAT graph at ≥0.80; keep the lexicographically first SampleID in each component, including singletons.
- SNP ALT metrics: unequal single-base A/C/G/T REF/ALT states only. Exclude `*` individually; retain eligible SNP states at the same multiallelic site. Current AC is the copy weight; Core presence determines whether that weight is covered.
- Selection: complete diploid target genotypes, low-frequency marginal-contribution weights 1/2/2, then total marginal gain, call rate, and fixed Wild list order. Plan B extends Plan A.
- Direct portfolio measurement: called alleles in partial GT count. This is distinct from the complete-genotype selection convention.
- Random controls: 1,000 sets per plan, seed 20260620, preserved sample/region order and historical random-stream consumption.

## Verified research anchors

| Portfolio | Samples | Missing wild targets | Retention (%) |
|---|---:|---:|---:|
| Current | 334 | 180,796 | 99.2604 |
| Core | 286 | 247,641 | 98.9870 |
| Plan A | 329 | 140,023 | 99.4272 |
| Plan B | 397 | 64,729 | 99.7352 |

Target universe: 24,445,213. Core-induced loss: 66,845; Plans A/B rescue 28,283/49,608 of these, in addition to recovering 79,335/133,304 originally missing targets.

Within-collection SNP ALT states: 62,594,768 Current; 62,112,843 retained; 481,925 lost. Presence retention is 99.2301%; copy-weighted coverage is 99.9807%. The two metrics use different denominators and are separate from the wild-baseline analysis.

## Validation and provenance

[Validation report](docs/VALIDATION.md) summarizes the packaged tests and the completed full research-data reproduction. [Release notes](CHANGELOG.md) summarize this initial public release. Before public release, the consolidated workflow was rerun on the complete research inputs and reproduced the configured manuscript anchors across allele representation, de-redundancy, supplementation, random controls, selected-sample QC, SNP ALT metrics and direct final-portfolio measurement.

Zenodo DOI: not yet assigned to this release. Do not describe this release as permanently archived until publication of the matching Zenodo record. This package does not assign a license or infer the author list.
