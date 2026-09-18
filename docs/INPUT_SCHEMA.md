# Input schema

All identifiers are case-sensitive. Coordinates are matched by exact `(CHROM, POS)` pairs.

## AC/AN tables

Configured group keys are `wild`, `ex_situ`, `pingshun`, `tongguan` and `zhendong`, plus one regional file for each Wild region. Files may be plain text or gzip-compressed and contain:

```text
CHROM    POS    ALT_AC    AN
```

`0 <= ALT_AC <= AN` is required. The five genome-wide group files must have the same number of rows in the same coordinate order. Regional files must contain every resource-garden-missing target exactly once and no non-target coordinate.

## Metadata

Required columns:

```text
SampleID    Material_type    Ex_situ_site    Region
```

Accepted values are `Wild` or `Ex_situ` for `Material_type`; `Pingshun`, `Tongguan` or `Zhendong` for resource-garden accessions; and one of the six configured Wild regions for field-collected samples. Optional provenance columns retained in outputs are `Original_province`, `Original_city` and `Original_county`.

## Sample lists

Each file contains one SampleID per non-empty line. Wild and resource-garden lists must be disjoint. Pingshun, Tongguan and Zhendong must be disjoint and must partition the combined resource-garden list. The fixed Wild sample-list order is used as the final Greedy tie-breaker.

## Missing-site VCF

Requirements:

- all resource-garden-missing target coordinates;
- all configured Wild samples;
- `GT` for Greedy analysis;
- `GT`, `DP` and `GQ` for strict carrier-quality checking;
- biallelic genotype codes `0` and `1`.

Additional samples are allowed; only configured Wild samples are used for target-carrier calculations.

## Step 06 cohort-level PASS SNP VCF

Configured key: `pass_snp_vcf` (internally aliased to `full_snp_vcf`).

The verified project path is:

```text
01_qc/04.samples_filtered.PASS.snp.vcf.gz
```

Requirements:

- the exact unified cohort-level PASS SNP set used for the within-collection de-redundancy ALT analysis;
- all 601 samples, including every current resource-garden accession;
- `GT`;
- 78,184,349 records for the verified manuscript run.

This input is not interchangeable with the deposited 25,502,888-site diversity-analysis VCF.

## PI_HAT pair table

Required columns, with accepted aliases:

```text
IID1    IID2    PI_HAT
```

Every unordered pair may appear at most once. Values must lie in `[0,1]`.

## ROH summary

Required column: `SampleID`. The selected-accession >2 Mb check additionally requires either length-class count fields such as `NSEG_2_5Mb` and `NSEG_gt5Mb`, or a recognized maximum-segment-length field. Missing selected-accession ROH records cause failure.

## Regional quotas and random order

`config/regional_quotas.tsv` contains:

```text
Plan    Region    Sample_N    Random_sampling_order
```

Every plan must specify all six regions, including zero quotas. Plan B quotas must be greater than or equal to Plan A quotas region by region. `Random_sampling_order` is part of exact random-stream reproducibility; zero-quota regions are recorded but do not consume random numbers.

## Step 10 diversity SNP VCF

Key `diversity_snp_vcf`: `02_filter/03.diversity.miss0.8.bi.noMAF.polymorphic.vcf.gz`.
This is the 25,502,888-record biallelic diversity input, with all Wild and Current samples and GT; records must be coordinate-sorted and unique. It is distinct from the 78,184,349-record PASS input, which can include multiallelic and spanning-deletion states. Never rename one input to satisfy the other's path.

The AC/AN tables for Steps 01–03 must be derived from the diversity input and its declared sample lists. The VCF minor-allele target and AC/AN target must therefore refer to the same sites and alleles. No upstream filtering or historical AC/AN-export command is reconstructed by this release.

For Step 09 use a ROH summary with >2 Mb counts or maximum lengths. This release reports the threshold checks and does not regenerate full ROH tracts or the precision-corrected total-length columns of Supplementary Table 14. Exact length values for the paper must be taken from the reviewed Results3 outputs.

ROH/PI_HAT are external derived inputs: the package does not run PLINK to estimate them. Sample lists and order must remain unchanged for exact selection and random-control reproducibility. Historical serialized fields `Material_type=Ex_situ` and `Ex_situ_site` remain required even though prose terminology is resource-garden.
