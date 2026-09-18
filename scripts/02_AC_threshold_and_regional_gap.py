#!/usr/bin/env python3
"""Regional occurrence and normalized gaps of ex situ-missing wild alleles."""
from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

from analysis_utils import (
    prepare_step_output,
    MAF_BIN_ORDER,
    THRESHOLD_ORDER,
    THRESHOLDS,
    canonical_region,
    default_config_path,
    input_path,
    load_config,
    load_missing_targets,
    iter_acan_rows,
    outputs_root,
    read_sample_list,
    read_tsv,
    step_output_dir,
    target_allele_count,
    write_key_value,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--missing-table", default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    out = prepare_step_output(cfg, "02_regional_occurrence_and_gap")
    missing_path = Path(args.missing_table) if args.missing_table else (
        outputs_root(cfg) / "01_allele_representation" / "05.ex_situ_missing_wild_minor_alleles.tsv.gz"
    )
    targets, order = load_missing_targets(missing_path)

    metadata_rows = read_tsv(input_path(cfg, "metadata"))
    metadata: dict[str, dict[str, str]] = {}
    for row in metadata_rows:
        sid = row.get("SampleID", "").strip()
        if not sid:
            continue
        if sid in metadata:
            raise ValueError(f"Duplicated SampleID in metadata: {sid}")
        metadata[sid] = row
    wild_samples = read_sample_list(input_path(cfg, "wild_samples"))
    absent = [sid for sid in wild_samples if sid not in metadata]
    if absent:
        raise ValueError(f"Wild sample IDs absent from metadata: {absent[:10]}")
    region_sample_n = Counter(canonical_region(metadata[sid].get("Region")) for sid in wild_samples)
    configured_counts = cfg.get("parameters", {}).get("region_sample_counts", {}) or {}
    for raw_region, expected_n in configured_counts.items():
        region = canonical_region(raw_region)
        observed_n = region_sample_n.get(region, 0)
        if observed_n != int(expected_n):
            raise ValueError(
                f"Wild region sample-count mismatch for {region}: "
                f"metadata={observed_n}, configured={expected_n}"
            )

    configured_region_files = cfg["inputs"].get("regional_acan", {})
    if not configured_region_files:
        raise ValueError("No regional_acan mapping is present in the config")
    # Regional paths are stored in a nested mapping. Canonicalization must not
    # collapse two configured aliases onto the same biological region.
    region_files: dict[str, Path] = {}
    for region, path_value in configured_region_files.items():
        canonical = canonical_region(region)
        if canonical in region_files:
            raise ValueError(f"Duplicate regional AC/AN mapping after name normalization: {canonical}")
        path = Path(path_value)
        region_files[canonical] = path if path.is_absolute() else Path(cfg["_project_root"]) / path

    missing_region_counts = [region for region in region_files if region_sample_n[region] <= 0]
    if missing_region_counts:
        raise ValueError(f"No wild metadata samples for regions: {missing_region_counts}")
    unconfigured_regions = sorted(set(region_sample_n) - set(region_files))
    if unconfigured_regions:
        raise ValueError(f"Wild metadata contains regions without AC/AN inputs: {unconfigured_regions}")

    presence: dict[tuple[str, str], set[str]] = {key: set() for key in order}
    regional_counts: dict[tuple[str, str], dict[str, int]] = {key: {} for key in order}
    region_rows_scanned = Counter()
    unmatched_rows = Counter()

    for region, path in region_files.items():
        seen_region_sites: set[tuple[str, str]] = set()
        for row in iter_acan_rows(path, region):
            key = (row.chrom, row.pos)
            region_rows_scanned[region] += 1
            if key not in targets:
                unmatched_rows[region] += 1
                continue
            if key in seen_region_sites:
                raise ValueError(f"Duplicate {region} AC/AN site: {row.chrom}:{row.pos}")
            seen_region_sites.add(key)
            count = target_allele_count(row, targets[key]["minor"])
            regional_counts[key][region] = count
            if count > 0:
                presence[key].add(region)
        missing_for_region = set(targets) - seen_region_sites
        if missing_for_region:
            preview = ", ".join(f"{c}:{pos}" for c, pos in sorted(missing_for_region)[:10])
            raise ValueError(
                f"Regional AC/AN table for {region} lacks {len(missing_for_region)} missing targets; "
                f"examples: {preview}"
            )
        if unmatched_rows[region]:
            raise ValueError(f"{region} AC/AN table contains {unmatched_rows[region]} non-target rows")

    occurrence = defaultdict(int)
    occurrence_maf = defaultdict(int)
    region_specific = defaultdict(int)
    site_path = out / "01.missing_allele_regional_occurrence_by_site.tsv.gz"
    site_fields = [
        "CHROM", "POS", "Wild_minor_allele", "Wild_minor_count", "Wild_minor_MAF",
        "Wild_MAF_bin", "Presence_region_count", "Presence_regions", "Region_specific_region",
    ] + [f"{region}_minor_count" for region in region_files]
    with gzip.open(site_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=site_fields)
        writer.writeheader()
        for key in order:
            target = targets[key]
            regions = sorted(presence[key])
            n_region = len(regions)
            specific = regions[0] if n_region == 1 else "NA"
            for threshold_name in THRESHOLD_ORDER:
                if target["minor_count"] < THRESHOLDS[threshold_name]:
                    continue
                occurrence[(threshold_name, n_region)] += 1
                occurrence_maf[(threshold_name, target["maf_bin"], n_region)] += 1
                if n_region == 1:
                    region_specific[(threshold_name, specific)] += 1
            row = {
                "CHROM": key[0],
                "POS": key[1],
                "Wild_minor_allele": target["minor"],
                "Wild_minor_count": target["minor_count"],
                "Wild_minor_MAF": target["maf"],
                "Wild_MAF_bin": target["maf_bin"],
                "Presence_region_count": n_region,
                "Presence_regions": ";".join(regions) if regions else "None",
                "Region_specific_region": specific,
            }
            for region in region_files:
                row[f"{region}_minor_count"] = regional_counts[key].get(region, 0)
            writer.writerow(row)

    occurrence_rows = []
    for threshold_name in THRESHOLD_ORDER:
        for n_region in range(0, len(region_files) + 1):
            occurrence_rows.append({
                "AC_threshold": threshold_name,
                "Wild_minor_count_min": THRESHOLDS[threshold_name],
                "Presence_region_count": n_region,
                "N_missing_Wild_minor_alleles": occurrence[(threshold_name, n_region)],
            })
    write_tsv(
        out / "02.region_presence_by_AC_threshold.tsv",
        occurrence_rows,
        ["AC_threshold", "Wild_minor_count_min", "Presence_region_count", "N_missing_Wild_minor_alleles"],
    )

    occurrence_maf_rows = []
    for threshold_name in THRESHOLD_ORDER:
        for mbin in MAF_BIN_ORDER:
            for n_region in range(0, len(region_files) + 1):
                occurrence_maf_rows.append({
                    "AC_threshold": threshold_name,
                    "Wild_minor_count_min": THRESHOLDS[threshold_name],
                    "Wild_MAF_bin": mbin,
                    "Presence_region_count": n_region,
                    "N_missing_Wild_minor_alleles": occurrence_maf[(threshold_name, mbin, n_region)],
                })
    write_tsv(
        out / "03.region_presence_by_AC_threshold_and_MAF.tsv",
        occurrence_maf_rows,
        ["AC_threshold", "Wild_minor_count_min", "Wild_MAF_bin", "Presence_region_count", "N_missing_Wild_minor_alleles"],
    )

    gap_rows = []
    for threshold_name in THRESHOLD_ORDER:
        for region in region_files:
            count = region_specific[(threshold_name, region)]
            sample_n = region_sample_n[region]
            gap_rows.append({
                "AC_threshold": threshold_name,
                "Wild_minor_count_min": THRESHOLDS[threshold_name],
                "Region": region,
                "Wild_sample_N": sample_n,
                "Region_specific_missing_Wild_minor_alleles": count,
                "Region_specific_gap_per_sample": count / sample_n,
            })
    write_tsv(
        out / "04.region_specific_gap_normalized_by_AC_threshold.tsv",
        gap_rows,
        ["AC_threshold", "Wild_minor_count_min", "Region", "Wild_sample_N",
         "Region_specific_missing_Wild_minor_alleles", "Region_specific_gap_per_sample"],
    )

    zero_region_targets = occurrence[("all", 0)]
    if zero_region_targets != 0:
        raise ValueError(f"{zero_region_targets} missing targets were absent from all Wild region tables")

    check = [
        ("Missing_targets_loaded", len(targets)),
        ("Wild_metadata_samples", sum(region_sample_n.values())),
        ("Regions", len(region_files)),
        ("Targets_detected_in_zero_regions_all", occurrence[("all", 0)]),
        ("Region_specific_targets_all", sum(region_specific[("all", r)] for r in region_files)),
    ]
    for region in region_files:
        check.extend([
            (f"{region}_sample_N", region_sample_n[region]),
            (f"{region}_ACAN_rows_scanned", region_rows_scanned[region]),
            (f"{region}_unmatched_rows", unmatched_rows[region]),
        ])
    write_key_value(out / "05.internal_check.tsv", check)
    print(f"[OK] Regional occurrence analysis complete: {out}")


if __name__ == "__main__":
    main()
