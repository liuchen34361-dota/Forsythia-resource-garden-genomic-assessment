#!/usr/bin/env python3
"""Quota-constrained greedy Plan A/Plan B selection and recovery robustness.

Candidate ranking follows the manuscript exactly: weighted marginal recovery,
then unweighted total marginal recovery, then target-site call rate, and finally
the fixed Wild sample-list order.  Plan B is initialized with Plan A and then
continued to the expanded quotas.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path
from typing import Any

from analysis_utils import (
    prepare_step_output,
    MAF_BIN_ORDER,
    THRESHOLD_ORDER,
    THRESHOLDS,
    build_target_bitsets_from_vcf,
    canonical_region,
    default_config_path,
    input_path,
    load_config,
    load_missing_targets,
    package_path,
    read_sample_list,
    read_tsv,
    step_output_dir,
    write_key_value,
    write_tsv,
)

PLAN_A = "Plan_A_conservative"
PLAN_B = "Plan_B_expanded"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--missing-table", default=None)
    return parser.parse_args()


def load_quotas(path: Path) -> dict[str, OrderedDict[str, int]]:
    rows = read_tsv(path)
    quotas: dict[str, OrderedDict[str, int]] = defaultdict(OrderedDict)
    for row in rows:
        plan = row["Plan"]
        region = canonical_region(row["Region"])
        n = int(row.get("Sample_N", row.get("Quota", "0")))
        if region in quotas[plan]:
            raise ValueError(f"Duplicate quota for {plan}, {region}")
        quotas[plan][region] = n
    for plan in (PLAN_A, PLAN_B):
        if plan not in quotas:
            raise ValueError(f"Quota file lacks {plan}")
    return dict(quotas)


def build_masks(targets: dict[tuple[str, str], dict[str, Any]]) -> dict[tuple[str, str], int]:
    masks: dict[tuple[str, str], int] = defaultdict(int)
    for info in targets.values():
        bit = 1 << int(info["idx"])
        for threshold in THRESHOLD_ORDER:
            if int(info["minor_count"]) < THRESHOLDS[threshold]:
                continue
            masks[(threshold, "All")] |= bit
            masks[(threshold, str(info["maf_bin"]))] |= bit
    return dict(masks)


def recovery_rows(
    plan: str,
    samples: list[str],
    bits: int,
    masks: dict[tuple[str, str], int],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for threshold in THRESHOLD_ORDER:
        for bin_label in ("All", *MAF_BIN_ORDER):
            mask = masks.get((threshold, bin_label), 0)
            total = mask.bit_count()
            recovered = (bits & mask).bit_count()
            rows.append({
                "Plan": plan,
                "Wild_supplement_N": len(samples),
                "AC_threshold": threshold,
                "Wild_minor_count_min": THRESHOLDS[threshold],
                "Wild_MAF_bin": bin_label,
                "Missing_Wild_minor_alleles_total": total,
                "Recovered_missing_Wild_minor_alleles": recovered,
                "Recovery_rate": recovered / total if total else 0.0,
            })
    return rows


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    out = prepare_step_output(cfg, "07_greedy_plan_selection")
    if args.missing_table:
        missing_path = Path(args.missing_table)
    else:
        missing_path = step_output_dir(cfg, "01_allele_representation") / "05.ex_situ_missing_wild_minor_alleles.tsv.gz"
    targets, target_order = load_missing_targets(missing_path)
    masks = build_masks(targets)

    wild_samples = read_sample_list(input_path(cfg, "wild_samples"))
    metadata_rows = read_tsv(input_path(cfg, "metadata"))
    metadata: dict[str, dict[str, str]] = {}
    for row in metadata_rows:
        sample = row.get("SampleID", "").strip()
        if not sample:
            continue
        if sample in metadata:
            raise ValueError(f"Duplicated SampleID in metadata: {sample}")
        metadata[sample] = row
    missing_metadata = [sample for sample in wild_samples if sample not in metadata]
    if missing_metadata:
        raise ValueError(f"Wild samples missing from metadata: {missing_metadata[:10]}")
    sample_region = {sample: canonical_region(metadata[sample].get("Region")) for sample in wild_samples}
    if any(region == "NA" for region in sample_region.values()):
        raise ValueError("Every Wild sample must have a region")
    region_to_samples: dict[str, list[str]] = defaultdict(list)
    for sample in wild_samples:
        region_to_samples[sample_region[sample]].append(sample)

    quotas = load_quotas(package_path(cfg, cfg["package_files"]["regional_quotas"]))
    wild_regions = set(region_to_samples)
    for plan, quota in quotas.items():
        quota_regions = set(quota)
        if quota_regions != wild_regions:
            missing = sorted(wild_regions - quota_regions)
            extra = sorted(quota_regions - wild_regions)
            raise ValueError(
                f"Quota regions do not match Wild metadata for {plan}: "
                f"missing={missing}, extra={extra}"
            )
        for region, n in quota.items():
            if n < 0 or n > len(region_to_samples.get(region, [])):
                raise ValueError(f"Invalid quota for {plan}, {region}: {n}")
    for region, plan_a_n in quotas[PLAN_A].items():
        if quotas[PLAN_B].get(region, 0) < plan_a_n:
            raise ValueError(f"Plan B quota is smaller than Plan A for {region}")

    vcf = input_path(cfg, "missing_site_vcf")
    bitset_result = build_target_bitsets_from_vcf(vcf, wild_samples, targets)
    if bitset_result.missing_targets:
        preview = ", ".join(
            f"{chrom}:{pos}" for chrom, pos in bitset_result.missing_targets[:10]
        )
        raise ValueError(
            f"Missing-site VCF lacks {len(bitset_result.missing_targets)} targets; examples: {preview}"
        )
    sample_bits = bitset_result.sample_bits
    called_sites = bitset_result.called_sites
    matched_n = bitset_result.processed_targets

    target_n = len(targets)
    call_rate = {sample: called_sites[sample] / target_n for sample in wild_samples}
    sample_order_index = {sample: index + 1 for index, sample in enumerate(wild_samples)}
    call_rows = [{
        "SampleID": sample,
        "Region": sample_region[sample],
        "Called_target_sites": called_sites[sample],
        "Total_target_sites": target_n,
        "Target_site_call_rate": call_rate[sample],
    } for sample in wild_samples]
    write_tsv(out / "00.Wild_target_site_call_rate.tsv", call_rows, list(call_rows[0]))

    weights = {key: float(value) for key, value in cfg.get("parameters", {}).get("greedy_weights", {}).items()}
    for key in ("MAF_le_0.01", "0.01_lt_MAF_le_0.05", "0.05_lt_MAF_le_0.10"):
        if key not in weights:
            raise ValueError(f"Missing greedy weight: {key}")

    def marginal(sample: str, covered: int) -> tuple[float, int, dict[str, int]]:
        new_bits = sample_bits[sample] & ~covered
        counts = {"All": new_bits.bit_count()}
        for bin_label in MAF_BIN_ORDER:
            counts[bin_label] = (new_bits & masks.get(("all", bin_label), 0)).bit_count()
        score = sum(counts[key] * weights.get(key, 0.0) for key in MAF_BIN_ORDER)
        return score, counts["All"], counts

    def select(
        quota: OrderedDict[str, int],
        initial_samples: list[str],
        initial_covered: int,
        stage: str,
    ) -> tuple[list[str], int, list[dict[str, Any]]]:
        selected = list(initial_samples)
        selected_set = set(selected)
        covered = initial_covered
        region_counts = Counter(sample_region[sample] for sample in selected)
        rows: list[dict[str, Any]] = []
        while any(region_counts[region] < required_n for region, required_n in quota.items()):
            best_sample: str | None = None
            best_key: tuple[float, int, float] | None = None
            best_counts: dict[str, int] | None = None
            for sample in wild_samples:  # fixed input order resolves the final tie
                if sample in selected_set:
                    continue
                region = sample_region[sample]
                if region not in quota or region_counts[region] >= quota[region]:
                    continue
                score, total_new, counts = marginal(sample, covered)
                key = (score, total_new, call_rate[sample])
                if best_key is None or key > best_key:
                    best_sample = sample
                    best_key = key
                    best_counts = counts
            if best_sample is None or best_key is None or best_counts is None:
                raise RuntimeError(f"No candidate available while filling {stage}")
            covered |= sample_bits[best_sample]
            selected.append(best_sample)
            selected_set.add(best_sample)
            region = sample_region[best_sample]
            region_counts[region] += 1
            meta = metadata[best_sample]
            rows.append({
                "Selection_stage": stage,
                "Plan_step": len(selected),
                "SampleID": best_sample,
                "Region": region,
                "Original_province": meta.get("Original_province", "NA"),
                "Original_city": meta.get("Original_city", "NA"),
                "Original_county": meta.get("Original_county", "NA"),
                "Region_selected_N_after_step": region_counts[region],
                "Region_quota": quota[region],
                "Marginal_weighted_score": best_key[0],
                "Marginal_new_All": best_counts["All"],
                "Marginal_new_MAF_le_0.01": best_counts["MAF_le_0.01"],
                "Marginal_new_0.01_lt_MAF_le_0.05": best_counts["0.01_lt_MAF_le_0.05"],
                "Marginal_new_0.05_lt_MAF_le_0.10": best_counts["0.05_lt_MAF_le_0.10"],
                "Target_site_call_rate": call_rate[best_sample],
                "Tie_break_order_index": sample_order_index[best_sample],
            })
        return selected, covered, rows

    plan_a_samples, plan_a_bits, plan_a_rows = select(quotas[PLAN_A], [], 0, "Greedy_Plan_A")
    plan_b_samples, plan_b_bits, plan_b_added_rows = select(
        quotas[PLAN_B], plan_a_samples, plan_a_bits, "Greedy_Plan_B_additional"
    )
    plan_b_rows: list[dict[str, Any]] = []
    a_lookup = {row["SampleID"]: row for row in plan_a_rows}
    for step, sample in enumerate(plan_a_samples, start=1):
        inherited = dict(a_lookup[sample])
        inherited["Selection_stage"] = "Inherited_from_Greedy_Plan_A"
        inherited["Plan_step"] = step
        inherited["Region_quota"] = quotas[PLAN_B][sample_region[sample]]
        plan_b_rows.append(inherited)
    plan_b_rows.extend(plan_b_added_rows)

    selection_fields = list(plan_a_rows[0])
    write_tsv(out / "01.Greedy_Plan_A_Wild_supplements.tsv", plan_a_rows, selection_fields)
    write_tsv(out / "02.Greedy_Plan_B_Wild_supplements.tsv", plan_b_rows, selection_fields)

    all_recovery_rows = (
        recovery_rows("Greedy_Plan_A", plan_a_samples, plan_a_bits, masks)
        + recovery_rows("Greedy_Plan_B", plan_b_samples, plan_b_bits, masks)
    )
    write_tsv(
        out / "03.Greedy_plan_recovery_by_AC_threshold_and_MAF.tsv",
        all_recovery_rows,
        list(all_recovery_rows[0]),
    )
    write_tsv(
        out / "04.Greedy_plan_recovery_overall_and_MAF.tsv",
        (row for row in all_recovery_rows if row["AC_threshold"] == "all"),
        list(all_recovery_rows[0]),
    )

    cumulative_rows: list[dict[str, Any]] = []
    for plan, samples in (("Greedy_Plan_A", plan_a_samples), ("Greedy_Plan_B", plan_b_samples)):
        covered = 0
        for step, sample in enumerate(samples, start=1):
            covered |= sample_bits[sample]
            for bin_label in ("All", *MAF_BIN_ORDER):
                mask = masks.get(("all", bin_label), 0)
                cumulative_rows.append({
                    "Plan": plan,
                    "Step": step,
                    "Selected_sample": sample,
                    "Wild_MAF_bin": bin_label,
                    "Recovered_missing_Wild_minor_alleles": (covered & mask).bit_count(),
                    "Missing_Wild_minor_alleles_total": mask.bit_count(),
                    "Recovery_rate": (covered & mask).bit_count() / mask.bit_count() if mask else 0.0,
                })
    write_tsv(out / "05.Greedy_cumulative_recovery.tsv", cumulative_rows, list(cumulative_rows[0]))

    region_rows: list[dict[str, Any]] = []
    for plan, samples in (("Greedy_Plan_A", plan_a_samples), ("Greedy_Plan_B", plan_b_samples)):
        counts = Counter(sample_region[sample] for sample in samples)
        quota = quotas[PLAN_A if plan.endswith("A") else PLAN_B]
        for region, expected_n in quota.items():
            region_rows.append({
                "Plan": plan,
                "Region": region,
                "Selected_sample_N": counts[region],
                "Quota": expected_n,
                "Quota_match": counts[region] == expected_n,
            })
    write_tsv(out / "06.Greedy_plan_region_composition.tsv", region_rows, list(region_rows[0]))

    cache_payload = {
        "format_version": 1,
        "missing_table": str(missing_path),
        "target_count": target_n,
        "wild_sample_order": wild_samples,
        "sample_bits_hex": {sample: hex(sample_bits[sample]) for sample in wild_samples},
        "called_sites": called_sites,
    }
    with gzip.open(out / "07.Wild_target_carrier_bitsets.json.gz", "wt", encoding="utf-8") as handle:
        json.dump(cache_payload, handle, sort_keys=True)

    lookup = {(row["Plan"], row["AC_threshold"], row["Wild_MAF_bin"]): row for row in all_recovery_rows}
    missing_common_target_n = masks.get(("all", "MAF_gt_0.10"), 0).bit_count()
    observed = {
        "Greedy_Plan_A_size": len(plan_a_samples),
        "Greedy_Plan_B_size": len(plan_b_samples),
        "Greedy_Plan_A_recovered_all": lookup[("Greedy_Plan_A", "all", "All")]["Recovered_missing_Wild_minor_alleles"],
        "Greedy_Plan_B_recovered_all": lookup[("Greedy_Plan_B", "all", "All")]["Recovered_missing_Wild_minor_alleles"],
        "Greedy_Plan_A_recovered_AC_ge_2": lookup[("Greedy_Plan_A", "AC_ge_2", "All")]["Recovered_missing_Wild_minor_alleles"],
        "Greedy_Plan_A_recovered_AC_ge_3": lookup[("Greedy_Plan_A", "AC_ge_3", "All")]["Recovered_missing_Wild_minor_alleles"],
    }
    write_key_value(out / "08.internal_check.tsv", [
        ("Missing_target_N", target_n),
        ("Wild_sample_N", len(wild_samples)),
        ("VCF_target_sites_matched", matched_n),
        ("Missing_targets_MAF_gt_0.10", missing_common_target_n),
        ("Plan_A_quota_total", sum(quotas[PLAN_A].values())),
        ("Plan_B_quota_total", sum(quotas[PLAN_B].values())),
        ("Plan_A_is_subset_of_Plan_B", set(plan_a_samples).issubset(plan_b_samples)),
        *observed.items(),
    ])
    print(f"[OK] Greedy Plan A/Plan B selection complete: {out}")


if __name__ == "__main__":
    main()
