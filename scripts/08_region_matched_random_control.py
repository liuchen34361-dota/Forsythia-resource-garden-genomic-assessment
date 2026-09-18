#!/usr/bin/env python3
"""Region-matched random control for Greedy Plan A and Plan B.

Historical-replay compatibility
-------------------------------
The manuscript's random-control results were generated with Python's standard
``random`` module after seeding once with 20260620.  For each plan, the
historical workflow first consumed 1,000 all-Wild draws and then sampled the
positive-quota regions in a specific insertion order.  Because successive
``random.sample`` calls consume the pseudorandom stream, changing that region
order changes every downstream portfolio even when the seed, pools and quotas
are unchanged.

This cleaned implementation therefore treats ``Random_sampling_order`` in
``config/regional_quotas.tsv`` as part of the reproducibility specification.
Configured quotas remain the authoritative composition check.  Zero-quota
regions are recorded for completeness but do not consume random numbers.
"""
from __future__ import annotations

import argparse
import gzip
import json
import random
import statistics
import sys
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path
from typing import Any

from analysis_utils import (
    prepare_step_output,
    MAF_BIN_ORDER,
    canonical_region,
    default_config_path,
    input_path,
    load_config,
    load_missing_targets,
    package_path,
    quantile,
    read_sample_list,
    read_tsv,
    step_output_dir,
    write_key_value,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--missing-table", default=None)
    return parser.parse_args()


def load_quotas(path: Path) -> tuple[dict[str, OrderedDict[str, int]], str]:
    """Load regional quotas while preserving an explicit RNG sampling order.

    ``Random_sampling_order`` is preferred.  A file without that column remains
    usable for generic/synthetic datasets; in that case, row order is used.
    """
    raw_by_plan: dict[str, list[tuple[int, int, str, int]]] = defaultdict(list)
    has_explicit_order = False

    for row_index, row in enumerate(read_tsv(path), start=1):
        plan = row["Plan"]
        region = canonical_region(row["Region"])
        sample_n = int(row.get("Sample_N", row.get("Quota", "0")))
        if sample_n < 0:
            raise ValueError(f"Negative quota for {plan}, {region}: {sample_n}")

        order_text = str(row.get("Random_sampling_order", "")).strip()
        if order_text:
            order = int(order_text)
            if order <= 0:
                raise ValueError(
                    f"Random_sampling_order must be positive for {plan}, {region}: {order}"
                )
            has_explicit_order = True
        else:
            order = row_index

        raw_by_plan[plan].append((order, row_index, region, sample_n))

    result: dict[str, OrderedDict[str, int]] = {}
    for plan, records in raw_by_plan.items():
        seen_regions: set[str] = set()
        explicit_orders = [order for order, _, _, _ in records]
        if len(explicit_orders) != len(set(explicit_orders)):
            raise ValueError(f"Duplicated Random_sampling_order within {plan}")

        ordered: OrderedDict[str, int] = OrderedDict()
        for _, _, region, sample_n in sorted(records, key=lambda item: (item[0], item[1])):
            if region in seen_regions:
                raise ValueError(f"Duplicate quota for {plan}, {region}")
            seen_regions.add(region)
            ordered[region] = sample_n
        result[plan] = ordered

    source = (
        "Random_sampling_order column in config/regional_quotas.tsv"
        if has_explicit_order
        else "row order in regional_quotas.tsv (no Random_sampling_order column)"
    )
    return result, source


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    out = prepare_step_output(cfg, "08_region_matched_random_control")
    greedy_dir = step_output_dir(cfg, "07_greedy_plan_selection")
    missing_path = Path(args.missing_table) if args.missing_table else (
        step_output_dir(cfg, "01_allele_representation")
        / "05.ex_situ_missing_wild_minor_alleles.tsv.gz"
    )
    targets, _ = load_missing_targets(missing_path)

    masks = {"All": 0, **{bin_label: 0 for bin_label in MAF_BIN_ORDER}}
    for info in targets.values():
        bit = 1 << int(info["idx"])
        masks["All"] |= bit
        masks[str(info["maf_bin"])] |= bit

    cache_path = greedy_dir / "07.Wild_target_carrier_bitsets.json.gz"
    with gzip.open(cache_path, "rt", encoding="utf-8") as handle:
        cache = json.load(handle)

    wild_samples = cache["wild_sample_order"]
    configured_wild_samples = read_sample_list(input_path(cfg, "wild_samples"))
    if wild_samples != configured_wild_samples:
        raise ValueError(
            "Carrier cache Wild sample order does not match the configured sample list"
        )

    sample_bits = {
        sample: int(value, 16)
        for sample, value in cache["sample_bits_hex"].items()
    }
    if set(sample_bits) != set(wild_samples):
        raise ValueError("Carrier cache sample keys do not match its Wild sample order")
    if int(cache["target_count"]) != len(targets):
        raise ValueError("Carrier cache target count differs from the missing table")

    metadata_rows = read_tsv(input_path(cfg, "metadata"))
    metadata: dict[str, dict[str, str]] = {}
    for row in metadata_rows:
        sample = row.get("SampleID", "").strip()
        if not sample:
            continue
        if sample in metadata:
            raise ValueError(f"Duplicated SampleID in metadata: {sample}")
        metadata[sample] = row

    absent_metadata = [sample for sample in wild_samples if sample not in metadata]
    if absent_metadata:
        raise ValueError(f"Wild samples missing from metadata: {absent_metadata[:10]}")

    sample_region = {
        sample: canonical_region(metadata[sample].get("Region"))
        for sample in wild_samples
    }
    region_to_samples: dict[str, list[str]] = defaultdict(list)
    for sample in wild_samples:
        region_to_samples[sample_region[sample]].append(sample)

    def plan_samples(path: Path) -> list[str]:
        samples: list[str] = []
        seen: set[str] = set()
        for row in read_tsv(path):
            sample = row["SampleID"]
            if sample not in seen:
                samples.append(sample)
                seen.add(sample)
        return samples

    plans = OrderedDict(
        [
            (
                "Greedy_Plan_A",
                plan_samples(greedy_dir / "01.Greedy_Plan_A_Wild_supplements.tsv"),
            ),
            (
                "Greedy_Plan_B",
                plan_samples(greedy_dir / "02.Greedy_Plan_B_Wild_supplements.tsv"),
            ),
        ]
    )

    quotas, region_order_source = load_quotas(
        package_path(cfg, cfg["package_files"]["regional_quotas"])
    )
    quota_key = {
        "Greedy_Plan_A": "Plan_A_conservative",
        "Greedy_Plan_B": "Plan_B_expanded",
    }

    parameters = cfg.get("parameters", {})
    random_parameters = parameters.get("random_control", {})
    seed = int(random_parameters.get("seed", parameters.get("random_seed", 20260620)))
    permutations = int(
        random_parameters.get(
            "n_permutations",
            parameters.get("random_permutations", 1000),
        )
    )
    burn_legacy = bool(
        random_parameters.get(
            "legacy_rng_compatibility",
            parameters.get("legacy_random_stream_burn_all_wild", True),
        )
    )
    if permutations <= 0:
        raise ValueError("random_permutations must be positive")

    # A dedicated Random instance has the same stream as random.seed(seed)
    # while avoiding accidental interference from unrelated code.
    rng = random.Random(seed)

    def merged_bits(samples: list[str]) -> int:
        bits = 0
        for sample in samples:
            bits |= sample_bits[sample]
        return bits

    def recovery(bits: int, bin_label: str) -> tuple[int, int, float]:
        mask = masks[bin_label]
        total = mask.bit_count()
        recovered = (bits & mask).bit_count()
        return recovered, total, recovered / total if total else 0.0

    observed: dict[tuple[str, str], tuple[int, int, float]] = {}
    for plan, samples in plans.items():
        bits = merged_bits(samples)
        for bin_label in ("All", *MAF_BIN_ORDER):
            observed[(plan, bin_label)] = recovery(bits, bin_label)

    distribution_rows: list[dict[str, Any]] = []
    rates: dict[tuple[str, str], list[float]] = defaultdict(list)
    sampling_order_rows: list[dict[str, Any]] = []

    for plan, selected in plans.items():
        plan_quota = quotas[quota_key[plan]]
        selected_counts = Counter(sample_region[sample] for sample in selected)

        unknown_regions = set(selected_counts) - set(plan_quota)
        if unknown_regions:
            raise ValueError(
                f"{plan} contains regions absent from configured quotas: "
                f"{sorted(unknown_regions)}"
            )

        mismatches = {
            region: (selected_counts.get(region, 0), quota)
            for region, quota in plan_quota.items()
            if selected_counts.get(region, 0) != quota
        }
        if mismatches:
            raise ValueError(
                f"{plan} region composition does not match configured quotas: "
                f"{mismatches}"
            )

        for order_index, (region, sample_n) in enumerate(plan_quota.items(), start=1):
            available_n = len(region_to_samples.get(region, []))
            if sample_n > available_n:
                raise ValueError(
                    f"Cannot sample {sample_n} accessions from {region}; "
                    f"available={available_n}"
                )
            sampling_order_rows.append(
                {
                    "Plan": plan,
                    "Region_sampling_order": order_index,
                    "Region": region,
                    "Quota": sample_n,
                    "Available_Wild_samples": available_n,
                    "Included_in_random_sampling": sample_n > 0,
                }
            )

        # The historical script generated 1,000 all-Wild portfolios before the
        # region-matched block for each plan.  The recovery calculation itself
        # consumes no random numbers, so consuming the draws is sufficient.
        if burn_legacy:
            for _ in range(permutations):
                rng.sample(wild_samples, len(selected))

        for permutation in range(1, permutations + 1):
            picked: list[str] = []
            for region, n in plan_quota.items():
                if n <= 0:
                    continue
                picked.extend(rng.sample(region_to_samples[region], n))

            bits = merged_bits(picked)
            for bin_label in ("All", *MAF_BIN_ORDER):
                recovered, total, rate = recovery(bits, bin_label)
                distribution_rows.append(
                    {
                        "Plan": plan,
                        "Control_type": "random_region_matched",
                        "Permutation": permutation,
                        "Wild_supplement_N": len(selected),
                        "Wild_MAF_bin": bin_label,
                        "Recovered_missing_Wild_minor_alleles": recovered,
                        "Missing_Wild_minor_alleles_total": total,
                        "Recovery_rate": rate,
                    }
                )
                rates[(plan, bin_label)].append(rate)

    write_tsv(
        out / "01.region_matched_random_recovery_distribution.tsv",
        distribution_rows,
        list(distribution_rows[0]),
    )

    summary_rows: list[dict[str, Any]] = []
    for plan in plans:
        for bin_label in ("All", *MAF_BIN_ORDER):
            values = rates[(plan, bin_label)]
            obs_recovered, total, obs_rate = observed[(plan, bin_label)]
            q975 = quantile(values, 0.975)
            sd = statistics.stdev(values) if len(values) > 1 else 0.0
            summary_rows.append(
                {
                    "Plan": plan,
                    "Control_type": "random_region_matched",
                    "Wild_supplement_N": len(plans[plan]),
                    "Wild_MAF_bin": bin_label,
                    "Missing_Wild_minor_alleles_total": total,
                    "Greedy_recovered": obs_recovered,
                    "Greedy_recovery_rate": obs_rate,
                    "Random_N": len(values),
                    "Random_mean": statistics.mean(values),
                    "Random_sd": sd,
                    "Random_min": min(values),
                    "Random_q025": quantile(values, 0.025),
                    "Random_q05": quantile(values, 0.05),
                    "Random_q50": quantile(values, 0.50),
                    "Random_q95": quantile(values, 0.95),
                    "Random_q975": q975,
                    "Random_max": max(values),
                    "Greedy_minus_random_q975_percentage_points": (
                        obs_rate - q975
                    )
                    * 100.0,
                    "Empirical_p_random_ge_greedy": (
                        sum(value >= obs_rate for value in values) + 1
                    )
                    / (len(values) + 1),
                    "Greedy_percentile_vs_random": sum(
                        value <= obs_rate for value in values
                    )
                    / len(values),
                    "Greedy_z_score": (
                        (obs_rate - statistics.mean(values)) / sd if sd else 0.0
                    ),
                }
            )

    write_tsv(
        out / "02.Greedy_vs_region_matched_random_summary.tsv",
        summary_rows,
        list(summary_rows[0]),
    )
    write_tsv(
        out / "03.region_sampling_order_and_quotas.tsv",
        sampling_order_rows,
        [
            "Plan",
            "Region_sampling_order",
            "Region",
            "Quota",
            "Available_Wild_samples",
            "Included_in_random_sampling",
        ],
    )

    validation_rows: list[dict[str, Any]] = []
    expected_gains = parameters.get("expected_random_gain_pp", {})
    tolerance = float(parameters.get("random_gain_validation_tolerance_pp", 0.02))
    for row in summary_rows:
        expected = expected_gains.get(row["Plan"], {}).get(row["Wild_MAF_bin"])
        if expected is None:
            continue
        observed_gain = float(row["Greedy_minus_random_q975_percentage_points"])
        difference = abs(observed_gain - float(expected))
        status = "PASS" if difference <= tolerance else "FAIL"
        validation_rows.append(
            {
                "Plan": row["Plan"],
                "Wild_MAF_bin": row["Wild_MAF_bin"],
                "Observed_gain_pp": observed_gain,
                "Expected_manuscript_gain_pp": expected,
                "Tolerance_pp": tolerance,
                "Absolute_difference_pp": difference,
                "Status": status,
            }
        )

    if validation_rows:
        write_tsv(
            out / "04.manuscript_random_gain_validation.tsv",
            validation_rows,
            list(validation_rows[0]),
        )
        if parameters.get("validate_known_results", False) and any(
            row["Status"] == "FAIL" for row in validation_rows
        ):
            failed = [row for row in validation_rows if row["Status"] == "FAIL"]
            raise AssertionError(
                f"Random-control manuscript validation failed: {failed}"
            )

    write_key_value(
        out / "05.internal_check.tsv",
        [
            ("Random_seed", seed),
            ("Random_permutations", permutations),
            ("Legacy_all_wild_burn_enabled", burn_legacy),
            ("Python_version", sys.version.split()[0]),
            ("Region_order_source", region_order_source),
            ("Wild_sample_N", len(wild_samples)),
            ("Plan_A_sample_N", len(plans["Greedy_Plan_A"])),
            ("Plan_B_sample_N", len(plans["Greedy_Plan_B"])),
            ("Distribution_rows", len(distribution_rows)),
            ("Summary_rows", len(summary_rows)),
        ],
    )
    print(f"[OK] Region-matched random control complete: {out}")


if __name__ == "__main__":
    main()
