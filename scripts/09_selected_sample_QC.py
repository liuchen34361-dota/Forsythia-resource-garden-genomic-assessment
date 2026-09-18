#!/usr/bin/env python3
"""Quality and redundancy checks for Greedy-selected Wild accessions."""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from analysis_utils import (
    prepare_step_output,
    default_config_path,
    input_path,
    load_config,
    quantile,
    read_tsv,
    step_output_dir,
    write_key_value,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(default_config_path()))
    return parser.parse_args()


def first_column(row: dict[str, str], candidates: Iterable[str]) -> str | None:
    lower = {key.lower(): key for key in row}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def number(row: dict[str, str], candidates: Iterable[str]) -> float | None:
    column = first_column(row, candidates)
    if column is None:
        return None
    text = row.get(column, "").strip()
    if text in {"", ".", "NA", "NaN", "nan"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def roh_gt_2mb(row: dict[str, str]) -> tuple[bool, str]:
    count_2_5 = number(row, ["NSEG_2_5Mb", "N_ROH_2_5Mb", "ROH_2_5Mb_N"])
    count_gt5 = number(row, ["NSEG_gt5Mb", "NSEG_5Mb_plus", "N_ROH_gt5Mb"])
    if count_2_5 is not None or count_gt5 is not None:
        total = (count_2_5 or 0.0) + (count_gt5 or 0.0)
        return total > 0, "NSEG_2_5Mb_plus_NSEG_gt5Mb"
    max_kb = number(row, [
        "Max_ROH_KB", "Max_ROH_length_KB", "Maximum_ROH_KB",
        "Max_total_ROH_KB_observed", "Max_ROH_segment_KB",
    ])
    if max_kb is not None:
        return max_kb > 2000.0, "maximum_ROH_length_KB"
    max_bp = number(row, ["Max_ROH_BP", "Max_ROH_length_bp", "Maximum_ROH_BP"])
    if max_bp is not None:
        return max_bp > 2_000_000.0, "maximum_ROH_length_BP"
    raise ValueError("ROH summary lacks recognized >2 Mb segment columns")


def unique_plan_samples(path: Path) -> list[str]:
    samples: list[str] = []
    seen: set[str] = set()
    for row in read_tsv(path):
        sample = row["SampleID"]
        if sample not in seen:
            samples.append(sample)
            seen.add(sample)
    return samples


def main() -> None:
    cfg = load_config(parse_args().config)
    out = prepare_step_output(cfg, "09_selected_sample_QC")
    greedy = step_output_dir(cfg, "07_greedy_plan_selection")
    pihat = step_output_dir(cfg, "05_PIHAT_components_and_core")

    plans = {
        "Greedy_Plan_A": unique_plan_samples(greedy / "01.Greedy_Plan_A_Wild_supplements.tsv"),
        "Greedy_Plan_B": unique_plan_samples(greedy / "02.Greedy_Plan_B_Wild_supplements.tsv"),
    }
    call_rows = read_tsv(greedy / "00.Wild_target_site_call_rate.tsv")
    call_rate: dict[str, float] = {}
    for row in call_rows:
        sample = row["SampleID"]
        if sample in call_rate:
            raise ValueError(f"Duplicated SampleID in target-site call-rate table: {sample}")
        call_rate[sample] = float(row["Target_site_call_rate"])
    wild_rates = sorted(call_rate.values())
    if not wild_rates:
        raise ValueError("Wild target-site call-rate table is empty")

    status_rows = read_tsv(pihat / "08.operational_global_sample_status.tsv")
    component: dict[str, dict[str, str]] = {}
    for row in status_rows:
        sample = row["SampleID"]
        if sample in component:
            raise ValueError(f"Duplicated SampleID in operational PI_HAT status: {sample}")
        component[sample] = row
    roh_rows = read_tsv(input_path(cfg, "roh_summary"))
    roh: dict[str, dict[str, str]] = {}
    for row in roh_rows:
        sample = row.get("SampleID", "").strip()
        if not sample:
            continue
        if sample in roh:
            raise ValueError(f"Duplicated SampleID in ROH summary: {sample}")
        roh[sample] = row

    if not set(plans["Greedy_Plan_A"]).issubset(plans["Greedy_Plan_B"]):
        raise ValueError("Greedy Plan A is not a subset of Greedy Plan B")
    selected_all = set(plans["Greedy_Plan_B"])
    absent_call = sorted(selected_all - set(call_rate))
    absent_component = sorted(selected_all - set(component))
    absent_roh = sorted(selected_all - set(roh))
    if absent_call or absent_component or absent_roh:
        raise ValueError(
            f"Missing selected-sample annotations: call={absent_call[:5]}, "
            f"component={absent_component[:5]}, ROH={absent_roh[:5]}"
        )
    operational = float(cfg.get("parameters", {}).get("operational_pihat", 0.80))
    wrong_threshold = [
        sample for sample in selected_all
        if abs(float(component[sample]["PI_HAT_min"]) - operational) > 1e-12
    ]
    if wrong_threshold:
        raise ValueError(f"Selected samples were joined to the wrong PI_HAT threshold: {wrong_threshold[:10]}")

    threshold = float(cfg.get("parameters", {}).get("selected_low_call_rate_threshold", 0.60))
    sample_rows: list[dict[str, Any]] = []
    duplicate_rows: list[dict[str, Any]] = []
    plan_rows: list[dict[str, Any]] = []

    for plan, samples in plans.items():
        component_counts = Counter(component[sample]["Global_component_ID"] for sample in samples)
        duplicate_components = {cid: n for cid, n in component_counts.items() if n > 1}
        for cid, n in sorted(duplicate_components.items()):
            members = [sample for sample in samples if component[sample]["Global_component_ID"] == cid]
            duplicate_rows.append({
                "Plan": plan,
                "Global_component_ID": cid,
                "Selected_sample_N": n,
                "Selected_samples": ";".join(members),
            })

        long_roh_n = 0
        low_call_n = 0
        selected_rates: list[float] = []
        roh_method: str | None = None
        for sample in samples:
            has_long, method = roh_gt_2mb(roh[sample])
            if roh_method is None:
                roh_method = method
            long_roh_n += has_long
            rate = call_rate[sample]
            low_call_n += rate < threshold
            selected_rates.append(rate)
            percentile = sum(value <= rate for value in wild_rates) / len(wild_rates)
            sample_rows.append({
                "Plan": plan,
                "SampleID": sample,
                "Target_site_call_rate": rate,
                "Call_rate_percentile_within_all_Wild": percentile,
                "Below_selected_call_rate_threshold": rate < threshold,
                "Global_component_ID_PIHAT_ge_0.80": component[sample]["Global_component_ID"],
                "Global_component_size_PIHAT_ge_0.80": component[sample]["Global_component_size"],
                "ROH_gt_2Mb": has_long,
                "ROH_gt_2Mb_detection_method": method,
            })

        plan_rows.append({
            "Plan": plan,
            "Selected_sample_N": len(samples),
            "Call_rate_min": min(selected_rates),
            "Call_rate_p05": quantile(selected_rates, 0.05),
            "Call_rate_median": quantile(selected_rates, 0.50),
            "Call_rate_mean": sum(selected_rates) / len(selected_rates),
            "Call_rate_max": max(selected_rates),
            "Below_call_rate_threshold_N": low_call_n,
            "Call_rate_threshold": threshold,
            "Selected_with_ROH_gt_2Mb_N": long_roh_n,
            "Repeated_PIHAT_component_N": len(duplicate_components),
            "Repeated_selected_accession_N": sum(duplicate_components.values()),
            "ROH_detection_method": roh_method,
        })

    write_tsv(out / "01.selected_Wild_accession_QC.tsv", sample_rows, list(sample_rows[0]))
    write_tsv(out / "02.selected_Wild_QC_summary_by_plan.tsv", plan_rows, list(plan_rows[0]))
    duplicate_fields = ["Plan", "Global_component_ID", "Selected_sample_N", "Selected_samples"]
    write_tsv(out / "03.repeated_PIHAT_component_check.tsv", duplicate_rows, duplicate_fields)

    failed_long = sum(int(row["Selected_with_ROH_gt_2Mb_N"]) for row in plan_rows)
    failed_duplicates = sum(int(row["Repeated_PIHAT_component_N"]) for row in plan_rows)
    plan_b_summary = next(row for row in plan_rows if row["Plan"] == "Greedy_Plan_B")
    write_key_value(out / "04.internal_check.tsv", [
        ("All_Wild_call_rate_N", len(wild_rates)),
        ("All_Wild_call_rate_min", min(wild_rates)),
        ("All_Wild_call_rate_median", quantile(wild_rates, 0.50)),
        ("All_Wild_call_rate_max", max(wild_rates)),
        ("Plan_A_sample_N", len(plans["Greedy_Plan_A"])),
        ("Plan_B_sample_N", len(plans["Greedy_Plan_B"])),
        ("Selected_with_ROH_gt_2Mb_across_plan_checks", failed_long),
        ("Repeated_PIHAT_components_across_plan_checks", failed_duplicates),
        ("Plan_B_selected_with_ROH_gt_2Mb", plan_b_summary["Selected_with_ROH_gt_2Mb_N"]),
        ("Plan_B_repeated_PIHAT_components", plan_b_summary["Repeated_PIHAT_component_N"]),
        ("Plan_A_is_subset_of_Plan_B", True),
    ])
    if cfg.get("parameters", {}).get("validate_known_results", False):
        if failed_long != 0:
            raise AssertionError("Selected-accession QC found ROH >2 Mb")
        if failed_duplicates != 0:
            raise AssertionError("Selected-accession QC found repeated PI_HAT components")
    print(f"[OK] Selected-accession QC complete: {out}")


if __name__ == "__main__":
    main()
