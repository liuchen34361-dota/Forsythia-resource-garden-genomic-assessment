#!/usr/bin/env python3
"""Validate cross-step invariants and manuscript-level numerical anchors."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from analysis_utils import (
    default_config_path,
    load_config,
    outputs_root,
    read_sample_list,
    read_tsv,
    write_key_value,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(default_config_path()))
    return parser.parse_args()


def key_values(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    result: dict[str, str] = {}
    for row in rows:
        key = row.get("Item", "")
        if not key:
            raise ValueError(f"Missing Item value in {path}")
        if key in result:
            raise ValueError(f"Duplicated Item in {path}: {key}")
        result[key] = row.get("Value", row.get("Count", ""))
    return result


def as_float(value: Any) -> float:
    if isinstance(value, bool):
        return float(int(value))
    return float(str(value))


def unique_samples(path: Path) -> list[str]:
    rows = read_tsv(path)
    samples: list[str] = []
    seen: set[str] = set()
    for row in rows:
        sample = row.get("SampleID", "").strip()
        if not sample:
            raise ValueError(f"Blank SampleID in {path}")
        if sample in seen:
            continue
        seen.add(sample)
        samples.append(sample)
    return samples


def main() -> None:
    cfg = load_config(parse_args().config)
    root = outputs_root(cfg)
    out = root / "validation"
    out.mkdir(parents=True, exist_ok=True)

    s1 = key_values(root / "01_allele_representation" / "09.internal_check.tsv")
    s1_overall = read_tsv(root / "01_allele_representation" / "01.wild_minor_allele_retention_overall.tsv")
    s1_complement = read_tsv(root / "01_allele_representation" / "06.resource_garden_complementarity_summary.tsv")
    s2_occurrence = read_tsv(root / "02_regional_occurrence_and_gap" / "02.region_presence_by_AC_threshold.tsv")
    s3_overall = read_tsv(root / "03_sample_size_standardized_retention" / "01.sample_size_standardized_retention_overall.tsv")
    s4_rows = read_tsv(root / "04_callability_and_carrier_quality" / "03.Wild_carrier_DP_GQ_summary.tsv")
    s5 = key_values(root / "05_PIHAT_components_and_core" / "09.internal_check.tsv")
    s6 = key_values(root / "06_core_SNP_ALT_coverage" / "04.internal_check.tsv")
    s6_summary = read_tsv(root / "06_core_SNP_ALT_coverage" / "01.deredundancy_ALT_coverage_summary.tsv")
    s7 = key_values(root / "07_greedy_plan_selection" / "08.internal_check.tsv")
    s7_recovery = read_tsv(root / "07_greedy_plan_selection" / "03.Greedy_plan_recovery_by_AC_threshold_and_MAF.tsv")
    s8 = key_values(root / "08_region_matched_random_control" / "05.internal_check.tsv")
    s8_summary_rows = read_tsv(root / "08_region_matched_random_control" / "02.Greedy_vs_region_matched_random_summary.tsv")
    s9_rows = read_tsv(root / "09_selected_sample_QC" / "02.selected_Wild_QC_summary_by_plan.tsv")

    carrier_all = next(
        (
            row for row in s4_rows
            if row.get("AC_threshold") == "all" and row.get("Wild_MAF_bin") == "All"
        ),
        None,
    )
    if carrier_all is None:
        raise ValueError("Cannot identify all-site Wild-carrier DP/GQ summary row")
    qc_by_plan = {row["Plan"]: row for row in s9_rows}
    if "Greedy_Plan_B" not in qc_by_plan:
        raise ValueError("Selected-accession QC lacks Greedy_Plan_B")
    retention_by_target = {row["Target_group"]: row for row in s1_overall}
    complement_all = next((row for row in s1_complement if row.get("Wild_MAF_bin") == "All"), None)
    if complement_all is None:
        raise ValueError("Resource-garden complementarity summary lacks the All row")
    rarefaction_by_target = {row["Target_set"]: row for row in s3_overall}
    alt_all = next((row for row in s6_summary if row.get("Current_ALT_freq_bin") == "All"), None)
    if alt_all is None:
        raise ValueError("ALT coverage summary lacks the All row")
    recovery_lookup = {
        (row["Plan"], row["AC_threshold"], row["Wild_MAF_bin"]): row
        for row in s7_recovery
    }
    random_summary_lookup = {
        (row["Plan"], row["Wild_MAF_bin"]): row
        for row in s8_summary_rows
    }

    plan_a = unique_samples(root / "07_greedy_plan_selection" / "01.Greedy_Plan_A_Wild_supplements.tsv")
    plan_b = unique_samples(root / "07_greedy_plan_selection" / "02.Greedy_Plan_B_Wild_supplements.tsv")
    core = unique_samples(root / "05_PIHAT_components_and_core" / "06.deredundant_exsitu_core.tsv")
    low = unique_samples(root / "05_PIHAT_components_and_core" / "07.lower_priority_redundant_exsitu.tsv")

    observed: dict[str, float | int] = {
        "total_samples": int(float(s5["Total_samples"])),
        "wild_sample_n": len(read_sample_list(Path(cfg["_project_root"]) / cfg["inputs"]["wild_samples"])),
        "ex_situ_total": len(read_sample_list(Path(cfg["_project_root"]) / cfg["inputs"]["exsitu_samples"])),
        "pingshun_sample_n": len(read_sample_list(Path(cfg["_project_root"]) / cfg["inputs"]["pingshun_samples"])),
        "tongguan_sample_n": len(read_sample_list(Path(cfg["_project_root"]) / cfg["inputs"]["tongguan_samples"])),
        "zhendong_sample_n": len(read_sample_list(Path(cfg["_project_root"]) / cfg["inputs"]["zhendong_samples"])),
        "wild_minor_total": int(float(s1["Total_wild_SNP_minor_alleles"])),
        "combined_ex_situ_retained": int(float(s1["Combined_exsitu_retained"])),
        "combined_ex_situ_missing": int(float(s1["Combined_exsitu_missing"])),
        "combined_ex_situ_retention_rate": int(float(s1["Combined_exsitu_retained"]))
        / int(float(s1["Total_wild_SNP_minor_alleles"])),
        "pingshun_retained": int(float(s1["Pingshun_retained"])),
        "tongguan_retained": int(float(s1["Tongguan_retained"])),
        "zhendong_retained": int(float(s1["Zhendong_retained"])),
        "missing_ac_ge_2": int(float(s1["AC_ge_2_missing"])),
        "missing_ac_ge_3": int(float(s1["AC_ge_3_missing"])),
        "wild_carrier_genotype_n": int(float(carrier_all["Carrier_genotype_N"])),
        "wild_carrier_dp_median": as_float(carrier_all["DP_Median"]),
        "wild_carrier_gq_median": as_float(carrier_all["GQ_Median"]),
        "de_redundant_core_n": len(core),
        "lower_priority_redundant_n": len(low),
        "pingshun_core_n": int(float(s5["Pingshun_core_size"])),
        "tongguan_core_n": int(float(s5["Tongguan_core_size"])),
        "zhendong_core_n": int(float(s5["Zhendong_core_size"])),
        "reported_alt_states_total": int(float(s6["SNP_state_denominator"])),
        "reported_alt_states_retained": int(float(s6["SNP_state_numerator"])),
        "reported_alt_states_lost": int(float(s6["SNP_state_denominator"]))
        - int(float(s6["SNP_state_numerator"])),
        "reported_alt_state_retention_rate": int(float(s6["SNP_state_numerator"]))
        / int(float(s6["SNP_state_denominator"])),
        "alt_coverage_vcf_site_n": int(float(s6["VCF_sites_processed"])),
        "true_alt_copy_total": int(float(alt_all["Current_ex_situ_ALT_copies"])),
        "true_alt_copy_covered": int(float(alt_all["Current_ALT_copies_covered_by_core_presence"])),
        "true_alt_copy_uncovered": int(float(alt_all["Current_ALT_copies_not_covered_by_core_presence"])),
        "true_alt_copy_weighted_coverage_rate": as_float(alt_all["True_ALT_copy_weighted_coverage_rate"]),
        "random_seed": int(float(s8["Random_seed"])),
        "random_permutations": int(float(s8["Random_permutations"])),
        "plan_a_random_gain_all_pp_rounded2": round(as_float(random_summary_lookup[("Greedy_Plan_A", "All")]["Greedy_minus_random_q975_percentage_points"]), 2),
        "plan_a_random_gain_very_rare_pp_rounded2": round(as_float(random_summary_lookup[("Greedy_Plan_A", "MAF_le_0.01")]["Greedy_minus_random_q975_percentage_points"]), 2),
        "plan_a_random_gain_low_frequency_pp_rounded2": round(as_float(random_summary_lookup[("Greedy_Plan_A", "0.01_lt_MAF_le_0.05")]["Greedy_minus_random_q975_percentage_points"]), 2),
        "plan_b_random_gain_all_pp_rounded2": round(as_float(random_summary_lookup[("Greedy_Plan_B", "All")]["Greedy_minus_random_q975_percentage_points"]), 2),
        "plan_b_random_gain_very_rare_pp_rounded2": round(as_float(random_summary_lookup[("Greedy_Plan_B", "MAF_le_0.01")]["Greedy_minus_random_q975_percentage_points"]), 2),
        "plan_a_wild_n": len(plan_a),
        "priority_2_wild_n": len(plan_b) - len(plan_a),
        "plan_b_wild_n": len(plan_b),
        "plan_a_recovered": int(float(s7["Greedy_Plan_A_recovered_all"])),
        "plan_b_recovered": int(float(s7["Greedy_Plan_B_recovered_all"])),
        "plan_a_recovery_rate": int(float(s7["Greedy_Plan_A_recovered_all"]))
        / int(float(s1["Combined_exsitu_missing"])),
        "plan_b_recovery_rate": int(float(s7["Greedy_Plan_B_recovered_all"]))
        / int(float(s1["Combined_exsitu_missing"])),
        "selected_qc_roh_gt2mb_n": int(float(qc_by_plan["Greedy_Plan_B"]["Selected_with_ROH_gt_2Mb_N"])),
        "selected_qc_duplicate_component_n": int(float(qc_by_plan["Greedy_Plan_B"]["Repeated_PIHAT_component_N"])),
    }

    invariant_rows: list[dict[str, Any]] = []

    def invariant(name: str, passed: bool, detail: str) -> None:
        invariant_rows.append({"Invariant": name, "Status": "PASS" if passed else "FAIL", "Detail": detail})

    invariant(
        "wild_retained_plus_missing_equals_total",
        observed["combined_ex_situ_retained"] + observed["combined_ex_situ_missing"] == observed["wild_minor_total"],
        f"{observed['combined_ex_situ_retained']} + {observed['combined_ex_situ_missing']} vs {observed['wild_minor_total']}",
    )
    invariant(
        "core_and_low_partition_current_ex_situ",
        set(core).isdisjoint(low) and len(core) + len(low) == observed["ex_situ_total"],
        f"core={len(core)}, low={len(low)}, ex_situ={observed['ex_situ_total']}",
    )
    invariant(
        "plan_A_is_subset_of_plan_B",
        set(plan_a).issubset(plan_b),
        f"Plan_A={len(plan_a)}, Plan_B={len(plan_b)}",
    )
    invariant(
        "ALT_state_arithmetic",
        observed["reported_alt_states_retained"] + observed["reported_alt_states_lost"]
        == observed["reported_alt_states_total"],
        f"retained={observed['reported_alt_states_retained']}, lost={observed['reported_alt_states_lost']}, total={observed['reported_alt_states_total']}",
    )
    invariant(
        "complementarity_still_missing_equals_combined_exsitu_missing",
        int(float(complement_all["Still_missing_from_all_ex_situ"]))
        == observed["combined_ex_situ_missing"],
        f"complement={complement_all['Still_missing_from_all_ex_situ']}, "
        f"combined_missing={observed['combined_ex_situ_missing']}",
    )
    regional_all_total = sum(
        int(float(row["N_missing_Wild_minor_alleles"]))
        for row in s2_occurrence if row.get("AC_threshold") == "all"
    )
    invariant(
        "regional_occurrence_partitions_missing_targets",
        regional_all_total == observed["combined_ex_situ_missing"],
        f"regional_total={regional_all_total}, missing={observed['combined_ex_situ_missing']}",
    )
    garden_target_map = {
        "PS_RG_observed": "Pingshun",
        "TG_RG_observed": "Tongguan",
        "ZD_RG_observed": "Zhendong",
    }
    for rarefaction_name, target_name in garden_target_map.items():
        rarefied_observed = float(rarefaction_by_target[rarefaction_name]["Retained_Wild_SNP_minor_alleles"])
        audit_observed = float(retention_by_target[target_name]["Retained_Wild_minor_alleles"])
        invariant(
            f"rarefaction_observed_matches_allele_audit_{target_name}",
            abs(rarefied_observed - audit_observed) < 1e-8,
            f"rarefaction={rarefied_observed}, audit={audit_observed}",
        )
    true_copy_total = int(float(alt_all["Current_ex_situ_ALT_copies"]))
    true_copy_covered = int(float(alt_all["Current_ALT_copies_covered_by_core_presence"]))
    true_copy_uncovered = int(float(alt_all["Current_ALT_copies_not_covered_by_core_presence"]))
    invariant(
        "true_ALT_copy_weighted_arithmetic",
        true_copy_covered + true_copy_uncovered == true_copy_total,
        f"covered={true_copy_covered}, uncovered={true_copy_uncovered}, total={true_copy_total}",
    )
    invariant(
        "greedy_full_target_recovery_matches_internal_counts",
        int(float(recovery_lookup[("Greedy_Plan_A", "all", "All")]["Recovered_missing_Wild_minor_alleles"]))
        == observed["plan_a_recovered"]
        and int(float(recovery_lookup[("Greedy_Plan_B", "all", "All")]["Recovered_missing_Wild_minor_alleles"]))
        == observed["plan_b_recovered"],
        f"PlanA={observed['plan_a_recovered']}, PlanB={observed['plan_b_recovered']}",
    )
    invariant(
        "random_distribution_row_count",
        int(float(s8["Distribution_rows"]))
        == int(float(s8["Random_permutations"])) * 2 * 5,
        f"rows={s8['Distribution_rows']}, permutations={s8['Random_permutations']}",
    )
    invariant(
        "selected_QC_no_missing_join_silent_fallback",
        len(plan_b) == int(float(qc_by_plan["Greedy_Plan_B"]["Selected_sample_N"])),
        f"Plan_B={len(plan_b)}, QC={qc_by_plan['Greedy_Plan_B']['Selected_sample_N']}",
    )

    portfolios = {r["Portfolio"]: r for r in read_tsv(root / "10_final_wild_portfolios" / "01.final_portfolios.tsv")}
    invariant("direct_current_matches_ACAN_baseline",
        int(portfolios["Current"]["Target_total"]) == observed["wild_minor_total"]
        and int(portfolios["Current"]["Missing"]) == observed["combined_ex_situ_missing"],
        "Direct GT universe and Current missing count must agree with Step 01")
    invariant("direct_core_size_matches_components",
        int(portfolios["Core"]["Sample_N"]) == observed["de_redundant_core_n"],
        "Step 10 core size must agree with Step 05")

    expected = cfg.get("expected_results", {}) or {}
    validation_rows: list[dict[str, Any]] = []
    for item, expected_value in expected.items():
        if item not in observed:
            validation_rows.append({
                "Item": item,
                "Observed": "NA",
                "Expected": expected_value,
                "Absolute_difference": "NA",
                "Tolerance": "NA",
                "Status": "NOT_IMPLEMENTED",
            })
            continue
        observed_value = observed[item]
        tolerance = 1e-8 if isinstance(expected_value, float) and not float(expected_value).is_integer() else 0.0
        difference = abs(float(observed_value) - float(expected_value))
        validation_rows.append({
            "Item": item,
            "Observed": observed_value,
            "Expected": expected_value,
            "Absolute_difference": difference,
            "Tolerance": tolerance,
            "Status": "PASS" if difference <= tolerance else "FAIL",
        })

    write_tsv(
        out / "01.manuscript_anchor_validation.tsv",
        validation_rows,
        ["Item", "Observed", "Expected", "Absolute_difference", "Tolerance", "Status"],
    )
    write_tsv(out / "02.cross_step_invariants.tsv", invariant_rows, ["Invariant", "Status", "Detail"])
    write_key_value(out / "03.observed_key_results.tsv", observed.items())

    failures = [row for row in validation_rows if row["Status"] in {"FAIL", "NOT_IMPLEMENTED"}]
    failures += [row for row in invariant_rows if row["Status"] == "FAIL"]
    if failures:
        preview = "\n".join(str(row) for row in failures[:20])
        raise AssertionError(f"Final workflow validation failed:\n{preview}")
    print(f"[OK] Final cross-step and manuscript-anchor validation complete: {out}")


if __name__ == "__main__":
    main()
