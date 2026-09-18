#!/usr/bin/env python3
"""Run the complete cleaned workflow on the packaged synthetic dataset."""
from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def key_values(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    return {row["Item"]: row["Value"] for row in rows}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    package_root = Path(__file__).resolve().parents[1]
    source_project = package_root / "tests" / "synthetic_project"
    with tempfile.TemporaryDirectory(prefix="forsythia_audit_smoke_") as tmp:
        project = Path(tmp) / "synthetic_project"
        shutil.copytree(source_project, project, ignore=shutil.ignore_patterns("results", "__pycache__", "*.pyc"))
        config = project / "config" / "analysis_config.json"
        subprocess.run(
            ["bash", str(package_root / "run_pipeline.sh"), str(config), str(project)],
            check=True,
            cwd=project,
        )
        results = project / "results"

        step1 = key_values(results / "01_allele_representation" / "09.internal_check.tsv")
        require(step1["Total_wild_SNP_minor_alleles"] == "6", "Synthetic Wild minor total should be 6")
        require(step1["Combined_exsitu_retained"] == "4", "Synthetic ex situ retained should be 4")
        require(step1["Combined_exsitu_missing"] == "2", "Synthetic ex situ missing should be 2")
        require(step1["AC_ge_2_missing"] == "2", "Synthetic AC>=2 missing should be 2")
        require(step1["AC_ge_3_missing"] == "2", "Synthetic AC>=3 missing should be 2")

        complement = {
            row["Wild_MAF_bin"]: row
            for row in read_tsv(results / "01_allele_representation" / "06.resource_garden_complementarity_summary.tsv")
        }
        require(complement["All"]["Still_missing_from_all_ex_situ"] == "2", "Complementarity target mismatch")
        require(complement["All"]["Recovered_by_TG_RG_or_ZD_RG"] == "3", "Complementarity recovery mismatch")

        pihat = key_values(results / "05_PIHAT_components_and_core" / "09.internal_check.tsv")
        require(pihat["Exsitu_core_size"] == "3", "Synthetic core size should be 3")
        require(pihat["Exsitu_low_priority_size"] == "1", "Synthetic lower-priority size should be 1")

        alt = {
            row["Current_ALT_freq_bin"]: row
            for row in read_tsv(results / "06_core_SNP_ALT_coverage" / "01.deredundancy_ALT_coverage_summary.tsv")
        }["All"]
        require(alt["Current_ex_situ_observed_ALT_states"] == "4", "Synthetic ALT-state denominator should be 4")
        require(alt["ALT_states_retained_in_core"] == "3", "Synthetic ALT-state numerator should be 3")
        require(alt["Current_ex_situ_ALT_copies"] == "9", "Synthetic ALT-copy denominator should be 9")
        require(alt["Current_ALT_copies_covered_by_core_presence"] == "8", "Synthetic copy-weighted numerator should be 8")

        greedy = key_values(results / "07_greedy_plan_selection" / "08.internal_check.tsv")
        require(greedy["Greedy_Plan_A_size"] == "2", "Synthetic Plan A size should be 2")
        require(greedy["Greedy_Plan_B_size"] == "4", "Synthetic Plan B size should be 4")
        require(greedy["Plan_A_is_subset_of_Plan_B"] == "True", "Synthetic plans must be nested")

        random_check = key_values(results / "08_region_matched_random_control" / "05.internal_check.tsv")
        require(random_check["Distribution_rows"] == "200", "Synthetic random distribution should have 200 rows")

        qc = {
            row["Plan"]: row
            for row in read_tsv(results / "09_selected_sample_QC" / "02.selected_Wild_QC_summary_by_plan.tsv")
        }
        require(qc["Greedy_Plan_B"]["Selected_with_ROH_gt_2Mb_N"] == "0", "Synthetic Plan B should have no >2 Mb ROH")
        require(qc["Greedy_Plan_B"]["Repeated_PIHAT_component_N"] == "0", "Synthetic Plan B should have no repeated PI_HAT component")

        invariants = read_tsv(results / "validation" / "02.cross_step_invariants.tsv")
        failed = [row for row in invariants if row["Status"] != "PASS"]
        require(not failed, f"Synthetic cross-step invariants failed: {failed}")

    print("[PASS] Synthetic end-to-end smoke test completed successfully.")


if __name__ == "__main__":
    main()
