#!/usr/bin/env python3
"""Resource-garden SNP ALT-state presence retention and copy-weighted coverage.
Only unequal single-base A/C/G/T REF/ALT states qualify. Other ALT states
are excluded individually; eligible SNP states in the same record remain.
Weights are Current AC, never Core AC. Called partial GT alleles count.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from analysis_utils import (
    prepare_step_output,
    default_config_path,
    input_path,
    load_config,
    open_text,
    read_sample_list,
    read_vcf_header_samples,
    read_tsv,
    step_output_dir,
    write_key_value,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(default_config_path()))
    return parser.parse_args()


def parse_gt_alt_coverage(gt: str, allow_partial: bool) -> list[int]:
    if not gt or gt in {".", "./.", ".|."}:
        return []
    alleles: list[int] = []
    for token in gt.replace("|", "/").split("/"):
        if token == ".":
            if allow_partial:
                continue
            return []
        try:
            alleles.append(int(token))
        except ValueError:
            return []
    return alleles


def alt_frequency_bin(ac: int, an: int) -> str:
    if ac <= 0 or an <= 0:
        return "ALT_absent"
    af = ac / an
    if ac == 1:
        return "AC_eq_1"
    if ac <= 5:
        return "AC_2_to_5"
    if af <= 0.01:
        return "AF_le_0.01"
    if af <= 0.05:
        return "0.01_lt_AF_le_0.05"
    if af <= 0.10:
        return "0.05_lt_AF_le_0.10"
    return "AF_gt_0.10"


def main() -> None:
    cfg = load_config(parse_args().config)
    out = prepare_step_output(cfg, "06_core_SNP_ALT_coverage")
    current_samples = read_sample_list(input_path(cfg, "exsitu_samples"))
    current_set = set(current_samples)

    core_path = step_output_dir(cfg, "05_PIHAT_components_and_core") / "06.deredundant_exsitu_core.tsv"
    core_rows = read_tsv(core_path)
    core_samples = [row["SampleID"] for row in core_rows]
    core_set = set(core_samples)
    if len(core_set) != len(core_samples):
        raise ValueError("Duplicated SampleID in de-redundant core")
    if not core_set:
        raise ValueError("De-redundant core is empty")
    if not core_set.issubset(current_set):
        extra = sorted(core_set - current_set)
        raise ValueError(f"Core contains samples absent from current ex situ: {extra[:10]}")
    low_set = current_set - core_set
    low_path = step_output_dir(cfg, "05_PIHAT_components_and_core") / "07.lower_priority_redundant_exsitu.tsv"
    low_rows = read_tsv(low_path)
    low_output_samples = [row["SampleID"] for row in low_rows]
    if len(low_output_samples) != len(set(low_output_samples)):
        raise ValueError("Duplicated SampleID in lower-priority ex situ output")
    if set(low_output_samples) != low_set:
        missing = sorted(low_set - set(low_output_samples))
        extra = sorted(set(low_output_samples) - low_set)
        raise ValueError(
            "Core/lower-priority partition disagrees with the configured current ex situ set: "
            f"missing={missing[:10]}, extra={extra[:10]}"
        )

    metadata_rows = read_tsv(input_path(cfg, "metadata"))
    metadata: dict[str, dict[str, str]] = {}
    for row in metadata_rows:
        sample = row.get("SampleID", "").strip()
        if not sample:
            continue
        if sample in metadata:
            raise ValueError(f"Duplicated SampleID in metadata: {sample}")
        metadata[sample] = row
    missing_meta = sorted(current_set - set(metadata))
    if missing_meta:
        raise ValueError(f"Current ex situ samples missing from metadata: {missing_meta[:10]}")

    sample_file = out / "00.current_exsitu.samples.txt"
    sample_file.write_text("\n".join(current_samples) + "\n", encoding="utf-8")

    vcf = input_path(cfg, "full_snp_vcf")
    tools = cfg.get("tools", {})
    tool = tools.get("bcftools", "bcftools")
    executable = shutil.which(str(tool)) or (str(tool) if Path(str(tool)).exists() else None)
    prefer_bcftools = bool(tools.get("prefer_bcftools", True))
    allow_direct = bool(tools.get("allow_direct_vcf_fallback", False))
    use_bcftools = prefer_bcftools and executable is not None
    if prefer_bcftools and executable is None and not allow_direct:
        raise RuntimeError(
            "bcftools was requested but is unavailable. Install bcftools or set "
            "tools.allow_direct_vcf_fallback=true (the direct mode is much slower on the full VCF)."
        )

    process: subprocess.Popen[str] | None = None
    if use_bcftools:
        header = subprocess.check_output(
            [str(executable), "view", "-h", "-S", str(sample_file), str(vcf)],
            text=True,
        )
        vcf_samples: list[str] | None = None
        for line in header.splitlines():
            if line.startswith("#CHROM"):
                vcf_samples = line.split("\t")[9:]
                break
        if vcf_samples is None:
            raise ValueError("bcftools did not return a #CHROM header")
        if set(vcf_samples) != current_set:
            raise ValueError("bcftools current-ex-situ subset does not match the configured sample list")
        process = subprocess.Popen(
            [str(executable), "view", "-H", "-S", str(sample_file), str(vcf)],
            stdout=subprocess.PIPE,
            stderr=(out / "bcftools.stderr.log").open("w"),
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        line_iter = process.stdout
        input_mode = "bcftools_subset"
    else:
        full_samples = read_vcf_header_samples(vcf)
        full_index = {sample: index for index, sample in enumerate(full_samples)}
        absent = sorted(current_set - set(full_index))
        if absent:
            raise ValueError(f"Current ex situ samples absent from full VCF: {absent[:10]}")
        selected_indices = [index for index, sample in enumerate(full_samples) if sample in current_set]
        vcf_samples = [full_samples[index] for index in selected_indices]
        if set(vcf_samples) != current_set:
            raise ValueError("Direct VCF subset does not match current ex situ sample list")

        def direct_subset_lines():
            with open_text(vcf, "rt") as handle:
                for raw_line in handle:
                    if raw_line.startswith("#"):
                        continue
                    full_fields = raw_line.rstrip("\n").split("\t")
                    if len(full_fields) != 9 + len(full_samples):
                        raise ValueError(
                            f"Malformed full VCF row near {full_fields[0] if full_fields else '?'}"
                        )
                    subset_fields = full_fields[:9] + [full_fields[9 + index] for index in selected_indices]
                    yield "\t".join(subset_fields) + "\n"

        line_iter = direct_subset_lines()
        input_mode = "direct_python_subset"

    group_flags: dict[str, dict[str, Any]] = {}
    for sample in vcf_samples:
        site = metadata[sample].get("Ex_situ_site", "NA") or "NA"
        group_flags[sample] = {
            "core": sample in core_set,
            "low": sample in low_set,
            "site": site,
        }

    bins = (
        "All", "AC_eq_1", "AC_2_to_5", "AF_le_0.01",
        "0.01_lt_AF_le_0.05", "0.05_lt_AF_le_0.10", "AF_gt_0.10",
    )
    state_total = Counter()
    state_retained = Counter()
    copy_total = Counter()
    copy_covered = Counter()
    core_copy_sum = Counter()
    lost_source_states = defaultdict(Counter)
    lost_source_copies = defaultdict(Counter)

    lost_path = out / "03.lost_ALT_states_after_deredundancy.tsv.gz"
    lost_fields = [
        "CHROM", "POS", "REF", "ALT_index", "ALT", "Current_AN", "Current_ALT_AC",
        "Current_ALT_AF", "Core_AN", "Core_ALT_AC", "Current_ALT_freq_bin",
        "Lost_current_ALT_copy_weight", "Loss_source_gardens",
    ]
    allow_partial = bool(cfg.get("parameters", {}).get("alt_coverage_allow_partial_genotypes", True))

    site_n = 0
    total_alt_records = 0
    bad_column_lines = 0

    with gzip.open(lost_path, "wt", encoding="utf-8", newline="") as lost_handle:
        lost_writer = csv.DictWriter(lost_handle, delimiter="\t", fieldnames=lost_fields, lineterminator="\n")
        lost_writer.writeheader()
        for line in line_iter:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 + len(vcf_samples):
                bad_column_lines += 1
                continue
            site_n += 1
            if site_n % 1_000_000 == 0:
                print(f"[INFO] ALT coverage: processed {site_n:,} VCF sites", flush=True)
            chrom, pos, ref = fields[0], fields[1], fields[3]
            alts = fields[4].split(",")
            n_alt = len(alts)
            total_alt_records += n_alt
            fmt = fields[8].split(":")
            if "GT" not in fmt:
                raise ValueError(f"GT is absent at {chrom}:{pos}")
            gt_index = fmt.index("GT")

            current_ac = [0] * n_alt
            core_ac = [0] * n_alt
            low_site_ac: dict[str, list[int]] = defaultdict(lambda: [0] * n_alt)
            current_an = 0
            core_an = 0

            for sample, sample_field in zip(vcf_samples, fields[9:]):
                values = sample_field.split(":")
                if gt_index >= len(values):
                    continue
                alleles = parse_gt_alt_coverage(values[gt_index], allow_partial)
                if not alleles:
                    continue
                flags = group_flags[sample]
                for allele in alleles:
                    if allele < 0 or allele > n_alt:
                        raise ValueError(f"GT index outside ALT range at {chrom}:{pos}")
                    current_an += 1
                    if 1 <= allele <= n_alt:
                        current_ac[allele - 1] += 1
                    if flags["core"]:
                        core_an += 1
                        if 1 <= allele <= n_alt:
                            core_ac[allele - 1] += 1
                    else:
                        if 1 <= allele <= n_alt:
                            low_site_ac[str(flags["site"])][allele - 1] += 1

            for index, alt in enumerate(alts):
                if len(ref) != 1 or len(alt) != 1 or ref.upper() not in "ACGT" or alt.upper() not in "ACGT" or ref.upper() == alt.upper():
                    continue
                ac_current = current_ac[index]
                if ac_current <= 0:
                    continue
                ac_core = core_ac[index]
                bin_label = alt_frequency_bin(ac_current, current_an)
                for key in ("All", bin_label):
                    state_total[key] += 1
                    copy_total[key] += ac_current
                    core_copy_sum[key] += ac_core
                    if ac_core > 0:
                        state_retained[key] += 1
                        copy_covered[key] += ac_current

                if ac_core <= 0:
                    low_copy_sum = sum(counts[index] for counts in low_site_ac.values())
                    if low_copy_sum != ac_current:
                        raise ValueError(
                            f"Lost ALT-copy partition mismatch at {chrom}:{pos}, ALT index {index + 1}: "
                            f"current={ac_current}, lower_priority={low_copy_sum}"
                        )
                    source_sites = sorted(site for site, counts in low_site_ac.items() if counts[index] > 0)
                    source = ";".join(source_sites) if source_sites else "Unexpected_no_low_priority_carrier"
                    for key in ("All", bin_label):
                        lost_source_states[source][key] += 1
                        lost_source_copies[source][key] += ac_current
                    lost_writer.writerow({
                        "CHROM": chrom,
                        "POS": pos,
                        "REF": ref,
                        "ALT_index": index + 1,
                        "ALT": alt,
                        "Current_AN": current_an,
                        "Current_ALT_AC": ac_current,
                        "Current_ALT_AF": ac_current / current_an if current_an else 0.0,
                        "Core_AN": core_an,
                        "Core_ALT_AC": ac_core,
                        "Current_ALT_freq_bin": bin_label,
                        "Lost_current_ALT_copy_weight": ac_current,
                        "Loss_source_gardens": source,
                    })

    if process is not None:
        stderr = "See bcftools.stderr.log"
        return_code = process.wait()
        if return_code != 0:
            raise RuntimeError(f"bcftools view failed:\n{stderr}")
    if bad_column_lines:
        raise ValueError(f"Encountered {bad_column_lines} malformed VCF rows")

    summary_rows: list[dict[str, Any]] = []
    for bin_label in bins:
        total_states = state_total[bin_label]
        retained_states = state_retained[bin_label]
        total_copies = copy_total[bin_label]
        covered_copies = copy_covered[bin_label]
        summary_rows.append({
            "Current_ALT_freq_bin": bin_label,
            "Current_ex_situ_observed_ALT_states": total_states,
            "ALT_states_retained_in_core": retained_states,
            "ALT_states_lost_after_deredundancy": total_states - retained_states,
            "ALT_state_presence_retention_rate": retained_states / total_states if total_states else 0.0,
            "Current_ex_situ_ALT_copies": total_copies,
            "Current_ALT_copies_covered_by_core_presence": covered_copies,
            "Current_ALT_copies_not_covered_by_core_presence": total_copies - covered_copies,
            "True_ALT_copy_weighted_coverage_rate": covered_copies / total_copies if total_copies else 0.0,
            "Core_ALT_copy_sum_auxiliary": core_copy_sum[bin_label],
        })
    write_tsv(out / "01.deredundancy_ALT_coverage_summary.tsv", summary_rows, list(summary_rows[0]))

    source_rows: list[dict[str, Any]] = []
    for source in sorted(lost_source_states):
        for bin_label in bins:
            source_rows.append({
                "Loss_source_gardens": source,
                "Current_ALT_freq_bin": bin_label,
                "Lost_ALT_states": lost_source_states[source][bin_label],
                "Lost_current_ALT_copies": lost_source_copies[source][bin_label],
            })
    source_fields = [
        "Loss_source_gardens", "Current_ALT_freq_bin",
        "Lost_ALT_states", "Lost_current_ALT_copies",
    ]
    write_tsv(out / "02.lost_ALT_source_summary.tsv", source_rows, source_fields)

    observed = {
        "SNP_observed_ALT_states": state_total["All"],
        "SNP_retained_ALT_states": state_retained["All"],
    }
    write_key_value(out / "04.internal_check.tsv", [
        ("VCF_sites_processed", site_n),
        ("Total_ALT_records", total_alt_records),
        ("Current_ex_situ_sample_N", len(current_set)),
        ("Core_sample_N", len(core_set)),
        ("Lower_priority_sample_N", len(low_set)),
        ("Core_and_lower_priority_partition_verified", True),
        ("Metric_definition", "presence_based_ALT_state_retention"),
        ("True_copy_weighted_metric_reported_separately", True),
        ("SNP_state_denominator", state_total["All"]),
        ("SNP_state_numerator", state_retained["All"]),
        ("SNP_state_retention_rate", state_retained["All"] / state_total["All"] if state_total["All"] else 0.0),
        ("True_copy_weight_denominator", copy_total["All"]),
        ("True_copy_weight_numerator", copy_covered["All"]),
        ("True_copy_weighted_coverage_rate", copy_covered["All"] / copy_total["All"] if copy_total["All"] else 0.0),
        ("Partial_genotypes_allowed_for_called_allele_counting", allow_partial),
        ("VCF_subset_mode", input_mode),
    ])
    print(f"[OK] ALT coverage metrics complete: {out}")


if __name__ == "__main__":
    main()
