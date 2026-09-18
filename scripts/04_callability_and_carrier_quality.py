#!/usr/bin/env python3
"""Callability and Wild-carrier DP/GQ checks for ex situ-missing alleles.

Site call rates are calculated from AC/AN tables.  Wild-carrier genotype depth
and genotype quality are summarized from the missing-site VCF.  The script is
strict about sample membership and about agreement between the reconstructed
missing status and the Step 01 missing-allele table.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from typing import Any

from analysis_utils import (
    prepare_step_output,
    MAF_BIN_ORDER,
    THRESHOLD_ORDER,
    THRESHOLDS,
    default_config_path,
    input_path,
    iter_aligned_acan,
    load_config,
    load_missing_targets,
    maf_bin,
    open_text,
    parse_gt,
    read_sample_list,
    step_output_dir,
    target_allele_count,
    wild_minor_from_acan,
    write_key_value,
    write_tsv,
)


class HistogramStats:
    """Memory-bounded call-rate summaries at 0.001 resolution."""

    def __init__(self) -> None:
        self.n = 0
        self.total = 0.0
        self.minimum: float | None = None
        self.maximum: float | None = None
        self.hist = [0] * 1001
        self.lt_050 = 0
        self.lt_080 = 0
        self.lt_090 = 0
        self.lt_095 = 0

    def add(self, value: float) -> None:
        value = max(0.0, min(1.0, value))
        self.n += 1
        self.total += value
        self.minimum = value if self.minimum is None else min(self.minimum, value)
        self.maximum = value if self.maximum is None else max(self.maximum, value)
        self.hist[int(round(value * 1000))] += 1
        self.lt_050 += value < 0.50
        self.lt_080 += value < 0.80
        self.lt_090 += value < 0.90
        self.lt_095 += value < 0.95

    def mean(self) -> float:
        return self.total / self.n if self.n else math.nan

    def quantile(self, q: float) -> float:
        if not self.n:
            return math.nan
        target = q * (self.n - 1)
        cumulative = 0
        for index, count in enumerate(self.hist):
            cumulative += count
            if cumulative - 1 >= target:
                return index / 1000.0
        return 1.0


class ValueStats:
    def __init__(self) -> None:
        self.values: list[float] = []

    def add(self, value: str | None) -> None:
        if value in {None, "", ".", "NA"}:
            return
        try:
            self.values.append(float(value))
        except ValueError:
            return

    def summary(self) -> dict[str, Any]:
        if not self.values:
            return {key: math.nan for key in ("Mean", "Min", "P05", "P25", "Median", "P75", "P95", "Max")} | {"N": 0}
        values = sorted(self.values)

        def q(probability: float) -> float:
            if len(values) == 1:
                return values[0]
            position = probability * (len(values) - 1)
            low = math.floor(position)
            high = math.ceil(position)
            fraction = position - low
            return values[low] * (1.0 - fraction) + values[high] * fraction

        return {
            "N": len(values),
            "Mean": sum(values) / len(values),
            "Min": values[0],
            "P05": q(0.05),
            "P25": q(0.25),
            "Median": q(0.50),
            "P75": q(0.75),
            "P95": q(0.95),
            "Max": values[-1],
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--missing-table", default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    out = prepare_step_output(cfg, "04_callability_and_carrier_quality")
    missing_path = (
        input_path(cfg, "missing_table_override")
        if args.missing_table is None and cfg["inputs"].get("missing_table_override")
        else None
    )
    if args.missing_table:
        from pathlib import Path
        missing_path = Path(args.missing_table)
    if missing_path is None:
        missing_path = step_output_dir(cfg, "01_allele_representation") / "05.ex_situ_missing_wild_minor_alleles.tsv.gz"

    targets, target_order = load_missing_targets(missing_path)
    target_set = set(targets)
    wild_samples = read_sample_list(input_path(cfg, "wild_samples"))
    exsitu_samples = read_sample_list(input_path(cfg, "exsitu_samples"))
    ploidy = int(cfg.get("parameters", {}).get("ploidy", 2))
    n_wild = len(wild_samples)
    n_exsitu = len(exsitu_samples)

    call_stats: dict[tuple[str, str, str, str], HistogramStats] = defaultdict(HistogramStats)
    # keys: status, threshold, MAF/All, group(Wild/Ex_situ)
    matched_missing = 0
    reconstructed_disagreement = 0
    rows_scanned = 0

    site_path = out / "04.ex_situ_missing_site_callability_by_site.tsv.gz"
    site_fields = [
        "CHROM", "POS", "Wild_minor_allele", "Wild_minor_count", "Wild_minor_MAF",
        "Wild_MAF_bin", "Wild_AN", "Ex_situ_AN", "Wild_call_rate", "Ex_situ_call_rate",
    ]
    with gzip.open(site_path, "wt", encoding="utf-8", newline="") as site_handle:
        site_writer = csv.DictWriter(site_handle, delimiter="\t", fieldnames=site_fields, lineterminator="\n")
        site_writer.writeheader()
        for rows in iter_aligned_acan({
            "Wild": input_path(cfg, "wild_acan"),
            "Ex_situ": input_path(cfg, "exsitu_acan"),
        }):
            rows_scanned += 1
            definition = wild_minor_from_acan(rows["Wild"])
            if definition is None:
                continue
            minor, minor_count, minor_maf = definition
            mbin = maf_bin(minor_maf)
            key = (rows["Wild"].chrom, rows["Wild"].pos)
            exsitu_count = target_allele_count(rows["Ex_situ"], minor)
            reconstructed_missing = exsitu_count <= 0
            listed_missing = key in target_set
            if reconstructed_missing != listed_missing:
                reconstructed_disagreement += 1
                continue
            status = "Ex_situ_missing" if reconstructed_missing else "Ex_situ_retained"
            if reconstructed_missing:
                matched_missing += 1

            wild_denominator = ploidy * n_wild
            exsitu_denominator = ploidy * n_exsitu
            if rows["Wild"].an > wild_denominator or rows["Ex_situ"].an > exsitu_denominator:
                raise ValueError(
                    f"AN exceeds ploidy × sample count at {key[0]}:{key[1]}: "
                    f"Wild_AN={rows['Wild'].an}/{wild_denominator}, "
                    f"Ex_situ_AN={rows['Ex_situ'].an}/{exsitu_denominator}"
                )
            wild_call = rows["Wild"].an / wild_denominator
            exsitu_call = rows["Ex_situ"].an / exsitu_denominator
            for threshold in THRESHOLD_ORDER:
                if minor_count < THRESHOLDS[threshold]:
                    continue
                for bin_label in (mbin, "All"):
                    call_stats[(status, threshold, bin_label, "Wild")].add(wild_call)
                    call_stats[(status, threshold, bin_label, "Ex_situ")].add(exsitu_call)

            if reconstructed_missing:
                site_writer.writerow({
                    "CHROM": key[0], "POS": key[1], "Wild_minor_allele": minor,
                    "Wild_minor_count": minor_count, "Wild_minor_MAF": minor_maf,
                    "Wild_MAF_bin": mbin, "Wild_AN": rows["Wild"].an,
                    "Ex_situ_AN": rows["Ex_situ"].an,
                    "Wild_call_rate": wild_call, "Ex_situ_call_rate": exsitu_call,
                })

    if reconstructed_disagreement:
        raise ValueError(
            f"{reconstructed_disagreement} sites disagreed between AC/AN reconstruction and the missing table"
        )
    if matched_missing != len(targets):
        raise ValueError(f"Matched {matched_missing} missing sites, expected {len(targets)}")

    summary_rows: list[dict[str, Any]] = []
    for status in ("Ex_situ_missing", "Ex_situ_retained"):
        for threshold in THRESHOLD_ORDER:
            for bin_label in (*MAF_BIN_ORDER, "All"):
                wild = call_stats[(status, threshold, bin_label, "Wild")]
                exsitu = call_stats[(status, threshold, bin_label, "Ex_situ")]
                if wild.n != exsitu.n:
                    raise ValueError("Wild and ex situ callability denominators differ")
                summary_rows.append({
                    "Allele_status": status,
                    "AC_threshold": threshold,
                    "Wild_minor_count_min": THRESHOLDS[threshold],
                    "Wild_MAF_bin": bin_label,
                    "Site_N": wild.n,
                    "Wild_call_rate_mean": wild.mean(),
                    "Wild_call_rate_min": wild.minimum,
                    "Wild_call_rate_p01_approx": wild.quantile(0.01),
                    "Wild_call_rate_p05_approx": wild.quantile(0.05),
                    "Wild_call_rate_median_approx": wild.quantile(0.50),
                    "Wild_call_rate_p95_approx": wild.quantile(0.95),
                    "Wild_call_rate_max": wild.maximum,
                    "Wild_call_rate_lt_0.50_N": wild.lt_050,
                    "Wild_call_rate_lt_0.80_N": wild.lt_080,
                    "Wild_call_rate_lt_0.90_N": wild.lt_090,
                    "Wild_call_rate_lt_0.95_N": wild.lt_095,
                    "Wild_call_rate_lt_0.50_fraction": wild.lt_050 / wild.n if wild.n else math.nan,
                    "Wild_call_rate_lt_0.80_fraction": wild.lt_080 / wild.n if wild.n else math.nan,
                    "Wild_call_rate_lt_0.90_fraction": wild.lt_090 / wild.n if wild.n else math.nan,
                    "Wild_call_rate_lt_0.95_fraction": wild.lt_095 / wild.n if wild.n else math.nan,
                    "Ex_situ_call_rate_mean": exsitu.mean(),
                    "Ex_situ_call_rate_min": exsitu.minimum,
                    "Ex_situ_call_rate_p01_approx": exsitu.quantile(0.01),
                    "Ex_situ_call_rate_p05_approx": exsitu.quantile(0.05),
                    "Ex_situ_call_rate_median_approx": exsitu.quantile(0.50),
                    "Ex_situ_call_rate_p95_approx": exsitu.quantile(0.95),
                    "Ex_situ_call_rate_max": exsitu.maximum,
                    "Ex_situ_call_rate_lt_0.50_N": exsitu.lt_050,
                    "Ex_situ_call_rate_lt_0.80_N": exsitu.lt_080,
                    "Ex_situ_call_rate_lt_0.90_N": exsitu.lt_090,
                    "Ex_situ_call_rate_lt_0.95_N": exsitu.lt_095,
                    "Ex_situ_call_rate_lt_0.50_fraction": exsitu.lt_050 / exsitu.n if exsitu.n else math.nan,
                    "Ex_situ_call_rate_lt_0.80_fraction": exsitu.lt_080 / exsitu.n if exsitu.n else math.nan,
                    "Ex_situ_call_rate_lt_0.90_fraction": exsitu.lt_090 / exsitu.n if exsitu.n else math.nan,
                    "Ex_situ_call_rate_lt_0.95_fraction": exsitu.lt_095 / exsitu.n if exsitu.n else math.nan,
                })
    call_fields = list(summary_rows[0])
    write_tsv(out / "01.missing_vs_retained_site_callability_summary.tsv", summary_rows, call_fields)
    write_tsv(
        out / "02.ex_situ_missing_site_callability_summary.tsv",
        (row for row in summary_rows if row["Allele_status"] == "Ex_situ_missing"),
        call_fields,
    )

    # Wild-carrier genotype DP/GQ at missing sites.
    carrier_dp: dict[tuple[str, str], ValueStats] = defaultdict(ValueStats)
    carrier_gq: dict[tuple[str, str], ValueStats] = defaultdict(ValueStats)
    carrier_genotypes = defaultdict(int)
    carrier_copies = defaultdict(int)
    matched_vcf: set[tuple[str, str]] = set()
    vcf_rows = 0
    format_dp_sites = 0
    format_gq_sites = 0

    with open_text(input_path(cfg, "missing_site_vcf"), "rt") as handle:
        sample_names: list[str] | None = None
        wild_indices: list[int] = []
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_names = line.rstrip("\n").split("\t")[9:]
                index = {sample: i for i, sample in enumerate(sample_names)}
                absent = [sample for sample in wild_samples if sample not in index]
                if absent:
                    raise ValueError(f"Wild samples absent from missing-site VCF: {absent[:10]}")
                wild_indices = [index[sample] for sample in wild_samples]
                continue
            if sample_names is None:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 + len(sample_names):
                raise ValueError(f"Malformed VCF row near {fields[0] if fields else '?'}")
            key = (fields[0], fields[1])
            if key not in targets:
                continue
            if key in matched_vcf:
                raise ValueError(f"Duplicate target in VCF: {key[0]}:{key[1]}")
            matched_vcf.add(key)
            vcf_rows += 1
            info = targets[key]
            fmt = fields[8].split(":")
            if "GT" not in fmt:
                raise ValueError(f"GT is absent at {key[0]}:{key[1]}")
            fmt_index = {name: i for i, name in enumerate(fmt)}
            format_dp_sites += "DP" in fmt_index
            format_gq_sites += "GQ" in fmt_index
            genotypes = fields[9:]
            target_code = 1 if info["minor"] == "ALT" else 0
            for index in wild_indices:
                values = genotypes[index].split(":")
                gt = values[fmt_index["GT"]] if fmt_index["GT"] < len(values) else "."
                alleles = parse_gt(gt)
                if not alleles:
                    continue
                dosage = sum(1 for allele in alleles if allele == target_code)
                if dosage <= 0:
                    continue
                dp = values[fmt_index["DP"]] if "DP" in fmt_index and fmt_index["DP"] < len(values) else None
                gq = values[fmt_index["GQ"]] if "GQ" in fmt_index and fmt_index["GQ"] < len(values) else None
                for threshold in THRESHOLD_ORDER:
                    if info["minor_count"] < THRESHOLDS[threshold]:
                        continue
                    for bin_label in (info["maf_bin"], "All"):
                        carrier_genotypes[(threshold, bin_label)] += 1
                        carrier_copies[(threshold, bin_label)] += dosage
                        carrier_dp[(threshold, bin_label)].add(dp)
                        carrier_gq[(threshold, bin_label)].add(gq)

    missing_vcf_targets = target_set - matched_vcf
    if missing_vcf_targets:
        preview = ", ".join(f"{c}:{p}" for c, p in sorted(missing_vcf_targets)[:10])
        raise ValueError(f"Missing-site VCF lacks {len(missing_vcf_targets)} targets; examples: {preview}")

    carrier_rows: list[dict[str, Any]] = []
    for threshold in THRESHOLD_ORDER:
        for bin_label in (*MAF_BIN_ORDER, "All"):
            dp = carrier_dp[(threshold, bin_label)].summary()
            gq = carrier_gq[(threshold, bin_label)].summary()
            row: dict[str, Any] = {
                "AC_threshold": threshold,
                "Wild_minor_count_min": THRESHOLDS[threshold],
                "Wild_MAF_bin": bin_label,
                "Carrier_genotype_N": carrier_genotypes[(threshold, bin_label)],
                "Carrier_target_allele_copies": carrier_copies[(threshold, bin_label)],
            }
            row.update({f"DP_{key}": value for key, value in dp.items()})
            row.update({f"GQ_{key}": value for key, value in gq.items()})
            carrier_rows.append(row)
    strict_carrier_qc = bool(cfg.get("parameters", {}).get("strict_carrier_DP_GQ", True))
    all_carriers = carrier_genotypes[("all", "All")]
    all_dp_n = len(carrier_dp[("all", "All")].values)
    all_gq_n = len(carrier_gq[("all", "All")].values)
    if strict_carrier_qc and (format_dp_sites != len(targets) or format_gq_sites != len(targets)):
        raise ValueError(
            "DP/GQ FORMAT fields are not present at every missing target site: "
            f"DP={format_dp_sites}/{len(targets)}, GQ={format_gq_sites}/{len(targets)}"
        )
    if strict_carrier_qc and (all_dp_n != all_carriers or all_gq_n != all_carriers):
        raise ValueError(
            "Some Wild carrier genotypes lack numeric DP/GQ values: "
            f"carriers={all_carriers}, DP={all_dp_n}, GQ={all_gq_n}"
        )

    write_tsv(out / "03.Wild_carrier_DP_GQ_summary.tsv", carrier_rows, list(carrier_rows[0]))

    write_key_value(out / "05.internal_check.tsv", [
        ("Wild_sample_N", n_wild),
        ("Ex_situ_sample_N", n_exsitu),
        ("ACAN_rows_scanned", rows_scanned),
        ("Missing_targets_loaded", len(targets)),
        ("Missing_targets_matched_in_ACAN", matched_missing),
        ("Missing_targets_matched_in_VCF", len(matched_vcf)),
        ("VCF_target_rows", vcf_rows),
        ("VCF_target_rows_with_DP", format_dp_sites),
        ("VCF_target_rows_with_GQ", format_gq_sites),
        ("Wild_carrier_genotypes_all", all_carriers),
        ("Wild_carrier_genotypes_with_DP_all", all_dp_n),
        ("Wild_carrier_genotypes_with_GQ_all", all_gq_n),
        ("Strict_carrier_DP_GQ_validation", strict_carrier_qc),
    ])
    print(f"[OK] Callability and carrier-quality checks complete: {out}")


if __name__ == "__main__":
    main()
