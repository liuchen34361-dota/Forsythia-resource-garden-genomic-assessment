#!/usr/bin/env python3
"""Deterministic sample-size standardization of Wild minor-allele retention.

PS-RG retention is rarefied at the allele-copy level to the TG-RG and ZD-RG
sample sizes.  The expected probability of detecting at least one copy of each
Wild SNP minor allele is calculated from the hypergeometric distribution.
"""
from __future__ import annotations

import argparse
import math
from collections import defaultdict
from functools import lru_cache
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
    maf_bin,
    read_sample_list,
    step_output_dir,
    target_allele_count,
    wild_minor_from_acan,
    write_key_value,
    write_tsv,
)


@lru_cache(maxsize=None)
def detection_probability(an: int, target_ac: int, draw_n: int) -> float:
    """Return P(X >= 1) for a hypergeometric draw without replacement."""
    if an <= 0 or target_ac <= 0 or draw_n <= 0:
        return 0.0
    draw_n = min(draw_n, an)
    if target_ac >= an or draw_n >= an or an - target_ac < draw_n:
        return 1.0
    # P(X=0) = C(an-target_ac, draw_n) / C(an, draw_n).
    log_p0 = (
        math.lgamma(an - target_ac + 1)
        - math.lgamma(draw_n + 1)
        - math.lgamma(an - target_ac - draw_n + 1)
        - math.lgamma(an + 1)
        + math.lgamma(draw_n + 1)
        + math.lgamma(an - draw_n + 1)
    )
    return max(0.0, min(1.0, 1.0 - math.exp(log_p0)))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(default_config_path()))
    return parser.parse_args()


def main() -> None:
    cfg = load_config(parse_args().config)
    out = prepare_step_output(cfg, "03_sample_size_standardized_retention")

    paths = {
        "Wild": input_path(cfg, "wild_acan"),
        "Pingshun": input_path(cfg, "pingshun_acan"),
        "Tongguan": input_path(cfg, "tongguan_acan"),
        "Zhendong": input_path(cfg, "zhendong_acan"),
    }
    n_ps = len(read_sample_list(input_path(cfg, "pingshun_samples")))
    n_tg = len(read_sample_list(input_path(cfg, "tongguan_samples")))
    n_zd = len(read_sample_list(input_path(cfg, "zhendong_samples")))
    if min(n_ps, n_tg, n_zd) <= 0:
        raise ValueError("Resource-garden sample lists must be non-empty")

    targets = (
        ("PS_RG_observed", "observed", n_ps),
        ("TG_RG_observed", "observed", n_tg),
        ("ZD_RG_observed", "observed", n_zd),
        (f"PS_RG_rarefied_to_TG_RG_sample_size_n{n_tg}_expected", "allele_copy_rarefied_expected", n_tg),
        (f"PS_RG_rarefied_to_ZD_RG_sample_size_n{n_zd}_expected", "allele_copy_rarefied_expected", n_zd),
    )

    totals = defaultdict(int)
    totals_maf = defaultdict(int)
    retained = defaultdict(float)
    retained_maf = defaultdict(float)
    rows_scanned = 0
    wild_polymorphic = 0

    for rows in iter_aligned_acan(paths):
        rows_scanned += 1
        definition = wild_minor_from_acan(rows["Wild"])
        if definition is None:
            continue
        wild_polymorphic += 1
        minor, minor_count, minor_maf = definition
        mbin = maf_bin(minor_maf)

        ps_count = target_allele_count(rows["Pingshun"], minor)
        tg_count = target_allele_count(rows["Tongguan"], minor)
        zd_count = target_allele_count(rows["Zhendong"], minor)
        ps_an = rows["Pingshun"].an

        # Preserve the historical deterministic implementation: site-level
        # called allele copies are scaled in proportion to accession number.
        draw_to_tg = min(ps_an, max(0, int(round(ps_an * n_tg / n_ps))))
        draw_to_zd = min(ps_an, max(0, int(round(ps_an * n_zd / n_ps))))

        values = {
            "PS_RG_observed": 1.0 if ps_count > 0 else 0.0,
            "TG_RG_observed": 1.0 if tg_count > 0 else 0.0,
            "ZD_RG_observed": 1.0 if zd_count > 0 else 0.0,
            f"PS_RG_rarefied_to_TG_RG_sample_size_n{n_tg}_expected": detection_probability(ps_an, ps_count, draw_to_tg),
            f"PS_RG_rarefied_to_ZD_RG_sample_size_n{n_zd}_expected": detection_probability(ps_an, ps_count, draw_to_zd),
        }

        for threshold in THRESHOLD_ORDER:
            if minor_count < THRESHOLDS[threshold]:
                continue
            totals[threshold] += 1
            totals_maf[(threshold, mbin)] += 1
            for target_name, _, _ in targets:
                retained[(threshold, target_name)] += values[target_name]
                retained_maf[(threshold, target_name, mbin)] += values[target_name]

    overall_rows: list[dict[str, Any]] = []
    for target_name, comparison_type, nominal_n in targets:
        total = totals["all"]
        value = retained[("all", target_name)]
        overall_rows.append({
            "Target_set": target_name,
            "Comparison_type": comparison_type,
            "Nominal_sample_size": nominal_n,
            "Total_Wild_SNP_minor_alleles": total,
            "Retained_Wild_SNP_minor_alleles": value,
            "Retention_rate": value / total if total else math.nan,
            "Missing_or_not_retained": total - value,
            "Missing_rate": (total - value) / total if total else math.nan,
        })
    write_tsv(
        out / "01.sample_size_standardized_retention_overall.tsv",
        overall_rows,
        ["Target_set", "Comparison_type", "Nominal_sample_size", "Total_Wild_SNP_minor_alleles",
         "Retained_Wild_SNP_minor_alleles", "Retention_rate", "Missing_or_not_retained", "Missing_rate"],
    )

    maf_rows: list[dict[str, Any]] = []
    for target_name, comparison_type, nominal_n in targets:
        for mbin in MAF_BIN_ORDER:
            total = totals_maf[("all", mbin)]
            value = retained_maf[("all", target_name, mbin)]
            maf_rows.append({
                "Target_set": target_name,
                "Comparison_type": comparison_type,
                "Nominal_sample_size": nominal_n,
                "Wild_MAF_bin": mbin,
                "Total_Wild_SNP_minor_alleles": total,
                "Retained_Wild_SNP_minor_alleles": value,
                "Retention_rate": value / total if total else math.nan,
                "Missing_or_not_retained": total - value,
                "Missing_rate": (total - value) / total if total else math.nan,
            })
    write_tsv(
        out / "02.sample_size_standardized_retention_by_MAF_bin.tsv",
        maf_rows,
        ["Target_set", "Comparison_type", "Nominal_sample_size", "Wild_MAF_bin",
         "Total_Wild_SNP_minor_alleles", "Retained_Wild_SNP_minor_alleles",
         "Retention_rate", "Missing_or_not_retained", "Missing_rate"],
    )

    threshold_rows: list[dict[str, Any]] = []
    for threshold in THRESHOLD_ORDER:
        for target_name, comparison_type, nominal_n in targets:
            total = totals[threshold]
            value = retained[(threshold, target_name)]
            threshold_rows.append({
                "AC_threshold": threshold,
                "Wild_minor_count_min": THRESHOLDS[threshold],
                "Target_set": target_name,
                "Comparison_type": comparison_type,
                "Nominal_sample_size": nominal_n,
                "Total_Wild_SNP_minor_alleles": total,
                "Retained_Wild_SNP_minor_alleles": value,
                "Retention_rate": value / total if total else math.nan,
                "Missing_or_not_retained": total - value,
                "Missing_rate": (total - value) / total if total else math.nan,
            })
    write_tsv(
        out / "03.sample_size_standardized_retention_by_AC_threshold.tsv",
        threshold_rows,
        ["AC_threshold", "Wild_minor_count_min", "Target_set", "Comparison_type",
         "Nominal_sample_size", "Total_Wild_SNP_minor_alleles",
         "Retained_Wild_SNP_minor_alleles", "Retention_rate",
         "Missing_or_not_retained", "Missing_rate"],
    )

    write_key_value(out / "04.internal_check.tsv", [
        ("ACAN_rows_scanned", rows_scanned),
        ("Wild_polymorphic_sites", wild_polymorphic),
        ("PS_RG_sample_N", n_ps),
        ("TG_RG_sample_N", n_tg),
        ("ZD_RG_sample_N", n_zd),
        ("Hypergeometric_cache_size", detection_probability.cache_info().currsize),
    ])
    print(f"[OK] Sample-size standardization complete: {out}")


if __name__ == "__main__":
    main()
