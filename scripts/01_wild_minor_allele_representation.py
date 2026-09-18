#!/usr/bin/env python3
"""Wild minor-allele retention, AC robustness, and garden complementarity.

This script reconstructs the core allele-audit tables directly from the five
aligned AC/AN files.  It replaces the missing historical script that produced
the original retention/missing tables.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

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
    step_output_dir,
    target_allele_count,
    wild_minor_from_acan,
    write_key_value,
    write_tsv,
)

TARGETS = ("Ex_situ", "Pingshun", "Tongguan", "Zhendong")
ACAN_KEYS = {
    "Wild": "wild_acan",
    "Ex_situ": "exsitu_acan",
    "Pingshun": "pingshun_acan",
    "Tongguan": "tongguan_acan",
    "Zhendong": "zhendong_acan",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(default_config_path()))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    out = prepare_step_output(cfg, "01_allele_representation")

    paths = {label: input_path(cfg, key) for label, key in ACAN_KEYS.items()}

    total_by_target = defaultdict(int)
    retained_by_target = defaultdict(int)
    total_by_target_maf = defaultdict(int)
    retained_by_target_maf = defaultdict(int)
    total_by_threshold = defaultdict(int)
    missing_by_threshold = defaultdict(int)
    total_by_threshold_maf = defaultdict(int)
    missing_by_threshold_maf = defaultdict(int)
    complement = defaultdict(int)

    missing_table_path = out / "05.ex_situ_missing_wild_minor_alleles.tsv.gz"
    complement_site_path = out / "07.resource_garden_complementarity_by_site.tsv.gz"
    missing_sites_path = out / "08.ex_situ_missing_sites.tsv"

    missing_fields = [
        "CHROM", "POS", "Wild_AC", "Wild_AN", "Wild_REF_count", "Wild_ALT_count",
        "Wild_minor_allele", "Wild_minor_count", "Wild_minor_MAF", "Wild_MAF_bin",
        "Ex_situ_AC", "Ex_situ_AN", "Ex_situ_REF_count", "Ex_situ_ALT_count",
    ]
    complement_fields = [
        "CHROM", "POS", "Wild_minor_allele", "Wild_minor_count", "Wild_minor_MAF",
        "Wild_MAF_bin", "Pingshun_minor_count", "Tongguan_minor_count",
        "Zhendong_minor_count", "Complementarity_class",
    ]

    rows_scanned = 0
    wild_polymorphic = 0
    missing_all_exsitu = 0
    ps_missing_total = 0

    with gzip.open(missing_table_path, "wt", encoding="utf-8", newline="") as missing_handle, \
         gzip.open(complement_site_path, "wt", encoding="utf-8", newline="") as comp_handle, \
         missing_sites_path.open("w", encoding="utf-8", newline="") as site_handle:
        missing_writer = csv.DictWriter(missing_handle, delimiter="\t", fieldnames=missing_fields)
        comp_writer = csv.DictWriter(comp_handle, delimiter="\t", fieldnames=complement_fields)
        site_writer = csv.writer(site_handle, delimiter="\t")
        missing_writer.writeheader()
        comp_writer.writeheader()
        site_writer.writerow(["CHROM", "POS"])

        for rows in iter_aligned_acan(paths):
            rows_scanned += 1
            wild = rows["Wild"]
            definition = wild_minor_from_acan(wild)
            if definition is None:
                continue
            wild_polymorphic += 1
            minor, minor_count, maf = definition
            mbin = maf_bin(maf)

            counts = {target: target_allele_count(rows[target], minor) for target in TARGETS}

            for target in TARGETS:
                total_by_target[target] += 1
                total_by_target_maf[(target, mbin)] += 1
                if counts[target] > 0:
                    retained_by_target[target] += 1
                    retained_by_target_maf[(target, mbin)] += 1

            exsitu_missing = counts["Ex_situ"] <= 0
            for threshold_name in THRESHOLD_ORDER:
                if minor_count < THRESHOLDS[threshold_name]:
                    continue
                total_by_threshold[threshold_name] += 1
                total_by_threshold_maf[(threshold_name, mbin)] += 1
                if exsitu_missing:
                    missing_by_threshold[threshold_name] += 1
                    missing_by_threshold_maf[(threshold_name, mbin)] += 1

            if exsitu_missing:
                missing_all_exsitu += 1
                missing_writer.writerow({
                    "CHROM": wild.chrom,
                    "POS": wild.pos,
                    "Wild_AC": wild.alt_ac,
                    "Wild_AN": wild.an,
                    "Wild_REF_count": wild.an - wild.alt_ac,
                    "Wild_ALT_count": wild.alt_ac,
                    "Wild_minor_allele": minor,
                    "Wild_minor_count": minor_count,
                    "Wild_minor_MAF": f"{maf:.17g}",
                    "Wild_MAF_bin": mbin,
                    "Ex_situ_AC": rows["Ex_situ"].alt_ac,
                    "Ex_situ_AN": rows["Ex_situ"].an,
                    "Ex_situ_REF_count": rows["Ex_situ"].an - rows["Ex_situ"].alt_ac,
                    "Ex_situ_ALT_count": rows["Ex_situ"].alt_ac,
                })
                site_writer.writerow([wild.chrom, wild.pos])

            if counts["Pingshun"] <= 0:
                ps_missing_total += 1
                tg_present = counts["Tongguan"] > 0
                zd_present = counts["Zhendong"] > 0
                if tg_present and zd_present:
                    category = "Both_TG_RG_and_ZD_RG"
                elif tg_present:
                    category = "TG_RG_only"
                elif zd_present:
                    category = "ZD_RG_only"
                else:
                    category = "Still_missing_from_all_ex_situ"
                complement[("All", category)] += 1
                complement[(mbin, category)] += 1
                comp_writer.writerow({
                    "CHROM": wild.chrom,
                    "POS": wild.pos,
                    "Wild_minor_allele": minor,
                    "Wild_minor_count": minor_count,
                    "Wild_minor_MAF": f"{maf:.17g}",
                    "Wild_MAF_bin": mbin,
                    "Pingshun_minor_count": counts["Pingshun"],
                    "Tongguan_minor_count": counts["Tongguan"],
                    "Zhendong_minor_count": counts["Zhendong"],
                    "Complementarity_class": category,
                })

            if rows_scanned % 5_000_000 == 0:
                print(
                    f"[INFO] rows={rows_scanned:,}; wild-polymorphic={wild_polymorphic:,}; "
                    f"ex-situ-missing={missing_all_exsitu:,}",
                    file=sys.stderr,
                )

    overall_rows = []
    for target in TARGETS:
        total = total_by_target[target]
        retained = retained_by_target[target]
        missing = total - retained
        overall_rows.append({
            "Target_group": target,
            "Wild_minor_alleles": total,
            "Retained_Wild_minor_alleles": retained,
            "Missing_Wild_minor_alleles": missing,
            "Minor_allele_retention_rate": retained / total if total else float("nan"),
            "Missing_rate": missing / total if total else float("nan"),
        })
    write_tsv(
        out / "01.wild_minor_allele_retention_overall.tsv",
        overall_rows,
        ["Target_group", "Wild_minor_alleles", "Retained_Wild_minor_alleles",
         "Missing_Wild_minor_alleles", "Minor_allele_retention_rate", "Missing_rate"],
    )

    maf_rows = []
    for target in TARGETS:
        for mbin in MAF_BIN_ORDER:
            total = total_by_target_maf[(target, mbin)]
            retained = retained_by_target_maf[(target, mbin)]
            missing = total - retained
            maf_rows.append({
                "Target_group": target,
                "Wild_MAF_bin": mbin,
                "Wild_minor_alleles": total,
                "Retained_Wild_minor_alleles": retained,
                "Missing_Wild_minor_alleles": missing,
                "Minor_allele_retention_rate": retained / total if total else float("nan"),
                "Missing_rate": missing / total if total else float("nan"),
            })
    write_tsv(
        out / "02.wild_minor_allele_retention_by_MAF_bin.tsv",
        maf_rows,
        ["Target_group", "Wild_MAF_bin", "Wild_minor_alleles",
         "Retained_Wild_minor_alleles", "Missing_Wild_minor_alleles",
         "Minor_allele_retention_rate", "Missing_rate"],
    )

    threshold_rows = []
    threshold_maf_rows = []
    for threshold_name in THRESHOLD_ORDER:
        total = total_by_threshold[threshold_name]
        missing = missing_by_threshold[threshold_name]
        retained = total - missing
        threshold_rows.append({
            "AC_threshold": threshold_name,
            "Wild_minor_count_min": THRESHOLDS[threshold_name],
            "Total_Wild_SNP_minor_alleles": total,
            "Missing_in_all_Ex_situ": missing,
            "Retained_in_Ex_situ": retained,
            "Retention_rate": retained / total if total else float("nan"),
            "Missing_rate": missing / total if total else float("nan"),
        })
        for mbin in MAF_BIN_ORDER:
            total_bin = total_by_threshold_maf[(threshold_name, mbin)]
            missing_bin = missing_by_threshold_maf[(threshold_name, mbin)]
            retained_bin = total_bin - missing_bin
            threshold_maf_rows.append({
                "AC_threshold": threshold_name,
                "Wild_minor_count_min": THRESHOLDS[threshold_name],
                "Wild_MAF_bin": mbin,
                "Total_Wild_SNP_minor_alleles": total_bin,
                "Missing_in_all_Ex_situ": missing_bin,
                "Retained_in_Ex_situ": retained_bin,
                "Retention_rate": retained_bin / total_bin if total_bin else float("nan"),
                "Missing_rate": missing_bin / total_bin if total_bin else float("nan"),
            })
    write_tsv(
        out / "03.retention_AC_threshold_robustness.tsv",
        threshold_rows,
        ["AC_threshold", "Wild_minor_count_min", "Total_Wild_SNP_minor_alleles",
         "Missing_in_all_Ex_situ", "Retained_in_Ex_situ", "Retention_rate", "Missing_rate"],
    )
    write_tsv(
        out / "04.retention_AC_threshold_by_MAF_bin.tsv",
        threshold_maf_rows,
        ["AC_threshold", "Wild_minor_count_min", "Wild_MAF_bin",
         "Total_Wild_SNP_minor_alleles", "Missing_in_all_Ex_situ",
         "Retained_in_Ex_situ", "Retention_rate", "Missing_rate"],
    )

    categories = (
        "TG_RG_only", "ZD_RG_only", "Both_TG_RG_and_ZD_RG",
        "Still_missing_from_all_ex_situ",
    )
    complement_rows = []
    for mbin in ("All", *MAF_BIN_ORDER):
        counts = {category: complement[(mbin, category)] for category in categories}
        total = sum(counts.values())
        recovered = total - counts["Still_missing_from_all_ex_situ"]
        complement_rows.append({
            "Wild_MAF_bin": mbin,
            "PS_RG_missing_Wild_minor_alleles": total,
            **counts,
            "Recovered_by_TG_RG_or_ZD_RG": recovered,
            "Recovery_rate_by_Tongguan_or_Zhendong": recovered / total if total else float("nan"),
            "Still_missing_rate_in_all_exsitu": counts["Still_missing_from_all_ex_situ"] / total if total else float("nan"),
        })
    write_tsv(
        out / "06.resource_garden_complementarity_summary.tsv",
        complement_rows,
        ["Wild_MAF_bin", "PS_RG_missing_Wild_minor_alleles", *categories,
         "Recovered_by_TG_RG_or_ZD_RG", "Recovery_rate_by_Tongguan_or_Zhendong",
         "Still_missing_rate_in_all_exsitu"],
    )

    observed = {
        "Total_wild_SNP_minor_alleles": wild_polymorphic,
        "Combined_exsitu_retained": retained_by_target["Ex_situ"],
        "Combined_exsitu_missing": missing_all_exsitu,
        "Pingshun_retained": retained_by_target["Pingshun"],
        "Tongguan_retained": retained_by_target["Tongguan"],
        "Zhendong_retained": retained_by_target["Zhendong"],
        "AC_ge_2_missing": missing_by_threshold["AC_ge_2"],
        "AC_ge_3_missing": missing_by_threshold["AC_ge_3"],
    }
    check_items = [
        ("AC_AN_rows_scanned", rows_scanned),
        ("Wild_polymorphic_SNPs", wild_polymorphic),
        ("Ex_situ_missing_Wild_minor_alleles", missing_all_exsitu),
        ("PS_RG_missing_Wild_minor_alleles", ps_missing_total),
    ] + list(observed.items())
    write_key_value(out / "09.internal_check.tsv", check_items)

    print(f"[OK] Allele-gap audit complete: {out}")


if __name__ == "__main__":
    main()
