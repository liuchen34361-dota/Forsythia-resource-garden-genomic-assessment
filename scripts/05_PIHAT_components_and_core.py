#!/usr/bin/env python3
"""PI_HAT threshold sensitivity and deterministic ex situ de-redundancy.

Global connected components are retained for manuscript-wide threshold
sensitivity and for selected-Wild-accession QC.  The operational 334-accession
ex situ core is constructed explicitly on the induced ex situ graph at
PI_HAT >= 0.80.  Singletons are retained; within every multi-accession
component the lexicographically first SampleID is retained.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

from analysis_utils import (
    prepare_step_output,
    canonical_region,
    connected_components,
    default_config_path,
    ensure_unique,
    input_path,
    load_config,
    read_sample_list,
    read_tsv,
    step_output_dir,
    write_key_value,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(default_config_path()))
    return parser.parse_args()


def clean(value: Any) -> str:
    text = "" if value is None else str(value).strip().replace("\r", "")
    return "NA" if not text or text.lower() in {"na", "nan", "none"} else text


def normalized_material(value: Any) -> str:
    text = clean(value).lower().replace("-", "_").replace(" ", "_")
    if text == "wild":
        return "Wild"
    if text in {"ex_situ", "exsitu"}:
        return "Ex_situ"
    return clean(value)


def find_column(names: Iterable[str], candidates: Iterable[str], label: str) -> str:
    names = list(names)
    lower = {name.lower(): name for name in names}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    raise ValueError(f"Cannot identify {label}; columns={names}")


def threshold_label(value: float) -> str:
    return f"PIHAT_ge_{value:.2f}"


def assign_components(
    nodes: list[str],
    edges: list[tuple[str, str]],
    prefix: str,
) -> tuple[list[list[str]], dict[str, str], dict[str, str], dict[str, int]]:
    components = connected_components(nodes, edges)
    sample_component: dict[str, str] = {}
    representative: dict[str, str] = {}
    size: dict[str, int] = {}
    for index, members in enumerate(components, start=1):
        component_id = f"{prefix}_C{index:04d}"
        rep = members[0]
        for sample in members:
            sample_component[sample] = component_id
            representative[sample] = rep
            size[sample] = len(members)
    return components, sample_component, representative, size


def main() -> None:
    cfg = load_config(parse_args().config)
    out = prepare_step_output(cfg, "05_PIHAT_components_and_core")

    metadata_rows = read_tsv(input_path(cfg, "pihat_metadata"))
    if not metadata_rows:
        raise ValueError("PI_HAT metadata is empty")
    sample_col = find_column(metadata_rows[0], ["SampleID", "IID"], "metadata SampleID")
    sample_order: list[str] = []
    metadata: dict[str, dict[str, str]] = {}
    for row in metadata_rows:
        sample = clean(row.get(sample_col))
        if sample == "NA":
            continue
        if sample in metadata:
            raise ValueError(f"Duplicated SampleID in metadata: {sample}")
        material = normalized_material(row.get("Material_type"))
        site = clean(row.get("Ex_situ_site"))
        if material == "Wild":
            site = "NA"
        metadata[sample] = {
            **{key: clean(value) for key, value in row.items()},
            "SampleID": sample,
            "Material_type": material,
            "Ex_situ_site": site,
            "Region": canonical_region(row.get("Region") or row.get("Registered_region")),
        }
        sample_order.append(sample)
    ensure_unique(sample_order, "PI_HAT metadata SampleID")

    pair_path = input_path(cfg, "pihat_pairwise")
    with pair_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header in {pair_path}")
        iid1_col = find_column(reader.fieldnames, ["IID1", "Sample1", "ID1"], "IID1")
        iid2_col = find_column(reader.fieldnames, ["IID2", "Sample2", "ID2"], "IID2")
        pihat_col = find_column(reader.fieldnames, ["PI_HAT", "PIHAT"], "PI_HAT")
        pairs: list[tuple[str, str, float]] = []
        pair_keys: set[tuple[str, str]] = set()
        for row_number, row in enumerate(reader, start=2):
            left = clean(row.get(iid1_col))
            right = clean(row.get(iid2_col))
            if left not in metadata or right not in metadata:
                raise ValueError(f"Pair row {row_number} contains sample absent from metadata: {left}, {right}")
            if left == right:
                continue
            try:
                pihat = float(row[pihat_col])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid PI_HAT at row {row_number}") from exc
            if not 0.0 <= pihat <= 1.0:
                raise ValueError(f"PI_HAT outside [0, 1] at row {row_number}: {pihat}")
            key = tuple(sorted((left, right)))
            if key in pair_keys:
                raise ValueError(f"Duplicated PI_HAT pair: {key[0]}, {key[1]}")
            pair_keys.add(key)
            pairs.append((left, right, pihat))

    parameters = cfg.get("parameters", {})
    thresholds = [float(value) for value in parameters.get("pihat_thresholds", [0.90, 0.80, 0.70, 0.50])]
    operational = float(parameters.get("operational_pihat", 0.80))
    if operational not in thresholds:
        thresholds.append(operational)
    thresholds = sorted(set(thresholds), reverse=True)

    groups: dict[str, list[str]] = {
        "All": list(sample_order),
        "Wild": [s for s in sample_order if metadata[s]["Material_type"] == "Wild"],
        "Ex_situ": [s for s in sample_order if metadata[s]["Material_type"] == "Ex_situ"],
        "Pingshun": [s for s in sample_order if metadata[s]["Ex_situ_site"] == "Pingshun"],
        "Tongguan": [s for s in sample_order if metadata[s]["Ex_situ_site"] == "Tongguan"],
        "Zhendong": [s for s in sample_order if metadata[s]["Ex_situ_site"] == "Zhendong"],
    }

    unsupported_material = sorted(
        sample for sample in sample_order
        if metadata[sample]["Material_type"] not in {"Wild", "Ex_situ"}
    )
    if unsupported_material:
        raise ValueError(f"Samples have unsupported Material_type values: {unsupported_material[:10]}")
    unsupported_sites = sorted(
        sample for sample in groups["Ex_situ"]
        if metadata[sample]["Ex_situ_site"] not in {"Pingshun", "Tongguan", "Zhendong"}
    )
    if unsupported_sites:
        raise ValueError(f"Ex situ samples have unsupported Ex_situ_site values: {unsupported_sites[:10]}")

    configured_groups = {
        "Wild": read_sample_list(input_path(cfg, "wild_samples")),
        "Ex_situ": read_sample_list(input_path(cfg, "exsitu_samples")),
        "Pingshun": read_sample_list(input_path(cfg, "pingshun_samples")),
        "Tongguan": read_sample_list(input_path(cfg, "tongguan_samples")),
        "Zhendong": read_sample_list(input_path(cfg, "zhendong_samples")),
    }
    for group_name, configured_samples in configured_groups.items():
        observed_set = set(groups[group_name])
        configured_set = set(configured_samples)
        if observed_set != configured_set:
            missing = sorted(configured_set - observed_set)
            extra = sorted(observed_set - configured_set)
            raise ValueError(
                f"PI_HAT metadata/sample-list mismatch for {group_name}: "
                f"missing={missing[:10]}, extra={extra[:10]}"
            )
    if set(groups["Wild"]) & set(groups["Ex_situ"]):
        raise ValueError("Wild and ex situ sample sets overlap")
    garden_union = set(groups["Pingshun"]) | set(groups["Tongguan"]) | set(groups["Zhendong"])
    if garden_union != set(groups["Ex_situ"]):
        raise ValueError("Resource-garden sample lists do not exactly partition the ex situ collection")

    global_summary: list[dict[str, Any]] = []
    group_summary: list[dict[str, Any]] = []
    component_rows: list[dict[str, Any]] = []
    status_rows: list[dict[str, Any]] = []
    operational_global_status: dict[str, dict[str, Any]] = {}
    operational_exsitu: tuple[list[list[str]], dict[str, str], dict[str, str], dict[str, int]] | None = None

    for threshold in thresholds:
        label = threshold_label(threshold)
        selected_edges = [(left, right) for left, right, pihat in pairs if pihat >= threshold]
        components, sample_component, representative, component_size = assign_components(
            sample_order, selected_edges, label
        )
        global_summary.append({
            "Threshold": label,
            "PI_HAT_min": threshold,
            "Total_samples": len(sample_order),
            "High_similarity_pair_count": len(selected_edges),
            "Effective_independent_N": len(components),
            "Redundant_component_N": sum(len(c) >= 2 for c in components),
            "Samples_in_redundant_components": sum(len(c) for c in components if len(c) >= 2),
            "Redundant_nonrepresentative_N": len(sample_order) - len(components),
            "Redundancy_rate": (len(sample_order) - len(components)) / len(sample_order),
        })

        for members in components:
            component_id = sample_component[members[0]]
            material_counts = Counter(metadata[s]["Material_type"] for s in members)
            site_counts = Counter(metadata[s]["Ex_situ_site"] for s in members if metadata[s]["Ex_situ_site"] != "NA")
            component_rows.append({
                "Threshold": label,
                "PI_HAT_min": threshold,
                "Component_ID": component_id,
                "Component_size": len(members),
                "Representative_sample_global": representative[members[0]],
                "N_Wild": material_counts.get("Wild", 0),
                "N_Ex_situ": material_counts.get("Ex_situ", 0),
                "N_Pingshun": site_counts.get("Pingshun", 0),
                "N_Tongguan": site_counts.get("Tongguan", 0),
                "N_Zhendong": site_counts.get("Zhendong", 0),
                "Samples": ";".join(members),
            })
        for sample in sample_order:
            row = {
                "Threshold": label,
                "PI_HAT_min": threshold,
                "SampleID": sample,
                "Global_component_ID": sample_component[sample],
                "Global_component_size": component_size[sample],
                "Global_representative_sample": representative[sample],
                "Is_global_representative": sample == representative[sample],
                "Material_type": metadata[sample]["Material_type"],
                "Ex_situ_site": metadata[sample]["Ex_situ_site"],
                "Region": metadata[sample]["Region"],
            }
            status_rows.append(row)
            if math_is_close(threshold, operational):
                operational_global_status[sample] = row

        # Historical group summary based on global component membership plus a
        # transparent induced-subgraph summary for each group.
        for group_name, members in groups.items():
            global_effective = len({sample_component[sample] for sample in members})
            member_set = set(members)
            induced_edges = [(l, r) for l, r in selected_edges if l in member_set and r in member_set]
            induced_components = connected_components(members, induced_edges) if members else []
            for scope, effective in (
                ("global_component_membership", global_effective),
                ("induced_subgraph", len(induced_components)),
            ):
                group_summary.append({
                    "Threshold": label,
                    "PI_HAT_min": threshold,
                    "Graph_scope": scope,
                    "Group": group_name,
                    "Nominal_N": len(members),
                    "Effective_independent_N": effective,
                    "Redundant_reduction_N": len(members) - effective,
                    "Redundancy_rate": (len(members) - effective) / len(members) if members else 0.0,
                })

        exsitu_members = groups["Ex_situ"]
        exsitu_set = set(exsitu_members)
        exsitu_edges = [(l, r) for l, r in selected_edges if l in exsitu_set and r in exsitu_set]
        exsitu_assignment = assign_components(exsitu_members, exsitu_edges, f"EXSITU_{label}")
        if math_is_close(threshold, operational):
            operational_exsitu = exsitu_assignment

    if operational_exsitu is None:
        raise RuntimeError("Operational PI_HAT threshold was not processed")
    ex_components, ex_component, ex_representative, ex_size = operational_exsitu

    core_rows: list[dict[str, Any]] = []
    low_rows: list[dict[str, Any]] = []
    status_exsitu_rows: list[dict[str, Any]] = []
    for sample in groups["Ex_situ"]:
        keep = sample == ex_representative[sample]
        row = {
            "SampleID": sample,
            "Operational_threshold": operational,
            "Ex_situ_component_ID": ex_component[sample],
            "Ex_situ_component_size": ex_size[sample],
            "Representative_sample": ex_representative[sample],
            "Is_representative": keep,
            "Management_class": "Core_keep" if keep else "Lower_priority_redundant",
            "Representative_rule": "First_accession_in_ascending_lexicographic_SampleID_order",
            "Material_type": "Ex_situ",
            "Ex_situ_site": metadata[sample]["Ex_situ_site"],
            "Region": metadata[sample]["Region"],
        }
        status_exsitu_rows.append(row)
        (core_rows if keep else low_rows).append(row)

    write_tsv(
        out / "01.PIHAT_global_threshold_summary.tsv",
        global_summary,
        list(global_summary[0]),
    )
    write_tsv(out / "02.PIHAT_group_threshold_summary.tsv", group_summary, list(group_summary[0]))
    write_tsv(out / "03.global_component_summary.tsv.gz", component_rows, list(component_rows[0]))
    write_tsv(out / "04.global_sample_component_status.tsv.gz", status_rows, list(status_rows[0]))
    write_tsv(out / "05.operational_exsitu_component_status.tsv", status_exsitu_rows, list(status_exsitu_rows[0]))
    write_tsv(out / "06.deredundant_exsitu_core.tsv", core_rows, list(core_rows[0]))
    write_tsv(out / "07.lower_priority_redundant_exsitu.tsv", low_rows, list(low_rows[0]))
    write_tsv(
        out / "08.operational_global_sample_status.tsv",
        (operational_global_status[sample] for sample in sample_order),
        list(next(iter(operational_global_status.values()))),
    )

    site_core = Counter(row["Ex_situ_site"] for row in core_rows)
    site_low = Counter(row["Ex_situ_site"] for row in low_rows)
    observed = {
        "Total_samples": len(sample_order),
        "Wild_samples": len(groups["Wild"]),
        "Ex_situ_samples": len(groups["Ex_situ"]),
        "Exsitu_core_size": len(core_rows),
        "Exsitu_low_priority_size": len(low_rows),
        "Pingshun_core_size": site_core["Pingshun"],
        "Tongguan_core_size": site_core["Tongguan"],
        "Zhendong_core_size": site_core["Zhendong"],
    }
    write_key_value(out / "09.internal_check.tsv", [
        ("Pairwise_rows_loaded", len(pairs)),
        ("Configured_sample_lists_match_metadata", True),
        ("Operational_threshold", operational),
        ("Operational_exsitu_component_N", len(ex_components)),
        *observed.items(),
        ("Pingshun_lower_priority_size", site_low["Pingshun"]),
        ("Tongguan_lower_priority_size", site_low["Tongguan"]),
        ("Zhendong_lower_priority_size", site_low["Zhendong"]),
    ])
    print(f"[OK] PI_HAT components and ex situ core complete: {out}")


def math_is_close(left: float, right: float) -> bool:
    return abs(left - right) < 1e-12


if __name__ == "__main__":
    main()
