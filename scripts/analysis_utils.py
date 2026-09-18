#!/usr/bin/env python3
"""Shared utilities for the Forsythia resource-garden genomic assessment workflow.

The workflow intentionally uses only the Python standard library.  Large VCF
subsetting in the ALT-coverage step requires a working ``bcftools`` executable.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from itertools import zip_longest
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

MAF_BIN_ORDER = (
    "MAF_le_0.01",
    "0.01_lt_MAF_le_0.05",
    "0.05_lt_MAF_le_0.10",
    "MAF_gt_0.10",
)
THRESHOLD_ORDER = ("all", "AC_ge_2", "AC_ge_3")
THRESHOLDS = {"all": 1, "AC_ge_2": 2, "AC_ge_3": 3}

REGION_ALIASES = {
    "qinling": "Qinling Mountains",
    "qinling mountain": "Qinling Mountains",
    "qinling mountains": "Qinling Mountains",
    "southern loess": "Southern Loess Plateau",
    "southern loess plateau": "Southern Loess Plateau",
    "lvliang": "Lvliang Mountain",
    "lvliang mountain": "Lvliang Mountain",
    "lvliang mountains": "Lvliang Mountain",
    "taiyue": "Taiyue Mountains",
    "taiyue mountain": "Taiyue Mountains",
    "taiyue mountains": "Taiyue Mountains",
    "taihang": "Taihang Mountains",
    "taihang mountain": "Taihang Mountains",
    "taihang mountains": "Taihang Mountains",
    "zhongtiao": "Zhongtiao Mountain",
    "zhongtiaoshan": "Zhongtiao Mountain",
    "zhongtiao mountain": "Zhongtiao Mountain",
    "zhongtiao mountains": "Zhongtiao Mountain",
}


def package_root() -> Path:
    return Path(__file__).resolve().parents[1]


def default_config_path() -> Path:
    return package_root() / "config" / "analysis_config.json"


def load_config(path: str | Path) -> dict[str, Any]:
    """Load and normalize the public JSON configuration.

    The cleaned scripts use a flat set of historical input aliases internally,
    while the public configuration groups sample lists and AC/AN tables for
    readability.  This function is the single compatibility layer between the
    two representations.  Relative project paths are resolved against the
    current working directory; ``run_pipeline.sh`` changes to the requested
    project root before invoking any step.
    """
    path = Path(path).expanduser().resolve()
    with path.open("r", encoding="utf-8") as handle:
        cfg = json.load(handle)
    if not isinstance(cfg, dict):
        raise ValueError(f"Configuration must be a JSON object: {path}")

    cfg["_config_path"] = str(path)
    cfg["_config_dir"] = str(path.parent)
    cfg["_package_root"] = str(package_root())
    root = Path(cfg.get("project_root", ".")).expanduser()
    cfg["_project_root"] = str(root.resolve())

    # Public key is singular; retain the historical plural alias used by the
    # consolidated analysis scripts.
    if "outputs_root" not in cfg:
        cfg["outputs_root"] = cfg.get("output_root", "genomic_audit_reproducibility_results")

    inputs = cfg.setdefault("inputs", {})
    if not isinstance(inputs, dict):
        raise ValueError("config.inputs must be a JSON object")
    sample_lists = inputs.get("sample_lists", {}) or {}
    acan = inputs.get("acan", {}) or {}
    aliases = {
        "wild_samples": sample_lists.get("wild"),
        "exsitu_samples": sample_lists.get("ex_situ"),
        "pingshun_samples": sample_lists.get("pingshun"),
        "tongguan_samples": sample_lists.get("tongguan"),
        "zhendong_samples": sample_lists.get("zhendong"),
        "wild_acan": acan.get("wild"),
        "exsitu_acan": acan.get("ex_situ"),
        "pingshun_acan": acan.get("pingshun"),
        "tongguan_acan": acan.get("tongguan"),
        "zhendong_acan": acan.get("zhendong"),
        "regional_acan": inputs.get("region_acan"),
        "pihat_pairwise": inputs.get("pairwise_pihat"),
        "pihat_metadata": inputs.get("pihat_metadata") or inputs.get("metadata"),
        "roh_summary": inputs.get("roh_sample_summary"),
        "full_snp_vcf": inputs.get("pass_snp_vcf", inputs.get("full_filtered_biallelic_snp_vcf")),
    }
    for key, value in aliases.items():
        if key not in inputs and value is not None:
            inputs[key] = value

    parameters = cfg.setdefault("parameters", {})
    if not isinstance(parameters, dict):
        raise ValueError("config.parameters must be a JSON object")
    if "operational_pihat" not in parameters:
        parameters["operational_pihat"] = parameters.get("operational_pihat_threshold", 0.80)
    if "greedy_weights" not in parameters:
        parameters["greedy_weights"] = parameters.get("greedy_maf_weights", {})
    random_control = parameters.get("random_control", {}) or {}
    parameters.setdefault("random_seed", random_control.get("seed", 20260620))
    parameters.setdefault("random_permutations", random_control.get("n_permutations", 1000))
    parameters.setdefault(
        "legacy_random_stream_burn_all_wild",
        random_control.get("legacy_rng_compatibility", True),
    )
    parameters.setdefault("ploidy", 2)

    package_files = cfg.setdefault("package_files", {})
    local_quota = path.parent / "regional_quotas.tsv"
    package_files.setdefault(
        "regional_quotas",
        str(local_quota.resolve()) if local_quota.is_file() else "config/regional_quotas.tsv",
    )
    return cfg


def project_path(cfg: Mapping[str, Any], value: str | Path) -> Path:
    value = Path(value)
    if value.is_absolute():
        return value
    return Path(cfg["_project_root"]) / value


def input_path(cfg: Mapping[str, Any], key: str) -> Path:
    try:
        value = cfg["inputs"][key]
    except KeyError as exc:
        raise KeyError(f"Missing config input key: {key}") from exc
    return project_path(cfg, value)


def package_path(cfg: Mapping[str, Any], value: str | Path) -> Path:
    value = Path(value)
    if value.is_absolute():
        return value
    return Path(cfg["_package_root"]) / value


def outputs_root(cfg: Mapping[str, Any]) -> Path:
    root = project_path(
        cfg,
        cfg.get("outputs_root", cfg.get("output_root", "genomic_audit_reproducibility_results")),
    )
    root.mkdir(parents=True, exist_ok=True)
    return root


def step_output_dir(cfg: Mapping[str, Any], step_name: str) -> Path:
    out = outputs_root(cfg) / step_name
    out.mkdir(parents=True, exist_ok=True)
    return out


def ensure_file(path: str | Path, label: str | None = None) -> Path:
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"Missing {label or 'file'}: {p}")
    return p


def ensure_unique(values: Sequence[str], label: str) -> None:
    seen: set[str] = set()
    duplicates: list[str] = []
    for value in values:
        if value in seen:
            duplicates.append(value)
        seen.add(value)
    if duplicates:
        preview = ", ".join(sorted(set(duplicates))[:10])
        raise ValueError(f"Duplicated {label}: {preview}")


def open_text(path: str | Path, mode: str = "rt"):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, mode, encoding=None if "b" in mode else "utf-8", newline="")
    return p.open(mode, encoding=None if "b" in mode else "utf-8", newline="")


def read_sample_list(path: str | Path) -> list[str]:
    ensure_file(path, "sample list")
    with Path(path).open("r", encoding="utf-8") as handle:
        samples = [line.strip() for line in handle if line.strip()]
    ensure_unique(samples, f"sample IDs in {path}")
    if not samples:
        raise ValueError(f"No samples found in {path}")
    return samples


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    ensure_file(path)
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")
        return list(reader)


def iter_tsv(path: str | Path) -> Iterator[dict[str, str]]:
    ensure_file(path)
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")
        yield from reader


def write_tsv(path: str | Path, rows: Iterable[Mapping[str, Any]], fieldnames: Sequence[str]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open_text(p, "wt") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_value(row.get(key, "")) for key in fieldnames})


def write_key_value(path: str | Path, items: Iterable[tuple[str, Any]]) -> None:
    write_tsv(path, ({"Item": k, "Value": v} for k, v in items), ["Item", "Value"])


def format_value(value: Any) -> Any:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.17g}"
    return value


def require_columns(fieldnames: Sequence[str] | None, required: Sequence[str], label: str) -> None:
    names = set(fieldnames or [])
    missing = [name for name in required if name not in names]
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")


def canonical_region(value: str | None) -> str:
    text = "" if value is None else str(value)
    text = " ".join(text.replace("_", " ").split())
    if not text:
        return "NA"
    return REGION_ALIASES.get(text.lower(), text)


def maf_bin(maf: float) -> str:
    if maf <= 0.01:
        return MAF_BIN_ORDER[0]
    if maf <= 0.05:
        return MAF_BIN_ORDER[1]
    if maf <= 0.10:
        return MAF_BIN_ORDER[2]
    return MAF_BIN_ORDER[3]


@dataclass(frozen=True)
class ACANRow:
    chrom: str
    pos: str
    alt_ac: int
    an: int


def parse_acan_line(line: str, label: str, line_number: int) -> ACANRow:
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 4:
        raise ValueError(f"Malformed {label} AC/AN row {line_number}: expected >=4 columns")
    try:
        return ACANRow(parts[0], parts[1], int(parts[2]), int(parts[3]))
    except ValueError as exc:
        raise ValueError(f"Non-integer AC/AN at {label} row {line_number}: {line.rstrip()}") from exc


def _looks_like_acan_header(line: str) -> bool:
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 4:
        return False
    try:
        int(parts[2])
        int(parts[3])
        return False
    except ValueError:
        tokens = {part.strip().upper().lstrip("#") for part in parts[:4]}
        return bool(tokens & {"CHROM", "POS", "AC", "ALT_AC", "AN"})


def iter_acan_rows(path: str | Path, label: str) -> Iterator[ACANRow]:
    """Stream one AC/AN table, permitting one optional header and blank lines."""
    seen_data = False
    with open_text(ensure_file(path, f"{label} AC/AN"), "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            if not seen_data and _looks_like_acan_header(line):
                seen_data = True  # header has now been consumed; a second one is invalid
                continue
            row = parse_acan_line(line, label, line_number)
            if row.alt_ac < 0 or row.an < 0 or row.alt_ac > row.an:
                raise ValueError(
                    f"Invalid AC/AN at {label} row {line_number}: AC={row.alt_ac}, AN={row.an}"
                )
            seen_data = True
            yield row


def iter_aligned_acan(paths: Mapping[str, str | Path]) -> Iterator[dict[str, ACANRow]]:
    """Stream AC/AN files in lockstep and fail on length/coordinate mismatch.

    Each file may contain a single optional header.  Alignment is checked after
    header removal rather than by physical line number.
    """
    labels = list(paths)
    iterators = [iter_acan_rows(paths[label], label) for label in labels]
    for row_number, rows_tuple in enumerate(zip_longest(*iterators), start=1):
        if any(row is None for row in rows_tuple):
            ended = [labels[i] for i, row in enumerate(rows_tuple) if row is None]
            raise ValueError(f"AC/AN files have different data-row counts near row {row_number}; ended: {ended}")
        rows = {label: row for label, row in zip(labels, rows_tuple) if row is not None}
        coordinates = {(row.chrom, row.pos) for row in rows.values()}
        if len(coordinates) != 1:
            detail = ", ".join(f"{label}={row.chrom}:{row.pos}" for label, row in rows.items())
            raise ValueError(f"AC/AN coordinate mismatch at data row {row_number}: {detail}")
        yield rows


def wild_minor_from_acan(wild: ACANRow) -> tuple[str, int, float] | None:
    if wild.an <= 0:
        return None
    ref_count = wild.an - wild.alt_ac
    if wild.alt_ac <= 0 or ref_count <= 0:
        return None
    # Ties are assigned to ALT, matching the original workflow.
    if wild.alt_ac <= ref_count:
        minor = "ALT"
        count = wild.alt_ac
    else:
        minor = "REF"
        count = ref_count
    return minor, count, count / wild.an


def target_allele_count(row: ACANRow, target: str) -> int:
    if target == "ALT":
        return row.alt_ac
    if target == "REF":
        return row.an - row.alt_ac
    raise ValueError(f"Unknown target allele label: {target}")


def parse_gt(gt: str) -> list[int]:
    if not gt or gt in {".", "./.", ".|."}:
        return []
    alleles: list[int] = []
    for token in gt.replace("|", "/").split("/"):
        if token == ".":
            return []
        try:
            alleles.append(int(token))
        except ValueError:
            return []
    return alleles


def target_dosage(gt: str, target: str) -> int:
    alleles = parse_gt(gt)
    if len(alleles) != 2:
        return 0
    code = 1 if target == "ALT" else 0
    return sum(1 for allele in alleles if allele == code)


def load_missing_targets(path: str | Path) -> tuple[dict[tuple[str, str], dict[str, Any]], list[tuple[str, str]]]:
    targets: dict[tuple[str, str], dict[str, Any]] = {}
    order: list[tuple[str, str]] = []
    with open_text(ensure_file(path, "missing-allele table"), "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(
            reader.fieldnames,
            ["CHROM", "POS", "Wild_minor_allele", "Wild_minor_count", "Wild_minor_MAF"],
            "missing-allele table",
        )
        for index, row in enumerate(reader):
            key = (row["CHROM"], row["POS"])
            if key in targets:
                raise ValueError(f"Duplicate missing target: {key[0]}:{key[1]}")
            minor_count = int(float(row["Wild_minor_count"]))
            maf = float(row["Wild_minor_MAF"])
            targets[key] = {
                "idx": index,
                "minor": row["Wild_minor_allele"],
                "minor_count": minor_count,
                "maf": maf,
                "maf_bin": maf_bin(maf),
                "row": row,
            }
            order.append(key)
    if not targets:
        raise ValueError(f"No missing targets loaded from {path}")
    return targets, order


def read_vcf_header_samples(vcf_path: str | Path) -> list[str]:
    with open_text(ensure_file(vcf_path, "VCF"), "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                ensure_unique(samples, f"VCF samples in {vcf_path}")
                return samples
    raise ValueError(f"No #CHROM header found in {vcf_path}")


@dataclass
class TargetBitsetResult:
    sample_bits: dict[str, int]
    called_sites: dict[str, int]
    processed_targets: int
    unmatched_vcf_sites: int
    missing_targets: list[tuple[str, str]]


def build_target_bitsets_from_vcf(
    vcf_path: str | Path,
    selected_samples: Sequence[str],
    targets: Mapping[tuple[str, str], Mapping[str, Any]],
) -> TargetBitsetResult:
    """Build one integer carrier bitset per sample from a missing-site VCF.

    The VCF may contain all samples.  Only ``selected_samples`` are examined.
    ``called_sites`` counts non-missing GT fields at target sites and is used as
    the Greedy call-rate tie breaker.
    """
    selected_samples = list(selected_samples)
    ensure_unique(selected_samples, "selected sample IDs")
    sample_bits = {sample: 0 for sample in selected_samples}
    called_sites = {sample: 0 for sample in selected_samples}
    matched: set[tuple[str, str]] = set()
    unmatched = 0

    with open_text(ensure_file(vcf_path, "missing-site VCF"), "rt") as handle:
        sample_names: list[str] | None = None
        selected_indices: list[tuple[str, int]] = []
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_names = line.rstrip("\n").split("\t")[9:]
                sample_to_index = {sample: i for i, sample in enumerate(sample_names)}
                missing_samples = [sample for sample in selected_samples if sample not in sample_to_index]
                if missing_samples:
                    raise ValueError(
                        f"{len(missing_samples)} requested samples are absent from VCF; "
                        f"examples: {', '.join(missing_samples[:10])}"
                    )
                selected_indices = [(sample, sample_to_index[sample]) for sample in selected_samples]
                continue
            if sample_names is None:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 + len(sample_names):
                raise ValueError(f"Malformed VCF row at {fields[0] if fields else '?'}")
            key = (fields[0], fields[1])
            if key not in targets:
                unmatched += 1
                continue
            if key in matched:
                raise ValueError(f"Duplicate target site in VCF: {key[0]}:{key[1]}")
            matched.add(key)
            target = targets[key]
            fmt = fields[8].split(":")
            if "GT" not in fmt:
                raise ValueError(f"GT is absent at {key[0]}:{key[1]}")
            gt_index = fmt.index("GT")
            genotypes = fields[9:]
            bit = 1 << int(target["idx"])
            for sample, index in selected_indices:
                values = genotypes[index].split(":")
                if gt_index >= len(values):
                    continue
                gt = values[gt_index]
                alleles = parse_gt(gt)
                if len(alleles) != 2:
                    continue
                if any(a not in (0, 1) for a in alleles):
                    raise ValueError(f"Non-biallelic GT at {key}")
                called_sites[sample] += 1
                code = 1 if target["minor"] == "ALT" else 0
                if any(allele == code for allele in alleles):
                    sample_bits[sample] |= bit

    missing_targets = sorted(set(targets) - matched)
    return TargetBitsetResult(
        sample_bits=sample_bits,
        called_sites=called_sites,
        processed_targets=len(matched),
        unmatched_vcf_sites=unmatched,
        missing_targets=missing_targets,
    )


def target_masks(targets: Mapping[tuple[str, str], Mapping[str, Any]]) -> dict[str, int]:
    masks = {"All": 0, **{name: 0 for name in MAF_BIN_ORDER}}
    for target in targets.values():
        bit = 1 << int(target["idx"])
        masks["All"] |= bit
        masks[str(target["maf_bin"])] |= bit
    return masks


def recovery_from_bits(bits: int, masks: Mapping[str, int]) -> dict[str, dict[str, float | int]]:
    result: dict[str, dict[str, float | int]] = {}
    for label, mask in masks.items():
        total = mask.bit_count()
        recovered = (bits & mask).bit_count()
        result[label] = {
            "Total": total,
            "Recovered": recovered,
            "Recovery_rate": recovered / total if total else math.nan,
        }
    return result


def merge_sample_bits(samples: Iterable[str], sample_bits: Mapping[str, int]) -> int:
    merged = 0
    for sample in samples:
        merged |= sample_bits[sample]
    return merged


def quantile(values: Sequence[float], q: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = q * (len(ordered) - 1)
    low = int(math.floor(position))
    high = int(math.ceil(position))
    fraction = position - low
    return ordered[low] * (1.0 - fraction) + ordered[high] * fraction


class UnionFind:
    def __init__(self, items: Iterable[str]):
        self.parent = {item: item for item in items}
        self.rank = {item: 0 for item in items}

    def find(self, item: str) -> str:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left == root_right:
            return
        if self.rank[root_left] < self.rank[root_right]:
            self.parent[root_left] = root_right
        elif self.rank[root_left] > self.rank[root_right]:
            self.parent[root_right] = root_left
        else:
            self.parent[root_right] = root_left
            self.rank[root_left] += 1


def connected_components(nodes: Sequence[str], edges: Iterable[tuple[str, str]]) -> list[list[str]]:
    nodes = list(nodes)
    ensure_unique(nodes, "graph nodes")
    node_set = set(nodes)
    uf = UnionFind(nodes)
    for left, right in edges:
        if left in node_set and right in node_set and left != right:
            uf.union(left, right)
    grouped: dict[str, list[str]] = defaultdict(list)
    for node in nodes:
        grouped[uf.find(node)].append(node)
    components = [sorted(component) for component in grouped.values()]
    return sorted(components, key=lambda component: component[0])


def run_bcftools_view(vcf: str | Path, sample_file: str | Path) -> tuple[subprocess.Popen[str], list[str]]:
    """Return a streaming ``bcftools view`` process and subset sample order."""
    ensure_file(vcf, "VCF")
    ensure_file(sample_file, "sample list")
    header = subprocess.check_output(
        ["bcftools", "view", "-h", "-S", str(sample_file), str(vcf)],
        text=True,
    )
    samples: list[str] | None = None
    for line in header.splitlines():
        if line.startswith("#CHROM"):
            samples = line.split("\t")[9:]
            break
    if samples is None:
        raise ValueError("Cannot find #CHROM header in bcftools output")
    process = subprocess.Popen(
        ["bcftools", "view", "-H", "-S", str(sample_file), str(vcf)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    return process, samples


def prepare_step_output(cfg, step_name):
    """Fail on an existing nonempty step; never silently overwrite results."""
    out = step_output_dir(cfg, step_name)
    if any(out.iterdir()):
        raise FileExistsError(f"Output step is not empty: {out}; select a new output_root")
    (out / "configuration_used.json").write_text(json.dumps(cfg, indent=2) + "\n")
    return out
