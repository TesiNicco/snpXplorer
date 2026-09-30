#!/usr/bin/env python3
"""Extract an all-variant LD matrix for a set of SNPXplorer haplotypes.

This is an offline utility. It does not modify or participate in the web app.
Haplotype membership is read from the SNPXplorer haplotype table and LD is read
from the chromosome-specific SQLite databases used by SNPXplorer.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import re
import shutil
import sqlite3
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence


LOGGER = logging.getLogger("haplotype_region_ld")
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_HAPLOTYPE_TABLE = (
    REPO_ROOT
    / "Data/databases/haplotypes/ld_clusters_with_gwas_AI_latest.tsv.gz"
)
DEFAULT_LD_DIRECTORY = REPO_ROOT / "Data/databases/LD_db"

HAPLOTYPE_RE = re.compile(
    r"^(?:chr)?(?P<chrom>[0-9]+|X|Y|MT|M):(?P<start>[0-9,]+)-(?P<end>[0-9,]+)$",
    re.IGNORECASE,
)
VARIANT_RE = re.compile(
    r"^(?:chr)?(?P<chrom>[0-9]+|X|Y|MT|M):(?P<pos>[0-9]+):"
    r"(?P<ref>[^:]+):(?P<alt>[^:]+)$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class HaplotypeRequest:
    label: str
    haplotype_id: str
    chrom: str
    start: int
    end: int


@dataclass(frozen=True)
class HaplotypeRecord:
    request: HaplotypeRequest
    kb: str
    declared_n_snps: int
    variants: tuple[str, ...]


@dataclass
class VariantRecord:
    variant_id: str
    chrom: str
    position: int
    ref: str
    alt: str
    haplotype_labels: set[str]
    haplotype_ids: set[str]
    db_id: int | None = None
    db_uniq: str | None = None
    mapping_status: str = "not_checked"


def normalize_chromosome(value: str) -> str:
    chrom = re.sub(r"^chr", "", value.strip(), flags=re.IGNORECASE).upper()
    return "MT" if chrom == "M" else chrom


def chromosome_sort_key(chrom: str) -> tuple[int, str]:
    if chrom.isdigit():
        return (int(chrom), "")
    return ({"X": 23, "Y": 24, "MT": 25}.get(chrom, 10_000), chrom)


def parse_haplotype_id(value: str) -> tuple[str, int, int]:
    match = HAPLOTYPE_RE.fullmatch(value.strip())
    if not match:
        raise ValueError(
            f"Invalid haplotype ID {value!r}; expected CHR:START-END, "
            "for example 7:12194394-12246599"
        )
    chrom = normalize_chromosome(match.group("chrom"))
    start = int(match.group("start").replace(",", ""))
    end = int(match.group("end").replace(",", ""))
    if start > end:
        raise ValueError(f"Haplotype start exceeds end in {value!r}")
    return chrom, start, end


def make_request(label: str, haplotype_id: str) -> HaplotypeRequest:
    clean_label = label.strip()
    clean_id = haplotype_id.strip()
    if not clean_label:
        raise ValueError("Haplotype labels cannot be empty")
    chrom, start, end = parse_haplotype_id(clean_id)
    return HaplotypeRequest(clean_label, clean_id, chrom, start, end)


def load_requests(
    cli_values: Sequence[str] | None, manifest: Path | None
) -> list[HaplotypeRequest]:
    requests: list[HaplotypeRequest] = []

    for value in cli_values or []:
        if "=" not in value:
            raise ValueError(
                f"Invalid --haplotype value {value!r}; expected LABEL=CHR:START-END"
            )
        label, haplotype_id = value.split("=", 1)
        requests.append(make_request(label, haplotype_id))

    if manifest is not None:
        with manifest.open("rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"label", "haplotype_id"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValueError(
                    f"Manifest {manifest} must contain tab-separated columns: "
                    "label and haplotype_id"
                )
            for row_number, row in enumerate(reader, start=2):
                try:
                    requests.append(make_request(row["label"], row["haplotype_id"]))
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"Manifest line {row_number}: {exc}") from exc

    if not requests:
        raise ValueError("Provide at least one --haplotype or --manifest")

    labels = [request.label for request in requests]
    if len(labels) != len(set(labels)):
        raise ValueError("Haplotype labels must be unique")
    ids = [request.haplotype_id for request in requests]
    if len(ids) != len(set(ids)):
        raise ValueError("Haplotype IDs must be unique")
    chromosomes = {request.chrom for request in requests}
    if len(chromosomes) != 1:
        raise ValueError("All haplotypes must be on the same chromosome")
    return requests


def iter_haplotype_lines(path: Path, request: HaplotypeRequest) -> Iterator[str]:
    """Yield tabix-region lines, falling back to a sequential gzip scan."""
    index_path = Path(f"{path}.tbi")
    tabix = shutil.which("tabix")
    if tabix and index_path.exists():
        region = f"{request.chrom}:{request.start}-{request.end}"
        result = subprocess.run(
            [tabix, str(path), region],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            yield from result.stdout.splitlines()
            return
        LOGGER.warning(
            "Tabix lookup failed for %s (%s); scanning the file instead",
            request.haplotype_id,
            result.stderr.strip(),
        )

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t", 3)
            if len(fields) < 3:
                continue
            try:
                chrom = normalize_chromosome(fields[0])
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            if chrom == request.chrom and start <= request.end and end >= request.start:
                yield line.rstrip("\n")


def load_haplotype_record(path: Path, request: HaplotypeRequest) -> HaplotypeRecord:
    matches: list[list[str]] = []
    for line in iter_haplotype_lines(path, request):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        try:
            chrom = normalize_chromosome(fields[0])
            start = int(fields[1])
            end = int(fields[2])
        except ValueError:
            continue
        if (chrom, start, end) == (request.chrom, request.start, request.end):
            matches.append(fields)

    if not matches:
        raise ValueError(
            f"Haplotype {request.label} ({request.haplotype_id}) was not found in {path}"
        )
    if len(matches) > 1:
        raise ValueError(
            f"Haplotype {request.haplotype_id} has {len(matches)} exact records in {path}"
        )

    fields = matches[0]
    variants = tuple(value for value in fields[5].split("|") if value)
    try:
        declared_n_snps = int(fields[4])
    except ValueError as exc:
        raise ValueError(f"Invalid NSNPS for haplotype {request.haplotype_id}") from exc
    if declared_n_snps != len(variants):
        LOGGER.warning(
            "%s declares %d SNPs but contains %d variant identifiers",
            request.haplotype_id,
            declared_n_snps,
            len(variants),
        )
    return HaplotypeRecord(request, fields[3], declared_n_snps, variants)


def parse_variant_id(variant_id: str) -> tuple[str, int, str, str]:
    match = VARIANT_RE.fullmatch(variant_id.strip())
    if not match:
        raise ValueError(
            f"Invalid variant identifier {variant_id!r}; expected CHR:POS:REF:ALT"
        )
    return (
        normalize_chromosome(match.group("chrom")),
        int(match.group("pos")),
        match.group("ref").upper(),
        match.group("alt").upper(),
    )


def combine_variants(records: Sequence[HaplotypeRecord]) -> list[VariantRecord]:
    variants: dict[str, VariantRecord] = {}
    for record in records:
        for variant_id in record.variants:
            chrom, position, ref, alt = parse_variant_id(variant_id)
            if chrom != record.request.chrom:
                raise ValueError(
                    f"Variant {variant_id} is not on haplotype chromosome "
                    f"{record.request.chrom}"
                )
            if variant_id not in variants:
                variants[variant_id] = VariantRecord(
                    variant_id=variant_id,
                    chrom=chrom,
                    position=position,
                    ref=ref,
                    alt=alt,
                    haplotype_labels=set(),
                    haplotype_ids=set(),
                )
            variants[variant_id].haplotype_labels.add(record.request.label)
            variants[variant_id].haplotype_ids.add(record.request.haplotype_id)

    return sorted(
        variants.values(),
        key=lambda item: (
            chromosome_sort_key(item.chrom),
            item.position,
            item.ref,
            item.alt,
        ),
    )


def chunked(values: Sequence[int], size: int = 900) -> Iterator[Sequence[int]]:
    for offset in range(0, len(values), size):
        yield values[offset : offset + size]


def open_read_only_database(path: Path) -> sqlite3.Connection:
    # Read-only mode prevents accidental changes while still allowing SQLite to
    # see committed data in an existing WAL file. ``immutable=1`` is avoided
    # because it can ignore WAL contents.
    uri = f"file:{path.resolve()}?mode=ro"
    connection = sqlite3.connect(uri, uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def map_variants_to_database(
    connection: sqlite3.Connection, variants: Sequence[VariantRecord]
) -> None:
    by_position: dict[int, list[sqlite3.Row]] = defaultdict(list)
    positions = sorted({variant.position for variant in variants})
    for position_chunk in chunked(positions):
        placeholders = ",".join("?" for _ in position_chunk)
        query = f"SELECT id, uniq, pos FROM variants WHERE pos IN ({placeholders})"
        for row in connection.execute(query, tuple(position_chunk)):
            by_position[int(row["pos"])].append(row)

    for variant in variants:
        candidates: list[sqlite3.Row] = []
        for row in by_position.get(variant.position, []):
            parts = str(row["uniq"]).split(":")
            if len(parts) != 3:
                continue
            _, allele1, allele2 = parts
            observed = (allele1.upper(), allele2.upper())
            if observed in {(variant.ref, variant.alt), (variant.alt, variant.ref)}:
                candidates.append(row)

        if not candidates:
            variant.mapping_status = "not_found"
            continue
        candidates.sort(key=lambda row: int(row["id"]))
        selected = candidates[0]
        variant.db_id = int(selected["id"])
        variant.db_uniq = str(selected["uniq"])
        variant.mapping_status = "matched" if len(candidates) == 1 else "ambiguous_first_id"


def query_stored_ld(
    connection: sqlite3.Connection, variants: Sequence[VariantRecord]
) -> dict[tuple[int, int], tuple[float, int]]:
    selected_ids = {variant.db_id for variant in variants if variant.db_id is not None}
    pairs: dict[tuple[int, int], tuple[float, int]] = {}
    for db_id in sorted(selected_ids):
        rows = connection.execute(
            "SELECT v1, v2, r2_milli, dist_bp FROM ld WHERE v1 = ? "
            "UNION ALL "
            "SELECT v1, v2, r2_milli, dist_bp FROM ld WHERE v2 = ?",
            (db_id, db_id),
        )
        for row in rows:
            v1 = int(row["v1"])
            v2 = int(row["v2"])
            if v1 not in selected_ids or v2 not in selected_ids or v1 == v2:
                continue
            key = (min(v1, v2), max(v1, v2))
            value = (int(row["r2_milli"]) / 1000.0, int(row["dist_bp"]))
            previous = pairs.get(key)
            if previous is None or value[0] > previous[0]:
                pairs[key] = value
    return pairs


def output_paths(prefix: Path) -> dict[str, Path]:
    return {
        "matrix": Path(f"{prefix}.ld_matrix.tsv.gz"),
        "mask": Path(f"{prefix}.ld_observed_mask.tsv.gz"),
        "pairs": Path(f"{prefix}.ld_observed_pairs.tsv.gz"),
        "heatmap": Path(f"{prefix}.ld_heatmap.png"),
        "haplotype_heatmap": Path(f"{prefix}.ld_heatmap_by_haplotype.png"),
        "position_heatmap": Path(f"{prefix}.ld_heatmap_genomic_position.png"),
        "variants": Path(f"{prefix}.variants.tsv"),
        "haplotypes": Path(f"{prefix}.haplotypes.tsv"),
        "qc": Path(f"{prefix}.qc.json"),
    }


def write_square_matrix(
    path: Path,
    variants: Sequence[VariantRecord],
    pairs: dict[tuple[int, int], tuple[float, int]],
    missing_value: str,
    observed_mask: bool,
) -> None:
    missing = "0" if missing_value == "zero" else "NA"
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["variant_id", *(variant.variant_id for variant in variants)])
        for row_index, row_variant in enumerate(variants):
            values: list[str] = []
            for column_index, column_variant in enumerate(variants):
                if row_index == column_index:
                    values.append("1")
                    continue
                if row_variant.db_id is None or column_variant.db_id is None:
                    values.append("0" if observed_mask else missing)
                    continue
                key = tuple(sorted((row_variant.db_id, column_variant.db_id)))
                if key in pairs:
                    values.append("1" if observed_mask else f"{pairs[key][0]:.3f}")
                else:
                    values.append("0" if observed_mask else missing)
            writer.writerow([row_variant.variant_id, *values])


def matrix_value(
    row_variant: VariantRecord,
    column_variant: VariantRecord,
    pairs: dict[tuple[int, int], tuple[float, int]],
    missing_value: str,
) -> float:
    if row_variant.variant_id == column_variant.variant_id:
        return 1.0
    if row_variant.db_id is None or column_variant.db_id is None:
        return 0.0 if missing_value == "zero" else float("nan")
    key = tuple(sorted((row_variant.db_id, column_variant.db_id)))
    if key in pairs:
        return pairs[key][0]
    return 0.0 if missing_value == "zero" else float("nan")


def heatmap_ticks(n_variants: int, maximum: int = 60) -> list[int]:
    """Return readable, evenly spaced tick indices while retaining both ends."""
    if n_variants <= maximum:
        return list(range(n_variants))
    step = max(1, (n_variants - 1 + maximum - 2) // (maximum - 1))
    ticks = list(range(0, n_variants, step))
    if ticks[-1] != n_variants - 1:
        ticks.append(n_variants - 1)
    return ticks


def natural_label_key(label: str) -> tuple[tuple[int, object], ...]:
    """Sort A, B, ..., H2, H10 naturally and case-insensitively."""
    return tuple(
        (0, int(part)) if part.isdigit() else (1, part.casefold())
        for part in re.split(r"([0-9]+)", label)
        if part
    )


def variants_in_haplotype_order(
    variants: Sequence[VariantRecord], haplotypes: Sequence[HaplotypeRecord]
) -> tuple[list[VariantRecord], list[str], list[str]]:
    labels = sorted(
        (record.request.label for record in haplotypes), key=natural_label_key
    )
    label_rank = {label: index for index, label in enumerate(labels)}

    primary_labels: dict[str, str] = {}
    for variant in variants:
        known_labels = [
            label for label in variant.haplotype_labels if label in label_rank
        ]
        if not known_labels:
            raise ValueError(f"Variant {variant.variant_id} has no supplied haplotype label")
        primary_labels[variant.variant_id] = min(
            known_labels, key=lambda label: label_rank[label]
        )

    ordered = sorted(
        variants,
        key=lambda variant: (
            label_rank[primary_labels[variant.variant_id]],
            variant.position,
            variant.ref,
            variant.alt,
        ),
    )
    ordered_primary_labels = [
        primary_labels[variant.variant_id] for variant in ordered
    ]
    return ordered, ordered_primary_labels, labels


def write_ld_heatmap(
    path: Path,
    variants: Sequence[VariantRecord],
    pairs: dict[tuple[int, int], tuple[float, int]],
    missing_value: str,
    dpi: int,
    annotation_labels: Sequence[str] | None = None,
    haplotype_labels: Sequence[str] | None = None,
) -> None:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import BoundaryNorm, ListedColormap
        import numpy as np
    except ImportError as exc:
        raise ValueError(
            "Matplotlib is required to create the LD heatmap"
        ) from exc

    values = np.asarray(
        [
            [matrix_value(row, column, pairs, missing_value) for column in variants]
            for row in variants
        ],
        dtype=float,
    )
    # Display the diagonal and lower half only. Missing values remain masked too.
    triangle_mask = np.triu(np.ones(values.shape, dtype=bool), k=1)
    plotted_values = np.ma.masked_where(triangle_mask | np.isnan(values), values)

    # These are the same five LD bands used by the SNPXplorer haplotype heatmap.
    colors = ["navy", "lightblue", "green", "orange", "red"]
    boundaries = [-0.0000001, 0.2000001, 0.4000001, 0.6000001, 0.8000001, 1.0000001]
    colormap = ListedColormap(colors, name="snpxplorer_ld")
    colormap.set_bad("white")
    normalization = BoundaryNorm(boundaries, colormap.N)

    n_variants = len(variants)
    figure_side = max(8.0, min(20.0, 4.0 + 0.28 * n_variants))
    if annotation_labels is None:
        figure = plt.figure(figsize=(figure_side, figure_side + 1.0))
        grid = figure.add_gridspec(
            2, 1, height_ratios=[1.2, figure_side], hspace=0.08
        )
        legend_grid = grid[0].subgridspec(
            3,
            3,
            height_ratios=[1, 1, 1],
            width_ratios=[1, 2, 1],
            hspace=0,
        )
        legend_axis = figure.add_subplot(legend_grid[1, 1])
        axis = figure.add_subplot(grid[1])
        annotation_axis = None
        row_annotation_axis = None
    else:
        figure = plt.figure(figsize=(figure_side, figure_side + 1.5))
        grid = figure.add_gridspec(
            3,
            2,
            height_ratios=[1.2, figure_side, 0.65],
            width_ratios=[0.65, figure_side],
            hspace=0.04,
            wspace=0.03,
        )
        legend_grid = grid[0, 1].subgridspec(
            3,
            3,
            height_ratios=[1, 1, 1],
            width_ratios=[1, 2, 1],
            hspace=0,
        )
        legend_axis = figure.add_subplot(legend_grid[1, 1])
        row_annotation_axis = figure.add_subplot(grid[1, 0])
        axis = figure.add_subplot(grid[1, 1], sharey=row_annotation_axis)
        annotation_axis = figure.add_subplot(grid[2, 1], sharex=axis)

    image = axis.imshow(
        plotted_values,
        cmap=colormap,
        norm=normalization,
        interpolation="nearest",
        origin="upper",
        aspect="equal",
        rasterized=True,
    )
    if annotation_labels is not None:
        axis.set_anchor("N")

    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)
    if annotation_axis is None:
        axis.set_xlabel("SNPs (physical-position order)")
        axis.set_ylabel("SNPs (physical-position order)")

    if annotation_axis is not None:
        if haplotype_labels is None:
            raise ValueError("Haplotype labels are required for the annotation strip")
        label_rank = {
            label: index for index, label in enumerate(haplotype_labels)
        }
        annotation_values = np.asarray(
            [[label_rank[label] for label in annotation_labels]], dtype=float
        )
        if len(haplotype_labels) <= 20:
            annotation_colors = plt.get_cmap("tab20").colors[
                : len(haplotype_labels)
            ]
        else:
            annotation_colors = [
                plt.get_cmap("hsv")(index / len(haplotype_labels))
                for index in range(len(haplotype_labels))
            ]
        haplotype_colormap = ListedColormap(annotation_colors, name="haplotypes")
        annotation_axis.imshow(
            annotation_values,
            cmap=haplotype_colormap,
            vmin=-0.5,
            vmax=max(0.5, len(haplotype_labels) - 0.5),
            interpolation="nearest",
            aspect="auto",
        )
        row_annotation_axis.imshow(
            annotation_values.T,
            cmap=haplotype_colormap,
            vmin=-0.5,
            vmax=max(0.5, len(haplotype_labels) - 0.5),
            interpolation="nearest",
            aspect="auto",
        )
        annotation_axis.set_yticks([])

        groups: list[tuple[str, int, int]] = []
        group_start = 0
        for index in range(1, len(annotation_labels) + 1):
            if (
                index == len(annotation_labels)
                or annotation_labels[index] != annotation_labels[group_start]
            ):
                groups.append(
                    (annotation_labels[group_start], group_start, index - 1)
                )
                group_start = index
        centers = [(start + end) / 2 for _, start, end in groups]
        annotation_axis.set_xticks(
            centers, labels=[label for label, _, _ in groups]
        )
        annotation_axis.tick_params(
            axis="x",
            top=False,
            labeltop=False,
            bottom=True,
            labelbottom=True,
            labelsize=22,
            length=0,
        )
        for label in annotation_axis.get_xticklabels():
            label.set_fontweight("bold")
        row_annotation_axis.set_yticks(
            centers, labels=[label for label, _, _ in groups]
        )
        row_annotation_axis.tick_params(
            axis="y",
            left=True,
            labelleft=True,
            right=False,
            labelright=False,
            labelsize=22,
            length=0,
        )
        row_annotation_axis.set_xticks([])
        for label in row_annotation_axis.get_yticklabels():
            label.set_fontweight("bold")
        for _, _, end in groups[:-1]:
            boundary = end + 0.5
            annotation_axis.axvline(boundary, color="white", linewidth=1.5)
            row_annotation_axis.axhline(boundary, color="white", linewidth=1.5)
            axis.axvline(boundary, color="white", linewidth=0.8)
            axis.axhline(boundary, color="white", linewidth=0.8)
        annotation_axis.set_xlim(-0.5, n_variants - 0.5)
        row_annotation_axis.set_ylim(n_variants - 0.5, -0.5)
        axis.tick_params(
            axis="both",
            left=False,
            labelleft=False,
            bottom=False,
            labelbottom=False,
        )
        axis.set_xlabel("")
        axis.set_ylabel("")
        for annotation_spine in (
            *annotation_axis.spines.values(),
            *row_annotation_axis.spines.values(),
        ):
            annotation_spine.set_visible(False)

    colorbar = figure.colorbar(
        image,
        cax=legend_axis,
        orientation="horizontal",
        boundaries=boundaries,
        ticks=[0.1, 0.3, 0.5, 0.7, 0.9],
    )
    colorbar.ax.set_xticklabels(
        ["0–0.2", ">0.2–0.4", ">0.4–0.6", ">0.6–0.8", ">0.8–1.0"]
    )
    colorbar.ax.tick_params(axis="x", labelsize=18, length=0, pad=5)
    colorbar.set_label(r"LD $r^2$", fontsize=20, labelpad=8)
    colorbar.ax.xaxis.set_label_position("top")

    if annotation_axis is None:
        figure.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.95)
    else:
        figure.subplots_adjust(left=0.10, right=0.98, bottom=0.08, top=0.95)
        # The equal-aspect heatmap can be slightly smaller than its GridSpec cell.
        # Bind both annotation strips to its final rendered bounds so their ends
        # align exactly, then use the same small gap on the left and bottom.
        figure.canvas.draw()
        heatmap_position = axis.get_position()
        row_position = row_annotation_axis.get_position()
        column_position = annotation_axis.get_position()
        annotation_gap = 0.006
        row_annotation_axis.set_position(
            [
                heatmap_position.x0 - annotation_gap - row_position.width,
                heatmap_position.y0,
                row_position.width,
                heatmap_position.height,
            ]
        )
        annotation_axis.set_position(
            [
                heatmap_position.x0,
                heatmap_position.y0 - annotation_gap - column_position.height,
                heatmap_position.width,
                column_position.height,
            ]
        )
    figure.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def write_position_constrained_ld_heatmap(
    path: Path,
    variants: Sequence[VariantRecord],
    pairs: dict[tuple[int, int], tuple[float, int]],
    dpi: int,
) -> None:
    """Draw observed LD as a roof above a true genomic-position axis."""
    try:
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
        from matplotlib.colors import BoundaryNorm, ListedColormap
        from matplotlib.ticker import FuncFormatter, MaxNLocator
        import numpy as np
    except ImportError as exc:
        raise ValueError(
            "Matplotlib and NumPy are required to create the genomic LD heatmap"
        ) from exc

    if len(variants) < 2:
        raise ValueError("At least two variants are required for the genomic LD heatmap")

    positions = np.asarray([variant.position for variant in variants], dtype=float)
    genomic_start = float(positions.min())
    genomic_end = float(positions.max())
    genomic_span = max(1.0, genomic_end - genomic_start)
    positive_gaps = np.diff(np.unique(positions))
    positive_gaps = positive_gaps[positive_gaps > 0]
    if len(positive_gaps) == 0:
        raise ValueError("The genomic LD heatmap requires variants at distinct positions")

    # Keep cells visible at overview scale without stretching them across large
    # genomic gaps. Very dense variants may overlap slightly at this resolution.
    minimum_display_half_width = genomic_span / 1_200.0
    maximum_half_width = max(
        minimum_display_half_width,
        float(np.median(positive_gaps)) * 0.45,
    )
    half_widths: list[float] = []
    for index, position in enumerate(positions):
        neighbor_half_gaps = [maximum_half_width]
        if index > 0 and position > positions[index - 1]:
            neighbor_half_gaps.append((position - positions[index - 1]) * 0.45)
        if index + 1 < len(positions) and positions[index + 1] > position:
            neighbor_half_gaps.append((positions[index + 1] - position) * 0.45)
        half_widths.append(
            max(minimum_display_half_width, min(neighbor_half_gaps))
        )

    polygons: list[list[tuple[float, float]]] = []
    ld_values: list[float] = []
    for row_index in range(1, len(variants)):
        row_variant = variants[row_index]
        if row_variant.db_id is None:
            continue
        row_low = positions[row_index] - half_widths[row_index]
        row_high = positions[row_index] + half_widths[row_index]
        for column_index in range(row_index):
            column_variant = variants[column_index]
            if column_variant.db_id is None:
                continue
            key = tuple(sorted((row_variant.db_id, column_variant.db_id)))
            if key not in pairs:
                continue
            column_low = positions[column_index] - half_widths[column_index]
            column_high = positions[column_index] + half_widths[column_index]

            # Transform a matrix cell into a rotated polygon:
            # genomic X = (column coordinate + row coordinate) / 2
            # height    = (row coordinate - column coordinate) / 2
            rectangle = [
                (column_low, row_low),
                (column_high, row_low),
                (column_high, row_high),
                (column_low, row_high),
            ]
            polygons.append(
                [
                    ((column + row) / 2.0, (row - column) / 2.0)
                    for column, row in rectangle
                ]
            )
            ld_values.append(pairs[key][0])

    colors = ["navy", "lightblue", "green", "orange", "red"]
    boundaries = [-0.0000001, 0.2000001, 0.4000001, 0.6000001, 0.8000001, 1.0000001]
    colormap = ListedColormap(colors, name="snpxplorer_ld")
    normalization = BoundaryNorm(boundaries, colormap.N)

    figure = plt.figure(figsize=(16, 9))
    grid = figure.add_gridspec(
        2, 1, height_ratios=[1.3, 8.0], hspace=0.10
    )
    legend_grid = grid[0].subgridspec(
        3,
        3,
        height_ratios=[1, 1, 1],
        width_ratios=[1, 2, 1],
        hspace=0,
    )
    legend_axis = figure.add_subplot(legend_grid[1, 1])
    axis = figure.add_subplot(grid[1])

    # Match the conventional matrix view: pairs absent from the SQLite LD
    # table share the lowest (navy) display category with observed r² values
    # from 0 to 0.2. Draw the complete position-constrained roof first, then
    # overlay every observed LD cell below.
    roof_midpoint = (genomic_start + genomic_end) / 2.0
    axis.fill(
        [genomic_start, roof_midpoint, genomic_end],
        [0.0, genomic_span / 2.0, 0.0],
        color="navy",
        edgecolor="none",
        zorder=0,
    )

    collection = PolyCollection(
        polygons,
        array=np.asarray(ld_values),
        cmap=colormap,
        norm=normalization,
        edgecolors="none",
        antialiaseds=False,
        rasterized=True,
        zorder=2,
    )
    axis.add_collection(collection)

    axis.set_xlim(genomic_start - genomic_span * 0.01, genomic_end + genomic_span * 0.01)
    axis.set_ylim(0, genomic_span * 0.52)
    axis.set_aspect("equal", adjustable="box", anchor="S")
    axis.scatter(
        positions,
        np.zeros_like(positions),
        marker="|",
        s=18,
        linewidths=0.5,
        color="black",
        clip_on=False,
        zorder=3,
    )
    axis.xaxis.set_major_locator(MaxNLocator(nbins=8))
    axis.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1e6:.3f}"))
    axis.tick_params(axis="x", labelsize=14)
    axis.set_xlabel(
        f"Genomic position on chromosome {variants[0].chrom} (Mb)", fontsize=17
    )
    axis.set_yticks([])
    for name, spine in axis.spines.items():
        spine.set_visible(name == "bottom")

    colorbar = figure.colorbar(
        collection,
        cax=legend_axis,
        orientation="horizontal",
        boundaries=boundaries,
        ticks=[0.1, 0.3, 0.5, 0.7, 0.9],
    )
    colorbar.ax.set_xticklabels(
        ["NA and 0–0.2", ">0.2–0.4", ">0.4–0.6", ">0.6–0.8", ">0.8–1.0"]
    )
    colorbar.ax.tick_params(axis="x", labelsize=15, length=0, pad=5)
    colorbar.set_label(r"LD $r^2$", fontsize=18, labelpad=7)
    colorbar.ax.xaxis.set_label_position("top")

    figure.subplots_adjust(left=0.06, right=0.99, bottom=0.12, top=0.96)
    figure.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def write_outputs(
    prefix: Path,
    haplotypes: Sequence[HaplotypeRecord],
    variants: Sequence[VariantRecord],
    pairs: dict[tuple[int, int], tuple[float, int]],
    missing_value: str,
    heatmap_dpi: int,
    haplotype_table: Path,
    ld_database: Path,
) -> dict[str, Path]:
    paths = output_paths(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    write_square_matrix(paths["matrix"], variants, pairs, missing_value, False)
    write_square_matrix(paths["mask"], variants, pairs, missing_value, True)
    write_ld_heatmap(
        paths["heatmap"],
        variants,
        pairs,
        missing_value,
        heatmap_dpi,
    )
    haplotype_ordered_variants, primary_labels, ordered_haplotype_labels = (
        variants_in_haplotype_order(variants, haplotypes)
    )
    write_ld_heatmap(
        paths["haplotype_heatmap"],
        haplotype_ordered_variants,
        pairs,
        missing_value,
        heatmap_dpi,
        primary_labels,
        ordered_haplotype_labels,
    )
    write_position_constrained_ld_heatmap(
        paths["position_heatmap"], variants, pairs, heatmap_dpi
    )
    haplotype_plot_index = {
        variant.variant_id: index
        for index, variant in enumerate(haplotype_ordered_variants, start=1)
    }
    primary_plot_haplotype = {
        variant.variant_id: label
        for variant, label in zip(haplotype_ordered_variants, primary_labels)
    }

    db_id_to_variant = {
        variant.db_id: variant for variant in variants if variant.db_id is not None
    }
    with gzip.open(paths["pairs"], "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["variant_1", "variant_2", "r2", "distance_bp"])
        sorted_pairs = sorted(
            pairs.items(),
            key=lambda item: (
                db_id_to_variant[item[0][0]].position,
                db_id_to_variant[item[0][1]].position,
            ),
        )
        for (v1, v2), (r2, distance) in sorted_pairs:
            writer.writerow(
                [
                    db_id_to_variant[v1].variant_id,
                    db_id_to_variant[v2].variant_id,
                    f"{r2:.3f}",
                    distance,
                ]
            )

    with paths["variants"].open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "matrix_index",
                "variant_id",
                "chromosome",
                "position",
                "ref",
                "alt",
                "haplotype_labels",
                "haplotype_ids",
                "in_multiple_haplotypes",
                "haplotype_plot_index",
                "primary_plot_haplotype",
                "ld_db_id",
                "ld_db_uniq",
                "mapping_status",
            ]
        )
        for index, variant in enumerate(variants, start=1):
            writer.writerow(
                [
                    index,
                    variant.variant_id,
                    variant.chrom,
                    variant.position,
                    variant.ref,
                    variant.alt,
                    "|".join(sorted(variant.haplotype_labels)),
                    "|".join(sorted(variant.haplotype_ids)),
                    len(variant.haplotype_ids) > 1,
                    haplotype_plot_index[variant.variant_id],
                    primary_plot_haplotype[variant.variant_id],
                    "" if variant.db_id is None else variant.db_id,
                    variant.db_uniq or "",
                    variant.mapping_status,
                ]
            )

    with paths["haplotypes"].open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "label",
                "haplotype_id",
                "chromosome",
                "start",
                "end",
                "span_kb",
                "declared_n_snps",
                "retrieved_n_variant_ids",
                "unique_n_variant_ids",
            ]
        )
        for record in haplotypes:
            writer.writerow(
                [
                    record.request.label,
                    record.request.haplotype_id,
                    record.request.chrom,
                    record.request.start,
                    record.request.end,
                    record.kb,
                    record.declared_n_snps,
                    len(record.variants),
                    len(set(record.variants)),
                ]
            )

    mapped_count = sum(variant.db_id is not None for variant in variants)
    possible_pairs = len(variants) * (len(variants) - 1) // 2
    possible_mapped_pairs = mapped_count * (mapped_count - 1) // 2
    qc = {
        "haplotype_source": str(haplotype_table.resolve()),
        "ld_source": str(ld_database.resolve()),
        "ld_source_note": (
            "All LD edges stored for the selected variants were retrieved from SQLite. "
            "No additional r2 threshold was applied by this script. An absent edge is "
            "not distinguishable from an edge omitted when the database was built."
        ),
        "missing_edge_matrix_value": "0" if missing_value == "zero" else "NA",
        "heatmap_dpi": heatmap_dpi,
        "heatmap_shape": "lower_triangle_including_diagonal",
        "position_heatmap_shape": "roof_with_genomic_position_axis",
        "position_heatmap_missing_pairs": "navy (combined with observed r2 0-0.2)",
        "haplotype_heatmap_order": ordered_haplotype_labels,
        "multi_haplotype_variant_plot_assignment": (
            "Earliest natural-sorted haplotype label; each variant is plotted once"
        ),
        "heatmap_ld_bands": {
            "0.0-0.2": "navy",
            ">0.2-0.4": "lightblue",
            ">0.4-0.6": "green",
            ">0.6-0.8": "orange",
            ">0.8-1.0": "red",
        },
        "n_haplotypes": len(haplotypes),
        "n_unique_variants": len(variants),
        "n_variants_in_multiple_haplotypes": sum(
            len(variant.haplotype_ids) > 1 for variant in variants
        ),
        "n_variants_mapped_to_ld_database": mapped_count,
        "n_variants_not_mapped_to_ld_database": len(variants) - mapped_count,
        "n_ambiguous_variant_mappings": sum(
            variant.mapping_status == "ambiguous_first_id" for variant in variants
        ),
        "n_all_possible_variant_pairs": possible_pairs,
        "n_possible_pairs_among_mapped_variants": possible_mapped_pairs,
        "n_pairs_observed_in_ld_database": len(pairs),
        "n_pairs_absent_among_mapped_variants": possible_mapped_pairs - len(pairs),
        "outputs": {name: str(path.resolve()) for name, path in paths.items()},
    }
    with paths["qc"].open("wt", encoding="utf-8") as handle:
        json.dump(qc, handle, indent=2)
        handle.write("\n")
    return paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build a physical-position-ordered LD matrix for all variants in one or "
            "more SNPXplorer haplotypes, using the SNPXplorer SQLite LD database."
        )
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "--haplotype",
        action="append",
        metavar="LABEL=CHR:START-END",
        help="Labeled haplotype; repeat this option for every haplotype",
    )
    inputs.add_argument(
        "--manifest",
        type=Path,
        help="Tab-separated file with label and haplotype_id columns",
    )
    parser.add_argument(
        "--haplotype-table",
        type=Path,
        default=DEFAULT_HAPLOTYPE_TABLE,
        help=f"Haplotype membership table (default: {DEFAULT_HAPLOTYPE_TABLE})",
    )
    parser.add_argument(
        "--ld-directory",
        type=Path,
        default=DEFAULT_LD_DIRECTORY,
        help=f"Directory containing ld_chr*.sqlite (default: {DEFAULT_LD_DIRECTORY})",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("haplotype_region_ld"),
        help="Output path prefix (default: ./haplotype_region_ld)",
    )
    parser.add_argument(
        "--missing-value",
        choices=("zero", "nan"),
        default="zero",
        help=(
            "Value used in the matrix when an LD edge is absent from SQLite: zero "
            "or NA (default: zero). The observed-mask output always identifies these."
        ),
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Enable more detailed logging"
    )
    parser.add_argument(
        "--heatmap-dpi",
        type=int,
        default=300,
        help="Resolution of the PNG heatmap (default: 300)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    try:
        requests = load_requests(args.haplotype, args.manifest)
        if args.heatmap_dpi < 72:
            raise ValueError("--heatmap-dpi must be at least 72")
        if not args.haplotype_table.is_file():
            raise FileNotFoundError(f"Haplotype table not found: {args.haplotype_table}")
        LOGGER.info("Loading %d haplotype records", len(requests))
        haplotypes = [
            load_haplotype_record(args.haplotype_table, request)
            for request in requests
        ]
        variants = combine_variants(haplotypes)
        chromosome = requests[0].chrom
        ld_database = args.ld_directory / f"ld_chr{chromosome}.sqlite"
        if not ld_database.is_file():
            raise FileNotFoundError(f"LD database not found: {ld_database}")

        LOGGER.info(
            "Combined haplotypes contain %d unique variants on chromosome %s",
            len(variants),
            chromosome,
        )
        with open_read_only_database(ld_database) as connection:
            map_variants_to_database(connection, variants)
            mapped = sum(variant.db_id is not None for variant in variants)
            LOGGER.info("Mapped %d/%d variants to SQLite", mapped, len(variants))
            pairs = query_stored_ld(connection, variants)
            LOGGER.info("Retrieved %d stored within-set LD pairs", len(pairs))

        paths = write_outputs(
            args.output_prefix,
            haplotypes,
            variants,
            pairs,
            args.missing_value,
            args.heatmap_dpi,
            args.haplotype_table,
            ld_database,
        )
        LOGGER.info("Finished. QC summary: %s", paths["qc"])
        return 0
    except (OSError, sqlite3.Error, ValueError) as exc:
        LOGGER.error("%s", exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
