#!/usr/bin/env python3
"""P1/P2/P5 virtual-Visium-55 reconstruction from measured 16 um Visium HD.

This is a geometry-only degradation followed by a fully independent fit of the
reference-free cell_resolved_v3 mathematical engine.  It never reuses fitted
16 um profiles/programs/cell expression, Xenium, AI annotations, scRNA-seq, or
the upstream fixed 44 HE programs.
"""

from __future__ import annotations

import argparse
import gc
import json
import math
import os
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

import cell_resolved_v3 as core


SOURCE_BIN_UM = 16.0
SPOT_DIAMETER_UM = 55.0
SPOT_RADIUS_UM = SPOT_DIAMETER_UM / 2.0
SPOT_PITCH_UM = 100.0
DEFAULT_COVERAGE = 0.95


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", choices=("P1", "P2", "P5"), required=True)
    parser.add_argument(
        "--data-root", type=Path, default=Path('DATA/Visium_HD')
    )
    parser.add_argument("--positions-tsv", type=Path, required=True)
    parser.add_argument(
        "--eligible-cell-h5ad", type=Path,
        help=(
            "Optional full measured-domain CellViT result. When provided, every "
            "listed real nucleus is assigned to virtual spots by within-component "
            "Voronoi contour overlap, including nuclei between capture circles."
        ),
    )
    parser.add_argument(
        "--marker-catalog", type=Path,
        default=Path(
            'DATA/marker_catalog.json'
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allocation-output", type=Path, required=True)
    parser.add_argument("--qc-json", type=Path, required=True)
    parser.add_argument("--spot-positions-output", type=Path, required=True)
    parser.add_argument("--coverage-threshold", type=float, default=DEFAULT_COVERAGE)
    parser.add_argument("--embedding-pcs", type=int, default=16)
    parser.add_argument("--pca-sample-per-class", type=int, default=5000)
    parser.add_argument("--program-genes", type=int, default=3000)
    parser.add_argument("--residual-genes", type=int, default=768)
    parser.add_argument("--residual-pcs", type=int, default=16)
    parser.add_argument("--graph-neighbors", type=int, default=24)
    parser.add_argument("--graph-radius-spots", type=float, default=6.0)
    parser.add_argument("--class-prior-strength", type=float, default=0.08)
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--permutations", type=int, default=100)
    parser.add_argument("--bootstrap-replicates", type=int, default=20)
    parser.add_argument("--cell-batch-size", type=int, default=1024)
    parser.add_argument("--calibration-cell-h5ad", type=Path)
    parser.add_argument("--calibration-max-cells", type=int, default=12000)
    parser.add_argument("--subspot-neighbors", type=int, default=9)
    parser.add_argument("--subspot-radius-spots", type=float, default=2.0)
    parser.add_argument("--subspot-spatial-scale", type=float, default=0.5)
    parser.add_argument("--morphology-offset-scale", type=float, default=0.5)
    parser.add_argument("--within-type-temperature", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=20260826)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--skip-ablations", action="store_true")
    args = parser.parse_args()
    if not 0 < args.coverage_threshold <= 1:
        parser.error("--coverage-threshold must be in (0, 1]")
    # Core spatial coordinates remain in 16 um array units.  Six virtual spots
    # therefore correspond to this radius in those coordinates.
    args.graph_radius_bins = (
        args.graph_radius_spots * SPOT_PITCH_UM / SOURCE_BIN_UM
    )
    args.subspot_radius_bins = (
        args.subspot_radius_spots * SPOT_PITCH_UM / SOURCE_BIN_UM
    )
    return args


def source_paths(args: argparse.Namespace) -> dict[str, Path]:
    stem = f"Visium_HD_Human_Colon_Cancer_{args.sample}_tissue_image"
    root = args.data_root / args.sample / "binned_outputs" / "square_016um"
    return {
        "matrix": root / "filtered_feature_bc_matrix.h5",
        "positions_parquet": root / "spatial" / "tissue_positions.parquet",
        "cell_features": (
            args.data_root / "cellvit_pp_outputs" / args.sample
            / f"{stem}_cell_features.h5"
        ),
    }


def quadrant_disk_rectangle_area(a: float, b: float, radius: float) -> float:
    """Area of [0,a] x [0,b] inside a radius-r disk centered at the origin."""
    a = min(max(float(a), 0.0), radius)
    b = min(max(float(b), 0.0), radius)
    if a <= 0 or b <= 0:
        return 0.0
    r2 = radius * radius
    if a * a + b * b <= r2 * (1.0 + 1e-14):
        return a * b
    xa = math.sqrt(max(r2 - a * a, 0.0))
    xb = math.sqrt(max(r2 - b * b, 0.0))
    value = 0.5 * (
        a * xa + b * xb
        + r2 * (
            math.asin(min(a / radius, 1.0))
            + math.asin(min(b / radius, 1.0))
            - math.pi / 2.0
        )
    )
    return max(value, 0.0)


def oriented_disk_primitive(x: float, y: float, radius: float) -> float:
    sx = -1.0 if x < 0 else 1.0
    sy = -1.0 if y < 0 else 1.0
    return sx * sy * quadrant_disk_rectangle_area(abs(x), abs(y), radius)


def circle_rectangle_area(
    center_x: float, center_y: float, radius: float, col: int, row: int
) -> float:
    """Exact area of a circle intersected with a unit square centered at col,row."""
    x0, x1 = col - 0.5 - center_x, col + 0.5 - center_x
    y0, y1 = row - 0.5 - center_y, row + 0.5 - center_y
    area = (
        oriented_disk_primitive(x1, y1, radius)
        - oriented_disk_primitive(x0, y1, radius)
        - oriented_disk_primitive(x1, y0, radius)
        + oriented_disk_primitive(x0, y0, radius)
    )
    return float(np.clip(area, 0.0, 1.0))


def circle_polygon_intersection_area(
    polygon: np.ndarray, center: np.ndarray, radius: float
) -> float:
    """Exact disk/simple-polygon intersection by edge triangle decomposition."""
    if len(polygon) < 3:
        return 0.0
    points = np.asarray(polygon, np.float64) - np.asarray(center, np.float64)
    r2 = radius * radius
    total = 0.0
    for a, b in zip(points, np.roll(points, -1, axis=0)):
        direction = b - a
        qa = float(np.dot(direction, direction))
        cuts = [0.0, 1.0]
        if qa > 1e-20:
            qb = 2.0 * float(np.dot(a, direction))
            qc = float(np.dot(a, a)) - r2
            discriminant = qb * qb - 4.0 * qa * qc
            if discriminant > 0:
                root = math.sqrt(discriminant)
                for value in ((-qb - root) / (2.0 * qa), (-qb + root) / (2.0 * qa)):
                    if 1e-12 < value < 1.0 - 1e-12:
                        cuts.append(value)
        cuts.sort()
        for left, right in zip(cuts[:-1], cuts[1:]):
            p = a + left * direction
            q = a + right * direction
            midpoint = 0.5 * (p + q)
            cross = float(p[0] * q[1] - p[1] * q[0])
            if float(np.dot(midpoint, midpoint)) <= r2 * (1.0 + 1e-12):
                total += 0.5 * cross
            else:
                dot = float(np.dot(p, q))
                total += 0.5 * r2 * math.atan2(cross, dot)
    return abs(total)


def phase_offsets(pitch_units: float) -> list[tuple[str, float, float]]:
    step = pitch_units / 3.0
    result = [("main", 0.0, 0.0)]
    for k in range(6):
        angle = k * math.pi / 3.0
        result.append(
            (
                f"shift_{k}",
                step * math.cos(angle),
                step * math.sin(angle),
            )
        )
    return result


def build_phase_overlap(
    positions: pd.DataFrame,
    phase_x: float,
    phase_y: float,
    coverage_threshold: float,
) -> tuple[sparse.csr_matrix, pd.DataFrame, dict]:
    """Build accepted virtual spots and exact spot-by-positive-16um-bin overlap."""
    col = positions["array_col"].to_numpy(np.int32)
    row = positions["array_row"].to_numpy(np.int32)
    lookup = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(row, col))}
    radius = SPOT_RADIUS_UM / SOURCE_BIN_UM
    pitch = SPOT_PITCH_UM / SOURCE_BIN_UM
    row_pitch = math.sqrt(3.0) * pitch / 2.0
    xmin, xmax = float(col.min()), float(col.max())
    ymin, ymax = float(row.min()), float(row.max())
    m0 = math.floor((ymin - radius - phase_y) / row_pitch) - 1
    m1 = math.ceil((ymax + radius - phase_y) / row_pitch) + 1
    circle_area = math.pi * radius * radius
    spot_meta: list[dict] = []
    overlap_rows: list[int] = []
    overlap_cols: list[int] = []
    overlap_values: list[float] = []
    candidate_count = 0
    for m in range(m0, m1 + 1):
        center_y = phase_y + row_pitch * m
        parity = m % 2
        n0 = math.floor(
            (xmin - radius - phase_x) / pitch - 0.5 * parity
        ) - 1
        n1 = math.ceil(
            (xmax + radius - phase_x) / pitch - 0.5 * parity
        ) + 1
        for n in range(n0, n1 + 1):
            candidate_count += 1
            center_x = phase_x + pitch * (n + 0.5 * parity)
            local: list[tuple[int, float]] = []
            r0 = math.ceil(center_y - radius - 0.5)
            r1 = math.floor(center_y + radius + 0.5)
            c0 = math.ceil(center_x - radius - 0.5)
            c1 = math.floor(center_x + radius + 0.5)
            for rr in range(r0, r1 + 1):
                for cc in range(c0, c1 + 1):
                    source_bin = lookup.get((rr, cc))
                    if source_bin is None:
                        continue
                    area = circle_rectangle_area(center_x, center_y, radius, cc, rr)
                    if area > 1e-12:
                        local.append((source_bin, area))
            covered = sum(value for _, value in local)
            coverage = covered / circle_area
            if coverage + 1e-12 < coverage_threshold or not local:
                continue
            spot_index = len(spot_meta)
            spot_meta.append(
                {
                    "hex_row": int(m),
                    "hex_col": int(n),
                    "array_row": center_y,
                    "array_col": center_x,
                    "coverage_fraction": coverage,
                }
            )
            for source_bin, value in local:
                overlap_rows.append(spot_index)
                overlap_cols.append(source_bin)
                overlap_values.append(value)
    overlap = sparse.coo_matrix(
        (
            np.asarray(overlap_values, np.float64),
            (
                np.asarray(overlap_rows, np.int32),
                np.asarray(overlap_cols, np.int32),
            ),
        ),
        shape=(len(spot_meta), len(positions)),
    ).tocsr()
    overlap.sum_duplicates()
    bin_sum = np.asarray(overlap.sum(axis=0)).ravel()
    if len(bin_sum) and float(bin_sum.max()) > 1.0 + 1e-10:
        raise RuntimeError(
            f"Virtual spot circles overlap in source-bin mass; max={bin_sum.max():.12g}"
        )
    meta = pd.DataFrame(spot_meta)
    coverages = meta["coverage_fraction"].to_numpy() if len(meta) else np.asarray([])
    audit = {
        "candidate_spots": candidate_count,
        "accepted_spots": int(overlap.shape[0]),
        "coverage_threshold": float(coverage_threshold),
        "coverage_min": float(coverages.min()) if len(coverages) else None,
        "coverage_median": float(np.median(coverages)) if len(coverages) else None,
        "coverage_max": float(coverages.max()) if len(coverages) else None,
        "overlap_links": int(overlap.nnz),
        "max_source_bin_captured_fraction": (
            float(bin_sum.max()) if len(bin_sum) else 0.0
        ),
    }
    return overlap, meta, audit


def build_virtual_spots(
    matrix: sparse.csc_matrix,
    positions: pd.DataFrame,
    affine: dict,
    sample: str,
    coverage_threshold: float,
) -> tuple[sparse.csc_matrix, sparse.csr_matrix, pd.DataFrame, np.ndarray, dict]:
    pitch_units = SPOT_PITCH_UM / SOURCE_BIN_UM
    audits: dict[str, dict] = {}
    main_overlap = None
    main_meta = None
    for name, phase_x, phase_y in phase_offsets(pitch_units):
        overlap, meta, audit = build_phase_overlap(
            positions, phase_x, phase_y, coverage_threshold
        )
        audit["phase_x_16um_units"] = phase_x
        audit["phase_y_16um_units"] = phase_y
        audits[name] = audit
        if name == "main":
            main_overlap, main_meta = overlap, meta
    assert main_overlap is not None and main_meta is not None
    if main_overlap.shape[0] == 0:
        raise RuntimeError("No virtual Visium spots passed the coverage gate")
    core.log("Aggregating positive 16 um bins into accepted virtual 55 um spots")
    virtual = (matrix @ main_overlap.T.tocsc()).tocsc()
    virtual.sum_duplicates()
    virtual.eliminate_zeros()
    library = np.asarray(virtual.sum(axis=0)).ravel()
    keep = library > 0
    if not np.all(keep):
        main_overlap = main_overlap[keep].tocsr()
        main_meta = main_meta.loc[keep].reset_index(drop=True)
        virtual = virtual[:, keep].tocsc()
        library = library[keep]
    barcodes = np.asarray(
        [
            f"{sample}_VIS55_m{int(m):+06d}_n{int(n):+06d}"
            for m, n in zip(main_meta["hex_row"], main_meta["hex_col"])
        ],
        dtype="U",
    )
    transform = affine["transform"]
    intercept = affine["intercept"]
    uv = main_meta[["array_col", "array_row"]].to_numpy(np.float64)
    xy = intercept[None, :] + uv @ transform.T
    main_meta = main_meta.copy()
    main_meta.insert(0, "barcode", barcodes)
    main_meta.insert(1, "in_tissue", 1)
    main_meta["pxl_col_in_fullres"] = xy[:, 0]
    main_meta["pxl_row_in_fullres"] = xy[:, 1]
    captured = np.asarray(main_overlap.sum(axis=0)).ravel()
    gap_weight = np.maximum(1.0 - captured, 0.0)
    gap_gene = np.asarray(matrix @ gap_weight).ravel().astype(np.float64)
    original_gene = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64)
    spot_gene = np.asarray(virtual.sum(axis=1)).ravel().astype(np.float64)
    closure_error = float(np.max(np.abs(original_gene - spot_gene - gap_gene)))
    geometry_audit = {
        "schema": "virtual_visium55_geometry_v1",
        "source_bin_um": SOURCE_BIN_UM,
        "spot_diameter_um": SPOT_DIAMETER_UM,
        "spot_radius_um": SPOT_RADIUS_UM,
        "spot_center_pitch_um": SPOT_PITCH_UM,
        "lattice": "hexagonal_odd_row_offset",
        "main_phase": "array_origin_(0,0)",
        "phase_sensitivity": audits,
        "main_spots_after_positive_mass_gate": int(virtual.shape[1]),
        "main_spot_matrix_nnz": int(virtual.nnz),
        "source_total_mass": float(matrix.sum()),
        "virtual_spot_total_mass": float(virtual.sum()),
        "gap_total_mass": float(gap_gene.sum()),
        "max_gene_mass_closure_error": closure_error,
    }
    return virtual, main_overlap, main_meta, gap_gene, geometry_audit


def virtual_cell_ownership(
    centroid_uv: np.ndarray, centers_uv: np.ndarray, radius: float
) -> tuple[np.ndarray, np.ndarray]:
    tree = cKDTree(centers_uv)
    distance, owner = tree.query(centroid_uv, k=1, distance_upper_bound=radius + 1e-10)
    eligible = np.isfinite(distance) & (owner < len(centers_uv))
    owner_all = np.full(len(centroid_uv), -1, np.int32)
    owner_all[eligible] = owner[eligible].astype(np.int32)
    return owner_all, eligible


def spot_source_components(
    overlap: sparse.csr_matrix, source_component: np.ndarray
) -> np.ndarray:
    """Assign each virtual spot to its dominant measured 16 um component."""
    result = np.full(overlap.shape[0], -1, np.int32)
    for spot in range(overlap.shape[0]):
        start, end = int(overlap.indptr[spot]), int(overlap.indptr[spot + 1])
        bins = overlap.indices[start:end]
        weights = overlap.data[start:end].astype(np.float64)
        if not len(bins):
            continue
        components = source_component[bins]
        mass = np.bincount(
            components, weights=weights,
            minlength=int(source_component.max(initial=-1)) + 1,
        )
        result[spot] = int(np.argmax(mass))
    if np.any(result < 0):
        raise RuntimeError("A virtual spot has no measured-domain component")
    return result


def full_real_nucleus_domain(
    selected_h5ad: Path,
    all_cell_id: np.ndarray,
    source_component: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Map the frozen measured-domain real nuclei back to the CellViT table."""
    selected = ad.read_h5ad(selected_h5ad, backed="r")
    try:
        if "cell_id" not in selected.obs or "owner_bin_index" not in selected.obs:
            raise RuntimeError(
                "--eligible-cell-h5ad must contain cell_id and owner_bin_index"
            )
        selected_id = selected.obs["cell_id"].to_numpy(np.int64)
        selected_owner = selected.obs["owner_bin_index"].to_numpy(np.int64)
    finally:
        selected.file.close()
    if len(np.unique(selected_id)) != len(selected_id):
        raise RuntimeError("Eligible-cell H5AD contains duplicate cell_id")
    if np.any((selected_owner < 0) | (selected_owner >= len(source_component))):
        raise RuntimeError("Eligible-cell H5AD has invalid measured-bin ownership")
    lookup = pd.Series(selected_owner, index=selected_id)
    owner_source_bin = lookup.reindex(all_cell_id).fillna(-1).to_numpy(np.int64)
    eligible = owner_source_bin >= 0
    component_all = np.full(len(all_cell_id), -1, np.int32)
    component_all[eligible] = source_component[owner_source_bin[eligible]]
    if int(eligible.sum()) != len(selected_id):
        raise RuntimeError(
            "Eligible-cell H5AD cell IDs do not map one-to-one to CellViT features"
        )
    return eligible, component_all


def nearest_spot_within_component(
    centroid_uv: np.ndarray,
    eligible: np.ndarray,
    cell_component: np.ndarray,
    centers_uv: np.ndarray,
    spot_component: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Nearest virtual-spot owner without crossing measured tissue components."""
    owner = np.full(len(centroid_uv), -1, np.int32)
    cross_component = np.zeros(len(centroid_uv), dtype=bool)
    for component in np.unique(cell_component[eligible]):
        cells = np.flatnonzero(eligible & (cell_component == component))
        spots = np.flatnonzero(spot_component == component)
        if not len(spots):
            distance, nearest = cKDTree(centers_uv).query(centroid_uv[cells], k=1)
            del distance
            owner[cells] = nearest.astype(np.int32)
            cross_component[cells] = True
            continue
        _, nearest = cKDTree(centers_uv[spots]).query(centroid_uv[cells], k=1)
        owner[cells] = spots[np.asarray(nearest, np.int64)].astype(np.int32)
    if np.any(owner[eligible] < 0):
        raise RuntimeError("A measured-domain real nucleus has no virtual-spot owner")
    return owner, cross_component


def clip_general_halfplane(
    polygon: np.ndarray, normal: np.ndarray, bound: float
) -> np.ndarray:
    """Clip a polygon to dot(normal, point) <= bound."""
    if len(polygon) == 0:
        return polygon
    result: list[np.ndarray] = []
    previous = polygon[-1]
    previous_value = float(np.dot(normal, previous) - bound)
    previous_inside = previous_value <= 1e-10
    for current in polygon:
        current_value = float(np.dot(normal, current) - bound)
        current_inside = current_value <= 1e-10
        if current_inside != previous_inside:
            denominator = previous_value - current_value
            if abs(denominator) > 1e-15:
                fraction = previous_value / denominator
                result.append(previous + fraction * (current - previous))
        if current_inside:
            result.append(current)
        previous = current
        previous_value = current_value
        previous_inside = current_inside
    return np.asarray(result, np.float64).reshape((-1, 2))


def local_voronoi_overlap(
    polygon: np.ndarray, candidate_centers: np.ndarray
) -> list[tuple[int, float]]:
    """Exact polygon fractions in the local nearest-centre Voronoi partition."""
    total_area = core.polygon_area(polygon)
    if total_area <= 1e-8:
        return []
    entries: list[tuple[int, float]] = []
    for local_site, center in enumerate(candidate_centers):
        clipped = polygon
        for other_site, other in enumerate(candidate_centers):
            if local_site == other_site:
                continue
            normal = other - center
            bound = 0.5 * float(np.dot(other, other) - np.dot(center, center))
            clipped = clip_general_halfplane(clipped, normal, bound)
            if len(clipped) < 3:
                break
        area = core.polygon_area(clipped)
        if area > max(1e-8, total_area * 1e-7):
            entries.append((local_site, area / total_area))
    total = sum(value for _, value in entries)
    if total <= 0:
        return []
    return [(site, value / total) for site, value in entries]


def build_full_domain_voronoi_overlap(
    nuclei: dict,
    affine: dict,
    centers_uv: np.ndarray,
    spot_component: np.ndarray,
    cross_component_all: np.ndarray,
    radius: float,
) -> tuple[sparse.csr_matrix, dict]:
    """Allocate every real measured-domain nucleus over local spot Voronoi cells."""
    inverse = affine["inverse"]
    intercept = affine["intercept"]
    offsets = nuclei.pop("all_contour_offsets")
    contours = nuclei.pop("all_contour_xy")
    source = nuclei["source_index"]
    owner = nuclei["owner_bin"]
    component = nuclei["domain_component"]
    trees = {
        int(value): (
            np.flatnonzero(spot_component == value),
            cKDTree(centers_uv[spot_component == value]),
        )
        for value in np.unique(spot_component)
    }
    all_tree = cKDTree(centers_uv)
    link_spot: list[int] = []
    link_cell: list[int] = []
    link_fraction: list[float] = []
    fallback = np.zeros(len(source), dtype=bool)
    direct_fraction = np.zeros(len(source), np.float32)
    evidence = np.ones(len(source), np.int8)
    evidence[cross_component_all[source]] = 2
    last = time.time()
    for local_i, source_i in enumerate(source):
        start, end = int(offsets[source_i]), int(offsets[source_i + 1])
        polygon_xy = contours[start:end]
        valid = len(polygon_xy) >= 3 and np.isfinite(polygon_xy).all()
        polygon = np.empty((0, 2), np.float64)
        if valid:
            polygon = (polygon_xy - intercept) @ inverse.T
            valid = core.polygon_area(polygon) > 1e-8

        if int(component[local_i]) in trees:
            spots, tree = trees[int(component[local_i])]
        else:
            spots = np.arange(len(centers_uv), dtype=np.int64)
            tree = all_tree
        k = min(12, len(spots))
        _, nearest_local = tree.query(nuclei["centroid_uv"][local_i], k=k)
        nearest_local = np.atleast_1d(nearest_local).astype(np.int64)
        candidates = spots[nearest_local]
        entries: list[tuple[int, float]] = []
        if valid:
            local_entries = local_voronoi_overlap(polygon, centers_uv[candidates])
            entries = [(int(candidates[index]), value) for index, value in local_entries]
            centroid = nuclei["centroid_uv"][local_i].astype(np.float64)
            extent = float(np.linalg.norm(polygon - centroid, axis=1).max(initial=0))
            direct_candidates = all_tree.query_ball_point(
                centroid, radius + extent + 1e-10
            )
            total_area = core.polygon_area(polygon)
            captured = 0.0
            for spot in direct_candidates:
                captured += circle_polygon_intersection_area(
                    polygon, centers_uv[int(spot)], radius
                ) / total_area
            direct_fraction[local_i] = np.float32(min(max(captured, 0.0), 1.0))
        if not entries:
            entries = [(int(owner[local_i]), 1.0)]
            fallback[local_i] = True
            distance = float(
                np.linalg.norm(
                    nuclei["centroid_uv"][local_i] - centers_uv[int(owner[local_i])]
                )
            )
            direct_fraction[local_i] = np.float32(1.0 if distance <= radius else 0.0)
        if direct_fraction[local_i] > 1e-7 and evidence[local_i] < 2:
            evidence[local_i] = 0
        for spot, value in entries:
            link_spot.append(spot)
            link_cell.append(local_i)
            link_fraction.append(value)
        if time.time() - last > 60:
            core.log(
                f"Computed full-domain Voronoi contour overlap for "
                f"{local_i + 1:,}/{len(source):,} real nuclei"
            )
            last = time.time()

    order = np.lexsort((np.asarray(link_cell), np.asarray(link_spot)))
    spots = np.asarray(link_spot, np.int32)[order]
    cells = np.asarray(link_cell, np.int32)[order]
    fractions = np.asarray(link_fraction, np.float32)[order]
    indptr = np.zeros(len(centers_uv) + 1, np.int64)
    np.cumsum(np.bincount(spots, minlength=len(centers_uv)), out=indptr[1:])
    geometry = sparse.csr_matrix(
        (fractions, cells, indptr), shape=(len(centers_uv), len(source))
    )
    csc = geometry.tocsc()
    allocation_sum = np.asarray(geometry.sum(axis=0)).ravel()
    nuclei["geometry_fallback"] = fallback
    nuclei["measured_overlap_fraction_sum"] = direct_fraction
    nuclei["direct_capture_overlap_fraction_sum"] = direct_fraction
    nuclei["overlap_bin_count"] = np.diff(csc.indptr).astype(np.int16)
    nuclei["spatial_evidence_level"] = evidence
    nuclei["spatial_evidence_categories"] = [
        "direct_55um_capture_overlap",
        "within_component_voronoi_interpolation",
        "cross_component_nearest_spot_fallback",
    ]
    audit = {
        "geometry": (
            "all-real-nucleus contour overlap with within-component local "
            "Voronoi cells of accepted virtual 55um spots"
        ),
        "allocation_weight_sum_max_error": float(
            np.max(np.abs(allocation_sum - 1.0), initial=0)
        ),
        "n_overlap_links": int(geometry.nnz),
        "n_multispot_cells": int((np.diff(csc.indptr) > 1).sum()),
        "max_spots_per_cell": int(np.diff(csc.indptr).max(initial=0)),
        "n_invalid_contour_fallback_cells": int(fallback.sum()),
        "n_direct_capture_overlap_cells": int((evidence == 0).sum()),
        "n_within_component_interpolated_cells": int((evidence == 1).sum()),
        "n_cross_component_fallback_cells": int((evidence == 2).sum()),
        "direct_capture_fraction_min": float(direct_fraction.min(initial=0)),
        "direct_capture_fraction_median": float(np.median(direct_fraction)),
        "direct_capture_fraction_max": float(direct_fraction.max(initial=0)),
    }
    return geometry, audit


def build_circle_contour_overlap(
    nuclei: dict,
    affine: dict,
    centers_uv: np.ndarray,
    radius: float,
) -> tuple[sparse.csr_matrix, dict]:
    inverse = affine["inverse"]
    intercept = affine["intercept"]
    offsets = nuclei.pop("all_contour_offsets")
    contours = nuclei.pop("all_contour_xy")
    source = nuclei["source_index"]
    owner = nuclei["owner_bin"]
    tree = cKDTree(centers_uv)
    link_spot: list[int] = []
    link_cell: list[int] = []
    link_fraction: list[float] = []
    fallback = np.zeros(len(source), dtype=bool)
    last = time.time()
    for local_i, source_i in enumerate(source):
        start, end = int(offsets[source_i]), int(offsets[source_i + 1])
        poly_xy = contours[start:end]
        valid_poly = len(poly_xy) >= 3 and np.isfinite(poly_xy).all()
        entries: list[tuple[int, float]] = []
        if valid_poly:
            poly = (poly_xy - intercept) @ inverse.T
            total_area = core.polygon_area(poly)
            if total_area > 1e-10:
                centroid = nuclei["centroid_uv"][local_i].astype(np.float64)
                extent = float(np.linalg.norm(poly - centroid, axis=1).max(initial=0))
                candidates = tree.query_ball_point(centroid, radius + extent + 1e-10)
                for spot in candidates:
                    area = circle_polygon_intersection_area(
                        poly, centers_uv[int(spot)], radius
                    )
                    if area > max(1e-10, total_area * 1e-7):
                        entries.append((int(spot), area / total_area))
        if not entries:
            entries = [(int(owner[local_i]), 1.0)]
            fallback[local_i] = True
        fraction_sum = sum(value for _, value in entries)
        if fraction_sum > 1.0 + 1e-7:
            entries = [(spot, value / fraction_sum) for spot, value in entries]
        for spot, value in entries:
            link_spot.append(spot)
            link_cell.append(local_i)
            link_fraction.append(value)
        if time.time() - last > 60:
            core.log(
                f"Computed virtual-circle contour overlap for "
                f"{local_i + 1:,}/{len(source):,} eligible nuclei"
            )
            last = time.time()
    order = np.lexsort((np.asarray(link_cell), np.asarray(link_spot)))
    spots = np.asarray(link_spot, np.int32)[order]
    cells = np.asarray(link_cell, np.int32)[order]
    fractions = np.asarray(link_fraction, np.float32)[order]
    indptr = np.zeros(len(centers_uv) + 1, np.int64)
    np.cumsum(np.bincount(spots, minlength=len(centers_uv)), out=indptr[1:])
    geometry = sparse.csr_matrix(
        (fractions, cells, indptr), shape=(len(centers_uv), len(source))
    )
    csc = geometry.tocsc()
    overlap_sum = np.asarray(geometry.sum(axis=0)).ravel()
    nuclei["geometry_fallback"] = fallback
    nuclei["measured_overlap_fraction_sum"] = overlap_sum.astype(np.float32)
    nuclei["overlap_bin_count"] = np.diff(csc.indptr).astype(np.int16)
    audit = {
        "geometry": "nucleus_contour_area_overlap_with_accepted_55um_circles",
        "n_overlap_links": int(geometry.nnz),
        "n_multispot_cells": int((np.diff(csc.indptr) > 1).sum()),
        "max_spots_per_cell": int(np.diff(csc.indptr).max(initial=0)),
        "n_geometry_fallback_cells": int(fallback.sum()),
        "min_measured_overlap_fraction_sum": (
            float(overlap_sum.min()) if len(overlap_sum) else 0.0
        ),
        "median_measured_overlap_fraction_sum": (
            float(np.median(overlap_sum)) if len(overlap_sum) else 0.0
        ),
        "max_measured_overlap_fraction_sum": (
            float(overlap_sum.max(initial=0)) if len(overlap_sum) else 0.0
        ),
    }
    return geometry, audit


def hex_components(positions: pd.DataFrame) -> np.ndarray:
    rows = positions["hex_row"].to_numpy(np.int32)
    cols = positions["hex_col"].to_numpy(np.int32)
    lookup = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(rows, cols))}
    edge_r: list[int] = []
    edge_c: list[int] = []
    for i, (row, col) in enumerate(zip(rows, cols)):
        if row % 2 == 0:
            neighbors = (
                (row, col - 1), (row, col + 1),
                (row - 1, col - 1), (row - 1, col),
                (row + 1, col - 1), (row + 1, col),
            )
        else:
            neighbors = (
                (row, col - 1), (row, col + 1),
                (row - 1, col), (row - 1, col + 1),
                (row + 1, col), (row + 1, col + 1),
            )
        for neighbor in neighbors:
            j = lookup.get((int(neighbor[0]), int(neighbor[1])))
            if j is not None:
                edge_r.append(i)
                edge_c.append(j)
    graph = sparse.coo_matrix(
        (np.ones(len(edge_r), np.uint8), (edge_r, edge_c)),
        shape=(len(rows), len(rows)),
    ).tocsr()
    _, labels = sparse.csgraph.connected_components(graph, directed=False)
    return labels.astype(np.int32)


def atomic_spot_positions(path: Path, positions: pd.DataFrame) -> None:
    partial = path.with_suffix(path.suffix + ".partial")
    if path.exists() or partial.exists():
        raise FileExistsError(f"Refusing to overwrite spot positions: {path}")
    positions.to_csv(partial, sep="\t", index=False, compression="gzip")
    os.replace(partial, path)


def append_degradation_provenance(
    path: Path,
    source_matrix: sparse.csc_matrix,
    source_barcodes: np.ndarray,
    source_bin_index: np.ndarray,
    overlap: sparse.csr_matrix,
    spot_positions: pd.DataFrame,
    gap_gene: np.ndarray,
    geometry_audit: dict,
) -> None:
    captured = np.asarray(overlap.sum(axis=0)).ravel()
    gap_weight = np.maximum(1.0 - captured, 0.0)
    with h5py.File(path, "r+") as h:
        h.attrs["schema"] = "visium_hd_virtual_visium55_cell_allocation_v4"
        h.attrs["observation_geometry"] = (
            "55um_diameter_circles_100um_pitch_hexagonal"
        )
        h.attrs["source_resolution_um"] = SOURCE_BIN_UM
        degradation = h.create_group("source_016um_to_virtual_visium55")
        degradation.attrs["format"] = "spot_by_positive_source_bin_csr"
        degradation.attrs["geometry_json"] = json.dumps(
            geometry_audit, ensure_ascii=False
        )
        degradation.create_dataset(
            "shape", data=np.asarray(overlap.shape, np.int64)
        )
        degradation.create_dataset(
            "indptr", data=overlap.indptr, compression="lzf"
        )
        degradation.create_dataset(
            "source_bin_index", data=overlap.indices, compression="lzf"
        )
        degradation.create_dataset(
            "area_fraction", data=overlap.data, compression="lzf"
        )
        degradation.create_dataset(
            "positive_source_bin_index", data=source_bin_index, compression="lzf"
        )
        core.write_text_dataset(degradation, "source_barcode", source_barcodes)
        degradation.create_dataset(
            "spot_coverage_fraction",
            data=spot_positions["coverage_fraction"].to_numpy(np.float64),
            compression="lzf",
        )
        degradation.create_dataset(
            "spot_hex_row",
            data=spot_positions["hex_row"].to_numpy(np.int32),
            compression="lzf",
        )
        degradation.create_dataset(
            "spot_hex_col",
            data=spot_positions["hex_col"].to_numpy(np.int32),
            compression="lzf",
        )
        degradation.create_dataset(
            "spot_center_uv",
            data=spot_positions[["array_col", "array_row"]].to_numpy(np.float64),
            compression="lzf",
        )
        degradation.create_dataset(
            "source_bin_gap_fraction", data=gap_weight, compression="lzf"
        )
        degradation.create_dataset(
            "gap_gene_mass", data=gap_gene, compression="lzf"
        )
        source_gene = np.asarray(source_matrix.sum(axis=1)).ravel().astype(np.float64)
        observed = h["observed"]
        virtual_shape = tuple(int(x) for x in observed["shape"][:])
        virtual_matrix = sparse.csc_matrix(
            (
                observed["count"][:],
                observed["gene_index"][:],
                observed["bin_indptr"][:],
            ),
            shape=virtual_shape,
        )
        virtual_gene = np.asarray(virtual_matrix.sum(axis=1)).ravel().astype(np.float64)
        degradation.attrs["max_gene_mass_closure_error"] = float(
            np.max(np.abs(source_gene - virtual_gene - gap_gene))
        )
        h.flush()


def update_output_schema(path: Path) -> None:
    with h5py.File(path, "r+") as h:
        h.attrs["schema"] = "virtual_visium55_spot2cell_mv2_v1"
        h.attrs["X_semantics"] = (
            "fractional counts constrained by virtual 55um spot mass"
        )
        h.flush()


def main() -> None:
    args = parse_args()
    paths = source_paths(args)
    required = {
        **paths,
        "positions_tsv": args.positions_tsv,
        "marker_catalog": args.marker_catalog,
    }
    if args.eligible_cell_h5ad is not None:
        required["eligible_cell_h5ad"] = args.eligible_cell_h5ad
    if args.calibration_cell_h5ad is not None:
        required["calibration_cell_h5ad"] = args.calibration_cell_h5ad
    for name, path in required.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing {name}: {path}")
    output_partial = args.output.with_suffix(args.output.suffix + ".partial")
    allocation_partial = args.allocation_output.with_suffix(
        args.allocation_output.suffix + ".partial"
    )
    formal_paths = (
        args.output, args.allocation_output, args.qc_json,
        args.spot_positions_output, output_partial, allocation_partial,
    )
    if not args.audit_only:
        for path in formal_paths:
            if path.exists():
                raise FileExistsError(f"Refusing to overwrite existing output: {path}")
        for path in (
            args.output, args.allocation_output, args.qc_json,
            args.spot_positions_output,
        ):
            path.parent.mkdir(parents=True, exist_ok=True)

    (
        source_matrix, source_library, source_barcodes,
        gene_id, gene_name, feature_type, source_bin_index, zero_count,
    ) = core.load_matrix(paths["matrix"])
    with h5py.File(paths["cell_features"], "r") as handle:
        all_centroid = handle["centroid_xy"][:].astype(np.float32)
        all_cell_id = handle["cell_id"][:].astype(np.int64)
    (
        source_positions, _, _, centroid_uv_all, affine, spatial_audit,
    ) = core.align_positions(args.positions_tsv, source_barcodes, all_centroid)
    del all_centroid

    (
        matrix, source_to_spot, spot_positions, gap_gene, degradation_audit,
    ) = build_virtual_spots(
        source_matrix, source_positions, affine, args.sample,
        args.coverage_threshold,
    )
    barcodes = spot_positions["barcode"].to_numpy(dtype="U")
    library = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    centers_uv = spot_positions[["array_col", "array_row"]].to_numpy(np.float64)
    if args.eligible_cell_h5ad is None:
        owner_all, eligible = virtual_cell_ownership(
            centroid_uv_all, centers_uv, SPOT_RADIUS_UM / SOURCE_BIN_UM
        )
        component_all = np.full(len(eligible), -1, np.int32)
        cross_component_all = np.zeros(len(eligible), dtype=bool)
        spot_component = np.zeros(len(centers_uv), np.int32)
    else:
        source_component = core.measured_components(source_positions)
        spot_component = spot_source_components(source_to_spot, source_component)
        eligible, component_all = full_real_nucleus_domain(
            args.eligible_cell_h5ad, all_cell_id, source_component
        )
        owner_all, cross_component_all = nearest_spot_within_component(
            centroid_uv_all, eligible, component_all, centers_uv, spot_component
        )
    nuclei = core.load_eligible_nuclei(
        paths["cell_features"], owner_all, eligible, centroid_uv_all
    )
    nuclei["domain_component"] = component_all[eligible]
    if args.eligible_cell_h5ad is None:
        geometry, geometry_audit = build_circle_contour_overlap(
            nuclei, affine, centers_uv, SPOT_RADIUS_UM / SOURCE_BIN_UM
        )
    else:
        geometry, geometry_audit = build_full_domain_voronoi_overlap(
            nuclei, affine, centers_uv, spot_component,
            cross_component_all, SPOT_RADIUS_UM / SOURCE_BIN_UM,
        )
    del (
        centroid_uv_all, owner_all, eligible, component_all,
        cross_component_all, all_cell_id,
    )
    geometry_audit["estimated_allocation_responsibility_records"] = int(
        np.dot(
            np.diff(matrix.indptr).astype(np.int64),
            np.diff(geometry.indptr).astype(np.int64),
        )
    )
    base = core.base_capacity(nuclei)
    base_composition, base_mass, n_cells_per_spot = core.census_from_geometry(
        geometry, nuclei["class_id"], base
    )
    occupied = n_cells_per_spot > 0
    summary = {
        "schema": "virtual_visium55_spot2cell_mv2_summary_v1",
        "algorithm": (
            "measured 16um mass -> exact 55um circle integration -> "
            "independent CellViT census/composition/local-state fit -> continuous "
            "subspot interpolation -> morphology residual -> two-stage KL allocation -> "
            + (
                "all-real-nucleus within-component Voronoi allocation -> "
                "exact occupied-spot gene projection"
                if args.eligible_cell_h5ad is not None else
                "exact occupied-spot gene projection"
            )
        ),
        "sample": args.sample,
        "independent_from_16um_fitted_model": True,
        "independent_sample_fit": True,
        "no_scrna_reference": True,
        "xenium_used_for_fit": False,
        "ai_annotation_used_for_fit": False,
        "fixed_44_programs_used_downstream": False,
        "real_nucleus_domain_rule": (
            "every CellViT nucleus in the positive-count 16um measured domain "
            "receives expression; no synthetic or outside-domain nuclei"
            if args.eligible_cell_h5ad is not None else
            "only nuclei with physical overlap to accepted 55um circles"
        ),
        "spot_gap_interpolation": (
            "within-component local Voronoi contour overlap"
            if args.eligible_cell_h5ad is not None else "disabled"
        ),
        "eligible_cell_h5ad": (
            str(args.eligible_cell_h5ad)
            if args.eligible_cell_h5ad is not None else None
        ),
        "occupied_spot_background_component": False,
        "within_type_temperature": float(args.within_type_temperature),
        "subspot_radius_spots": float(args.subspot_radius_spots),
        "source_positive_16um_bins": int(source_matrix.shape[1]),
        "source_zero_count_filtered_bins_excluded": int(zero_count),
        "n_genes": int(matrix.shape[0]),
        "n_accepted_virtual_spots": int(matrix.shape[1]),
        "n_eligible_cells": int(geometry.shape[1]),
        "n_spots_with_overlapping_eligible_cells": int(occupied.sum()),
        "n_spots_without_overlapping_eligible_cells": int((~occupied).sum()),
        "spatial_audit": spatial_audit,
        "degradation_audit": degradation_audit,
        "geometry_audit": geometry_audit,
        "class_cell_counts": {
            core.CLASS_NAMES[int(value)]: int(count)
            for value, count in zip(
                *np.unique(nuclei["class_id"], return_counts=True)
            )
        },
    }
    core.log(json.dumps(summary, ensure_ascii=False))
    if args.audit_only:
        return

    marker_sets, catalog = core.load_marker_sets(args.marker_catalog, gene_name)
    occupied_gene = np.asarray(matrix @ occupied.astype(np.float64)).ravel()
    global_mean = occupied_gene / max(float(occupied_gene.sum()), 1.0)
    prior = core.marker_prior(global_mean, gene_name, marker_sets)
    fit_genes = core.select_fit_genes(matrix, gene_name, marker_sets)
    blocks = core.spatial_blocks(spot_positions)
    alpha, independent_profiles, alpha_cv = core.fit_alpha_and_profiles(
        matrix, library, base_mass, prior, fit_genes, blocks,
        args.class_prior_strength, args.cv_folds, args.seed,
    )
    capacity = (
        base.astype(np.float64) * alpha[nuclei["class_id"]]
    ).astype(np.float32)
    composition, census_mass, n_cells_per_spot = core.census_from_geometry(
        geometry, nuclei["class_id"], capacity
    )
    shrinkage, shrinkage_audit = core.hierarchy_shrinkage_cv(
        matrix, library, composition, prior, fit_genes, blocks,
        args.class_prior_strength, args.cv_folds,
    )
    profiles = core.apply_hierarchy_shrinkage(
        independent_profiles, shrinkage, composition, global_mean
    )

    core.log("Learning class-internal CellViT features for the 55 um branch")
    cell_features, embedding_scores, embedding_pca = core.learn_cell_features(
        paths["cell_features"], nuclei, args.embedding_pcs,
        args.pca_sample_per_class, args.seed,
    )
    groups = core.aggregate_bin_type_features(
        geometry, nuclei["class_id"], capacity, cell_features
    )
    class_audit = core.profile_audit(
        composition, groups, profiles, gene_name, marker_sets, alpha, shrinkage
    )
    no_marker_prior = np.tile(global_mean[None, :], (core.N_CLASSES, 1))
    no_marker_profiles = core.fit_profiles(
        matrix, composition, library, no_marker_prior,
        args.class_prior_strength,
    ).astype(np.float32)
    no_marker_audit = core.profile_audit(
        composition, groups, no_marker_profiles, gene_name, marker_sets,
        alpha, np.ones(core.N_CLASSES, np.float32),
    )
    summary.update(
        {
            "marker_catalog_schema_version": catalog.get("schema_version"),
            "marker_catalog_workbook_sha256": catalog.get("source_workbook_sha256"),
            "marker_genes_present": {
                name: len(genes) for name, genes in marker_sets.items()
            },
            "alpha": alpha.tolist(),
            "alpha_cv": alpha_cv,
            "hierarchy_shrinkage_cv": shrinkage_audit,
            "class_fit_audit": class_audit,
        }
    )

    program_genes, residual_genes, variability = core.select_program_genes(
        matrix, library, occupied, gene_name, marker_sets,
        args.program_genes, args.residual_genes, args.seed,
    )
    graph_genes = residual_genes[::2]
    state_fit_genes = residual_genes[1::2]
    program_genes = np.setdiff1d(program_genes, graph_genes, assume_unique=False)
    summary["program_gene_count"] = int(len(program_genes))
    summary["residual_graph_gene_count"] = int(len(graph_genes))
    summary["state_selection_gene_count"] = int(len(state_fit_genes))
    summary["gene_split_crossfit"] = (
        "deterministic disjoint graph-gene and state-selection-gene halves"
    )
    residual_scores, residual_pca = core.compute_residual_pcs(
        matrix, library, composition, profiles, graph_genes,
        occupied, args.residual_pcs, args.seed,
    )
    component_id = hex_components(spot_positions)
    summary["measured_domain_connected_components"] = int(
        len(np.unique(component_id))
    )
    group_scores, loadings, program_count, local_audit = core.learn_local_programs(
        matrix, library, spot_positions, component_id, composition, profiles,
        residual_scores, state_fit_genes, program_genes, groups,
        args.graph_neighbors, args.graph_radius_bins, args.cv_folds,
        args.bootstrap_replicates, args.seed,
    )
    morphology_link_scores, feature_coef, reliability, feature_audit = (
        core.learn_feature_program_mapping(
            nuclei, capacity, cell_features, groups, composition, residual_scores,
            spot_positions, group_scores, program_count, args.permutations,
            args.cv_folds, args.seed + 7000,
        )
    )
    local_state_reliability = core.local_program_reliability(local_audit)
    group_link_class = nuclei["class_id"][groups["link_cell"]]
    raw_group_link_scores = group_scores[groups["link_group"]].astype(
        np.float32, copy=False
    )
    group_link_scores = (
        raw_group_link_scores
        * local_state_reliability[group_link_class, None]
    ).astype(np.float32)
    morphology_offsets = (
        morphology_link_scores.astype(np.float32, copy=False) - raw_group_link_scores
    ).astype(np.float32)
    spatial_offsets, subspot_audit = core.interpolate_subspot_offsets(
        nuclei, geometry, capacity, groups, spot_positions, component_id,
        composition, residual_scores, group_scores, program_count, local_audit,
        args.subspot_neighbors, args.subspot_radius_bins,
    )
    spatial_scale, morphology_scale, calibration_audit = core.calibrate_offset_scales(
        args.calibration_cell_h5ad, args.calibration_max_cells,
        nuclei, geometry, capacity, spot_positions, gene_name, profiles, loadings,
        program_count, group_link_scores, spatial_offsets, morphology_offsets,
        args.subspot_spatial_scale, args.morphology_offset_scale,
        args.seed + 17000,
    )
    link_scores = (
        group_link_scores
        + spatial_scale * spatial_offsets
        + morphology_scale * morphology_offsets
    ).astype(np.float32)
    summary["local_type_program_audit"] = local_audit
    summary["nucleus_feature_program_audit"] = feature_audit
    summary["subspot_spatial_interpolation_audit"] = subspot_audit
    summary["heterogeneity_calibration"] = calibration_audit
    summary["selected_subspot_spatial_scale"] = float(spatial_scale)
    summary["selected_morphology_offset_scale"] = float(morphology_scale)
    summary["global_profile_bootstrap"] = core.global_profile_bootstrap(
        matrix, library, composition, prior, profiles, fit_genes, blocks,
        args.class_prior_strength, args.bootstrap_replicates, args.seed + 9000,
    )
    if args.skip_ablations:
        ablations = []
        summary["ablations_skipped"] = True
    else:
        ablations = core.run_ablations(
            matrix, library, spot_positions, component_id, composition, profiles,
            residual_scores, state_fit_genes, groups, class_audit,
            no_marker_audit, feature_audit, args.graph_neighbors,
            args.graph_radius_bins, args.cv_folds, args.seed + 11000,
        )
        summary["ablations"] = ablations

    core.log(f"Hashing immutable {args.sample} inputs for virtual-Visium provenance")
    provenance_paths = {
        "source_016um_matrix": paths["matrix"],
        "source_016um_positions": args.positions_tsv,
        "cell_features": paths["cell_features"],
        "marker_catalog": args.marker_catalog,
    }
    if args.calibration_cell_h5ad is not None:
        provenance_paths["heterogeneity_calibration_reference"] = (
            args.calibration_cell_h5ad
        )
    summary["input_provenance"] = {
        name: {
            "path": str(path),
            "bytes": int(path.stat().st_size),
            "sha256": core.sha256_file(path),
        }
        for name, path in provenance_paths.items()
    }
    summary["outputs"] = {
        "cell_resolved_h5ad": str(args.output),
        "allocation_h5": str(args.allocation_output),
        "qc_json": str(args.qc_json),
        "virtual_spot_positions": str(args.spot_positions_output),
    }

    cell_total, cell_program, conservation = core.write_allocation_file(
        allocation_partial, args, paths, matrix, library, barcodes,
        np.arange(matrix.shape[1], dtype=np.int32),
        gene_id, gene_name, geometry, nuclei, capacity,
        link_scores, profiles, loadings, program_count, summary,
    )
    summary["count_conservation"] = conservation
    append_degradation_provenance(
        allocation_partial, source_matrix, source_barcodes, source_bin_index,
        source_to_spot, spot_positions, gap_gene, degradation_audit,
    )
    nnz = core.write_h5ad(
        output_partial, allocation_partial, args, matrix, barcodes,
        gene_id, gene_name, feature_type, spot_positions, geometry, nuclei,
        capacity, cell_features, embedding_scores, cell_total, cell_program,
        profiles, loadings, program_count, reliability, shrinkage,
        marker_sets, program_genes, graph_genes, groups, group_scores,
        embedding_pca, residual_pca, summary,
    )
    update_output_schema(output_partial)
    summary["cell_matrix_nnz"] = nnz
    validation = core.validate_output_files(
        output_partial, allocation_partial, geometry.shape[1], matrix.shape[0]
    )
    validation["virtual_spot_gene_closure_error"] = (
        degradation_audit["max_gene_mass_closure_error"]
    )
    validation["outside_accepted_circle_cells"] = 0
    summary["output_validation"] = validation
    core.write_qc_files(
        args.qc_json, summary, class_audit, local_audit, feature_audit, ablations
    )
    atomic_spot_positions(args.spot_positions_output, spot_positions)
    os.replace(allocation_partial, args.allocation_output)
    os.replace(output_partial, args.output)
    core.log(f"Completed virtual Visium 55 um {args.sample}: {args.output}")


if __name__ == "__main__":
    main()
