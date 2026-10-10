#!/usr/bin/env python3
"""CellViT-aware reconstruction of Visium HD bins at nucleus resolution.

The model is deliberately reference-free: measured bin counts, the fixed
CellViT hierarchy/nuclear features, and COAD plus universal-TME positive marker
sets are the only inputs.  Every observed gene in a bin with an eligible
overlapping nucleus is assigned to those nuclei; zero-nucleus measured bins are
kept unassigned.  No expression is extrapolated outside the measured domain.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import os
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import optimize, sparse
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold


CLASS_NAMES = [
    "Neoplastic", "Dead", "Epithelial", "Connective_tissue",
    "Neutrophil", "Lymphocyte", "Plasma", "Eosinophil",
    "Inflammatory_unresolved",
]
N_CLASSES = len(CLASS_NAMES)
MAX_PROGRAMS = 4
EPS = 1e-15
BIN_SIZE_UM = 16
SPATIAL_BLOCK_SIZE_BINS = 64

PARENT_MEMBERS = {
    0: (0, 2),       # epithelial lineage
    1: (),           # damaged: global measured profile
    2: (0, 2),
    3: (3,),
    4: (4, 7),       # granulocyte
    5: (5, 6),       # lymphoid
    6: (5, 6),
    7: (4, 7),
    8: (4, 5, 6, 7, 8),  # immune parent
}


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--sample", choices=("P1", "P2", "P5"), required=True)
    p.add_argument("--data-root", type=Path, default=Path('DATA/Visium_HD'))
    p.add_argument("--bin-size-um", type=int, choices=(8, 16), default=16)
    p.add_argument("--positions-tsv", type=Path, required=True)
    p.add_argument(
        "--marker-catalog", type=Path,
        default=Path('DATA/marker_catalog.json'),
    )
    p.add_argument("--output", type=Path, required=True, help="Output .cell_resolved.h5ad")
    p.add_argument("--allocation-output", type=Path)
    p.add_argument("--qc-json", type=Path)
    p.add_argument("--embedding-pcs", type=int, default=16)
    p.add_argument("--pca-sample-per-class", type=int, default=5000)
    p.add_argument("--program-genes", type=int, default=3000)
    p.add_argument("--residual-genes", type=int, default=768)
    p.add_argument("--residual-pcs", type=int, default=16)
    p.add_argument("--graph-neighbors", type=int, default=24)
    p.add_argument(
        "--graph-radius-bins", type=float,
        help="Maximum graph radius in bins; default preserves the original 96 um radius",
    )
    p.add_argument("--class-prior-strength", type=float, default=0.08)
    p.add_argument("--cv-folds", type=int, default=5)
    p.add_argument("--permutations", type=int, default=100)
    p.add_argument("--bootstrap-replicates", type=int, default=20)
    p.add_argument("--cell-batch-size", type=int, default=1024)
    p.add_argument(
        "--calibration-cell-h5ad", type=Path,
        help=(
            "Optional same-section 2 um nucleus-indexed quasi-reference used only "
            "to calibrate the magnitude of centered sub-spot spatial and morphology "
            "offsets. It is never used to fit class profiles or local state axes."
        ),
    )
    p.add_argument("--calibration-max-cells", type=int, default=12000)
    p.add_argument("--subspot-neighbors", type=int, default=9)
    p.add_argument(
        "--subspot-radius-bins", type=float,
        help="Maximum radius for continuous cell-position interpolation; defaults to the local graph radius.",
    )
    p.add_argument("--subspot-spatial-scale", type=float, default=0.5)
    p.add_argument("--morphology-offset-scale", type=float, default=0.5)
    p.add_argument("--within-type-temperature", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=20260826)
    p.add_argument("--audit-only", action="store_true")
    p.add_argument("--skip-ablations", action="store_true")
    args = p.parse_args()
    if args.subspot_neighbors < 1:
        p.error("--subspot-neighbors must be positive")
    if args.within_type_temperature <= 0:
        p.error("--within-type-temperature must be positive")
    return args


def input_paths(args: argparse.Namespace) -> dict[str, Path]:
    stem = f"Visium_HD_Human_Colon_Cancer_{args.sample}_tissue_image"
    root = (
        args.data_root / args.sample / "binned_outputs"
        / f"square_{args.bin_size_um:03d}um"
    )
    return {
        "matrix": root / "filtered_feature_bc_matrix.h5",
        "positions_parquet": root / "spatial" / "tissue_positions.parquet",
        "cell_features": args.data_root / "cellvit_pp_outputs" / args.sample / f"{stem}_cell_features.h5",
    }


def text_array(values) -> np.ndarray:
    a = np.asarray(values)
    if a.dtype.kind == "S":
        return a.astype("U")
    if a.dtype.kind == "O":
        return np.asarray([
            x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in a
        ], dtype="U")
    return a.astype("U")


def sha256_file(path: Path, block: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            data = handle.read(block)
            if not data:
                break
            digest.update(data)
    return digest.hexdigest()


def load_matrix(path: Path):
    log(f"Loading measured matrix: {path}")
    with h5py.File(path, "r") as h:
        g = h["matrix"]
        data = g["data"][:]
        indices = g["indices"][:].astype(np.int32, copy=False)
        indptr = g["indptr"][:].astype(np.int64, copy=False)
        shape = tuple(int(x) for x in g["shape"][:])
        barcodes = text_array(g["barcodes"][:])
        f = g["features"]
        gene_id = text_array(f["id"][:])
        gene_name = text_array(f["name"][:])
        feature_type = text_array(f["feature_type"][:])
    matrix = sparse.csc_matrix((data, indices, indptr), shape=shape)
    library = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    positive = library > 0
    source_bin_index = np.flatnonzero(positive).astype(np.int32)
    zero_count = int((~positive).sum())
    if zero_count:
        log(f"Excluding {zero_count} filtered barcodes with zero measured counts")
        matrix = matrix[:, positive].tocsc()
        barcodes = barcodes[positive]
        library = library[positive]
        gc.collect()
    return (
        matrix, library, barcodes, gene_id, gene_name, feature_type,
        source_bin_index, zero_count,
    )


def align_positions(path: Path, barcodes: np.ndarray, centroid_xy: np.ndarray):
    pos = pd.read_csv(path, sep="\t", compression="gzip")
    required = [
        "barcode", "in_tissue", "array_row", "array_col",
        "pxl_row_in_fullres", "pxl_col_in_fullres",
    ]
    missing = sorted(set(required) - set(pos.columns))
    if missing:
        raise RuntimeError(f"Missing position columns: {missing}")
    if pos["barcode"].duplicated().any():
        raise RuntimeError("Duplicate barcodes in position table")
    pos = pos.set_index("barcode", drop=False)
    absent = [bc for bc in barcodes if bc not in pos.index]
    if absent:
        raise RuntimeError(f"{len(absent)} measured barcodes lack spatial positions")
    pos = pos.loc[barcodes].reset_index(drop=True)
    if not np.all(pos["in_tissue"].to_numpy() == 1):
        raise RuntimeError("A positive measured barcode is not in_tissue=1")

    col = pos["array_col"].to_numpy(np.float64)
    row = pos["array_row"].to_numpy(np.float64)
    design = np.column_stack([np.ones(len(pos)), col, row])
    cx = np.linalg.lstsq(
        design, pos["pxl_col_in_fullres"].to_numpy(np.float64), rcond=None
    )[0]
    cy = np.linalg.lstsq(
        design, pos["pxl_row_in_fullres"].to_numpy(np.float64), rcond=None
    )[0]
    rx = pos["pxl_col_in_fullres"].to_numpy() - design @ cx
    ry = pos["pxl_row_in_fullres"].to_numpy() - design @ cy
    max_residual = float(np.max(np.abs(np.concatenate([rx, ry]))))
    transform = np.array([[cx[1], cx[2]], [cy[1], cy[2]]], np.float64)
    intercept = np.array([cx[0], cy[0]], np.float64)
    if abs(np.linalg.det(transform)) < 1e-6 or max_residual > 2.0:
        raise RuntimeError(f"Invalid array/image affine; max residual={max_residual:.3f}px")

    inverse = np.linalg.inv(transform)
    centroid_uv = (centroid_xy.astype(np.float64) - intercept) @ inverse.T
    nearest_col = np.rint(centroid_uv[:, 0]).astype(np.int32)
    nearest_row = np.rint(centroid_uv[:, 1]).astype(np.int32)
    inside_square = (
        (np.abs(centroid_uv[:, 0] - nearest_col) <= 0.5 + 1e-7)
        & (np.abs(centroid_uv[:, 1] - nearest_row) <= 0.5 + 1e-7)
    )
    pr = pos["array_row"].to_numpy(np.int32)
    pc = pos["array_col"].to_numpy(np.int32)
    min_row = min(int(pr.min()), int(nearest_row.min(initial=0)))
    min_col = min(int(pc.min()), int(nearest_col.min(initial=0)))
    row_shift = -min(0, min_row)
    col_shift = -min(0, min_col)
    max_row = int(max(pr.max(), nearest_row.max(initial=0))) + row_shift
    max_col = int(max(pc.max(), nearest_col.max(initial=0))) + col_shift
    lookup = np.full((max_row + 1, max_col + 1), -1, np.int32)
    lookup[pr + row_shift, pc + col_shift] = np.arange(len(pos), dtype=np.int32)
    rr = nearest_row + row_shift
    cc = nearest_col + col_shift
    valid = (
        inside_square & (rr >= 0) & (cc >= 0)
        & (rr < lookup.shape[0]) & (cc < lookup.shape[1])
    )
    owner = np.full(len(centroid_xy), -1, np.int32)
    owner[valid] = lookup[rr[valid], cc[valid]]
    eligible = owner >= 0
    audit = {
        "max_affine_residual_px": max_residual,
        "grid_vector_array_col_px": transform[:, 0].tolist(),
        "grid_vector_array_row_px": transform[:, 1].tolist(),
        "bin_edge_px_mean": float(
            0.5 * (np.linalg.norm(transform[:, 0]) + np.linalg.norm(transform[:, 1]))
        ),
        "n_nuclei_total": int(len(centroid_xy)),
        "n_nuclei_inside_positive_measured_bins": int(eligible.sum()),
        "n_nuclei_excluded_outside_positive_measured_bins": int((~eligible).sum()),
    }
    affine = {
        "transform": transform,
        "inverse": inverse,
        "intercept": intercept,
        "lookup": lookup,
        "row_shift": row_shift,
        "col_shift": col_shift,
    }
    return pos, owner, eligible, centroid_uv, affine, audit


def load_eligible_nuclei(path: Path, owner_all: np.ndarray, eligible: np.ndarray,
                         centroid_uv_all: np.ndarray):
    idx = np.flatnonzero(eligible).astype(np.int64)
    with h5py.File(path, "r") as h:
        labels = json.loads(h.attrs["hierarchical_labels_json"])
        expected = {str(i): name for i, name in enumerate(CLASS_NAMES)}
        if labels != expected:
            raise RuntimeError(f"Unexpected hierarchy: {labels}")
        pannuke_conf = h["pannuke_confidence"][idx].astype(np.float32)
        class_id = h["hierarchical_class_id"][idx].astype(np.uint8)
        lizard = h["lizard_probabilities"][idx].astype(np.float32)
        conditional = lizard[:, [0, 2, 3, 4]]
        conditional /= np.maximum(conditional.sum(1, keepdims=True), 1e-8)
        hierarchical_conf = pannuke_conf.copy()
        for class_value, local_col in ((4, 0), (5, 1), (6, 2), (7, 3)):
            sel = class_id == class_value
            hierarchical_conf[sel] *= conditional[sel, local_col]
        unresolved = class_id == 8
        hierarchical_conf[unresolved] *= np.maximum(
            1.0 - conditional[unresolved].max(1, initial=0), 0.05
        )
        offsets = h["contour_offsets"][:].astype(np.int64)
        contour_xy = h["contour_xy"][:].astype(np.float64)
        return {
            "source_index": idx,
            "owner_bin": owner_all[idx],
            "centroid_xy": h["centroid_xy"][idx].astype(np.float32),
            "centroid_uv": centroid_uv_all[idx].astype(np.float32),
            "cell_id": h["cell_id"][idx].astype(np.int64),
            "class_id": class_id,
            "morphology": h["morphology"][idx].astype(np.float32),
            "pannuke_confidence": pannuke_conf,
            "hierarchical_confidence": hierarchical_conf.astype(np.float32),
            "lizard_probabilities": lizard.astype(np.float16),
            "morphology_names": json.loads(h.attrs["morphology_feature_names_json"]),
            "all_contour_offsets": offsets,
            "all_contour_xy": contour_xy,
        }


def polygon_area(poly: np.ndarray) -> float:
    if len(poly) < 3:
        return 0.0
    x = poly[:, 0]
    y = poly[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))))


def clip_halfplane(poly: np.ndarray, axis: int, bound: float,
                   keep_greater: bool) -> np.ndarray:
    if len(poly) == 0:
        return poly
    result = []
    previous = poly[-1]
    previous_inside = previous[axis] >= bound if keep_greater else previous[axis] <= bound
    for current in poly:
        current_inside = current[axis] >= bound if keep_greater else current[axis] <= bound
        if current_inside != previous_inside:
            delta = current[axis] - previous[axis]
            if abs(delta) > 1e-15:
                t = (bound - previous[axis]) / delta
                result.append(previous + t * (current - previous))
        if current_inside:
            result.append(current)
        previous = current
        previous_inside = current_inside
    return np.asarray(result, dtype=np.float64).reshape((-1, 2))


def rectangle_intersection_area(poly: np.ndarray, col: int, row: int) -> float:
    clipped = poly
    clipped = clip_halfplane(clipped, 0, col - 0.5, True)
    clipped = clip_halfplane(clipped, 0, col + 0.5, False)
    clipped = clip_halfplane(clipped, 1, row - 0.5, True)
    clipped = clip_halfplane(clipped, 1, row + 0.5, False)
    return polygon_area(clipped)


def build_contour_overlap(nuclei: dict, affine: dict, n_bins: int):
    inverse = affine["inverse"]
    intercept = affine["intercept"]
    lookup = affine["lookup"]
    row_shift = affine["row_shift"]
    col_shift = affine["col_shift"]
    offsets = nuclei.pop("all_contour_offsets")
    contours = nuclei.pop("all_contour_xy")
    source = nuclei["source_index"]
    owner = nuclei["owner_bin"]
    link_bin: list[int] = []
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
            total_area = polygon_area(poly)
            if total_area > 1e-8:
                c0 = int(math.ceil(float(poly[:, 0].min()) - 0.5))
                c1 = int(math.floor(float(poly[:, 0].max()) + 0.5))
                r0 = int(math.ceil(float(poly[:, 1].min()) - 0.5))
                r1 = int(math.floor(float(poly[:, 1].max()) + 0.5))
                for row in range(r0, r1 + 1):
                    rr = row + row_shift
                    if rr < 0 or rr >= lookup.shape[0]:
                        continue
                    for col in range(c0, c1 + 1):
                        cc = col + col_shift
                        if cc < 0 or cc >= lookup.shape[1]:
                            continue
                        b = int(lookup[rr, cc])
                        if b < 0:
                            continue
                        area = rectangle_intersection_area(poly, col, row)
                        if area > max(1e-8, total_area * 1e-6):
                            entries.append((b, area / total_area))
        if not entries:
            entries = [(int(owner[local_i]), 1.0)]
            fallback[local_i] = True
        total_fraction = sum(value for _, value in entries)
        if total_fraction > 1.0 + 1e-5:
            entries = [(b, value / total_fraction) for b, value in entries]
        for b, value in entries:
            link_bin.append(b)
            link_cell.append(local_i)
            link_fraction.append(value)
        if time.time() - last > 60:
            log(f"Computed contour overlap for {local_i + 1:,}/{len(source):,} eligible nuclei")
            last = time.time()
    order = np.lexsort((np.asarray(link_cell), np.asarray(link_bin)))
    bins = np.asarray(link_bin, np.int32)[order]
    cells = np.asarray(link_cell, np.int32)[order]
    fractions = np.asarray(link_fraction, np.float32)[order]
    indptr = np.zeros(n_bins + 1, np.int64)
    np.cumsum(np.bincount(bins, minlength=n_bins), out=indptr[1:])
    geometry = sparse.csr_matrix((fractions, cells, indptr), shape=(n_bins, len(source)))
    csc = geometry.tocsc()
    overlap_sum = np.asarray(geometry.sum(axis=0)).ravel()
    audit = {
        "geometry": f"nucleus_contour_area_overlap_with_measured_{BIN_SIZE_UM}um_bins",
        "n_overlap_links": int(geometry.nnz),
        "n_multibin_cells": int((np.diff(csc.indptr) > 1).sum()),
        "max_bins_per_cell": int(np.diff(csc.indptr).max(initial=0)),
        "n_geometry_fallback_cells": int(fallback.sum()),
        "min_measured_overlap_fraction_sum": float(overlap_sum.min()) if len(overlap_sum) else 0.0,
        "median_measured_overlap_fraction_sum": float(np.median(overlap_sum)),
        "max_measured_overlap_fraction_sum": float(overlap_sum.max(initial=0)),
    }
    nuclei["geometry_fallback"] = fallback
    nuclei["measured_overlap_fraction_sum"] = overlap_sum.astype(np.float32)
    nuclei["overlap_bin_count"] = np.diff(csc.indptr).astype(np.int16)
    return geometry, audit


def base_capacity(nuclei: dict) -> np.ndarray:
    class_id = nuclei["class_id"]
    area = nuclei["morphology"][:, 0].astype(np.float64)
    result = np.ones(len(area), np.float64)
    for c in range(N_CLASSES):
        sel = class_id == c
        valid = sel & np.isfinite(area) & (area > 0)
        median = float(np.median(area[valid])) if valid.any() else 1.0
        result[sel] = np.sqrt(np.maximum(area[sel], 1.0) / max(median, 1.0))
    result = np.clip(result, 0.25, 4.0)
    confidence = np.clip(
        np.nan_to_num(nuclei["hierarchical_confidence"], nan=0.05), 0.05, 1.0
    )
    return (result * confidence).astype(np.float32)


def census_from_geometry(geometry: sparse.csr_matrix, class_id: np.ndarray,
                         capacity: np.ndarray):
    link_bins = np.repeat(
        np.arange(geometry.shape[0], dtype=np.int32), np.diff(geometry.indptr)
    )
    link_cells = geometry.indices.astype(np.int32, copy=False)
    link_weight = geometry.data.astype(np.float64) * capacity[link_cells]
    mass = sparse.coo_matrix(
        (link_weight, (link_bins, class_id[link_cells])),
        shape=(geometry.shape[0], N_CLASSES),
    ).tocsr().toarray().astype(np.float64)
    total = mass.sum(1)
    composition = np.divide(
        mass, total[:, None], out=np.zeros_like(mass), where=total[:, None] > 0
    )
    n_cells = np.diff(geometry.indptr).astype(np.int32)
    return composition, mass, n_cells


def robust_z(values: np.ndarray) -> np.ndarray:
    median = np.nanmedian(values, axis=0)
    scale = np.nanmedian(np.abs(values - median), axis=0) * 1.4826
    scale[~np.isfinite(scale) | (scale < 1e-6)] = 1.0
    return np.clip(
        np.nan_to_num((values - median) / scale, nan=0.0, posinf=5.0, neginf=-5.0),
        -5, 5,
    )


def neighborhood_features(nuclei: dict, radius_bins: float = 3.0) -> np.ndarray:
    uv = nuclei["centroid_uv"].astype(np.float64)
    area = np.log1p(np.maximum(nuclei["morphology"][:, 0].astype(np.float64), 0))
    conf = nuclei["hierarchical_confidence"].astype(np.float64)
    classes = nuclei["class_id"]
    result = np.zeros((len(uv), 3), np.float32)
    for c in range(N_CLASSES):
        members = np.flatnonzero(classes == c)
        if not len(members):
            continue
        tree = cKDTree(uv[members])
        k = min(9, len(members))
        distance, index = tree.query(uv[members], k=k, workers=-1)
        if k == 1:
            index = index[:, None]
        neighbors = index[:, 1:] if k > 1 else index[:, :0]
        if neighbors.shape[1]:
            result[members, 1] = np.mean(area[members][neighbors], axis=1)
            result[members, 2] = np.mean(conf[members][neighbors], axis=1)
        try:
            count = tree.query_ball_point(uv[members], radius_bins, return_length=True)
        except TypeError:
            count = np.asarray([len(x) for x in tree.query_ball_point(uv[members], radius_bins)])
        result[members, 0] = np.log1p(np.maximum(np.asarray(count) - 1, 0))
    return result


def learn_cell_features(feature_h5: Path, nuclei: dict, n_components: int,
                        sample_per_class: int, seed: int):
    rng = np.random.default_rng(seed)
    source = nuclei["source_index"]
    classes = nuclei["class_id"]
    embedding_scores = np.zeros((len(source), n_components), np.float32)
    pca_models: list[PCA | None] = [None] * N_CLASSES
    with h5py.File(feature_h5, "r") as h:
        ds = h["cellvit_embedding_float16"]
        for c in range(N_CLASSES):
            members = np.flatnonzero(classes == c)
            if len(members) < 2:
                continue
            fit_members = members
            if len(fit_members) > sample_per_class:
                fit_members = np.sort(rng.choice(fit_members, sample_per_class, replace=False))
            sample = ds[source[fit_members]].astype(np.float32)
            k = min(n_components, sample.shape[0] - 1, sample.shape[1])
            pca = PCA(n_components=k, svd_solver="randomized", random_state=seed + c)
            pca.fit(sample)
            pca_models[c] = pca
            for start in range(0, len(members), 8192):
                part = members[start:start + 8192]
                embedding_scores[part, :k] = pca.transform(
                    ds[source[part]].astype(np.float32)
                ).astype(np.float32)
            del sample

    morph = nuclei["morphology"].astype(np.float64)
    confidence = nuclei["hierarchical_confidence"].astype(np.float64)[:, None]
    neighbor = neighborhood_features(nuclei).astype(np.float64)
    raw = np.concatenate([embedding_scores.astype(np.float64), morph, confidence, neighbor], axis=1)
    features = np.zeros_like(raw, dtype=np.float32)
    log_area = np.log1p(np.maximum(morph[:, 0], 0))
    for c in range(N_CLASSES):
        members = np.flatnonzero(classes == c)
        if not len(members):
            continue
        nuisance = np.column_stack([
            np.ones(len(members)), log_area[members], confidence[members, 0]
        ])
        coef = np.linalg.lstsq(nuisance, raw[members], rcond=None)[0]
        residual = raw[members] - nuisance @ coef
        # Keep confidence itself as an explicit, gated feature after residualising
        # all morphology/embedding/neighbourhood channels against confidence.
        residual[:, n_components + morph.shape[1]] = confidence[members, 0]
        features[members] = robust_z(residual).astype(np.float32)
    return features, embedding_scores, pca_models


def aggregate_bin_type_features(geometry: sparse.csr_matrix, class_id: np.ndarray,
                                capacity: np.ndarray, features: np.ndarray):
    link_bin = np.repeat(
        np.arange(geometry.shape[0], dtype=np.int32), np.diff(geometry.indptr)
    )
    link_cell = geometry.indices.astype(np.int32, copy=False)
    link_class = class_id[link_cell]
    link_weight = geometry.data.astype(np.float64) * capacity[link_cell]
    key = link_bin.astype(np.int64) * N_CLASSES + link_class.astype(np.int64)
    unique, inverse = np.unique(key, return_inverse=True)
    mass = np.bincount(inverse, weights=link_weight, minlength=len(unique)).astype(np.float64)
    count = np.bincount(inverse, minlength=len(unique)).astype(np.int32)
    normalized = link_weight / np.maximum(mass[inverse], EPS)
    aggregator = sparse.coo_matrix(
        (normalized, (inverse, link_cell)), shape=(len(unique), geometry.shape[1])
    ).tocsr()
    mean = np.asarray(aggregator @ features).astype(np.float32)
    second = np.asarray(aggregator @ np.square(features)).astype(np.float32)
    variance = np.maximum(second - np.square(mean), 0).astype(np.float32)
    return {
        "bin_index": (unique // N_CLASSES).astype(np.int32),
        "class_id": (unique % N_CLASSES).astype(np.uint8),
        "mass": mass.astype(np.float32),
        "n_cells": count,
        "feature_mean": mean,
        "feature_variance": variance,
        "link_group": inverse.astype(np.int32),
        "link_bin": link_bin,
        "link_cell": link_cell,
        "link_weight": link_weight.astype(np.float32),
        "aggregator": aggregator,
    }


def load_marker_sets(path: Path, gene_names: np.ndarray):
    catalog = json.loads(path.read_text(encoding="utf-8"))
    sets = {name: set() for name in CLASS_NAMES}
    for record in catalog["tissues"]["COAD"]["organ_marker_records"]:
        category = str(record.get("marker_category", ""))
        genes = set(record.get("positive_markers", []))
        if "谱系" in category:
            sets["Neoplastic"].update(genes)
            sets["Epithelial"].update(genes)
        elif "正常" in category:
            sets["Epithelial"].update(genes)
        elif "肿瘤" in category or "侵袭" in category:
            sets["Neoplastic"].update(genes)
    immune = ["Neutrophil", "Lymphocyte", "Plasma", "Eosinophil", "Inflammatory_unresolved"]
    for record in catalog["universal_tme_records"]:
        rc = str(record.get("class", ""))
        state = str(record.get("cell_type_or_state", ""))
        desc = rc + " " + state
        genes = set(record.get("positive_markers", []))
        targets: list[str] = []
        if rc == "免疫" or "泛白细胞" in desc:
            targets = immune
        elif "中性粒" in desc:
            targets = ["Neutrophil"]
        elif "嗜酸" in desc:
            targets = ["Eosinophil"]
        elif "浆细胞" in desc:
            targets = ["Plasma"]
        elif any(x in rc for x in ("T细胞", "NK细胞", "B细胞")):
            targets = ["Lymphocyte"]
        elif any(x in desc for x in ("髓系", "单核", "巨噬", "树突", "肥大细胞")):
            targets = ["Inflammatory_unresolved"]
        elif any(x in desc for x in ("成纤维", "内皮", "周细胞", "平滑肌", "脂肪", "施旺", "基质")):
            targets = ["Connective_tissue"]
        for target in targets:
            sets[target].update(genes)
    present = set(gene_names.tolist())
    return {name: sorted(genes & present) for name, genes in sets.items()}, catalog


def marker_indices(gene_names: np.ndarray, marker_sets: dict[str, list[str]]) -> np.ndarray:
    names = {gene for genes in marker_sets.values() for gene in genes}
    return np.asarray([i for i, name in enumerate(gene_names) if name in names], np.int32)


def marker_prior(global_mean: np.ndarray, gene_names: np.ndarray,
                 marker_sets: dict[str, list[str]]) -> np.ndarray:
    prior = np.tile(global_mean[None, :], (N_CLASSES, 1)).astype(np.float64)
    gene_classes: dict[str, set[int]] = {}
    for c, name in enumerate(CLASS_NAMES):
        for gene in marker_sets[name]:
            gene_classes.setdefault(gene, set()).add(c)
    locations: dict[str, list[int]] = {}
    for i, name in enumerate(gene_names):
        locations.setdefault(str(name), []).append(i)
    for gene, expected in gene_classes.items():
        weights = np.full(N_CLASSES, 0.5)
        weights[list(expected)] = 4.0
        weights /= weights.mean()
        for i in locations.get(gene, []):
            prior[:, i] = global_mean[i] * weights
    prior = np.maximum(prior, EPS)
    prior /= prior.sum(1, keepdims=True)
    return prior


def spatial_blocks(positions: pd.DataFrame, block_size: int | None = None) -> np.ndarray:
    if block_size is None:
        block_size = SPATIAL_BLOCK_SIZE_BINS
    row = positions["array_row"].to_numpy(np.int64)
    col = positions["array_col"].to_numpy(np.int64)
    br = (row - row.min()) // block_size
    bc = (col - col.min()) // block_size
    return (br * (int(bc.max(initial=0)) + 1) + bc).astype(np.int64)


def composition_from_alpha(base_mass: np.ndarray, alpha: np.ndarray) -> np.ndarray:
    adjusted = base_mass * alpha[None, :]
    total = adjusted.sum(1)
    return np.divide(
        adjusted, total[:, None], out=np.zeros_like(adjusted), where=total[:, None] > 0
    )


def fit_profiles(matrix: sparse.csc_matrix, composition: np.ndarray,
                 library: np.ndarray, prior: np.ndarray, strength: float,
                 bins: np.ndarray | None = None,
                 genes: np.ndarray | None = None) -> np.ndarray:
    if bins is None:
        bins = np.flatnonzero(composition.sum(1) > 0)
    fit_all_genes = genes is None
    if genes is None:
        genes = np.arange(matrix.shape[0], dtype=np.int32)
    comp = composition[bins].astype(np.float64)
    lib = library[bins].astype(np.float64)
    median = float(np.median(lib)) if len(lib) else 1.0
    weight = np.clip(lib / max(median, 1.0), 0.25, 4.0)
    gram = comp.T @ (comp * weight[:, None])
    multiplier = comp * (weight / np.maximum(lib, 1.0))[:, None]
    cross = np.asarray(matrix[genes, :][:, bins] @ multiplier).T.astype(np.float64)
    positive_diag = np.diag(gram)
    positive_diag = positive_diag[positive_diag > 0]
    scale = float(np.median(positive_diag)) if len(positive_diag) else 1.0
    penalty = max(strength * scale, 1e-8)
    system = gram + np.eye(N_CLASSES) * penalty
    result = np.linalg.solve(system, cross + penalty * prior[:, genes])
    result = np.maximum(result, np.maximum(prior[:, genes] * 1e-4, EPS))
    if fit_all_genes or len(genes) == matrix.shape[0]:
        result /= np.maximum(result.sum(1, keepdims=True), EPS)
    return result


def observed_probability(matrix: sparse.csc_matrix, library: np.ndarray,
                         genes: np.ndarray, bins: np.ndarray) -> np.ndarray:
    return matrix[genes, :][:, bins].T.toarray().astype(np.float64) / library[bins, None]


def select_fit_genes(matrix: sparse.csc_matrix, gene_names: np.ndarray,
                     marker_sets: dict[str, list[str]], limit: int = 512) -> np.ndarray:
    marker = marker_indices(gene_names, marker_sets).tolist()
    total = np.asarray(matrix.sum(1)).ravel()
    ranked = np.argsort(-total, kind="stable")
    chosen = list(marker)
    seen = set(chosen)
    for gene in ranked:
        if len(chosen) >= limit:
            break
        if int(gene) not in seen and total[gene] > 0:
            chosen.append(int(gene))
            seen.add(int(gene))
    return np.asarray(sorted(chosen), np.int32)


def fit_alpha_and_profiles(matrix: sparse.csc_matrix, library: np.ndarray,
                           base_mass: np.ndarray, prior: np.ndarray,
                           fit_genes: np.ndarray, blocks: np.ndarray,
                           strength: float, folds: int, seed: int):
    occupied = np.flatnonzero(base_mass.sum(1) > 0)
    rng = np.random.default_rng(seed)
    optimize_bins = occupied
    if len(optimize_bins) > 20000:
        optimize_bins = np.sort(rng.choice(optimize_bins, 20000, replace=False))
    observed = observed_probability(matrix, library, fit_genes, optimize_bins)
    weight = np.clip(
        library[optimize_bins] / np.median(library[optimize_bins]), 0.25, 4.0
    )
    alpha = np.ones(N_CLASSES, np.float64)
    profiles = fit_profiles(matrix, composition_from_alpha(base_mass, alpha), library, prior, strength)
    for _ in range(3):
        def objective(raw: np.ndarray) -> float:
            centered = raw - raw.mean()
            candidate = np.exp(np.clip(centered, -2.0, 2.0))
            comp = composition_from_alpha(base_mass[optimize_bins], candidate)
            predicted = np.maximum(comp @ profiles[:, fit_genes], EPS)
            loss = np.sum(weight[:, None] * np.square(observed - predicted))
            return float(loss + 0.01 * np.square(centered).sum())
        result = optimize.minimize(
            objective, np.log(alpha), method="L-BFGS-B",
            bounds=[(-2.0, 2.0)] * N_CLASSES,
            options={"maxiter": 60, "ftol": 1e-10},
        )
        raw = result.x - result.x.mean()
        alpha = np.exp(raw)
        profiles = fit_profiles(
            matrix, composition_from_alpha(base_mass, alpha), library,
            prior, strength,
        )

    unique_blocks = np.unique(blocks[occupied])
    n_splits = min(folds, len(unique_blocks))
    alpha_kept = True
    cv_audit = {"n_splits": int(n_splits), "alpha_loss": [], "unit_loss": []}
    if n_splits >= 2:
        cv = GroupKFold(n_splits=n_splits)
        for train_local, test_local in cv.split(occupied, groups=blocks[occupied]):
            train = occupied[train_local]
            test = occupied[test_local]
            for candidate, key in ((alpha, "alpha_loss"), (np.ones(N_CLASSES), "unit_loss")):
                comp = composition_from_alpha(base_mass, candidate)
                fold_profile = fit_profiles(
                    matrix, comp, library, prior, strength,
                    bins=train, genes=fit_genes,
                )
                truth = observed_probability(matrix, library, fit_genes, test)
                prediction = comp[test] @ fold_profile
                cv_audit[key].append(float(np.mean(np.square(truth - prediction))))
        alpha_loss = np.asarray(cv_audit["alpha_loss"])
        unit_loss = np.asarray(cv_audit["unit_loss"])
        best = min(alpha_loss.mean(), unit_loss.mean())
        best_values = alpha_loss if alpha_loss.mean() <= unit_loss.mean() else unit_loss
        one_se = float(best_values.std(ddof=1) / math.sqrt(len(best_values))) if len(best_values) > 1 else 0.0
        if unit_loss.mean() <= best + one_se:
            alpha = np.ones(N_CLASSES)
            alpha_kept = False
            profiles = fit_profiles(
                matrix, composition_from_alpha(base_mass, alpha), library,
                prior, strength,
            )
        cv_audit.update({
            "alpha_mean": float(alpha_loss.mean()),
            "unit_mean": float(unit_loss.mean()),
            "one_se": one_se,
            "alpha_kept": alpha_kept,
        })
    return alpha.astype(np.float32), profiles.astype(np.float32), cv_audit


def parent_profile(profiles: np.ndarray, class_index: int,
                   abundance: np.ndarray, global_profile: np.ndarray) -> np.ndarray:
    members = PARENT_MEMBERS[class_index]
    if not members:
        return global_profile
    if len(members) == 1 and members[0] == class_index:
        return profiles[class_index]
    weights = abundance[list(members)].astype(np.float64)
    if weights.sum() <= 0:
        weights[:] = 1.0
    return np.average(profiles[list(members)], axis=0, weights=weights)


def hierarchy_shrinkage_cv(matrix: sparse.csc_matrix, library: np.ndarray,
                           composition: np.ndarray, prior: np.ndarray,
                           fit_genes: np.ndarray, blocks: np.ndarray,
                           strength: float, folds: int):
    occupied = np.flatnonzero(composition.sum(1) > 0)
    abundance = composition[occupied].sum(0)
    global_total = float(matrix[:, occupied].sum())
    global_fit = np.asarray(matrix[fit_genes, :][:, occupied].sum(1)).ravel().astype(np.float64)
    global_fit = np.maximum(global_fit / max(global_total, EPS), EPS)
    lambdas = np.ones(N_CLASSES, np.float64)
    audit: dict[str, dict] = {}
    n_splits = min(folds, len(np.unique(blocks[occupied])))
    if n_splits < 2:
        return lambdas.astype(np.float32), {"reason": "insufficient spatial blocks"}
    grid = np.asarray([0.0, 0.25, 0.5, 0.75, 1.0])
    losses = np.zeros((N_CLASSES, len(grid), n_splits), np.float64)
    cv = GroupKFold(n_splits=n_splits)
    for fold, (train_local, test_local) in enumerate(cv.split(occupied, groups=blocks[occupied])):
        train = occupied[train_local]
        test = occupied[test_local]
        fold_profiles = fit_profiles(
            matrix, composition, library, prior, strength,
            bins=train, genes=fit_genes,
        )
        truth = observed_probability(matrix, library, fit_genes, test)
        for c in range(N_CLASSES):
            parent = parent_profile(fold_profiles, c, abundance, global_fit)
            for j, value in enumerate(grid):
                candidate = fold_profiles.copy()
                candidate[c] = value * fold_profiles[c] + (1.0 - value) * parent
                prediction = composition[test] @ candidate
                losses[c, j, fold] = np.mean(np.square(truth - prediction))
    for c, name in enumerate(CLASS_NAMES):
        if PARENT_MEMBERS[c] == (c,):
            lambdas[c] = 1.0
            audit[name] = {
                "independent_fraction": 1.0,
                "reason": "class is its own terminal parent; no shrinkage alternative",
                "grid": grid.tolist(),
            }
            continue
        mean = losses[c].mean(1)
        se = losses[c].std(1, ddof=1) / math.sqrt(n_splits) if n_splits > 1 else np.zeros(len(grid))
        best = int(np.argmin(mean))
        eligible = np.flatnonzero(mean <= mean[best] + se[best])
        chosen = int(eligible[0])
        lambdas[c] = grid[chosen]
        audit[name] = {
            "independent_fraction": float(grid[chosen]),
            "cv_mean_loss": mean.tolist(),
            "cv_se": se.tolist(),
            "grid": grid.tolist(),
        }
    return lambdas.astype(np.float32), audit


def apply_hierarchy_shrinkage(profiles: np.ndarray, lambdas: np.ndarray,
                              composition: np.ndarray,
                              global_profile: np.ndarray) -> np.ndarray:
    original = profiles.astype(np.float64)
    abundance = composition.sum(0)
    result = original.copy()
    for c in range(N_CLASSES):
        parent = parent_profile(original, c, abundance, global_profile)
        result[c] = lambdas[c] * original[c] + (1.0 - lambdas[c]) * parent
        result[c] = np.maximum(result[c], EPS)
        result[c] /= result[c].sum()
    return result.astype(np.float32)


def profile_audit(composition: np.ndarray, groups: dict, profiles: np.ndarray,
                  gene_names: np.ndarray, marker_sets: dict[str, list[str]],
                  alpha: np.ndarray, shrinkage: np.ndarray):
    occupied = composition.sum(1) > 0
    design = composition[occupied]
    correlation = np.corrcoef(design.T)
    correlation = np.nan_to_num(correlation, nan=0.0)
    lookup: dict[str, list[int]] = {}
    for i, name in enumerate(gene_names):
        lookup.setdefault(str(name), []).append(i)
    marker_qc = {}
    coverage = {}
    for c, name in enumerate(CLASS_NAMES):
        idx = [i for gene in marker_sets[name] for i in lookup.get(gene, [])]
        if idx:
            other = np.delete(profiles[:, idx], c, axis=0).max(0)
            ratio = profiles[c, idx] / np.maximum(other, EPS)
            marker_ratio = float(np.median(ratio))
        else:
            marker_ratio = None
        gi = np.flatnonzero(groups["class_id"] == c)
        masses = groups["mass"][gi].astype(np.float64)
        effective = float(np.square(masses.sum()) / max(np.square(masses).sum(), EPS))
        coverage[name] = {
            "n_bin_type_groups": int(len(gi)),
            "effective_bin_support": effective,
            "marker_count": int(len(idx)),
            "marker_expected_to_best_other_median_ratio": marker_ratio,
            "alpha": float(alpha[c]),
            "independent_profile_fraction": float(shrinkage[c]),
        }
        marker_qc[name] = marker_ratio
    gram = design.T @ design
    return {
        "design_rank": int(np.linalg.matrix_rank(design)),
        "condition_number": float(np.linalg.cond(gram + np.eye(N_CLASSES) * 1e-12)),
        "composition_correlation": correlation.tolist(),
        "class_coverage": coverage,
        "marker_ratios": marker_qc,
    }


def select_program_genes(matrix: sparse.csc_matrix, library: np.ndarray,
                         occupied: np.ndarray, gene_names: np.ndarray,
                         marker_sets: dict[str, list[str]], n_program: int,
                         n_residual: int, seed: int):
    rng = np.random.default_rng(seed)
    bins = np.flatnonzero(occupied)
    if len(bins) > 12000:
        bins = np.sort(rng.choice(bins, 12000, replace=False))
    mean = np.zeros(matrix.shape[0], np.float64)
    variance = np.zeros(matrix.shape[0], np.float64)
    scale = (10000.0 / library[bins]).astype(np.float32)
    for start in range(0, matrix.shape[0], 384):
        end = min(start + 384, matrix.shape[0])
        block = matrix[start:end, :][:, bins].T.toarray().astype(np.float32)
        block *= scale[:, None]
        np.log1p(block, out=block)
        mean[start:end] = block.mean(0, dtype=np.float64)
        variance[start:end] = block.var(0, dtype=np.float64)
    score = variance / np.maximum(mean, 0.05)
    score[~np.isfinite(score)] = -np.inf
    score[np.asarray(matrix.sum(1)).ravel() <= 0] = -np.inf
    ranked = np.argsort(-score, kind="stable")
    markers = marker_indices(gene_names, marker_sets)

    def build(limit: int) -> np.ndarray:
        result = list(markers.tolist())
        seen = set(result)
        for gene in ranked:
            if len(result) >= limit:
                break
            if int(gene) not in seen and np.isfinite(score[gene]):
                result.append(int(gene))
                seen.add(int(gene))
        return np.asarray(sorted(result), np.int32)

    return build(n_program), build(n_residual), score


def compute_residual_pcs(matrix: sparse.csc_matrix, library: np.ndarray,
                         composition: np.ndarray, profiles: np.ndarray,
                         genes: np.ndarray, occupied: np.ndarray,
                         n_components: int, seed: int):
    observed = matrix[genes, :].T.toarray().astype(np.float32)
    observed /= library[:, None].astype(np.float32)
    expected = composition @ profiles[:, genes]
    residual = np.log1p(10000.0 * observed) - np.log1p(10000.0 * expected)
    indices = np.flatnonzero(occupied)
    rng = np.random.default_rng(seed)
    fit = indices
    if len(fit) > 30000:
        fit = np.sort(rng.choice(fit, 30000, replace=False))
    k = min(n_components, residual.shape[1], len(fit) - 1)
    pca = PCA(n_components=k, svd_solver="randomized", random_state=seed)
    pca.fit(residual[fit])
    scores = np.zeros((matrix.shape[1], k), np.float32)
    for start in range(0, len(indices), 8192):
        part = indices[start:start + 8192]
        scores[part] = pca.transform(residual[part]).astype(np.float32)
    return scores, pca


def build_type_graph(coords: np.ndarray, components: np.ndarray,
                     composition: np.ndarray, residual: np.ndarray,
                     target_fraction: np.ndarray, neighbors: int,
                     radius: float, use_space: bool = True,
                     use_composition: bool = True,
                     use_residual: bool = True):
    n = len(coords)
    if n == 0:
        return sparse.csr_matrix((0, 0)), {}
    k = min(max(1, neighbors + 1), n)
    tree = cKDTree(coords)
    distance, neighbor = tree.query(
        coords, k=k, distance_upper_bound=radius, workers=-1
    )
    if k == 1:
        distance = distance[:, None]
        neighbor = neighbor[:, None]
    valid = np.isfinite(distance) & (neighbor < n)
    rows = np.repeat(np.arange(n), k)[valid.ravel()]
    cols = neighbor.ravel()[valid.ravel()].astype(np.int64)
    ds = distance.ravel()[valid.ravel()]
    same_component = components[rows] == components[cols]
    rows, cols, ds = rows[same_component], cols[same_component], ds[same_component]
    dc = np.linalg.norm(composition[rows] - composition[cols], axis=1)
    de = np.linalg.norm(residual[rows] - residual[cols], axis=1)
    nonself = rows != cols
    spatial_values = ds[nonself & (ds > 0)]
    composition_values = dc[nonself & (dc > 0)]
    residual_values = de[nonself & (de > 0)]
    sigma_space = float(np.median(spatial_values)) if len(spatial_values) else 1.0
    sigma_comp = float(np.median(composition_values)) if len(composition_values) else 1.0
    sigma_residual = float(np.median(residual_values)) if len(residual_values) else 1.0
    weight = np.ones(len(rows), np.float64)
    if use_space:
        weight *= np.exp(-0.5 * np.square(ds / max(sigma_space, 1e-6)))
    if use_composition:
        weight *= np.exp(-0.5 * np.square(dc / max(sigma_comp, 1e-6)))
    if use_residual:
        weight *= np.exp(-0.5 * np.square(de / max(sigma_residual, 1e-6)))
    weight *= np.sqrt(np.maximum(target_fraction[rows] * target_fraction[cols], 0))
    weight[rows == cols] = 1.0
    graph = sparse.coo_matrix((weight, (rows, cols)), shape=(n, n)).tocsr()
    graph = 0.5 * (graph + graph.T)
    graph.setdiag(np.maximum(graph.diagonal(), 1.0))
    graph.eliminate_zeros()
    return graph, {
        "n_nodes": int(n),
        "n_edges_directed": int(graph.nnz),
        "sigma_spatial_bins": sigma_space,
        "sigma_composition": sigma_comp,
        "sigma_residual_pcs": sigma_residual,
        "use_space": bool(use_space),
        "use_composition": bool(use_composition),
        "use_residual": bool(use_residual),
    }


def measured_components(positions: pd.DataFrame) -> np.ndarray:
    row = positions["array_row"].to_numpy(np.int32)
    col = positions["array_col"].to_numpy(np.int32)
    location = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(row, col))}
    edge_r: list[int] = []
    edge_c: list[int] = []
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            for i, (r, c) in enumerate(zip(row, col)):
                j = location.get((int(r + dr), int(c + dc)))
                if j is not None:
                    edge_r.append(i)
                    edge_c.append(j)
    graph = sparse.coo_matrix(
        (np.ones(len(edge_r), np.uint8), (edge_r, edge_c)),
        shape=(len(row), len(row)),
    ).tocsr()
    _, labels = sparse.csgraph.connected_components(graph, directed=False)
    return labels.astype(np.int32)


def type_signal(matrix: sparse.csc_matrix, genes: np.ndarray, bins: np.ndarray,
                target: int, composition: np.ndarray, profiles: np.ndarray,
                library: np.ndarray):
    fraction = composition[bins, target].astype(np.float64)
    observed = matrix[genes, :][:, bins].T.toarray().astype(np.float64)
    observed /= library[bins, None]
    mixture = composition[bins] @ profiles[:, genes]
    signal = observed - (mixture - fraction[:, None] * profiles[target, genes])
    return fraction, signal


def graph_local_logratio(signal: np.ndarray, fraction: np.ndarray,
                         baseline: np.ndarray, graph: sparse.csr_matrix):
    denominator = np.asarray(graph @ np.square(fraction)).ravel()
    positive = denominator[denominator > 0]
    ridge = 0.5 * (float(np.median(positive)) if len(positive) else 1.0)
    numerator = graph @ (fraction[:, None] * signal) + ridge * baseline[None, :]
    local = numerator / (denominator[:, None] + ridge)
    local = np.maximum(local, np.maximum(baseline[None, :] * 1e-3, EPS))
    logratio = np.log(local) - np.log(np.maximum(baseline[None, :], EPS))
    return np.clip(logratio, -5, 5), ridge


def choose_program_count(logratio: np.ndarray, blocks: np.ndarray,
                         max_programs: int, folds: int, seed: int):
    n_splits = min(folds, len(np.unique(blocks)))
    maximum = min(max_programs, logratio.shape[1], len(logratio) - 1)
    if n_splits < 2 or maximum < 1 or len(logratio) < 80:
        return 0, {"reason": "insufficient groups or spatial blocks"}
    candidate = np.arange(maximum + 1)
    fold_loss = np.zeros((len(candidate), n_splits), np.float64)
    cv = GroupKFold(n_splits=n_splits)
    for fold, (train, test) in enumerate(cv.split(logratio, groups=blocks)):
        center = logratio[train].mean(0)
        fold_loss[0, fold] = np.mean(np.square(logratio[test] - center))
        pca = PCA(
            n_components=maximum, svd_solver="randomized", random_state=seed + fold
        ).fit(logratio[train])
        transformed = pca.transform(logratio[test])
        for k in candidate[1:]:
            reconstruction = transformed[:, :k] @ pca.components_[:k] + pca.mean_
            fold_loss[k, fold] = np.mean(np.square(logratio[test] - reconstruction))
    mean = fold_loss.mean(1)
    se = fold_loss.std(1, ddof=1) / math.sqrt(n_splits)
    best = int(np.argmin(mean))
    eligible = np.flatnonzero(mean <= mean[best] + se[best])
    chosen = int(eligible[0])
    return chosen, {
        "candidate_programs": candidate.tolist(),
        "cv_mean_reconstruction_mse": mean.tolist(),
        "cv_se": se.tolist(),
        "one_se_selected_programs": chosen,
        "improvement_over_k0": float(mean[0] - mean[chosen]),
    }


def program_bootstrap_stability(logratio: np.ndarray, components: np.ndarray,
                                k: int, blocks: np.ndarray, replicates: int,
                                seed: int):
    if k == 0 or replicates <= 0:
        return {"replicates": 0, "loading_subspace_cosine_median": None,
                "z_reconstruction_correlation_median": None}
    rng = np.random.default_rng(seed)
    unique = np.unique(blocks)
    cosines: list[float] = []
    correlations: list[float] = []
    full_reconstruction = (logratio @ components[:k].T) @ components[:k]
    for replicate in range(replicates):
        sampled_blocks = rng.choice(unique, len(unique), replace=True)
        rows = np.concatenate([np.flatnonzero(blocks == block) for block in sampled_blocks])
        if len(rows) > 10000:
            rows = rng.choice(rows, 10000, replace=False)
        if len(rows) <= k:
            continue
        pca = PCA(
            n_components=k, svd_solver="randomized", random_state=seed + replicate
        ).fit(logratio[rows])
        singular = np.linalg.svd(components[:k] @ pca.components_.T, compute_uv=False)
        cosines.append(float(np.mean(singular)))
        reconstruction = (logratio @ pca.components_.T) @ pca.components_
        a = full_reconstruction.ravel()
        b = reconstruction.ravel()
        if a.std() > 0 and b.std() > 0:
            correlations.append(float(np.corrcoef(a, b)[0, 1]))
    return {
        "replicates": int(len(cosines)),
        "loading_subspace_cosine_median": float(np.median(cosines)) if cosines else None,
        "loading_subspace_cosine_q05": float(np.quantile(cosines, 0.05)) if cosines else None,
        "z_reconstruction_correlation_median": float(np.median(correlations)) if correlations else None,
        "z_reconstruction_correlation_q05": float(np.quantile(correlations, 0.05)) if correlations else None,
    }


def learn_local_programs(matrix: sparse.csc_matrix, library: np.ndarray,
                         positions: pd.DataFrame, component_id: np.ndarray,
                         composition: np.ndarray, profiles: np.ndarray,
                         residual_scores: np.ndarray,
                         residual_genes: np.ndarray, program_genes: np.ndarray,
                         groups: dict, neighbors: int, radius: float,
                         folds: int, bootstrap_replicates: int, seed: int):
    group_scores = np.zeros((len(groups["bin_index"]), MAX_PROGRAMS), np.float32)
    loadings = np.zeros((N_CLASSES, MAX_PROGRAMS, matrix.shape[0]), np.float32)
    counts = np.zeros(N_CLASSES, np.uint8)
    audits: dict[str, dict] = {}
    coords = positions[["array_row", "array_col"]].to_numpy(np.float64)
    blocks_all = spatial_blocks(positions)
    for c, name in enumerate(CLASS_NAMES):
        group_index = np.flatnonzero(groups["class_id"] == c)
        if len(group_index) < 80:
            audits[name] = {
                "n_bin_type_groups": int(len(group_index)), "n_programs": 0,
                "reason": "fewer than 80 bin-type groups",
            }
            continue
        bins = groups["bin_index"][group_index]
        fraction = composition[bins, c]
        graph, graph_audit = build_type_graph(
            coords[bins], component_id[bins], composition[bins],
            residual_scores[bins], fraction, neighbors, radius,
        )
        fraction_r, signal_r = type_signal(
            matrix, residual_genes, bins, c, composition, profiles, library
        )
        prediction_loss, global_loss = graph_prediction_fold_losses(
            signal_r, fraction_r, profiles[c, residual_genes], graph,
            blocks_all[bins], folds,
        )
        graph_gate = {
            "fold_prediction_mse": prediction_loss.tolist(),
            "fold_global_baseline_mse": global_loss.tolist(),
            "prediction_mse_mean": float(prediction_loss.mean()) if len(prediction_loss) else None,
            "global_baseline_mse_mean": float(global_loss.mean()) if len(global_loss) else None,
            "passed": bool(
                len(prediction_loss)
                and prediction_loss.mean() < global_loss.mean()
            ),
        }
        if not graph_gate["passed"]:
            audits[name] = {
                "n_bin_type_groups": int(len(group_index)), "n_programs": 0,
                "reason": "local graph did not improve held-out prediction over global type baseline",
                "graph": graph_audit, "local_graph_gate": graph_gate,
            }
            continue
        lr, local_ridge = graph_local_logratio(
            signal_r, fraction_r, profiles[c, residual_genes], graph
        )
        k, cv_audit = choose_program_count(
            lr, blocks_all[bins], MAX_PROGRAMS, folds, seed + c * 101
        )
        counts[c] = k
        if k == 0:
            audits[name] = {
                "n_bin_type_groups": int(len(group_index)), "n_programs": 0,
                "graph": graph_audit, "local_graph_gate": graph_gate,
                "program_cv": cv_audit,
            }
            continue
        weight = np.clip(
            fraction * np.sqrt(library[bins] / np.median(library[bins])), 0.05, 4.0
        ).astype(np.float64)
        weighted_mean = np.average(lr, axis=0, weights=weight)
        lr_centered = lr - weighted_mean
        pca = PCA(n_components=k, svd_solver="randomized", random_state=seed + c)
        z = pca.fit_transform(lr_centered).astype(np.float64)
        z -= np.average(z, axis=0, weights=weight)
        group_scores[group_index, :k] = z.astype(np.float32)
        gram = z.T @ (z * weight[:, None])
        ridge = max(0.02 * float(np.median(np.diag(gram))), 1e-8)
        projector = np.linalg.solve(
            gram + np.eye(k) * ridge, z.T * weight[None, :]
        )
        for start in range(0, len(program_genes), 256):
            genes = program_genes[start:start + 256]
            fraction_g, signal_g = type_signal(
                matrix, genes, bins, c, composition, profiles, library
            )
            local_lr, _ = graph_local_logratio(
                signal_g, fraction_g, profiles[c, genes], graph
            )
            local_lr -= np.average(local_lr, axis=0, weights=weight)
            coefficient_block = (projector @ local_lr).astype(np.float32)
            for program in range(k):
                loadings[c, program, genes] = coefficient_block[program]
        for program in range(k):
            norm = float(np.linalg.norm(loadings[c, program, program_genes]))
            if norm <= 1e-12:
                group_scores[group_index, program] = 0
                continue
            loadings[c, program] /= norm
            group_scores[group_index, program] *= norm
            values = loadings[c, program, program_genes]
            pivot = int(np.argmax(np.abs(values)))
            if values[pivot] < 0:
                loadings[c, program] *= -1
                group_scores[group_index, program] *= -1
        stability = program_bootstrap_stability(
            lr_centered, pca.components_, k, blocks_all[bins],
            bootstrap_replicates, seed + 1000 + c * 100,
        )
        audits[name] = {
            "n_bin_type_groups": int(len(group_index)),
            "n_programs": int(k),
            "local_ridge": float(local_ridge),
            "program_gene_count": int(len(program_genes)),
            "residual_gene_count": int(len(residual_genes)),
            "graph": graph_audit,
            "local_graph_gate": graph_gate,
            "program_cv": cv_audit,
            "bootstrap": stability,
        }
        log(f"{name}: selected {k} local programs from {len(group_index):,} bin-type groups")
    return group_scores, loadings, counts, audits


def program_r2(y: np.ndarray, prediction: np.ndarray) -> np.ndarray:
    denominator = np.square(y - y.mean(0)).sum(0)
    return 1.0 - np.square(y - prediction).sum(0) / np.maximum(denominator, EPS)


def ridge_oof(mean_features: np.ndarray, nuisance: np.ndarray, y: np.ndarray,
              blocks: np.ndarray, alpha: float, folds: int):
    y = np.asarray(y, dtype=np.float64)
    if y.ndim == 1:
        y = y[:, None]
    n_splits = min(folds, len(np.unique(blocks)))
    if n_splits < 2:
        return np.tile(y.mean(0), (len(y), 1)), np.tile(y.mean(0), (len(y), 1))
    pred_nuisance = np.zeros_like(y, dtype=np.float64)
    pred_full = np.zeros_like(y, dtype=np.float64)
    cv = GroupKFold(n_splits=n_splits)
    for train, test in cv.split(y, groups=blocks):
        nm = nuisance[train].mean(0)
        ns = nuisance[train].std(0)
        ns[ns < 1e-6] = 1.0
        fm = mean_features[train].mean(0)
        fs = mean_features[train].std(0)
        fs[fs < 1e-6] = 1.0
        nuisance_train = (nuisance[train] - nm) / ns
        nuisance_test = (nuisance[test] - nm) / ns
        mean_train = (mean_features[train] - fm) / fs
        mean_test = (mean_features[test] - fm) / fs
        nuisance_model = Ridge(alpha=alpha).fit(nuisance_train, y[train])
        full_model = Ridge(alpha=alpha).fit(
            np.concatenate([mean_train, nuisance_train], axis=1), y[train]
        )
        nuisance_prediction = np.asarray(
            nuisance_model.predict(nuisance_test), dtype=np.float64
        )
        full_prediction = np.asarray(
            full_model.predict(
                np.concatenate([mean_test, nuisance_test], axis=1)
            ),
            dtype=np.float64,
        )
        if nuisance_prediction.ndim == 1:
            nuisance_prediction = nuisance_prediction[:, None]
        if full_prediction.ndim == 1:
            full_prediction = full_prediction[:, None]
        pred_nuisance[test] = nuisance_prediction
        pred_full[test] = full_prediction
    return pred_nuisance, pred_full


def ridge_coefficient_matrix(model: Ridge) -> np.ndarray:
    """Return sklearn Ridge coefficients as target-by-feature in all cases."""
    coefficient = np.asarray(model.coef_, dtype=np.float64)
    return coefficient[None, :] if coefficient.ndim == 1 else coefficient


def learn_feature_program_mapping(nuclei: dict, capacity: np.ndarray,
                                  cell_features: np.ndarray, groups: dict,
                                  composition: np.ndarray,
                                  residual_scores: np.ndarray,
                                  positions: pd.DataFrame,
                                  group_scores: np.ndarray,
                                  program_count: np.ndarray,
                                  permutations: int, folds: int, seed: int):
    n_features = cell_features.shape[1]
    coefficient = np.zeros((N_CLASSES, n_features, MAX_PROGRAMS), np.float32)
    reliability = np.zeros((N_CLASSES, MAX_PROGRAMS), np.float32)
    link_scores = np.zeros((len(groups["link_cell"]), MAX_PROGRAMS), np.float32)
    audits: dict[str, dict] = {}
    rng = np.random.default_rng(seed)
    block_all = spatial_blocks(positions)
    eta_grid = np.asarray([0.0, 0.25, 0.5, 0.75, 1.0])
    for c, name in enumerate(CLASS_NAMES):
        k = int(program_count[c])
        group_index = np.flatnonzero(groups["class_id"] == c)
        if k == 0 or len(group_index) < 80:
            audits[name] = {
                "n_groups": int(len(group_index)), "n_programs": k,
                "programs": [], "reason": "no validated local programs",
            }
            continue
        bins = groups["bin_index"][group_index]
        mean_all = groups["feature_mean"][group_index].astype(np.float64)
        nuisance_all = np.concatenate([
            groups["feature_variance"][group_index],
            composition[bins], residual_scores[bins],
        ], axis=1).astype(np.float64)
        y_all = group_scores[group_index, :k].astype(np.float64)
        sample_local = np.arange(len(group_index))
        if len(sample_local) > 12000:
            sample_local = np.sort(rng.choice(sample_local, 12000, replace=False))
        mean = mean_all[sample_local]
        nuisance = nuisance_all[sample_local]
        y = y_all[sample_local]
        blocks = block_all[bins[sample_local]]
        best_alpha = 100.0
        best_score = -np.inf
        best_nuisance = best_full = None
        for alpha in (1.0, 10.0, 100.0):
            pred_nuisance, pred_full = ridge_oof(
                mean, nuisance, y, blocks, alpha, folds
            )
            score = float(np.mean(program_r2(y, pred_full)))
            if score > best_score:
                best_score = score
                best_alpha = alpha
                best_nuisance = pred_nuisance
                best_full = pred_full
        assert best_nuisance is not None and best_full is not None
        full_r2 = program_r2(y, best_full)
        nuisance_r2 = program_r2(y, best_nuisance)
        incremental = full_r2 - nuisance_r2

        cells = np.flatnonzero(nuclei["class_id"] == c)
        sampled_groups = group_index[sample_local]
        aggregator = groups["aggregator"][sampled_groups][:, cells]
        local_cell_features = cell_features[cells]
        null_increment = np.zeros((permutations, k), np.float32)
        for permutation in range(permutations):
            permuted = rng.permutation(len(cells))
            permuted_mean = np.asarray(aggregator @ local_cell_features[permuted]).astype(np.float64)
            null_nuisance, null_full = ridge_oof(
                permuted_mean, nuisance, y, blocks, best_alpha, folds
            )
            null_increment[permutation] = (
                program_r2(y, null_full) - program_r2(y, null_nuisance)
            ).astype(np.float32)
        null_q95 = (
            np.quantile(null_increment, 0.95, axis=0)
            if permutations else np.full(k, np.inf)
        )
        effect_prediction = best_full - best_nuisance
        program_audit = []
        for program in range(k):
            losses = [
                float(np.square(
                    y[:, program]
                    - (best_nuisance[:, program] + eta * effect_prediction[:, program])
                ).sum())
                for eta in eta_grid
            ]
            selected_eta = float(eta_grid[int(np.argmin(losses))])
            gate = bool(full_r2[program] > 0 and incremental[program] > null_q95[program])
            if not gate:
                selected_eta = 0.0
            reliability[c, program] = selected_eta
            program_audit.append({
                "program": int(program),
                "heldout_r2_full": float(full_r2[program]),
                "heldout_r2_nuisance": float(nuisance_r2[program]),
                "incremental_r2": float(incremental[program]),
                "permutation_q95_incremental_r2": float(null_q95[program]),
                "permutation_pass": gate,
                "eta": selected_eta,
            })

        mean_center = mean_all.mean(0)
        mean_scale = mean_all.std(0)
        mean_scale[mean_scale < 1e-6] = 1.0
        nuisance_center = nuisance_all.mean(0)
        nuisance_scale = nuisance_all.std(0)
        nuisance_scale[nuisance_scale < 1e-6] = 1.0
        full_design = np.concatenate([
            (mean_all - mean_center) / mean_scale,
            (nuisance_all - nuisance_center) / nuisance_scale,
        ], axis=1)
        final_model = Ridge(alpha=best_alpha).fit(full_design, y_all)
        final_coef = ridge_coefficient_matrix(final_model)
        mean_coef = (
            final_coef[:, :n_features] / mean_scale[None, :]
        ).T
        mean_coef *= reliability[c, :k][None, :]
        coefficient[c, :, :k] = mean_coef.astype(np.float32)

        class_links = np.flatnonzero(nuclei["class_id"][groups["link_cell"]] == c)
        link_group = groups["link_group"][class_links]
        link_cell = groups["link_cell"][class_links]
        deviation = (
            cell_features[link_cell] - groups["feature_mean"][link_group]
        ) @ mean_coef
        scale = np.std(y_all, axis=0)
        deviation = np.clip(
            deviation, -2.0 * np.maximum(scale, 1e-6), 2.0 * np.maximum(scale, 1e-6)
        )
        weights = groups["link_weight"][class_links].astype(np.float64)
        for program in range(k):
            weighted = np.bincount(
                link_group, weights=weights * deviation[:, program],
                minlength=len(groups["bin_index"]),
            )
            deviation[:, program] -= weighted[link_group] / np.maximum(
                groups["mass"][link_group], EPS
            )
        link_scores[class_links, :k] = (
            group_scores[link_group, :k] + deviation
        ).astype(np.float32)
        zero_mean_error = 0.0
        if k:
            check = np.bincount(
                link_group, weights=weights * deviation[:, 0],
                minlength=len(groups["bin_index"]),
            )
            zero_mean_error = float(np.max(np.abs(check[group_index])))
        audits[name] = {
            "n_groups": int(len(group_index)),
            "n_sampled_groups_for_cv": int(len(sample_local)),
            "n_programs": k,
            "spatial_cv_blocks": int(len(np.unique(blocks))),
            "best_ridge_alpha": best_alpha,
            "permutations": int(permutations),
            "max_abs_weighted_zero_mean_error": zero_mean_error,
            "programs": program_audit,
        }
        passed = int(np.count_nonzero(reliability[c, :k]))
        log(f"{name}: morphology gate passed {passed}/{k} programs")
    return link_scores, coefficient, reliability, audits


def global_profile_bootstrap(matrix: sparse.csc_matrix, library: np.ndarray,
                             composition: np.ndarray, prior: np.ndarray,
                             profiles: np.ndarray, genes: np.ndarray,
                             blocks: np.ndarray, strength: float,
                             replicates: int, seed: int):
    if replicates <= 0:
        return {"replicates": 0}
    occupied = np.flatnonzero(composition.sum(1) > 0)
    unique = np.unique(blocks[occupied])
    rng = np.random.default_rng(seed)
    reference = profiles[:, genes].astype(np.float64)
    correlations = np.full((replicates, N_CLASSES), np.nan, np.float64)
    for replicate in range(replicates):
        sampled_blocks = rng.choice(unique, len(unique), replace=True)
        sampled = np.concatenate([
            occupied[blocks[occupied] == block] for block in sampled_blocks
        ])
        if len(sampled) > 20000:
            sampled = rng.choice(sampled, 20000, replace=False)
        estimate = fit_profiles(
            matrix, composition, library, prior, strength,
            bins=sampled, genes=genes,
        )
        for c in range(N_CLASSES):
            if estimate[c].std() > 0 and reference[c].std() > 0:
                correlations[replicate, c] = np.corrcoef(estimate[c], reference[c])[0, 1]
    return {
        "replicates": int(replicates),
        "genes": int(len(genes)),
        "per_class": {
            CLASS_NAMES[c]: {
                "pcc_median": float(np.nanmedian(correlations[:, c])),
                "pcc_q05": float(np.nanquantile(correlations[:, c], 0.05)),
            }
            for c in range(N_CLASSES)
        },
    }


def graph_prediction_mse(signal: np.ndarray, fraction: np.ndarray,
                         baseline: np.ndarray, graph: sparse.csr_matrix,
                         blocks: np.ndarray, folds: int) -> float:
    prediction_loss, _ = graph_prediction_fold_losses(
        signal, fraction, baseline, graph, blocks, folds
    )
    return float(np.mean(prediction_loss)) if len(prediction_loss) else float("nan")


def graph_prediction_fold_losses(signal: np.ndarray, fraction: np.ndarray,
                                 baseline: np.ndarray, graph: sparse.csr_matrix,
                                 blocks: np.ndarray, folds: int):
    n_splits = min(folds, len(np.unique(blocks)))
    if n_splits < 2:
        return np.asarray([]), np.asarray([])
    target = np.maximum(
        signal / np.maximum(fraction[:, None], 1e-6),
        np.maximum(baseline[None, :] * 1e-3, EPS),
    )
    target = np.clip(
        np.log(target) - np.log(np.maximum(baseline[None, :], EPS)), -5, 5
    )
    prediction_loss = []
    baseline_loss = []
    cv = GroupKFold(n_splits=n_splits)
    for train, test in cv.split(target, groups=blocks):
        subgraph = graph[test][:, train]
        denominator = np.asarray(subgraph @ np.square(fraction[train])).ravel()
        positive = denominator[denominator > 0]
        ridge = 0.5 * (float(np.median(positive)) if len(positive) else 1.0)
        numerator = subgraph @ (fraction[train, None] * signal[train])
        local = (numerator + ridge * baseline[None, :]) / (denominator[:, None] + ridge)
        local = np.maximum(local, np.maximum(baseline[None, :] * 1e-3, EPS))
        prediction = np.clip(
            np.log(local) - np.log(np.maximum(baseline[None, :], EPS)), -5, 5
        )
        prediction_loss.append(float(np.mean(np.square(target[test] - prediction))))
        baseline_loss.append(float(np.mean(np.square(target[test]))))
    return np.asarray(prediction_loss), np.asarray(baseline_loss)


def local_program_reliability(local_audit: dict) -> np.ndarray:
    """Convert spatial held-out improvement into a bounded type reliability.

    This reliability is intentionally based only on the already frozen spatial
    cross-validation losses.  It must not be tuned against CellViT/Xenium
    labels, cluster counts, or an embedding appearance.
    """
    reliability = np.zeros(N_CLASSES, np.float32)
    for class_index, class_name in enumerate(CLASS_NAMES):
        gate = local_audit.get(class_name, {}).get("local_graph_gate", {})
        prediction = gate.get("prediction_mse_mean")
        baseline = gate.get("global_baseline_mse_mean")
        if (
            not gate.get("passed", False)
            or prediction is None
            or baseline is None
            or not np.isfinite(prediction)
            or not np.isfinite(baseline)
            or baseline <= 0
        ):
            continue
        reliability[class_index] = np.float32(
            np.clip(1.0 - float(prediction) / float(baseline), 0.0, 1.0)
        )
    return reliability


def interpolate_subspot_offsets(
    nuclei: dict,
    geometry: sparse.csr_matrix,
    capacity: np.ndarray,
    groups: dict,
    positions: pd.DataFrame,
    component_id: np.ndarray,
    composition: np.ndarray,
    residual_scores: np.ndarray,
    group_scores: np.ndarray,
    program_count: np.ndarray,
    local_audit: dict,
    neighbors: int,
    radius: float,
):
    """Interpolate a continuous state at each nucleus position and centre it.

    Neighbouring bins provide context only.  The interpolated state is centred
    with the exact A*u weights inside every bin-by-class group, so it can only
    redistribute that group's mass among its real nuclei.
    """
    n_links = geometry.nnz
    offsets = np.zeros((n_links, MAX_PROGRAMS), np.float32)
    link_cell = groups["link_cell"]
    link_bin = groups["link_bin"]
    link_group = groups["link_group"]
    link_weight = groups["link_weight"].astype(np.float64)
    link_class = nuclei["class_id"][link_cell]
    cell_xy = nuclei["centroid_uv"][:, [0, 1]].astype(np.float64)
    bin_xy = positions[["array_col", "array_row"]].to_numpy(np.float64)
    reliability = local_program_reliability(local_audit)
    lookup = np.full((geometry.shape[0], N_CLASSES), -1, np.int64)
    lookup[groups["bin_index"], groups["class_id"]] = np.arange(
        len(groups["bin_index"]), dtype=np.int64
    )
    audits: dict[str, dict] = {}

    for c, class_name in enumerate(CLASS_NAMES):
        k = int(program_count[c])
        class_links = np.flatnonzero(link_class == c)
        class_groups = np.flatnonzero(groups["class_id"] == c)
        if k == 0 or not len(class_links) or len(class_groups) < 2:
            audits[class_name] = {
                "n_links": int(len(class_links)),
                "n_support_bins": int(len(class_groups)),
                "n_programs": k,
                "reliability": float(reliability[c]),
                "reason": "insufficient validated local state support",
            }
            continue

        support_bins = groups["bin_index"][class_groups]
        tree = cKDTree(bin_xy[support_bins])
        query_k = min(max(1, int(neighbors)), len(support_bins))
        distance, local_index = tree.query(
            cell_xy[link_cell[class_links]],
            k=query_k,
            distance_upper_bound=float(radius),
            workers=-1,
        )
        if query_k == 1:
            distance = distance[:, None]
            local_index = local_index[:, None]
        valid = np.isfinite(distance) & (local_index < len(support_bins))
        safe_index = np.minimum(local_index, len(support_bins) - 1)
        neighbour_bins = support_bins[safe_index]
        owner_bins = link_bin[class_links]
        valid &= component_id[neighbour_bins] == component_id[owner_bins, None]

        positive_distance = distance[valid & (distance > 0)]
        sigma_space = (
            float(np.median(positive_distance))
            if len(positive_distance) else max(float(radius) / 2.0, 1.0)
        )
        comp_delta = composition[neighbour_bins] - composition[owner_bins, None, :]
        comp_d2 = np.square(comp_delta).sum(2)
        residual_delta = residual_scores[neighbour_bins] - residual_scores[owner_bins, None, :]
        residual_d2 = np.square(residual_delta).sum(2)
        comp_positive = comp_d2[valid & (comp_d2 > 0)]
        residual_positive = residual_d2[valid & (residual_d2 > 0)]
        sigma_comp2 = float(np.median(comp_positive)) if len(comp_positive) else 1.0
        sigma_residual2 = (
            float(np.median(residual_positive)) if len(residual_positive) else 1.0
        )
        spatial = np.exp(-0.5 * np.square(distance) / max(sigma_space ** 2, EPS))
        boundary = 0.05 + 0.95 * np.exp(
            -0.5 * comp_d2 / max(sigma_comp2, EPS)
            -0.5 * residual_d2 / max(sigma_residual2, EPS)
        )
        type_support = np.sqrt(np.maximum(composition[neighbour_bins, c], 0.0))
        weight = np.where(valid, spatial * boundary * type_support, 0.0)
        denominator = weight.sum(1)
        neighbour_group = lookup[neighbour_bins, c]
        neighbour_group = np.maximum(neighbour_group, 0)
        field = np.einsum(
            "ln,lnk->lk", weight,
            group_scores[neighbour_group, :k].astype(np.float64),
            optimize=True,
        )
        field = np.divide(
            field, denominator[:, None], out=np.zeros_like(field),
            where=denominator[:, None] > 0,
        )
        missing = denominator <= 0
        if np.any(missing):
            owner_group = lookup[owner_bins[missing], c]
            field[missing] = group_scores[owner_group, :k]

        local_group = link_group[class_links]
        centred = field.copy()
        for axis in range(k):
            weighted_sum = np.bincount(
                local_group,
                weights=link_weight[class_links] * field[:, axis],
                minlength=len(groups["bin_index"]),
            )
            centre = weighted_sum / np.maximum(groups["mass"], EPS)
            centred[:, axis] -= centre[local_group]
            scale = float(np.std(group_scores[class_groups, axis]))
            if scale > 1e-8:
                centred[:, axis] = np.clip(centred[:, axis], -2.0 * scale, 2.0 * scale)
            # Clipping can perturb the exact centre; restore it once.
            correction = np.bincount(
                local_group,
                weights=link_weight[class_links] * centred[:, axis],
                minlength=len(groups["bin_index"]),
            ) / np.maximum(groups["mass"], EPS)
            centred[:, axis] -= correction[local_group]
        centred *= float(reliability[c])
        offsets[class_links, :k] = centred.astype(np.float32)

        check = np.zeros(k, np.float64)
        for axis in range(k):
            sums = np.bincount(
                local_group,
                weights=link_weight[class_links] * centred[:, axis],
                minlength=len(groups["bin_index"]),
            )
            check[axis] = np.max(np.abs(sums[class_groups]), initial=0.0)
        audits[class_name] = {
            "n_links": int(len(class_links)),
            "n_support_bins": int(len(class_groups)),
            "n_programs": k,
            "neighbors": int(query_k),
            "radius": float(radius),
            "sigma_space": sigma_space,
            "sigma_composition_squared": sigma_comp2,
            "sigma_residual_squared": sigma_residual2,
            "reliability": float(reliability[c]),
            "missing_interpolation_links": int(np.count_nonzero(missing)),
            "max_abs_weighted_zero_mean_error": float(check.max(initial=0.0)),
        }
    return offsets, audits


def aggregate_link_states_to_cells(
    geometry: sparse.csr_matrix,
    capacity: np.ndarray,
    *states: np.ndarray,
) -> list[np.ndarray]:
    link_cell = geometry.indices
    weight = geometry.data.astype(np.float64) * capacity[link_cell].astype(np.float64)
    denominator = np.bincount(
        link_cell, weights=weight, minlength=geometry.shape[1]
    ).astype(np.float64)
    outputs = []
    for state in states:
        numerator = np.zeros((geometry.shape[1], state.shape[1]), np.float64)
        np.add.at(numerator, link_cell, weight[:, None] * state.astype(np.float64))
        outputs.append(np.divide(
            numerator, denominator[:, None], out=np.zeros_like(numerator),
            where=denominator[:, None] > 0,
        ).astype(np.float32))
    return outputs


def _row_correlation(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left = left - left.mean(1, keepdims=True)
    right = right - right.mean(1, keepdims=True)
    denominator = np.sqrt(np.square(left).sum(1) * np.square(right).sum(1))
    return np.divide(
        (left * right).sum(1), denominator,
        out=np.full(len(left), np.nan, np.float64), where=denominator > 1e-12,
    )


def calibrate_offset_scales(
    reference_path: Path | None,
    max_cells: int,
    nuclei: dict,
    geometry: sparse.csr_matrix,
    capacity: np.ndarray,
    positions: pd.DataFrame,
    gene_names: np.ndarray,
    profiles: np.ndarray,
    loadings: np.ndarray,
    program_count: np.ndarray,
    group_link_scores: np.ndarray,
    spatial_offsets: np.ndarray,
    morphology_offsets: np.ndarray,
    default_spatial: float,
    default_morphology: float,
    seed: int,
):
    if reference_path is None:
        return float(default_spatial), float(default_morphology), {
            "calibrated": False,
            "reason": "no high-resolution calibration reference supplied",
            "spatial_scale": float(default_spatial),
            "morphology_scale": float(default_morphology),
        }
    if not reference_path.exists():
        raise FileNotFoundError(reference_path)

    import anndata as ad

    reference = ad.read_h5ad(reference_path)
    if "cell_id" not in reference.obs:
        raise ValueError("Calibration H5AD lacks obs['cell_id']")
    current_lookup = {
        int(cell_id): index for index, cell_id in enumerate(nuclei["cell_id"])
    }
    pairs = [
        (ref_index, current_lookup.get(int(cell_id), -1))
        for ref_index, cell_id in enumerate(reference.obs["cell_id"].to_numpy())
    ]
    pairs = np.asarray([pair for pair in pairs if pair[1] >= 0], np.int64)
    current_gene_lookup = {}
    for index, gene in enumerate(text_array(gene_names)):
        current_gene_lookup.setdefault(str(gene), index)
    common = [
        (ref_index, current_gene_lookup[str(gene)])
        for ref_index, gene in enumerate(reference.var_names.astype(str))
        if str(gene) in current_gene_lookup
    ]
    if len(common) < 30 or len(pairs) < 200:
        return float(default_spatial), float(default_morphology), {
            "calibrated": False,
            "reason": "insufficient matched cells or common genes",
            "matched_cells": int(len(pairs)),
            "common_genes": int(len(common)),
            "spatial_scale": float(default_spatial),
            "morphology_scale": float(default_morphology),
        }
    ref_gene = np.asarray([value[0] for value in common], np.int64)
    model_gene = np.asarray([value[1] for value in common], np.int64)
    ref_rows = pairs[:, 0]
    current_rows = pairs[:, 1]
    ref_counts = reference.X[ref_rows][:, ref_gene]
    if sparse.issparse(ref_counts):
        ref_total = np.asarray(ref_counts.sum(1)).ravel()
    else:
        ref_total = np.asarray(ref_counts).sum(1)
    keep = ref_total >= 3
    ref_rows = ref_rows[keep]
    current_rows = current_rows[keep]
    ref_total = ref_total[keep]
    if len(ref_rows) > max_cells:
        rng = np.random.default_rng(seed)
        chosen = np.sort(rng.choice(len(ref_rows), max_cells, replace=False))
        ref_rows = ref_rows[chosen]
        current_rows = current_rows[chosen]
    reference_values = reference.X[ref_rows][:, ref_gene]
    if sparse.issparse(reference_values):
        reference_values = reference_values.toarray()
    reference_values = np.asarray(reference_values, np.float64)
    row_total = reference_values.sum(1)
    reference_log = np.log1p(
        reference_values * np.divide(
            10000.0, row_total, out=np.zeros_like(row_total), where=row_total > 0
        )[:, None]
    )

    group_cell, spatial_cell, morphology_cell = aggregate_link_states_to_cells(
        geometry, capacity, group_link_scores, spatial_offsets, morphology_offsets
    )
    geometry_csc = geometry.tocsc()
    primary_bin = np.zeros(len(current_rows), np.int32)
    for out_index, cell in enumerate(current_rows):
        start, end = int(geometry_csc.indptr[cell]), int(geometry_csc.indptr[cell + 1])
        bins = geometry_csc.indices[start:end]
        values = geometry_csc.data[start:end]
        primary_bin[out_index] = int(bins[int(np.argmax(values))]) if len(bins) else 0
    blocks = spatial_blocks(positions)[primary_bin]
    unique_blocks = np.unique(blocks)
    holdout_blocks = unique_blocks[np.arange(len(unique_blocks)) % 5 == 0]
    test_mask = np.isin(blocks, holdout_blocks)
    train_mask = ~test_mask
    class_id = nuclei["class_id"][current_rows]
    group_key = primary_bin.astype(np.int64) * N_CLASSES + class_id.astype(np.int64)

    def centre_by_group(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        _, inverse = np.unique(group_key, return_inverse=True)
        count = np.bincount(inverse).astype(np.float64)
        total = np.zeros((len(count), values.shape[1]), np.float64)
        np.add.at(total, inverse, values)
        centred = values - total[inverse] / count[inverse, None]
        return centred, count[inverse] >= 2

    reference_residual, replicated = centre_by_group(reference_log)
    rng = np.random.default_rng(seed + 91)
    permuted_reference = reference_residual.copy()
    for key in np.unique(group_key):
        members = np.flatnonzero(group_key == key)
        if len(members) > 1:
            permuted_reference[members] = reference_residual[np.roll(members, 1)]

    def prediction(spatial_scale: float, morphology_scale: float) -> np.ndarray:
        state = (
            group_cell[current_rows].astype(np.float64)
            + spatial_scale * spatial_cell[current_rows].astype(np.float64)
            + morphology_scale * morphology_cell[current_rows].astype(np.float64)
        )
        result = np.empty((len(current_rows), len(model_gene)), np.float64)
        for c in np.unique(class_id):
            selected = class_id == c
            k = int(program_count[c])
            logits = np.log(np.maximum(profiles[c, model_gene], EPS))[None, :]
            if k:
                logits = logits + state[selected, :k] @ loadings[c, :k][:, model_gene]
            logits -= logits.max(1, keepdims=True)
            probability = np.exp(np.clip(logits, -30, 0))
            probability /= np.maximum(probability.sum(1, keepdims=True), EPS)
            result[selected] = np.log1p(10000.0 * probability)
        return result

    def metrics(predicted: np.ndarray, mask: np.ndarray) -> dict:
        pred_residual, _ = centre_by_group(predicted)
        valid = mask & replicated
        matched = _row_correlation(pred_residual[valid], reference_residual[valid])
        random_corr = _row_correlation(pred_residual[valid], permuted_reference[valid])
        total_corr = _row_correlation(predicted[mask], reference_log[mask])
        pred_var = np.var(pred_residual[valid], axis=0)
        ref_var = np.var(reference_residual[valid], axis=0)
        variance_corr = (
            float(np.corrcoef(pred_var, ref_var)[0, 1])
            if np.std(pred_var) > 0 and np.std(ref_var) > 0 else 0.0
        )
        ratio = float(np.median(
            (pred_var + 1e-8) / np.maximum(ref_var + 1e-8, 1e-8)
        ))
        matched_median = float(np.nanmedian(matched)) if len(matched) else -1.0
        random_median = float(np.nanmedian(random_corr)) if len(random_corr) else 0.0
        return {
            "n_cells": int(np.count_nonzero(mask)),
            "n_within_group_cells": int(np.count_nonzero(valid)),
            "within_group_matched_median": matched_median,
            "within_group_random_median": random_median,
            "within_group_advantage": matched_median - random_median,
            "total_profile_median": float(np.nanmedian(total_corr)) if len(total_corr) else -1.0,
            "within_group_gene_variance_correlation": variance_corr,
            "within_group_variance_ratio_median": ratio,
        }

    candidates = []
    best = None
    for spatial_scale in (0.0, 0.5, 1.0, 2.0, 4.0, 8.0):
        for morphology_scale in (0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
            predicted = prediction(spatial_scale, morphology_scale)
            train = metrics(predicted, train_mask)
            variance_ratio = max(
                train["within_group_variance_ratio_median"], 1e-8
            )
            if variance_ratio < 0.10:
                variance_penalty = abs(math.log(variance_ratio / 0.10))
            elif variance_ratio > 1.0:
                variance_penalty = abs(math.log(variance_ratio))
            else:
                variance_penalty = 0.0
            objective = (
                train["within_group_advantage"]
                + 0.10 * train["total_profile_median"]
                + 0.15 * train["within_group_gene_variance_correlation"]
                - 0.08 * variance_penalty
            )
            candidate = {
                "spatial_scale": spatial_scale,
                "morphology_scale": morphology_scale,
                "training_objective": float(objective),
                "training": train,
            }
            candidates.append(candidate)
            if best is None or objective > best[0]:
                best = (objective, spatial_scale, morphology_scale, predicted, train)
    assert best is not None
    test = metrics(best[3], test_mask)
    gate_passed = bool(
        best[4]["within_group_advantage"] > 0
        and test["within_group_advantage"] > 0
        and test["within_group_gene_variance_correlation"] > 0
    )
    final_spatial = float(best[1]) if gate_passed else 0.0
    final_morphology = float(best[2]) if gate_passed else 0.0
    return final_spatial, final_morphology, {
        "calibrated": True,
        "reference": str(reference_path),
        "matched_cells_after_count_gate": int(len(current_rows)),
        "common_genes": int(len(common)),
        "train_cells": int(np.count_nonzero(train_mask)),
        "heldout_cells": int(np.count_nonzero(test_mask)),
        "pre_gate_spatial_scale": float(best[1]),
        "pre_gate_morphology_scale": float(best[2]),
        "heldout_cell_specific_gate_passed": gate_passed,
        "selected_spatial_scale": final_spatial,
        "selected_morphology_scale": final_morphology,
        "training": best[4],
        "heldout": test,
        "candidates": candidates,
    }


def run_ablations(matrix: sparse.csc_matrix, library: np.ndarray,
                  positions: pd.DataFrame, component_id: np.ndarray,
                  composition: np.ndarray, profiles: np.ndarray,
                  residual_scores: np.ndarray, residual_genes: np.ndarray,
                  groups: dict, marker_audit: dict, no_marker_audit: dict,
                  feature_audit: dict,
                  neighbors: int, radius: float, folds: int, seed: int):
    rng = np.random.default_rng(seed)
    coords = positions[["array_row", "array_col"]].to_numpy(np.float64)
    blocks = spatial_blocks(positions)
    genes = residual_genes[:min(128, len(residual_genes))]
    variants = {
        "full": (True, True, True),
        "without_spatial": (False, True, True),
        "without_composition": (True, False, True),
        "without_residual": (True, True, False),
    }
    graph_scores = {key: [] for key in variants}
    for c in range(N_CLASSES):
        gi = np.flatnonzero(groups["class_id"] == c)
        if len(gi) < 80:
            continue
        if len(gi) > 4000:
            gi = np.sort(rng.choice(gi, 4000, replace=False))
        bins = groups["bin_index"][gi]
        fraction, signal = type_signal(
            matrix, genes, bins, c, composition, profiles, library
        )
        for key, flags in variants.items():
            graph, _ = build_type_graph(
                coords[bins], component_id[bins], composition[bins],
                residual_scores[bins], fraction, neighbors, radius,
                use_space=flags[0], use_composition=flags[1], use_residual=flags[2],
            )
            score = graph_prediction_mse(
                signal, fraction, profiles[c, genes], graph, blocks[bins], folds
            )
            if np.isfinite(score):
                graph_scores[key].append(score)
    morphology_observed = []
    morphology_null = []
    for audit in feature_audit.values():
        for program in audit.get("programs", []):
            morphology_observed.append(program["incremental_r2"])
            morphology_null.append(program["permutation_q95_incremental_r2"])
    class_marker = [
        value for value in marker_audit.get("marker_ratios", {}).values()
        if value is not None
    ]
    class_marker_without = [
        value for value in no_marker_audit.get("marker_ratios", {}).values()
        if value is not None
    ]
    global_mean = np.asarray(matrix.sum(1)).ravel().astype(np.float64)
    global_mean /= max(global_mean.sum(), EPS)
    occupied = np.flatnonzero(composition.sum(1) > 0)
    eval_bins = occupied
    if len(eval_bins) > 10000:
        eval_bins = np.sort(rng.choice(eval_bins, 10000, replace=False))
    eval_genes = genes
    truth = observed_probability(matrix, library, eval_genes, eval_bins)
    global_prediction = composition[eval_bins] @ profiles[:, eval_genes]
    capacity_only_prediction = np.tile(global_mean[eval_genes], (len(eval_bins), 1))
    return [
        {"ablation": "without_marker_prior", "metric": "median_marker_enrichment",
         "value": float(np.median(class_marker_without)) if class_marker_without else None,
         "main": float(np.median(class_marker)) if class_marker else None},
        {"ablation": "without_spatial_kernel", "metric": "heldout_logratio_mse",
         "value": float(np.mean(graph_scores["without_spatial"])) if graph_scores["without_spatial"] else None,
         "main": float(np.mean(graph_scores["full"])) if graph_scores["full"] else None},
        {"ablation": "without_composition_kernel", "metric": "heldout_logratio_mse",
         "value": float(np.mean(graph_scores["without_composition"])) if graph_scores["without_composition"] else None,
         "main": float(np.mean(graph_scores["full"])) if graph_scores["full"] else None},
        {"ablation": "without_residual_kernel", "metric": "heldout_logratio_mse",
         "value": float(np.mean(graph_scores["without_residual"])) if graph_scores["without_residual"] else None,
         "main": float(np.mean(graph_scores["full"])) if graph_scores["full"] else None},
        {"ablation": "within_class_nucleus_feature_permutation",
         "metric": "median_incremental_r2_vs_null_q95",
         "value": float(np.median(morphology_null)) if morphology_null else None,
         "main": float(np.median(morphology_observed)) if morphology_observed else None},
        {"ablation": "shared_bin_type_state", "metric": "morphology_incremental_r2_removed",
         "value": 0.0, "main": float(np.median(morphology_observed)) if morphology_observed else None},
        {"ablation": "contour_capacity_only", "metric": "bin_expression_mse",
         "value": float(np.mean(np.square(truth - capacity_only_prediction))),
         "main": float(np.mean(np.square(truth - global_prediction)))},
    ]


def allocate_bin(matrix: sparse.csc_matrix, bin_index: int,
                 geometry: sparse.csr_matrix, capacity: np.ndarray,
                 class_id: np.ndarray, link_scores: np.ndarray,
                 profiles: np.ndarray, loadings: np.ndarray,
                 program_count: np.ndarray,
                 within_type_temperature: float = 1.0):
    cell_start, cell_end = int(geometry.indptr[bin_index]), int(geometry.indptr[bin_index + 1])
    matrix_start, matrix_end = int(matrix.indptr[bin_index]), int(matrix.indptr[bin_index + 1])
    cells = geometry.indices[cell_start:cell_end]
    overlap = geometry.data[cell_start:cell_end].astype(np.float64)
    genes = matrix.indices[matrix_start:matrix_end].astype(np.int32, copy=False)
    counts = matrix.data[matrix_start:matrix_end].astype(np.float64)
    if not len(cells) or not len(genes):
        return cells, genes, counts, np.empty((len(cells), len(genes)), np.float32)
    local_class = class_id[cells]
    local_scores = link_scores[cell_start:cell_end]
    link_capacity = overlap * capacity[cells].astype(np.float64)
    classes = np.unique(local_class)
    type_propensity = np.zeros((len(classes), len(genes)), np.float64)
    conditional = np.zeros((len(cells), len(genes)), np.float64)
    for class_row, c in enumerate(classes):
        selected = local_class == c
        k = int(program_count[c])
        class_capacity = link_capacity[selected]
        class_mass = float(class_capacity.sum())
        if class_mass <= 0:
            continue
        mean_score = np.average(
            local_scores[selected, :k].astype(np.float64),
            axis=0, weights=class_capacity,
        ) if k else np.zeros(0, np.float64)
        type_modifier = np.zeros(len(genes), np.float64)
        cell_modifier = np.zeros((int(selected.sum()), len(genes)), np.float64)
        if k:
            type_modifier = (
                mean_score[None, :]
                @ loadings[c, :k][:, genes].astype(np.float64)
            ).ravel()
            centred = local_scores[selected, :k].astype(np.float64) - mean_score[None, :]
            cell_modifier = (
                centred @ loadings[c, :k][:, genes].astype(np.float64)
            ) / max(float(within_type_temperature), EPS)
        type_propensity[class_row] = (
            class_mass * np.maximum(profiles[c, genes], EPS)
            * np.exp(np.clip(type_modifier, -5, 5))
        )
        within = class_capacity[:, None] * np.exp(np.clip(cell_modifier, -5, 5))
        within_denominator = within.sum(0)
        missing_within = within_denominator <= 0
        if np.any(missing_within):
            within[:, missing_within] = class_capacity[:, None]
            within_denominator[missing_within] = class_mass
        conditional[selected] = within / np.maximum(within_denominator[None, :], EPS)

    type_denominator = type_propensity.sum(0)
    missing_type = type_denominator <= 0
    if np.any(missing_type):
        for class_row, c in enumerate(classes):
            selected = local_class == c
            type_propensity[class_row, missing_type] = link_capacity[selected].sum()
        type_denominator[missing_type] = max(link_capacity.sum(), EPS)
    type_probability = type_propensity / np.maximum(type_denominator[None, :], EPS)
    responsibility = np.zeros_like(conditional)
    for class_row, c in enumerate(classes):
        selected = local_class == c
        responsibility[selected] = conditional[selected] * type_probability[class_row][None, :]
    # Correct the float32 representation so per-gene sums remain as close to one
    # as possible after storage.
    q = responsibility.astype(np.float32)
    if len(cells):
        pivot = np.argmax(responsibility, axis=0)
        correction = 1.0 - q.sum(0, dtype=np.float64)
        q[pivot, np.arange(len(genes))] = (
            q[pivot, np.arange(len(genes))].astype(np.float64) + correction
        ).astype(np.float32)
    return cells, genes, counts, q


def write_text_dataset(group: h5py.Group, name: str, values) -> h5py.Dataset:
    dataset = group.create_dataset(
        name, data=np.asarray(values, dtype=object), dtype=h5py.string_dtype("utf-8")
    )
    dataset.attrs["encoding-type"] = "string-array"
    dataset.attrs["encoding-version"] = "0.2.0"
    return dataset


def encoded_array(group: h5py.Group, name: str, values, **kwargs) -> h5py.Dataset:
    dataset = group.create_dataset(name, data=values, **kwargs)
    dataset.attrs["encoding-type"] = "array"
    dataset.attrs["encoding-version"] = "0.2.0"
    return dataset


def encoded_scalar_string(group: h5py.Group, name: str, value: str) -> h5py.Dataset:
    dataset = group.create_dataset(name, data=np.asarray(value, dtype=h5py.string_dtype("utf-8")))
    dataset.attrs["encoding-type"] = "string"
    dataset.attrs["encoding-version"] = "0.2.0"
    return dataset


def init_dict_group(parent: h5py.Group, name: str) -> h5py.Group:
    group = parent.create_group(name)
    group.attrs["encoding-type"] = "dict"
    group.attrs["encoding-version"] = "0.1.0"
    return group


def init_dataframe(parent: h5py.Group, name: str, index: np.ndarray,
                   columns: list[str]) -> h5py.Group:
    group = parent.create_group(name)
    group.attrs["encoding-type"] = "dataframe"
    group.attrs["encoding-version"] = "0.2.0"
    group.attrs["_index"] = "_index"
    group.attrs["column-order"] = np.asarray(columns, dtype=h5py.string_dtype("utf-8"))
    write_text_dataset(group, "_index", index)
    return group


def write_categorical(group: h5py.Group, name: str, codes: np.ndarray,
                      categories: list[str]) -> None:
    cat = group.create_group(name)
    cat.attrs["encoding-type"] = "categorical"
    cat.attrs["encoding-version"] = "0.2.0"
    cat.attrs["ordered"] = False
    encoded_array(cat, "codes", codes.astype(np.int8), compression="lzf")
    write_text_dataset(cat, "categories", categories)


def append_dataset(dataset: h5py.Dataset, values: np.ndarray) -> None:
    start = dataset.shape[0]
    dataset.resize((start + len(values),))
    dataset[start:] = values


def output_resolution_by_class(shrinkage: np.ndarray, program_count: np.ndarray,
                               reliability: np.ndarray):
    names = [
        "capacity_only", "class_global", "multiview_type_state",
        "morphology_refined",
    ]
    code = np.zeros(N_CLASSES, np.int8)
    for c in range(N_CLASSES):
        if shrinkage[c] > 0:
            code[c] = 1
        if program_count[c] > 0:
            code[c] = 2
        if np.any(reliability[c, :int(program_count[c])] > 0):
            code[c] = 3
    return code, names


def write_allocation_file(path: Path, args: argparse.Namespace,
                          paths: dict[str, Path], matrix: sparse.csc_matrix,
                          library: np.ndarray, barcodes: np.ndarray,
                          source_bin_index: np.ndarray, gene_id: np.ndarray,
                          gene_name: np.ndarray, geometry: sparse.csr_matrix,
                          nuclei: dict, capacity: np.ndarray,
                          link_scores: np.ndarray, profiles: np.ndarray,
                          loadings: np.ndarray, program_count: np.ndarray,
                          summary: dict):
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    n_links = np.diff(geometry.indptr).astype(np.int64)
    bin_nnz = np.diff(matrix.indptr).astype(np.int64)
    record_count_by_bin = n_links * bin_nnz
    record_indptr = np.zeros(matrix.shape[1] + 1, np.int64)
    np.cumsum(record_count_by_bin, out=record_indptr[1:])
    total_records = int(record_indptr[-1])
    log(
        f"Writing allocation provenance: {matrix.nnz:,} observed bin-gene values, "
        f"{total_records:,} cell responsibilities"
    )
    occupied = n_links > 0
    unassigned = ~occupied
    unassigned_gene = np.asarray(matrix @ unassigned.astype(np.float64)).ravel()
    occupied_gene = np.asarray(matrix @ occupied.astype(np.float64)).ravel()
    assigned_gene = np.zeros(matrix.shape[0], np.float64)
    cell_total = np.zeros(geometry.shape[1], np.float64)
    program_numerator = np.zeros((geometry.shape[1], MAX_PROGRAMS), np.float64)
    max_bin_gene_error = 0.0
    with h5py.File(path, "w") as h:
        h.attrs.update({
            "schema": f"visium_hd_{BIN_SIZE_UM:03d}um_spot2cell_mv2_allocation_v1",
            "complete": 0,
            "sample": args.sample,
            "expression_units": "fractional_captured_counts",
            "no_spatial_extrapolation": True,
            "occupied_bin_background_component": False,
            "geometry": "nucleus_contour_area_overlap",
            "allocation_model": "two_stage_spot_gene_to_type_then_KL_regularized_within_type",
            "within_type_temperature": float(
                getattr(args, "within_type_temperature", 1.0)
            ),
            "matrix_source": str(paths["matrix"]),
            "cell_feature_source": str(paths["cell_features"]),
            "summary_json": json.dumps(summary, ensure_ascii=False),
        })
        observed = h.create_group("observed")
        observed.attrs["format"] = "gene_by_bin_csc"
        observed.create_dataset("shape", data=np.asarray(matrix.shape, np.int64))
        observed.create_dataset("bin_indptr", data=matrix.indptr, compression="lzf")
        observed.create_dataset("gene_index", data=matrix.indices, compression="lzf")
        observed.create_dataset("count", data=matrix.data, compression="lzf")
        observed.create_dataset("positive_source_bin_index", data=source_bin_index, compression="lzf")
        write_text_dataset(observed, "barcode", barcodes)
        write_text_dataset(observed, "gene_id", gene_id)
        write_text_dataset(observed, "gene_name", gene_name)

        gg = h.create_group("geometry")
        gg.attrs["format"] = "bin_by_cell_csr"
        gg.create_dataset("shape", data=np.asarray(geometry.shape, np.int64))
        gg.create_dataset("indptr", data=geometry.indptr, compression="lzf")
        gg.create_dataset("cell_index", data=geometry.indices, compression="lzf")
        gg.create_dataset("overlap_fraction", data=geometry.data, compression="lzf")
        gg.create_dataset("geometry_fallback", data=nuclei["geometry_fallback"], compression="lzf")

        allocation = h.create_group("allocation")
        allocation.attrs["layout"] = (
            "For bin b, reshape responsibility[bin_record_indptr[b]:bin_record_indptr[b+1]] "
            "as (observed genes in bin, geometry cells in bin)."
        )
        allocation.create_dataset("bin_record_indptr", data=record_indptr, compression="lzf")
        chunk = (min(1_000_000, max(total_records, 1)),)
        responsibility_ds = allocation.create_dataset(
            "responsibility", shape=(total_records,), dtype=np.float32,
            chunks=chunk, compression="lzf",
        )
        last = time.time()
        for b in np.flatnonzero(occupied):
            cells, genes, counts, responsibility = allocate_bin(
                matrix, int(b), geometry, capacity, nuclei["class_id"],
                link_scores, profiles, loadings, program_count,
                getattr(args, "within_type_temperature", 1.0),
            )
            start, end = int(record_indptr[b]), int(record_indptr[b + 1])
            responsibility_ds[start:end] = responsibility.T.ravel()
            allocation_mass = responsibility.astype(np.float64) * counts[None, :]
            reconstructed = allocation_mass.sum(0)
            if len(reconstructed):
                max_bin_gene_error = max(
                    max_bin_gene_error, float(np.max(np.abs(reconstructed - counts)))
                )
            np.add.at(assigned_gene, genes, reconstructed)
            link_total = allocation_mass.sum(1)
            np.add.at(cell_total, cells, link_total)
            program_numerator[cells] += link_total[:, None] * link_scores[
                geometry.indptr[b]:geometry.indptr[b + 1]
            ]
            if time.time() - last > 60:
                log(f"Allocated {b + 1:,}/{matrix.shape[1]:,} bins")
                last = time.time()
        audit = h.create_group("audit")
        audit.create_dataset("observed_gene_counts_occupied_bins", data=occupied_gene, compression="lzf")
        audit.create_dataset("assigned_gene_counts", data=assigned_gene, compression="lzf")
        audit.create_dataset("unassigned_gene_counts_zero_cell_bins", data=unassigned_gene, compression="lzf")
        audit.create_dataset("unassigned_bin_index", data=np.flatnonzero(unassigned), compression="lzf")
        audit.attrs.update({
            "max_abs_gene_conservation_error": float(np.max(np.abs(assigned_gene - occupied_gene))),
            "max_abs_bin_gene_float32_error": max_bin_gene_error,
            "observed_total_counts": float(library.sum()),
            "observed_counts_occupied_bins": float(library[occupied].sum()),
            "unassigned_counts_zero_cell_bins": float(library[unassigned].sum()),
            "assigned_cell_counts": float(cell_total.sum()),
        })
        h.attrs["complete"] = 1
        h.flush()
    cell_program = np.divide(
        program_numerator, cell_total[:, None],
        out=np.zeros_like(program_numerator), where=cell_total[:, None] > 0,
    ).astype(np.float32)
    audit_summary = {
        "observed_total_counts": float(library.sum()),
        "observed_counts_occupied_bins": float(library[occupied].sum()),
        "unassigned_counts_zero_cell_bins": float(library[unassigned].sum()),
        "assigned_cell_counts": float(cell_total.sum()),
        "max_abs_gene_error": float(np.max(np.abs(assigned_gene - occupied_gene))),
        "max_abs_bin_gene_float32_error": max_bin_gene_error,
        "allocation_records": total_records,
    }
    return cell_total, cell_program, audit_summary


def marker_class_strings(gene_names: np.ndarray,
                         marker_sets: dict[str, list[str]]) -> np.ndarray:
    mapping: dict[str, list[str]] = {}
    for class_name, genes in marker_sets.items():
        for gene in genes:
            mapping.setdefault(gene, []).append(class_name)
    return np.asarray([",".join(mapping.get(str(gene), [])) for gene in gene_names], dtype="U")


def write_h5ad(path: Path, allocation_path: Path, args: argparse.Namespace,
               matrix: sparse.csc_matrix, barcodes: np.ndarray,
               gene_id: np.ndarray, gene_name: np.ndarray,
               feature_type: np.ndarray, positions: pd.DataFrame,
               geometry: sparse.csr_matrix, nuclei: dict,
               capacity: np.ndarray, cell_features: np.ndarray,
               embedding_scores: np.ndarray, cell_total: np.ndarray,
               cell_program: np.ndarray, profiles: np.ndarray,
               loadings: np.ndarray, program_count: np.ndarray,
               reliability: np.ndarray, shrinkage: np.ndarray,
               marker_sets: dict[str, list[str]], program_genes: np.ndarray,
               residual_genes: np.ndarray, groups: dict,
               group_scores: np.ndarray, pca_models: list[PCA | None],
               residual_pca: PCA, summary: dict):
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    n_cells, n_genes = geometry.shape[1], matrix.shape[0]
    resolution_by_class, resolution_names = output_resolution_by_class(
        shrinkage, program_count, reliability
    )
    with h5py.File(path, "w") as h, h5py.File(allocation_path, "r") as allocation:
        h.attrs["encoding-type"] = "anndata"
        h.attrs["encoding-version"] = "0.1.0"
        h.attrs["complete"] = 0
        h.attrs["schema"] = f"visium_hd_{BIN_SIZE_UM:03d}um_spot2cell_mv2_v1"

        obs_columns = [
            "source_cell_index", "cell_id", "owner_bin_index",
            "cellvit_label", "cellvit_class_id", "cellvit_confidence",
            "hierarchical_class", "hierarchical_class_id",
            "pannuke_confidence", "hierarchical_confidence", "rna_capacity",
            "allocated_count_total", "rna_total",
            "centroid_x_px", "centroid_y_px",
            "centroid_array_col", "centroid_array_row",
            "expression_resolution", "resolution_level",
            "geometry_fallback", "overlap_bin_count",
            "measured_overlap_fraction_sum",
        ]
        if "direct_capture_overlap_fraction_sum" in nuclei:
            obs_columns.append("direct_capture_overlap_fraction_sum")
        if "spatial_evidence_level" in nuclei:
            obs_columns.append("spatial_evidence_level")
        obs_index = np.asarray([
            f"{args.sample}_{int(cell_id)}" for cell_id in nuclei["cell_id"]
        ], dtype="U")
        obs = init_dataframe(h, "obs", obs_index, obs_columns)
        obs_values = [
            ("source_cell_index", nuclei["source_index"]),
            ("cell_id", nuclei["cell_id"]),
            ("owner_bin_index", nuclei["owner_bin"]),
            ("cellvit_class_id", nuclei["class_id"]),
            ("cellvit_confidence", nuclei["hierarchical_confidence"]),
            ("hierarchical_class_id", nuclei["class_id"]),
            ("pannuke_confidence", nuclei["pannuke_confidence"]),
            ("hierarchical_confidence", nuclei["hierarchical_confidence"]),
            ("rna_capacity", capacity),
            ("allocated_count_total", cell_total),
            ("rna_total", cell_total),
            ("centroid_x_px", nuclei["centroid_xy"][:, 0]),
            ("centroid_y_px", nuclei["centroid_xy"][:, 1]),
            ("centroid_array_col", nuclei["centroid_uv"][:, 0]),
            ("centroid_array_row", nuclei["centroid_uv"][:, 1]),
            ("geometry_fallback", nuclei["geometry_fallback"]),
            ("overlap_bin_count", nuclei["overlap_bin_count"]),
            ("measured_overlap_fraction_sum", nuclei["measured_overlap_fraction_sum"]),
        ]
        if "direct_capture_overlap_fraction_sum" in nuclei:
            obs_values.append(
                (
                    "direct_capture_overlap_fraction_sum",
                    nuclei["direct_capture_overlap_fraction_sum"],
                )
            )
        for name, values in obs_values:
            encoded_array(obs, name, values, compression="lzf")
        write_categorical(obs, "cellvit_label", nuclei["class_id"], CLASS_NAMES)
        write_categorical(obs, "hierarchical_class", nuclei["class_id"], CLASS_NAMES)
        write_categorical(
            obs, "expression_resolution",
            resolution_by_class[nuclei["class_id"]], resolution_names,
        )
        write_categorical(
            obs, "resolution_level",
            resolution_by_class[nuclei["class_id"]], resolution_names,
        )
        if "spatial_evidence_level" in nuclei:
            write_categorical(
                obs, "spatial_evidence_level",
                nuclei["spatial_evidence_level"],
                nuclei["spatial_evidence_categories"],
            )

        var_columns = [
            "gene_name", "feature_type", "marker_classes",
            "is_marker", "is_program_gene", "is_residual_graph_gene",
            "cell_state_modeled",
        ]
        var = init_dataframe(h, "var", gene_id, var_columns)
        write_text_dataset(var, "gene_name", gene_name)
        write_text_dataset(var, "feature_type", feature_type)
        marker_strings = marker_class_strings(gene_name, marker_sets)
        write_text_dataset(var, "marker_classes", marker_strings)
        encoded_array(var, "is_marker", marker_strings != "", compression="lzf")
        encoded_array(var, "is_program_gene", np.isin(np.arange(n_genes), program_genes), compression="lzf")
        encoded_array(var, "is_residual_graph_gene", np.isin(np.arange(n_genes), residual_genes), compression="lzf")
        encoded_array(var, "cell_state_modeled", np.any(loadings != 0, axis=(0, 1)), compression="lzf")

        obsm = init_dict_group(h, "obsm")
        encoded_array(obsm, "spatial", nuclei["centroid_xy"], compression="lzf")
        encoded_array(obsm, "spatial_array_uv", nuclei["centroid_uv"], compression="lzf")
        encoded_array(obsm, "morphology", nuclei["morphology"], compression="lzf")
        encoded_array(obsm, "cellvit_embedding_pca", embedding_scores, compression="lzf")
        encoded_array(obsm, "standardized_nucleus_features", cell_features, compression="lzf")
        encoded_array(obsm, "program_activity", cell_program, compression="lzf")
        varm = init_dict_group(h, "varm")
        encoded_array(
            varm, "program_loadings",
            np.transpose(loadings, (2, 0, 1)).reshape(n_genes, N_CLASSES * MAX_PROGRAMS),
            compression="lzf",
        )
        init_dict_group(h, "obsp")
        init_dict_group(h, "varp")

        uns = init_dict_group(h, "uns")
        write_text_dataset(uns, "class_names", CLASS_NAMES)
        write_text_dataset(uns, "program_axis", [f"program_{i}" for i in range(MAX_PROGRAMS)])
        encoded_array(uns, "class_gene_probability", profiles, compression="lzf")
        encoded_array(uns, "program_count_by_class", program_count)
        encoded_array(uns, "morphology_reliability", reliability, compression="lzf")
        encoded_array(uns, "independent_profile_fraction", shrinkage, compression="lzf")
        encoded_array(uns, "bin_type_bin_index", groups["bin_index"], compression="lzf")
        encoded_array(uns, "bin_type_class_id", groups["class_id"], compression="lzf")
        encoded_array(uns, "bin_type_program_activity", group_scores, compression="lzf")
        encoded_scalar_string(uns, "summary_json", json.dumps(summary, ensure_ascii=False))
        encoded_scalar_string(uns, "marker_sets_json", json.dumps(marker_sets, ensure_ascii=False))
        pca_group = init_dict_group(uns, "cellvit_pca_by_class")
        for c, model in enumerate(pca_models):
            class_group = init_dict_group(pca_group, CLASS_NAMES[c])
            if model is not None:
                encoded_array(class_group, "components", model.components_.astype(np.float32), compression="lzf")
                encoded_array(class_group, "mean", model.mean_.astype(np.float32), compression="lzf")
                encoded_array(class_group, "explained_variance_ratio", model.explained_variance_ratio_.astype(np.float32))
        residual_group = init_dict_group(uns, "residual_pca")
        encoded_array(residual_group, "components", residual_pca.components_.astype(np.float32), compression="lzf")
        encoded_array(residual_group, "mean", residual_pca.mean_.astype(np.float32), compression="lzf")
        encoded_array(residual_group, "explained_variance_ratio", residual_pca.explained_variance_ratio_.astype(np.float32))

        x = h.create_group("X")
        x.attrs["encoding-type"] = "csr_matrix"
        x.attrs["encoding-version"] = "0.1.0"
        x.attrs["shape"] = np.asarray([n_cells, n_genes], np.int64)
        x_data = x.create_dataset(
            "data", shape=(0,), maxshape=(None,), dtype=np.float32,
            chunks=(1_000_000,), compression="lzf",
        )
        x_indices = x.create_dataset(
            "indices", shape=(0,), maxshape=(None,), dtype=np.int32,
            chunks=(1_000_000,), compression="lzf",
        )
        x_indptr = x.create_dataset("indptr", shape=(n_cells + 1,), dtype=np.int64, compression="lzf")
        x_indptr[0] = 0

        layers = init_dict_group(h, "layers")
        normalized = layers.create_group("log1p_cp10k")
        normalized.attrs["encoding-type"] = "csr_matrix"
        normalized.attrs["encoding-version"] = "0.1.0"
        normalized.attrs["shape"] = np.asarray([n_cells, n_genes], np.int64)
        normalized_data = normalized.create_dataset(
            "data", shape=(0,), maxshape=(None,), dtype=np.float32,
            chunks=(1_000_000,), compression="lzf",
        )

        responsibility_ds = allocation["allocation/responsibility"]
        record_indptr = allocation["allocation/bin_record_indptr"][:]
        geometry_csc = geometry.tocsc()
        cursor = 0
        last = time.time()
        for cell_start in range(0, n_cells, args.cell_batch_size):
            cell_end = min(cell_start + args.cell_batch_size, n_cells)
            candidate_bins = np.unique(geometry_csc.indices[
                geometry_csc.indptr[cell_start]:geometry_csc.indptr[cell_end]
            ])
            row_parts = []
            col_parts = []
            data_parts = []
            for b in candidate_bins:
                gs, ge = int(geometry.indptr[b]), int(geometry.indptr[b + 1])
                cells = geometry.indices[gs:ge]
                selected = (cells >= cell_start) & (cells < cell_end)
                if not selected.any():
                    continue
                ms, me = int(matrix.indptr[b]), int(matrix.indptr[b + 1])
                genes = matrix.indices[ms:me]
                counts = matrix.data[ms:me].astype(np.float64)
                rs, re = int(record_indptr[b]), int(record_indptr[b + 1])
                q = responsibility_ds[rs:re].reshape((len(genes), len(cells)))[:, selected]
                mass = q.astype(np.float64) * counts[:, None]
                selected_cells = cells[selected] - cell_start
                row_parts.append(np.repeat(selected_cells, len(genes)))
                col_parts.append(np.tile(genes, len(selected_cells)))
                data_parts.append(mass.T.ravel().astype(np.float32))
            if data_parts:
                batch = sparse.coo_matrix(
                    (np.concatenate(data_parts), (np.concatenate(row_parts), np.concatenate(col_parts))),
                    shape=(cell_end - cell_start, n_genes),
                ).tocsr()
                batch.sum_duplicates()
                batch.sort_indices()
            else:
                batch = sparse.csr_matrix((cell_end - cell_start, n_genes), dtype=np.float32)
            append_dataset(x_data, batch.data.astype(np.float32, copy=False))
            append_dataset(x_indices, batch.indices.astype(np.int32, copy=False))
            x_indptr[cell_start + 1:cell_end + 1] = cursor + batch.indptr[1:]
            row_total = np.asarray(batch.sum(1)).ravel().astype(np.float64)
            row_scale = np.divide(
                10000.0, row_total, out=np.zeros_like(row_total), where=row_total > 0
            )
            repeated_scale = np.repeat(row_scale, np.diff(batch.indptr))
            norm_data = np.log1p(batch.data.astype(np.float64) * repeated_scale).astype(np.float32)
            append_dataset(normalized_data, norm_data)
            cursor += batch.nnz
            if time.time() - last > 60:
                log(f"Built cell CSR through {cell_end:,}/{n_cells:,}; nnz={cursor:,}")
                last = time.time()
        normalized["indices"] = x_indices
        normalized["indptr"] = x_indptr
        h.attrs["complete"] = 1
        h.flush()
    return int(cursor)


def validate_output_files(h5ad_path: Path, allocation_path: Path,
                          expected_cells: int, expected_genes: int) -> dict:
    result = {}
    with h5py.File(allocation_path, "r") as h:
        if int(h.attrs.get("complete", 0)) != 1:
            raise RuntimeError("Allocation output lacks complete=1")
        observed = float(h["audit"].attrs["observed_total_counts"])
        assigned = float(h["audit"].attrs["assigned_cell_counts"])
        unassigned = float(h["audit"].attrs["unassigned_counts_zero_cell_bins"])
        if abs(observed - assigned - unassigned) > max(1e-4, observed * 1e-6):
            raise RuntimeError("Allocation total-count conservation failed")
        result["allocation"] = {
            "complete": True,
            "observed_total": observed,
            "assigned_total": assigned,
            "unassigned_total": unassigned,
            "max_abs_gene_error": float(h["audit"].attrs["max_abs_gene_conservation_error"]),
            "max_abs_bin_gene_error": float(h["audit"].attrs["max_abs_bin_gene_float32_error"]),
        }
    with h5py.File(h5ad_path, "r") as h:
        if int(h.attrs.get("complete", 0)) != 1:
            raise RuntimeError("H5AD output lacks complete=1")
        shape = tuple(int(x) for x in h["X"].attrs["shape"])
        if shape != (expected_cells, expected_genes):
            raise RuntimeError(f"Unexpected H5AD shape: {shape}")
        indptr = h["X/indptr"]
        if len(indptr) != expected_cells + 1:
            raise RuntimeError("H5AD CSR indptr length mismatch")
        if int(indptr[-1]) != len(h["X/data"]) or len(h["X/data"]) != len(h["X/indices"]):
            raise RuntimeError("H5AD CSR data/indices/indptr mismatch")
        result["h5ad"] = {
            "complete": True,
            "shape": list(shape),
            "nnz": int(indptr[-1]),
        }
    return result


def atomic_tsv(path: Path, frame: pd.DataFrame) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    partial = path.with_suffix(path.suffix + ".partial")
    if partial.exists():
        raise FileExistsError(f"Preserved incomplete output exists: {partial}")
    frame.to_csv(partial, sep="\t", index=False)
    os.replace(partial, path)


def write_qc_files(path: Path, summary: dict, class_audit: dict,
                   local_audit: dict, feature_audit: dict,
                   ablations: list[dict]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    class_rows = []
    for name in CLASS_NAMES:
        row = {"class": name}
        row.update(class_audit["class_coverage"][name])
        row["n_programs"] = local_audit.get(name, {}).get("n_programs", 0)
        bootstrap = local_audit.get(name, {}).get("bootstrap", {})
        row["program_loading_subspace_cosine_median"] = bootstrap.get(
            "loading_subspace_cosine_median"
        )
        row["program_z_reconstruction_correlation_median"] = bootstrap.get(
            "z_reconstruction_correlation_median"
        )
        class_rows.append(row)
    program_rows = []
    for name, audit in feature_audit.items():
        for item in audit.get("programs", []):
            program_rows.append({"class": name, **item})
    prefix = path.name[:-len(".qc.json")] if path.name.endswith(".qc.json") else path.stem
    atomic_tsv(path.parent / f"{prefix}.class_qc.tsv", pd.DataFrame(class_rows))
    atomic_tsv(path.parent / f"{prefix}.program_qc.tsv", pd.DataFrame(program_rows))
    atomic_tsv(path.parent / f"{prefix}.ablation_qc.tsv", pd.DataFrame(ablations))
    partial = path.with_suffix(path.suffix + ".partial")
    partial.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    os.replace(partial, path)


def main() -> None:
    global BIN_SIZE_UM, SPATIAL_BLOCK_SIZE_BINS
    args = parse_args()
    BIN_SIZE_UM = int(args.bin_size_um)
    SPATIAL_BLOCK_SIZE_BINS = max(1, int(round(1024.0 / BIN_SIZE_UM)))
    if args.graph_radius_bins is None:
        args.graph_radius_bins = 96.0 / BIN_SIZE_UM
    if args.subspot_radius_bins is None:
        args.subspot_radius_bins = args.graph_radius_bins
    paths = input_paths(args)
    for name, path in {
        **paths, "positions_tsv": args.positions_tsv,
        "marker_catalog": args.marker_catalog,
    }.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing {name}: {path}")
    if args.calibration_cell_h5ad is not None and not args.calibration_cell_h5ad.exists():
        raise FileNotFoundError(
            f"Missing calibration_cell_h5ad: {args.calibration_cell_h5ad}"
        )
    allocation_output = args.allocation_output or args.output.with_name(
        f"{args.sample}.allocation.h5"
    )
    qc_json = args.qc_json or args.output.with_name(f"{args.sample}.qc.json")
    output_partial = args.output.with_suffix(args.output.suffix + ".partial")
    allocation_partial = allocation_output.with_suffix(allocation_output.suffix + ".partial")
    for path in (args.output, allocation_output, qc_json, output_partial, allocation_partial):
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    allocation_output.parent.mkdir(parents=True, exist_ok=True)

    (
        matrix, library, barcodes, gene_id, gene_name, feature_type,
        source_bin_index, zero_count,
    ) = load_matrix(paths["matrix"])
    with h5py.File(paths["cell_features"], "r") as h:
        all_centroid = h["centroid_xy"][:].astype(np.float32)
    positions, owner_all, eligible, centroid_uv_all, affine, spatial_audit = align_positions(
        args.positions_tsv, barcodes, all_centroid
    )
    del all_centroid
    nuclei = load_eligible_nuclei(
        paths["cell_features"], owner_all, eligible, centroid_uv_all
    )
    del centroid_uv_all, owner_all, eligible
    geometry, geometry_audit = build_contour_overlap(nuclei, affine, matrix.shape[1])
    geometry_audit["estimated_allocation_responsibility_records"] = int(
        np.dot(
            np.diff(matrix.indptr).astype(np.int64),
            np.diff(geometry.indptr).astype(np.int64),
        )
    )
    base = base_capacity(nuclei)
    base_composition, base_mass, n_cells_per_bin = census_from_geometry(
        geometry, nuclei["class_id"], base
    )
    occupied = n_cells_per_bin > 0
    summary = {
        "schema": f"visium_hd_{BIN_SIZE_UM:03d}um_spot2cell_mv2_summary_v1",
        "algorithm": (
            "CellViT fixed census -> contour overlap -> marker-anchored independent "
            "class profiles -> edge-aware local type states -> continuous subspot "
            "state interpolation -> permutation-gated morphology residual -> "
            "two-stage KL-regularized exact occupied-bin gene projection"
        ),
        "sample": args.sample,
        "bin_size_um": BIN_SIZE_UM,
        "spatial_block_size_bins": SPATIAL_BLOCK_SIZE_BINS,
        "spatial_block_size_um": SPATIAL_BLOCK_SIZE_BINS * BIN_SIZE_UM,
        "graph_radius_bins": float(args.graph_radius_bins),
        "graph_radius_um": float(args.graph_radius_bins * BIN_SIZE_UM),
        "subspot_radius_bins": float(args.subspot_radius_bins),
        "subspot_radius_um": float(args.subspot_radius_bins * BIN_SIZE_UM),
        "within_type_temperature": float(args.within_type_temperature),
        "independent_sample_fit": True,
        "no_scrna_reference": True,
        "no_spatial_extrapolation": True,
        "occupied_bin_background_component": False,
        "spatial_rule": (
            f"centroid inside positive-count in-tissue measured {BIN_SIZE_UM} um domain; "
            "evidence weighted by nucleus-contour overlap with measured bins"
        ),
        "n_genes": int(matrix.shape[0]),
        "n_positive_measured_bins": int(matrix.shape[1]),
        "n_zero_count_filtered_bins_excluded": int(zero_count),
        "n_eligible_cells": int(geometry.shape[1]),
        "n_bins_with_overlapping_eligible_cells": int(occupied.sum()),
        "n_bins_without_overlapping_eligible_cells": int((~occupied).sum()),
        "spatial_audit": spatial_audit,
        "geometry_audit": geometry_audit,
        "class_cell_counts": {
            CLASS_NAMES[int(value)]: int(count)
            for value, count in zip(*np.unique(nuclei["class_id"], return_counts=True))
        },
    }
    log(json.dumps(summary, ensure_ascii=False))
    if args.audit_only:
        return

    marker_sets, catalog = load_marker_sets(args.marker_catalog, gene_name)
    occupied_gene = np.asarray(matrix @ occupied.astype(np.float64)).ravel()
    global_mean = occupied_gene / max(float(occupied_gene.sum()), 1.0)
    prior = marker_prior(global_mean, gene_name, marker_sets)
    fit_genes = select_fit_genes(matrix, gene_name, marker_sets)
    blocks = spatial_blocks(positions)
    alpha, independent_profiles, alpha_cv = fit_alpha_and_profiles(
        matrix, library, base_mass, prior, fit_genes, blocks,
        args.class_prior_strength, args.cv_folds, args.seed,
    )
    capacity = (base.astype(np.float64) * alpha[nuclei["class_id"]]).astype(np.float32)
    composition, census_mass, n_cells_per_bin = census_from_geometry(
        geometry, nuclei["class_id"], capacity
    )
    shrinkage, shrinkage_audit = hierarchy_shrinkage_cv(
        matrix, library, composition, prior, fit_genes, blocks,
        args.class_prior_strength, args.cv_folds,
    )
    profiles = apply_hierarchy_shrinkage(
        independent_profiles, shrinkage, composition, global_mean
    )

    log("Learning class-internal CellViT embedding PCs and residualized morphology features")
    cell_features, embedding_scores, embedding_pca = learn_cell_features(
        paths["cell_features"], nuclei, args.embedding_pcs,
        args.pca_sample_per_class, args.seed,
    )
    groups = aggregate_bin_type_features(
        geometry, nuclei["class_id"], capacity, cell_features
    )
    class_audit = profile_audit(
        composition, groups, profiles, gene_name, marker_sets, alpha, shrinkage
    )
    no_marker_prior = np.tile(global_mean[None, :], (N_CLASSES, 1))
    no_marker_profiles = fit_profiles(
        matrix, composition, library, no_marker_prior,
        args.class_prior_strength,
    ).astype(np.float32)
    no_marker_audit = profile_audit(
        composition, groups, no_marker_profiles, gene_name, marker_sets,
        alpha, np.ones(N_CLASSES, np.float32),
    )
    summary.update({
        "marker_catalog_schema_version": catalog.get("schema_version"),
        "marker_catalog_workbook_sha256": catalog.get("source_workbook_sha256"),
        "marker_genes_present": {name: len(genes) for name, genes in marker_sets.items()},
        "alpha": alpha.tolist(),
        "alpha_cv": alpha_cv,
        "hierarchy_shrinkage_cv": shrinkage_audit,
        "class_fit_audit": class_audit,
    })

    program_genes, residual_genes, variability = select_program_genes(
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
    residual_scores, residual_pca = compute_residual_pcs(
        matrix, library, composition, profiles, graph_genes,
        occupied, args.residual_pcs, args.seed,
    )
    component_id = measured_components(positions)
    summary["measured_domain_connected_components"] = int(len(np.unique(component_id)))
    group_scores, loadings, program_count, local_audit = learn_local_programs(
        matrix, library, positions, component_id, composition, profiles,
        residual_scores, state_fit_genes, program_genes, groups,
        args.graph_neighbors, args.graph_radius_bins, args.cv_folds,
        args.bootstrap_replicates, args.seed,
    )
    morphology_link_scores, feature_coef, reliability, feature_audit = learn_feature_program_mapping(
        nuclei, capacity, cell_features, groups, composition, residual_scores,
        positions, group_scores, program_count, args.permutations,
        args.cv_folds, args.seed + 7000,
    )
    local_state_reliability = local_program_reliability(local_audit)
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
    spatial_offsets, subspot_audit = interpolate_subspot_offsets(
        nuclei, geometry, capacity, groups, positions, component_id,
        composition, residual_scores, group_scores, program_count, local_audit,
        args.subspot_neighbors, args.subspot_radius_bins,
    )
    spatial_scale, morphology_scale, calibration_audit = calibrate_offset_scales(
        args.calibration_cell_h5ad, args.calibration_max_cells,
        nuclei, geometry, capacity, positions, gene_name, profiles, loadings,
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
    summary["global_profile_bootstrap"] = global_profile_bootstrap(
        matrix, library, composition, prior, profiles, fit_genes, blocks,
        args.class_prior_strength, args.bootstrap_replicates, args.seed + 9000,
    )
    if args.skip_ablations:
        ablations = []
        summary["ablations_skipped"] = True
    else:
        ablations = run_ablations(
            matrix, library, positions, component_id, composition, profiles,
            residual_scores, state_fit_genes, groups, class_audit,
            no_marker_audit, feature_audit, args.graph_neighbors,
            args.graph_radius_bins, args.cv_folds, args.seed + 11000,
        )
        summary["ablations"] = ablations

    log("Hashing immutable inputs for provenance")
    provenance_paths = {
        "matrix": paths["matrix"],
        "cell_features": paths["cell_features"],
        "marker_catalog": args.marker_catalog,
    }
    if args.calibration_cell_h5ad is not None:
        provenance_paths["heterogeneity_calibration_reference"] = (
            args.calibration_cell_h5ad
        )
    summary["input_provenance"] = {
        name: {
            "path": str(path), "bytes": int(path.stat().st_size),
            "sha256": sha256_file(path),
        }
        for name, path in provenance_paths.items()
    }
    summary["outputs"] = {
        "cell_resolved_h5ad": str(args.output),
        "allocation_h5": str(allocation_output),
        "qc_json": str(qc_json),
    }

    cell_total, cell_program, conservation = write_allocation_file(
        allocation_partial, args, paths, matrix, library, barcodes,
        source_bin_index, gene_id, gene_name, geometry, nuclei, capacity,
        link_scores, profiles, loadings, program_count, summary,
    )
    summary["count_conservation"] = conservation
    nnz = write_h5ad(
        output_partial, allocation_partial, args, matrix, barcodes,
        gene_id, gene_name, feature_type, positions, geometry, nuclei,
        capacity, cell_features, embedding_scores, cell_total, cell_program,
        profiles, loadings, program_count, reliability, shrinkage,
        marker_sets, program_genes, graph_genes, groups, group_scores,
        embedding_pca, residual_pca, summary,
    )
    summary["cell_matrix_nnz"] = nnz
    validation = validate_output_files(
        output_partial, allocation_partial, geometry.shape[1], matrix.shape[0]
    )
    summary["output_validation"] = validation
    write_qc_files(
        qc_json, summary, class_audit, local_audit, feature_audit, ablations
    )
    os.replace(allocation_partial, allocation_output)
    os.replace(output_partial, args.output)
    log(f"Completed {args.sample}: {args.output}")


if __name__ == "__main__":
    main()
