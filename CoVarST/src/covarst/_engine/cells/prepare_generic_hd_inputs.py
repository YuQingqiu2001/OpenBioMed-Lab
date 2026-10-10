#!/usr/bin/env python3
"""Prepare genuine PanNuke-only HD inputs without changing source data.

Modes: metadata (no segmentation needed), coarse (native16-only census/geometry/
profiles/features), direct2um (independent measured evaluation reference), test.
Dependencies: existing ST_GJH numpy/scipy/pandas/h5py/pyarrow/anndata/sklearn.
The coarse path has no 2um argument and never reads a fine-resolution reference.
CPU affinity and BLAS are limited to eight cores; no GPU package is imported.
"""
from __future__ import annotations

import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "8"
if hasattr(os, "sched_setaffinity"):
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:8])

import argparse
import hashlib
import importlib.util
import json
import sys
import tempfile
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import sparse
from scipy.spatial import cKDTree

SCHEMA = "generic_hd_pannuke_native16_inputs_v1"
CLASSES = ["Neoplastic", "Inflammatory", "Connective", "Dead", "Epithelial"]
RAW_PANNUKE = {0: "Background", 1: "Neoplastic", 2: "Inflammatory", 3: "Connective", 4: "Dead", 5: "Epithelial"}
# Broad lineage anchors only. Neoplastic/normal epithelial status comes from
# the fixed image census, not a transcriptomic marker claim.
DEFAULT_MARKERS = {
    "Neoplastic": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Inflammatory": ["PTPRC", "CD3D", "CD3E", "CD79A", "MS4A1", "LST1", "TYROBP", "FCER1G"],
    "Connective": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PECAM1", "VWF"],
    "Dead": [],
}
POSITION_COLUMNS = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]


def decode(values):
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in values])


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            result.update(block)
    return result.hexdigest()


def fingerprint(path, do_hash=True):
    path = Path(path)
    answer = {"path": str(path.resolve()), "bytes": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns}
    if do_hash:
        answer["sha256"] = digest(path)
    return answer


def write_json(path, payload):
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    partial = path.with_name(path.name + ".partial")
    if partial.exists():
        raise FileExistsError(partial)
    partial.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    partial.replace(path)


def load_core(path):
    spec = importlib.util.spec_from_file_location("generic_hd_original_core", path / "cell_resolved_v3.py")
    core = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(core)
    # Runtime schema parameters only; on-disk mathematical source is untouched.
    core.CLASS_NAMES = CLASSES.copy()
    core.N_CLASSES = len(CLASSES)
    core.PARENT_MEMBERS = {0: (0, 4), 1: (1,), 2: (2,), 3: (), 4: (0, 4)}
    core.BIN_SIZE_UM = 16
    return core


def check_input_resolution(args, expected=None):
    factors = json.loads(args.scalefactors.read_text(encoding="utf-8"))
    resolution = float(factors["bin_size_um"])
    if expected is not None and resolution != expected:
        raise ValueError(f"{args.mode} requires {expected}um input, got {resolution}")
    if resolution not in (2.0, 16.0):
        raise ValueError("Only native2/native16 resolutions are supported")
    with h5py.File(args.matrix, "r") as h:
        first = decode(h["matrix/barcodes"][:min(10, h["matrix/barcodes"].shape[0])])
    if not all(value.startswith(f"s_{int(resolution):03d}um_") for value in first):
        raise ValueError("Barcode resolution disagrees with scalefactors")
    return factors, int(resolution)


def positions_for_barcodes(path, barcodes):
    positions = pd.read_parquet(path, columns=POSITION_COLUMNS)
    if positions.barcode.duplicated().any():
        raise ValueError("Duplicate position barcodes")
    positions = positions.set_index("barcode", drop=False)
    missing = ~pd.Index(barcodes).isin(positions.index)
    if missing.any():
        raise ValueError(f"{int(missing.sum())} measured barcodes lack positions")
    return positions.loc[barcodes].reset_index(drop=True)


def affine_from_positions(positions, mpp, resolution):
    design = np.column_stack([np.ones(len(positions)), positions.array_col, positions.array_row])
    observed = positions[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(float)
    coefficient = np.linalg.lstsq(design, observed, rcond=None)[0]
    transform = coefficient[1:].T
    residual = observed - design @ coefficient
    maximum = float(np.abs(residual).max())
    edge_um = np.linalg.norm(transform, axis=0) * mpp
    if abs(np.linalg.det(transform)) < 1e-8 or maximum > 2:
        raise ValueError(f"Invalid image affine: max residual {maximum} px")
    if not np.allclose(edge_um, resolution, rtol=.02, atol=.05):
        raise ValueError(f"Image MPP and array grid disagree: {edge_um}um versus {resolution}um")
    return transform, coefficient[0], {"transform_array_to_pixel": transform.tolist(), "intercept_pixel": coefficient[0].tolist(), "max_affine_residual_px": maximum, "grid_edge_um": edge_um.tolist()}


def metadata(args):
    factors, resolution = check_input_resolution(args)
    with h5py.File(args.matrix, "r") as h:
        group = h["matrix"]
        barcodes = decode(group["barcodes"][:])
        shape = [int(x) for x in group["shape"][:]]
        nnz = int(group["data"].shape[0])
        features = decode(group["features/name"][:])
        feature_types = sorted(set(decode(group["features/feature_type"][:])))
        pointer = group["indptr"][:]
    positions = positions_for_barcodes(args.positions, barcodes)
    _, _, affine = affine_from_positions(positions, float(factors["microns_per_pixel"]), resolution)
    result = {"schema": SCHEMA, "mode": "metadata", "sample": args.sample, "complete": True, "bin_size_um": resolution, "raw_shape_genes_by_bins": shape, "raw_nnz": nnz, "n_unique_symbols": len(set(features)), "feature_types": feature_types, "n_filtered_barcodes_outside_tissue": int((positions.in_tissue != 1).sum()), "n_structurally_empty_columns": int(np.count_nonzero(np.diff(pointer) == 0)), "positive_library_not_claimed_from_metadata": True, "microns_per_pixel": float(factors["microns_per_pixel"]), "affine": affine, "features_available": bool(args.features and args.features.exists()), "no_cell_expression_created": True, "inputs": {name: fingerprint(getattr(args, name), False) for name in ("matrix", "positions", "scalefactors")}}
    path = args.output_dir / f"{args.sample}.{resolution:03d}um.metadata.json"
    write_json(path, result)
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)
    return result


def load_positive_tissue(args, core):
    values = core.load_matrix(args.matrix)
    matrix, library, barcodes, gene_id, gene_name, feature_type, original_columns, n_zero = values
    positions = positions_for_barcodes(args.positions, barcodes)
    keep = positions.in_tissue.to_numpy() == 1
    excluded = int((~keep).sum())
    if excluded:
        matrix = matrix[:, keep].tocsc()
        library = library[keep]; barcodes = barcodes[keep]; original_columns = original_columns[keep]
        positions = positions.loc[keep].reset_index(drop=True)
    if not matrix.shape[1] or np.any(library <= 0):
        raise ValueError("No positive measured tissue bins")
    if np.any(matrix.data < 0) or not np.isfinite(matrix.data).all() or np.any(matrix.data != np.floor(matrix.data)):
        raise ValueError("Native matrix is not nonnegative integer measured counts")
    audit = {"raw_zero_library_columns_excluded": int(n_zero), "positive_non_tissue_columns_excluded": excluded, "positive_in_tissue_bins": matrix.shape[1], "raw_gene_features": matrix.shape[0], "total_counts_in_domain": int(matrix.sum()), "domain_rule": "filtered matrix membership AND in_tissue=1 AND strictly positive native all-gene library"}
    return matrix, library, barcodes, gene_id, gene_name, feature_type, original_columns, positions, audit


def feature_metadata(path, mpp):
    if path is None or not path.is_file():
        raise FileNotFoundError("A real completed PanNuke CellViT feature H5 is required; metadata mode can run before segmentation")
    with h5py.File(path, "r") as h:
        if int(h.attrs.get("feature_complete", 0)) != 1:
            raise ValueError("Feature H5 is incomplete")
        if bool(h.attrs.get("lizard_applied", False)):
            raise ValueError("This adapter requires PanNuke-only features")
        declared = {int(k): v for k, v in json.loads(h.attrs["pannuke_labels_json"]).items()}
        if declared != RAW_PANNUKE:
            raise ValueError(f"Unexpected PanNuke coding: {declared}")
        if not np.isclose(float(h.attrs["microns_per_pixel"]), mpp, rtol=1e-5):
            raise ValueError("Feature MPP and Space Ranger MPP disagree")
        ids = h["cell_id"][:].astype(np.int64)
        xy = h["centroid_xy"][:].astype(np.float64)
        classes = h["pannuke_class_id"][:].astype(np.int64)
        if len(np.unique(ids)) != len(ids) or xy.shape != (len(ids), 2) or len(classes) != len(ids):
            raise ValueError("Feature rows/IDs disagree")
        if not np.isfinite(xy).all() or not np.isin(classes, np.arange(6)).all():
            raise ValueError("Invalid centroid/class metadata")
        if h["cellvit_embedding_float16"].shape != (len(ids), 1280):
            raise ValueError("A complete 1280-D feature bank is required")
    return ids, xy, classes


def load_nuclei(path, owner_all, eligible, uv_all):
    indices = np.flatnonzero(eligible).astype(np.int64)
    with h5py.File(path, "r") as h:
        confidence = h["pannuke_confidence"][indices].astype(np.float32)
        raw_class = h["pannuke_class_id"][indices].astype(np.int64)
        return {"source_index": indices, "owner_bin": owner_all[indices], "centroid_xy": h["centroid_xy"][indices].astype(np.float32), "centroid_uv": uv_all[indices].astype(np.float32), "cell_id": h["cell_id"][indices].astype(np.int64), "class_id": (raw_class - 1).astype(np.uint8), "pannuke_class_id": raw_class.astype(np.uint8), "morphology": h["morphology"][indices].astype(np.float32), "pannuke_confidence": confidence, "hierarchical_confidence": confidence, "morphology_names": json.loads(h.attrs["morphology_feature_names_json"]), "all_contour_offsets": h["contour_offsets"][:].astype(np.int64), "all_contour_xy": h["contour_xy"][:].astype(np.float64)}


def write_h5ad_atomic(data, path):
    partial = path.with_name(path.name + ".partial")
    if path.exists() or partial.exists():
        raise FileExistsError(path)
    data.write_h5ad(partial, compression="gzip")
    with h5py.File(partial, "r+") as h:
        # AnnData 0.10 promotes repeated symbols to categoricals. The original
        # factor loader explicitly requires a flat string array at these keys.
        # Rewrite only newly created temporary outputs, preserving H5AD encoding.
        for name in ("gene_name", "gene_symbol"):
            key = f"var/{name}"
            if key in h and isinstance(h[key], h5py.Group):
                group = h[key]
                if "codes" not in group or "categories" not in group:
                    raise ValueError(f"Unexpected variable encoding at {key}")
                codes = group["codes"][:]
                if np.any(codes < 0):
                    raise ValueError(f"Missing gene identifier at {key}")
                values = decode(group["categories"][:])[codes]
                del h[key]
                dataset = h.create_dataset(key, data=np.asarray(values, object), dtype=h5py.string_dtype("utf-8"))
                dataset.attrs.update({"encoding-type": "string-array", "encoding-version": "0.2.0"})
        h.attrs["complete"] = 1
        h.attrs["schema"] = SCHEMA
    partial.replace(path)


def write_ledger(path, matrix, barcodes, gene_id, gene_name, original_columns, geometry):
    partial = path.with_name(path.name + ".partial")
    if path.exists() or partial.exists():
        raise FileExistsError(path)
    text = h5py.string_dtype("utf-8")
    with h5py.File(partial, "x") as h:
        h.attrs.update(schema=SCHEMA, complete=0, raw_measured_counts=True, source_resolution_um=16, no_fine_expression_used=True)
        observed = h.create_group("observed")
        for name, values in {"shape": np.asarray(matrix.shape), "positive_source_bin_index": original_columns, "bin_indptr": matrix.indptr, "gene_index": matrix.indices, "count": matrix.data}.items():
            observed.create_dataset(name, data=values, compression="lzf")
        for name, values in {"barcode": barcodes, "gene_id": gene_id, "gene_name": gene_name}.items():
            observed.create_dataset(name, data=np.asarray(values, object), dtype=text)
        group = h.create_group("geometry")
        for name, values in {"shape": np.asarray(geometry.shape), "indptr": geometry.indptr, "cell_index": geometry.indices, "overlap_fraction": geometry.data}.items():
            group.create_dataset(name, data=values, compression="lzf")
        h.attrs["complete"] = 1
    partial.replace(path)


def coarse(args):
    started = time.time()
    factors, _ = check_input_resolution(args, 16)
    mpp = float(factors["microns_per_pixel"])
    _, xy, raw_class = feature_metadata(args.features, mpp)
    core = load_core(args.pipeline_dir)
    matrix, library, barcodes, gene_id, genes, feature_type, original_columns, positions, domain_audit = load_positive_tissue(args, core)
    _, _, affine_audit = affine_from_positions(positions, mpp, 16)
    positions_path = args.output_dir / f"{args.sample}.native16.positions.tsv.gz"
    if positions_path.exists():
        raise FileExistsError(positions_path)
    positions.to_csv(positions_path, sep="\t", index=False, compression="gzip")
    positions, owner, eligible, uv, affine, spatial_audit = core.align_positions(positions_path, barcodes, xy)
    eligible &= raw_class != 0
    if not eligible.any():
        raise ValueError("No non-background nuclei in native16 positive tissue domain")
    nuclei = load_nuclei(args.features, owner, eligible, uv)
    geometry, geometry_audit = core.build_contour_overlap(nuclei, affine, matrix.shape[1])
    base = core.base_capacity(nuclei)
    _, mass, n_cells_per_bin = core.census_from_geometry(geometry, nuclei["class_id"], base)
    occupied = n_cells_per_bin > 0
    marker_sets = DEFAULT_MARKERS.copy() if args.marker_map is None else json.loads(args.marker_map.read_text(encoding="utf-8"))
    if set(marker_sets) != set(CLASSES):
        raise ValueError("Marker map must specify exactly the five PanNuke cell classes")
    marker_sets = {key: sorted(set(value).intersection(genes)) for key, value in marker_sets.items()}
    occupied_counts = np.asarray(matrix @ occupied.astype(np.float64)).ravel()
    global_mean = occupied_counts / occupied_counts.sum()
    prior = core.marker_prior(global_mean, genes, marker_sets)
    fit_genes = core.select_fit_genes(matrix, genes, marker_sets)
    blocks = core.spatial_blocks(positions)
    alpha, profiles, alpha_audit = core.fit_alpha_and_profiles(matrix, library, mass, prior, fit_genes, blocks, args.class_prior_strength, args.cv_folds, args.seed)
    capacity = (base * alpha[nuclei["class_id"]]).astype(np.float32)
    composition, _, _ = core.census_from_geometry(geometry, nuclei["class_id"], capacity)
    shrinkage, shrinkage_audit = core.hierarchy_shrinkage_cv(matrix, library, composition, prior, fit_genes, blocks, args.class_prior_strength, args.cv_folds)
    profiles = core.apply_hierarchy_shrinkage(profiles, shrinkage, composition, global_mean)
    if not np.isfinite(profiles).all() or np.any(profiles < 0) or not np.allclose(profiles.sum(1), 1, atol=2e-6):
        raise ValueError("Global class probabilities fail full-gene normalization")
    features, embedding_pcs, _ = core.learn_cell_features(args.features, nuclei, args.embedding_pcs, args.pca_sample_per_class, args.seed)
    if not np.isfinite(features).all():
        raise ValueError("Nucleus feature standardization is not finite")
    # The panel uses only native16 measured abundance plus broad lineage anchors.
    selected_features = core.select_fit_genes(matrix, genes, marker_sets, limit=args.panel_size)
    panel = list(dict.fromkeys(genes[selected_features].tolist()))
    if len(panel) < 20:
        raise ValueError("Fewer than 20 usable native16 panel genes")
    inputs = {name: fingerprint(getattr(args, name)) for name in ("matrix", "positions", "scalefactors", "features")}
    summary = {"schema": SCHEMA, "sample": args.sample, "class_names": CLASSES, "cell_label_source": "PanNuke raw IDs1..5 mapped to0..4; Background excluded; no Lizard", "class_counts": {name: int(np.count_nonzero(nuclei["class_id"] == j)) for j, name in enumerate(CLASSES)}, "domain": domain_audit, "spatial_audit": spatial_audit, "geometry_audit": geometry_audit, "affine_audit": affine_audit, "n_cells": len(nuclei["cell_id"]), "n_genes": len(genes), "n_bins_without_nuclei": int((~occupied).sum()), "unassigned_count_mass_zero_nucleus_bins": int(matrix[:, ~occupied].sum()), "marker_sets": marker_sets, "marker_policy": "explicit broad pan-lineage anchors only; no CRC organ markers; same anchors for neoplastic and normal epithelial image classes", "panel_genes": panel, "panel_selection": "native16 all-gene measured abundance plus broad anchors; no 2um read", "alpha": alpha.tolist(), "alpha_cv": alpha_audit, "hierarchy_shrinkage": shrinkage.tolist(), "hierarchy_shrinkage_cv": shrinkage_audit, "full_gene_profile_sum_max_error": float(np.abs(profiles.sum(1) - 1).max()), "heterogeneity_calibration": {"calibrated": False, "reason": "high-resolution inputs structurally absent from coarse preparation"}, "2um_expression_used_for_fit": False, "no_cell_expression_matrix_created": True, "input_provenance": inputs, "source_core_sha256": digest(args.pipeline_dir / "cell_resolved_v3.py"), "adapter_sha256": digest(__file__), "features_basis": "same original per-class embedding PCA plus nuisance-residualized morphology/confidence/neighborhood; CPU-only", "CV_interpretation": "conditional/transductive coarse-only profile estimation, not independent held-out validation", "elapsed_seconds": time.time() - started}
    obs = pd.DataFrame({"cell_id": nuclei["cell_id"], "source_cell_index": nuclei["source_index"], "owner_bin_index": nuclei["owner_bin"], "cellvit_class_id": nuclei["class_id"], "cellvit_label": pd.Categorical([CLASSES[int(x)] for x in nuclei["class_id"]], categories=CLASSES), "pannuke_class_id": nuclei["pannuke_class_id"], "rna_capacity": capacity, "pannuke_confidence": nuclei["pannuke_confidence"], "centroid_x_px": nuclei["centroid_xy"][:, 0], "centroid_y_px": nuclei["centroid_xy"][:, 1]}, index=pd.Index([f"{args.sample}_nucleus_{x}" for x in nuclei["cell_id"]], name="cell"))
    var = pd.DataFrame({"gene_name": genes, "feature_type": feature_type}, index=pd.Index(gene_id, name="gene_id"))
    census = ad.AnnData(obs=obs, var=var, shape=(len(obs), len(var)))
    census.obsm["spatial"] = nuclei["centroid_xy"]
    census.obsm["spatial_array_uv"] = nuclei["centroid_uv"]
    census.obsm["standardized_nucleus_features"] = features
    census.obsm["cellvit_embedding_pca"] = embedding_pcs
    census.obsm["morphology"] = nuclei["morphology"]
    census.uns["class_names"] = np.asarray(CLASSES, object)
    census.uns["class_gene_probability"] = profiles
    census.uns["summary_json"] = json.dumps(summary, ensure_ascii=False)
    census.uns["marker_sets_json"] = json.dumps(marker_sets)
    cell_path = args.output_dir / f"{args.sample}.native16.census_profiles.h5ad"
    ledger_path = args.output_dir / f"{args.sample}.native16.raw_ledger.h5"
    panel_path = args.output_dir / f"{args.sample}.native16.panel_genes.h5ad"
    write_h5ad_atomic(census, cell_path)
    write_ledger(ledger_path, matrix, barcodes, gene_id, genes, original_columns, geometry)
    panel_data = ad.AnnData(shape=(0, len(panel)), obs=pd.DataFrame(index=pd.Index([], dtype=str)), var=pd.DataFrame({"gene_name": np.asarray(panel)}, index=pd.Index(panel, name="gene")))
    panel_data.uns["selection"] = "native16-only gene identifiers; no expression matrix"
    write_h5ad_atomic(panel_data, panel_path)
    summary["outputs"] = {"cell_h5ad": str(cell_path), "allocation": str(ledger_path), "panel_gene_source": str(panel_path), "positions_tsv": str(positions_path)}
    summary["fit_command"] = f"python fit_native16_bayesgraph.py --sample {args.sample} --class-schema input --cell-h5ad '{cell_path}' --allocation '{ledger_path}' --raw16 '{args.matrix}' --positions '{args.positions}' --scalefactors '{args.scalefactors}' --panel-gene-source '{panel_path}' --output-dir '<FRESH_FACTOR_DIR>'"
    write_json(args.output_dir / f"{args.sample}.native16.preparation.qc.json", summary)
    print(json.dumps({"complete": True, "sample": args.sample, "n_cells": len(obs), "outputs": summary["outputs"], "fit_command": summary["fit_command"]}, ensure_ascii=False, indent=2), flush=True)


def finite_assignment(rows, cols, uv, radius_bins, neighbors=3, scale=.75, quantile=.99, chunk=250000):
    """Original finite Voronoi/octagonal rule, with <=8 KDTree workers."""
    if len(uv) < 2:
        raise ValueError("At least two genuine nuclei are required")
    tree = cKDTree(uv)
    k = min(neighbors + 1, len(uv))
    distance = tree.query(uv, k=k, workers=8)[0]
    spacing = np.median(distance[:, 1:], axis=1)
    raw_radius = np.maximum(radius_bins, scale * spacing)
    cap = float(np.quantile(raw_radius, quantile))
    radius = np.minimum(raw_radius, cap)
    assignment = np.full(len(rows), -1, np.int32)
    for start in range(0, len(rows), chunk):
        end = min(start + chunk, len(rows))
        point = np.column_stack([cols[start:end], rows[start:end]])
        _, nearest = tree.query(point, k=1, workers=8)
        delta = point - uv[nearest]
        octagonal = np.maximum.reduce([abs(delta[:, 0]), abs(delta[:, 1]), abs(delta[:, 0] + delta[:, 1]) / np.sqrt(2), abs(delta[:, 0] - delta[:, 1]) / np.sqrt(2)])
        valid = octagonal <= radius[nearest] + 1e-7
        assignment[start:end][valid] = nearest[valid]
    return assignment, spacing, radius, cap


def direct(args):
    started = time.time()
    factors, _ = check_input_resolution(args, 2)
    mpp = float(factors["microns_per_pixel"])
    all_ids, xy, raw_classes = feature_metadata(args.features, mpp)
    core = load_core(args.pipeline_dir)
    matrix, library, barcodes, gene_id, genes, feature_type, original_columns, positions, domain_audit = load_positive_tissue(args, core)
    transform, intercept, affine = affine_from_positions(positions, mpp, 2)
    uv_all = (xy - intercept) @ np.linalg.inv(transform).T
    # Exactly the positive native2 tissue mask. No 16um expression/profile is read.
    rr = np.rint(uv_all[:, 1]).astype(np.int64); cc = np.rint(uv_all[:, 0]).astype(np.int64)
    keys = set(zip(positions.array_row.astype(int), positions.array_col.astype(int)))
    eligible = np.fromiter(((int(r), int(c)) in keys for r, c in zip(rr, cc)), dtype=bool, count=len(rr)) & (raw_classes > 0)
    source_index = np.flatnonzero(eligible)
    if not len(source_index):
        raise ValueError("No nucleus centers in positive native2 tissue domain")
    with h5py.File(args.features, "r") as h:
        area = h["morphology"][source_index, 0].astype(np.float64)
        confidence = h["pannuke_confidence"][source_index].astype(np.float32)
    if not np.isfinite(area).all() or np.any(area <= 0):
        raise ValueError("Finite-domain reference requires valid positive nucleus areas")
    edge_px = .5 * np.linalg.norm(transform, axis=0).sum()
    radius_bins = np.sqrt(area / np.pi) / edge_px
    assignment, spacing, radius, cap = finite_assignment(positions.array_row.to_numpy(float), positions.array_col.to_numpy(float), uv_all[source_index], radius_bins, args.local_neighbors, args.radius_scale, args.radius_cap_quantile)
    assigned = np.flatnonzero(assignment >= 0)
    ownership = sparse.csc_matrix((np.ones(len(assigned), np.int64), (assigned, assignment[assigned])), shape=(matrix.shape[1], len(source_index)))
    cell_counts = (matrix.astype(np.int64) @ ownership).T.tocsr()
    cell_counts.sum_duplicates(); cell_counts.sort_indices(); cell_counts.eliminate_zeros()
    expected_gene = np.asarray(matrix[:, assigned].sum(axis=1)).ravel().astype(np.int64)
    actual_gene = np.asarray(cell_counts.sum(axis=0)).ravel().astype(np.int64)
    if not np.array_equal(expected_gene, actual_gene):
        raise ValueError("Direct reference all-gene integer conservation failure")
    all_counts = np.asarray(cell_counts.sum(axis=1)).ravel().astype(np.int64)
    ids = all_ids[source_index]
    obs = pd.DataFrame({"cell_id": ids, "source_cell_index": source_index, "pannuke_class_id": raw_classes[source_index], "cellvit_class_id": raw_classes[source_index] - 1, "pannuke_confidence": confidence, "centroid_x_px": xy[source_index, 0], "centroid_y_px": xy[source_index, 1], "centroid_array_col_2um": uv_all[source_index, 0], "centroid_array_row_2um": uv_all[source_index, 1], "assigned_all_gene_umi": all_counts, "assigned_2um_bin_count": np.bincount(assignment[assigned], minlength=len(source_index)), "finite_domain_radius_um": radius * 2, "local_spacing_um": spacing * 2}, index=pd.Index([f"{args.sample}_nucleus_{x}" for x in ids], name="cell"))
    var = pd.DataFrame({"gene_symbol": genes, "feature_type": feature_type}, index=pd.Index(gene_id, name="gene_id"))
    result = ad.AnnData(X=cell_counts, obs=obs, var=var)
    result.obsm["spatial"] = xy[source_index].astype(np.float32)
    result.uns["X_semantics"] = "integer measured2um counts summed inside finite nucleus domains; unsmoothed"
    result.uns["normalization"] = "full native-gene library before duplicate symbol collapse or panel subset"
    result.uns["not_independent_single_cell_ground_truth"] = True
    output = args.output_dir / f"{args.sample}.direct2um.all_gene_counts.h5ad"
    write_h5ad_atomic(result, output)
    assignment_path = args.output_dir / f"{args.sample}.direct2um.assignment.h5"
    with h5py.File(assignment_path, "x") as h:
        h.attrs.update(complete=0, schema=SCHEMA, reference_only=True, no_fitting=True)
        h.create_dataset("cell_id", data=ids)
        h.create_dataset("positive_source_bin_index", data=original_columns, compression="lzf")
        h.create_dataset("bin_cell_index", data=assignment, compression="lzf")
        h.attrs["complete"] = 1
    qc = {"schema": SCHEMA, "complete": True, "sample": args.sample, "mode": "direct2um", "reference_only": True, "no_expression_model": True, "no_marker_prior": True, "not_independent_single_cell_ground_truth": True, "domain": domain_audit, "affine": affine, "n_cells": len(ids), "n_positive_library_cells": int((all_counts > 0).sum()), "n_zero_library_cells_retained": int((all_counts == 0).sum()), "n_assigned_bins": len(assigned), "assigned_gene_conservation_max_error": 0, "assigned_count_mass": int(actual_gene.sum()), "unassigned_count_mass": int(matrix.sum()) - int(actual_gene.sum()), "finite_domain": {"local_neighbors": args.local_neighbors, "radius_scale": args.radius_scale, "radius_cap_quantile": args.radius_cap_quantile, "cap_um": cap * 2, "rule": "unique nearest-nucleus Euclidean Voronoi intersected with finite octagonal cap and native2 positive in-tissue domain"}, "inputs": {key: fingerprint(getattr(args, key)) for key in ("matrix", "positions", "scalefactors", "features")}, "source_math": {"path": str(args.pipeline_dir / "build_visium2_voronoi_cell_reference.py"), "sha256": digest(args.pipeline_dir / "build_visium2_voronoi_cell_reference.py")}, "outputs": {"counts": str(output), "assignment": str(assignment_path)}, "elapsed_seconds": time.time() - started}
    write_json(args.output_dir / f"{args.sample}.direct2um.qc.json", qc)
    print(json.dumps({"complete": True, "sample": args.sample, "n_cells": len(ids), "outputs": qc["outputs"]}, ensure_ascii=False, indent=2), flush=True)


def self_test(pipeline):
    core = load_core(pipeline)
    polygon = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]])
    assert abs(core.rectangle_intersection_area(polygon, 0, 0) - 1) < 1e-12
    assert core.rectangle_intersection_area(polygon, 2, 0) == 0
    grid = pd.DataFrame({"array_row": [0, 0, 1, 1], "array_col": [0, 1, 0, 1], "pxl_col_in_fullres": [10, 74, 10, 74], "pxl_row_in_fullres": [20, 20, 84, 84]})
    transform, intercept, audit = affine_from_positions(grid, .25, 16)
    np.testing.assert_allclose(transform, np.eye(2) * 64, atol=1e-10)
    np.testing.assert_allclose(intercept, [10, 20], atol=1e-10)
    probability = np.array([[1., 2., 7.]]) / 10
    np.testing.assert_allclose(probability[:, [0, 2]], [[.1, .7]])
    assert not np.isclose(probability[:, [0, 2]].sum(), 1)
    profile = core.marker_prior(np.array([.1, .2, .7]), np.array(["EPCAM", "PTPRC", "COL1A1"]), DEFAULT_MARKERS)
    np.testing.assert_allclose(profile.sum(1), 1, atol=1e-12)
    assignment, _, _, _ = finite_assignment(np.array([0, 0, 0]), np.array([0, 4, 20]), np.array([[0., 0.], [4., 0.]]), np.array([.5, .5]))
    assert np.array_equal(assignment, [0, 1, -1])
    # Disposable unit fixtures, never represented as biological outputs.
    with tempfile.TemporaryDirectory(prefix="generic_hd_unit_") as temporary:
        path = Path(temporary) / "duplicate_symbols.h5ad"
        fixture = ad.AnnData(shape=(2, 4), obs=pd.DataFrame({"cell_id": [7, 9]}, index=["a", "b"]), var=pd.DataFrame({"gene_name": ["A", "A", "B", "B"]}, index=["g1", "g2", "g3", "g4"]))
        write_h5ad_atomic(fixture, path)
        with h5py.File(path, "r") as h:
            assert isinstance(h["var/gene_name"], h5py.Dataset)
            assert decode(h["var/gene_name"][:]).tolist() == ["A", "A", "B", "B"]
            assert "X" not in h
        restored = ad.read_h5ad(path)
        assert restored.shape == (2, 4) and restored.X is None
        assert restored.var.gene_name.tolist() == ["A", "A", "B", "B"]
        # A measurement can be in the raw matrix but excluded by tissue mask.
        pospath = Path(temporary) / "positions.tsv.gz"
        tiny_positions = pd.DataFrame({"barcode": ["u", "v", "z", "w"], "in_tissue": [1, 1, 1, 0], "array_row": [0, 0, 1, 1], "array_col": [0, 1, 0, 1], "pxl_row_in_fullres": [0., 0., 64., 64.], "pxl_col_in_fullres": [0., 64., 0., 64.]})
        tiny_positions.to_csv(pospath, sep="\t", index=False, compression="gzip")
        _, owner, eligible, _, _, _ = core.align_positions(pospath, np.asarray(["u", "v", "z"]), np.asarray([[0., 0.], [64., 0.], [0., 64.], [64., 64.], [500., 500.]]))
        assert np.array_equal(eligible, [True, True, True, False, False])
        assert np.array_equal(owner[:2], [0, 1])
    print("SELF_TEST_PASS: original contour clipping, physical affine, full-gene normalization, five-class prior dimensions, finite-domain non-extrapolation, duplicate-symbol H5AD roundtrip, raw/tissue-domain exclusion", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=["metadata", "coarse", "direct2um", "test"], required=True)
    p.add_argument("--sample", default="TEST")
    p.add_argument("--matrix", type=Path)
    p.add_argument("--positions", type=Path)
    p.add_argument("--scalefactors", type=Path)
    p.add_argument("--features", type=Path)
    p.add_argument("--output-dir", type=Path)
    p.add_argument("--pipeline-dir", type=Path, default=Path(__file__).resolve().parent)
    p.add_argument("--marker-map", type=Path)
    p.add_argument("--panel-size", type=int, default=278)
    p.add_argument("--embedding-pcs", type=int, default=16)
    p.add_argument("--pca-sample-per-class", type=int, default=5000)
    p.add_argument("--class-prior-strength", type=float, default=.08)
    p.add_argument("--cv-folds", type=int, default=5)
    p.add_argument("--seed", type=int, default=20260828)
    p.add_argument("--local-neighbors", type=int, default=3)
    p.add_argument("--radius-scale", type=float, default=.75)
    p.add_argument("--radius-cap-quantile", type=float, default=.99)
    args = p.parse_args()
    if args.mode == "test":
        self_test(args.pipeline_dir); return
    for key in ("matrix", "positions", "scalefactors", "output_dir"):
        if getattr(args, key) is None:
            p.error(f"--{key.replace('_', '-')} is required")
    if args.mode == "coarse" and any("002um" in str(getattr(args, key)) or "direct2um" in str(getattr(args, key)) for key in ("matrix", "positions", "scalefactors", "features")):
        raise ValueError("High-resolution paths are forbidden in coarse preparation")
    if not 20 <= args.panel_size <= 1000 or not 0 < args.radius_cap_quantile <= 1:
        p.error("Invalid panel size or radius cap quantile")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    {"metadata": metadata, "coarse": coarse, "direct2um": direct}[args.mode](args)


if __name__ == "__main__":
    main()
