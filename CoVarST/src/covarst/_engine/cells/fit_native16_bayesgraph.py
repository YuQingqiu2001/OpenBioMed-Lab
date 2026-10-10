#!/usr/bin/env python3
"""Adapt frozen BayesGraph functions to native 16 um CRC bins.

Inputs: native 16 um count ledger/contour geometry; a completed cell H5AD's
        census, image features and global class profiles; fixed panel gene IDs.
Outputs: all-gene factor HDF5 plus provenance/preflight JSON in a NEW directory.
Dependencies: existing Ubuntu-20.04 ST_GJH environment and original pipeline.
Example (read-only metadata/64-column check, no fitting):
  python fit_native16_bayesgraph.py --sample P1 --preflight-only
Example (run only after preflight/throughput review):
  python fit_native16_bayesgraph.py --sample P1 --output-dir /new/run/P1

No original source is edited. Old X, allocation responsibilities, old Programs,
and calibration offsets are never read. This is a factor model, not a new
mass-conserving allocation. Inherited spatial CV uses globally estimated native
16 um class profiles and is therefore conditional/transductive model selection.
"""

from __future__ import annotations

import argparse
import copy
import gc
from functools import partial
import hashlib
import importlib
import inspect
import json
import os
from pathlib import Path
import re
import sys
import time

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components
import torch
import scipy

SCHEMA = "native16_BayesGraph_frozen_function_adapter_v1"
EPS = 1e-10
READ_LOG: list[dict] = []
ROOT = Path(__file__).resolve().parent
SOURCE = Path(__file__).resolve().parent


def decode(x):
    return np.asarray([v.decode("utf-8") if isinstance(v, bytes) else str(v)
                       for v in np.asarray(x)], dtype=str)


def safe_input(path: Path) -> Path:
    path = path.resolve()
    # Forbid high-resolution expression/reference paths even if passed explicitly.
    if re.search(r"(?:002um|(?<!\d)2um|Visium2|calibration)", str(path), re.I):
        raise ValueError(f"Forbidden high-resolution/calibration input: {path}")
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def read(h, key, selection=()):
    READ_LOG.append({"file": str(Path(h.filename).resolve()), "dataset": key,
                     "selection": str(selection)[:120]})
    return h[key][selection]


def read_obs(h, key):
    p = "obs/" + key
    if isinstance(h[p], h5py.Dataset):
        return read(h, p)
    return read(h, p + "/categories")[read(h, p + "/codes")]


def fingerprint(path):
    stat = path.stat()
    return {"path": str(path.resolve()), "bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def json_default(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def write_json(path, value):
    temporary = Path(str(path) + ".partial")
    if path.exists() or temporary.exists():
        raise FileExistsError(path)
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    default=json_default), encoding="utf-8")
    os.replace(temporary, path)


def install_runtime_compatibility(mv4):
    """Map the renamed SciPy CG keyword only; preserve its numeric tolerance.

    SciPy 1.10 exposes relative tolerance as ``tol`` while the original source
    calls the newer name ``rtol``. This process-local wrapper neither edits the
    original source nor changes the objective, solver, tolerances or maxiter.
    """
    original = mv4.splinalg.cg
    parameters = inspect.signature(original).parameters
    audit = {"scipy": scipy.__version__, "cg_signature_before": str(inspect.signature(original)),
             "cg_relative_tolerance_keyword_alias": "not_needed"}
    if "rtol" not in parameters and "tol" in parameters:
        def cg_compat(*args, rtol=None, **kwargs):
            if rtol is not None:
                if "tol" in kwargs:
                    raise TypeError("Specify either rtol or tol, not both")
                kwargs["tol"] = rtol
            return original(*args, **kwargs)
        mv4.splinalg.cg = cg_compat
        audit["cg_relative_tolerance_keyword_alias"] = "rtol_to_tol_same_numeric_value"
    elif "rtol" not in parameters:
        raise RuntimeError("Unsupported SciPy CG signature")
    system = sparse.csr_matrix([[2., .25], [.25, 3.]])
    rhs = np.asarray([1., 2.])
    solution, status = mv4.splinalg.cg(system, rhs, rtol=1e-7, atol=0., maxiter=800)
    if status != 0 or not np.allclose(system @ solution, rhs, rtol=1e-7, atol=0.):
        raise RuntimeError("CG compatibility smoke test failed")
    audit["cg_smoke_passed"] = True
    if audit["cg_relative_tolerance_keyword_alias"] != "not_needed":
        reference, reference_status = original(system, rhs, tol=1e-7, atol=0., maxiter=800)
        if reference_status != status or not np.array_equal(solution, reference):
            raise RuntimeError("CG keyword alias is not exactly equivalent")
        audit["cg_alias_direct_tol_bitwise_equal"] = True
    lsmr_parameters = inspect.signature(mv4.splinalg.lsmr).parameters
    if not all(key in lsmr_parameters for key in ["atol", "btol", "maxiter"]):
        raise RuntimeError("Unsupported LSMR signature")
    lsmr_result = mv4.splinalg.lsmr(system, rhs, atol=1e-9, btol=1e-9, maxiter=200)
    if not np.allclose(system @ lsmr_result[0], rhs, rtol=1e-8, atol=1e-8):
        raise RuntimeError("LSMR compatibility smoke failed")
    audit["lsmr_signature"] = str(inspect.signature(mv4.splinalg.lsmr))
    audit["lsmr_smoke_passed"] = True
    return audit


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True)
    p.add_argument("--class-schema", choices=["crc9", "input"], default="crc9",
                   help="input: class dimensions/names come only from input metadata, e.g. 5 PanNuke classes")
    p.add_argument("--pipeline-dir", type=Path, default=SOURCE)
    p.add_argument("--cell-h5ad", type=Path)
    p.add_argument("--allocation", type=Path)
    p.add_argument("--raw16", type=Path)
    p.add_argument("--positions", type=Path)
    p.add_argument("--scalefactors", type=Path)
    p.add_argument("--panel-gene-source", type=Path,
                   help="H5AD var identifiers only; X/obsm/state are never read")
    p.add_argument("--output-dir", type=Path)
    p.add_argument("--checkpoint-dir", type=Path,
                   help="Optional exact numerical cache; defaults to output-dir/numerical_checkpoints")
    p.add_argument("--preflight-only", action="store_true")
    p.add_argument("--preflight-json", type=Path)
    p.add_argument("--ranks", default="0,8")
    p.add_argument("--trend-lambdas", default="0.05,0.2")
    p.add_argument("--spatial-folds", type=int, default=5)
    p.add_argument("--spatial-buffer-um", type=float, default=150.0)
    p.add_argument("--coordinate-unit-um", type=float, default=100.0,
                   help="Original refine_cells radius 3.5 corresponds to 350 um at default")
    p.add_argument("--cv-steps-v55", type=int, default=100)
    p.add_argument("--final-steps-v55", type=int, default=320)
    p.add_argument("--fit-batch-size", type=int, default=384)
    p.add_argument("--trend-batch-size", type=int, default=4096)
    p.add_argument("--cpu-threads", type=int, default=8)
    p.add_argument("--max-cuda-gib", type=float, default=16.0)
    p.add_argument("--permutations", type=int, default=100)
    p.add_argument("--gene-block", type=int, default=512)
    p.add_argument("--seed", type=int, default=20260828)
    p.add_argument("--device", default="cuda")
    a = p.parse_args()
    s = a.sample
    if not re.fullmatch(r"[A-Za-z0-9_.+-]+", s):
        p.error("Sample ID must be a safe filename identifier")
    input_keys = ["cell_h5ad", "allocation", "raw16", "positions", "scalefactors", "panel_gene_source"]
    if s not in ["P1", "P2", "P5"] and any(getattr(a, k) is None for k in input_keys):
        p.error("For non-CRC samples supply all six input paths explicitly")
    native = ROOT / s / "binned_outputs/square_016um"
    defaults = {
        "cell_h5ad": ROOT / f"cell_resolved_016um/{s}/{s}.cell_resolved.h5ad",
        "allocation": ROOT / f"cell_resolved_016um/{s}/{s}.allocation.h5",
        "raw16": native / "filtered_feature_bc_matrix.h5",
        "positions": native / "spatial/tissue_positions.parquet",
        "scalefactors": native / "spatial/scalefactors_json.json",
        "panel_gene_source": ROOT / f"cell_resolved_spot2cell_mv3/{s}/visium55_v3/{s}.visium55.spot2cell_mv3.P_state_278.h5ad",
    }
    for k, v in defaults.items():
        setattr(a, k, safe_input(getattr(a, k) or v))
    a.spatial_buffer = a.spatial_buffer_um / a.coordinate_unit_um
    ranks = [int(x) for x in a.ranks.split(",")]
    if not ranks or any(x < 0 or x > 8 for x in ranks) or len(ranks) != len(set(ranks)):
        p.error("Ranks must be unique integers in [0,8]")
    if a.coordinate_unit_um <= 0 or a.spatial_buffer_um < 0 or a.spatial_folds < 2:
        p.error("Invalid physical scale, buffer, or fold count")
    if min(a.fit_batch_size, a.trend_batch_size, a.cpu_threads) < 1 or a.cpu_threads > 8:
        p.error("Batch sizes must be positive and CPU threads must be in [1,8]")
    if not 0 < a.max_cuda_gib <= 16:
        p.error("This adapter permits at most 16 GiB of additional Torch CUDA memory")
    if not a.preflight_only and a.output_dir is None:
        p.error("A fresh --output-dir is required for fitting")
    return a


def square_components(row, col):
    """4-neighbour components of the positive measured native square domain."""
    lookup = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(row, col))}
    if len(lookup) != len(row):
        raise ValueError("Duplicated native square coordinates")
    ii, jj = [], []
    for (r, c), i in lookup.items():
        for key in ((r + 1, c), (r, c + 1)):
            j = lookup.get(key)
            if j is not None:
                ii.extend((i, j)); jj.extend((j, i))
    adjacency = sparse.csr_matrix((np.ones(len(ii), np.uint8), (ii, jj)),
                                  shape=(len(row), len(row)))
    n, labels = connected_components(adjacency, directed=False)
    return labels.astype(np.int32), int(n)


def panel_genes(path):
    with h5py.File(path, "r") as h:
        if "var" not in h:
            raise ValueError("Panel source must have a fixed small var gene list, not an all-gene factor")
        key = "gene_name" if "gene_name" in h["var"] else h["var"].attrs.get("_index", "_index")
        if isinstance(key, bytes):
            key = key.decode()
        genes = decode(read(h, "var/" + key))
    if len(genes) > 1000 or len(genes) < 20 or len(set(genes)) != len(genes):
        raise ValueError("Expected a fixed unique panel (20-1000 gene identifiers)")
    return genes


def inspect_inputs(a):
    """Read metadata and small vectors; verify 64 complete raw-ledger columns."""
    with h5py.File(a.cell_h5ad, "r") as h:
        summary = json.loads(read(h, "uns/summary_json"))
        if (summary.get("heterogeneity_calibration") or {}).get("calibrated", False):
            raise ValueError("Use an uncalibrated cell source; calibrated sources are excluded")
        if "heterogeneity_calibration_reference" in summary.get("input_provenance", {}):
            raise ValueError("Cell source provenance includes a high-resolution calibration reference")
        key = "gene_name" if "gene_name" in h["var"] else h["var"].attrs.get("_index", "_index")
        if isinstance(key, bytes):
            key = key.decode()
        genes = decode(read(h, "var/" + key))
        cell_id = read_obs(h, "cell_id").astype(np.int64)
        owner = read_obs(h, "owner_bin_index").astype(np.int64)
        class_id = read_obs(h, "cellvit_class_id").astype(np.int32)
        cell_uv = read(h, "obsm/spatial_array_uv").astype(np.float64)
        feature_shape = list(h["obsm/standardized_nucleus_features"].shape)
        profiles = read(h, "uns/class_gene_probability")
        names = decode(read(h, "uns/class_names"))
        if not np.all(np.isfinite(profiles)) or np.any(profiles < 0):
            raise ValueError("Invalid global class profiles")
        if profiles.shape != (len(names), len(genes)):
            raise ValueError("Class profile shape mismatch")
        if len(names) != len(set(names)) or len(names) < 1 or class_id.min() < 0 or class_id.max() >= len(names):
            raise ValueError("Class names or numeric class IDs are inconsistent")
    fixed_panel = panel_genes(a.panel_gene_source)
    lookup = {g: i for i, g in enumerate(genes)}
    missing = [g for g in fixed_panel if g not in lookup]
    if missing:
        raise ValueError(f"Missing fixed panel genes: {missing[:12]}")
    scalefactors = json.loads(a.scalefactors.read_text(encoding="utf-8"))
    READ_LOG.append({"file": str(a.scalefactors), "dataset": "JSON metadata"})
    if abs(float(scalefactors["bin_size_um"]) - 16.0) > 1e-8:
        raise ValueError("Only native 16 um bins are supported")
    with h5py.File(a.allocation, "r") as h, h5py.File(a.raw16, "r") as raw:
        if int(h.attrs.get("complete", 0)) != 1:
            raise ValueError("Incomplete allocation ledger")
        shape = tuple(int(x) for x in read(h, "observed/shape"))
        geometry_shape = tuple(int(x) for x in read(h, "geometry/shape"))
        barcodes = decode(read(h, "observed/barcode"))
        source_bin = read(h, "observed/positive_source_bin_index").astype(np.int64)
        raw_barcode = decode(read(raw, "matrix/barcodes"))
        if not np.array_equal(raw_barcode[source_bin], barcodes):
            raise ValueError("Native raw matrix order disagrees with ledger")
        if not np.array_equal(decode(read(raw, "matrix/features/name")), genes):
            raise ValueError("Native gene names/order disagree with cell source")
        if not np.array_equal(decode(read(h, "observed/gene_name")), genes):
            raise ValueError("Ledger genes/order disagree with cell source")
        if geometry_shape != (shape[1], len(cell_id)):
            raise ValueError("Geometry and census dimensions disagree")
        if len(set(cell_id)) != len(cell_id) or owner.min() < 0 or owner.max() >= shape[1]:
            raise ValueError("Invalid cell IDs or owner bins")
        old_ptr = read(h, "observed/bin_indptr")
        raw_ptr = read(raw, "matrix/indptr")
        sampled = np.unique(np.linspace(0, shape[1] - 1, min(64, shape[1])).astype(int))
        for b in sampled:
            t = source_bin[b]
            aa, ab = int(old_ptr[b]), int(old_ptr[b + 1])
            ra, rb = int(raw_ptr[t]), int(raw_ptr[t + 1])
            if not np.array_equal(read(h, "observed/gene_index", slice(aa, ab)),
                                  read(raw, "matrix/indices", slice(ra, rb))):
                raise ValueError(f"Gene-index ledger mismatch in bin {b}")
            if not np.array_equal(read(h, "observed/count", slice(aa, ab)),
                                  read(raw, "matrix/data", slice(ra, rb))):
                raise ValueError(f"Count ledger mismatch in bin {b}")
        observed_nnz = int(h["observed/count"].size)
    positions = pd.read_parquet(a.positions, columns=["barcode", "in_tissue", "array_row", "array_col",
                                                       "pxl_row_in_fullres", "pxl_col_in_fullres"])
    READ_LOG.append({"file": str(a.positions), "dataset": "positions-only parquet"})
    if positions["barcode"].duplicated().any():
        raise ValueError("Duplicate positions barcodes")
    positions = positions.set_index("barcode").loc[barcodes].reset_index()
    if not positions["in_tissue"].eq(1).all():
        raise ValueError("Ledger includes non-tissue bins")
    owner_xy = positions[["array_col", "array_row"]].to_numpy(float)[owner]
    owner_offset = float(np.max(np.abs(cell_uv - owner_xy)))
    if not np.isfinite(cell_uv).all() or owner_offset > 0.501:
        raise ValueError(f"Cell coordinates disagree with native owner-square coordinates: {owner_offset}")
    component, n_components = square_components(positions.array_row, positions.array_col)
    # Compatibility names only: these are scaled square coordinates, never a hex lattice.
    positions["hex_col"] = positions.array_col.to_numpy(float) * 16.0 / a.coordinate_unit_um
    positions["hex_row"] = positions.array_row.to_numpy(float) * 16.0 / a.coordinate_unit_um
    report = {
        "schema": SCHEMA, "sample": a.sample, "preflight_passed": True,
        "cell_source_schema": summary.get("schema"), "n_cells": len(cell_id),
        "n_positive_native16_bins": shape[1], "n_genes": len(genes),
        "n_panel_genes": len(fixed_panel), "class_names": names.tolist(),
        "nucleus_features_shape": feature_shape, "ledger_nnz": observed_nnz,
        "max_abs_cell_to_owner_square_coordinate_offset": owner_offset,
        "raw_ledger_check": {"all_barcodes_and_genes_match": True,
                             "complete_columns_exactly_checked": len(sampled),
                             "all_nonzero_values_checked": False},
        "square_4nn_components": n_components, "coordinate_unit_um": a.coordinate_unit_um,
        "cv_buffer_um": a.spatial_buffer_um, "cell_graph_radius_um": 3.5 * a.coordinate_unit_um,
        "global_profiles_estimated_from_all_native16_bins": True,
        "CV_interpretation": "conditional/transductive native16 model selection; not independent profile validation",
        "2um_expression_used_for_fit": False, "old_cell_expression_read": False,
        "old_programs_or_offsets_read": False, "calibration_source_read": False,
        "source_calibration": summary.get("heterogeneity_calibration"),
        "profile_provenance": summary.get("input_provenance", {}),
        "panel_usage": "gene identifiers only; no panel X, expression, coordinates, or states read",
        "approx_bin_draws_per_positive_bin": {
            "per_cv_fit": a.fit_batch_size * a.cv_steps_v55 / shape[1],
            "per_final_fit": a.fit_batch_size * a.final_steps_v55 / shape[1]},
        "inputs": {k: fingerprint(getattr(a, k)) for k in
                   ["cell_h5ad", "allocation", "raw16", "positions", "scalefactors", "panel_gene_source"]},
        "environment": {"python": sys.version, "numpy": np.__version__, "torch": torch.__version__,
                        "h5py": h5py.__version__, "cuda_available": torch.cuda.is_available()},
    }
    small = {"gene_names": genes, "panel_genes": fixed_panel,
             "panel_index": np.asarray([lookup[g] for g in fixed_panel], np.int32),
             "cell_id": cell_id, "owner_bin": owner.astype(np.int32), "class_id": class_id,
             "profiles": profiles.astype(np.float64), "spots": positions,
             "spot_component": component, "cell_component": component[owner]}
    return report, small


def load_large_inputs(a, data):
    with h5py.File(a.cell_h5ad, "r") as h:
        data["capacity"] = read_obs(h, "rna_capacity").astype(np.float64)
        data["coords"] = read(h, "obsm/spatial_array_uv").astype(np.float64) * 16.0 / a.coordinate_unit_um
        data["spatial_px"] = read(h, "obsm/spatial").astype(np.float64)
        data["features"] = np.clip(np.nan_to_num(read(h, "obsm/standardized_nucleus_features")), -6., 6.).astype(np.float64)
    if not np.isfinite(data["coords"]).all() or not np.isfinite(data["capacity"]).all() or np.any(data["capacity"] <= 0):
        raise ValueError("Invalid census geometry or capacity")
    with h5py.File(a.allocation, "r") as h:
        data["geometry"] = sparse.csr_matrix((read(h, "geometry/overlap_fraction").astype(np.float64),
                                               read(h, "geometry/cell_index").astype(np.int32),
                                               read(h, "geometry/indptr").astype(np.int64)),
                                              shape=tuple(read(h, "geometry/shape")))
        counts = read(h, "observed/count")
        if not np.isfinite(counts).all() or np.any(counts < 0) or np.any(counts != np.floor(counts)):
            raise ValueError("Native observed ledger must contain finite integer UMI counts")
        data["matrix"] = sparse.csc_matrix((counts.astype(np.float32),
                                             read(h, "observed/gene_index").astype(np.int32),
                                             read(h, "observed/bin_indptr").astype(np.int64)),
                                            shape=tuple(read(h, "observed/shape")))
    if np.any(np.asarray(data["matrix"].sum(axis=0)).ravel() <= 0):
        raise ValueError("Ledger contains a zero-count bin")
    data["profiles"] /= np.maximum(data["profiles"].sum(1, keepdims=True), EPS)
    mass = np.column_stack([np.asarray(data["geometry"] @ (data["capacity"] * (data["class_id"] == c))).ravel()
                            for c in range(len(data["profiles"]))])
    total = mass.sum(1)
    composition = np.divide(mass, total[:, None], out=np.zeros_like(mass), where=total[:, None] > 0).astype(np.float32)
    panel_counts = data["matrix"][data["panel_index"]].T.toarray().astype(np.float32)
    panel_profiles = np.maximum(data["profiles"][:, data["panel_index"]], EPS)
    panel_profiles /= panel_profiles.sum(1, keepdims=True)
    return composition, total, panel_counts, panel_profiles.astype(np.float32)


def main():
    a = parse_args()
    torch.set_num_threads(a.cpu_threads)
    torch.set_num_interop_threads(1)
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:a.cpu_threads])
    if not a.preflight_only and a.device.startswith("cuda"):
        torch.cuda.set_per_process_memory_fraction(
            a.max_cuda_gib * (1 << 30) / torch.cuda.get_device_properties(0).total_memory)
    sys.path.insert(0, str(a.pipeline_dir))
    dual = importlib.import_module("p1_dual_bayesgraph_panel")
    sweep = importlib.import_module("p1_bayesgraph_fullgene_rank_sweep")
    bg = importlib.import_module("spot2cell_bayesgraph_core")
    runtime_compatibility = install_runtime_compatibility(dual.mv4)
    # Explicit parameter forwarding to the original function, with no source or
    # objective changes. dual.bg is this same module object.
    bg.fit_spot_program = partial(bg.fit_spot_program,
                                 batch_size=a.fit_batch_size,
                                 trend_batch_size=a.trend_batch_size)
    report, data = inspect_inputs(a)
    report["runtime_compatibility"] = runtime_compatibility
    if a.class_schema == "input":
        # Only category metadata and the resulting class dimension change.
        # All called fitting/inverse functions already operate on this dimension.
        names = tuple(report["class_names"])
        dual.CLASS_NAMES = names
        dual.mv3.CLASS_NAMES = names
        dual.proto.CLASS_NAMES = names
        sweep.mv3.CLASS_NAMES = names
        report["class_schema"] = "input_metadata_dynamic"
        report["CRC_Lizard_labels_or_organ_markers_imported_by_adapter"] = False
    elif list(dual.CLASS_NAMES) != report["class_names"]:
        raise ValueError("Original BayesGraph class order does not match native16 census")
    else:
        report["class_schema"] = "crc9_fixed_input"
    report["original_sources_sha256"] = {name: sha256(a.pipeline_dir / name) for name in
        ["p1_dual_bayesgraph_panel.py", "p1_bayesgraph_fullgene_rank_sweep.py",
         "p1_visium55_mv4_panel_prototype.py", "cell_resolved_visium55_mv3.py",
         "spot2cell_bayesgraph_core.py", "spot2cell_mv4_core.py"]}
    report["adapter_sha256"] = sha256(Path(__file__))
    report["read_audit"] = READ_LOG.copy()
    if a.preflight_json:
        write_json(a.preflight_json, report)
    if a.preflight_only:
        print(json.dumps({k: v for k, v in report.items() if k != "read_audit"},
                         ensure_ascii=False, indent=2, default=json_default))
        return
    if a.output_dir.exists():
        raise FileExistsError(f"Output must be a new directory: {a.output_dir}")
    a.output_dir.mkdir(parents=True)
    import bayesgraph_runtime
    bayesgraph_runtime.install_checkpoints(dual.mv4,
        a.checkpoint_dir or a.output_dir / "numerical_checkpoints", report["original_sources_sha256"])
    report["checkpoint_helper_sha256"] = sha256(Path(bayesgraph_runtime.__file__))
    write_json(a.output_dir / f"{a.sample}.native16.preflight.json", report)
    print(f"[{a.sample}] Loading raw native16 ledger and image-only nucleus features", flush=True)
    composition, total, panel_counts, panel_profiles = load_large_inputs(a, data)
    coords = data["spots"][["hex_col", "hex_row"]].to_numpy(np.float64)
    context = np.column_stack((np.sqrt(np.maximum(composition, 0)), np.log1p(total)))
    graph = bg.build_spatial_graph(coords, component=data["spot_component"], context=context,
                                   neighbors=6, radius_factor=2.5)
    started = time.time()
    results = {}
    for rank in dual.parse_numbers(a.ranks, int):
        stage_started = time.time()
        local = copy.copy(a)
        local.ranks = str(rank)
        local.force_rank_v55 = rank
        final, _, cv, selected, folds = dual.fit_v55(local, panel_counts, composition,
            panel_profiles, coords, data["spot_component"], graph)
        fit_finished = time.time()
        print(f"[{a.sample}] rank{rank} CV/final complete; starting cell refinement", flush=True)
        state, contrast, morph, edges, audits, program_count, support = dual.refine_cells(
            local, data, data["panel_index"], composition, final, panel_counts.sum(axis=1))
        refine_finished = time.time()
        print(f"[{a.sample}] rank{rank} refinement complete; extending all genes", flush=True)
        loadings = sweep.extend_v55_loadings(data, composition, final, data["panel_index"], a.gene_block)
        loadings[program_count == 0] = 0
        summary = {k: v for k, v in report.items() if k != "read_audit"}
        summary.update({"rank": rank, "selected": selected, "cv": cv, "spot_graph": graph.audit,
            "program_count_by_class": program_count, "support_spots_by_class": support,
            "class_refinement": audits, "parameters": vars(a),
            "factor_formula": "softmax(log(class_gene_probability)+cell_state@program_loadings)",
            "factor_is_count_conserving_allocation": False,
            "fixed_44_programs_used": False, "external_reference_used": False,
            "stage_seconds": {"CV_and_final_fit": fit_finished - stage_started,
                              "cell_refinement": refine_finished - fit_finished,
                              "all_gene_extension": time.time() - refine_finished}})
        output = a.output_dir / f"{a.sample}.16um.BayesGraph_rank{rank}_all_gene.factorized.h5"
        sweep.write_factor(output, modality=f"native16_BG_rank{rank}", genes=data["gene_names"],
            cell_id=data["cell_id"], class_id=data["class_id"], spatial=data["spatial_px"],
            profiles=data["profiles"], loadings=loadings, state=state, program_count=program_count,
            summary=summary, sample=a.sample, extra={"spot_program_activity": final.w.astype(np.float32),
                "graph_contrast": contrast, "morph_contrast": morph, "spatial_fold": folds.astype(np.int8),
                "cell_graph_edges": edges.astype(np.float32), "panel_gene_index": data["panel_index"]})
        results[str(rank)] = {"path": str(output), "selected": selected}
        print(f"[{a.sample}] Completed rank{rank}: {output}", flush=True)
        if final.model is not None:
            final.model.to("cpu")
        del final, state, contrast, morph, loadings
        torch.cuda.empty_cache(); gc.collect()
    write_json(a.output_dir / f"{a.sample}.native16.completed.json",
        {"schema": SCHEMA, "status": "complete", "results": results,
         "elapsed_seconds": time.time() - started, "parameters": vars(a), "read_audit": READ_LOG,
         "2um_expression_used_for_fit": False, "old_cell_expression_read": False})


if __name__ == "__main__":
    main()
