#!/usr/bin/env python3
"""Build sample-specific full-gene BayesGraph factors for 2 um and V55 rank0-rank8.

The expensive rank sweep is stored as an exact all-gene factorization rather
than nine duplicated dense cell-by-gene matrices.  Downstream analysis
regenerates log1p-CP10K expression in cell blocks.  V55 never reads 2 um
expression during fitting.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import time

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import torch

import cell_resolved_visium55_mv3 as mv3
import p1_dual_bayesgraph_panel as dual
import spot2cell_bayesgraph_core as bg


EPS = 1.0e-10
SCHEMA = "BayesGraph_fullgene_rank_sweep_v1"


def log(message: str) -> None:
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), message, flush=True)


def parse_args() -> argparse.Namespace:
    root = Path('DATA/Visium_HD')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", choices=("P1", "P2", "P5"), default="P1")
    parser.add_argument(
        "--reference-2um-panel", type=Path,
        default=None,
    )
    parser.add_argument(
        "--assignment-2um", type=Path,
        default=None,
    )
    parser.add_argument(
        "--matrix-2um", type=Path,
        default=None,
    )
    parser.add_argument("--mv2-h5ad", type=Path, default=None)
    parser.add_argument("--mv2-allocation", type=Path, default=None)
    parser.add_argument("--spot-positions", type=Path, default=None)
    parser.add_argument("--mv3-state", type=Path, default=None)
    parser.add_argument("--panel-template", type=Path, default=None)
    parser.add_argument(
        "--output-dir", type=Path,
        default=None,
    )
    parser.add_argument("--ranks", default="0,1,2,3,4,5,6,7,8")
    parser.add_argument("--trend-lambdas", default="0.05,0.2")
    parser.add_argument("--tau-grid", default="0,2,5,10,20,50")
    parser.add_argument("--spatial-folds", type=int, default=5)
    parser.add_argument("--spatial-buffer", type=float, default=1.5)
    parser.add_argument("--cv-steps-v55", type=int, default=100)
    parser.add_argument("--final-steps-v55", type=int, default=320)
    parser.add_argument("--cv-steps-2um", type=int, default=160)
    parser.add_argument("--final-steps-2um", type=int, default=360)
    parser.add_argument("--permutations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260828)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--gene-block", type=int, default=512)
    parser.add_argument("--skip-2um", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    sample = args.sample
    mv2 = root / f"cell_resolved_spot2cell_mv2/{sample}/visium55_v2"
    mv3_root = root / f"cell_resolved_spot2cell_mv3/{sample}/visium55_v3"
    defaults = {
        "reference_2um_panel": root / f"cell_resolved_002um_voronoi/{sample}/{sample}.Visium2_direct.full_support_input.h5ad",
        "assignment_2um": root / f"cell_resolved_002um_voronoi/{sample}/{sample}.Visium2_direct.bin_assignment.h5",
        "matrix_2um": Path('DATA/Visium_HD'),
        "mv2_h5ad": mv2 / f"{sample}.visium55.spot2cell_mv2.h5ad",
        "mv2_allocation": mv2 / f"{sample}.visium55.spot2cell_mv2.allocation.h5",
        "spot_positions": mv2 / f"{sample}.visium55.spots.tsv.gz",
        "mv3_state": mv3_root / f"{sample}.visium55.spot2cell_mv3.state.h5",
        "panel_template": mv3_root / f"{sample}.visium55.spot2cell_mv3.P_state_278.h5ad",
        "output_dir": root / f"cell_resolved_spot2cell_rank_sweep_fullgene_v1/{sample}",
    }
    for name, value in defaults.items():
        if getattr(args, name) is None:
            setattr(args, name, value)
    return args


def decode(values) -> np.ndarray:
    values = np.asarray(values)
    if values.dtype.kind == "S":
        return values.astype("U")
    return np.asarray([
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ], dtype="U")


def sha256(path: Path, block: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(block):
            digest.update(chunk)
    return digest.hexdigest()


def reserve(path: Path) -> Path:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    partial = Path(str(path) + ".partial")
    if partial.exists():
        raise FileExistsError(f"Preserved partial exists: {partial}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return partial


def atomic_json(value: dict, path: Path) -> None:
    partial = reserve(path)
    partial.write_text(json.dumps(value, ensure_ascii=False, indent=2, default=json_default), encoding="utf-8")
    os.replace(partial, path)


def json_default(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    raise TypeError(type(value).__name__)


def atomic_h5ad(data: ad.AnnData, path: Path) -> None:
    partial = reserve(path)
    data.write_h5ad(partial, compression="lzf")
    os.replace(partial, path)


def write_text(group: h5py.Group, name: str, values) -> None:
    group.create_dataset(name, data=np.asarray(values, dtype=h5py.string_dtype("utf-8")))


def aggregate_full_2um(args: argparse.Namespace, panel: ad.AnnData) -> tuple[sparse.csr_matrix, np.ndarray, Path]:
    output = args.output_dir / f"{args.sample}.2um.direct_all_gene_counts.h5ad"
    if output.exists():
        log(f"Reusing complete all-gene 2 um aggregation: {output}")
        data = ad.read_h5ad(output)
        symbols = (
            data.var["gene_symbol"].astype(str).to_numpy()
            if "gene_symbol" in data.var else data.var_names.astype(str).to_numpy()
        )
        return sparse.csr_matrix(data.X, dtype=np.int32), symbols, output

    log("Reading full filtered 2 um matrix and frozen bin-to-cell assignment")
    with h5py.File(args.matrix_2um, "r") as handle:
        group = handle["matrix"]
        shape = tuple(int(v) for v in group["shape"][:])
        matrix = sparse.csc_matrix(
            (
                group["data"][:].astype(np.int32),
                group["indices"][:].astype(np.int32),
                group["indptr"][:].astype(np.int64),
            ),
            shape=shape,
        )
        genes = decode(group["features/name"][:])
        gene_ids = decode(group["features/id"][:])
    with h5py.File(args.assignment_2um, "r") as handle:
        bin_cell = handle["bin_cell_index"][:].astype(np.int32)
        cell_id = handle["cell_id"][:].astype(np.int64)
    if matrix.shape[1] != len(bin_cell):
        raise RuntimeError(f"2 um bin order mismatch: {matrix.shape[1]} != {len(bin_cell)}")
    if len(cell_id) != panel.n_obs:
        raise RuntimeError("2 um assignment and panel cell count differ")
    valid = np.flatnonzero(bin_cell >= 0)
    assignment = sparse.csc_matrix(
        (np.ones(len(valid), np.int8), (valid, bin_cell[valid])),
        shape=(matrix.shape[1], panel.n_obs),
    )
    log(f"Aggregating {len(valid):,} assigned 2 um bins across {len(genes):,} genes")
    expected_total = int(matrix[:, valid].sum())
    counts = (matrix @ assignment).T.tocsr().astype(np.int32)
    counts.sum_duplicates(); counts.sort_indices(); counts.eliminate_zeros()
    del matrix, assignment
    if int(counts.sum()) != expected_total:
        raise RuntimeError(
            f"Assigned 2 um UMI conservation failure: {int(counts.sum())} != {expected_total}"
        )
    obs = panel.obs.copy()
    if len(set(gene_ids)) != len(gene_ids):
        raise RuntimeError("2 um Ensembl feature IDs are not unique")
    var = pd.DataFrame(
        {"gene_symbol": genes}, index=pd.Index(gene_ids, name="gene_id")
    )
    result = ad.AnnData(X=counts, obs=obs, var=var)
    for key in panel.obsm:
        result.obsm[key] = np.asarray(panel.obsm[key]).copy()
    result.uns["schema"] = SCHEMA
    result.uns["X_semantics"] = "integer counts from frozen finite cell domains over measured positive 2 um bins"
    result.uns["not_independent_single_cell_ground_truth"] = True
    result.uns["source_matrix_sha256"] = sha256(args.matrix_2um)
    result.uns["assignment_sha256"] = sha256(args.assignment_2um)
    atomic_h5ad(result, output)
    return counts, genes, output


def extend_direct_loadings(
    counts: sparse.csr_matrix,
    class_id: np.ndarray,
    profiles: np.ndarray,
    state: np.ndarray,
    rank: int,
    block_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    centered = np.zeros_like(state, dtype=np.float32)
    loadings = np.zeros((len(mv3.CLASS_NAMES), rank, counts.shape[1]), np.float32)
    if rank == 0:
        return centered, loadings
    library = np.asarray(counts.sum(1)).ravel().astype(np.float64)
    for c, name in enumerate(mv3.CLASS_NAMES):
        rows = np.flatnonzero((class_id == c) & (library > 0))
        if len(rows) < max(40, rank + 2):
            continue
        z = state[rows, :rank].astype(np.float64)
        weight = np.sqrt(np.maximum(library[rows], 1.0))
        mean = np.average(z, axis=0, weights=weight)
        z -= mean
        centered[rows, :rank] = z.astype(np.float32)
        gram = z.T @ (z * weight[:, None])
        ridge = max(0.02 * float(np.median(np.diag(gram))), 1.0e-6)
        projector = np.linalg.solve(gram + ridge * np.eye(rank), z.T * weight[None])
        base = profiles[c].astype(np.float64)
        for start in range(0, counts.shape[1], block_size):
            end = min(start + block_size, counts.shape[1])
            observed = counts[rows, start:end].toarray().astype(np.float64)
            observed *= np.divide(1.0e4, library[rows], out=np.zeros(len(rows)), where=library[rows] > 0)[:, None]
            target = np.log1p(observed) - np.log1p(1.0e4 * base[start:end])[None]
            target -= np.average(target, axis=0, weights=weight)
            loadings[c, :, start:end] = np.clip(projector @ target, -4.0, 4.0).astype(np.float32)
        log(f"2 um {name}: extended rank {rank} to {counts.shape[1]:,} genes")
    return centered, loadings


def extend_v55_loadings(
    data: dict,
    composition: np.ndarray,
    fit: bg.ProgramFit,
    panel_index: np.ndarray,
    block_size: int,
) -> np.ndarray:
    rank = fit.v.shape[1]
    profiles = np.asarray(data["profiles"], np.float64)
    profiles /= np.maximum(profiles.sum(1, keepdims=True), EPS)
    loadings = np.zeros((len(mv3.CLASS_NAMES), rank, len(data["gene_names"])), np.float32)
    if rank == 0:
        return loadings
    matrix = data["matrix"]
    library = np.asarray(matrix.sum(0)).ravel().astype(np.float64)
    for c, name in enumerate(mv3.CLASS_NAMES):
        bins = np.flatnonzero((composition[:, c] > 1.0e-4) & (library > 0))
        if len(bins) < max(40, rank + 2):
            continue
        fraction = composition[bins, c].astype(np.float64)
        z = fit.w[bins, c, :rank].astype(np.float64)
        weight = np.maximum(fraction * np.sqrt(library[bins]), 0.05)
        gram = z.T @ (z * weight[:, None])
        ridge = max(0.02 * float(np.median(np.diag(gram))), 1.0e-6)
        projector = np.linalg.solve(gram + ridge * np.eye(rank), z.T * weight[None])
        for start in range(0, matrix.shape[0], block_size):
            end = min(start + block_size, matrix.shape[0])
            observed = matrix[start:end, bins].T.toarray().astype(np.float64)
            observed *= np.divide(1.0, library[bins], out=np.zeros(len(bins)), where=library[bins] > 0)[:, None]
            profile_block = profiles[:, start:end]
            mixture = composition[bins].astype(np.float64) @ profile_block
            signal = observed - (mixture - fraction[:, None] * profile_block[c][None])
            signal = np.maximum(signal, 0.0) / np.maximum(fraction[:, None], 1.0e-4)
            target = np.log1p(1.0e4 * signal) - np.log1p(1.0e4 * profile_block[c])[None]
            target = np.clip(target, -6.0, 6.0)
            target -= np.average(target, axis=0, weights=weight)
            loadings[c, :, start:end] = np.clip(projector @ target, -4.0, 4.0).astype(np.float32)
        loadings[c][:, panel_index] = fit.v[c, :rank]
        log(f"V55 {name}: extended rank {rank} to {matrix.shape[0]:,} genes")
    return loadings


def write_factor(
    path: Path,
    *,
    modality: str,
    genes: np.ndarray,
    cell_id: np.ndarray,
    class_id: np.ndarray,
    spatial: np.ndarray,
    profiles: np.ndarray,
    loadings: np.ndarray,
    state: np.ndarray,
    program_count: np.ndarray,
    summary: dict,
    sample: str,
    extra: dict[str, np.ndarray] | None = None,
) -> None:
    partial = reserve(path)
    with h5py.File(partial, "w") as handle:
        handle.attrs.update({
            "schema": SCHEMA,
            "complete": 0,
            "sample": sample,
            "modality": modality,
            "expression_semantics": "all-gene factorized inferred cell probability",
            "summary_json": json.dumps(summary, ensure_ascii=False, default=json_default),
        })
        write_text(handle, "gene_name", genes)
        cells = handle.create_group("cells")
        cells.create_dataset("cell_id", data=cell_id, compression="lzf")
        cells.create_dataset("class_id", data=class_id.astype(np.uint8), compression="lzf")
        cells.create_dataset("spatial", data=spatial.astype(np.float32), compression="lzf")
        model = handle.create_group("model")
        model.create_dataset("class_gene_probability", data=profiles.astype(np.float32), compression="lzf")
        model.create_dataset("program_loadings", data=loadings.astype(np.float32), compression="lzf")
        model.create_dataset("program_count_by_class", data=program_count.astype(np.uint8))
        model.create_dataset("cell_program_activity", data=state.astype(np.float32), compression="lzf")
        if extra:
            for name, values in extra.items():
                model.create_dataset(name, data=values, compression="lzf")
        handle.attrs["complete"] = 1
        handle.flush()
    os.replace(partial, path)


def run_2um(args: argparse.Namespace) -> dict:
    output = args.output_dir / f"{args.sample}.2um.BayesGraph_corrected_all_gene.factorized.h5"
    if args.resume and output.exists():
        with h5py.File(output, "r") as handle:
            if int(handle.attrs.get("complete", 0)) != 1:
                raise RuntimeError(f"Incomplete resume target: {output}")
            saved = json.loads(handle.attrs["summary_json"])
        log(f"Reusing completed 2 um factor: {output}")
        return {"factor": str(output), "reused": True, **saved}
    panel = ad.read_h5ad(args.reference_2um_panel)
    counts, genes, counts_path = aggregate_full_2um(args, panel)
    class_id = panel.obs["cellvit_class_id"].to_numpy(np.int64)
    coords = np.asarray(panel.obsm["spatial"], np.float64)
    panel_counts = sparse.csr_matrix(panel.layers["measured_2um_sum_counts"], dtype=np.int32)
    graph = bg.build_spatial_graph(coords, group=class_id, neighbors=6, radius_factor=4.0)
    local = copy.copy(args)
    final, _, cv, selected = dual.fit_2um(local, panel_counts, coords, class_id, graph)
    profiles = bg.class_profiles(counts, class_id, len(mv3.CLASS_NAMES)).astype(np.float32)
    state, loadings = extend_direct_loadings(
        counts, class_id, profiles, final.w, int(selected["rank"]), args.gene_block,
    )
    summary = {
        "selected": selected,
        "cv": cv,
        "graph": graph.audit,
        "n_cells": int(len(class_id)),
        "n_genes": int(len(genes)),
        "counts_path": str(counts_path),
        "not_independent_single_cell_ground_truth": True,
        "external_reference_used": False,
    }
    write_factor(
        output, modality="2um_BG_corrected", genes=genes,
        cell_id=panel.obs["cell_id"].to_numpy(np.int64), class_id=class_id,
        spatial=coords, profiles=profiles, loadings=loadings, state=state,
        program_count=np.full(len(mv3.CLASS_NAMES), int(selected["rank"]), np.int32),
        summary=summary, sample=args.sample,
        extra={"posterior_tau": np.asarray([float(selected["selected_tau"])], np.float32)},
    )
    if final.model is not None:
        final.model.to("cpu")
    torch.cuda.empty_cache()
    return {"factor": str(output), **summary}


def run_v55(args: argparse.Namespace) -> dict:
    data, _, panel_genes, panel_index, panel_counts, panel_profiles, composition, type_mass = dual.load_v55(args)
    coords_spot = data["spots"][["hex_col", "hex_row"]].to_numpy(np.float64)
    graph_context = np.column_stack((np.sqrt(np.maximum(composition, 0.0)), np.log1p(type_mass)))
    graph = bg.build_spatial_graph(
        coords_spot, component=data["spot_component"], context=graph_context,
        neighbors=6, radius_factor=2.5,
    )
    _, _, obsm = dual.v55_obs_var(args.panel_template)
    spatial_cell = np.asarray(obsm["spatial"], np.float64)
    full_profiles = np.asarray(data["profiles"], np.float64)
    full_profiles /= np.maximum(full_profiles.sum(1, keepdims=True), EPS)
    results = {}
    for rank in dual.parse_numbers(args.ranks, int):
        output = args.output_dir / f"{args.sample}.V55.BayesGraph_rank{rank}_all_gene.factorized.h5"
        if args.resume and output.exists():
            with h5py.File(output, "r") as handle:
                if int(handle.attrs.get("complete", 0)) != 1:
                    raise RuntimeError(f"Incomplete resume target: {output}")
                saved = json.loads(handle.attrs["summary_json"])
            log(f"Reusing completed V55 rank {rank}: {output}")
            results[str(rank)] = {"factor": str(output), "reused": True, **saved}
            continue
        local = copy.copy(args)
        local.ranks = str(rank)
        local.force_rank_v55 = int(rank)
        log(f"Starting V55 full-gene rank {rank}")
        final, _, cv, selected, folds = dual.fit_v55(
            local, panel_counts, composition, panel_profiles, coords_spot,
            data["spot_component"], graph,
        )
        state, graph_contrast, morph_contrast, cell_edges, class_audits, program_count, support = dual.refine_cells(
            local, data, panel_index, composition, final, panel_counts.sum(axis=1),
        )
        full_loadings = extend_v55_loadings(data, composition, final, panel_index, args.gene_block)
        full_loadings[program_count == 0] = 0
        summary = {
            "rank": int(rank),
            "selected": selected,
            "cv": cv,
            "spot_graph": graph.audit,
            "program_count_by_class": program_count,
            "support_spots_by_class": support,
            "class_refinement": class_audits,
            "n_cells": int(len(data["class_id"])),
            "n_spots": int(len(panel_counts)),
            "n_genes": int(len(data["gene_names"])),
            "fixed_44_programs_used": False,
            "2um_expression_used_for_fit": False,
            "external_reference_used": False,
            "factor_formula": "softmax(log(class_gene_probability)+cell_state@program_loadings)",
        }
        write_factor(
            output, modality=f"V55_BG_rank{rank}", genes=np.asarray(data["gene_names"], str),
            cell_id=data["cell_id"], class_id=data["class_id"], spatial=spatial_cell,
            profiles=full_profiles, loadings=full_loadings, state=state,
            program_count=program_count, summary=summary, sample=args.sample,
            extra={
                "spot_program_activity": final.w.astype(np.float32),
                "graph_contrast": graph_contrast.astype(np.float32),
                "morph_contrast": morph_contrast.astype(np.float32),
                "spatial_fold": folds.astype(np.int8),
                "cell_graph_edges": cell_edges.astype(np.int32),
            },
        )
        results[str(rank)] = {"factor": str(output), **summary}
        if final.model is not None:
            final.model.to("cpu")
        del final, state, graph_contrast, morph_contrast, full_loadings
        torch.cuda.empty_cache()
    return results


def main() -> None:
    args = parse_args()
    required = [
        args.reference_2um_panel, args.assignment_2um, args.matrix_2um,
        args.mv2_h5ad, args.mv2_allocation, args.spot_positions,
        args.mv3_state, args.panel_template,
    ]
    for path in required:
        if not path.exists():
            raise FileNotFoundError(path)
    if args.output_dir.exists() and not args.resume:
        raise FileExistsError(f"Refusing to overwrite output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=args.resume)
    started = time.time()
    summary = {
        "schema": SCHEMA,
        "status": "running",
        "parameters": vars(args),
        "inputs": {str(path): sha256(path) for path in required},
    }
    if not args.skip_2um:
        summary["2um"] = run_2um(args)
    summary["V55"] = run_v55(args)
    summary["elapsed_seconds"] = time.time() - started
    summary["status"] = "complete"
    atomic_json(summary, args.output_dir / f"{args.sample}.fullgene_rank_sweep.qc.json")
    log(f"Completed full-gene rank sweep: {args.output_dir}")


if __name__ == "__main__":
    main()
