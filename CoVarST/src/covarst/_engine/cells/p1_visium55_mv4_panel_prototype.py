#!/usr/bin/env python3
"""P1 V55 Spot2Cell-MV4 278-gene development prototype.

The fit reads only the P1 V55 observation, CellViT geometry/features and the
independently learned V55 type/program model.  The existing 278-gene H5AD is
used as a schema and gene-name template; neither 2 um expression nor any
external single-cell reference is read during fitting or allocation.

This is a development-panel prototype.  It writes several predeclared
graph/morphology-strength candidates without selecting one from 2 um results.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
from scipy import sparse

import cell_resolved_visium55_mv3 as mv3
import spot2cell_mv4_core as mv4


EPS = 1.0e-12
CLASS_NAMES = mv3.CLASS_NAMES
MAX_PROGRAMS = mv3.MAX_PROGRAMS


def parse_beta_pairs(value: str) -> list[tuple[float, float]]:
    result = []
    for item in value.split(";"):
        graph, morph = item.split(",")
        pair = (float(graph), float(morph))
        if pair[0] < 0 or pair[1] < 0:
            raise argparse.ArgumentTypeError("beta values must be non-negative")
        result.append(pair)
    if not result:
        raise argparse.ArgumentTypeError("at least one beta pair is required")
    return result


def parse_args() -> argparse.Namespace:
    root = Path('DATA/Visium_HD')
    mv2 = root / "cell_resolved_spot2cell_mv2/P1/visium55_v2"
    mv3_root = root / "cell_resolved_spot2cell_mv3/P1/visium55_v3"
    out = root / "cell_resolved_spot2cell_mv4/P1/visium55_panel_prototype_v1"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mv2-h5ad", type=Path, default=mv2 / "P1.visium55.spot2cell_mv2.h5ad")
    parser.add_argument("--mv2-allocation", type=Path, default=mv2 / "P1.visium55.spot2cell_mv2.allocation.h5")
    parser.add_argument("--spot-positions", type=Path, default=mv2 / "P1.visium55.spots.tsv.gz")
    parser.add_argument("--mv3-state", type=Path, default=mv3_root / "P1.visium55.spot2cell_mv3.state.h5")
    parser.add_argument("--panel-template", type=Path, default=mv3_root / "P1.visium55.spot2cell_mv3.P_state_278.h5ad")
    parser.add_argument("--output-dir", type=Path, default=out)
    parser.add_argument("--permutations", type=int, default=20)
    parser.add_argument("--spatial-folds", type=int, default=5)
    parser.add_argument("--spatial-buffer", type=float, default=1.5)
    parser.add_argument("--graph-neighbors", type=int, default=10)
    parser.add_argument("--graph-radius", type=float, default=3.5)
    parser.add_argument("--boundary-strength", type=float, default=3.0)
    parser.add_argument("--lambda-prior", type=float, default=0.25)
    parser.add_argument("--lambda-graph", type=float, default=0.50)
    parser.add_argument("--program-genes", type=int, default=3000)
    parser.add_argument("--residual-genes", type=int, default=768)
    parser.add_argument("--residual-pcs", type=int, default=16)
    parser.add_argument(
        "--beta-pairs", type=parse_beta_pairs,
        default=parse_beta_pairs("0,0;0.25,0.25;0.5,0.25;0.5,0.5;1,0.5;1,1"),
        help="semicolon-separated graph,morph pairs",
    )
    parser.add_argument("--seed", type=int, default=20260827)
    return parser.parse_args()


def log(message: str) -> None:
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), message, flush=True)


def robust_standardize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    median = np.nanmedian(values, axis=0)
    mad = 1.4826 * np.nanmedian(np.abs(values - median), axis=0)
    std = np.nanstd(values, axis=0)
    scale = np.where(mad > 1.0e-6, mad, np.where(std > 1.0e-6, std, 1.0))
    result = (values - median) / scale
    result = np.clip(np.nan_to_num(result), -6.0, 6.0)
    keep = np.std(result, axis=0) > 1.0e-8
    return result[:, keep]


def load_data(args: argparse.Namespace) -> tuple[dict, dict]:
    for path in (
        args.mv2_h5ad, args.mv2_allocation, args.spot_positions,
        args.mv3_state, args.panel_template,
    ):
        if not path.exists():
            raise FileNotFoundError(path)
    bridge = SimpleNamespace(
        mv2_h5ad=args.mv2_h5ad,
        mv2_allocation=args.mv2_allocation,
        spot_positions=args.spot_positions,
        panel_reference=args.panel_template,
        state_output=args.mv3_state,
    )
    data = mv3.load_inputs(bridge)
    _, mv3_summary = mv3.load_state(bridge, data)
    # MV4 deliberately disables MV3's direct per-gene feature map.  A low-rank
    # weak state must pass first; gene residuals are a later, separate gate.
    data.pop("feature_gene_coef", None)
    data.pop("feature_gene_clip", None)
    data.pop("cell_feature_basis", None)
    return data, mv3_summary


def graph_gene_context(args: argparse.Namespace, data: dict) -> tuple[np.ndarray, dict]:
    composition, _, n_cells = mv3.core.census_from_geometry(
        data["geometry"], data["class_id"], data["capacity"],
    )
    library = np.asarray(data["matrix"].sum(axis=0)).ravel().astype(np.float64)
    occupied = n_cells > 0
    _, residual_genes, selection_audit = mv3.core.select_program_genes(
        data["matrix"], library, occupied, data["gene_names"], data["marker_sets"],
        args.program_genes, args.residual_genes, args.seed + 41000,
    )
    graph_genes = residual_genes[::2]
    residual_scores, pca = mv3.core.compute_residual_pcs(
        data["matrix"], library, composition, data["profiles"], graph_genes,
        occupied, args.residual_pcs, args.seed + 42000,
    )
    spot_context = robust_standardize(np.column_stack((
        composition,
        residual_scores,
        np.log1p(n_cells)[:, None],
        np.log1p(library)[:, None],
    )))
    audit = {
        "graph_gene_count": int(len(graph_genes)),
        "graph_gene_indices_sha256": mv3.hashlib.sha256(
            np.asarray(graph_genes, dtype=np.int32).tobytes()
        ).hexdigest(),
        "graph_context_dimensions": int(spot_context.shape[1]),
        "residual_pca_explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "gene_selection_audit": selection_audit,
    }
    return spot_context, audit


def group_context(
    data: dict,
    group_bins: np.ndarray,
    measurement: sparse.csr_matrix,
    spot_context: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    coords = data["spots"].iloc[group_bins][["hex_col", "hex_row"]].to_numpy(np.float64)
    spatial = robust_standardize(coords)
    spatial_basis = np.column_stack((
        spatial,
        np.square(spatial),
        spatial[:, 0] * spatial[:, 1],
    ))
    row_square = np.asarray(measurement.multiply(measurement).sum(1)).ravel()
    effective_n = np.divide(1.0, row_square, out=np.ones_like(row_square), where=row_square > EPS)
    context = robust_standardize(np.column_stack((
        spatial_basis,
        spot_context[group_bins],
        np.log1p(effective_n)[:, None],
    )))
    return coords, context


def fit_states(args: argparse.Namespace, data: dict, spot_context: np.ndarray):
    n_cells = len(data["class_id"])
    graph_contrast = np.zeros((n_cells, MAX_PROGRAMS), np.float32)
    morph_contrast = np.zeros((n_cells, MAX_PROGRAMS), np.float32)
    all_edges = []
    audits = []
    safe_owner = np.clip(data["owner_bin"], 0, len(data["spots"]) - 1)
    for class_value, class_name in enumerate(CLASS_NAMES):
        cells, groups, measurement = mv3.class_measurement_matrix(data, class_value)
        k = int(data["program_count"][class_value])
        class_audit = {
            "class": class_name,
            "n_cells": int(len(cells)),
            "n_groups": int(len(groups)),
            "program_count": k,
        }
        if k <= 0 or len(groups) < max(40, args.spatial_folds * 6):
            class_audit["gate"] = "class_global"
            class_audit["reason"] = "no stable programs or insufficient groups"
            audits.append(class_audit)
            continue
        bins = data["group_bin"][groups]
        target = data["group_scores"][groups, :k]
        coords, context = group_context(data, bins, measurement, spot_context)
        folds = mv4.contiguous_spatial_folds(coords, args.spatial_folds)
        bin_to_fold = np.full(data["geometry"].shape[0], -1, np.int32)
        bin_to_fold[bins] = folds
        local_cell_fold = bin_to_fold[safe_owner[cells]]
        # A rare boundary-crossing cell whose owner lacks a same-class group is
        # assigned the fold of its strongest measurement row.
        missing_fold = local_cell_fold < 0
        if np.any(missing_fold):
            membership = measurement[:, missing_fold].tocsc()
            local_cell_fold[missing_fold] = folds[np.asarray(membership.argmax(axis=0)).ravel()]
        log(f"{class_name}: fitting cross-fitted slice-internal weak model")
        weak = mv4.fit_slice_internal_weak_state(
            measurement=measurement,
            cell_features=data["features"][cells],
            group_target=target,
            group_context=context,
            group_coords=coords,
            cell_fold=local_cell_fold,
            n_folds=args.spatial_folds,
            permutations=args.permutations,
            buffer_distance=args.spatial_buffer,
            seed=args.seed + class_value * 1000,
        )
        local_boundary_context = spot_context[safe_owner[cells]]
        graph = mv4.build_boundary_aware_graph(
            coords=data["coords"][cells],
            features=data["features"][cells],
            boundary_context=local_boundary_context,
            component=data["cell_component"][cells],
            owner_group=safe_owner[cells],
            neighbors=args.graph_neighbors,
            radius=args.graph_radius,
            boundary_strength=args.boundary_strength,
        )
        log(f"{class_name}: solving boundary-aware graph inverse")
        _, local_graph, inverse_audit = mv4.solve_graph_inverse(
            measurement=measurement,
            group_target=target,
            laplacian=graph.laplacian,
            weak_prior=weak.morph_contrast,
            lambda_prior=args.lambda_prior,
            lambda_graph=args.lambda_graph,
        )
        local_morph, orthogonal_audit = mv4.orthogonalize_nullspace_contrasts(
            measurement, local_graph, weak.morph_contrast,
        )
        graph_contrast[cells, :k] = local_graph[:, :k]
        morph_contrast[cells, :k] = local_morph[:, :k]
        if len(graph.edges):
            global_edge = graph.edges.copy()
            global_edge[:, 0] = cells[global_edge[:, 0].astype(np.int64)]
            global_edge[:, 1] = cells[global_edge[:, 1].astype(np.int64)]
            all_edges.append(np.column_stack((
                global_edge,
                np.full(len(global_edge), class_value, np.float64),
            )))
        class_audit.update({
            "gate": "morphology_refined" if np.any(weak.eta > 0) else "graph_refined",
            "weak_model": weak.audit,
            "graph": graph.audit,
            "inverse": inverse_audit,
            "orthogonalization": orthogonal_audit,
        })
        audits.append(class_audit)
        log(
            f"{class_name}: eta>0 {np.count_nonzero(weak.eta)}/{k}; "
            f"edges={graph.audit['n_edges']:,}"
        )
    edges = np.concatenate(all_edges) if all_edges else np.empty((0, 6), np.float64)
    return graph_contrast, morph_contrast, edges, audits


def spot_class_state(data: dict) -> np.ndarray:
    result = np.zeros(
        (data["geometry"].shape[0], len(CLASS_NAMES), MAX_PROGRAMS), np.float32,
    )
    result[data["group_bin"], data["group_class"]] = data["group_scores"].astype(np.float32)
    return result


def panel_indices(args: argparse.Namespace, data: dict):
    with h5py.File(args.panel_template, "r") as h:
        panel_names = mv3.read_var_names(h)
    lookup = {str(gene): index for index, gene in enumerate(data["gene_names"])}
    index = np.asarray([lookup.get(str(gene), -1) for gene in panel_names], np.int32)
    if np.any(index < 0):
        raise RuntimeError("V55 gene space lacks one or more prototype panel genes")
    reverse = np.full(len(data["gene_names"]), -1, np.int32)
    reverse[index] = np.arange(len(index), dtype=np.int32)
    return panel_names, index, reverse


def candidate_name(beta_graph: float, beta_morph: float) -> str:
    return f"G{int(round(beta_graph * 100)):03d}_M{int(round(beta_morph * 100)):03d}"


def allocate_candidate(
    data: dict,
    panel_index: np.ndarray,
    panel_reverse: np.ndarray,
    local_state: np.ndarray,
) -> tuple[np.ndarray, dict]:
    n_cells = len(data["class_id"])
    mass = np.zeros((n_cells, len(panel_index)), np.float32)
    geometry = data["geometry"]
    matrix = data["matrix"]
    occupied = np.diff(geometry.indptr) > 0
    observed = np.asarray(matrix[panel_index] @ occupied.astype(np.float64)).ravel()
    assigned = np.zeros(len(panel_index), np.float64)
    max_bin_gene_error = 0.0
    last = time.time()
    for b in np.flatnonzero(occupied):
        gs, ge = int(geometry.indptr[b]), int(geometry.indptr[b + 1])
        ms, me = int(matrix.indptr[b]), int(matrix.indptr[b + 1])
        cells = geometry.indices[gs:ge]
        overlap = geometry.data[gs:ge].astype(np.float64)
        genes = matrix.indices[ms:me]
        counts = matrix.data[ms:me].astype(np.float64)
        local_panel = panel_reverse[genes]
        keep = local_panel >= 0
        if not np.any(keep):
            continue
        genes = genes[keep]
        counts = counts[keep]
        local_panel = local_panel[keep]
        link_capacity = overlap * data["capacity"][cells]
        log_propensity = np.empty((len(cells), len(genes)), np.float64)
        classes = data["class_id"][cells]
        for class_value in np.unique(classes):
            selected = classes == class_value
            k = int(data["program_count"][class_value])
            modifier = np.zeros((int(np.sum(selected)), len(genes)), np.float64)
            if k > 0:
                modifier = (
                    local_state[cells[selected], :k]
                    @ data["loadings"][class_value, :k][:, genes]
                )
            log_propensity[selected] = (
                np.log(np.maximum(link_capacity[selected, None], EPS))
                + np.log(np.maximum(data["profiles"][class_value, genes][None, :], EPS))
                + np.clip(modifier, -6.0, 6.0)
            )
        _, block = mv4.exact_mass_projection(
            counts, log_propensity, link_capacity,
            library_ipf_strength=0.0, library_ipf_iterations=0,
        )
        mass[np.ix_(cells, local_panel)] += block.astype(np.float32)
        reconstructed = np.sum(block, axis=0)
        np.add.at(assigned, local_panel, reconstructed)
        max_bin_gene_error = max(
            max_bin_gene_error,
            float(np.max(np.abs(reconstructed - counts), initial=0.0)),
        )
        if time.time() - last > 60:
            log(f"projected {b + 1:,}/{matrix.shape[1]:,} V55 spots")
            last = time.time()
    return mass, {
        "panel_observed_mass": float(np.sum(observed)),
        "panel_assigned_mass": float(np.sum(assigned)),
        "max_abs_panel_gene_conservation_error": float(np.max(np.abs(observed - assigned))),
        "max_abs_bin_gene_error": max_bin_gene_error,
        "library_ipf_strength": 0.0,
    }


def write_candidate(
    args: argparse.Namespace,
    name: str,
    mass: np.ndarray,
    cell_state: np.ndarray,
    beta_graph: float,
    beta_morph: float,
    audit: dict,
) -> Path:
    output = args.output_dir / f"P1.visium55.spot2cell_mv4.{name}.X_recon_278.h5ad"
    partial = output.with_suffix(output.suffix + ".partial")
    if output.exists() or partial.exists():
        raise FileExistsError(f"Refusing to overwrite: {output}")
    with h5py.File(args.panel_template, "r") as old, h5py.File(partial, "w") as h:
        for key, value in old.attrs.items():
            h.attrs[key] = value
        h.attrs.update({
            "schema": "visium55_spot2cell_mv4_X_recon_panel_prototype_v1",
            "complete": 0,
            "X_semantics": "exact V55 spot-gene-conserved fractional reconstructed mass",
            "external_reference_used_for_fit": False,
            "beta_graph": beta_graph,
            "beta_morph": beta_morph,
            "library_ipf_strength": 0.0,
        })
        for key in old.keys():
            if key != "X":
                old.copy(key, h)
        if "program_activity" in h["obsm"]:
            del h["obsm/program_activity"]
        mv3.encoded_array(
            h["obsm"], "program_activity", cell_state.astype(np.float32), compression="lzf",
        )
        x = h.create_dataset(
            "X", data=mass, dtype=np.float32,
            chunks=(min(2048, len(mass)), mass.shape[1]), compression="lzf",
        )
        x.attrs["encoding-type"] = "array"
        x.attrs["encoding-version"] = "0.2.0"
        if "P_state_semantics" in h["uns"]:
            del h["uns/P_state_semantics"]
        mv3.encoded_scalar_string(
            h["uns"], "X_recon_semantics",
            "MV4 responsibility-level graph+morphology model with exact V55 spot-gene projection",
        )
        mv3.encoded_scalar_string(
            h["uns"], "mv4_candidate_audit_json", json.dumps(audit, ensure_ascii=False),
        )
        h.attrs["complete"] = 1
        h.flush()
    os.replace(partial, output)
    return output


def atomic_json(value: dict, path: Path) -> None:
    def json_default(item):
        if isinstance(item, np.ndarray):
            return item.tolist()
        if isinstance(item, np.generic):
            return item.item()
        if isinstance(item, Path):
            return str(item)
        raise TypeError(f"Object of type {type(item).__name__} is not JSON serializable")

    partial = path.with_suffix(path.suffix + ".partial")
    with partial.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, default=json_default)
        handle.write("\n")
    os.replace(partial, path)


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite prototype directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True)
    log("Loading frozen V55 observation, geometry and independent type/program model")
    data, mv3_summary = load_data(args)
    log("Computing disjoint graph-gene boundary context from V55 only")
    spot_context, graph_gene_audit = graph_gene_context(args, data)
    graph_contrast, morph_contrast, edges, class_audits = fit_states(
        args, data, spot_context,
    )
    base_lookup = spot_class_state(data)
    safe_owner = np.clip(data["owner_bin"], 0, base_lookup.shape[0] - 1)
    base_cell_state = base_lookup[safe_owner, data["class_id"]]
    panel_names, panel_index, panel_reverse = panel_indices(args, data)
    state_output = args.output_dir / "P1.visium55.spot2cell_mv4.state.h5"
    with h5py.File(state_output.with_suffix(".h5.partial"), "w") as h:
        h.attrs.update({
            "schema": "P1_visium55_spot2cell_mv4_state_panel_prototype_v1",
            "complete": 0,
            "external_reference_used_for_fit": False,
            "fixed_44_programs_used": False,
        })
        h.create_dataset("cell_id", data=data["cell_id"], compression="lzf")
        h.create_dataset("class_id", data=data["class_id"].astype(np.uint8), compression="lzf")
        h.create_dataset("graph_contrast", data=graph_contrast, compression="lzf")
        h.create_dataset("morph_contrast", data=morph_contrast, compression="lzf")
        edge_group = h.create_group("boundary_aware_graph")
        for column, name, dtype in (
            (0, "source_cell_index", np.int32),
            (1, "target_cell_index", np.int32),
            (2, "weight", np.float32),
            (3, "biological_boundary_score", np.float32),
            (4, "cross_technical_spot", np.uint8),
            (5, "class_id", np.uint8),
        ):
            edge_group.create_dataset(name, data=edges[:, column].astype(dtype), compression="lzf")
        h.attrs["complete"] = 1
        h.flush()
    os.replace(state_output.with_suffix(".h5.partial"), state_output)

    candidates = []
    for beta_graph, beta_morph in args.beta_pairs:
        name = candidate_name(beta_graph, beta_morph)
        log(f"Allocating candidate {name}")
        cell_state = (
            base_cell_state
            + beta_graph * graph_contrast
            + beta_morph * morph_contrast
        ).astype(np.float32)
        mass, allocation_audit = allocate_candidate(
            data, panel_index, panel_reverse, cell_state,
        )
        candidate_audit = {
            "name": name,
            "beta_graph": beta_graph,
            "beta_morph": beta_morph,
            "allocation": allocation_audit,
            "selection_status": "predeclared_unselected_development_candidate",
        }
        output = write_candidate(
            args, name, mass, cell_state, beta_graph, beta_morph, candidate_audit,
        )
        candidate_audit["path"] = str(output)
        candidates.append(candidate_audit)
        log(
            f"{name}: conservation error "
            f"{allocation_audit['max_abs_panel_gene_conservation_error']:.3g}"
        )
        del mass, cell_state
    summary = {
        "schema": "P1_visium55_spot2cell_mv4_panel_prototype_summary_v1",
        "status": "development_candidates_not_frozen",
        "n_cells": int(len(data["cell_id"])),
        "n_genes": int(len(panel_names)),
        "n_graph_edges": int(len(edges)),
        "external_reference_used_for_fit": False,
        "fixed_44_programs_used": False,
        "mv3_type_program_model_provenance": mv3_summary.get("schema", "unknown"),
        "graph_gene_context": graph_gene_audit,
        "classes": class_audits,
        "candidates": candidates,
        "parameters": {
            "spatial_folds": args.spatial_folds,
            "spatial_buffer": args.spatial_buffer,
            "permutations": args.permutations,
            "graph_neighbors": args.graph_neighbors,
            "graph_radius": args.graph_radius,
            "boundary_strength": args.boundary_strength,
            "lambda_prior": args.lambda_prior,
            "lambda_graph": args.lambda_graph,
            "seed": args.seed,
        },
    }
    atomic_json(summary, args.output_dir / "P1.visium55.spot2cell_mv4.prototype.qc.json")
    log(f"Completed MV4 panel prototype: {args.output_dir}")


if __name__ == "__main__":
    main()
