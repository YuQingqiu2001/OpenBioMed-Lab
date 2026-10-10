#!/usr/bin/env python3
"""Sample-specific virtual-Visium-55 MV3 spot-to-cell reconstruction.

MV3 keeps the independently fitted MV2 spot/type expression model, but replaces
the under-identified within-spot allocation layer with:

1. a same-CellViT-class cell graph;
2. a global graph-regularised inverse for continuous cell state ``z``;
3. a factorised, library-free cell expression state ``P_state``; and
4. a joint per-spot KL/IPF mass projection that exactly conserves every observed
   spot-by-gene count while softly regularising cell library sizes.

The 2 um direct data are never read by this script.  They are reserved for the
separate post-hoc validation script.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import linalg as splinalg
from scipy.spatial import cKDTree

import cell_resolved_v3 as core
from cell_resolved_visium55 import hex_components


EPS = 1.0e-12
OLD_MAX_PROGRAMS = 4
MAX_PROGRAMS = 8
CLASS_NAMES = list(core.CLASS_NAMES)


def parse_args() -> argparse.Namespace:
    root = Path('DATA/Visium_HD')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", choices=("P1", "P2", "P5"), default="P1")
    parser.add_argument(
        "--mv2-h5ad", type=Path,
        default=None,
    )
    parser.add_argument(
        "--mv2-allocation", type=Path,
        default=None,
    )
    parser.add_argument(
        "--spot-positions", type=Path,
        default=None,
    )
    parser.add_argument(
        "--state-output", type=Path,
        default=None,
    )
    parser.add_argument(
        "--output", type=Path,
        default=None,
    )
    parser.add_argument(
        "--allocation-output", type=Path,
        default=None,
    )
    parser.add_argument(
        "--pstate-panel-output", type=Path,
        default=None,
    )
    parser.add_argument(
        "--xmass-panel-output", type=Path,
        default=None,
    )
    parser.add_argument(
        "--qc-json", type=Path,
        default=None,
    )
    parser.add_argument(
        "--panel-reference", type=Path,
        default=None,
        help="Used only to obtain the common 278-gene names; expression is never read.",
    )
    parser.add_argument(
        "--mode", choices=("state", "panel", "allocate", "all"), default="all",
    )
    parser.add_argument("--graph-neighbors", type=int, default=10)
    parser.add_argument("--graph-radius-array-units", type=float, default=3.5)
    parser.add_argument("--ridge-permutations", type=int, default=100)
    parser.add_argument(
        "--reuse-mv2-programs", action="store_true",
        help="Diagnostic fallback only. Default MV3 refits 0-8 type programs from V55.",
    )
    parser.add_argument("--program-genes", type=int, default=3000)
    parser.add_argument("--residual-genes", type=int, default=768)
    parser.add_argument("--residual-pcs", type=int, default=16)
    parser.add_argument("--program-cv-folds", type=int, default=5)
    parser.add_argument("--program-bootstrap-replicates", type=int, default=20)
    parser.add_argument("--feature-gene-permutations", type=int, default=20)
    parser.add_argument("--lambda-feature-multiplier", type=float, default=0.25)
    parser.add_argument("--lambda-graph-multiplier", type=float, default=0.50)
    parser.add_argument("--library-ipf-strength", type=float, default=0.25)
    parser.add_argument("--library-ipf-iterations", type=int, default=4)
    parser.add_argument("--cell-batch-size", type=int, default=1024)
    parser.add_argument("--seed", type=int, default=20260827)
    args = parser.parse_args()
    sample = args.sample
    mv2 = root / f"cell_resolved_spot2cell_mv2/{sample}/visium55_v2"
    out = root / f"cell_resolved_spot2cell_mv3/{sample}/visium55_v3"
    defaults = {
        "mv2_h5ad": mv2 / f"{sample}.visium55.spot2cell_mv2.h5ad",
        "mv2_allocation": mv2 / f"{sample}.visium55.spot2cell_mv2.allocation.h5",
        "spot_positions": mv2 / f"{sample}.visium55.spots.tsv.gz",
        "state_output": out / f"{sample}.visium55.spot2cell_mv3.state.h5",
        "output": out / f"{sample}.visium55.spot2cell_mv3.h5ad",
        "allocation_output": out / f"{sample}.visium55.spot2cell_mv3.allocation.h5",
        "pstate_panel_output": out / f"{sample}.visium55.spot2cell_mv3.P_state_278.h5ad",
        "xmass_panel_output": out / f"{sample}.visium55.spot2cell_mv3.X_mass_278.precheck.h5ad",
        "qc_json": out / f"{sample}.visium55.spot2cell_mv3.qc.json",
        "panel_reference": root / f"cell_resolved_002um_voronoi/{sample}/{sample}.Visium2_direct.full_support_input.h5ad",
    }
    for name, value in defaults.items():
        if getattr(args, name) is None:
            setattr(args, name, value)
    if not 0.0 <= args.library_ipf_strength <= 1.0:
        parser.error("--library-ipf-strength must be in [0,1]")
    if args.graph_neighbors < 2:
        parser.error("--graph-neighbors must be at least 2")
    return args


def log(message: str) -> None:
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), message, flush=True)


def decode(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values)
    if values.dtype.kind in {"S", "O"}:
        return np.asarray([
            item.decode("utf-8") if isinstance(item, bytes) else str(item)
            for item in values
        ], dtype="U")
    return values.astype("U")


def sha256(path: Path, block: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(block)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def read_obs(handle: h5py.File, name: str) -> np.ndarray:
    item = handle[f"obs/{name}"]
    if isinstance(item, h5py.Group):
        codes = item["codes"][:]
        categories = decode(item["categories"][:])
        result = np.full(len(codes), "", dtype="U128")
        valid = codes >= 0
        result[valid] = categories[codes[valid]]
        return result
    return item[:]


def read_var_names(handle: h5py.File) -> np.ndarray:
    if "gene_name" in handle["var"]:
        return decode(handle["var/gene_name"][:])
    index_name = handle["var"].attrs.get("_index", "_index")
    if isinstance(index_name, bytes):
        index_name = index_name.decode("utf-8")
    return decode(handle[f"var/{index_name}"][:])


def load_inputs(args: argparse.Namespace) -> dict:
    for path in (
        args.mv2_h5ad, args.mv2_allocation, args.spot_positions,
        args.panel_reference,
    ):
        if not path.exists():
            raise FileNotFoundError(path)
    with h5py.File(args.mv2_h5ad, "r") as h:
        class_id = read_obs(h, "cellvit_class_id").astype(np.int32)
        owner_bin = read_obs(h, "owner_bin_index").astype(np.int32)
        cell_id = read_obs(h, "cell_id").astype(np.int64)
        capacity = read_obs(h, "rna_capacity").astype(np.float64)
        coords = h["obsm/spatial_array_uv"][:].astype(np.float64)
        features = h["obsm/standardized_nucleus_features"][:].astype(np.float64)
        features = np.clip(np.nan_to_num(features), -6.0, 6.0)
        profiles = h["uns/class_gene_probability"][:].astype(np.float64)
        flat_loadings = h["varm/program_loadings"][:].astype(np.float64)
        n_genes = profiles.shape[1]
        old_loadings = flat_loadings.reshape(
            n_genes, len(CLASS_NAMES), OLD_MAX_PROGRAMS
        )
        old_loadings = np.transpose(old_loadings, (1, 2, 0))
        loadings = np.zeros((len(CLASS_NAMES), MAX_PROGRAMS, n_genes), np.float64)
        loadings[:, :OLD_MAX_PROGRAMS] = old_loadings
        old_program_count = h["uns/program_count_by_class"][:].astype(np.int32)
        program_count = old_program_count.copy()
        group_bin = h["uns/bin_type_bin_index"][:].astype(np.int32)
        group_class = h["uns/bin_type_class_id"][:].astype(np.int32)
        old_group_scores = h["uns/bin_type_program_activity"][:].astype(np.float64)
        group_scores = np.zeros((len(old_group_scores), MAX_PROGRAMS), np.float64)
        group_scores[:, :OLD_MAX_PROGRAMS] = old_group_scores
        gene_names = read_var_names(h)
        marker_json = h["uns/marker_sets_json"][()]
        if isinstance(marker_json, bytes):
            marker_json = marker_json.decode("utf-8")
        marker_sets = json.loads(str(marker_json))
    with h5py.File(args.mv2_allocation, "r") as h:
        geometry_shape = tuple(h["geometry/shape"][:].astype(np.int64))
        geometry = sparse.csr_matrix(
            (
                h["geometry/overlap_fraction"][:].astype(np.float64),
                h["geometry/cell_index"][:].astype(np.int32),
                h["geometry/indptr"][:].astype(np.int32),
            ),
            shape=geometry_shape,
        )
        matrix_shape = tuple(h["observed/shape"][:].astype(np.int64))
        matrix = sparse.csc_matrix(
            (
                h["observed/count"][:].astype(np.float64),
                h["observed/gene_index"][:].astype(np.int32),
                h["observed/bin_indptr"][:].astype(np.int32),
            ),
            shape=matrix_shape,
        )
    spots = pd.read_csv(args.spot_positions, sep="\t")
    if len(spots) != geometry.shape[0]:
        raise RuntimeError("Spot positions and allocation geometry are not aligned")
    spot_component = hex_components(spots)
    safe_owner = np.clip(owner_bin, 0, len(spots) - 1)
    cell_component = spot_component[safe_owner]
    return {
        "class_id": class_id,
        "owner_bin": owner_bin,
        "cell_id": cell_id,
        "capacity": capacity,
        "coords": coords,
        "features": features,
        "profiles": profiles,
        "loadings": loadings,
        "program_count": program_count,
        "group_bin": group_bin,
        "group_class": group_class,
        "group_scores": group_scores,
        "marker_sets": marker_sets,
        "gene_names": gene_names,
        "geometry": geometry,
        "matrix": matrix,
        "spots": spots,
        "spot_component": spot_component,
        "cell_component": cell_component,
    }


def refit_mv3_programs(args: argparse.Namespace, data: dict) -> None:
    """Refit 0-8 type-specific axes from V55 only; never reuse the fixed 44 axes."""
    log("Refitting independent 0-8 type-specific Programs from V55 residuals")
    composition, _, n_cells_per_spot = core.census_from_geometry(
        data["geometry"], data["class_id"], data["capacity"],
    )
    groups = core.aggregate_bin_type_features(
        data["geometry"], data["class_id"], data["capacity"], data["features"],
    )
    library = np.asarray(data["matrix"].sum(axis=0)).ravel().astype(np.float64)
    occupied = n_cells_per_spot > 0
    program_genes, residual_genes, _ = core.select_program_genes(
        data["matrix"], library, occupied,
        data["gene_names"], data["marker_sets"],
        args.program_genes, args.residual_genes, args.seed + 31000,
    )
    graph_genes = residual_genes[::2]
    state_selection_genes = residual_genes[1::2]
    program_genes = np.setdiff1d(
        program_genes, graph_genes, assume_unique=False,
    )
    residual_scores, residual_pca = core.compute_residual_pcs(
        data["matrix"], library, composition, data["profiles"], graph_genes,
        occupied, args.residual_pcs, args.seed + 32000,
    )
    previous_max = core.MAX_PROGRAMS
    core.MAX_PROGRAMS = MAX_PROGRAMS
    try:
        group_scores, loadings, program_count, audit = core.learn_local_programs(
            data["matrix"], library, data["spots"], data["spot_component"],
            composition, data["profiles"], residual_scores,
            state_selection_genes, program_genes, groups,
            neighbors=24, radius=37.5,
            folds=args.program_cv_folds,
            bootstrap_replicates=args.program_bootstrap_replicates,
            seed=args.seed + 33000,
        )
    finally:
        core.MAX_PROGRAMS = previous_max
    data["group_bin"] = groups["bin_index"].astype(np.int32)
    data["group_class"] = groups["class_id"].astype(np.int32)
    data["group_scores"] = group_scores.astype(np.float64)
    data["loadings"] = loadings.astype(np.float64)
    data["program_count"] = program_count.astype(np.int32)
    feature_gene_coef, feature_gene_clip, feature_gene_audit = fit_direct_gene_feature_map(
        args, data, composition, groups, library, program_genes,
        group_scores, loadings, program_count,
    )
    data["feature_gene_coef"] = feature_gene_coef
    data["feature_gene_clip"] = feature_gene_clip
    data["cell_feature_basis"] = centred_cell_feature_basis(data)
    data["program_refit"] = {
            "schema": "V55_MV3_independent_type_programs_v1",
        "max_programs": MAX_PROGRAMS,
        "selected_program_count": program_count.astype(int).tolist(),
        "program_gene_count": int(len(program_genes)),
        "residual_graph_gene_count": int(len(graph_genes)),
        "state_selection_gene_count": int(len(state_selection_genes)),
        "fixed_44_programs_used": False,
        "audit": audit,
        "direct_feature_gene_mapping": feature_gene_audit,
        "residual_pca_explained_variance_ratio": (
            residual_pca.explained_variance_ratio_.astype(float).tolist()
        ),
    }


def class_measurement_matrix(data: dict, class_value: int):
    class_cells = np.flatnonzero(data["class_id"] == class_value).astype(np.int32)
    class_groups = np.flatnonzero(data["group_class"] == class_value).astype(np.int32)
    if not len(class_cells) or not len(class_groups):
        return class_cells, class_groups, sparse.csr_matrix((0, len(class_cells)))
    cell_local = np.full(len(data["class_id"]), -1, np.int32)
    cell_local[class_cells] = np.arange(len(class_cells), dtype=np.int32)
    bin_to_group = np.full(data["geometry"].shape[0], -1, np.int32)
    bin_to_group[data["group_bin"][class_groups]] = np.arange(len(class_groups), dtype=np.int32)
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []
    geometry = data["geometry"]
    for b in data["group_bin"][class_groups]:
        start, end = int(geometry.indptr[b]), int(geometry.indptr[b + 1])
        cells = geometry.indices[start:end]
        selected = data["class_id"][cells] == class_value
        if not np.any(selected):
            continue
        chosen = cells[selected]
        weight = geometry.data[start:end][selected] * data["capacity"][chosen]
        total = float(weight.sum())
        if total <= 0:
            weight = np.ones(len(chosen), np.float64)
            total = float(len(chosen))
        rows.append(np.full(len(chosen), bin_to_group[b], np.int32))
        cols.append(cell_local[chosen])
        vals.append(weight / total)
    if not rows:
        return class_cells, class_groups, sparse.csr_matrix((len(class_groups), len(class_cells)))
    measurement = sparse.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(len(class_groups), len(class_cells)),
    ).tocsr()
    measurement.sum_duplicates()
    return class_cells, class_groups, measurement


def spatial_fold_ids(spots: pd.DataFrame, bins: np.ndarray, n_folds: int = 5) -> np.ndarray:
    row = spots.iloc[bins]["hex_row"].to_numpy(np.int64)
    col = spots.iloc[bins]["hex_col"].to_numpy(np.int64)
    return np.mod(3 * row + 2 * col, n_folds).astype(np.int32)


def ridge_fit(x: np.ndarray, y: np.ndarray, alpha: float):
    mean_x = x.mean(0)
    scale_x = x.std(0)
    scale_x[scale_x < 1.0e-6] = 1.0
    mean_y = y.mean(0)
    z = (x - mean_x) / scale_x
    lhs = z.T @ z + float(alpha) * np.eye(z.shape[1])
    coef_standard = np.linalg.solve(lhs, z.T @ (y - mean_y))
    coef = coef_standard / scale_x[:, None]
    intercept = mean_y - mean_x @ coef
    return coef, intercept


def ridge_cv_predictions(
    x: np.ndarray, y: np.ndarray, folds: np.ndarray, alpha: float,
) -> np.ndarray:
    prediction = np.zeros_like(y)
    for fold in np.unique(folds):
        test = folds == fold
        train = ~test
        if int(train.sum()) <= x.shape[1] + 2 or not np.any(test):
            prediction[test] = y[train].mean(0) if np.any(train) else 0.0
            continue
        coef, intercept = ridge_fit(x[train], y[train], alpha)
        prediction[test] = x[test] @ coef + intercept
    return prediction


def r2_columns(truth: np.ndarray, prediction: np.ndarray) -> np.ndarray:
    denominator = np.sum(np.square(truth - truth.mean(0)), axis=0)
    numerator = np.sum(np.square(truth - prediction), axis=0)
    return np.divide(
        1.0 - numerator / np.maximum(denominator, EPS),
        1.0,
        out=np.full(truth.shape[1], -np.inf),
        where=denominator > EPS,
    )


def centred_cell_feature_basis(data: dict) -> np.ndarray:
    """Class-internal feature contrasts with local measurement means removed."""
    basis = np.zeros_like(data["features"], dtype=np.float32)
    for class_value in range(len(CLASS_NAMES)):
        cells, _, measurement = class_measurement_matrix(data, class_value)
        if not len(cells) or measurement.shape[0] == 0:
            continue
        local = data["features"][cells]
        group_mean = measurement @ local
        back_weight = np.asarray(measurement.T.sum(1)).ravel()
        back = np.asarray(measurement.T @ group_mean)
        back = np.divide(
            back, back_weight[:, None], out=np.zeros_like(back),
            where=back_weight[:, None] > 0,
        )
        basis[cells] = np.clip(local - back, -5.0, 5.0).astype(np.float32)
    return basis


def fit_direct_gene_feature_map(
    args: argparse.Namespace,
    data: dict,
    composition: np.ndarray,
    groups: dict,
    library: np.ndarray,
    program_genes: np.ndarray,
    group_scores: np.ndarray,
    loadings: np.ndarray,
    program_count: np.ndarray,
):
    """Map group-mean nuclear features to residual genes, then apply contrasts.

    Coefficients are learned only when spatially held-out prediction beats a
    permutation-derived null.  At cell level the coefficient acts on locally
    centred feature contrasts, so it cannot simply copy an ecological group
    mean into every nucleus.
    """
    n_features = data["features"].shape[1]
    n_genes = data["matrix"].shape[0]
    coefficient = np.zeros((len(CLASS_NAMES), n_features, n_genes), np.float32)
    clip_value = np.zeros((len(CLASS_NAMES), n_genes), np.float32)
    audits: dict[str, dict] = {}
    rng = np.random.default_rng(args.seed + 44000)
    alphas = (10.0, 100.0, 1000.0)
    for class_value, class_name in enumerate(CLASS_NAMES):
        group_index = np.flatnonzero(groups["class_id"] == class_value)
        if len(group_index) < 80:
            audits[class_name] = {
                "n_groups": int(len(group_index)), "passed_genes": 0,
                "reason": "fewer than 80 bin-type groups",
            }
            continue
        bins = groups["bin_index"][group_index]
        fraction_g, signal_g = core.type_signal(
            data["matrix"], program_genes, bins, class_value,
            composition, data["profiles"], library,
        )
        baseline = np.maximum(data["profiles"][class_value, program_genes], EPS)
        target = np.maximum(
            signal_g / np.maximum(fraction_g[:, None], 1.0e-6),
            baseline[None, :] * 1.0e-3,
        )
        target = np.clip(
            np.log(target) - np.log(baseline[None, :]), -5.0, 5.0,
        )
        k = int(program_count[class_value])
        if k > 0:
            target -= (
                group_scores[group_index, :k]
                @ loadings[class_value, :k][:, program_genes]
            )
        x = groups["feature_mean"][group_index].astype(np.float64)
        folds = spatial_fold_ids(data["spots"], bins)
        sample_size = min(512, target.shape[1])
        sample_gene = np.sort(rng.choice(target.shape[1], sample_size, replace=False))
        alpha_score = []
        for alpha in alphas:
            prediction = ridge_cv_predictions(x, target[:, sample_gene], folds, alpha)
            alpha_score.append(float(np.nanmedian(r2_columns(target[:, sample_gene], prediction))))
        selected_alpha = float(alphas[int(np.argmax(alpha_score))])
        prediction = ridge_cv_predictions(x, target, folds, selected_alpha)
        heldout_r2 = r2_columns(target, prediction)
        null_values = []
        null_gene = sample_gene[:min(128, len(sample_gene))]
        for _ in range(args.feature_gene_permutations):
            permuted = target[rng.permutation(len(target))][:, null_gene]
            null_prediction = ridge_cv_predictions(x, permuted, folds, selected_alpha)
            null_values.append(r2_columns(permuted, null_prediction))
        null_threshold = float(
            np.quantile(np.concatenate(null_values), 0.95)
        ) if null_values else 0.0
        excess = heldout_r2 - max(null_threshold, 0.0)
        eta = np.zeros_like(excess)
        eta[excess > 0.0] = 0.25
        eta[excess > 0.005] = 0.50
        eta[excess > 0.020] = 0.75
        eta[excess > 0.050] = 1.00
        coef, _ = ridge_fit(x, target, selected_alpha)
        coef *= eta[None, :]

        # Calibrate the permitted cell contrast from V55 target variability,
        # never from 2 um.  The cap prevents nuclear features from creating more
        # within-cell log-ratio spread than supported across measured spots.
        class_cells, _, measurement = class_measurement_matrix(data, class_value)
        local_features = data["features"][class_cells]
        group_mean = measurement @ local_features
        back_weight = np.asarray(measurement.T.sum(1)).ravel()
        back = np.asarray(measurement.T @ group_mean)
        back = np.divide(
            back, back_weight[:, None], out=np.zeros_like(back),
            where=back_weight[:, None] > 0,
        )
        contrast = np.clip(local_features - back, -5.0, 5.0)
        calibration_gene = np.flatnonzero(eta > 0)
        scale = np.ones(target.shape[1], np.float64)
        for start in range(0, len(calibration_gene), 256):
            local_gene = calibration_gene[start:start + 256]
            effect_sd = np.std(contrast @ coef[:, local_gene], axis=0)
            target_sd = np.std(target[:, local_gene], axis=0)
            scale[local_gene] = np.minimum(
                np.divide(
                    0.75 * target_sd, np.maximum(effect_sd, EPS),
                    out=np.ones_like(target_sd), where=effect_sd > EPS,
                ),
                3.0,
            )
        coef *= scale[None, :]
        coefficient[class_value, :, program_genes] = coef.T.astype(np.float32)
        clip = 2.0 * np.std(target, axis=0)
        clip_value[class_value, program_genes] = np.clip(clip, 0.0, 5.0).astype(np.float32)
        audits[class_name] = {
            "n_groups": int(len(group_index)),
            "candidate_genes": int(len(program_genes)),
            "selected_alpha": selected_alpha,
            "alpha_cv_median_r2": alpha_score,
            "permutation_replicates": int(args.feature_gene_permutations),
            "permutation_r2_q95": null_threshold,
            "passed_genes": int(np.count_nonzero(eta > 0)),
            "eta_counts": {
                str(value): int(np.count_nonzero(eta == value))
                for value in (0.0, 0.25, 0.5, 0.75, 1.0)
            },
            "heldout_r2_median_passed": (
                float(np.median(heldout_r2[eta > 0])) if np.any(eta > 0) else None
            ),
        }
        log(
            f"{class_name}: direct nucleus-feature gene gate passed "
            f"{np.count_nonzero(eta > 0):,}/{len(program_genes):,} genes"
        )
    return coefficient, clip_value, audits


def feature_prior(
    measurement: sparse.csr_matrix,
    group_features: np.ndarray,
    cell_features: np.ndarray,
    target: np.ndarray,
    folds: np.ndarray,
    permutations: int,
    seed: int,
):
    alphas = (1.0, 10.0, 100.0, 1000.0)
    candidate = []
    for alpha in alphas:
        prediction = ridge_cv_predictions(group_features, target, folds, alpha)
        candidate.append(r2_columns(target, prediction))
    candidate = np.asarray(candidate)
    selected_alpha = float(alphas[int(np.argmax(np.nanmedian(candidate, axis=1)))])
    observed = r2_columns(
        target,
        ridge_cv_predictions(group_features, target, folds, selected_alpha),
    )
    rng = np.random.default_rng(seed)
    null = np.zeros((permutations, target.shape[1]), np.float64)
    for replicate in range(permutations):
        permuted = target[rng.permutation(len(target))]
        null[replicate] = r2_columns(
            permuted,
            ridge_cv_predictions(group_features, permuted, folds, selected_alpha),
        )
    null_q95 = np.quantile(null, 0.95, axis=0)
    excess = observed - np.maximum(null_q95, 0.0)
    eta = np.zeros(target.shape[1], np.float64)
    eta[excess > 0.0] = 0.25
    eta[excess > 0.005] = 0.50
    eta[excess > 0.020] = 0.75
    eta[excess > 0.050] = 1.00
    coef, intercept = ridge_fit(group_features, target, selected_alpha)
    raw_cell = cell_features @ coef + intercept
    predicted_group = measurement @ raw_cell
    # The group-level mapping is used only for within-group contrasts.  This
    # prevents ecological group differences from being copied uncritically to
    # individual cells.
    back_weight = np.asarray(measurement.T.sum(1)).ravel()
    back_pred = np.asarray(measurement.T @ predicted_group)
    back_pred = np.divide(
        back_pred, back_weight[:, None], out=np.zeros_like(back_pred),
        where=back_weight[:, None] > 0,
    )
    contrast = raw_cell - back_pred
    target_sd = np.std(target, axis=0)
    contrast_sd = np.std(contrast, axis=0)
    scale = np.minimum(
        np.divide(
            0.35 * target_sd, np.maximum(contrast_sd, EPS),
            out=np.zeros_like(target_sd), where=contrast_sd > EPS,
        ),
        3.0,
    )
    contrast *= scale[None, :] * eta[None, :]
    return contrast, {
        "selected_alpha": selected_alpha,
        "heldout_r2": observed.tolist(),
        "permutation_q95": null_q95.tolist(),
        "eta": eta.tolist(),
        "contrast_scale": scale.tolist(),
    }


def cell_baseline(measurement: sparse.csr_matrix, target: np.ndarray) -> np.ndarray:
    back_weight = np.asarray(measurement.T.sum(1)).ravel()
    value = np.asarray(measurement.T @ target)
    value = np.divide(
        value, back_weight[:, None], out=np.zeros_like(value),
        where=back_weight[:, None] > 0,
    )
    return value


def build_state_graph(
    coords: np.ndarray,
    features: np.ndarray,
    context: np.ndarray,
    component: np.ndarray,
    neighbors: int,
    radius: float,
):
    n_cells = len(coords)
    if n_cells <= 1:
        return sparse.csr_matrix((n_cells, n_cells)), np.empty((0, 3), np.float64)
    query_k = min(neighbors + 1, n_cells)
    distance, index = cKDTree(coords).query(
        coords, k=query_k, distance_upper_bound=radius, workers=-1,
    )
    if query_k == 1:
        distance = distance[:, None]
        index = index[:, None]
    source = np.repeat(np.arange(n_cells, dtype=np.int32), query_k - 1)
    target = index[:, 1:].ravel().astype(np.int64)
    spatial_distance = distance[:, 1:].ravel()
    valid = (
        np.isfinite(spatial_distance) & (target >= 0) & (target < n_cells)
        & (target != source)
    )
    source = source[valid]
    target = target[valid].astype(np.int32)
    spatial_distance = spatial_distance[valid]
    same_component = component[source] == component[target]
    source = source[same_component]
    target = target[same_component]
    spatial_distance = spatial_distance[same_component]
    if not len(source):
        return sparse.csr_matrix((n_cells, n_cells)), np.empty((0, 3), np.float64)
    feature_distance = np.sqrt(np.mean(np.square(features[source] - features[target]), axis=1))
    if context.shape[1]:
        context_distance = np.sqrt(np.mean(np.square(context[source] - context[target]), axis=1))
    else:
        context_distance = np.zeros(len(source), np.float64)

    def positive_median(values: np.ndarray, fallback: float = 1.0) -> float:
        values = values[np.isfinite(values) & (values > 0)]
        return float(np.median(values)) if len(values) else fallback

    spatial_scale = positive_median(spatial_distance)
    feature_scale = positive_median(feature_distance)
    context_scale = positive_median(context_distance)
    weight = np.exp(
        -0.5 * np.square(spatial_distance / max(spatial_scale, EPS))
        -0.25 * np.square(feature_distance / max(feature_scale, EPS))
        -0.50 * np.square(context_distance / max(context_scale, EPS))
    )
    directed = sparse.coo_matrix(
        (weight, (source, target)), shape=(n_cells, n_cells),
    ).tocsr()
    directed.sum_duplicates()
    mutual = directed.minimum(directed.T).tocsr()
    degree = np.asarray(mutual.sum(1)).ravel()
    isolated = degree <= 0
    if np.any(isolated):
        union = directed.maximum(directed.T).tocsr()
        mutual = (mutual + sparse.diags(isolated.astype(np.float64)) @ union).maximum(
            (mutual + sparse.diags(isolated.astype(np.float64)) @ union).T
        ).tocsr()
    mutual.setdiag(0)
    mutual.eliminate_zeros()
    degree = np.asarray(mutual.sum(1)).ravel()
    inv_sqrt = np.divide(
        1.0, np.sqrt(degree), out=np.zeros_like(degree), where=degree > 0,
    )
    normalized = sparse.diags(inv_sqrt) @ mutual @ sparse.diags(inv_sqrt)
    diagonal = (degree > 0).astype(np.float64)
    laplacian = sparse.diags(diagonal) - normalized
    upper = sparse.triu(mutual, k=1).tocoo()
    edges = np.column_stack((upper.row, upper.col, upper.data)).astype(np.float64)
    audit = {
        "n_edges": int(len(upper.data)),
        "n_isolated": int(np.sum(degree <= 0)),
        "degree_median": float(np.median(np.diff(mutual.indptr))),
        "spatial_scale": spatial_scale,
        "feature_scale": feature_scale,
        "context_scale": context_scale,
    }
    return laplacian.tocsr(), edges, audit


def solve_class_state(
    data: dict, class_value: int, args: argparse.Namespace,
):
    class_cells, class_groups, measurement = class_measurement_matrix(data, class_value)
    k = int(data["program_count"][class_value])
    if not len(class_cells) or not len(class_groups) or k <= 0:
        return class_cells, np.zeros((len(class_cells), MAX_PROGRAMS), np.float32), [], {
            "class": CLASS_NAMES[class_value],
            "n_cells": int(len(class_cells)),
            "n_groups": int(len(class_groups)),
            "program_count": k,
            "state_gate": "class_global",
        }
    target = data["group_scores"][class_groups, :k].astype(np.float64)
    group_features = measurement @ data["features"][class_cells]
    folds = spatial_fold_ids(data["spots"], data["group_bin"][class_groups])
    contrast, feature_audit = feature_prior(
        measurement, group_features, data["features"][class_cells], target, folds,
        args.ridge_permutations, args.seed + class_value * 1000,
    )
    baseline = cell_baseline(measurement, target)
    prior = baseline + contrast
    laplacian, local_edges, graph_audit = build_state_graph(
        data["coords"][class_cells], data["features"][class_cells], baseline,
        data["cell_component"][class_cells], args.graph_neighbors,
        args.graph_radius_array_units,
    )
    row_sq = np.asarray(measurement.multiply(measurement).sum(1)).ravel()
    effective_n = np.divide(1.0, row_sq, out=np.ones_like(row_sq), where=row_sq > EPS)
    effective_n = np.clip(effective_n, 1.0, 100.0)
    weighted_measurement = sparse.diags(effective_n)
    measurement_precision = measurement.T @ weighted_measurement @ measurement
    positive_diag = measurement_precision.diagonal()
    positive_diag = positive_diag[positive_diag > 0]
    base_scale = float(np.median(positive_diag)) if len(positive_diag) else 1.0
    lambda_feature = args.lambda_feature_multiplier * base_scale
    lambda_graph = args.lambda_graph_multiplier * base_scale
    system = (
        measurement_precision
        + lambda_feature * sparse.eye(len(class_cells), format="csr")
        + lambda_graph * laplacian
    ).tocsr()
    rhs = measurement.T @ (effective_n[:, None] * target) + lambda_feature * prior
    state = np.zeros((len(class_cells), k), np.float64)
    cg_status = []
    for program in range(k):
        solution, info = splinalg.cg(
            system, rhs[:, program], x0=prior[:, program],
            rtol=1.0e-6, atol=0.0, maxiter=600,
        )
        if info != 0:
            log(
                f"WARNING {CLASS_NAMES[class_value]} program {program}: "
                f"CG status {info}"
            )
        state[:, program] = solution
        cg_status.append(int(info))
    predicted_group = measurement @ state
    # Preserve each program's weighted global mean while retaining the graph
    # solution.  Extreme inverse values are winsorised only beyond the observed
    # group-state support.
    for program in range(k):
        offset = np.average(target[:, program], weights=effective_n) - np.average(
            predicted_group[:, program], weights=effective_n,
        )
        state[:, program] += offset
        lo, hi = np.quantile(target[:, program], [0.005, 0.995])
        width = max(float(hi - lo), 0.25)
        state[:, program] = np.clip(state[:, program], lo - width, hi + width)
    predicted_group = measurement @ state
    reconstruction_r2 = r2_columns(target, predicted_group)
    within_group_variance = []
    for program in range(k):
        variance = []
        for row in range(measurement.shape[0]):
            start, end = measurement.indptr[row], measurement.indptr[row + 1]
            cells = measurement.indices[start:end]
            weights = measurement.data[start:end]
            if len(cells) > 1:
                mean = float(np.dot(weights, state[cells, program]))
                variance.append(float(np.dot(weights, np.square(state[cells, program] - mean))))
        within_group_variance.append(float(np.median(variance)) if variance else 0.0)
    output = np.zeros((len(class_cells), MAX_PROGRAMS), np.float32)
    output[:, :k] = state.astype(np.float32)
    audit = {
        "class": CLASS_NAMES[class_value],
        "n_cells": int(len(class_cells)),
        "n_groups": int(len(class_groups)),
        "program_count": k,
        "effective_n_median": float(np.median(effective_n)),
        "lambda_base_scale": base_scale,
        "lambda_feature": lambda_feature,
        "lambda_graph": lambda_graph,
        "group_reconstruction_r2": reconstruction_r2.tolist(),
        "within_group_state_variance_median": within_group_variance,
        "cg_status": cg_status,
        "feature_mapping": feature_audit,
        "graph": graph_audit,
        "state_gate": (
            "morphology_refined"
            if np.any(np.asarray(feature_audit["eta"]) > 0)
            else "bin_type_state"
        ),
    }
    global_edges = []
    if len(local_edges):
        global_edges = np.column_stack((
            class_cells[local_edges[:, 0].astype(np.int64)],
            class_cells[local_edges[:, 1].astype(np.int64)],
            local_edges[:, 2],
            np.full(len(local_edges), class_value, np.float64),
        ))
    return class_cells, output, global_edges, audit


def write_state(args: argparse.Namespace, data: dict) -> dict:
    partial = args.state_output.with_suffix(args.state_output.suffix + ".partial")
    if args.state_output.exists() or partial.exists():
        raise FileExistsError(f"Refusing to overwrite state output: {args.state_output}")
    args.state_output.parent.mkdir(parents=True, exist_ok=True)
    state = np.zeros((len(data["class_id"]), MAX_PROGRAMS), np.float32)
    all_edges = []
    audits = []
    for class_value in range(len(CLASS_NAMES)):
        log(f"Solving global cell state for {CLASS_NAMES[class_value]}")
        cells, values, edges, audit = solve_class_state(data, class_value, args)
        state[cells] = values
        if len(edges):
            all_edges.append(edges)
        audits.append(audit)
        log(
            f"{CLASS_NAMES[class_value]}: cells={len(cells):,}, "
            f"groups={audit['n_groups']:,}, gate={audit['state_gate']}"
        )
    edge_values = (
        np.concatenate(all_edges, axis=0)
        if all_edges else np.empty((0, 4), np.float64)
    )
    summary = {
            "schema": "virtual_visium55_spot2cell_mv3_state_v1",
        "algorithm": (
            "same-class spatial-morphology-context graph + global "
            "graph-regularised cell-state inverse"
        ),
        "mv2_model_reused_only_for": [
            "independently fitted class_gene_probability",
            "independently fitted type-specific program_loadings",
            "spot-by-type program activity targets",
            "spot-cell measurement geometry",
        ],
        "fixed_44_programs_used": False,
        "external_labels_used": False,
        "scrna_reference_used": False,
        "visium2_expression_used_for_fit": False,
        "n_cells": int(len(state)),
        "n_state_graph_edges": int(len(edge_values)),
        "classes": audits,
        "type_program_refit": data.get(
            "program_refit",
            {
                "schema": "MV2_program_reuse_diagnostic",
                "max_programs": OLD_MAX_PROGRAMS,
                "selected_program_count": data["program_count"].astype(int).tolist(),
            },
        ),
        "parameters": {
            "graph_neighbors": args.graph_neighbors,
            "graph_radius_array_units": args.graph_radius_array_units,
            "ridge_permutations": args.ridge_permutations,
            "lambda_feature_multiplier": args.lambda_feature_multiplier,
            "lambda_graph_multiplier": args.lambda_graph_multiplier,
            "seed": args.seed,
        },
        "input_provenance": {
            "mv2_h5ad": {
                "path": str(args.mv2_h5ad), "sha256": sha256(args.mv2_h5ad),
            },
            "mv2_allocation": {
                "path": str(args.mv2_allocation), "sha256": sha256(args.mv2_allocation),
            },
        },
    }
    with h5py.File(partial, "w") as h:
        h.attrs["schema"] = summary["schema"]
        h.attrs["complete"] = 0
        h.attrs["summary_json"] = json.dumps(summary, ensure_ascii=False)
        h.create_dataset("cell_id", data=data["cell_id"], compression="lzf")
        h.create_dataset("class_id", data=data["class_id"].astype(np.uint8), compression="lzf")
        h.create_dataset("program_activity", data=state, compression="lzf")
        model = h.create_group("type_program_model")
        model.create_dataset(
            "program_count_by_class", data=data["program_count"].astype(np.int32),
        )
        model.create_dataset(
            "program_loadings", data=data["loadings"].astype(np.float32),
            compression="lzf",
        )
        model.create_dataset(
            "bin_type_bin_index", data=data["group_bin"].astype(np.int32),
            compression="lzf",
        )
        model.create_dataset(
            "bin_type_class_id", data=data["group_class"].astype(np.uint8),
            compression="lzf",
        )
        model.create_dataset(
            "bin_type_program_activity", data=data["group_scores"].astype(np.float32),
            compression="lzf",
        )
        if "feature_gene_coef" in data:
            model.create_dataset(
                "feature_gene_coef", data=data["feature_gene_coef"].astype(np.float32),
                compression="lzf",
            )
            model.create_dataset(
                "feature_gene_clip", data=data["feature_gene_clip"].astype(np.float32),
                compression="lzf",
            )
        graph = h.create_group("same_class_state_graph")
        graph.create_dataset("source_cell_index", data=edge_values[:, 0].astype(np.int32), compression="lzf")
        graph.create_dataset("target_cell_index", data=edge_values[:, 1].astype(np.int32), compression="lzf")
        graph.create_dataset("weight", data=edge_values[:, 2].astype(np.float32), compression="lzf")
        graph.create_dataset("class_id", data=edge_values[:, 3].astype(np.uint8), compression="lzf")
        h.attrs["complete"] = 1
        h.flush()
    os.replace(partial, args.state_output)
    return summary


def load_state(args: argparse.Namespace, data: dict):
    with h5py.File(args.state_output, "r") as h:
        if int(h.attrs.get("complete", 0)) != 1:
            raise RuntimeError("State file is incomplete")
        if not np.array_equal(h["cell_id"][:], data["cell_id"]):
            raise RuntimeError("State file cell order differs from MV2 cell order")
        state = h["program_activity"][:].astype(np.float64)
        if "type_program_model" in h:
            model = h["type_program_model"]
            data["program_count"] = model["program_count_by_class"][:].astype(np.int32)
            data["loadings"] = model["program_loadings"][:].astype(np.float64)
            data["group_bin"] = model["bin_type_bin_index"][:].astype(np.int32)
            data["group_class"] = model["bin_type_class_id"][:].astype(np.int32)
            data["group_scores"] = model["bin_type_program_activity"][:].astype(np.float64)
            if "feature_gene_coef" in model:
                data["feature_gene_coef"] = model["feature_gene_coef"][:].astype(np.float64)
                data["feature_gene_clip"] = model["feature_gene_clip"][:].astype(np.float64)
            else:
                data["feature_gene_coef"] = np.zeros(
                    (len(CLASS_NAMES), data["features"].shape[1], len(data["gene_names"])),
                    np.float64,
                )
                data["feature_gene_clip"] = np.zeros(
                    (len(CLASS_NAMES), len(data["gene_names"])), np.float64,
                )
            data["cell_feature_basis"] = centred_cell_feature_basis(data)
        summary = json.loads(h.attrs["summary_json"])
    return state, summary


def allocate_bin_mv3(
    matrix: sparse.csc_matrix,
    bin_index: int,
    geometry: sparse.csr_matrix,
    capacity: np.ndarray,
    class_id: np.ndarray,
    state: np.ndarray,
    profiles: np.ndarray,
    loadings: np.ndarray,
    program_count: np.ndarray,
    ipf_strength: float,
    ipf_iterations: int,
    cell_feature_basis: np.ndarray | None = None,
    feature_gene_coef: np.ndarray | None = None,
    feature_gene_clip: np.ndarray | None = None,
):
    gs, ge = int(geometry.indptr[bin_index]), int(geometry.indptr[bin_index + 1])
    ms, me = int(matrix.indptr[bin_index]), int(matrix.indptr[bin_index + 1])
    cells = geometry.indices[gs:ge]
    overlap = geometry.data[gs:ge].astype(np.float64)
    genes = matrix.indices[ms:me].astype(np.int32, copy=False)
    counts = matrix.data[ms:me].astype(np.float64)
    if not len(cells) or not len(genes):
        return cells, genes, counts, np.empty((len(cells), len(genes)), np.float32)
    link_capacity = overlap * capacity[cells]
    local_class = class_id[cells]
    propensity = np.zeros((len(cells), len(genes)), np.float64)
    for class_value in np.unique(local_class):
        selected = local_class == class_value
        k = int(program_count[class_value])
        modifier = np.zeros((int(selected.sum()), len(genes)), np.float64)
        if k > 0:
            modifier = state[cells[selected], :k] @ loadings[class_value, :k][:, genes]
        if cell_feature_basis is not None and feature_gene_coef is not None:
            feature_modifier = (
                cell_feature_basis[cells[selected]].astype(np.float64)
                @ feature_gene_coef[class_value][:, genes]
            )
            if feature_gene_clip is not None:
                limit = feature_gene_clip[class_value, genes][None, :]
                feature_modifier = np.clip(
                    feature_modifier, -limit, limit,
                )
            modifier += feature_modifier
        propensity[selected] = (
            link_capacity[selected, None]
            * np.maximum(profiles[class_value, genes][None, :], EPS)
            * np.exp(np.clip(modifier, -6.0, 6.0))
        )
    denominator = propensity.sum(0)
    missing = denominator <= EPS
    if np.any(missing):
        propensity[:, missing] = link_capacity[:, None]
        denominator[missing] = max(float(link_capacity.sum()), EPS)
    q = propensity / np.maximum(denominator[None, :], EPS)
    if ipf_strength > 0 and ipf_iterations > 0:
        target_library = counts.sum() * link_capacity / max(float(link_capacity.sum()), EPS)
        for _ in range(ipf_iterations):
            assigned = q @ counts
            scale = np.power(
                np.divide(
                    target_library, np.maximum(assigned, EPS),
                    out=np.ones_like(target_library), where=assigned > 0,
                ),
                ipf_strength,
            )
            scale = np.clip(scale, 0.25, 4.0)
            q *= scale[:, None]
            q /= np.maximum(q.sum(0, keepdims=True), EPS)
    q32 = q.astype(np.float32)
    pivot = np.argmax(q, axis=0)
    correction = 1.0 - q32.sum(0, dtype=np.float64)
    q32[pivot, np.arange(len(genes))] = (
        q32[pivot, np.arange(len(genes))].astype(np.float64) + correction
    ).astype(np.float32)
    return cells, genes, counts, q32


def write_allocation(args: argparse.Namespace, data: dict, state: np.ndarray, state_summary: dict):
    partial = args.allocation_output.with_suffix(args.allocation_output.suffix + ".partial")
    if args.allocation_output.exists() or partial.exists():
        raise FileExistsError(f"Refusing to overwrite allocation: {args.allocation_output}")
    args.allocation_output.parent.mkdir(parents=True, exist_ok=True)
    matrix = data["matrix"]
    geometry = data["geometry"]
    n_links = np.diff(geometry.indptr).astype(np.int64)
    bin_nnz = np.diff(matrix.indptr).astype(np.int64)
    record_indptr = np.zeros(matrix.shape[1] + 1, np.int64)
    np.cumsum(n_links * bin_nnz, out=record_indptr[1:])
    occupied = n_links > 0
    assigned_gene = np.zeros(matrix.shape[0], np.float64)
    occupied_gene = np.asarray(matrix @ occupied.astype(np.float64)).ravel()
    unassigned_gene = np.asarray(matrix @ (~occupied).astype(np.float64)).ravel()
    cell_total = np.zeros(geometry.shape[1], np.float64)
    max_bin_gene_error = 0.0
    max_library_log_ratio = 0.0
    summary = {
        **state_summary,
        "schema": "virtual_visium55_spot2cell_mv3_summary_v1",
        "allocation": (
            "joint cell-by-gene KL prior from A*u*P_state with exact spot-gene "
            "projection and soft cell-library IPF"
        ),
        "library_ipf_strength": args.library_ipf_strength,
        "library_ipf_iterations": args.library_ipf_iterations,
        "primary_expression_for_clustering": "P_state factorization / P_state_278 panel",
        "X_mass_resolution_gate": "conservation_only",
        "X_mass_warning": (
            "Exact spot-gene conservation necessarily retains V55 capture-spot imprint; "
            "do not use X_mass log1p_cp10k as the primary clustering representation."
        ),
    }
    with h5py.File(args.mv2_allocation, "r") as old, h5py.File(partial, "w") as h:
        for group_name in ("observed", "geometry", "source_016um_to_virtual_visium55"):
            if group_name in old:
                old.copy(group_name, h)
        h.attrs.update({
            "schema": "visium_hd_virtual_visium55_cell_allocation_mv3_v1",
            "complete": 0,
            "sample": args.sample,
            "expression_units": "fractional_captured_counts",
            "geometry": "all-real-nucleus within-component Voronoi overlap",
            "allocation_model": "joint_KL_P_state_exact_gene_projection_soft_library_IPF",
            "P_state_semantics": "continuous library-free expression probability",
            "X_mass_semantics": "exact spot-gene-conserved fractional count mass",
            "summary_json": json.dumps(summary, ensure_ascii=False),
        })
        allocation = h.create_group("allocation")
        allocation.attrs["layout"] = (
            "For bin b, reshape responsibility[bin_record_indptr[b]:"
            "bin_record_indptr[b+1]] as (observed genes in bin, geometry cells in bin)."
        )
        allocation.create_dataset("bin_record_indptr", data=record_indptr, compression="lzf")
        total_records = int(record_indptr[-1])
        responsibility_ds = allocation.create_dataset(
            "responsibility", shape=(total_records,), dtype=np.float32,
            chunks=(min(1_000_000, max(total_records, 1)),), compression="lzf",
        )
        last = time.time()
        log(f"Writing {total_records:,} MV3 allocation responsibilities")
        for b in np.flatnonzero(occupied):
            cells, genes, counts, responsibility = allocate_bin_mv3(
                matrix, int(b), geometry, data["capacity"], data["class_id"],
                state, data["profiles"], data["loadings"], data["program_count"],
                args.library_ipf_strength, args.library_ipf_iterations,
                data.get("cell_feature_basis"), data.get("feature_gene_coef"),
                data.get("feature_gene_clip"),
            )
            start, end = int(record_indptr[b]), int(record_indptr[b + 1])
            responsibility_ds[start:end] = responsibility.T.ravel()
            mass = responsibility.astype(np.float64) * counts[None, :]
            reconstructed = mass.sum(0)
            max_bin_gene_error = max(
                max_bin_gene_error,
                float(np.max(np.abs(reconstructed - counts))) if len(counts) else 0.0,
            )
            np.add.at(assigned_gene, genes, reconstructed)
            link_total = mass.sum(1)
            np.add.at(cell_total, cells, link_total)
            target_weight = (
                geometry.data[geometry.indptr[b]:geometry.indptr[b + 1]]
                * data["capacity"][cells]
            )
            target = (
                counts.sum() * target_weight
                / max(float(target_weight.sum()), EPS)
            )
            valid = (target > EPS) & (link_total > EPS)
            if np.any(valid):
                max_library_log_ratio = max(
                    max_library_log_ratio,
                    float(np.max(np.abs(np.log(link_total[valid] / target[valid])))),
                )
            if time.time() - last > 60:
                log(f"Allocated {b + 1:,}/{matrix.shape[1]:,} virtual spots")
                last = time.time()
        audit = h.create_group("audit")
        audit.create_dataset("observed_gene_counts_occupied_bins", data=occupied_gene, compression="lzf")
        audit.create_dataset("assigned_gene_counts", data=assigned_gene, compression="lzf")
        audit.create_dataset("unassigned_gene_counts_zero_cell_bins", data=unassigned_gene, compression="lzf")
        audit.create_dataset("unassigned_bin_index", data=np.flatnonzero(~occupied), compression="lzf")
        audit.attrs.update({
            "max_abs_gene_conservation_error": float(np.max(np.abs(assigned_gene - occupied_gene))),
            "max_abs_bin_gene_float32_error": max_bin_gene_error,
            "observed_total_counts": float(matrix.sum()),
            "observed_counts_occupied_bins": float(np.sum(occupied_gene)),
            "unassigned_counts_zero_cell_bins": float(np.sum(unassigned_gene)),
            "assigned_cell_counts": float(np.sum(cell_total)),
            "max_abs_cell_library_log_ratio_to_capacity_target": max_library_log_ratio,
        })
        state_group = h.create_group("cell_state")
        state_group.create_dataset("program_activity", data=state.astype(np.float32), compression="lzf")
        with h5py.File(args.state_output, "r") as state_file:
            state_file.copy("same_class_state_graph", state_group)
        h.attrs["complete"] = 1
        h.flush()
    os.replace(partial, args.allocation_output)
    return cell_total, {
        "observed_total_counts": float(matrix.sum()),
        "observed_counts_occupied_bins": float(np.sum(occupied_gene)),
        "unassigned_counts_zero_cell_bins": float(np.sum(unassigned_gene)),
        "assigned_cell_counts": float(np.sum(cell_total)),
        "max_abs_gene_error": float(np.max(np.abs(assigned_gene - occupied_gene))),
        "max_abs_bin_gene_float32_error": max_bin_gene_error,
        "allocation_records": int(record_indptr[-1]),
        "max_abs_cell_library_log_ratio_to_capacity_target": max_library_log_ratio,
    }, summary


def append_dataset(dataset: h5py.Dataset, values: np.ndarray) -> None:
    start = dataset.shape[0]
    dataset.resize((start + len(values),))
    dataset[start:] = values


def encoded_array(group: h5py.Group, name: str, values, **kwargs):
    dataset = group.create_dataset(name, data=values, **kwargs)
    dataset.attrs["encoding-type"] = "array"
    dataset.attrs["encoding-version"] = "0.2.0"
    return dataset


def encoded_scalar_string(group: h5py.Group, name: str, value: str):
    dataset = group.create_dataset(
        name, data=np.asarray(value, dtype=h5py.string_dtype("utf-8")),
    )
    dataset.attrs["encoding-type"] = "string"
    dataset.attrs["encoding-version"] = "0.2.0"
    return dataset


def write_text_dataset(group: h5py.Group, name: str, values):
    dataset = group.create_dataset(
        name, data=np.asarray(values, dtype=object), dtype=h5py.string_dtype("utf-8"),
    )
    dataset.attrs["encoding-type"] = "string-array"
    dataset.attrs["encoding-version"] = "0.2.0"
    return dataset


def write_categorical(group: h5py.Group, name: str, codes: np.ndarray, categories: list[str]):
    item = group.create_group(name)
    item.attrs["encoding-type"] = "categorical"
    item.attrs["encoding-version"] = "0.2.0"
    item.attrs["ordered"] = False
    encoded_array(item, "codes", codes.astype(np.int8), compression="lzf")
    categories_ds = item.create_dataset(
        "categories", data=np.asarray(categories, dtype=object),
        dtype=h5py.string_dtype("utf-8"),
    )
    categories_ds.attrs["encoding-type"] = "string-array"
    categories_ds.attrs["encoding-version"] = "0.2.0"


def initialise_h5ad(
    args: argparse.Namespace, data: dict, state: np.ndarray,
    cell_total: np.ndarray, summary: dict,
):
    partial = args.output.with_suffix(args.output.suffix + ".partial")
    if args.output.exists() or partial.exists():
        raise FileExistsError(f"Refusing to overwrite H5AD: {args.output}")
    with h5py.File(args.mv2_h5ad, "r") as old, h5py.File(partial, "w") as h:
        for key, value in old.attrs.items():
            h.attrs[key] = value
        h.attrs["complete"] = 0
        h.attrs["schema"] = "visium_hd_virtual_visium55_spot2cell_mv3_v1"
        h.attrs["X_semantics"] = "X_mass exact spot-gene-conserved fractional counts"
        h.attrs["P_state_semantics"] = (
            "factorised continuous expression probability in uns/mv3_state_factorization"
        )
        for group_name in ("obs", "var", "obsm", "varm", "uns", "obsp", "varp"):
            old.copy(group_name, h)
        for name in ("allocated_count_total", "rna_total"):
            del h[f"obs/{name}"]
            encoded_array(h["obs"], name, cell_total.astype(np.float64), compression="lzf")
        old_program = h["obsm/program_activity"][:]
        del h["obsm/program_activity"]
        encoded_array(h["obsm"], "program_activity_mv2", old_program, compression="lzf")
        encoded_array(h["obsm"], "program_activity", state.astype(np.float32), compression="lzf")
        encoded_array(h["obsm"], "P_state_program_activity", state.astype(np.float32), compression="lzf")
        if "program_loadings" in h["varm"]:
            del h["varm/program_loadings"]
        flattened_loadings = np.transpose(
            data["loadings"], (2, 0, 1)
        ).reshape(len(data["gene_names"]), len(CLASS_NAMES) * MAX_PROGRAMS)
        encoded_array(
            h["varm"], "program_loadings", flattened_loadings.astype(np.float32),
            compression="lzf",
        )
        if "program_count_by_class" in h["uns"]:
            del h["uns/program_count_by_class"]
        encoded_array(
            h["uns"], "program_count_by_class",
            data["program_count"].astype(np.int32),
        )
        for name, values in (
            ("bin_type_bin_index", data["group_bin"].astype(np.int32)),
            ("bin_type_class_id", data["group_class"].astype(np.uint8)),
            ("bin_type_program_activity", data["group_scores"].astype(np.float32)),
        ):
            if name in h["uns"]:
                del h[f"uns/{name}"]
            encoded_array(h["uns"], name, values, compression="lzf")
        if "program_axis" in h["uns"]:
            del h["uns/program_axis"]
        write_text_dataset(
            h["uns"], "program_axis", [f"program_{index}" for index in range(MAX_PROGRAMS)],
        )
        if "cell_state_modeled" in h["var"]:
            del h["var/cell_state_modeled"]
        encoded_array(
            h["var"], "cell_state_modeled",
            np.any(data["loadings"] != 0, axis=(0, 1)), compression="lzf",
        )
        for name in ("state_resolution", "post_allocation_resolution"):
            if name in h["obs"]:
                del h[f"obs/{name}"]
        # P_state and X_mass have deliberately separate evidence contracts.
        gates = []
        for item in summary["classes"]:
            gates.append(item["state_gate"])
        categories = ["class_global", "bin_type_state", "morphology_refined"]
        codes = np.asarray([categories.index(gates[c]) for c in data["class_id"]], np.int8)
        write_categorical(h["obs"], "state_resolution", codes, categories)
        write_categorical(
            h["obs"], "post_allocation_resolution",
            np.zeros(len(codes), np.int8), ["conservation_only"],
        )
        order = list(h["obs"].attrs.get("column-order", []))
        existing_order = [
            item.decode() if isinstance(item, bytes) else str(item) for item in order
        ]
        for name in ("state_resolution", "post_allocation_resolution"):
            if name not in existing_order:
                order.append(name)
        h["obs"].attrs["column-order"] = np.asarray(
            order, dtype=h5py.string_dtype("utf-8")
        )
        if "summary_json_mv2" not in h["uns"]:
            old_summary = h["uns/summary_json"][()]
            if isinstance(old_summary, bytes):
                old_summary = old_summary.decode("utf-8")
            encoded_scalar_string(h["uns"], "summary_json_mv2", str(old_summary))
        del h["uns/summary_json"]
        encoded_scalar_string(h["uns"], "summary_json", json.dumps(summary, ensure_ascii=False))
        factor = h["uns"].create_group("mv3_state_factorization")
        factor.attrs["encoding-type"] = "dict"
        factor.attrs["encoding-version"] = "0.1.0"
        factor.attrs["formula"] = "softmax(log(class_gene_probability)+z@program_loadings)"
        encoded_array(factor, "program_activity", state.astype(np.float32), compression="lzf")
        encoded_array(factor, "program_count_by_class", data["program_count"].astype(np.int32))
        if "feature_gene_coef" in data:
            encoded_array(
                factor, "feature_gene_coef",
                data["feature_gene_coef"].astype(np.float32), compression="lzf",
            )
            encoded_array(
                factor, "feature_gene_clip",
                data["feature_gene_clip"].astype(np.float32), compression="lzf",
            )
            factor.attrs["feature_formula"] = (
                "locally_centered_standardized_nucleus_features @ gated_gene_coefficients"
            )
        # X and the normalised layer are added in a second streaming pass.
        h.flush()
    return partial


def build_cell_csr(
    partial_h5ad: Path, allocation_path: Path, data: dict,
    cell_batch_size: int,
):
    matrix = data["matrix"]
    geometry = data["geometry"]
    n_cells = geometry.shape[1]
    n_genes = matrix.shape[0]
    with h5py.File(partial_h5ad, "r+") as h, h5py.File(allocation_path, "r") as allocation:
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
        x_indptr = x.create_dataset(
            "indptr", shape=(n_cells + 1,), dtype=np.int64, compression="lzf",
        )
        x_indptr[0] = 0
        layers = h.create_group("layers")
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
        for cell_start in range(0, n_cells, cell_batch_size):
            cell_end = min(cell_start + cell_batch_size, n_cells)
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
                if not np.any(selected):
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
                    (
                        np.concatenate(data_parts),
                        (np.concatenate(row_parts), np.concatenate(col_parts)),
                    ),
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
                10000.0, row_total, out=np.zeros_like(row_total), where=row_total > 0,
            )
            normalized_values = np.log1p(
                batch.data.astype(np.float64) * np.repeat(row_scale, np.diff(batch.indptr))
            ).astype(np.float32)
            append_dataset(normalized_data, normalized_values)
            cursor += batch.nnz
            if time.time() - last > 60:
                log(f"Built MV3 cell CSR through {cell_end:,}/{n_cells:,}; nnz={cursor:,}")
                last = time.time()
        normalized["indices"] = x_indices
        normalized["indptr"] = x_indptr
        h.attrs["complete"] = 1
        h.flush()
    return int(cursor)


def materialize_pstate_panel(args: argparse.Namespace, data: dict, state: np.ndarray):
    partial = args.pstate_panel_output.with_suffix(args.pstate_panel_output.suffix + ".partial")
    if args.pstate_panel_output.exists() or partial.exists():
        raise FileExistsError(f"Refusing to overwrite P_state panel: {args.pstate_panel_output}")
    with h5py.File(args.panel_reference, "r") as reference:
        panel_names = read_var_names(reference)
    gene_lookup = {str(gene): i for i, gene in enumerate(data["gene_names"])}
    panel_index = np.asarray([gene_lookup.get(str(gene), -1) for gene in panel_names], np.int32)
    present = panel_index >= 0
    if not np.all(present):
        missing = panel_names[~present]
        raise RuntimeError(f"MV2 full gene space lacks panel genes: {missing[:10].tolist()}")
    n_cells = len(data["class_id"])
    n_panel = len(panel_names)
    with h5py.File(partial, "w") as h:
        h.attrs["encoding-type"] = "anndata"
        h.attrs["encoding-version"] = "0.1.0"
        h.attrs["schema"] = "visium55_mv3_P_state_panel_v1"
        h.attrs["complete"] = 0
        obs = h.create_group("obs")
        obs.attrs["encoding-type"] = "dataframe"
        obs.attrs["encoding-version"] = "0.2.0"
        obs.attrs["_index"] = "_index"
        obs.attrs["column-order"] = np.asarray(
            ["cell_id", "cellvit_class_id", "cellvit_label", "owner_bin_index"],
            dtype=h5py.string_dtype("utf-8"),
        )
        obs_index = np.asarray([f"{args.sample}_{value}" for value in data["cell_id"]], dtype=object)
        index_ds = obs.create_dataset("_index", data=obs_index, dtype=h5py.string_dtype("utf-8"))
        index_ds.attrs["encoding-type"] = "string-array"
        index_ds.attrs["encoding-version"] = "0.2.0"
        encoded_array(obs, "cell_id", data["cell_id"], compression="lzf")
        encoded_array(obs, "cellvit_class_id", data["class_id"].astype(np.int32), compression="lzf")
        encoded_array(obs, "owner_bin_index", data["owner_bin"].astype(np.int32), compression="lzf")
        write_categorical(obs, "cellvit_label", data["class_id"].astype(np.int8), CLASS_NAMES)
        var = h.create_group("var")
        var.attrs["encoding-type"] = "dataframe"
        var.attrs["encoding-version"] = "0.2.0"
        var.attrs["_index"] = "_index"
        var.attrs["column-order"] = np.asarray(["gene_name"], dtype=h5py.string_dtype("utf-8"))
        var_index = var.create_dataset(
            "_index", data=np.asarray(panel_names, dtype=object), dtype=h5py.string_dtype("utf-8"),
        )
        var_index.attrs["encoding-type"] = "string-array"
        var_index.attrs["encoding-version"] = "0.2.0"
        gene_ds = var.create_dataset(
            "gene_name", data=np.asarray(panel_names, dtype=object), dtype=h5py.string_dtype("utf-8"),
        )
        gene_ds.attrs["encoding-type"] = "string-array"
        gene_ds.attrs["encoding-version"] = "0.2.0"
        obsm = h.create_group("obsm")
        obsm.attrs["encoding-type"] = "dict"
        obsm.attrs["encoding-version"] = "0.1.0"
        encoded_array(obsm, "spatial", data["coords"].astype(np.float32), compression="lzf")
        encoded_array(obsm, "program_activity", state.astype(np.float32), compression="lzf")
        for name in ("obsp", "varm", "varp", "layers"):
            group = h.create_group(name)
            group.attrs["encoding-type"] = "dict"
            group.attrs["encoding-version"] = "0.1.0"
        uns = h.create_group("uns")
        uns.attrs["encoding-type"] = "dict"
        uns.attrs["encoding-version"] = "0.1.0"
        encoded_scalar_string(
            uns, "P_state_semantics",
            "continuous library-free expression probability; X=log1p(1e4*P_state)",
        )
        x = h.create_dataset(
            "X", shape=(n_cells, n_panel), dtype=np.float32,
            chunks=(min(2048, n_cells), n_panel), compression="lzf",
        )
        x.attrs["encoding-type"] = "array"
        x.attrs["encoding-version"] = "0.2.0"
        for start in range(0, n_cells, 4096):
            end = min(start + 4096, n_cells)
            values = np.zeros((end - start, n_panel), np.float64)
            local_class = data["class_id"][start:end]
            for class_value in np.unique(local_class):
                selected = local_class == class_value
                k = int(data["program_count"][class_value])
                modifier = np.zeros((int(selected.sum()), n_panel), np.float64)
                if k > 0:
                    modifier = (
                        state[start:end][selected, :k]
                        @ data["loadings"][class_value, :k][:, panel_index]
                    )
                if "feature_gene_coef" in data:
                    feature_modifier = (
                        data["cell_feature_basis"][start:end][selected].astype(np.float64)
                        @ data["feature_gene_coef"][class_value][:, panel_index]
                    )
                    limit = data["feature_gene_clip"][class_value, panel_index][None, :]
                    modifier += np.clip(feature_modifier, -limit, limit)
                probability = (
                    np.maximum(data["profiles"][class_value, panel_index][None, :], EPS)
                    * np.exp(np.clip(modifier, -6.0, 6.0))
                )
                # Normalisation is performed over the full gene space, not over
                # the 278-gene panel.
                full_modifier = np.zeros((int(selected.sum()), data["profiles"].shape[1]), np.float64)
                if k > 0:
                    full_modifier = (
                        state[start:end][selected, :k]
                        @ data["loadings"][class_value, :k]
                    )
                if "feature_gene_coef" in data:
                    feature_modifier = (
                        data["cell_feature_basis"][start:end][selected].astype(np.float64)
                        @ data["feature_gene_coef"][class_value]
                    )
                    limit = data["feature_gene_clip"][class_value][None, :]
                    full_modifier += np.clip(feature_modifier, -limit, limit)
                full_denominator = np.sum(
                    np.maximum(data["profiles"][class_value][None, :], EPS)
                    * np.exp(np.clip(full_modifier, -6.0, 6.0)),
                    axis=1,
                )
                probability /= np.maximum(full_denominator[:, None], EPS)
                values[selected] = np.log1p(1.0e4 * probability)
            x[start:end] = values.astype(np.float32)
        h.attrs["complete"] = 1
        h.flush()
    os.replace(partial, args.pstate_panel_output)


def materialize_xmass_panel_precheck(
    args: argparse.Namespace, data: dict, state: np.ndarray,
) -> dict:
    """Run the full-gene MV3 projection but retain only the 278-gene cell mass.

    The IPF library factors are therefore identical to the eventual full output;
    this is a computational precheck, not a reduced-gene re-fit.
    """
    partial = args.xmass_panel_output.with_suffix(
        args.xmass_panel_output.suffix + ".partial"
    )
    if args.xmass_panel_output.exists() or partial.exists():
        raise FileExistsError(
            f"Refusing to overwrite X_mass precheck: {args.xmass_panel_output}"
        )
    if not args.pstate_panel_output.exists():
        materialize_pstate_panel(args, data, state)
    with h5py.File(args.panel_reference, "r") as reference:
        panel_names = read_var_names(reference)
    gene_lookup = {str(gene): index for index, gene in enumerate(data["gene_names"])}
    panel_index = np.asarray([gene_lookup.get(str(gene), -1) for gene in panel_names], np.int32)
    if np.any(panel_index < 0):
        raise RuntimeError("Full V55 gene space lacks one or more panel genes")
    panel_lookup = np.full(len(data["gene_names"]), -1, np.int32)
    panel_lookup[panel_index] = np.arange(len(panel_index), dtype=np.int32)
    mass = np.zeros((len(data["class_id"]), len(panel_index)), np.float32)
    assigned = np.zeros(len(panel_index), np.float64)
    geometry = data["geometry"]
    occupied = np.diff(geometry.indptr) > 0
    observed = np.asarray(
        data["matrix"][panel_index] @ occupied.astype(np.float64)
    ).ravel().astype(np.float64)
    last = time.time()
    for b in np.flatnonzero(occupied):
        cells, genes, counts, responsibility = allocate_bin_mv3(
            data["matrix"], int(b), geometry, data["capacity"], data["class_id"],
            state, data["profiles"], data["loadings"], data["program_count"],
            args.library_ipf_strength, args.library_ipf_iterations,
            data.get("cell_feature_basis"), data.get("feature_gene_coef"),
            data.get("feature_gene_clip"),
        )
        local_panel = panel_lookup[genes]
        keep = local_panel >= 0
        if np.any(keep):
            block = responsibility[:, keep].astype(np.float64) * counts[keep][None, :]
            mass[np.ix_(cells, local_panel[keep])] += block.astype(np.float32)
            np.add.at(assigned, local_panel[keep], block.sum(0))
        if time.time() - last > 60:
            log(f"MV3 full-gene precheck projected {b + 1:,}/{data['matrix'].shape[1]:,} spots")
            last = time.time()
    with h5py.File(args.pstate_panel_output, "r") as old, h5py.File(partial, "w") as h:
        for key, value in old.attrs.items():
            h.attrs[key] = value
        h.attrs["schema"] = "visium55_mv3_X_mass_panel_precheck_v1"
        h.attrs["X_semantics"] = (
            "fractional counts from full-gene MV3 allocation; only 278 genes retained"
        )
        h.attrs["complete"] = 0
        for name in old.keys():
            if name != "X":
                old.copy(name, h)
        x = h.create_dataset(
            "X", data=mass, dtype=np.float32,
            chunks=(min(2048, len(mass)), mass.shape[1]), compression="lzf",
        )
        x.attrs["encoding-type"] = "array"
        x.attrs["encoding-version"] = "0.2.0"
        if "P_state_semantics" in h["uns"]:
            del h["uns/P_state_semantics"]
        encoded_scalar_string(
            h["uns"], "X_mass_semantics",
            "full-gene joint KL/IPF allocation, exact spot-gene projection",
        )
        h.attrs["complete"] = 1
        h.flush()
    os.replace(partial, args.xmass_panel_output)
    return {
        "max_abs_panel_gene_conservation_error": float(np.max(np.abs(assigned - observed))),
        "panel_observed_mass": float(observed.sum()),
        "panel_assigned_mass": float(assigned.sum()),
    }


def atomic_json(value: dict, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    if path.exists() or partial.exists():
        raise FileExistsError(f"Refusing to overwrite QC: {path}")
    with partial.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    os.replace(partial, path)


def finalize_embedded_summary(
    h5ad_path: Path, allocation_path: Path, summary: dict,
) -> None:
    encoded = json.dumps(summary, ensure_ascii=False)
    with h5py.File(h5ad_path, "r+") as h:
        if "summary_json" in h["uns"]:
            del h["uns/summary_json"]
        encoded_scalar_string(h["uns"], "summary_json", encoded)
        h.flush()
    with h5py.File(allocation_path, "r+") as h:
        h.attrs["summary_json"] = encoded
        h.flush()


def main() -> None:
    args = parse_args()
    log(f"Loading frozen {args.sample} V55 MV2 model, geometry, and observed spot matrix")
    data = load_inputs(args)
    if args.mode in {"state", "all"} and not args.state_output.exists():
        if not args.reuse_mv2_programs:
            refit_mv3_programs(args, data)
        state_summary = write_state(args, data)
    else:
        _, state_summary = load_state(args, data)
    if args.mode == "state":
        log(f"Completed MV3 state inverse: {args.state_output}")
        return
    state, state_summary = load_state(args, data)
    if args.mode == "panel":
        if not args.pstate_panel_output.exists():
            materialize_pstate_panel(args, data, state)
        audit = materialize_xmass_panel_precheck(args, data, state)
        log(f"Completed MV3 X_mass precheck: {json.dumps(audit)}")
        return
    cell_total, conservation, summary = write_allocation(
        args, data, state, state_summary,
    )
    summary["count_conservation"] = conservation
    partial_h5ad = initialise_h5ad(args, data, state, cell_total, summary)
    nnz = build_cell_csr(
        partial_h5ad, args.allocation_output, data, args.cell_batch_size,
    )
    os.replace(partial_h5ad, args.output)
    if not args.pstate_panel_output.exists():
        materialize_pstate_panel(args, data, state)
    else:
        log(f"Reusing existing P_state panel: {args.pstate_panel_output}")
    summary["cell_matrix_nnz"] = nnz
    summary["outputs"] = {
        "cell_resolved_h5ad": str(args.output),
        "allocation_h5": str(args.allocation_output),
        "state_h5": str(args.state_output),
        "P_state_278_h5ad": str(args.pstate_panel_output),
    }
    summary["output_bytes"] = {
        str(path.name): int(path.stat().st_size)
        for path in (
            args.output, args.allocation_output, args.state_output,
            args.pstate_panel_output,
        )
    }
    finalize_embedded_summary(args.output, args.allocation_output, summary)
    atomic_json(summary, args.qc_json)
    log(f"Completed {args.sample} V55 MV3: {args.output}")


if __name__ == "__main__":
    main()
