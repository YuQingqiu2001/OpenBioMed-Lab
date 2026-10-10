#!/usr/bin/env python3
"""Core primitives for Spot2Cell-MV4.

MV4 combines a mathematical spot-to-cell inverse with a slice-internal,
multiple-instance weak model.  The weak model predicts only low-dimensional
within-lineage state from CellViT features; it never treats inferred cell
expression as a supervised target.  All gene mass is assigned by a final
per-spot softmax projection, so every observed spot-by-gene value is conserved.

This module is deliberately independent of P1 file paths.  It contains the
small, testable numerical contracts used by the V55 panel prototype and later
by the 16 um and relative-expression branches.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as splinalg
from scipy.spatial import cKDTree


EPS = 1.0e-12


@dataclass
class WeakStateResult:
    """Cross-fitted slice-internal morphology model result."""

    morph_contrast: np.ndarray
    oof_prediction: np.ndarray
    context_prediction: np.ndarray
    folds: np.ndarray
    eta: np.ndarray
    audit: dict


@dataclass
class GraphResult:
    """Boundary-aware graph and its audit information."""

    laplacian: sparse.csr_matrix
    edges: np.ndarray
    boundary_score: np.ndarray
    audit: dict


def _as_2d(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    if values.ndim == 1:
        values = values[:, None]
    if values.ndim != 2:
        raise ValueError("expected a one- or two-dimensional array")
    return values


def contiguous_spatial_folds(coords: np.ndarray, n_folds: int = 5) -> np.ndarray:
    """Split coordinates into deterministic contiguous bands.

    The MV3 modular hex-row fold assignment interleaved neighbouring spots and
    therefore did not provide a genuine spatial holdout.  MV4 projects the
    coordinates onto their first principal spatial axis and cuts the ordered
    positions into contiguous, nearly equal-sized bands.
    """

    coords = _as_2d(coords)
    if coords.shape[1] < 2:
        raise ValueError("spatial coordinates require at least two columns")
    if n_folds < 2:
        raise ValueError("n_folds must be at least two")
    if len(coords) < n_folds:
        raise ValueError("fewer spatial observations than folds")
    centred = coords[:, :2] - np.mean(coords[:, :2], axis=0, keepdims=True)
    if np.allclose(centred, 0.0):
        projection = np.arange(len(coords), dtype=np.float64)
    else:
        _, _, vh = np.linalg.svd(centred, full_matrices=False)
        axis = vh[0]
        # Fix the arbitrary SVD sign for reproducibility.
        pivot = int(np.argmax(np.abs(axis)))
        if axis[pivot] < 0:
            axis = -axis
        projection = centred @ axis
    order = np.argsort(projection, kind="mergesort")
    folds = np.empty(len(coords), np.int32)
    for fold, indices in enumerate(np.array_split(order, n_folds)):
        folds[indices] = fold
    return folds


def spline_feature_basis(
    features: np.ndarray,
    knots: Iterable[float] = (-1.0, 0.0, 1.0),
    clip: float = 6.0,
) -> np.ndarray:
    """Fixed low-complexity hinge basis for robust standardized features."""

    features = np.clip(_as_2d(features), -clip, clip)
    blocks = [features]
    for knot in knots:
        blocks.append(np.maximum(features - float(knot), 0.0))
    basis = np.concatenate(blocks, axis=1)
    # Constant or non-finite columns carry no cell-level information.
    basis = np.nan_to_num(basis, nan=0.0, posinf=clip, neginf=-clip)
    scale = np.std(basis, axis=0)
    return basis[:, scale > 1.0e-8]


def ridge_fit(x: np.ndarray, y: np.ndarray, alpha: float) -> tuple[np.ndarray, np.ndarray]:
    """Multi-output ridge with an unpenalised intercept."""

    x = _as_2d(x)
    y = _as_2d(y)
    if len(x) != len(y):
        raise ValueError("x and y row counts differ")
    if not len(x):
        raise ValueError("cannot fit an empty design")
    if x.shape[1] == 0:
        return np.zeros((0, y.shape[1])), np.mean(y, axis=0)
    mean_x = np.mean(x, axis=0)
    scale_x = np.std(x, axis=0)
    scale_x[scale_x < 1.0e-8] = 1.0
    mean_y = np.mean(y, axis=0)
    z = (x - mean_x) / scale_x
    lhs = z.T @ z + float(alpha) * np.eye(z.shape[1])
    coef_scaled = np.linalg.solve(lhs, z.T @ (y - mean_y))
    coef = coef_scaled / scale_x[:, None]
    intercept = mean_y - mean_x @ coef
    return coef, intercept


def r2_columns(truth: np.ndarray, prediction: np.ndarray) -> np.ndarray:
    truth = _as_2d(truth)
    prediction = _as_2d(prediction)
    denominator = np.sum(np.square(truth - np.mean(truth, axis=0)), axis=0)
    numerator = np.sum(np.square(truth - prediction), axis=0)
    return np.divide(
        denominator - numerator,
        np.maximum(denominator, EPS),
        out=np.full(truth.shape[1], -np.inf),
        where=denominator > EPS,
    )


def _training_mask_with_buffer(
    coords: np.ndarray,
    folds: np.ndarray,
    test_fold: int,
    buffer_distance: float,
) -> tuple[np.ndarray, np.ndarray]:
    test = folds == test_fold
    train = ~test
    if buffer_distance > 0 and np.any(test) and np.any(train):
        distance, _ = cKDTree(coords[test, :2]).query(coords[train, :2], k=1)
        train_indices = np.flatnonzero(train)
        train[train_indices[distance <= buffer_distance]] = False
    return train, test


def _crossfit_group_prediction(
    design: np.ndarray,
    target: np.ndarray,
    coords: np.ndarray,
    folds: np.ndarray,
    alpha: float,
    buffer_distance: float,
) -> tuple[np.ndarray, list[np.ndarray]]:
    prediction = np.zeros_like(target, dtype=np.float64)
    coefficients: list[np.ndarray] = []
    for fold in np.unique(folds):
        train, test = _training_mask_with_buffer(
            coords, folds, int(fold), buffer_distance,
        )
        if not np.any(test):
            continue
        if int(np.sum(train)) <= max(5, design.shape[1] // 4):
            prediction[test] = np.mean(target[train], axis=0) if np.any(train) else 0.0
            coefficients.append(np.zeros((design.shape[1], target.shape[1])))
            continue
        coef, intercept = ridge_fit(design[train], target[train], alpha)
        prediction[test] = design[test] @ coef + intercept
        coefficients.append(coef)
    return prediction, coefficients


def measurement_nullspace_center(
    measurement: sparse.spmatrix,
    values: np.ndarray,
    tolerance: float = 1.0e-9,
) -> tuple[np.ndarray, float]:
    """Project cell values into the null space of the group measurement.

    For every returned column ``v``, ``measurement @ v`` is numerically zero.
    This is stronger than subtracting a back-projected group mean when a nucleus
    overlaps more than one spot.
    """

    measurement = sparse.csr_matrix(measurement, dtype=np.float64)
    values = _as_2d(values)
    if measurement.shape[1] != len(values):
        raise ValueError("measurement and cell values are not aligned")
    if measurement.shape[0] == 0 or measurement.nnz == 0:
        return values.copy(), 0.0
    gram = (measurement @ measurement.T).tocsr()
    centred = values.copy()
    max_error = 0.0
    for column in range(values.shape[1]):
        rhs = np.asarray(measurement @ values[:, column]).ravel()
        solution = splinalg.lsmr(
            gram, rhs, atol=tolerance, btol=tolerance,
            maxiter=max(200, 4 * gram.shape[0]),
        )[0]
        centred[:, column] -= np.asarray(measurement.T @ solution).ravel()
        # One refinement removes accumulated sparse/iterative error.
        residual = np.asarray(measurement @ centred[:, column]).ravel()
        if np.max(np.abs(residual), initial=0.0) > tolerance:
            correction = splinalg.lsmr(
                gram, residual, atol=tolerance * 0.1, btol=tolerance * 0.1,
                maxiter=max(200, 4 * gram.shape[0]),
            )[0]
            centred[:, column] -= np.asarray(measurement.T @ correction).ravel()
        error = np.asarray(measurement @ centred[:, column]).ravel()
        max_error = max(max_error, float(np.max(np.abs(error), initial=0.0)))
    return centred, max_error


def _within_fold_permutation(
    values: np.ndarray, folds: np.ndarray, rng: np.random.Generator,
) -> np.ndarray:
    result = values.copy()
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        if len(rows) > 1:
            result[rows] = values[rng.permutation(rows)]
    return result


def _coefficient_stability(
    coefficients: list[np.ndarray], feature_start: int, n_outputs: int,
) -> np.ndarray:
    stability = np.zeros(n_outputs, np.float64)
    if len(coefficients) < 2:
        return stability
    for output in range(n_outputs):
        vectors = [coef[feature_start:, output] for coef in coefficients]
        similarities = []
        for left in range(len(vectors)):
            for right in range(left + 1, len(vectors)):
                norm = np.linalg.norm(vectors[left]) * np.linalg.norm(vectors[right])
                if norm > EPS:
                    similarities.append(float(np.dot(vectors[left], vectors[right]) / norm))
        stability[output] = float(np.median(similarities)) if similarities else 0.0
    return stability


def fit_slice_internal_weak_state(
    measurement: sparse.spmatrix,
    cell_features: np.ndarray,
    group_target: np.ndarray,
    group_context: np.ndarray,
    group_coords: np.ndarray,
    cell_fold: np.ndarray,
    *,
    n_folds: int = 5,
    alphas: tuple[float, ...] = (1.0, 10.0, 100.0, 1000.0),
    permutations: int = 100,
    buffer_distance: float = 0.0,
    seed: int = 20260827,
) -> WeakStateResult:
    """Fit a cross-fitted multiple-instance nucleus-to-state mapping.

    ``group_target`` is a spot-by-lineage low-dimensional state inferred only
    from the current slice.  The feature model is evaluated incrementally over
    ``group_context``.  Its cell-level output is projected into the exact null
    space of the measurement operator before it can affect allocation.
    """

    measurement = sparse.csr_matrix(measurement, dtype=np.float64)
    cell_features = _as_2d(cell_features)
    target = _as_2d(group_target)
    context = _as_2d(group_context)
    coords = _as_2d(group_coords)
    cell_fold = np.asarray(cell_fold, dtype=np.int32)
    if measurement.shape != (len(target), len(cell_features)):
        raise ValueError("measurement shape does not match groups and cells")
    if len(context) != len(target) or len(coords) != len(target):
        raise ValueError("group arrays are not aligned")
    if len(cell_fold) != len(cell_features):
        raise ValueError("cell_fold is not aligned to cells")

    folds = contiguous_spatial_folds(coords, n_folds=n_folds)
    basis = spline_feature_basis(cell_features)
    aggregate_basis = np.asarray(measurement @ basis)
    design = np.concatenate((context, aggregate_basis), axis=1)
    feature_start = context.shape[1]

    alpha_scores = []
    for alpha in alphas:
        prediction, _ = _crossfit_group_prediction(
            design, target, coords, folds, alpha, buffer_distance,
        )
        alpha_scores.append(float(np.nanmedian(r2_columns(target, prediction))))
    selected_alpha = float(alphas[int(np.argmax(alpha_scores))])
    full_prediction, fold_coefs = _crossfit_group_prediction(
        design, target, coords, folds, selected_alpha, buffer_distance,
    )
    context_prediction, _ = _crossfit_group_prediction(
        context, target, coords, folds, selected_alpha, buffer_distance,
    )
    full_r2 = r2_columns(target, full_prediction)
    context_r2 = r2_columns(target, context_prediction)
    increment = np.full_like(full_r2, -np.inf)
    finite_increment = np.isfinite(full_r2) & np.isfinite(context_r2)
    increment[finite_increment] = full_r2[finite_increment] - context_r2[finite_increment]
    stability = _coefficient_stability(fold_coefs, feature_start, target.shape[1])

    rng = np.random.default_rng(seed)
    null_increment = np.zeros((permutations, target.shape[1]), np.float64)
    for replicate in range(permutations):
        permuted_basis = _within_fold_permutation(aggregate_basis, folds, rng)
        null_design = np.concatenate((context, permuted_basis), axis=1)
        null_prediction, _ = _crossfit_group_prediction(
            null_design, target, coords, folds, selected_alpha, buffer_distance,
        )
        null_r2 = r2_columns(target, null_prediction)
        valid_null = np.isfinite(null_r2) & np.isfinite(context_r2)
        null_increment[replicate] = -np.inf
        null_increment[replicate, valid_null] = (
            null_r2[valid_null] - context_r2[valid_null]
        )
    null_q95 = (
        np.quantile(null_increment, 0.95, axis=0)
        if permutations > 0 else np.zeros(target.shape[1], np.float64)
    )
    excess = increment - np.maximum(null_q95, 0.0)
    eta = np.zeros(target.shape[1], np.float64)
    eligible = (full_r2 > 0.0) & (excess > 0.0) & (stability > 0.0)
    eta[eligible] = 0.25
    eta[eligible & (excess > 0.005)] = 0.50
    eta[eligible & (excess > 0.020)] = 0.75
    eta[eligible & (excess > 0.050)] = 1.00

    raw_cell = np.zeros((len(cell_features), target.shape[1]), np.float64)
    for fold, coef in zip(np.unique(folds), fold_coefs):
        chosen = cell_fold == fold
        if np.any(chosen) and coef.shape[0] >= feature_start + basis.shape[1]:
            raw_cell[chosen] = basis[chosen] @ coef[feature_start:, :]
    raw_cell *= eta[None, :]
    morph_contrast, centering_error = measurement_nullspace_center(
        measurement, raw_cell,
    )
    audit = {
        "schema": "spot2cell_mv4_slice_internal_weak_state_v1",
        "selected_alpha": selected_alpha,
        "alpha_cv_median_r2": alpha_scores,
        "full_heldout_r2": full_r2.tolist(),
        "context_heldout_r2": context_r2.tolist(),
        "incremental_heldout_r2": increment.tolist(),
        "permutation_increment_q95": null_q95.tolist(),
        "coefficient_stability_median_cosine": stability.tolist(),
        "eta": eta.tolist(),
        "n_groups": int(measurement.shape[0]),
        "n_cells": int(measurement.shape[1]),
        "n_basis": int(basis.shape[1]),
        "n_folds": int(n_folds),
        "buffer_distance": float(buffer_distance),
        "nullspace_centering_max_abs_error": centering_error,
        "external_reference_used": False,
    }
    return WeakStateResult(
        morph_contrast=morph_contrast.astype(np.float32),
        oof_prediction=full_prediction.astype(np.float32),
        context_prediction=context_prediction.astype(np.float32),
        folds=folds,
        eta=eta.astype(np.float32),
        audit=audit,
    )


def build_boundary_aware_graph(
    coords: np.ndarray,
    features: np.ndarray,
    boundary_context: np.ndarray,
    component: np.ndarray,
    owner_group: np.ndarray,
    *,
    neighbors: int = 10,
    radius: float = np.inf,
    boundary_strength: float = 3.0,
    technical_continuity_boost: float = 1.25,
) -> GraphResult:
    """Build a same-lineage graph that distinguishes seams from boundaries."""

    coords = _as_2d(coords)
    features = _as_2d(features)
    boundary_context = _as_2d(boundary_context)
    component = np.asarray(component)
    owner_group = np.asarray(owner_group)
    n_cells = len(coords)
    if not (len(features) == len(boundary_context) == len(component) == len(owner_group) == n_cells):
        raise ValueError("graph inputs are not aligned")
    if n_cells <= 1:
        return GraphResult(
            sparse.csr_matrix((n_cells, n_cells)), np.empty((0, 5)),
            np.empty(0), {"n_edges": 0, "n_isolated": n_cells},
        )
    query_k = min(max(2, neighbors + 1), n_cells)
    distance, index = cKDTree(coords[:, :2]).query(
        coords[:, :2], k=query_k, distance_upper_bound=radius, workers=-1,
    )
    if query_k == 1:
        distance = distance[:, None]
        index = index[:, None]
    source = np.repeat(np.arange(n_cells, dtype=np.int32), query_k - 1)
    target = index[:, 1:].reshape(-1)
    spatial_distance = distance[:, 1:].reshape(-1)
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
        return GraphResult(
            sparse.csr_matrix((n_cells, n_cells)), np.empty((0, 5)),
            np.empty(0), {"n_edges": 0, "n_isolated": n_cells},
        )

    feature_distance = np.sqrt(np.mean(np.square(features[source] - features[target]), axis=1))
    context_distance = np.sqrt(
        np.mean(np.square(boundary_context[source] - boundary_context[target]), axis=1)
    )

    def positive_median(values: np.ndarray) -> float:
        positive = values[np.isfinite(values) & (values > 0)]
        return float(np.median(positive)) if len(positive) else 1.0

    spatial_scale = positive_median(spatial_distance)
    feature_scale = positive_median(feature_distance)
    context_scale = positive_median(context_distance)
    boundary_score = 1.0 - np.exp(
        -0.5 * np.square(context_distance / max(context_scale, EPS))
    )
    base_weight = np.exp(
        -0.5 * np.square(spatial_distance / max(spatial_scale, EPS))
        -0.25 * np.square(feature_distance / max(feature_scale, EPS))
    )
    cross_technical = owner_group[source] != owner_group[target]
    seam_boost = np.ones(len(source), np.float64)
    seam_boost[cross_technical & (boundary_score < 0.5)] = technical_continuity_boost
    weight = base_weight * np.exp(-boundary_strength * boundary_score) * seam_boost
    weight = np.clip(weight, 0.0, 1.0)

    directed = sparse.coo_matrix(
        (weight, (source, target)), shape=(n_cells, n_cells),
    ).tocsr()
    directed.sum_duplicates()
    adjacency = directed.minimum(directed.T).tocsr()
    degree = np.asarray(adjacency.sum(1)).ravel()
    isolated = degree <= 0
    if np.any(isolated):
        union = directed.maximum(directed.T).tocsr()
        rescue = sparse.diags(isolated.astype(np.float64)) @ union
        adjacency = (adjacency + rescue).maximum((adjacency + rescue).T).tocsr()
    adjacency.setdiag(0.0)
    adjacency.eliminate_zeros()
    degree = np.asarray(adjacency.sum(1)).ravel()
    inv_sqrt = np.divide(
        1.0, np.sqrt(degree), out=np.zeros_like(degree), where=degree > 0,
    )
    normalized = sparse.diags(inv_sqrt) @ adjacency @ sparse.diags(inv_sqrt)
    laplacian = sparse.diags((degree > 0).astype(np.float64)) - normalized
    upper = sparse.triu(adjacency, k=1).tocoo()
    edge_source = upper.row.astype(np.int32)
    edge_target = upper.col.astype(np.int32)
    edge_context_distance = np.sqrt(np.mean(np.square(
        boundary_context[edge_source] - boundary_context[edge_target]
    ), axis=1))
    edge_boundary = 1.0 - np.exp(
        -0.5 * np.square(edge_context_distance / max(context_scale, EPS))
    )
    edge_cross = owner_group[edge_source] != owner_group[edge_target]
    edges = np.column_stack((
        edge_source, edge_target, upper.data, edge_boundary, edge_cross.astype(np.float64),
    )).astype(np.float64)
    audit = {
        "n_edges": int(len(upper.data)),
        "n_isolated": int(np.sum(degree <= 0)),
        "degree_median": float(np.median(np.diff(adjacency.indptr))),
        "spatial_scale": spatial_scale,
        "feature_scale": feature_scale,
        "boundary_context_scale": context_scale,
        "boundary_strength": float(boundary_strength),
        "cross_spot_low_boundary_edges": int(np.sum(edge_cross & (edge_boundary < 0.5))),
        "cross_spot_high_boundary_edges": int(np.sum(edge_cross & (edge_boundary >= 0.5))),
    }
    return GraphResult(laplacian.tocsr(), edges, edge_boundary, audit)


def solve_graph_inverse(
    measurement: sparse.spmatrix,
    group_target: np.ndarray,
    laplacian: sparse.spmatrix,
    weak_prior: np.ndarray,
    *,
    lambda_prior: float = 0.25,
    lambda_graph: float = 0.5,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Solve the graph-regularised cell-state inverse and return its contrast."""

    measurement = sparse.csr_matrix(measurement, dtype=np.float64)
    target = _as_2d(group_target)
    prior = _as_2d(weak_prior)
    laplacian = sparse.csr_matrix(laplacian, dtype=np.float64)
    if measurement.shape != (len(target), len(prior)):
        raise ValueError("inverse inputs are not aligned")
    row_square = np.asarray(measurement.multiply(measurement).sum(1)).ravel()
    effective_n = np.divide(1.0, row_square, out=np.ones_like(row_square), where=row_square > EPS)
    effective_n = np.clip(effective_n, 1.0, 100.0)
    precision = measurement.T @ sparse.diags(effective_n) @ measurement
    positive = precision.diagonal()
    positive = positive[positive > 0]
    scale = float(np.median(positive)) if len(positive) else 1.0
    prior_penalty = float(lambda_prior) * scale
    graph_penalty = float(lambda_graph) * scale
    system = (
        precision + prior_penalty * sparse.eye(measurement.shape[1], format="csr")
        + graph_penalty * laplacian
    ).tocsr()
    rhs = measurement.T @ (effective_n[:, None] * target) + prior_penalty * prior
    state = np.zeros_like(prior, dtype=np.float64)
    status = []
    for column in range(target.shape[1]):
        state[:, column], info = splinalg.cg(
            system, rhs[:, column], x0=prior[:, column],
            rtol=1.0e-7, atol=0.0, maxiter=800,
        )
        status.append(int(info))
    graph_contrast, centering_error = measurement_nullspace_center(measurement, state)
    predicted = np.asarray(measurement @ state)
    audit = {
        "group_reconstruction_r2": r2_columns(target, predicted).tolist(),
        "lambda_prior": prior_penalty,
        "lambda_graph": graph_penalty,
        "cg_status": status,
        "contrast_centering_max_abs_error": centering_error,
    }
    return state.astype(np.float32), graph_contrast.astype(np.float32), audit


def orthogonalize_nullspace_contrasts(
    measurement: sparse.spmatrix,
    primary: np.ndarray,
    secondary: np.ndarray,
) -> tuple[np.ndarray, dict]:
    """Remove the primary spatial direction from a secondary morphology term."""

    primary = _as_2d(primary)
    secondary = _as_2d(secondary)
    if primary.shape != secondary.shape:
        raise ValueError("contrast matrices must have identical shapes")
    result = secondary.copy()
    cosine_before = []
    for column in range(primary.shape[1]):
        denominator = float(np.dot(primary[:, column], primary[:, column]))
        cross = float(np.dot(primary[:, column], result[:, column]))
        norm = np.linalg.norm(primary[:, column]) * np.linalg.norm(result[:, column])
        cosine_before.append(cross / norm if norm > EPS else 0.0)
        if denominator > EPS:
            result[:, column] -= primary[:, column] * (cross / denominator)
    result, centering_error = measurement_nullspace_center(measurement, result)
    cosine_after = []
    for column in range(primary.shape[1]):
        norm = np.linalg.norm(primary[:, column]) * np.linalg.norm(result[:, column])
        cosine_after.append(
            float(np.dot(primary[:, column], result[:, column]) / norm) if norm > EPS else 0.0
        )
    return result.astype(np.float32), {
        "cosine_before": cosine_before,
        "cosine_after": cosine_after,
        "centering_max_abs_error": centering_error,
    }


def exact_mass_projection(
    observed: np.ndarray,
    log_propensity: np.ndarray,
    fallback_weight: np.ndarray,
    *,
    library_ipf_strength: float = 0.0,
    library_ipf_iterations: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Project cell logits to non-negative, exactly conserved gene mass."""

    observed = np.asarray(observed, dtype=np.float64)
    logits = _as_2d(log_propensity)
    fallback_weight = np.asarray(fallback_weight, dtype=np.float64)
    if logits.shape[1] != len(observed) or logits.shape[0] != len(fallback_weight):
        raise ValueError("projection inputs are not aligned")
    if np.any(observed < 0) or np.any(fallback_weight < 0):
        raise ValueError("observed mass and fallback weights must be non-negative")
    if not 0.0 <= library_ipf_strength <= 1.0:
        raise ValueError("library_ipf_strength must be in [0,1]")
    shifted = logits - np.max(logits, axis=0, keepdims=True)
    propensity = np.exp(np.clip(shifted, -80.0, 0.0))
    denominator = np.sum(propensity, axis=0)
    missing = ~np.isfinite(denominator) | (denominator <= EPS)
    if np.any(missing):
        fallback = fallback_weight / max(float(np.sum(fallback_weight)), EPS)
        propensity[:, missing] = fallback[:, None]
        denominator[missing] = 1.0
    responsibility = propensity / np.maximum(denominator[None, :], EPS)
    if library_ipf_strength > 0 and library_ipf_iterations > 0:
        target_library = (
            float(np.sum(observed)) * fallback_weight
            / max(float(np.sum(fallback_weight)), EPS)
        )
        for _ in range(library_ipf_iterations):
            assigned = responsibility @ observed
            scale = np.power(
                np.divide(
                    target_library, np.maximum(assigned, EPS),
                    out=np.ones_like(target_library), where=assigned > EPS,
                ),
                library_ipf_strength,
            )
            responsibility *= np.clip(scale, 0.25, 4.0)[:, None]
            responsibility /= np.maximum(np.sum(responsibility, axis=0, keepdims=True), EPS)
    responsibility = responsibility.astype(np.float32)
    pivot = np.argmax(responsibility, axis=0)
    correction = 1.0 - np.sum(responsibility, axis=0, dtype=np.float64)
    responsibility[pivot, np.arange(len(observed))] = (
        responsibility[pivot, np.arange(len(observed))].astype(np.float64) + correction
    ).astype(np.float32)
    mass = responsibility.astype(np.float64) * observed[None, :]
    return responsibility, mass
