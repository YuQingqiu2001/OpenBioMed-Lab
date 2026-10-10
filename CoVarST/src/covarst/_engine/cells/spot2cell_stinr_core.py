#!/usr/bin/env python3
"""Low-rank implicit spatial rate models for Spot2Cell-MV4.1-STINR.

This module adapts the coordinate-to-expression idea and sine activations from
the public STINR CVPR 2025 implementation to the P1 single-slice inverse
problem.  It deliberately does not contain DLPFC labels, an scRNA basis, or
file-system paths.  Observed counts are never overwritten: the model estimates
latent gene probabilities, while exact count allocation remains a separate
projection.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import io
import random
from typing import Iterable

import numpy as np
from scipy import sparse
import torch
from torch import nn
from torch.nn import functional as F


EPS = 1.0e-8


@dataclass
class FitResult:
    """A fitted INR plus compact training evidence."""

    model: "TypeConditionedINR"
    audit: dict


def set_deterministic_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True


def normalize_coordinates(
    coords: np.ndarray,
    center: np.ndarray | None = None,
    scale: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Robustly map two-dimensional coordinates to an order-one range."""

    coords = np.asarray(coords, dtype=np.float64)
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError("coords must have shape (n, >=2)")
    xy = coords[:, :2]
    if center is None:
        center = np.nanmedian(xy, axis=0)
    if scale is None:
        q05, q95 = np.nanquantile(xy, [0.05, 0.95], axis=0)
        scale = 0.5 * (q95 - q05)
    center = np.asarray(center, dtype=np.float64)
    scale = np.asarray(scale, dtype=np.float64)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-6), scale, 1.0)
    result = np.clip((xy - center) / scale, -3.0, 3.0)
    return result.astype(np.float32), center, scale


def split_sparse_counts(
    counts: sparse.spmatrix,
    fraction: float = 0.5,
    seed: int = 20260827,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix]:
    """Binomially thin every observed count into independent train/test parts."""

    if not 0.0 < fraction < 1.0:
        raise ValueError("fraction must lie strictly between zero and one")
    matrix = sparse.csr_matrix(counts)
    if np.any(matrix.data < 0):
        raise ValueError("counts must be non-negative")
    rounded = np.rint(matrix.data).astype(np.int64)
    if not np.allclose(matrix.data, rounded, atol=1.0e-6):
        raise ValueError("binomial thinning requires integer counts")
    rng = np.random.default_rng(seed)
    left_data = rng.binomial(rounded, fraction).astype(np.int32)
    right_data = (rounded - left_data).astype(np.int32)
    left = sparse.csr_matrix(
        (left_data, matrix.indices.copy(), matrix.indptr.copy()),
        shape=matrix.shape,
    )
    right = sparse.csr_matrix(
        (right_data, matrix.indices.copy(), matrix.indptr.copy()),
        shape=matrix.shape,
    )
    left.eliminate_zeros()
    right.eliminate_zeros()
    return left, right


def estimate_type_profiles(
    counts: sparse.spmatrix | np.ndarray,
    class_id: np.ndarray,
    n_classes: int,
    pseudocount: float = 0.5,
) -> np.ndarray:
    """Estimate modality-internal type offsets with weak parent-free shrinkage."""

    class_id = np.asarray(class_id, dtype=np.int64)
    if sparse.issparse(counts):
        matrix = sparse.csr_matrix(counts)
    else:
        matrix = np.asarray(counts, dtype=np.float64)
    if matrix.shape[0] != len(class_id):
        raise ValueError("counts and class_id are not aligned")
    global_sum = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    global_probability = (global_sum + pseudocount) / (
        np.sum(global_sum) + pseudocount * matrix.shape[1]
    )
    result = np.zeros((n_classes, matrix.shape[1]), dtype=np.float64)
    for value in range(n_classes):
        rows = class_id == value
        if not np.any(rows):
            result[value] = global_probability
            continue
        class_sum = np.asarray(matrix[rows].sum(axis=0)).ravel().astype(np.float64)
        # The global prior prevents rare classes from creating hard zeros.
        prior_mass = max(20.0, np.sqrt(max(float(np.sum(class_sum)), 1.0)))
        smoothed = class_sum + prior_mass * global_probability
        result[value] = smoothed / max(float(np.sum(smoothed)), EPS)
    return result.astype(np.float32)


class SineLayer(nn.Module):
    """Compact sine layer following the public STINR implementation."""

    def __init__(self, c_in: int, c_out: int, omega: float = 1.0):
        super().__init__()
        self.omega = float(omega)
        self.linear = nn.Linear(c_in, c_out)
        bound = np.sqrt(6.0 / float(c_in + c_out))
        nn.init.uniform_(self.linear.weight, -bound, bound)
        nn.init.zeros_(self.linear.bias)

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        return torch.sin(self.omega * self.linear(values))


class TypeConditionedINR(nn.Module):
    """Coordinate INR with independently learned type-specific gene loadings."""

    def __init__(
        self,
        baseline_probability: np.ndarray,
        rank: int = 8,
        hidden: int = 96,
        omega: float = 2.0,
        coord_center: np.ndarray | None = None,
        coord_scale: np.ndarray | None = None,
    ):
        super().__init__()
        baseline = np.asarray(baseline_probability, dtype=np.float32)
        if baseline.ndim != 2 or np.any(baseline < 0):
            raise ValueError("baseline_probability must be non-negative type-by-gene")
        baseline /= np.maximum(baseline.sum(axis=1, keepdims=True), EPS)
        self.n_classes, self.n_genes = baseline.shape
        self.rank = int(rank)
        self.hidden = int(hidden)
        self.omega = float(omega)
        self.register_buffer("baseline_log", torch.log(torch.from_numpy(baseline) + EPS))
        center = np.zeros(2, np.float32) if coord_center is None else np.asarray(coord_center, np.float32)
        scale = np.ones(2, np.float32) if coord_scale is None else np.asarray(coord_scale, np.float32)
        self.register_buffer("coord_center", torch.from_numpy(center))
        self.register_buffer("coord_scale", torch.from_numpy(scale))
        if self.rank > 0:
            self.encoder = nn.Sequential(
                SineLayer(2, hidden, omega=omega),
                SineLayer(hidden, hidden, omega=omega),
                nn.Linear(hidden, self.rank, bias=False),
            )
            nn.init.uniform_(
                self.encoder[-1].weight,
                -np.sqrt(6.0 / float(hidden + self.rank)),
                np.sqrt(6.0 / float(hidden + self.rank)),
            )
            self.loadings = nn.Parameter(
                torch.zeros(self.n_classes, self.rank, self.n_genes)
            )
            nn.init.normal_(self.loadings, mean=0.0, std=0.01)
        else:
            self.encoder = None
            self.register_parameter("loadings", None)

    def spatial_state(self, coords_normalized: torch.Tensor) -> torch.Tensor:
        if self.rank == 0:
            return torch.zeros(
                (len(coords_normalized), 0),
                dtype=coords_normalized.dtype,
                device=coords_normalized.device,
            )
        return self.encoder(coords_normalized)

    def logits(
        self,
        coords_normalized: torch.Tensor,
        class_id: torch.Tensor,
    ) -> torch.Tensor:
        baseline = self.baseline_log[class_id]
        if self.rank == 0:
            return baseline
        state = self.spatial_state(coords_normalized)
        chosen = self.loadings[class_id]
        residual = torch.einsum("bk,bkg->bg", state, chosen)
        return baseline + torch.clamp(residual, -8.0, 8.0)

    def probabilities(
        self,
        coords_normalized: torch.Tensor,
        class_id: torch.Tensor,
    ) -> torch.Tensor:
        return F.softmax(self.logits(coords_normalized, class_id), dim=1)


def model_sha256(model: nn.Module) -> str:
    buffer = io.BytesIO()
    torch.save(model.state_dict(), buffer)
    return hashlib.sha256(buffer.getvalue()).hexdigest()


def sparse_multinomial_nll(
    counts: sparse.spmatrix,
    probabilities: np.ndarray,
    rows: np.ndarray | None = None,
) -> float:
    """Mean negative log probability per held-out molecule."""

    matrix = sparse.csr_matrix(counts)
    probabilities = np.asarray(probabilities, dtype=np.float64)
    if probabilities.shape != matrix.shape:
        raise ValueError("probability and count shapes differ")
    if rows is None:
        rows = np.arange(matrix.shape[0], dtype=np.int64)
    rows = np.asarray(rows, dtype=np.int64)
    block = matrix[rows].tocoo()
    total = float(np.sum(block.data))
    if total <= 0:
        return float("nan")
    chosen = probabilities[rows[block.row], block.col]
    return float(-np.sum(block.data * np.log(np.maximum(chosen, EPS))) / total)


def dense_multinomial_nll(
    counts: np.ndarray,
    probabilities: np.ndarray,
    rows: np.ndarray | None = None,
) -> float:
    matrix = np.asarray(counts, dtype=np.float64)
    probabilities = np.asarray(probabilities, dtype=np.float64)
    if matrix.shape != probabilities.shape:
        raise ValueError("probability and count shapes differ")
    if rows is not None:
        matrix = matrix[np.asarray(rows, dtype=np.int64)]
        probabilities = probabilities[np.asarray(rows, dtype=np.int64)]
    total = float(np.sum(matrix))
    if total <= 0:
        return float("nan")
    return float(-np.sum(matrix * np.log(np.maximum(probabilities, EPS))) / total)


def _device(value: str | torch.device | None) -> torch.device:
    if value is None:
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(value)


def fit_direct_inr(
    counts: sparse.spmatrix,
    coords: np.ndarray,
    class_id: np.ndarray,
    baseline_probability: np.ndarray,
    *,
    rank: int = 8,
    hidden: int = 96,
    omega: float = 2.0,
    epochs: int = 12,
    batch_size: int = 4096,
    learning_rate: float = 2.0e-3,
    weight_decay: float = 1.0e-5,
    seed: int = 20260827,
    device: str | torch.device | None = None,
    eligible_rows: np.ndarray | None = None,
) -> FitResult:
    """Fit an INR to directly observed cell/location counts."""

    set_deterministic_seed(seed)
    matrix = sparse.csr_matrix(counts, dtype=np.float32)
    class_id = np.asarray(class_id, dtype=np.int64)
    coords_norm, center, scale = normalize_coordinates(coords)
    if matrix.shape[0] != len(coords_norm) or len(class_id) != len(coords_norm):
        raise ValueError("direct observations are not aligned")
    library = np.asarray(matrix.sum(axis=1)).ravel()
    if eligible_rows is None:
        eligible = np.flatnonzero(library > 0)
    else:
        eligible = np.asarray(eligible_rows, dtype=np.int64)
        eligible = eligible[library[eligible] > 0]
    if not len(eligible):
        raise ValueError("no positive-library observations are eligible")
    target_device = _device(device)
    model = TypeConditionedINR(
        baseline_probability,
        rank=rank,
        hidden=hidden,
        omega=omega,
        coord_center=center,
        coord_scale=scale,
    ).to(target_device)
    optimizer = torch.optim.Adamax(
        model.parameters(), lr=learning_rate, weight_decay=weight_decay,
    )
    rng = np.random.default_rng(seed + 1)
    loss_history: list[float] = []
    model.train()
    for _ in range(int(epochs)):
        shuffled = rng.permutation(eligible)
        epoch_loss = 0.0
        epoch_count = 0.0
        for start in range(0, len(shuffled), batch_size):
            rows = shuffled[start:start + batch_size]
            y_np = matrix[rows].toarray().astype(np.float32, copy=False)
            molecules = float(np.sum(y_np))
            if molecules <= 0:
                continue
            xy = torch.from_numpy(coords_norm[rows]).to(target_device)
            cls = torch.from_numpy(class_id[rows]).to(target_device)
            y = torch.from_numpy(y_np).to(target_device)
            log_probability = F.log_softmax(model.logits(xy, cls), dim=1)
            loss = -torch.sum(y * log_probability) / max(molecules, 1.0)
            if model.rank > 0:
                loss = loss + 1.0e-5 * torch.mean(torch.square(model.loadings))
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
            epoch_loss += float(loss.detach().cpu()) * molecules
            epoch_count += molecules
        loss_history.append(epoch_loss / max(epoch_count, 1.0))
    audit = {
        "schema": "mv41_stinr_direct_fit_v1",
        "rank": int(rank),
        "hidden": int(hidden),
        "omega": float(omega),
        "epochs": int(epochs),
        "batch_size": int(batch_size),
        "learning_rate": float(learning_rate),
        "weight_decay": float(weight_decay),
        "n_observations": int(matrix.shape[0]),
        "n_eligible": int(len(eligible)),
        "n_genes": int(matrix.shape[1]),
        "loss_history": loss_history,
        "model_sha256": model_sha256(model),
        "seed": int(seed),
        "device": str(target_device),
        "external_reference_used": False,
    }
    return FitResult(model=model, audit=audit)


def predict_cell_probabilities(
    model: TypeConditionedINR,
    coords: np.ndarray,
    class_id: np.ndarray,
    *,
    batch_size: int = 4096,
    device: str | torch.device | None = None,
) -> np.ndarray:
    coords_norm, _, _ = normalize_coordinates(
        coords,
        center=model.coord_center.detach().cpu().numpy(),
        scale=model.coord_scale.detach().cpu().numpy(),
    )
    class_id = np.asarray(class_id, dtype=np.int64)
    target_device = _device(device)
    model = model.to(target_device)
    model.eval()
    result = np.empty((len(coords_norm), model.n_genes), np.float32)
    with torch.inference_mode():
        for start in range(0, len(coords_norm), batch_size):
            end = min(start + batch_size, len(coords_norm))
            xy = torch.from_numpy(coords_norm[start:end]).to(target_device)
            cls = torch.from_numpy(class_id[start:end]).to(target_device)
            result[start:end] = model.probabilities(xy, cls).cpu().numpy()
    return result


def _spot_batch(
    measurement: sparse.csr_matrix,
    spots: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    block = measurement[np.asarray(spots, dtype=np.int64)].tocoo()
    unique_cells, inverse = np.unique(block.col, return_inverse=True)
    return (
        unique_cells.astype(np.int64),
        inverse.astype(np.int64),
        np.column_stack((block.row.astype(np.int64), block.data.astype(np.float32))),
    )


def _aggregate_spot_probability(
    cell_probability: torch.Tensor,
    local_cell: torch.Tensor,
    spot_row: torch.Tensor,
    weight: torch.Tensor,
    n_spots: int,
) -> torch.Tensor:
    weighted = cell_probability[local_cell] * weight[:, None]
    result = torch.zeros(
        (n_spots, cell_probability.shape[1]),
        dtype=cell_probability.dtype,
        device=cell_probability.device,
    )
    result.index_add_(0, spot_row, weighted)
    return result / torch.clamp(result.sum(dim=1, keepdim=True), min=EPS)


def fit_spot_inr(
    counts: np.ndarray,
    measurement: sparse.spmatrix,
    cell_coords: np.ndarray,
    cell_class_id: np.ndarray,
    baseline_probability: np.ndarray,
    *,
    rank: int = 8,
    hidden: int = 96,
    omega: float = 2.0,
    epochs: int = 30,
    spot_batch_size: int = 192,
    learning_rate: float = 2.0e-3,
    weight_decay: float = 1.0e-5,
    seed: int = 20260827,
    device: str | torch.device | None = None,
    eligible_spots: np.ndarray | None = None,
) -> FitResult:
    """Fit an INR through a sparse spot-to-cell measurement operator."""

    set_deterministic_seed(seed)
    y = np.asarray(counts, dtype=np.float32)
    measurement = sparse.csr_matrix(measurement, dtype=np.float32)
    if measurement.shape[0] != len(y) or measurement.shape[1] != len(cell_coords):
        raise ValueError("spot observations and measurement are not aligned")
    cell_class_id = np.asarray(cell_class_id, dtype=np.int64)
    coords_norm, center, scale = normalize_coordinates(cell_coords)
    library = y.sum(axis=1)
    if eligible_spots is None:
        eligible = np.flatnonzero((library > 0) & (np.diff(measurement.indptr) > 0))
    else:
        eligible = np.asarray(eligible_spots, dtype=np.int64)
        eligible = eligible[(library[eligible] > 0) & (np.diff(measurement.indptr)[eligible] > 0)]
    if not len(eligible):
        raise ValueError("no positive-library spots are eligible")
    target_device = _device(device)
    model = TypeConditionedINR(
        baseline_probability,
        rank=rank,
        hidden=hidden,
        omega=omega,
        coord_center=center,
        coord_scale=scale,
    ).to(target_device)
    optimizer = torch.optim.Adamax(
        model.parameters(), lr=learning_rate, weight_decay=weight_decay,
    )
    rng = np.random.default_rng(seed + 1)
    loss_history: list[float] = []
    model.train()
    for _ in range(int(epochs)):
        shuffled = rng.permutation(eligible)
        epoch_loss = 0.0
        epoch_count = 0.0
        for start in range(0, len(shuffled), spot_batch_size):
            spots = shuffled[start:start + spot_batch_size]
            unique_cells, inverse, links = _spot_batch(measurement, spots)
            if not len(unique_cells):
                continue
            xy = torch.from_numpy(coords_norm[unique_cells]).to(target_device)
            cls = torch.from_numpy(cell_class_id[unique_cells]).to(target_device)
            cell_probability = model.probabilities(xy, cls)
            local_cell = torch.from_numpy(inverse).to(target_device)
            spot_row = torch.from_numpy(links[:, 0].astype(np.int64)).to(target_device)
            weight = torch.from_numpy(links[:, 1].astype(np.float32)).to(target_device)
            predicted = _aggregate_spot_probability(
                cell_probability, local_cell, spot_row, weight, len(spots),
            )
            y_np = y[spots]
            molecules = float(np.sum(y_np))
            target = torch.from_numpy(y_np).to(target_device)
            loss = -torch.sum(target * torch.log(predicted + EPS)) / max(molecules, 1.0)
            if model.rank > 0:
                loss = loss + 1.0e-5 * torch.mean(torch.square(model.loadings))
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
            epoch_loss += float(loss.detach().cpu()) * molecules
            epoch_count += molecules
        loss_history.append(epoch_loss / max(epoch_count, 1.0))
    audit = {
        "schema": "mv41_stinr_spot_fit_v1",
        "rank": int(rank),
        "hidden": int(hidden),
        "omega": float(omega),
        "epochs": int(epochs),
        "spot_batch_size": int(spot_batch_size),
        "learning_rate": float(learning_rate),
        "weight_decay": float(weight_decay),
        "n_spots": int(y.shape[0]),
        "n_eligible": int(len(eligible)),
        "n_cells": int(measurement.shape[1]),
        "n_genes": int(y.shape[1]),
        "loss_history": loss_history,
        "model_sha256": model_sha256(model),
        "seed": int(seed),
        "device": str(target_device),
        "external_reference_used": False,
    }
    return FitResult(model=model, audit=audit)


def aggregate_cell_probabilities(
    measurement: sparse.spmatrix,
    cell_probability: np.ndarray,
    *,
    batch_size: int = 256,
) -> np.ndarray:
    measurement = sparse.csr_matrix(measurement, dtype=np.float64)
    values = np.asarray(cell_probability, dtype=np.float64)
    if measurement.shape[1] != len(values):
        raise ValueError("measurement and cell_probability are not aligned")
    result = np.empty((measurement.shape[0], values.shape[1]), np.float32)
    for start in range(0, measurement.shape[0], batch_size):
        end = min(start + batch_size, measurement.shape[0])
        block = np.asarray(measurement[start:end] @ values)
        block /= np.maximum(block.sum(axis=1, keepdims=True), EPS)
        result[start:end] = block.astype(np.float32)
    return result


def posterior_shrinkage(
    counts: sparse.spmatrix,
    prior_probability: np.ndarray,
    tau: float,
) -> np.ndarray:
    """Empirical-Bayes cell probability retaining each cell's observed counts."""

    matrix = sparse.csr_matrix(counts, dtype=np.float64)
    prior = np.asarray(prior_probability, dtype=np.float64)
    if matrix.shape != prior.shape:
        raise ValueError("counts and prior_probability are not aligned")
    if tau < 0:
        raise ValueError("tau must be non-negative")
    result = prior * float(tau)
    coo = matrix.tocoo()
    np.add.at(result, (coo.row, coo.col), coo.data)
    denominator = np.asarray(matrix.sum(axis=1)).ravel() + float(tau)
    zero = denominator <= 0
    denominator[zero] = 1.0
    result /= denominator[:, None]
    result[zero] = prior[zero]
    result /= np.maximum(result.sum(axis=1, keepdims=True), EPS)
    return result.astype(np.float32)


def select_posterior_tau(
    train_counts: sparse.spmatrix,
    test_counts: sparse.spmatrix,
    prior_probability: np.ndarray,
    candidates: Iterable[float] = (0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0),
) -> tuple[float, list[dict]]:
    audits: list[dict] = []
    best_tau = None
    best_loss = np.inf
    for value in candidates:
        tau = float(value)
        posterior = posterior_shrinkage(train_counts, prior_probability, tau)
        loss = sparse_multinomial_nll(test_counts, posterior)
        audits.append({"tau": tau, "heldout_nll_per_molecule": loss})
        if np.isfinite(loss) and loss < best_loss:
            best_loss = loss
            best_tau = tau
    if best_tau is None:
        best_tau = 0.0
    return float(best_tau), audits


def exact_panel_mass(
    counts_by_spot: np.ndarray,
    geometry: sparse.spmatrix,
    capacity: np.ndarray,
    cell_probability: np.ndarray,
) -> tuple[np.ndarray, dict]:
    """Allocate observed panel molecules by latent rates with exact conservation."""

    counts = np.asarray(counts_by_spot, dtype=np.float64)
    geometry = sparse.csr_matrix(geometry, dtype=np.float64)
    capacity = np.asarray(capacity, dtype=np.float64)
    probability = np.asarray(cell_probability, dtype=np.float64)
    if geometry.shape[0] != len(counts) or geometry.shape[1] != len(probability):
        raise ValueError("mass projection inputs are not aligned")
    if len(capacity) != len(probability):
        raise ValueError("capacity and cell_probability are not aligned")
    output = np.zeros_like(probability, dtype=np.float32)
    occupied = np.diff(geometry.indptr) > 0
    total_observed_by_gene = counts.sum(axis=0)
    unassigned_by_gene = counts[~occupied].sum(axis=0)
    observed = np.zeros(probability.shape[1], np.float64)
    assigned = np.zeros(probability.shape[1], np.float64)
    max_error = 0.0
    for spot in range(geometry.shape[0]):
        start, end = int(geometry.indptr[spot]), int(geometry.indptr[spot + 1])
        if start == end:
            continue
        cells = geometry.indices[start:end]
        link = geometry.data[start:end] * capacity[cells]
        prop = link[:, None] * probability[cells]
        denominator = prop.sum(axis=0)
        fallback = denominator <= EPS
        if np.any(fallback):
            prop[:, fallback] = link[:, None]
            denominator[fallback] = np.sum(link)
        responsibility = prop / np.maximum(denominator[None, :], EPS)
        block = responsibility * counts[spot][None, :]
        output[cells] += block.astype(np.float32)
        reconstructed = block.sum(axis=0)
        observed += counts[spot]
        assigned += reconstructed
        max_error = max(max_error, float(np.max(np.abs(reconstructed - counts[spot]))))
    audit = {
        "observed_panel_mass": float(np.sum(observed)),
        "observed_panel_mass_semantics": "positive-candidate-cell spots only",
        "total_observed_panel_mass": float(np.sum(total_observed_by_gene)),
        "unassigned_zero_cell_panel_mass": float(np.sum(unassigned_by_gene)),
        "unassigned_zero_cell_by_gene": unassigned_by_gene.astype(float).tolist(),
        "n_zero_cell_spots": int(np.count_nonzero(~occupied)),
        "assigned_panel_mass": float(np.sum(assigned)),
        "max_abs_gene_conservation_error": float(np.max(np.abs(observed - assigned))),
        "max_abs_spot_gene_conservation_error": max_error,
        "max_abs_total_balance_gene_error": float(np.max(np.abs(
            total_observed_by_gene - assigned - unassigned_by_gene
        ))),
    }
    return output, audit
