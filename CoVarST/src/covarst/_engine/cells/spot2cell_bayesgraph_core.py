#!/usr/bin/env python3
"""BayesTME-inspired sparse graph programs for slice-internal imputation.

This module implements only the mathematical ideas needed by the P1 MV4.2
prototype: low-rank type-specific expression programs, stochastic second-order
graph trend filtering, sparse gene loadings and count-likelihood fitting.  It is
not the BayesTME package and does not use external cell-type references.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import random

import numpy as np
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
import torch
from torch import nn


EPS = 1.0e-10


@dataclass
class GraphData:
    edges: np.ndarray
    edge_weight: np.ndarray
    triplets: np.ndarray
    triplet_weight: np.ndarray
    audit: dict


@dataclass
class ProgramFit:
    model: nn.Module | None
    w: np.ndarray
    v: np.ndarray
    probabilities: np.ndarray
    audit: dict


def set_deterministic(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(False)


def class_profiles(counts: sparse.spmatrix, class_id: np.ndarray, n_classes: int,
                   pseudocount: float = 0.25) -> np.ndarray:
    counts = sparse.csr_matrix(counts, dtype=np.float64)
    class_id = np.asarray(class_id, dtype=np.int64)
    result = np.full((n_classes, counts.shape[1]), pseudocount, np.float64)
    for value in range(n_classes):
        rows = np.flatnonzero(class_id == value)
        if len(rows):
            result[value] += np.asarray(counts[rows].sum(axis=0)).ravel()
    result /= np.maximum(result.sum(axis=1, keepdims=True), EPS)
    return result.astype(np.float32)


def _combined_groups(group: np.ndarray | None, component: np.ndarray | None,
                     n: int) -> np.ndarray:
    if group is None and component is None:
        return np.zeros(n, np.int64)
    if group is None:
        group = np.zeros(n, np.int64)
    if component is None:
        component = np.zeros(n, np.int64)
    group = np.asarray(group, dtype=np.int64)
    component = np.asarray(component, dtype=np.int64)
    pairs = np.column_stack((group, component))
    _, inverse = np.unique(pairs, axis=0, return_inverse=True)
    return inverse.astype(np.int64)


def build_spatial_graph(
    coords: np.ndarray,
    *,
    group: np.ndarray | None = None,
    component: np.ndarray | None = None,
    context: np.ndarray | None = None,
    neighbors: int = 6,
    radius_factor: float = 4.0,
) -> GraphData:
    """Build an undirected kNN graph and one second-difference triplet per node."""

    coords = np.asarray(coords, dtype=np.float64)
    n = len(coords)
    combined = _combined_groups(group, component, n)
    directed_src: list[np.ndarray] = []
    directed_dst: list[np.ndarray] = []
    directed_dist: list[np.ndarray] = []
    neighbor_table = np.full((n, neighbors), -1, np.int64)
    distance_table = np.full((n, neighbors), np.inf, np.float64)
    for label in np.unique(combined):
        rows = np.flatnonzero(combined == label)
        if len(rows) < 2:
            continue
        k = min(neighbors + 1, len(rows))
        model = NearestNeighbors(n_neighbors=k, metric="euclidean", n_jobs=-1)
        model.fit(coords[rows])
        distance, index = model.kneighbors(coords[rows], return_distance=True)
        distance = distance[:, 1:]
        index = rows[index[:, 1:]]
        first = distance[:, 0]
        local_scale = np.maximum(first, np.median(first[first > 0]) if np.any(first > 0) else 1.0)
        keep = distance <= radius_factor * local_scale[:, None]
        width = distance.shape[1]
        neighbor_table[np.ix_(rows, np.arange(width))] = np.where(keep, index, -1)
        distance_table[np.ix_(rows, np.arange(width))] = np.where(keep, distance, np.inf)
        rr, cc = np.where(keep)
        directed_src.append(rows[rr])
        directed_dst.append(index[rr, cc])
        directed_dist.append(distance[rr, cc])
    if not directed_src:
        return GraphData(
            edges=np.empty((0, 2), np.int64), edge_weight=np.empty(0, np.float32),
            triplets=np.empty((0, 3), np.int64), triplet_weight=np.empty(0, np.float32),
            audit={"n_nodes": n, "n_edges": 0, "n_triplets": 0},
        )
    src = np.concatenate(directed_src)
    dst = np.concatenate(directed_dst)
    dist = np.concatenate(directed_dist)
    pair = np.sort(np.column_stack((src, dst)), axis=1)
    keep = pair[:, 0] != pair[:, 1]
    pair = pair[keep]
    dist = dist[keep]
    order = np.lexsort((pair[:, 1], pair[:, 0]))
    pair = pair[order]
    dist = dist[order]
    unique = np.ones(len(pair), bool)
    unique[1:] = np.any(pair[1:] != pair[:-1], axis=1)
    edges = pair[unique]
    edge_dist = dist[unique]
    scale = np.median(edge_dist[edge_dist > 0]) if np.any(edge_dist > 0) else 1.0
    weight = np.exp(-0.5 * np.square(edge_dist / max(scale, EPS)))
    if context is not None and len(edges):
        context = np.asarray(context, dtype=np.float64)
        delta = np.linalg.norm(context[edges[:, 0]] - context[edges[:, 1]], axis=1)
        cscale = np.median(delta[delta > 0]) if np.any(delta > 0) else 1.0
        weight *= np.exp(-0.5 * np.square(delta / max(cscale, EPS)))
    weight = np.clip(weight, 1.0e-3, 1.0).astype(np.float32)

    # For each node choose the most nearly opposite neighbor to its nearest
    # neighbor.  This is a stochastic-friendly graph analogue of BayesTME's
    # three-node second-order trend-filtering penalty.
    triplets = []
    triplet_weight = []
    edge_lookup = {(int(a), int(b)): float(w) for (a, b), w in zip(edges, weight)}
    for center in range(n):
        nbr = neighbor_table[center]
        nbr = nbr[nbr >= 0]
        if len(nbr) < 2:
            continue
        vectors = coords[nbr] - coords[center]
        norm = np.linalg.norm(vectors, axis=1)
        valid = norm > EPS
        nbr = nbr[valid]
        vectors = vectors[valid]
        norm = norm[valid]
        if len(nbr) < 2:
            continue
        unit = vectors / norm[:, None]
        second = int(np.argmin(unit @ unit[0]))
        if second == 0:
            second = 1
        a, c = int(nbr[0]), int(nbr[second])
        wa = edge_lookup.get(tuple(sorted((a, center))), 1.0e-3)
        wc = edge_lookup.get(tuple(sorted((c, center))), 1.0e-3)
        triplets.append((a, center, c))
        triplet_weight.append(np.sqrt(wa * wc))
    triplets = np.asarray(triplets, np.int64).reshape(-1, 3)
    triplet_weight = np.asarray(triplet_weight, np.float32)
    audit = {
        "n_nodes": int(n),
        "n_edges": int(len(edges)),
        "n_triplets": int(len(triplets)),
        "neighbors": int(neighbors),
        "median_edge_distance": float(np.median(edge_dist)) if len(edge_dist) else None,
        "median_edge_weight": float(np.median(weight)) if len(weight) else None,
        "edge_sha256": hashlib.sha256(edges.astype(np.int64).tobytes()).hexdigest(),
    }
    return GraphData(edges, weight, triplets, triplet_weight, audit)


class SpotGraphProgram(nn.Module):
    def __init__(self, n_spots: int, n_types: int, n_genes: int, rank: int,
                 profiles: np.ndarray, composition: np.ndarray):
        super().__init__()
        self.w = nn.Parameter(torch.zeros(n_spots, n_types, rank))
        self.v = nn.Parameter(torch.randn(n_types, rank, n_genes) * 0.01)
        self.register_buffer("log_profiles", torch.log(torch.as_tensor(profiles) + EPS))
        self.register_buffer("composition", torch.as_tensor(composition, dtype=torch.float32))

    def centered_w(self) -> torch.Tensor:
        denom = self.composition.sum(dim=0).clamp_min(EPS)[:, None]
        mean = (self.w * self.composition[..., None]).sum(dim=0) / denom
        return self.w - mean[None]

    def forward(self, rows: torch.Tensor) -> torch.Tensor:
        local_w = self.centered_w()[rows]
        logits = self.log_profiles[None] + torch.einsum("btk,tkg->btg", local_w, self.v)
        type_probability = torch.softmax(logits, dim=-1)
        return torch.einsum("bt,btg->bg", self.composition[rows], type_probability).clamp_min(EPS)


class DirectGraphProgram(nn.Module):
    def __init__(self, n_nodes: int, n_types: int, n_genes: int, rank: int,
                 profiles: np.ndarray, class_id: np.ndarray):
        super().__init__()
        self.w = nn.Parameter(torch.zeros(n_nodes, rank))
        self.v = nn.Parameter(torch.randn(n_types, rank, n_genes) * 0.01)
        self.register_buffer("log_profiles", torch.log(torch.as_tensor(profiles) + EPS))
        self.register_buffer("class_id", torch.as_tensor(class_id, dtype=torch.long))

    def forward(self, rows: torch.Tensor) -> torch.Tensor:
        local_class = self.class_id[rows]
        logits = self.log_profiles[local_class] + torch.einsum(
            "bk,bkg->bg", self.w[rows], self.v[local_class],
        )
        return torch.softmax(logits, dim=-1).clamp_min(EPS)


def _trend_penalty(w: torch.Tensor, triplets: torch.Tensor,
                   weights: torch.Tensor, rng: np.random.Generator,
                   batch_size: int) -> torch.Tensor:
    if len(triplets) == 0:
        return w.sum() * 0.0
    take = min(batch_size, len(triplets))
    index = torch.as_tensor(rng.choice(len(triplets), take, replace=False),
                            dtype=torch.long, device=w.device)
    local = triplets[index]
    delta = w[local[:, 0]] - 2.0 * w[local[:, 1]] + w[local[:, 2]]
    shape = (take,) + (1,) * (delta.ndim - 1)
    return (delta.abs() * weights[index].reshape(shape)).mean()


def _finalize_program(w: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    w = np.asarray(w, np.float64)
    v = np.asarray(v, np.float64)
    if w.ndim == 2:
        # Direct-node W is shared as one tensor but its rows belong to different
        # types.  Per-type rescaling would require the class vector; keep the
        # fitted parameterization intact so W @ V remains exactly consistent.
        return w.astype(np.float32), v.astype(np.float32)
    for t in range(v.shape[0]):
        for k in range(v.shape[1]):
            norm = float(np.linalg.norm(v[t, k]))
            if norm <= EPS:
                continue
            v[t, k] /= norm
            if w.ndim == 3:
                w[:, t, k] *= norm
                pivot = int(np.argmax(np.abs(v[t, k])))
                if v[t, k, pivot] < 0:
                    v[t, k] *= -1
                    w[:, t, k] *= -1
    return w.astype(np.float32), v.astype(np.float32)


def fit_spot_program(
    counts: np.ndarray,
    composition: np.ndarray,
    profiles: np.ndarray,
    graph: GraphData,
    *,
    rank: int,
    lambda_trend: float,
    lambda_lasso: float = 0.02,
    lambda_l2: float = 1.0e-4,
    steps: int = 250,
    batch_size: int = 384,
    trend_batch_size: int = 4096,
    learning_rate: float = 0.03,
    eligible: np.ndarray | None = None,
    seed: int = 0,
    device: str = "cpu",
) -> ProgramFit:
    counts = np.asarray(counts, np.float32)
    composition = np.asarray(composition, np.float32)
    profiles = np.asarray(profiles, np.float32)
    if rank == 0:
        probability = composition @ profiles
        probability /= np.maximum(probability.sum(axis=1, keepdims=True), EPS)
        return ProgramFit(None, np.zeros((len(counts), composition.shape[1], 0), np.float32),
                          np.zeros((composition.shape[1], 0, counts.shape[1]), np.float32),
                          probability.astype(np.float32), {"rank": 0, "schema": "type_profile_only"})
    set_deterministic(seed)
    torch_device = torch.device(device)
    model = SpotGraphProgram(len(counts), composition.shape[1], counts.shape[1], rank,
                             profiles, composition).to(torch_device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    rng = np.random.default_rng(seed)
    library = counts.sum(axis=1)
    if eligible is None:
        eligible = np.flatnonzero(library > 0)
    else:
        eligible = np.asarray(eligible, np.int64)
        eligible = eligible[library[eligible] > 0]
    triplets = torch.as_tensor(graph.triplets, dtype=torch.long, device=torch_device)
    triplet_weight = torch.as_tensor(graph.triplet_weight, dtype=torch.float32, device=torch_device)
    history = []
    for step in range(steps):
        rows_np = rng.choice(eligible, min(batch_size, len(eligible)), replace=False)
        rows = torch.as_tensor(rows_np, dtype=torch.long, device=torch_device)
        y = torch.as_tensor(counts[rows_np], dtype=torch.float32, device=torch_device)
        probability = model(rows)
        nll = -(y * torch.log(probability)).sum() / y.sum().clamp_min(1.0)
        trend = _trend_penalty(model.w, triplets, triplet_weight, rng, trend_batch_size)
        penalty = lambda_trend * trend + lambda_lasso * model.v.abs().mean()
        penalty = penalty + lambda_l2 * (model.w.square().mean() + model.v.square().mean())
        loss = nll + penalty
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
        if step == 0 or (step + 1) % max(1, steps // 10) == 0:
            history.append({"step": int(step + 1), "loss": float(loss.detach().cpu()),
                            "nll": float(nll.detach().cpu()), "trend": float(trend.detach().cpu())})
    probabilities = predict_spot(model, batch_size=1024, device=device)
    w = model.centered_w().detach().cpu().numpy()
    v = model.v.detach().cpu().numpy()
    w, v = _finalize_program(w, v)
    audit = {"schema": "mv42_bayesgraph_spot_v1", "rank": int(rank),
             "lambda_trend": float(lambda_trend), "lambda_lasso": float(lambda_lasso),
             "steps": int(steps), "n_eligible": int(len(eligible)), "loss_history": history,
             "external_reference_used": False}
    return ProgramFit(model, w, v, probabilities, audit)


def predict_spot(model: SpotGraphProgram, batch_size: int = 1024,
                 device: str = "cpu") -> np.ndarray:
    model.eval()
    model.to(device)
    result = np.empty((len(model.w), model.v.shape[-1]), np.float32)
    with torch.no_grad():
        for start in range(0, len(result), batch_size):
            rows = torch.arange(start, min(start + batch_size, len(result)), device=device)
            result[start:start + len(rows)] = model(rows).cpu().numpy()
    return result


def fit_direct_program(
    counts: sparse.spmatrix,
    class_id: np.ndarray,
    profiles: np.ndarray,
    graph: GraphData,
    *,
    rank: int,
    lambda_trend: float,
    lambda_lasso: float = 0.02,
    lambda_l2: float = 1.0e-4,
    steps: int = 300,
    batch_size: int = 2048,
    trend_batch_size: int = 8192,
    learning_rate: float = 0.02,
    eligible: np.ndarray | None = None,
    seed: int = 0,
    device: str = "cpu",
) -> ProgramFit:
    counts = sparse.csr_matrix(counts, dtype=np.float32)
    class_id = np.asarray(class_id, np.int64)
    if rank == 0:
        probability = profiles[class_id].astype(np.float32)
        return ProgramFit(None, np.zeros((len(class_id), 0), np.float32),
                          np.zeros((profiles.shape[0], 0, profiles.shape[1]), np.float32),
                          probability, {"rank": 0, "schema": "class_profile_only"})
    set_deterministic(seed)
    torch_device = torch.device(device)
    model = DirectGraphProgram(len(class_id), profiles.shape[0], counts.shape[1], rank,
                               profiles, class_id).to(torch_device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    rng = np.random.default_rng(seed)
    library = np.asarray(counts.sum(axis=1)).ravel()
    if eligible is None:
        eligible = np.flatnonzero(library > 0)
    else:
        eligible = np.asarray(eligible, np.int64)
        eligible = eligible[library[eligible] > 0]
    triplets = torch.as_tensor(graph.triplets, dtype=torch.long, device=torch_device)
    triplet_weight = torch.as_tensor(graph.triplet_weight, dtype=torch.float32, device=torch_device)
    history = []
    for step in range(steps):
        rows_np = rng.choice(eligible, min(batch_size, len(eligible)), replace=False)
        rows = torch.as_tensor(rows_np, dtype=torch.long, device=torch_device)
        y = torch.as_tensor(counts[rows_np].toarray(), dtype=torch.float32, device=torch_device)
        probability = model(rows)
        nll = -(y * torch.log(probability)).sum() / y.sum().clamp_min(1.0)
        trend = _trend_penalty(model.w, triplets, triplet_weight, rng, trend_batch_size)
        loss = nll + lambda_trend * trend + lambda_lasso * model.v.abs().mean()
        loss = loss + lambda_l2 * (model.w.square().mean() + model.v.square().mean())
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
        if step == 0 or (step + 1) % max(1, steps // 10) == 0:
            history.append({"step": int(step + 1), "loss": float(loss.detach().cpu()),
                            "nll": float(nll.detach().cpu()), "trend": float(trend.detach().cpu())})
    probabilities = predict_direct(model, batch_size=4096, device=device)
    w = model.w.detach().cpu().numpy()
    v = model.v.detach().cpu().numpy()
    w, v = _finalize_program(w, v)
    audit = {"schema": "mv42_bayesgraph_direct_v1", "rank": int(rank),
             "lambda_trend": float(lambda_trend), "lambda_lasso": float(lambda_lasso),
             "steps": int(steps), "n_eligible": int(len(eligible)), "loss_history": history,
             "external_reference_used": False}
    return ProgramFit(model, w, v, probabilities, audit)


def predict_direct(model: DirectGraphProgram, batch_size: int = 4096,
                   device: str = "cpu") -> np.ndarray:
    model.eval()
    model.to(device)
    result = np.empty((len(model.w), model.v.shape[-1]), np.float32)
    with torch.no_grad():
        for start in range(0, len(result), batch_size):
            rows = torch.arange(start, min(start + batch_size, len(result)), device=device)
            result[start:start + len(rows)] = model(rows).cpu().numpy()
    return result


def dense_multinomial_nll(counts: np.ndarray, probability: np.ndarray,
                          rows: np.ndarray | None = None) -> float:
    counts = np.asarray(counts, np.float64)
    probability = np.asarray(probability, np.float64)
    if rows is not None:
        counts = counts[rows]
        probability = probability[rows]
    total = counts.sum()
    if total <= 0:
        return float("nan")
    return float(-(counts * np.log(np.maximum(probability, EPS))).sum() / total)


def sparse_multinomial_nll(counts: sparse.spmatrix, probability: np.ndarray) -> float:
    counts = sparse.coo_matrix(counts, dtype=np.float64)
    total = counts.data.sum()
    if total <= 0:
        return float("nan")
    return float(-(counts.data * np.log(np.maximum(probability[counts.row, counts.col], EPS))).sum() / total)


def split_sparse_counts(counts: sparse.spmatrix, seed: int) -> tuple[sparse.csr_matrix, sparse.csr_matrix]:
    counts = sparse.csr_matrix(counts)
    if np.any(counts.data != np.rint(counts.data)):
        raise ValueError("binomial thinning requires integer-valued counts")
    rng = np.random.default_rng(seed)
    left_data = rng.binomial(np.rint(counts.data).astype(np.int64), 0.5)
    left = sparse.csr_matrix((left_data, counts.indices.copy(), counts.indptr.copy()), shape=counts.shape)
    left.eliminate_zeros()
    right = counts.astype(np.int64) - left
    right.eliminate_zeros()
    return left.astype(np.int32), right.astype(np.int32)


def posterior_shrinkage(counts: sparse.spmatrix, prior: np.ndarray, tau: float) -> np.ndarray:
    counts = sparse.csr_matrix(counts, dtype=np.float64)
    prior = np.asarray(prior, np.float64)
    library = np.asarray(counts.sum(axis=1)).ravel()
    result = counts.toarray() + float(tau) * prior
    denominator = library + float(tau)
    zero = denominator <= 0
    denominator[zero] = 1.0
    result /= denominator[:, None]
    if np.any(zero):
        result[zero] = prior[zero]
    result /= np.maximum(result.sum(axis=1, keepdims=True), EPS)
    return result.astype(np.float32)
