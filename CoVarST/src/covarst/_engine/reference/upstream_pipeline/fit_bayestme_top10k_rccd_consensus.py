"""Fit a reconstruction-constrained consensus dictionary for 10K BayesTME.

The per-slide BayesTME factorization is frozen::

    P_local[s] = theta_rna[s] @ basis_probability[s]

This script learns a batch-neutral consensus reference W and sparse slide-local
maps A_s such that::

    P_calibrated[s] = theta_rna[s] @ A_s @ W_batch[s]

matches P_local[s].  W_batch is a strongly shrunk technology plus
source-by-technology calibration of the single shared biological dictionary W.
Patient-level holdout validation selects the smallest candidate dictionary that
meets reconstruction gates; the selected model is then refit on all slides.

No upstream BayesTME result is modified.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import torch
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA


EPS = 1.0e-8
METHOD = "bayestme_top10k_reconstruction_constrained_consensus_dictionary_v1"


def _atomic_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _normalise_rows(values: np.ndarray) -> np.ndarray:
    values = np.maximum(np.asarray(values, dtype=np.float64), 0.0)
    return values / np.maximum(values.sum(axis=-1, keepdims=True), EPS)


def _sha256(path: Path, block: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(block)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _quantiles(values: Iterable[float]) -> dict[str, float]:
    array = np.asarray(list(values), dtype=np.float64)
    return {
        "min": float(np.min(array)),
        "q05": float(np.quantile(array, 0.05)),
        "q25": float(np.quantile(array, 0.25)),
        "median": float(np.median(array)),
        "q75": float(np.quantile(array, 0.75)),
        "q95": float(np.quantile(array, 0.95)),
        "max": float(np.max(array)),
    }


@dataclass
class Bundle:
    basis: np.ndarray
    sqrt_basis: np.ndarray
    theta_gram: np.ndarray
    local_mass: np.ndarray
    stability: np.ndarray
    program_weight: np.ndarray
    core_mask: np.ndarray
    valid_mask: np.ndarray
    slide_weight: np.ndarray
    batch_index: np.ndarray
    technology_index: np.ndarray
    genes: np.ndarray
    gene_weight: np.ndarray
    slides: pd.DataFrame
    programs: pd.DataFrame
    theta: list[np.ndarray]
    library_size: list[np.ndarray]
    global_rows: list[np.ndarray]
    batch_names: list[str]
    technology_names: list[str]
    batch_to_technology: np.ndarray
    technology_prior: np.ndarray
    batch_prior_within_technology: np.ndarray

    @property
    def n_slides(self) -> int:
        return int(self.basis.shape[0])

    @property
    def k_max(self) -> int:
        return int(self.basis.shape[1])

    @property
    def n_genes(self) -> int:
        return int(self.basis.shape[2])


def _posterior_stability(slide_dir: Path, basis_mean: np.ndarray) -> np.ndarray:
    path = slide_dir / "expression_trace.npy"
    if not path.exists():
        return np.ones(basis_mean.shape[0], dtype=np.float32)
    draws = np.load(path, allow_pickle=False).astype(np.float32)
    draws /= np.maximum(draws.sum(axis=2, keepdims=True), EPS)
    root_draws = np.sqrt(np.maximum(draws, 0.0))
    root_mean = np.sqrt(np.maximum(basis_mean, 0.0))[None, :, :]
    hellinger_sq = 0.5 * np.square(root_draws - root_mean).sum(axis=2)
    uncertainty = hellinger_sq.mean(axis=0)
    stability = 1.0 / (1.0 + 20.0 * uncertainty)
    return np.clip(stability, 0.25, 1.0).astype(np.float32)


def _core_programs(mass: np.ndarray, cumulative: float = 0.99) -> np.ndarray:
    order = np.argsort(-mass, kind="stable")
    cumulative_mass = np.cumsum(mass[order])
    count = int(np.searchsorted(cumulative_mass, cumulative, side="left") + 1)
    result = np.zeros(len(mass), dtype=bool)
    result[order[:count]] = True
    result |= mass >= 0.02
    return result


def _balanced_slide_weights(slides: pd.DataFrame) -> np.ndarray:
    result = np.zeros(len(slides), dtype=np.float64)
    sources = sorted(slides["source"].astype(str).unique())
    for source in sources:
        source_rows = slides.index[slides["source"].astype(str).eq(source)].to_numpy()
        patients = sorted(slides.loc[source_rows, "patient"].astype(str).unique())
        for patient in patients:
            rows = slides.index[
                slides["source"].astype(str).eq(source)
                & slides["patient"].astype(str).eq(patient)
            ].to_numpy()
            result[rows] = 1.0 / (len(sources) * len(patients) * len(rows))
    result /= result.sum()
    return result.astype(np.float32)


def _gene_weights(
    basis: np.ndarray,
    valid: np.ndarray,
    program_weight: np.ndarray,
    slide_weight: np.ndarray,
    marker_genes: int = 100,
) -> np.ndarray:
    flat = basis[valid]
    flat_program_weight = program_weight[valid]
    repeated_slide_weight = np.repeat(slide_weight[:, None], valid.shape[1], axis=1)[valid]
    weights = np.maximum(flat_program_weight * repeated_slide_weight, EPS)
    weights /= weights.sum()
    mean_probability = np.sum(flat.astype(np.float64) * weights[:, None], axis=0)
    abundance_weight = 1.0 / np.sqrt(mean_probability + 1.0e-7)
    lower, upper = np.quantile(abundance_weight, [0.10, 0.90])
    abundance_weight = np.clip(abundance_weight, lower, upper)
    marker_genes = min(int(marker_genes), flat.shape[1])
    top = np.argpartition(flat, -marker_genes, axis=1)[:, -marker_genes:]
    document_frequency = np.bincount(top.ravel(), minlength=flat.shape[1])
    idf = np.log((len(flat) + 1.0) / (document_frequency + 1.0)) + 1.0
    idf = np.clip(idf, 1.0, np.quantile(idf, 0.95))
    result = abundance_weight * idf
    result /= result.mean()
    return result.astype(np.float32)


def load_bundle(
    queue_dir: Path,
    master_spot_index: Path,
    *,
    use_posterior_stability: bool,
    max_slides: int | None,
) -> Bundle:
    genes = np.load(queue_dir / "gene_order.npy", allow_pickle=False).astype(str)
    complete_paths = sorted((queue_dir / "slides").glob("*/complete.json"))
    if max_slides is not None:
        complete_paths = complete_paths[: int(max_slides)]
    if not complete_paths:
        raise RuntimeError("no completed per-slide BayesTME results")
    master = pd.read_csv(
        master_spot_index,
        usecols=[
            "slide_index",
            "sample_id",
            "source",
            "patient",
            "technology",
        ],
    ).drop_duplicates("slide_index")
    master = master.set_index("slide_index", drop=False)
    slide_records: list[dict[str, object]] = []
    program_records: list[pd.DataFrame] = []
    basis_blocks: list[np.ndarray] = []
    gram_blocks: list[np.ndarray] = []
    mass_blocks: list[np.ndarray] = []
    stability_blocks: list[np.ndarray] = []
    core_blocks: list[np.ndarray] = []
    theta_blocks: list[np.ndarray] = []
    library_blocks: list[np.ndarray] = []
    global_row_blocks: list[np.ndarray] = []
    for slide_position, complete_path in enumerate(complete_paths):
        slide_dir = complete_path.parent
        complete = json.loads(complete_path.read_text(encoding="utf-8"))
        slide_index = int(complete["slide_index"])
        if slide_index not in master.index:
            raise RuntimeError(f"slide {slide_index} absent from master spot index")
        metadata = master.loc[slide_index]
        basis = np.load(slide_dir / "basis_probability_mean.npy", allow_pickle=False).astype(np.float32)
        theta = np.load(slide_dir / "theta_rna_mean.npy", allow_pickle=False).astype(np.float32)
        local_index = pd.read_csv(slide_dir / "local_program_index.csv")
        spot_index = pd.read_csv(slide_dir / "spot_index.csv", usecols=["global_index"])
        library_size = np.load(slide_dir / "library_size.npy", allow_pickle=False).astype(np.float32)
        if basis.shape != (theta.shape[1], len(genes)):
            raise RuntimeError(f"{slide_dir.name}: basis/theta/gene mismatch")
        if len(local_index) != basis.shape[0] or len(library_size) != theta.shape[0]:
            raise RuntimeError(f"{slide_dir.name}: local metadata mismatch")
        if float(np.max(np.abs(basis.sum(axis=1) - 1.0))) > 1.0e-5:
            raise RuntimeError(f"{slide_dir.name}: basis rows do not sum to one")
        if float(np.max(np.abs(theta.sum(axis=1) - 1.0))) > 1.0e-5:
            raise RuntimeError(f"{slide_dir.name}: theta rows do not sum to one")
        mass = local_index["rna_mass_fraction"].to_numpy(dtype=np.float32)
        stability = (
            _posterior_stability(slide_dir, basis)
            if use_posterior_stability
            else np.ones(len(mass), dtype=np.float32)
        )
        local_index = local_index.copy()
        local_index["slide_position"] = slide_position
        local_index["slide_dir"] = str(slide_dir)
        local_index["patient"] = str(metadata["patient"])
        local_index["technology"] = str(metadata["technology"])
        local_index["batch"] = f"{metadata['source']}|{metadata['technology']}"
        local_index["posterior_stability"] = stability
        program_records.append(local_index)
        basis_blocks.append(basis)
        gram_blocks.append((theta.T @ theta / max(len(theta), 1)).astype(np.float32))
        mass_blocks.append(mass)
        stability_blocks.append(stability)
        core_blocks.append(_core_programs(mass))
        theta_blocks.append(theta)
        library_blocks.append(library_size)
        global_row_blocks.append(spot_index["global_index"].to_numpy(dtype=np.int64))
        slide_records.append(
            {
                "slide_position": slide_position,
                "slide_index": slide_index,
                "sample_id": str(metadata["sample_id"]),
                "source": str(metadata["source"]),
                "patient": str(metadata["patient"]),
                "technology": str(metadata["technology"]),
                "batch": f"{metadata['source']}|{metadata['technology']}",
                "slide_dir": str(slide_dir),
                "spots": int(theta.shape[0]),
                "local_k": int(theta.shape[1]),
            }
        )
    slides = pd.DataFrame(slide_records)
    programs = pd.concat(program_records, ignore_index=True)
    programs["global_local_program_row"] = np.arange(len(programs), dtype=np.int64)
    k_max = max(block.shape[0] for block in basis_blocks)
    n_slides = len(basis_blocks)
    n_genes = len(genes)
    basis_pad = np.zeros((n_slides, k_max, n_genes), dtype=np.float32)
    gram_pad = np.zeros((n_slides, k_max, k_max), dtype=np.float32)
    mass_pad = np.zeros((n_slides, k_max), dtype=np.float32)
    stability_pad = np.zeros_like(mass_pad)
    core_pad = np.zeros_like(mass_pad, dtype=bool)
    valid_pad = np.zeros_like(mass_pad, dtype=bool)
    for slide, basis in enumerate(basis_blocks):
        k = len(basis)
        basis_pad[slide, :k] = basis
        gram_pad[slide, :k, :k] = gram_blocks[slide]
        mass_pad[slide, :k] = mass_blocks[slide]
        stability_pad[slide, :k] = stability_blocks[slide]
        core_pad[slide, :k] = core_blocks[slide]
        valid_pad[slide, :k] = True
    slide_weight = _balanced_slide_weights(slides)
    program_weight = mass_pad * stability_pad
    gene_weight = _gene_weights(
        basis_pad,
        valid_pad,
        program_weight,
        slide_weight,
    )
    technology_names = sorted(slides["technology"].unique())
    technology_map = {name: index for index, name in enumerate(technology_names)}
    batch_names = sorted(slides["batch"].unique())
    batch_map = {name: index for index, name in enumerate(batch_names)}
    technology_index = slides["technology"].map(technology_map).to_numpy(dtype=np.int64)
    batch_index = slides["batch"].map(batch_map).to_numpy(dtype=np.int64)
    batch_to_technology = np.asarray(
        [technology_map[name.split("|", 1)[1]] for name in batch_names],
        dtype=np.int64,
    )
    technology_prior = np.zeros(len(technology_names), dtype=np.float64)
    batch_prior = np.zeros(len(batch_names), dtype=np.float64)
    for slide, weight in enumerate(slide_weight.astype(np.float64)):
        technology_prior[technology_index[slide]] += weight
        batch_prior[batch_index[slide]] += weight
    technology_prior /= technology_prior.sum()
    batch_prior_within = np.zeros_like(batch_prior)
    for technology in range(len(technology_names)):
        rows = np.flatnonzero(batch_to_technology == technology)
        total = batch_prior[rows].sum()
        if total > 0:
            batch_prior_within[rows] = batch_prior[rows] / total
    return Bundle(
        basis=basis_pad,
        sqrt_basis=np.sqrt(np.maximum(basis_pad, 0.0)).astype(np.float32),
        theta_gram=gram_pad,
        local_mass=mass_pad,
        stability=stability_pad,
        program_weight=program_weight,
        core_mask=core_pad,
        valid_mask=valid_pad,
        slide_weight=slide_weight,
        batch_index=batch_index,
        technology_index=technology_index,
        genes=genes,
        gene_weight=gene_weight,
        slides=slides,
        programs=programs,
        theta=theta_blocks,
        library_size=library_blocks,
        global_rows=global_row_blocks,
        batch_names=batch_names,
        technology_names=technology_names,
        batch_to_technology=batch_to_technology,
        technology_prior=technology_prior.astype(np.float32),
        batch_prior_within_technology=batch_prior_within.astype(np.float32),
    )


def _patient_holdout(slides: pd.DataFrame, fraction: float, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(int(seed))
    validation_patients: set[str] = set()
    for source in sorted(slides["source"].unique()):
        patients = sorted(slides.loc[slides["source"].eq(source), "patient"].astype(str).unique())
        shuffled = np.asarray(patients, dtype=object)
        rng.shuffle(shuffled)
        count = max(1, int(round(len(shuffled) * float(fraction))))
        count = min(count, max(len(shuffled) - 1, 1))
        validation_patients.update(str(value) for value in shuffled[:count])
    validation = slides["patient"].astype(str).isin(validation_patients).to_numpy()
    # A batch represented by a single patient must remain in training.
    for batch in sorted(slides["batch"].unique()):
        rows = slides.index[slides["batch"].eq(batch)].to_numpy()
        if np.all(validation[rows]):
            patient = str(slides.loc[rows[0], "patient"])
            validation[slides["patient"].astype(str).eq(patient).to_numpy()] = False
    train_rows = np.flatnonzero(~validation)
    validation_rows = np.flatnonzero(validation)
    if len(validation_rows) == 0:
        raise RuntimeError("patient holdout produced no validation slides")
    return train_rows.astype(np.int64), validation_rows.astype(np.int64)


def _initial_reference(
    bundle: Bundle,
    slide_rows: np.ndarray,
    programs: int,
    seed: int,
    fallback_slide_rows: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    initialization_rows = slide_rows
    selected_valid = bundle.valid_mask[initialization_rows]
    if int(selected_valid.sum()) < int(programs) and fallback_slide_rows is not None:
        fallback_slide_rows = np.asarray(fallback_slide_rows, dtype=np.int64)
        fallback_valid = bundle.valid_mask[fallback_slide_rows]
        if int(fallback_valid.sum()) >= int(programs):
            initialization_rows = fallback_slide_rows
            selected_valid = fallback_valid
            print(
                json.dumps(
                    {
                        "stage": "initial_reference_fallback",
                        "reason": "requested_programs_exceed_training_local_programs",
                        "requested_programs": int(programs),
                        "training_local_programs": int(bundle.valid_mask[slide_rows].sum()),
                        "initialization_local_programs": int(selected_valid.sum()),
                        "initialization_slides": int(len(initialization_rows)),
                    }
                ),
                flush=True,
            )
    if int(selected_valid.sum()) < int(programs):
        raise ValueError(
            f"cannot initialise {programs} consensus programs from "
            f"{int(selected_valid.sum())} available local programs"
        )
    selected_core = bundle.core_mask[initialization_rows] & selected_valid
    flat_basis = bundle.basis[initialization_rows][selected_valid]
    flat_core = selected_core[selected_valid]
    flat_program_weight = bundle.program_weight[initialization_rows][selected_valid]
    repeated_slide_weight = np.repeat(
        bundle.slide_weight[initialization_rows, None], bundle.k_max, axis=1
    )[selected_valid]
    fit_rows = np.flatnonzero(flat_core)
    if len(fit_rows) < programs:
        fit_rows = np.arange(len(flat_basis))
    fit_basis = flat_basis[fit_rows]
    fit_weight = np.maximum(flat_program_weight[fit_rows] * repeated_slide_weight[fit_rows], EPS)
    features = np.sqrt(np.maximum(fit_basis, 0.0)).astype(np.float32)
    components = min(50, features.shape[0] - 1, features.shape[1] - 1)
    embedding = PCA(
        n_components=components,
        svd_solver="randomized",
        random_state=int(seed),
    ).fit_transform(features)
    clusterer = MiniBatchKMeans(
        n_clusters=int(programs),
        random_state=int(seed),
        batch_size=min(512, len(embedding)),
        n_init=10,
        max_iter=300,
        reassignment_ratio=0.01,
    )
    labels = clusterer.fit_predict(embedding, sample_weight=fit_weight)
    reference = np.zeros((programs, bundle.n_genes), dtype=np.float64)
    reference_weight = np.zeros(programs, dtype=np.float64)
    for row, label in enumerate(labels):
        weight = float(fit_weight[row])
        reference[int(label)] += fit_basis[row] * weight
        reference_weight[int(label)] += weight
    rng = np.random.default_rng(int(seed) + 17)
    for label in range(programs):
        if reference_weight[label] <= 0:
            replacement = int(rng.choice(len(fit_basis), p=fit_weight / fit_weight.sum()))
            reference[label] = fit_basis[replacement]
            reference_weight[label] = 1.0
    reference /= np.maximum(reference.sum(axis=1, keepdims=True), EPS)
    assignment_logits = _mapping_logits(bundle.basis[slide_rows], bundle.valid_mask[slide_rows], reference)
    return reference.astype(np.float32), assignment_logits


def _mapping_logits(basis: np.ndarray, valid: np.ndarray, reference: np.ndarray) -> np.ndarray:
    result = np.full((basis.shape[0], basis.shape[1], len(reference)), -4.0, dtype=np.float32)
    root_reference = np.sqrt(np.maximum(reference, 0.0)).astype(np.float32)
    for slide in range(len(basis)):
        rows = np.flatnonzero(valid[slide])
        similarity = np.sqrt(np.maximum(basis[slide, rows], 0.0)) @ root_reference.T
        labels = np.argmax(similarity, axis=1)
        result[slide, rows, labels] = 4.0
    return result


class RCCD(torch.nn.Module):
    def __init__(
        self,
        reference: np.ndarray,
        mapping_logits: np.ndarray,
        technology_count: int,
        batch_count: int,
        batch_to_technology: np.ndarray,
        technology_prior: np.ndarray,
        batch_prior_within: np.ndarray,
        *,
        technology_bound: float,
        batch_bound: float,
    ) -> None:
        super().__init__()
        self.w_logits = torch.nn.Parameter(
            torch.from_numpy(np.log(np.maximum(reference, 1.0e-10))).float()
        )
        self.a_logits = torch.nn.Parameter(torch.from_numpy(mapping_logits).float())
        self.technology_raw = torch.nn.Parameter(torch.zeros(technology_count, reference.shape[1]))
        self.batch_raw = torch.nn.Parameter(torch.zeros(batch_count, reference.shape[1]))
        self.register_buffer("batch_to_technology", torch.from_numpy(batch_to_technology).long())
        self.register_buffer("technology_prior", torch.from_numpy(technology_prior).float())
        self.register_buffer("batch_prior_within", torch.from_numpy(batch_prior_within).float())
        self.technology_bound = float(technology_bound)
        self.batch_bound = float(batch_bound)

    def mapping(self, temperature: float) -> torch.Tensor:
        return torch.softmax(self.a_logits / float(temperature), dim=-1)

    def neutral_reference(self) -> torch.Tensor:
        return torch.softmax(self.w_logits, dim=-1)

    def deltas(self) -> tuple[torch.Tensor, torch.Tensor]:
        technology = self.technology_bound * torch.tanh(self.technology_raw)
        technology = technology - torch.sum(
            technology * self.technology_prior[:, None], dim=0, keepdim=True
        )
        batch = self.batch_bound * torch.tanh(self.batch_raw)
        centered = torch.zeros_like(batch)
        for technology_index in range(len(self.technology_prior)):
            rows = torch.nonzero(
                self.batch_to_technology.eq(technology_index),
                as_tuple=False,
            ).flatten()
            if len(rows) <= 1:
                centered[rows] = 0.0
                continue
            prior = self.batch_prior_within[rows]
            mean = torch.sum(batch[rows] * prior[:, None], dim=0, keepdim=True)
            centered[rows] = batch[rows] - mean
        return technology, centered

    def reference_by_batch(self) -> torch.Tensor:
        technology, batch = self.deltas()
        adapted = []
        for batch_index in range(len(self.batch_to_technology)):
            technology_index = int(self.batch_to_technology[batch_index].item())
            adapted.append(
                torch.softmax(
                    self.w_logits + technology[technology_index] + batch[batch_index],
                    dim=-1,
                )
            )
        return torch.stack(adapted, dim=0)


@dataclass
class TensorBundle:
    basis: torch.Tensor
    sqrt_basis: torch.Tensor
    log_basis: torch.Tensor
    gram: torch.Tensor
    program_weight: torch.Tensor
    core_mask: torch.Tensor
    valid_mask: torch.Tensor
    slide_weight: torch.Tensor
    batch_index: torch.Tensor
    gene_weight: torch.Tensor
    energy: torch.Tensor


def _tensor_bundle(bundle: Bundle, slide_rows: np.ndarray, device: torch.device) -> TensorBundle:
    basis = torch.from_numpy(bundle.basis[slide_rows]).to(device)
    gram = torch.from_numpy(bundle.theta_gram[slide_rows]).to(device)
    gene_weight = torch.from_numpy(bundle.gene_weight).to(device)
    q_basis = torch.einsum("sik,skg->sig", gram, basis)
    energy = torch.sum(basis * q_basis * gene_weight[None, None, :], dim=(1, 2))
    return TensorBundle(
        basis=basis,
        sqrt_basis=torch.from_numpy(bundle.sqrt_basis[slide_rows]).to(device),
        log_basis=torch.log1p(1.0e4 * basis),
        gram=gram,
        program_weight=torch.from_numpy(bundle.program_weight[slide_rows]).to(device),
        core_mask=torch.from_numpy(bundle.core_mask[slide_rows]).to(device),
        valid_mask=torch.from_numpy(bundle.valid_mask[slide_rows]).to(device),
        slide_weight=torch.from_numpy(bundle.slide_weight[slide_rows]).to(device),
        batch_index=torch.from_numpy(bundle.batch_index[slide_rows]).to(device),
        gene_weight=gene_weight,
        energy=torch.clamp(energy, min=EPS),
    )


def _mapping_regularizers(
    mapping: torch.Tensor,
    data: TensorBundle,
) -> tuple[torch.Tensor, torch.Tensor]:
    valid_weight = data.program_weight * data.valid_mask.float()
    valid_weight = valid_weight / torch.clamp(valid_weight.sum(dim=1, keepdim=True), min=EPS)
    entropy = -torch.sum(mapping * torch.log(torch.clamp(mapping, min=EPS)), dim=-1)
    entropy = entropy / math.log(mapping.shape[-1])
    entropy_slide = torch.sum(entropy * valid_weight, dim=1)
    entropy_loss = torch.sum(entropy_slide * data.slide_weight) / torch.clamp(data.slide_weight.sum(), min=EPS)
    core_weight = data.program_weight * data.core_mask.float()
    pair_weight = torch.sqrt(
        torch.clamp(core_weight[:, :, None] * core_weight[:, None, :], min=0.0)
    )
    identity = torch.eye(mapping.shape[1], device=mapping.device)[None, :, :]
    pair_weight = pair_weight * (1.0 - identity)
    overlap = torch.einsum("skg,slg->skl", mapping, mapping)
    collision_slide = torch.sum(overlap * pair_weight, dim=(1, 2)) / torch.clamp(
        torch.sum(pair_weight, dim=(1, 2)), min=EPS
    )
    collision_loss = torch.sum(collision_slide * data.slide_weight) / torch.clamp(
        data.slide_weight.sum(), min=EPS
    )
    return entropy_loss, collision_loss


def _reconstruction_loss(
    mapping: torch.Tensor,
    reference_by_batch: torch.Tensor,
    neutral_reference: torch.Tensor,
    data: TensorBundle,
    technology_delta: torch.Tensor | None,
    batch_delta: torch.Tensor | None,
    *,
    lambda_basis: float,
    lambda_log_basis: float,
    lambda_entropy: float,
    lambda_collision: float,
    lambda_diversity: float,
    lambda_technology: float,
    lambda_batch: float,
) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
    spot_slide = torch.zeros(len(data.basis), device=data.basis.device)
    basis_slide = torch.zeros_like(spot_slide)
    log_basis_slide = torch.zeros_like(spot_slide)
    for batch_index in torch.unique(data.batch_index).tolist():
        rows = torch.nonzero(
            data.batch_index.eq(int(batch_index)),
            as_tuple=False,
        ).flatten()
        predicted = torch.einsum(
            "skp,pg->skg", mapping[rows], reference_by_batch[int(batch_index)]
        )
        error = data.basis[rows] - predicted
        q_error = torch.einsum("sik,skg->sig", data.gram[rows], error)
        numerator = torch.sum(
            error * q_error * data.gene_weight[None, None, :], dim=(1, 2)
        )
        spot_slide[rows] = numerator / data.energy[rows]
        hellinger = 0.5 * torch.sum(
            torch.square(data.sqrt_basis[rows] - torch.sqrt(torch.clamp(predicted, min=0.0))),
            dim=2,
        )
        weights = data.program_weight[rows] * data.valid_mask[rows].float()
        weights = weights / torch.clamp(weights.sum(dim=1, keepdim=True), min=EPS)
        basis_slide[rows] = torch.sum(hellinger * weights, dim=1)
        log_error = torch.nn.functional.smooth_l1_loss(
            torch.log1p(1.0e4 * predicted),
            data.log_basis[rows],
            reduction="none",
            beta=0.25,
        )
        log_program = torch.sum(
            log_error * data.gene_weight[None, None, :],
            dim=2,
        ) / torch.clamp(torch.sum(data.gene_weight), min=EPS)
        log_basis_slide[rows] = torch.sum(log_program * weights, dim=1)
    slide_weight = data.slide_weight / torch.clamp(data.slide_weight.sum(), min=EPS)
    spot_loss = torch.sum(spot_slide * slide_weight)
    basis_loss = torch.sum(basis_slide * slide_weight)
    log_basis_loss = torch.sum(log_basis_slide * slide_weight)
    entropy_loss, collision_loss = _mapping_regularizers(mapping, data)
    root_reference = torch.sqrt(torch.clamp(neutral_reference, min=0.0))
    similarity = root_reference @ root_reference.T
    off_diagonal = ~torch.eye(len(similarity), dtype=bool, device=similarity.device)
    diversity_loss = torch.mean(torch.relu(similarity[off_diagonal] - 0.98) ** 2)
    technology_loss = (
        torch.mean(technology_delta**2)
        if technology_delta is not None
        else torch.zeros((), device=data.basis.device)
    )
    batch_loss = (
        torch.mean(batch_delta**2)
        if batch_delta is not None
        else torch.zeros((), device=data.basis.device)
    )
    total = (
        spot_loss
        + float(lambda_basis) * basis_loss
        + float(lambda_log_basis) * log_basis_loss
        + float(lambda_entropy) * entropy_loss
        + float(lambda_collision) * collision_loss
        + float(lambda_diversity) * diversity_loss
        + float(lambda_technology) * technology_loss
        + float(lambda_batch) * batch_loss
    )
    return total, {
        "spot": spot_loss,
        "basis": basis_loss,
        "log_basis": log_basis_loss,
        "entropy": entropy_loss,
        "collision": collision_loss,
        "diversity": diversity_loss,
        "technology": technology_loss,
        "batch": batch_loss,
    }


def _fit_model(
    bundle: Bundle,
    slide_rows: np.ndarray,
    reference_init: np.ndarray,
    mapping_init: np.ndarray,
    *,
    device: torch.device,
    epochs: int,
    warmup_epochs: int,
    learning_rate: float,
    seed: int,
    args: argparse.Namespace,
) -> tuple[RCCD, list[dict[str, float]]]:
    torch.manual_seed(int(seed))
    model = RCCD(
        reference_init,
        mapping_init,
        len(bundle.technology_names),
        len(bundle.batch_names),
        bundle.batch_to_technology,
        bundle.technology_prior,
        bundle.batch_prior_within_technology,
        technology_bound=float(args.technology_bound),
        batch_bound=float(args.batch_bound),
    ).to(device)
    data = _tensor_bundle(bundle, slide_rows, device)
    optimizer = torch.optim.Adam(model.parameters(), lr=float(learning_rate))
    history: list[dict[str, float]] = []
    for epoch in range(int(epochs)):
        optimizer.zero_grad(set_to_none=True)
        temperature = max(0.35, 1.0 - 0.65 * epoch / max(int(epochs) - 1, 1))
        mapping = model.mapping(temperature)
        neutral = model.neutral_reference()
        if epoch < int(warmup_epochs):
            references = neutral[None, :, :].repeat(len(bundle.batch_names), 1, 1)
            technology_delta = None
            batch_delta = None
        else:
            references = model.reference_by_batch()
            technology_delta, batch_delta = model.deltas()
        total, pieces = _reconstruction_loss(
            mapping,
            references,
            neutral,
            data,
            technology_delta,
            batch_delta,
            lambda_basis=float(args.lambda_basis),
            lambda_log_basis=float(args.lambda_log_basis),
            lambda_entropy=float(args.lambda_entropy),
            lambda_collision=float(args.lambda_collision),
            lambda_diversity=float(args.lambda_diversity),
            lambda_technology=float(args.lambda_technology),
            lambda_batch=float(args.lambda_batch),
        )
        if not torch.isfinite(total):
            raise FloatingPointError(f"non-finite loss at epoch {epoch}")
        total.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        optimizer.step()
        record = {"epoch": float(epoch + 1), "total": float(total.detach().cpu())}
        record.update({name: float(value.detach().cpu()) for name, value in pieces.items()})
        history.append(record)
        if epoch == 0 or (epoch + 1) % 25 == 0 or epoch + 1 == int(epochs):
            print(json.dumps({"fit": record}, ensure_ascii=False), flush=True)
    return model, history


def _fit_mapping_only(
    bundle: Bundle,
    slide_rows: np.ndarray,
    neutral_reference: np.ndarray,
    reference_by_batch: np.ndarray,
    *,
    device: torch.device,
    epochs: int,
    learning_rate: float,
    seed: int,
    args: argparse.Namespace,
) -> np.ndarray:
    torch.manual_seed(int(seed))
    mapping_init = _mapping_logits(
        bundle.basis[slide_rows], bundle.valid_mask[slide_rows], neutral_reference
    )
    logits = torch.nn.Parameter(torch.from_numpy(mapping_init).to(device))
    optimizer = torch.optim.Adam([logits], lr=float(learning_rate))
    data = _tensor_bundle(bundle, slide_rows, device)
    references = torch.from_numpy(reference_by_batch).to(device)
    neutral = torch.from_numpy(neutral_reference).to(device)
    for epoch in range(int(epochs)):
        optimizer.zero_grad(set_to_none=True)
        temperature = max(0.35, 1.0 - 0.65 * epoch / max(int(epochs) - 1, 1))
        mapping = torch.softmax(logits / temperature, dim=-1)
        total, _ = _reconstruction_loss(
            mapping,
            references,
            neutral,
            data,
            None,
            None,
            lambda_basis=float(args.lambda_basis),
            lambda_log_basis=float(args.lambda_log_basis),
            lambda_entropy=float(args.lambda_entropy),
            lambda_collision=float(args.lambda_collision),
            lambda_diversity=0.0,
            lambda_technology=0.0,
            lambda_batch=0.0,
        )
        total.backward()
        torch.nn.utils.clip_grad_norm_([logits], 5.0)
        optimizer.step()
    return torch.softmax(logits / 0.35, dim=-1).detach().cpu().numpy().astype(np.float32)


def _model_arrays(
    model: RCCD,
    temperature: float = 0.35,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    model.eval()
    with torch.no_grad():
        neutral = model.neutral_reference().cpu().numpy().astype(np.float32)
        by_batch = model.reference_by_batch().cpu().numpy().astype(np.float32)
        technology, batch = model.deltas()
        mapping = model.mapping(float(temperature)).cpu().numpy().astype(np.float32)
    return (
        neutral,
        by_batch,
        technology.cpu().numpy().astype(np.float32),
        batch.cpu().numpy().astype(np.float32),
        mapping,
    )


def _n95(values: np.ndarray) -> int:
    order = np.sort(np.asarray(values, dtype=np.float64))[::-1]
    return int(np.searchsorted(np.cumsum(order), 0.95, side="left") + 1)


def _effective(values: np.ndarray) -> np.ndarray:
    values = np.maximum(np.asarray(values, dtype=np.float64), EPS)
    values /= values.sum(axis=1, keepdims=True)
    return np.exp(-np.sum(values * np.log(values), axis=1))


def _evaluate_slide(
    theta: np.ndarray,
    basis: np.ndarray,
    mapping: np.ndarray,
    reference: np.ndarray,
    library_size: np.ndarray,
    gene_weight: np.ndarray,
    *,
    block: int,
) -> dict[str, float]:
    predicted_basis = mapping @ reference
    gram = theta.T @ theta / max(len(theta), 1)
    error = basis - predicted_basis
    numerator = np.sum((error * (gram @ error)) * gene_weight[None, :])
    denominator = np.sum((basis * (gram @ basis)) * gene_weight[None, :])
    relative = math.sqrt(max(float(numerator / max(denominator, EPS)), 0.0))
    spot_dot = np.zeros(len(theta), dtype=np.float64)
    spot_local_sq = np.zeros(len(theta), dtype=np.float64)
    spot_global_sq = np.zeros(len(theta), dtype=np.float64)
    gene_correlations: list[np.ndarray] = []
    count_correlations: list[np.ndarray] = []
    for start in range(0, basis.shape[1], int(block)):
        stop = min(start + int(block), basis.shape[1])
        local_probability = theta @ basis[:, start:stop]
        global_probability = theta @ predicted_basis[:, start:stop]
        local_log = np.log1p(1.0e4 * local_probability)
        global_log = np.log1p(1.0e4 * global_probability)
        spot_dot += np.sum(local_log * global_log, axis=1)
        spot_local_sq += np.sum(local_log**2, axis=1)
        spot_global_sq += np.sum(global_log**2, axis=1)
        local_center = local_log - local_log.mean(axis=0, keepdims=True)
        global_center = global_log - global_log.mean(axis=0, keepdims=True)
        denom = np.sqrt(np.sum(local_center**2, axis=0) * np.sum(global_center**2, axis=0))
        pcc = np.full(stop - start, np.nan, dtype=np.float64)
        valid = denom > 1.0e-12
        pcc[valid] = np.sum(local_center[:, valid] * global_center[:, valid], axis=0) / denom[valid]
        gene_correlations.append(pcc)
        local_count = local_probability * library_size[:, None]
        global_count = global_probability * library_size[:, None]
        local_count -= local_count.mean(axis=0, keepdims=True)
        global_count -= global_count.mean(axis=0, keepdims=True)
        count_denom = np.sqrt(np.sum(local_count**2, axis=0) * np.sum(global_count**2, axis=0))
        count_pcc = np.full(stop - start, np.nan, dtype=np.float64)
        count_valid = count_denom > 1.0e-12
        count_pcc[count_valid] = (
            np.sum(local_count[:, count_valid] * global_count[:, count_valid], axis=0)
            / count_denom[count_valid]
        )
        count_correlations.append(count_pcc)
    spot_cosine = spot_dot / np.maximum(np.sqrt(spot_local_sq * spot_global_sq), EPS)
    gene_pcc = np.concatenate(gene_correlations)
    count_pcc = np.concatenate(count_correlations)
    theta_global = theta @ mapping
    valid_gene = np.isfinite(gene_pcc)
    valid_count = np.isfinite(count_pcc)
    return {
        "relative_probability_error": float(relative),
        "median_spot_logcp10k_cosine": float(np.median(spot_cosine)),
        "q05_spot_logcp10k_cosine": float(np.quantile(spot_cosine, 0.05)),
        "median_gene_logcp10k_pcc": float(np.median(gene_pcc[valid_gene])),
        "fraction_gene_logcp10k_pcc_ge_0_95": float(np.mean(gene_pcc[valid_gene] >= 0.95)),
        "median_gene_count_pcc": float(np.median(count_pcc[valid_count])),
        "fraction_gene_count_pcc_ge_0_95": float(np.mean(count_pcc[valid_count] >= 0.95)),
        "local_n95": int(_n95(theta.mean(axis=0))),
        "global_n95": int(_n95(theta_global.mean(axis=0))),
        "local_effective_programs_mean": float(np.mean(_effective(theta))),
        "global_effective_programs_mean": float(np.mean(_effective(theta_global))),
        "mapping_primary_weight_mean": float(np.mean(np.max(mapping, axis=1))),
        "mapping_primary_weight_min": float(np.min(np.max(mapping, axis=1))),
    }


def evaluate_rows(
    bundle: Bundle,
    slide_rows: np.ndarray,
    mappings: np.ndarray,
    reference_by_batch: np.ndarray,
    *,
    block: int,
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for local_position, slide_row in enumerate(slide_rows):
        k = int(bundle.valid_mask[slide_row].sum())
        metrics = _evaluate_slide(
            bundle.theta[slide_row],
            bundle.basis[slide_row, :k],
            mappings[local_position, :k],
            reference_by_batch[int(bundle.batch_index[slide_row])],
            bundle.library_size[slide_row],
            bundle.gene_weight,
            block=int(block),
        )
        records.append(
            {
                **bundle.slides.iloc[int(slide_row)].to_dict(),
                **metrics,
            }
        )
    return pd.DataFrame(records)


def _summary(metrics: pd.DataFrame) -> dict[str, object]:
    columns = [
        "relative_probability_error",
        "median_spot_logcp10k_cosine",
        "q05_spot_logcp10k_cosine",
        "median_gene_logcp10k_pcc",
        "fraction_gene_logcp10k_pcc_ge_0_95",
        "median_gene_count_pcc",
        "fraction_gene_count_pcc_ge_0_95",
        "local_effective_programs_mean",
        "global_effective_programs_mean",
        "mapping_primary_weight_mean",
    ]
    result: dict[str, object] = {
        column: _quantiles(metrics[column].to_numpy(dtype=np.float64))
        for column in columns
    }
    result["local_n95_median"] = float(metrics["local_n95"].median())
    result["global_n95_median"] = float(metrics["global_n95"].median())
    return result


def _candidate_passes(summary: dict[str, object], args: argparse.Namespace) -> bool:
    relative = summary["relative_probability_error"]
    cosine = summary["median_spot_logcp10k_cosine"]
    gene = summary["median_gene_logcp10k_pcc"]
    return bool(
        relative["median"] <= float(args.gate_relative_median)
        and relative["q95"] <= float(args.gate_relative_q95)
        and cosine["median"] >= float(args.gate_spot_cosine)
        and gene["median"] >= float(args.gate_gene_pcc)
    )


def _save_final(
    bundle: Bundle,
    output_dir: Path,
    reference: np.ndarray,
    reference_by_batch: np.ndarray,
    technology_delta: np.ndarray,
    batch_delta: np.ndarray,
    mapping: np.ndarray,
    metrics: pd.DataFrame,
    history: list[dict[str, float]],
    candidate_table: pd.DataFrame,
    split: pd.DataFrame,
    selected_programs: int,
    selection_rule: str,
    args: argparse.Namespace,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=False)
    np.save(output_dir / "consensus_reference_probability_10k.npy", reference)
    np.save(
        output_dir / "consensus_reference_logcp10k_10k.npy",
        np.log1p(1.0e4 * reference).astype(np.float32),
    )
    np.save(output_dir / "consensus_reference_probability_by_batch_10k.npy", reference_by_batch)
    np.save(output_dir / "technology_gene_adapters.npy", technology_delta)
    np.save(output_dir / "batch_gene_adapters.npy", batch_delta)
    np.save(output_dir / "gene_order_10k.npy", bundle.genes)
    np.save(output_dir / "gene_weight_10k.npy", bundle.gene_weight)
    flat_mapping = mapping[bundle.valid_mask]
    np.save(output_dir / "local_to_consensus_mapping.npy", flat_mapping.astype(np.float32))
    total_spots = max(int(max(rows.max() for rows in bundle.global_rows)) + 1, sum(len(x) for x in bundle.theta))
    theta_consensus = np.zeros((total_spots, selected_programs), dtype=np.float32)
    theta_composition_consensus = np.zeros_like(theta_consensus)
    assignment_records: list[dict[str, object]] = []
    flat_offset = 0
    for slide in range(bundle.n_slides):
        k = int(bundle.valid_mask[slide].sum())
        local_mapping = mapping[slide, :k]
        rows = bundle.global_rows[slide]
        theta_consensus[rows] = bundle.theta[slide] @ local_mapping
        slide_dir = Path(str(bundle.slides.iloc[slide]["slide_dir"]))
        theta_composition = np.load(
            slide_dir / "theta_composition_mean.npy", allow_pickle=False
        ).astype(np.float32)
        theta_composition_consensus[rows] = theta_composition @ local_mapping
        for local_program in range(k):
            weights = local_mapping[local_program]
            order = np.argsort(-weights, kind="stable")
            assignment_records.append(
                {
                    **bundle.programs.iloc[flat_offset + local_program].to_dict(),
                    "consensus_primary_index": int(order[0]),
                    "consensus_primary": f"CP{int(order[0]):03d}",
                    "primary_weight": float(weights[order[0]]),
                    "secondary_index": int(order[1]),
                    "secondary_weight": float(weights[order[1]]),
                    "mapping_entropy": float(
                        -np.sum(weights * np.log(np.maximum(weights, EPS)))
                        / math.log(selected_programs)
                    ),
                }
            )
        flat_offset += k
    np.save(output_dir / "theta_rna_consensus.npy", theta_consensus)
    np.save(output_dir / "theta_composition_consensus.npy", theta_composition_consensus)
    assignment = pd.DataFrame(assignment_records)
    assignment.to_csv(output_dir / "local_program_assignment.csv", index=False)
    metrics.to_csv(output_dir / "reconstruction_metrics_by_slide.csv", index=False)
    candidate_table.to_csv(output_dir / "candidate_dictionary_selection.csv", index=False)
    split.to_csv(output_dir / "patient_holdout_split.csv", index=False)
    pd.DataFrame(history).to_csv(output_dir / "final_training_history.csv", index=False)
    batch_metadata = pd.DataFrame(
        {
            "batch_index": np.arange(len(bundle.batch_names)),
            "batch": bundle.batch_names,
            "technology_index": bundle.batch_to_technology,
            "technology": [
                bundle.technology_names[int(value)] for value in bundle.batch_to_technology
            ],
        }
    )
    batch_metadata.to_csv(output_dir / "batch_adapter_metadata.csv", index=False)
    program_mass = theta_consensus.astype(np.float64).sum(axis=0)
    program_mass /= np.maximum(program_mass.sum(), EPS)
    support = assignment.groupby("consensus_primary_index").agg(
        local_programs=("global_local_program_row", "count"),
        slides=("slide_index", "nunique"),
        patients=("patient", "nunique"),
        sources=("source", "nunique"),
        primary_weight_median=("primary_weight", "median"),
    )
    metadata_records = []
    for program in range(selected_programs):
        top = np.argsort(-reference[program], kind="stable")[:50]
        row = support.loc[program].to_dict() if program in support.index else {}
        metadata_records.append(
            {
                "consensus_program_index": program,
                "consensus_program": f"CP{program:03d}",
                "rna_mass_fraction": float(program_mass[program]),
                "top_genes": ";".join(bundle.genes[top].tolist()),
                **row,
            }
        )
    pd.DataFrame(metadata_records).to_csv(
        output_dir / "consensus_program_metadata.csv", index=False
    )
    summary = _summary(metrics)
    reference_path = output_dir / "consensus_reference_probability_10k.npy"
    theta_path = output_dir / "theta_rna_consensus.npy"
    contract = {
        "method": METHOD,
        "upstream_bayestme_modified": False,
        "slides": bundle.n_slides,
        "patients": int(bundle.slides["patient"].nunique()),
        "sources": bundle.slides["source"].value_counts().to_dict(),
        "technologies": bundle.slides["technology"].value_counts().to_dict(),
        "local_programs": int(bundle.valid_mask.sum()),
        "consensus_programs": int(selected_programs),
        "genes": bundle.n_genes,
        "objective": "P_local = theta_local @ B_local ~= theta_local @ A_slide @ W_batch",
        "neutral_reference": "one shared row-stochastic W",
        "batch_calibration": "technology plus centered source-by-technology logit adapters",
        "mapping": "nonnegative row-stochastic sparse soft A_slide",
        "candidate_selection_rule": selection_rule,
        "posterior_stability_used": bool(args.posterior_stability),
        "patient_holdout_fraction": float(args.holdout_fraction),
        "validation_gates": {
            "relative_error_median_max": float(args.gate_relative_median),
            "relative_error_q95_max": float(args.gate_relative_q95),
            "spot_cosine_median_min": float(args.gate_spot_cosine),
            "gene_pcc_median_min": float(args.gate_gene_pcc),
        },
        "full_fit_reconstruction_summary": summary,
        "reference_row_sum_maximum_error": float(
            np.max(np.abs(reference.sum(axis=1) - 1.0))
        ),
        "theta_row_sum_maximum_error": float(
            np.max(np.abs(theta_consensus.sum(axis=1) - 1.0))
        ),
        "mapping_row_sum_maximum_error": float(
            np.max(np.abs(flat_mapping.sum(axis=1) - 1.0))
        ),
        "reference_nonnegative": bool(np.all(reference >= 0)),
        "reference_sha256": _sha256(reference_path),
        "theta_sha256": _sha256(theta_path),
        "production_ready": bool(
            np.max(np.abs(reference.sum(axis=1) - 1.0)) <= 1.0e-5
            and np.max(np.abs(theta_consensus.sum(axis=1) - 1.0)) <= 1.0e-5
            and summary["median_spot_logcp10k_cosine"]["median"]
            >= float(args.gate_spot_cosine)
            and summary["median_gene_logcp10k_pcc"]["median"]
            >= float(args.gate_gene_pcc)
        ),
    }
    _atomic_json(output_dir / "consensus_contract.json", contract)
    return contract


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue-dir", type=Path, required=True)
    parser.add_argument("--master-spot-index", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--candidates", type=int, nargs="+", default=[64, 96, 128])
    parser.add_argument("--holdout-fraction", type=float, default=0.20)
    parser.add_argument("--candidate-epochs", type=int, default=140)
    parser.add_argument("--mapping-epochs", type=int, default=100)
    parser.add_argument("--final-epochs", type=int, default=240)
    parser.add_argument("--warmup-epochs", type=int, default=35)
    parser.add_argument("--learning-rate", type=float, default=0.025)
    parser.add_argument("--technology-bound", type=float, default=0.35)
    parser.add_argument("--batch-bound", type=float, default=0.20)
    parser.add_argument("--lambda-basis", type=float, default=0.20)
    parser.add_argument("--lambda-log-basis", type=float, default=0.0)
    parser.add_argument("--lambda-entropy", type=float, default=0.015)
    parser.add_argument("--lambda-collision", type=float, default=0.05)
    parser.add_argument("--lambda-diversity", type=float, default=0.02)
    parser.add_argument("--lambda-technology", type=float, default=0.05)
    parser.add_argument("--lambda-batch", type=float, default=0.10)
    parser.add_argument("--gate-relative-median", type=float, default=0.08)
    parser.add_argument("--gate-relative-q95", type=float, default=0.15)
    parser.add_argument("--gate-spot-cosine", type=float, default=0.99)
    parser.add_argument("--gate-gene-pcc", type=float, default=0.95)
    parser.add_argument("--evaluation-block", type=int, default=512)
    parser.add_argument("--seed", type=int, default=260802)
    parser.add_argument("--posterior-stability", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--max-slides", type=int)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.output_dir.exists():
        contract = args.output_dir / "consensus_contract.json"
        if contract.exists() and not args.dry_run:
            print(contract.read_text(encoding="utf-8"), flush=True)
            return
        if not args.dry_run:
            raise FileExistsError(f"refusing existing output directory: {args.output_dir}")
    random.seed(int(args.seed))
    np.random.seed(int(args.seed))
    torch.manual_seed(int(args.seed))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    bundle = load_bundle(
        args.queue_dir,
        args.master_spot_index,
        use_posterior_stability=bool(args.posterior_stability),
        max_slides=args.max_slides,
    )
    print(
        json.dumps(
            {
                "device": str(device),
                "slides": bundle.n_slides,
                "patients": int(bundle.slides["patient"].nunique()),
                "local_programs": int(bundle.valid_mask.sum()),
                "genes": bundle.n_genes,
                "batches": bundle.batch_names,
                "technologies": bundle.technology_names,
                "posterior_stability_quantiles": _quantiles(
                    bundle.stability[bundle.valid_mask]
                ),
            },
            ensure_ascii=False,
            indent=2,
        ),
        flush=True,
    )
    if args.dry_run:
        rows = np.arange(min(bundle.n_slides, 8), dtype=np.int64)
        programs = min(8, int(bundle.valid_mask[rows].sum()))
        reference, mapping = _initial_reference(bundle, rows, programs, int(args.seed))
        model, _ = _fit_model(
            bundle,
            rows,
            reference,
            mapping,
            device=device,
            epochs=2,
            warmup_epochs=1,
            learning_rate=float(args.learning_rate),
            seed=int(args.seed),
            args=args,
        )
        neutral, by_batch, _, _, fitted_mapping = _model_arrays(model)
        print(
            json.dumps(
                {
                    "dry_run": "passed",
                    "reference_shape": list(neutral.shape),
                    "batch_reference_shape": list(by_batch.shape),
                    "mapping_shape": list(fitted_mapping.shape),
                },
                indent=2,
            ),
            flush=True,
        )
        return
    train_rows, validation_rows = _patient_holdout(
        bundle.slides, float(args.holdout_fraction), int(args.seed)
    )
    split = bundle.slides.copy()
    split["split"] = "train"
    split.loc[validation_rows, "split"] = "validation"
    candidate_records: list[dict[str, object]] = []
    fitted_candidates: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
    for programs in sorted(set(int(value) for value in args.candidates)):
        print(json.dumps({"candidate": programs, "stage": "initialise"}), flush=True)
        reference_init, mapping_init = _initial_reference(
            bundle, train_rows, programs, int(args.seed) + programs
        )
        model, _ = _fit_model(
            bundle,
            train_rows,
            reference_init,
            mapping_init,
            device=device,
            epochs=int(args.candidate_epochs),
            warmup_epochs=min(int(args.warmup_epochs), int(args.candidate_epochs) // 2),
            learning_rate=float(args.learning_rate),
            seed=int(args.seed) + programs,
            args=args,
        )
        neutral, by_batch, technology_delta, batch_delta, _ = _model_arrays(model)
        validation_mapping = _fit_mapping_only(
            bundle,
            validation_rows,
            neutral,
            by_batch,
            device=device,
            epochs=int(args.mapping_epochs),
            learning_rate=float(args.learning_rate),
            seed=int(args.seed) + 10_000 + programs,
            args=args,
        )
        validation_metrics = evaluate_rows(
            bundle,
            validation_rows,
            validation_mapping,
            by_batch,
            block=int(args.evaluation_block),
        )
        summary = _summary(validation_metrics)
        passed = _candidate_passes(summary, args)
        record = {
            "consensus_programs": programs,
            "passed": passed,
            "relative_error_median": summary["relative_probability_error"]["median"],
            "relative_error_q95": summary["relative_probability_error"]["q95"],
            "spot_cosine_median": summary["median_spot_logcp10k_cosine"]["median"],
            "gene_pcc_median": summary["median_gene_logcp10k_pcc"]["median"],
            "count_pcc_median": summary["median_gene_count_pcc"]["median"],
            "global_n95_median": summary["global_n95_median"],
        }
        candidate_records.append(record)
        fitted_candidates[programs] = (neutral, by_batch, technology_delta, batch_delta)
        print(json.dumps({"candidate_result": record}, ensure_ascii=False), flush=True)
        del model
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
    candidate_table = pd.DataFrame(candidate_records).sort_values("consensus_programs")
    passing = candidate_table.loc[candidate_table["passed"]]
    if len(passing):
        selected_programs = int(passing.iloc[0]["consensus_programs"])
        selection_rule = "smallest_candidate_meeting_patient_holdout_reconstruction_gates"
    else:
        ranked = candidate_table.sort_values(
            ["gene_pcc_median", "spot_cosine_median", "relative_error_median"],
            ascending=[False, False, True],
        )
        selected_programs = int(ranked.iloc[0]["consensus_programs"])
        selection_rule = "best_patient_holdout_candidate_no_candidate_met_all_gates"
    selected_neutral, _, selected_technology, selected_batch = fitted_candidates[selected_programs]
    full_mapping_init = _mapping_logits(bundle.basis, bundle.valid_mask, selected_neutral)
    final_model = RCCD(
        selected_neutral,
        full_mapping_init,
        len(bundle.technology_names),
        len(bundle.batch_names),
        bundle.batch_to_technology,
        bundle.technology_prior,
        bundle.batch_prior_within_technology,
        technology_bound=float(args.technology_bound),
        batch_bound=float(args.batch_bound),
    )
    with torch.no_grad():
        final_model.technology_raw.copy_(
            torch.from_numpy(np.arctanh(np.clip(selected_technology / max(float(args.technology_bound), EPS), -0.95, 0.95))).float()
        )
        final_model.batch_raw.copy_(
            torch.from_numpy(np.arctanh(np.clip(selected_batch / max(float(args.batch_bound), EPS), -0.95, 0.95))).float()
        )
    # Reuse the common fitter after transferring the selected initial state.
    final_model = final_model.to(device)
    final_data = _tensor_bundle(bundle, np.arange(bundle.n_slides), device)
    optimizer = torch.optim.Adam(final_model.parameters(), lr=float(args.learning_rate))
    final_history: list[dict[str, float]] = []
    for epoch in range(int(args.final_epochs)):
        optimizer.zero_grad(set_to_none=True)
        temperature = max(
            0.35,
            1.0 - 0.65 * epoch / max(int(args.final_epochs) - 1, 1),
        )
        mapping = final_model.mapping(temperature)
        neutral = final_model.neutral_reference()
        references = final_model.reference_by_batch()
        technology_delta, batch_delta = final_model.deltas()
        total, pieces = _reconstruction_loss(
            mapping,
            references,
            neutral,
            final_data,
            technology_delta,
            batch_delta,
            lambda_basis=float(args.lambda_basis),
            lambda_log_basis=float(args.lambda_log_basis),
            lambda_entropy=float(args.lambda_entropy),
            lambda_collision=float(args.lambda_collision),
            lambda_diversity=float(args.lambda_diversity),
            lambda_technology=float(args.lambda_technology),
            lambda_batch=float(args.lambda_batch),
        )
        total.backward()
        torch.nn.utils.clip_grad_norm_(final_model.parameters(), 5.0)
        optimizer.step()
        record = {"epoch": float(epoch + 1), "total": float(total.detach().cpu())}
        record.update({name: float(value.detach().cpu()) for name, value in pieces.items()})
        final_history.append(record)
        if epoch == 0 or (epoch + 1) % 25 == 0 or epoch + 1 == int(args.final_epochs):
            print(json.dumps({"final_fit": record}), flush=True)
    reference, reference_by_batch, technology_delta, batch_delta, mapping = _model_arrays(final_model)
    metrics = evaluate_rows(
        bundle,
        np.arange(bundle.n_slides, dtype=np.int64),
        mapping,
        reference_by_batch,
        block=int(args.evaluation_block),
    )
    contract = _save_final(
        bundle,
        args.output_dir,
        reference,
        reference_by_batch,
        technology_delta,
        batch_delta,
        mapping,
        metrics,
        final_history,
        candidate_table,
        split,
        selected_programs,
        selection_rule,
        args,
    )
    print(json.dumps(contract, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
