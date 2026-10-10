"""Train an image-only mapper for the G384 pan-atlas consensus programs.

The shared Virchow2/morphology graph trunk has two compositional heads:

* cell-composition proportions -> ``theta_composition_consensus``;
* RNA-mixture proportions -> ``theta_rna_consensus``.

The frozen expression decoder is the slide-matched batch-adapted pan-atlas reference
``W``.  Crucially, expression PCC is never computed against ``theta @ W`` from
the consensus fit itself.  Its reference is the original per-slide 10K
BayesTME inverse ``theta_local @ basis_local``.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from openst_final.he_uni_mamba_bnn_inr import normalize_slide_coordinates
from run_virchow2_bayestme_global_mapper import _morphology_graph
from run_virchow2_fseg_k64_mapper import BoundedINR, FSegGATv2, _to_device_graph


EPS = 1.0e-8


class _GradientReverse(torch.autograd.Function):
    @staticmethod
    def forward(ctx, value: torch.Tensor, strength: float) -> torch.Tensor:
        ctx.strength = float(strength)
        return value.view_as(value)

    @staticmethod
    def backward(ctx, gradient: torch.Tensor) -> tuple[torch.Tensor, None]:
        return -ctx.strength * gradient, None


@dataclass
class PanAtlasSlide:
    key: str
    sample_id: str
    source: str
    patient: str
    split: str
    technology: str
    source_index: int
    technology_index: int
    reference_rows: np.ndarray
    local_rows: np.ndarray
    features: np.ndarray
    coordinates: np.ndarray
    theta_rna: np.ndarray
    theta_composition: np.ndarray
    local_theta_rna: np.ndarray
    local_basis_panel: np.ndarray
    consensus_basis_panel: np.ndarray
    consensus_reference: np.ndarray
    slide_dir: Path
    spot_index: pd.DataFrame
    graphs: tuple[
        tuple[np.ndarray, np.ndarray, np.ndarray],
        tuple[np.ndarray, np.ndarray, np.ndarray],
    ]


class PanAtlasDualMapper(nn.Module):
    def __init__(
        self,
        input_dim: int,
        *,
        hidden: int,
        programs: int,
        program_groups: np.ndarray,
        graph_layers: int,
        context_graph_layers: int,
        graph_heads: int,
        dropout: float,
        type_lift_alpha: float,
        sources: int,
        technologies: int,
        support_mass: float,
        hierarchical_logits: bool = False,
    ) -> None:
        super().__init__()
        self.programs = int(programs)
        self.type_lift_alpha = float(type_lift_alpha)
        self.support_mass = float(support_mass)
        self.hierarchical_logits = bool(hierarchical_logits)
        self.encoder = nn.Sequential(
            nn.LayerNorm(input_dim),
            nn.Linear(input_dim, hidden * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden * 2, hidden),
            nn.LayerNorm(hidden),
        )
        self.centered_encoder = nn.Sequential(
            nn.LayerNorm(input_dim),
            nn.Linear(input_dim, hidden * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden * 2, hidden),
            nn.LayerNorm(hidden),
        )
        self.centered_feature_gate = nn.Parameter(torch.tensor(-1.0))
        group_index = torch.as_tensor(program_groups, dtype=torch.long)
        if group_index.numel() != programs:
            raise ValueError("program-group mapping does not match programs")
        self.register_buffer("program_group_index", group_index)
        self.program_groups = int(group_index.max().item()) + 1
        self.local_graph = nn.ModuleList(
            [FSegGATv2(hidden, graph_heads, dropout) for _ in range(graph_layers)]
        )
        self.context_graph = nn.ModuleList(
            [
                FSegGATv2(hidden, graph_heads, dropout)
                for _ in range(context_graph_layers)
            ]
        )
        self.local_graph_gate = nn.Parameter(torch.tensor(-1.5))
        self.context_graph_gate = nn.Parameter(torch.tensor(-2.0))
        self.graph_fusion = nn.Sequential(
            nn.Linear(hidden * 3, hidden * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden * 2, hidden),
            nn.LayerNorm(hidden),
        )
        self.edge_predictor = nn.Sequential(
            nn.Linear(hidden * 2 + 1, hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, 1),
        )
        self.v4_graph_fusion = nn.Sequential(
            nn.Linear(hidden * 4, hidden * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden * 2, hidden),
            nn.LayerNorm(hidden),
        )
        self.slide_encoder = nn.Sequential(
            nn.Linear(hidden * 2, hidden),
            nn.GELU(),
            nn.LayerNorm(hidden),
        )
        fused = hidden * 3
        self.shared_local = nn.Sequential(
            nn.Linear(fused, hidden * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden * 2, hidden),
            nn.GELU(),
        )
        self.rna_local = nn.Linear(hidden, programs)
        self.composition_local = nn.Linear(hidden, programs)
        self.rna_slide = nn.Linear(hidden, programs)
        self.composition_slide = nn.Linear(hidden, programs)
        self.support_local = nn.Linear(hidden, programs)
        self.support_slide = nn.Linear(hidden, programs)
        self.support_metric = nn.Linear(hidden, hidden, bias=False)
        self.support_prototype = nn.Parameter(torch.empty(programs, hidden))
        nn.init.normal_(self.support_prototype, std=hidden**-0.5)
        self.support_type_gate = nn.Parameter(torch.tensor(-0.7))
        self.log_support_temperature = nn.Parameter(torch.tensor(math.log(0.20)))
        self.source_classifier = nn.Sequential(
            nn.Linear(hidden, hidden // 2),
            nn.GELU(),
            nn.Linear(hidden // 2, max(int(sources), 1)),
        )
        self.technology_classifier = nn.Sequential(
            nn.Linear(hidden, hidden // 2),
            nn.GELU(),
            nn.Linear(hidden // 2, max(int(technologies), 1)),
        )
        self.rna_group_local = nn.Linear(hidden, self.program_groups)
        self.composition_group_local = nn.Linear(hidden, self.program_groups)
        self.rna_group_slide = nn.Linear(hidden, self.program_groups)
        self.composition_group_slide = nn.Linear(hidden, self.program_groups)

        # A metric head learns which program identity is morphologically present;
        # the ordinary linear heads remain responsible for calibrated proportions.
        self.rna_metric = nn.Linear(hidden, hidden, bias=False)
        self.composition_metric = nn.Linear(hidden, hidden, bias=False)
        self.rna_prototype = nn.Parameter(torch.empty(programs, hidden))
        self.composition_prototype = nn.Parameter(torch.empty(programs, hidden))
        nn.init.normal_(self.rna_prototype, std=hidden**-0.5)
        nn.init.normal_(self.composition_prototype, std=hidden**-0.5)
        self.log_type_temperature = nn.Parameter(torch.tensor(math.log(0.20)))
        self.rna_type_gate = nn.Parameter(torch.tensor(-1.4))
        self.composition_type_gate = nn.Parameter(torch.tensor(-1.4))
        self.rna_group_gate = nn.Parameter(torch.tensor(0.0))
        self.composition_group_gate = nn.Parameter(torch.tensor(0.0))

        # Each program learns its own amount of local and niche-scale smoothing.
        # Boundary programs can keep these gates near zero, while broad tissue
        # programs can use more graph context.
        self.rna_local_refine = nn.Parameter(torch.full((programs,), -2.2))
        self.rna_context_refine = nn.Parameter(torch.full((programs,), -2.6))
        self.composition_local_refine = nn.Parameter(torch.full((programs,), -2.2))
        self.composition_context_refine = nn.Parameter(torch.full((programs,), -2.6))

    @staticmethod
    def _center(logits: torch.Tensor) -> torch.Tensor:
        return logits - logits.mean(dim=1, keepdim=True)

    def forward(
        self,
        features: torch.Tensor,
        graphs: tuple[
            tuple[torch.Tensor, torch.Tensor, torch.Tensor],
            tuple[torch.Tensor, torch.Tensor, torch.Tensor],
        ],
        *,
        support_gamma: float = 1.0,
        domain_strength: float = 0.0,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, dict[str, torch.Tensor]]:
        local_graph, context_graph = graphs
        raw = self.encoder(features)
        per_spot = F.layer_norm(features.float(), (features.shape[1],))
        centered_features = per_spot - per_spot.mean(dim=0, keepdim=True)
        centered = self.centered_encoder(centered_features)
        encoded = raw + torch.sigmoid(self.centered_feature_gate) * centered

        edge_source, edge_target, edge_prior = local_graph
        edge_input = torch.cat(
            [
                torch.abs(encoded[edge_source] - encoded[edge_target]),
                encoded[edge_source] * encoded[edge_target],
                torch.log(edge_prior.clamp_min(1.0e-6))[:, None],
            ],
            dim=1,
        )
        edge_probability = torch.sigmoid(self.edge_predictor(edge_input).squeeze(1))
        edge_probability = torch.where(
            edge_source == edge_target,
            torch.ones_like(edge_probability),
            edge_probability,
        )
        local_context = encoded
        if self.hierarchical_logits:
            for block in self.local_graph:
                local_context = block(local_context, local_graph)
            local_context = encoded + torch.sigmoid(self.local_graph_gate) * (
                local_context - encoded
            )
        homogeneous = self._edge_aggregate(
            local_context[edge_target],
            edge_source,
            edge_prior * edge_probability,
            len(encoded),
        )
        boundary = self._edge_aggregate(
            local_context[edge_source] - local_context[edge_target],
            edge_source,
            edge_prior * (1.0 - edge_probability),
            len(encoded),
        )
        niche_context = encoded
        for block in self.context_graph:
            niche_context = block(niche_context, context_graph)
        niche_context = encoded + torch.sigmoid(self.context_graph_gate) * (
            niche_context - encoded
        )
        spatial = self.v4_graph_fusion(
            torch.cat([encoded, homogeneous, boundary, niche_context], dim=1)
        )
        pooled = torch.cat(
            [spatial.mean(dim=0), spatial.std(dim=0, unbiased=False)], dim=0
        )
        slide = self.slide_encoder(pooled)
        repeated = slide[None, :].expand(len(spatial), -1)
        shared = self.shared_local(
            torch.cat([spatial, repeated, spatial * repeated], dim=1)
        )
        support_temperature = self.log_support_temperature.exp().clamp(0.05, 1.0)
        support_embedding = F.normalize(self.support_metric(shared), dim=1)
        support_prototype_score = (
            support_embedding @ F.normalize(self.support_prototype, dim=1).T
        ) / support_temperature
        support_logits = (
            self.support_local(shared)
            + self.support_slide(slide)[None, :]
            + torch.sigmoid(self.support_type_gate) * support_prototype_score
        )
        rna_group = self.rna_group_local(shared) + self.rna_group_slide(slide)[None, :]
        composition_group = (
            self.composition_group_local(shared)
            + self.composition_group_slide(slide)[None, :]
        )
        temperature = self.log_type_temperature.exp().clamp(0.05, 1.0)
        rna_embedding = F.normalize(self.rna_metric(shared), dim=1)
        composition_embedding = F.normalize(self.composition_metric(shared), dim=1)
        rna_type = rna_embedding @ F.normalize(self.rna_prototype, dim=1).T / temperature
        composition_type = (
            composition_embedding
            @ F.normalize(self.composition_prototype, dim=1).T
            / temperature
        )
        support_log_probability = F.logsigmoid(support_logits)
        rna = (
            self.rna_local(shared)
            + self.rna_slide(slide)[None, :]
            + float(support_gamma) * support_log_probability
        )
        composition = (
            self.composition_local(shared)
            + self.composition_slide(slide)[None, :]
            + float(support_gamma) * support_log_probability
        )
        if self.hierarchical_logits:
            rna_group_lift = rna_group[:, self.program_group_index]
            composition_group_lift = composition_group[:, self.program_group_index]
            rna = (
                rna
                + torch.sigmoid(self.rna_group_gate) * self._center(rna_group_lift)
                + torch.sigmoid(self.rna_type_gate) * self._center(rna_type)
            )
            composition = (
                composition
                + torch.sigmoid(self.composition_group_gate)
                * self._center(composition_group_lift)
                + torch.sigmoid(self.composition_type_gate)
                * self._center(composition_type)
            )
        reversed_slide = _GradientReverse.apply(slide, float(domain_strength))
        auxiliary = {
            "rna_group_logits": rna_group,
            "composition_group_logits": composition_group,
            "rna_type_logits": rna_type,
            "composition_type_logits": composition_type,
            "support_logits": support_logits,
            "edge_source": edge_source,
            "edge_target": edge_target,
            "edge_probability": edge_probability,
            "source_logits": self.source_classifier(reversed_slide[None, :]),
            "technology_logits": self.technology_classifier(reversed_slide[None, :]),
        }
        return self._center(rna), self._center(composition), spatial, auxiliary

    @staticmethod
    def _edge_aggregate(
        message: torch.Tensor,
        source: torch.Tensor,
        weight: torch.Tensor,
        spots: int,
    ) -> torch.Tensor:
        aggregate = message.new_zeros((int(spots), message.shape[1]))
        aggregate.index_add_(0, source, message * weight[:, None].to(message.dtype))
        denominator = weight.new_zeros(int(spots))
        denominator.index_add_(0, source, weight)
        return aggregate / denominator[:, None].clamp_min(EPS)

    @staticmethod
    def _neighbor_delta(
        logits: torch.Tensor,
        graph: tuple[torch.Tensor, torch.Tensor, torch.Tensor],
    ) -> torch.Tensor:
        source, target, weight = graph
        weighted = weight.to(logits.dtype)[:, None] * logits[target]
        aggregate = torch.zeros_like(logits)
        aggregate.index_add_(0, source, weighted)
        denominator = torch.zeros(
            logits.shape[0], dtype=logits.dtype, device=logits.device
        )
        denominator.index_add_(0, source, weight.to(logits.dtype))
        return aggregate / denominator[:, None].clamp_min(EPS) - logits

    @classmethod
    def _program_refine(
        cls,
        logits: torch.Tensor,
        local_graph: tuple[torch.Tensor, torch.Tensor, torch.Tensor],
        context_graph: tuple[torch.Tensor, torch.Tensor, torch.Tensor],
        local_gate: torch.Tensor,
        context_gate: torch.Tensor,
    ) -> torch.Tensor:
        return (
            logits
            + torch.sigmoid(local_gate)[None, :] * cls._neighbor_delta(logits, local_graph)
            + torch.sigmoid(context_gate)[None, :]
            * cls._neighbor_delta(logits, context_graph)
        )


def _canonical(sample_id: str) -> str:
    return str(sample_id).split(":", 1)[-1]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _resolve_feature_path(value: str, feature_root: Path | None) -> Path:
    path = Path(str(value))
    if path.is_file():
        return path
    if feature_root is not None:
        candidate = feature_root / path
        if candidate.is_file():
            return candidate
        candidate = feature_root / path.name
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(path)


def _normalise_rows(value: np.ndarray) -> np.ndarray:
    value = np.maximum(np.asarray(value, dtype=np.float32), 0.0)
    return value / np.maximum(value.sum(axis=1, keepdims=True), EPS)


def _select_loss_panel(
    consensus_dir: Path, *, count: int
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    genes = np.load(
        consensus_dir / "gene_order_10k.npy", allow_pickle=False
    ).astype(str)
    weights = np.load(
        consensus_dir / "gene_weight_10k.npy", allow_pickle=False
    ).astype(np.float64)
    if len(genes) != len(weights):
        raise RuntimeError("gene order and gene weights disagree")
    count = min(int(count), len(genes))
    order = np.argsort(-weights, kind="stable")[:count].astype(np.int64)
    selected_weight = np.sqrt(np.maximum(weights[order], EPS))
    selected_weight /= max(float(selected_weight.mean()), EPS)
    table = pd.DataFrame(
        {
            "loss_panel_rank": np.arange(1, len(order) + 1),
            "gene_index": order,
            "gene": genes[order],
            "consensus_gene_weight": weights[order],
            "loss_weight": selected_weight,
        }
    )
    return order, selected_weight.astype(np.float32), table


def _program_groups(reference: np.ndarray, *, count: int) -> np.ndarray:
    """Group fine programs by their full-10K expression geometry.

    The hierarchy supplies an easier coarse program-type task without changing
    the 384-row decoder or merging any final proportions.
    """
    reference = _normalise_rows(np.asarray(reference, dtype=np.float32))
    root = np.sqrt(np.maximum(reference, 0.0))
    similarity = np.clip(root @ root.T, 0.0, 1.0)
    distance = np.clip(1.0 - similarity, 0.0, 1.0)
    np.fill_diagonal(distance, 0.0)
    tree = linkage(squareform(distance, checks=False), method="average")
    labels = fcluster(tree, t=min(int(count), len(reference)), criterion="maxclust")
    # Canonicalise labels by first program index so repeated runs are stable.
    ordered = sorted(np.unique(labels), key=lambda value: int(np.flatnonzero(labels == value)[0]))
    remap = {int(value): index for index, value in enumerate(ordered)}
    return np.asarray([remap[int(value)] for value in labels], dtype=np.int64)


def _morphology_program_groups(
    slides: list[PanAtlasSlide],
    reference: np.ndarray,
    *,
    count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build train-patient-only program groups from residual HE directions.

    Each slide contributes equally within a patient, each patient contributes
    equally within a source, and each source contributes equally overall.  A
    small expression-geometry term stabilises programs with weak HE residuals
    without allowing the decoder alone to define the morphology hierarchy.
    """
    train = [slide for slide in slides if slide.split == "train"]
    if not train:
        raise RuntimeError("morphology grouping requires training slides")
    programs = int(reference.shape[0])
    patient_sum: dict[tuple[str, str], np.ndarray] = {}
    patient_count: dict[tuple[str, str], int] = {}
    compressed_dim: int | None = None
    for slide in train:
        features = np.asarray(slide.features, dtype=np.float32)
        features = features - features.mean(axis=1, keepdims=True)
        features /= np.sqrt(
            np.maximum(np.mean(np.square(features), axis=1, keepdims=True), EPS)
        )
        features -= features.mean(axis=0, keepdims=True)
        block = 8 if features.shape[1] % 8 == 0 else 1
        if block > 1:
            features = features.reshape(
                len(features), features.shape[1] // block, block
            ).mean(axis=2)
        if compressed_dim is None:
            compressed_dim = int(features.shape[1])
        elif int(features.shape[1]) != compressed_dim:
            raise RuntimeError("Virchow2 feature dimensions disagree across slides")

        theta = 0.5 * (
            np.asarray(slide.theta_rna, dtype=np.float32)
            + np.asarray(slide.theta_composition, dtype=np.float32)
        )
        prevalence = np.maximum(theta.mean(axis=0, keepdims=True), 1.0e-6)
        residual = theta / np.sqrt(prevalence)
        residual -= residual.mean(axis=0, keepdims=True)
        denominator = np.sqrt(np.maximum(np.square(residual).sum(axis=0), EPS))
        prototype = (residual.T @ features) / denominator[:, None]
        prototype /= np.sqrt(
            np.maximum(np.square(prototype).sum(axis=1, keepdims=True), EPS)
        )
        key = (slide.source, slide.patient)
        if key not in patient_sum:
            patient_sum[key] = np.zeros_like(prototype, dtype=np.float64)
            patient_count[key] = 0
        patient_sum[key] += prototype
        patient_count[key] += 1

    source_blocks: list[np.ndarray] = []
    for source in sorted({key[0] for key in patient_sum}):
        patient_blocks = [
            patient_sum[key] / max(patient_count[key], 1)
            for key in patient_sum
            if key[0] == source
        ]
        source_blocks.append(np.mean(patient_blocks, axis=0))
    morphology = np.mean(source_blocks, axis=0).astype(np.float32)
    morphology /= np.sqrt(
        np.maximum(np.square(morphology).sum(axis=1, keepdims=True), EPS)
    )
    morphology_similarity = np.clip(
        0.5 * (morphology @ morphology.T + 1.0), 0.0, 1.0
    )

    expression = _normalise_rows(np.asarray(reference, dtype=np.float32))
    expression_root = np.sqrt(np.maximum(expression, 0.0))
    expression_similarity = np.clip(expression_root @ expression_root.T, 0.0, 1.0)
    similarity = np.clip(
        0.75 * morphology_similarity + 0.25 * expression_similarity,
        0.0,
        1.0,
    )
    distance = np.clip(1.0 - similarity, 0.0, 1.0)
    np.fill_diagonal(distance, 0.0)
    tree = linkage(squareform(distance, checks=False), method="average")
    labels = fcluster(
        tree, t=min(int(count), programs), criterion="maxclust"
    )
    ordered = sorted(
        np.unique(labels), key=lambda value: int(np.flatnonzero(labels == value)[0])
    )
    remap = {int(value): index for index, value in enumerate(ordered)}
    groups = np.asarray([remap[int(value)] for value in labels], dtype=np.int64)
    return groups, morphology, similarity.astype(np.float32)


def _load_slides(
    *,
    consensus_dir: Path,
    per_slide_dir: Path,
    master_spot_index: Path,
    split_manifest: Path,
    panel: np.ndarray,
    feature_root: Path | None,
    neighbors: int,
    context_neighbors: int,
    sources: tuple[str, ...],
    require_all_splits: bool = True,
) -> tuple[list[PanAtlasSlide], np.ndarray, np.ndarray, pd.DataFrame]:
    theta_rna_all = np.load(
        consensus_dir / "theta_rna_consensus.npy", mmap_mode="r"
    )
    theta_composition_all = np.load(
        consensus_dir / "theta_composition_consensus.npy", mmap_mode="r"
    )
    reference = np.load(
        consensus_dir / "consensus_reference_probability_10k.npy",
        mmap_mode="r",
    )
    batch_reference_path = (
        consensus_dir / "consensus_reference_probability_by_batch_10k.npy"
    )
    batch_metadata_path = consensus_dir / "batch_adapter_metadata.csv"
    if batch_reference_path.exists() and batch_metadata_path.exists():
        reference_by_batch = np.load(batch_reference_path, mmap_mode="r")
        batch_metadata = pd.read_csv(batch_metadata_path)
        batch_lookup = {
            str(row.batch): int(row.batch_index)
            for row in batch_metadata.itertuples(index=False)
        }
        if reference_by_batch.shape[1:] != reference.shape:
            raise RuntimeError("batch-adapted and neutral consensus reference shapes differ")
        batch_panel = {
            key: np.asarray(reference_by_batch[index][:, panel], dtype=np.float32)
            for key, index in batch_lookup.items()
        }
        batch_full = {
            key: np.asarray(reference_by_batch[index], dtype=np.float32)
            for key, index in batch_lookup.items()
        }
    else:
        batch_lookup = {}
        batch_panel = {}
        batch_full = {}
    neutral_panel = np.asarray(reference[:, panel], dtype=np.float32)
    neutral_full = np.asarray(reference, dtype=np.float32)
    genes = np.load(
        consensus_dir / "gene_order_10k.npy", allow_pickle=False
    ).astype(str)
    if theta_rna_all.shape != theta_composition_all.shape:
        raise RuntimeError("RNA and composition consensus theta shapes differ")
    if theta_rna_all.shape[1] != reference.shape[0]:
        raise RuntimeError("consensus theta/reference dimensions differ")
    if reference.shape[1] != len(genes):
        raise RuntimeError("consensus reference/gene dimensions differ")

    master = pd.read_csv(master_spot_index)
    required_master = {
        "global_index",
        "slide_index",
        "local_index",
        "sample_id",
        "source",
        "patient",
        "technology",
        "barcode",
    }
    missing = required_master.difference(master.columns)
    if missing:
        raise RuntimeError(f"master spot index lacks columns: {sorted(missing)}")
    split = pd.read_csv(split_manifest)
    source_lookup = {value: index for index, value in enumerate(sorted(sources))}
    technology_lookup = {
        value: index
        for index, value in enumerate(sorted(master["technology"].astype(str).unique()))
    }
    assignments = pd.read_csv(consensus_dir / "local_program_assignment.csv")
    slide_directory = (
        assignments[["sample_id", "slide_dir"]]
        .drop_duplicates()
        .set_index("sample_id")["slide_dir"]
        .astype(str)
        .to_dict()
    )
    slides: list[PanAtlasSlide] = []
    coverage_records: list[dict[str, object]] = []
    for row in split.itertuples(index=False):
        source = str(row.source)
        if source not in sources:
            continue
        sample_id = str(row.reference_sample_id)
        block = master.loc[
            master["source"].astype(str).eq(source)
            & master["sample_id"].astype(str).eq(sample_id)
        ].sort_values("local_index", kind="stable")
        if not len(block):
            raise RuntimeError(f"master spots absent for {sample_id}")
        feature_path = _resolve_feature_path(str(row.feature_path), feature_root)
        with np.load(feature_path, allow_pickle=False) as archive:
            features = np.asarray(archive["features"], dtype=np.float16)
            barcodes = np.asarray(archive["barcode"]).astype(str)
            coordinates = np.asarray(archive["coords"], dtype=np.float32)
        lookup = {barcode: index for index, barcode in enumerate(barcodes)}
        matched = block["barcode"].astype(str).isin(lookup)
        coverage = float(matched.mean())
        if coverage < 0.90:
            raise RuntimeError(
                f"Virchow2/spot coverage below 90% for {sample_id}: {coverage:.3%}"
            )
        block = block.loc[matched].copy()
        feature_rows = np.asarray(
            [lookup[value] for value in block["barcode"].astype(str)],
            dtype=np.int64,
        )
        global_rows = block["global_index"].to_numpy(dtype=np.int64)
        local_rows = block["local_index"].to_numpy(dtype=np.int64)
        raw_slide_dir = Path(slide_directory[sample_id])
        slide_dir = raw_slide_dir
        if not slide_dir.is_dir():
            slide_dir = per_slide_dir / "slides" / raw_slide_dir.name
        if not slide_dir.is_dir():
            raise FileNotFoundError(raw_slide_dir)
        local_theta = np.load(
            slide_dir / "theta_rna_mean.npy", allow_pickle=False
        ).astype(np.float32)
        local_basis = np.load(
            slide_dir / "basis_probability_mean.npy", mmap_mode="r"
        )
        if local_theta.shape[1] != local_basis.shape[0]:
            raise RuntimeError(f"local theta/basis mismatch for {sample_id}")
        if int(local_rows.max()) >= len(local_theta):
            raise RuntimeError(f"local spot rows exceed theta for {sample_id}")
        local_feature = features[feature_rows]
        local_coordinates = coordinates[feature_rows]
        technology = str(block["technology"].iloc[0])
        consensus_batch = f"{source}|{technology}"
        consensus_basis_panel = batch_panel.get(consensus_batch, neutral_panel)
        consensus_reference = batch_full.get(consensus_batch, neutral_full)
        theta_rna = _normalise_rows(theta_rna_all[global_rows])
        theta_composition = _normalise_rows(theta_composition_all[global_rows])
        local_graph = _morphology_graph(
            local_coordinates,
            np.asarray(local_feature, dtype=np.float32),
            neighbors=int(neighbors),
        )
        context_graph = _morphology_graph(
            local_coordinates,
            np.asarray(local_feature, dtype=np.float32),
            neighbors=int(context_neighbors),
        )
        slides.append(
            PanAtlasSlide(
                key=f"{source}:{_canonical(sample_id)}",
                sample_id=sample_id,
                source=source,
                patient=str(row.patient),
                split=str(row.split),
                technology=technology,
                source_index=int(source_lookup[source]),
                technology_index=int(
                    technology_lookup[str(block["technology"].iloc[0])]
                ),
                reference_rows=global_rows,
                local_rows=local_rows,
                features=local_feature,
                coordinates=local_coordinates,
                theta_rna=theta_rna,
                theta_composition=theta_composition,
                local_theta_rna=local_theta[local_rows],
                local_basis_panel=np.asarray(
                    local_basis[:, panel], dtype=np.float32
                ),
                consensus_basis_panel=consensus_basis_panel,
                consensus_reference=consensus_reference,
                slide_dir=slide_dir,
                spot_index=block.reset_index(drop=True),
                graphs=(local_graph, context_graph),
            )
        )
        coverage_records.append(
            {
                "sample_id": sample_id,
                "source": source,
                "patient": str(row.patient),
                "split": str(row.split),
                "reference_spots": int(len(matched)),
                "matched_spots": int(matched.sum()),
                "coverage": coverage,
                "feature_path": str(feature_path),
                "slide_dir": str(slide_dir),
                "consensus_batch": consensus_batch,
                "batch_adapted_decoder": bool(consensus_batch in batch_lookup),
            }
        )
    if not slides:
        raise RuntimeError("no slides aligned")
    if require_all_splits and not {"train", "val", "test"}.issubset(
        {slide.split for slide in slides}
    ):
        raise RuntimeError("aligned slides do not contain train/val/test")
    dimensions = {int(slide.features.shape[1]) for slide in slides}
    if len(dimensions) != 1:
        raise RuntimeError(f"feature dimensions differ: {dimensions}")
    return slides, np.asarray(reference, dtype=np.float32), genes, pd.DataFrame(coverage_records)


def _training_spatial_reliability_gene_weights(
    slides: list[PanAtlasSlide],
    consensus_basis_panel: np.ndarray,
    base_weight: np.ndarray,
    *,
    block_size: int,
    floor: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Train-only full-panel weights aligned with spatial PCC recoverability."""
    records: list[tuple[str, str, np.ndarray, np.ndarray]] = []
    genes = int(consensus_basis_panel.shape[1])
    for slide in slides:
        if slide.split != "train":
            continue
        variance = np.zeros(genes, dtype=np.float32)
        reliability = np.full(genes, np.nan, dtype=np.float32)
        for start in range(0, genes, int(block_size)):
            stop = min(start + int(block_size), genes)
            local_probability = (
                slide.local_theta_rna @ slide.local_basis_panel[:, start:stop]
            )
            consensus_probability = (
                slide.theta_rna @ slide.consensus_basis_panel[:, start:stop]
            )
            local_log = np.log1p(1.0e4 * np.maximum(local_probability, 0.0))
            consensus_log = np.log1p(
                1.0e4 * np.maximum(consensus_probability, 0.0)
            )
            variance[start:stop] = np.var(local_log, axis=0).astype(np.float32)
            reliability[start:stop] = _correlation_vector(
                local_log, consensus_log
            )
        records.append((slide.source, slide.patient, variance, reliability))

    def balanced(position: int) -> np.ndarray:
        source_values = []
        for source in sorted({value[0] for value in records}):
            patient_values = []
            for patient in sorted(
                {value[1] for value in records if value[0] == source}
            ):
                patient_values.append(
                    np.nanmean(
                        np.stack(
                            [
                                value[position]
                                for value in records
                                if value[0] == source and value[1] == patient
                            ],
                            axis=0,
                        ),
                        axis=0,
                    )
                )
            source_values.append(np.nanmean(np.stack(patient_values), axis=0))
        return np.nanmean(np.stack(source_values), axis=0).astype(np.float32)

    spatial_variance = np.maximum(balanced(2), 0.0)
    oracle_reliability = np.nan_to_num(
        balanced(3), nan=0.0, posinf=1.0, neginf=-1.0
    )
    positive = spatial_variance[spatial_variance > 1.0e-8]
    scale = float(np.median(positive)) if len(positive) else 1.0
    spatial_score = np.sqrt(spatial_variance / max(scale, EPS))
    spatial_score = np.clip(spatial_score, 0.0, 4.0)
    reliability_score = np.clip(0.5 + 0.5 * oracle_reliability, 0.05, 1.0)
    base = np.sqrt(np.maximum(np.asarray(base_weight, dtype=np.float32), EPS))
    base /= max(float(base.mean()), EPS)
    weight = (float(floor) + reliability_score * spatial_score) * base
    upper = float(np.quantile(weight, 0.995))
    weight = np.clip(weight, max(float(floor) * 0.1, EPS), max(upper, 1.0))
    weight /= max(float(weight.mean()), EPS)
    return (
        weight.astype(np.float32),
        spatial_variance.astype(np.float32),
        oracle_reliability.astype(np.float32),
    )


def _estimate_rna_per_cell(slides: list[PanAtlasSlide]) -> np.ndarray:
    train = [slide for slide in slides if slide.split == "train"]
    rna = np.sum([slide.theta_rna.sum(axis=0) for slide in train], axis=0)
    composition = np.sum(
        [slide.theta_composition.sum(axis=0) for slide in train], axis=0
    )
    ratio = rna / np.maximum(composition, 1.0e-6)
    valid = (rna > 1.0e-4) & (composition > 1.0e-4) & np.isfinite(ratio)
    median = float(np.median(ratio[valid])) if np.any(valid) else 1.0
    ratio[~valid] = median
    ratio = np.clip(ratio / max(median, EPS), 0.05, 20.0)
    return ratio.astype(np.float32)


def _program_weights(consensus_dir: Path) -> np.ndarray:
    table = pd.read_csv(consensus_dir / "consensus_program_metadata.csv")
    # A few consensus rows receive only secondary soft assignments and thus
    # have no hard-primary support count in the metadata table.  They remain
    # valid rows of W; treat their missing hard support as zero rather than
    # allowing NaNs to poison the spatial and PCC losses.
    patient = np.minimum(
        table["patients"].fillna(0.0).to_numpy(dtype=np.float64) / 3.0,
        1.0,
    )
    source = np.minimum(
        table["sources"].fillna(0.0).to_numpy(dtype=np.float64) / 2.0,
        1.0,
    )
    mass = table["rna_mass_fraction"].fillna(0.0).to_numpy(dtype=np.float64)
    mass = np.sqrt(mass / max(float(np.mean(mass)), EPS))
    weight = (0.25 + 0.75 * np.sqrt(patient * source)) * np.clip(mass, 0.25, 4.0)
    weight /= max(float(np.mean(weight)), EPS)
    if not np.all(np.isfinite(weight)) or np.any(weight <= 0.0):
        raise RuntimeError("program loss weights are not finite and positive")
    return weight.astype(np.float32)


def _to_device_graph_pair(
    graphs: tuple[
        tuple[np.ndarray, np.ndarray, np.ndarray],
        tuple[np.ndarray, np.ndarray, np.ndarray],
    ],
    device: torch.device,
) -> tuple[
    tuple[torch.Tensor, torch.Tensor, torch.Tensor],
    tuple[torch.Tensor, torch.Tensor, torch.Tensor],
]:
    return (
        _to_device_graph(graphs[0], device),
        _to_device_graph(graphs[1], device),
    )


def _aggregate_program_groups(
    theta: torch.Tensor,
    group_index: torch.Tensor,
    groups: int,
) -> torch.Tensor:
    result = theta.new_zeros((theta.shape[0], int(groups)))
    result.index_add_(1, group_index, theta)
    return result


def _relative_program_type_target(
    theta: torch.Tensor, *, lift_alpha: float
) -> torch.Tensor:
    """Emphasise locally enriched identity separately from absolute abundance."""
    prevalence = theta.mean(dim=0, keepdim=True).clamp_min(1.0e-5)
    enriched = theta.clamp_min(EPS) / prevalence.pow(float(lift_alpha))
    return enriched / enriched.sum(dim=1, keepdim=True).clamp_min(EPS)


def _cumulative_support_target(
    rna_target: torch.Tensor,
    composition_target: torch.Tensor,
    *,
    mass: float,
) -> torch.Tensor:
    combined = 0.5 * (rna_target + composition_target)
    combined = combined / combined.sum(dim=1, keepdim=True).clamp_min(EPS)
    values, indices = torch.sort(combined, dim=1, descending=True)
    previous = torch.cumsum(values, dim=1) - values
    keep_sorted = previous < float(mass)
    support = torch.zeros_like(combined)
    support.scatter_(1, indices, keep_sorted.to(combined.dtype))
    return support


def _balanced_focal_support_loss(
    logits: torch.Tensor,
    target: torch.Tensor,
    program_weight: torch.Tensor,
    *,
    gamma_positive: float = 1.0,
    gamma_negative: float = 3.0,
) -> torch.Tensor:
    probability = torch.sigmoid(logits)
    positive = -target * torch.log(probability.clamp_min(EPS)) * (
        1.0 - probability
    ).pow(float(gamma_positive))
    negative = -(1.0 - target) * torch.log((1.0 - probability).clamp_min(EPS)) * (
        probability
    ).pow(float(gamma_negative))
    positive = positive * program_weight[None, :]
    positive_loss = positive.sum() / (
        target * program_weight[None, :]
    ).sum().clamp_min(EPS)
    negative_loss = negative.sum() / (1.0 - target).sum().clamp_min(EPS)
    return 0.5 * (positive_loss + negative_loss)


def _edge_similarity_loss(
    auxiliary: dict[str, torch.Tensor],
    support_target: torch.Tensor,
    composition_target: torch.Tensor,
) -> torch.Tensor:
    source = auxiliary["edge_source"]
    target = auxiliary["edge_target"]
    predicted = auxiliary["edge_probability"]
    nonself = source != target
    if not bool(nonself.any()):
        return predicted.new_zeros(())
    source = source[nonself]
    target = target[nonself]
    predicted = predicted[nonself]
    intersection = (
        support_target[source] * support_target[target]
    ).sum(dim=1)
    union = torch.maximum(
        support_target[source], support_target[target]
    ).sum(dim=1).clamp_min(1.0)
    support_overlap = intersection / union
    root_difference = torch.sqrt(composition_target[source].clamp_min(EPS)) - torch.sqrt(
        composition_target[target].clamp_min(EPS)
    )
    hellinger_sq = 0.5 * root_difference.square().sum(dim=1)
    theta_similarity = torch.exp(-hellinger_sq / 0.20)
    edge_target = (0.5 * support_overlap + 0.5 * theta_similarity).detach()
    return F.binary_cross_entropy(predicted.clamp(EPS, 1.0 - EPS), edge_target)


def _jsd(predicted: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    predicted = predicted.clamp_min(EPS)
    target = target.clamp_min(EPS)
    middle = 0.5 * (predicted + target)
    return 0.5 * (
        (predicted * (torch.log(predicted) - torch.log(middle))).sum(dim=1)
        + (target * (torch.log(target) - torch.log(middle))).sum(dim=1)
    ).mean()


def _hellinger(predicted: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    return 0.5 * (
        torch.sqrt(predicted.clamp_min(EPS))
        - torch.sqrt(target.clamp_min(EPS))
    ).square().sum(dim=1).mean()


def _weighted_correlation_loss(
    predicted: torch.Tensor,
    target: torch.Tensor,
    weight: torch.Tensor,
    *,
    target_norm_minimum: float,
) -> torch.Tensor:
    x = target - target.mean(dim=0, keepdim=True)
    y = predicted - predicted.mean(dim=0, keepdim=True)
    x_norm = torch.sqrt(x.square().sum(dim=0).clamp_min(1.0e-12))
    y_norm = torch.sqrt(y.square().sum(dim=0).clamp_min(1.0e-12))
    valid = (x_norm > target_norm_minimum) & (y_norm > 1.0e-5)
    if not bool(valid.any()):
        return predicted.new_tensor(1.0)
    correlation = (x[:, valid] * y[:, valid]).sum(dim=0) / (
        x_norm[valid] * y_norm[valid]
    ).clamp_min(EPS)
    local_weight = weight[valid]
    return 1.0 - (
        correlation.clamp(-1.0, 1.0) * local_weight
    ).sum() / local_weight.sum().clamp_min(EPS)


def _scheduled_weights(args: argparse.Namespace, epoch: int) -> dict[str, float]:
    ramp = float(
        np.clip(
            (int(epoch) - int(args.pcc_warmup_epochs) + 1)
            / max(int(args.pcc_ramp_epochs), 1),
            0.0,
            1.0,
        )
    )
    support_gamma = float(args.support_gamma_max) * (0.25 + 0.75 * float(
        np.clip((int(epoch) - 1) / max(int(args.support_ramp_epochs), 1), 0.0, 1.0)
    ))
    domain_ramp = float(
        np.clip(
            (int(epoch) - int(args.domain_warmup_epochs))
            / max(int(args.domain_ramp_epochs), 1),
            0.0,
            1.0,
        )
    )
    return {
        "rna_jsd": float(args.rna_jsd_weight),
        "composition_jsd": float(args.composition_jsd_weight),
        "theta_pcc": float(args.theta_pcc_weight) * ramp,
        "expression_log": float(args.expression_log_weight),
        "gene_pcc": float(args.gene_pcc_weight) * ramp,
        "bulk": float(args.bulk_weight),
        "gradient": float(args.gradient_weight) * max(ramp, 0.25),
        "consistency": float(args.consistency_weight),
        "entropy": float(args.entropy_weight),
        "group_jsd": float(args.group_jsd_weight),
        "program_type": float(args.program_type_weight) * max(ramp, 0.25),
        "support": float(args.support_weight),
        "edge": float(args.edge_weight) * max(ramp, 0.25),
        "domain": float(args.domain_weight) * domain_ramp,
        "support_gamma": support_gamma,
        "domain_ramp": domain_ramp,
        "pcc_ramp": ramp,
    }


def _loss(
    *,
    rna_logits: torch.Tensor,
    composition_logits: torch.Tensor,
    auxiliary: dict[str, torch.Tensor],
    rna_target: torch.Tensor,
    composition_target: torch.Tensor,
    local_expression_target: torch.Tensor,
    consensus_basis_panel: torch.Tensor,
    gene_weight: torch.Tensor,
    program_weight: torch.Tensor,
    rna_per_cell: torch.Tensor,
    program_group_index: torch.Tensor,
    graphs: tuple[
        tuple[torch.Tensor, torch.Tensor, torch.Tensor],
        tuple[torch.Tensor, torch.Tensor, torch.Tensor],
    ],
    type_lift_alpha: float,
    source_target: torch.Tensor,
    technology_target: torch.Tensor,
    support_mass: float,
    weights: dict[str, float],
) -> tuple[torch.Tensor, dict[str, float]]:
    rna = torch.softmax(rna_logits.float(), dim=1)
    composition = torch.softmax(composition_logits.float(), dim=1)
    predicted_expression = rna @ consensus_basis_panel
    predicted_log = torch.log1p(1.0e4 * predicted_expression)
    target_log = torch.log1p(1.0e4 * local_expression_target.clamp_min(0.0))

    rna_jsd = _jsd(rna, rna_target)
    composition_jsd = _jsd(composition, composition_target)
    rna_pcc = _weighted_correlation_loss(
        torch.sqrt(rna.clamp_min(EPS)),
        torch.sqrt(rna_target.clamp_min(EPS)),
        program_weight,
        target_norm_minimum=1.0e-3,
    )
    composition_pcc = _weighted_correlation_loss(
        torch.sqrt(composition.clamp_min(EPS)),
        torch.sqrt(composition_target.clamp_min(EPS)),
        program_weight,
        target_norm_minimum=1.0e-3,
    )
    theta_pcc = 0.5 * (rna_pcc + composition_pcc)
    support_target = _cumulative_support_target(
        rna_target,
        composition_target,
        mass=float(support_mass),
    )
    support = _balanced_focal_support_loss(
        auxiliary["support_logits"],
        support_target,
        program_weight,
    )
    edge = _edge_similarity_loss(auxiliary, support_target, composition_target)
    domain = 0.5 * (
        F.cross_entropy(auxiliary["source_logits"], source_target.reshape(1))
        + F.cross_entropy(
            auxiliary["technology_logits"], technology_target.reshape(1)
        )
    )

    groups = int(auxiliary["rna_group_logits"].shape[1])
    rna_group_target = _aggregate_program_groups(
        rna_target, program_group_index, groups
    )
    composition_group_target = _aggregate_program_groups(
        composition_target, program_group_index, groups
    )
    group_jsd = 0.5 * (
        _jsd(
            torch.softmax(auxiliary["rna_group_logits"].float(), dim=1),
            rna_group_target,
        )
        + _jsd(
            torch.softmax(auxiliary["composition_group_logits"].float(), dim=1),
            composition_group_target,
        )
    )
    rna_type_target = _relative_program_type_target(
        rna_target, lift_alpha=type_lift_alpha
    )
    composition_type_target = _relative_program_type_target(
        composition_target, lift_alpha=type_lift_alpha
    )
    program_type = 0.5 * (
        _jsd(
            torch.softmax(auxiliary["rna_type_logits"].float(), dim=1),
            rna_type_target,
        )
        + _jsd(
            torch.softmax(auxiliary["composition_type_logits"].float(), dim=1),
            composition_type_target,
        )
    )
    expression_per_gene = F.smooth_l1_loss(
        predicted_log, target_log, reduction="none", beta=0.25
    ).mean(dim=0)
    expression_log = (
        expression_per_gene * gene_weight
    ).sum() / gene_weight.sum().clamp_min(EPS)
    gene_pcc = _weighted_correlation_loss(
        predicted_log,
        target_log,
        gene_weight,
        target_norm_minimum=1.0e-3,
    )
    bulk = 0.5 * (
        _hellinger(rna.mean(dim=0, keepdim=True), rna_target.mean(dim=0, keepdim=True))
        + _hellinger(
            composition.mean(dim=0, keepdim=True),
            composition_target.mean(dim=0, keepdim=True),
        )
    )

    graph_gradients: list[torch.Tensor] = []
    for source, target_edge, graph_weight in graphs:
        nonself = source != target_edge
        if not bool(nonself.any()):
            continue
        source = source[nonself]
        target_edge = target_edge[nonself]
        edge_weight = graph_weight[nonself]
        predicted_gradient = torch.sqrt(rna[source].clamp_min(EPS)) - torch.sqrt(
            rna[target_edge].clamp_min(EPS)
        )
        target_gradient = torch.sqrt(rna_target[source].clamp_min(EPS)) - torch.sqrt(
            rna_target[target_edge].clamp_min(EPS)
        )
        gradient_per_program = F.smooth_l1_loss(
            predicted_gradient, target_gradient, reduction="none", beta=0.05
        )
        gradient_per_edge = (
            gradient_per_program * program_weight[None, :]
        ).sum(dim=1) / program_weight.sum().clamp_min(EPS)
        graph_gradients.append(
            (gradient_per_edge * edge_weight).sum()
            / edge_weight.sum().clamp_min(EPS)
        )
    if graph_gradients:
        gradient = 0.7 * graph_gradients[0]
        if len(graph_gradients) > 1:
            gradient = gradient + 0.3 * graph_gradients[1]
    else:
        gradient = rna_jsd.new_zeros(())

    rna_from_composition = composition * rna_per_cell[None, :]
    rna_from_composition /= rna_from_composition.sum(dim=1, keepdim=True).clamp_min(EPS)
    consistency = _hellinger(rna, rna_from_composition)
    log_k = math.log(max(rna.shape[1], 2))
    predicted_entropy = -(rna * torch.log(rna.clamp_min(EPS))).sum(dim=1) / log_k
    target_entropy = -(
        rna_target * torch.log(rna_target.clamp_min(EPS))
    ).sum(dim=1) / log_k
    entropy = F.smooth_l1_loss(predicted_entropy, target_entropy, beta=0.05)

    total = (
        weights["rna_jsd"] * rna_jsd
        + weights["composition_jsd"] * composition_jsd
        + weights["theta_pcc"] * theta_pcc
        + weights["expression_log"] * expression_log
        + weights["gene_pcc"] * gene_pcc
        + weights["bulk"] * bulk
        + weights["gradient"] * gradient
        + weights["consistency"] * consistency
        + weights["entropy"] * entropy
        + weights["group_jsd"] * group_jsd
        + weights["program_type"] * program_type
        + weights["support"] * support
        + weights["edge"] * edge
        + weights["domain"] * domain
    )
    parts = {
        "rna_jsd": float(rna_jsd.detach().cpu()),
        "composition_jsd": float(composition_jsd.detach().cpu()),
        "theta_pcc_loss": float(theta_pcc.detach().cpu()),
        "expression_log_huber": float(expression_log.detach().cpu()),
        "gene_pcc_loss": float(gene_pcc.detach().cpu()),
        "bulk_hellinger": float(bulk.detach().cpu()),
        "gradient": float(gradient.detach().cpu()),
        "rna_composition_consistency": float(consistency.detach().cpu()),
        "entropy_match": float(entropy.detach().cpu()),
        "program_group_jsd": float(group_jsd.detach().cpu()),
        "relative_program_type_jsd": float(program_type.detach().cpu()),
        "support_focal": float(support.detach().cpu()),
        "edge_similarity": float(edge.detach().cpu()),
        "domain_classification": float(domain.detach().cpu()),
    }
    return total, parts


def _correlation_vector(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    left -= left.mean(axis=0, keepdims=True)
    right -= right.mean(axis=0, keepdims=True)
    denominator = np.sqrt(
        np.square(left).sum(axis=0) * np.square(right).sum(axis=0)
    )
    result = np.full(left.shape[1], np.nan, dtype=np.float32)
    valid = denominator > 1.0e-12
    result[valid] = (
        (left[:, valid] * right[:, valid]).sum(axis=0) / denominator[valid]
    ).astype(np.float32)
    return result


def _balanced_pcc_summary(
    records: list[tuple[str, str, np.ndarray]]
) -> tuple[float, dict[str, float], np.ndarray]:
    source_vectors: list[np.ndarray] = []
    by_source: dict[str, float] = {}
    for source in sorted({value[0] for value in records}):
        patient_vectors: list[np.ndarray] = []
        for patient in sorted(
            {value[1] for value in records if value[0] == source}
        ):
            patient_vectors.append(
                np.nanmean(
                    np.stack(
                        [
                            vector
                            for current_source, current_patient, vector in records
                            if current_source == source and current_patient == patient
                        ],
                        axis=0,
                    ),
                    axis=0,
                )
            )
        source_vector = np.nanmean(np.stack(patient_vectors, axis=0), axis=0)
        source_vectors.append(source_vector)
        by_source[source] = float(np.nanmedian(source_vector))
    balanced = np.nanmean(np.stack(source_vectors, axis=0), axis=0)
    return float(np.nanmedian(balanced)), by_source, balanced.astype(np.float32)


@torch.no_grad()
def _evaluate(
    model: PanAtlasDualMapper,
    inr: BoundedINR | None,
    slides: list[PanAtlasSlide],
    *,
    split: str,
    consensus_basis_panel: torch.Tensor,
    gene_weight: torch.Tensor,
    program_weight: torch.Tensor,
    rna_per_cell: torch.Tensor,
    device: torch.device,
    weights: dict[str, float],
) -> tuple[float, dict[str, object], dict[str, tuple[np.ndarray, np.ndarray]]]:
    model.eval()
    if inr is not None:
        inr.eval()
    losses: list[float] = []
    slide_records: list[dict[str, object]] = []
    pcc_records: list[tuple[str, str, np.ndarray]] = []
    predictions: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for slide in slides:
        if slide.split != split:
            continue
        features = torch.from_numpy(np.asarray(slide.features, dtype=np.float32)).to(device)
        graphs = _to_device_graph_pair(slide.graphs, device)
        rna_logits, composition_logits, context, auxiliary = model(
            features,
            graphs,
            support_gamma=float(weights["support_gamma"]),
            domain_strength=0.0,
        )
        correction = None
        if inr is not None:
            coordinates = torch.from_numpy(
                normalize_slide_coordinates(slide.coordinates)
            ).to(device)
            correction = inr(coordinates, context)
            rna_correction, composition_correction = correction.split(
                model.programs, dim=1
            )
            rna_logits = rna_logits + rna_correction
            composition_logits = composition_logits + composition_correction
        rna_target = torch.from_numpy(slide.theta_rna).to(device)
        composition_target = torch.from_numpy(slide.theta_composition).to(device)
        local_expression = (
            torch.from_numpy(slide.local_theta_rna).to(device)
            @ torch.from_numpy(slide.local_basis_panel).to(device)
        )
        slide_consensus_basis_panel = torch.from_numpy(
            slide.consensus_basis_panel
        ).to(device)
        loss, parts = _loss(
            rna_logits=rna_logits,
            composition_logits=composition_logits,
            auxiliary=auxiliary,
            rna_target=rna_target,
            composition_target=composition_target,
            local_expression_target=local_expression,
            consensus_basis_panel=slide_consensus_basis_panel,
            gene_weight=gene_weight,
            program_weight=program_weight,
            rna_per_cell=rna_per_cell,
            program_group_index=model.program_group_index,
            graphs=graphs,
            type_lift_alpha=float(model.type_lift_alpha),
            source_target=torch.tensor(slide.source_index, device=device),
            technology_target=torch.tensor(slide.technology_index, device=device),
            support_mass=float(model.support_mass),
            weights=weights,
        )
        predicted_rna = torch.softmax(rna_logits, dim=1)
        predicted_composition = torch.softmax(composition_logits, dim=1)
        predicted_log = torch.log1p(
            1.0e4 * (predicted_rna @ slide_consensus_basis_panel)
        )
        target_log = torch.log1p(1.0e4 * local_expression)
        gene_pcc = _correlation_vector(
            target_log.cpu().numpy(), predicted_log.cpu().numpy()
        )
        theta_rna_pcc = _correlation_vector(
            np.sqrt(np.maximum(slide.theta_rna, EPS)),
            np.sqrt(np.maximum(predicted_rna.cpu().numpy(), EPS)),
        )
        theta_composition_pcc = _correlation_vector(
            np.sqrt(np.maximum(slide.theta_composition, EPS)),
            np.sqrt(np.maximum(predicted_composition.cpu().numpy(), EPS)),
        )
        record = {
            "slide_key": slide.key,
            "source": slide.source,
            "patient": slide.patient,
            "spots": int(len(slide.features)),
            "loss": float(loss.cpu()),
            "loss_panel_gene_pcc_median": float(np.nanmedian(gene_pcc)),
            "rna_program_pcc_median": float(np.nanmedian(theta_rna_pcc)),
            "composition_program_pcc_median": float(
                np.nanmedian(theta_composition_pcc)
            ),
            **parts,
            "inr_mean_abs_correction": (
                float(correction.abs().mean().cpu()) if correction is not None else 0.0
            ),
        }
        losses.append(float(loss.cpu()))
        slide_records.append(record)
        pcc_records.append((slide.source, slide.patient, gene_pcc))
        predictions[slide.key] = (
            predicted_rna.cpu().numpy().astype(np.float32),
            predicted_composition.cpu().numpy().astype(np.float32),
        )
    median, by_source, _ = _balanced_pcc_summary(pcc_records)
    summary = {
        "split": split,
        "slides": slide_records,
        "mean_loss": float(np.mean(losses)),
        "patient_source_balanced_loss_panel_gene_pcc_median": median,
        "by_source_patient_balanced_loss_panel_gene_pcc_median": by_source,
    }
    return float(np.mean(losses)), summary, predictions


def _full_10k_evaluation(
    *,
    slides: list[PanAtlasSlide],
    split: str,
    predictions: dict[str, tuple[np.ndarray, np.ndarray]],
    consensus_reference: np.ndarray,
    genes: np.ndarray,
    block_size: int,
) -> tuple[dict[str, object], pd.DataFrame]:
    predicted_records: list[tuple[str, str, np.ndarray]] = []
    oracle_records: list[tuple[str, str, np.ndarray]] = []
    slide_summaries: list[dict[str, object]] = []
    for slide in slides:
        if slide.split != split:
            continue
        predicted_rna = predictions[slide.key][0]
        local_basis = np.load(
            slide.slide_dir / "basis_probability_mean.npy", mmap_mode="r"
        )
        predicted_pcc = np.full(len(genes), np.nan, dtype=np.float32)
        oracle_pcc = np.full(len(genes), np.nan, dtype=np.float32)
        for start in range(0, len(genes), int(block_size)):
            stop = min(start + int(block_size), len(genes))
            target = slide.local_theta_rna @ np.asarray(
                local_basis[:, start:stop], dtype=np.float32
            )
            reference_block = np.asarray(
                slide.consensus_reference[:, start:stop], dtype=np.float32
            )
            predicted = predicted_rna @ reference_block
            oracle = slide.theta_rna @ reference_block
            target_log = np.log1p(1.0e4 * target)
            predicted_pcc[start:stop] = _correlation_vector(
                target_log, np.log1p(1.0e4 * predicted)
            )
            oracle_pcc[start:stop] = _correlation_vector(
                target_log, np.log1p(1.0e4 * oracle)
            )
        predicted_records.append((slide.source, slide.patient, predicted_pcc))
        oracle_records.append((slide.source, slide.patient, oracle_pcc))
        slide_summaries.append(
            {
                "slide_key": slide.key,
                "source": slide.source,
                "patient": slide.patient,
                "spots": int(len(slide.features)),
                "predicted_gene_pcc_median": float(np.nanmedian(predicted_pcc)),
                "oracle_consensus_gene_pcc_median": float(np.nanmedian(oracle_pcc)),
            }
        )
    predicted_median, predicted_by_source, predicted_vector = _balanced_pcc_summary(
        predicted_records
    )
    oracle_median, oracle_by_source, oracle_vector = _balanced_pcc_summary(
        oracle_records
    )
    gene_table = pd.DataFrame(
        {
            "gene_index": np.arange(len(genes), dtype=np.int64),
            "gene": genes,
            "predicted_patient_source_balanced_pcc": predicted_vector,
            "oracle_consensus_patient_source_balanced_pcc": oracle_vector,
        }
    )
    summary = {
        "split": split,
        "reference_for_pcc": "original_per_slide_top10k_theta_rna_times_basis_probability",
        "prediction_decoder": "predicted_theta_rna_times_slide_matched_batch_adapted_pan_atlas_W",
        "patient_source_balanced_gene_pcc_median": predicted_median,
        "by_source_patient_balanced_gene_pcc_median": predicted_by_source,
        "oracle_consensus_gene_pcc_median": oracle_median,
        "oracle_by_source_gene_pcc_median": oracle_by_source,
        "slides": slide_summaries,
    }
    return summary, gene_table


def _source_weights(
    slides: list[PanAtlasSlide], *, power: float
) -> dict[str, float]:
    train = [slide for slide in slides if slide.split == "train"]
    counts = {
        source: sum(slide.source == source for slide in train)
        for source in sorted({slide.source for slide in train})
    }
    raw = {
        source: (len(train) / (len(counts) * count)) ** float(power)
        for source, count in counts.items()
    }
    normalizer = sum(counts[key] * raw[key] for key in counts) / len(train)
    return {key: float(raw[key] / max(normalizer, EPS)) for key in raw}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--consensus-dir", type=Path, required=True)
    parser.add_argument("--per-slide-dir", type=Path, required=True)
    parser.add_argument("--master-spot-index", type=Path, required=True)
    parser.add_argument("--split-manifest", type=Path, required=True)
    parser.add_argument("--feature-root", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--sources",
        nargs="+",
        default=["HEST", "GSE210616", "GSE213688"],
    )
    parser.add_argument(
        "--loss-genes",
        type=int,
        default=10000,
        help="Genes used by expression-log and gene-PCC losses; default is full 10K.",
    )
    parser.add_argument("--epochs", type=int, default=55)
    parser.add_argument("--patience", type=int, default=10)
    parser.add_argument("--inr-epochs", type=int, default=8)
    parser.add_argument("--hidden", type=int, default=256)
    parser.add_argument("--graph-layers", type=int, default=2)
    parser.add_argument("--context-graph-layers", type=int, default=1)
    parser.add_argument("--graph-heads", type=int, default=4)
    parser.add_argument("--neighbors", type=int, default=6)
    parser.add_argument("--context-neighbors", type=int, default=18)
    parser.add_argument("--program-groups", type=int, default=64)
    parser.add_argument(
        "--program-group-method",
        choices=["expression", "morphology"],
        default="expression",
        help="Use train-patient-only residual HE directions for morphology groups.",
    )
    parser.add_argument(
        "--hierarchical-logits",
        action="store_true",
        help="Lift group/type logits into the final fine-program logits.",
    )
    parser.add_argument("--learning-rate", type=float, default=3.0e-4)
    parser.add_argument("--inr-learning-rate", type=float, default=5.0e-4)
    parser.add_argument("--source-balance-power", type=float, default=0.25)
    parser.add_argument("--rna-jsd-weight", type=float, default=0.70)
    parser.add_argument("--composition-jsd-weight", type=float, default=0.45)
    parser.add_argument("--theta-pcc-weight", type=float, default=1.25)
    parser.add_argument("--expression-log-weight", type=float, default=1.00)
    parser.add_argument("--gene-pcc-weight", type=float, default=3.50)
    parser.add_argument("--bulk-weight", type=float, default=0.30)
    parser.add_argument("--gradient-weight", type=float, default=0.10)
    parser.add_argument("--consistency-weight", type=float, default=0.20)
    parser.add_argument("--entropy-weight", type=float, default=0.10)
    parser.add_argument("--group-jsd-weight", type=float, default=0.0)
    parser.add_argument("--program-type-weight", type=float, default=0.0)
    parser.add_argument("--type-lift-alpha", type=float, default=0.50)
    parser.add_argument("--support-weight", type=float, default=0.80)
    parser.add_argument(
        "--support-gamma-max",
        type=float,
        default=1.0,
        help="Maximum support-logit gate applied to abundance logits; use 0 to disable gating.",
    )
    parser.add_argument("--edge-weight", type=float, default=0.40)
    parser.add_argument("--domain-weight", type=float, default=0.08)
    parser.add_argument("--support-mass", type=float, default=0.90)
    parser.add_argument("--support-ramp-epochs", type=int, default=8)
    parser.add_argument("--domain-warmup-epochs", type=int, default=5)
    parser.add_argument("--domain-ramp-epochs", type=int, default=12)
    parser.add_argument("--pcc-warmup-epochs", type=int, default=3)
    parser.add_argument("--pcc-ramp-epochs", type=int, default=8)
    parser.add_argument("--full-evaluation-block", type=int, default=512)
    parser.add_argument("--gene-weight-block", type=int, default=512)
    parser.add_argument("--gene-weight-floor", type=float, default=0.15)
    parser.add_argument("--disable-spatial-gene-weighting", action="store_true")
    parser.add_argument("--seed", type=int, default=260802)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--final-train-all",
        action="store_true",
        help=(
            "Train exactly --epochs passes on an all-train development manifest, "
            "without validation selection or test evaluation."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.time()
    np.random.seed(int(args.seed))
    torch.manual_seed(int(args.seed))
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is required")
    torch.set_float32_matmul_precision("high")
    device = torch.device("cuda")

    panel_np, gene_weight_np, panel_table = _select_loss_panel(
        args.consensus_dir, count=int(args.loss_genes)
    )
    slides, consensus_reference_np, genes, coverage = _load_slides(
        consensus_dir=args.consensus_dir,
        per_slide_dir=args.per_slide_dir,
        master_spot_index=args.master_spot_index,
        split_manifest=args.split_manifest,
        panel=panel_np,
        feature_root=args.feature_root,
        neighbors=int(args.neighbors),
        context_neighbors=int(args.context_neighbors),
        sources=tuple(str(value) for value in args.sources),
        require_all_splits=not bool(args.final_train_all),
    )
    if bool(args.final_train_all):
        if {slide.split for slide in slides} != {"train"}:
            raise RuntimeError("--final-train-all requires every aligned slide to be train")
        if int(args.inr_epochs) != 0:
            raise RuntimeError("--final-train-all requires --inr-epochs 0")
    if len(panel_np) == 10000 and not bool(args.disable_spatial_gene_weighting):
        (
            gene_weight_np,
            train_spatial_variance,
            train_oracle_reliability,
        ) = _training_spatial_reliability_gene_weights(
            slides,
            consensus_reference_np[:, panel_np],
            gene_weight_np,
            block_size=int(args.gene_weight_block),
            floor=float(args.gene_weight_floor),
        )
        panel_table["train_balanced_spatial_variance"] = train_spatial_variance
        panel_table["train_balanced_oracle_reliability"] = train_oracle_reliability
        panel_table["loss_weight"] = gene_weight_np
    rna_per_cell_np = _estimate_rna_per_cell(slides)
    program_weight_np = _program_weights(args.consensus_dir)
    source_weight = _source_weights(
        slides, power=float(args.source_balance_power)
    )
    programs = int(consensus_reference_np.shape[0])
    source_levels = [
        value
        for value, _ in sorted(
            {slide.source: slide.source_index for slide in slides}.items(),
            key=lambda item: item[1],
        )
    ]
    technology_levels = [
        value
        for value, _ in sorted(
            {slide.technology: slide.technology_index for slide in slides}.items(),
            key=lambda item: item[1],
        )
    ]
    morphology_prototypes_np: np.ndarray | None = None
    morphology_similarity_np: np.ndarray | None = None
    if str(args.program_group_method) == "morphology":
        (
            program_groups_np,
            morphology_prototypes_np,
            morphology_similarity_np,
        ) = _morphology_program_groups(
            slides,
            consensus_reference_np,
            count=int(args.program_groups),
        )
    else:
        program_groups_np = _program_groups(
            consensus_reference_np, count=int(args.program_groups)
        )
    input_dim = int(slides[0].features.shape[1])
    model = PanAtlasDualMapper(
        input_dim,
        hidden=int(args.hidden),
        programs=programs,
        program_groups=program_groups_np,
        graph_layers=int(args.graph_layers),
        context_graph_layers=int(args.context_graph_layers),
        graph_heads=int(args.graph_heads),
        dropout=0.10,
        type_lift_alpha=float(args.type_lift_alpha),
        sources=len(source_levels),
        technologies=len(technology_levels),
        support_mass=float(args.support_mass),
        hierarchical_logits=bool(args.hierarchical_logits),
    ).to(device)
    consensus_basis_panel = torch.from_numpy(
        consensus_reference_np[:, panel_np]
    ).to(device)
    gene_weight = torch.from_numpy(gene_weight_np).to(device)
    program_weight = torch.from_numpy(program_weight_np).to(device)
    rna_per_cell = torch.from_numpy(rna_per_cell_np).to(device)

    if args.dry_run:
        slide = next(value for value in slides if value.split == "train")
        model.train()
        features = torch.from_numpy(np.asarray(slide.features, dtype=np.float32)).to(device)
        graphs = _to_device_graph_pair(slide.graphs, device)
        dry_weights = _scheduled_weights(
            args, max(int(args.pcc_warmup_epochs), 1)
        )
        rna_logits, composition_logits, _, auxiliary = model(
            features,
            graphs,
            support_gamma=float(dry_weights["support_gamma"]),
            domain_strength=float(dry_weights["domain_ramp"]),
        )
        local_expression = (
            torch.from_numpy(slide.local_theta_rna).to(device)
            @ torch.from_numpy(slide.local_basis_panel).to(device)
        )
        slide_consensus_basis_panel = torch.from_numpy(
            slide.consensus_basis_panel
        ).to(device)
        loss, parts = _loss(
            rna_logits=rna_logits,
            composition_logits=composition_logits,
            auxiliary=auxiliary,
            rna_target=torch.from_numpy(slide.theta_rna).to(device),
            composition_target=torch.from_numpy(slide.theta_composition).to(device),
            local_expression_target=local_expression,
            consensus_basis_panel=slide_consensus_basis_panel,
            gene_weight=gene_weight,
            program_weight=program_weight,
            rna_per_cell=rna_per_cell,
            program_group_index=model.program_group_index,
            graphs=graphs,
            type_lift_alpha=float(args.type_lift_alpha),
            source_target=torch.tensor(slide.source_index, device=device),
            technology_target=torch.tensor(slide.technology_index, device=device),
            support_mass=float(args.support_mass),
            weights=dry_weights,
        )
        loss.backward()
        finite_gradients = all(
            parameter.grad is None or bool(torch.isfinite(parameter.grad).all())
            for parameter in model.parameters()
        )
        validation_slide = min(
            (value for value in slides if value.split == "val"),
            key=lambda value: len(value.features),
        )
        validation_loss, validation_summary, validation_prediction = _evaluate(
            model,
            None,
            [validation_slide],
            split="val",
            consensus_basis_panel=consensus_basis_panel,
            gene_weight=gene_weight,
            program_weight=program_weight,
            rna_per_cell=rna_per_cell,
            device=device,
            weights=_scheduled_weights(args, max(int(args.pcc_warmup_epochs), 1)),
        )
        full_summary, _ = _full_10k_evaluation(
            slides=[validation_slide],
            split="val",
            predictions=validation_prediction,
            consensus_reference=consensus_reference_np,
            genes=genes,
            block_size=int(args.full_evaluation_block),
        )
        print(
            json.dumps(
                {
                    "dry_run": True,
                    "slides": len(slides),
                    "programs": programs,
                    "loss_genes": len(panel_np),
                    "input_dim": input_dim,
                    "sample": slide.key,
                    "loss": float(loss.detach().cpu()),
                    "loss_parts": parts,
                    "finite_gradients": finite_gradients,
                    "validation_smoke_slide": validation_slide.key,
                    "validation_loss": validation_loss,
                    "validation_loss_panel_pcc": validation_summary[
                        "patient_source_balanced_loss_panel_gene_pcc_median"
                    ],
                    "validation_full10k_pcc": full_summary[
                        "patient_source_balanced_gene_pcc_median"
                    ],
                    "validation_full10k_oracle_pcc": full_summary[
                        "oracle_consensus_gene_pcc_median"
                    ],
                    "peak_cuda_memory_bytes": int(torch.cuda.max_memory_allocated()),
                },
                indent=2,
            ),
            flush=True,
        )
        return

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(f"refusing non-empty output {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    panel_table.to_csv(args.output_dir / "loss_gene_panel.csv", index=False)
    coverage.to_csv(args.output_dir / "aligned_slide_coverage.csv", index=False)
    np.save(args.output_dir / "loss_gene_indices.npy", panel_np)
    np.save(args.output_dir / "loss_gene_weight.npy", gene_weight_np)
    np.save(args.output_dir / "rna_per_cell_program_factor.npy", rna_per_cell_np)
    np.save(args.output_dir / "program_loss_weight.npy", program_weight_np)
    np.save(args.output_dir / "program_group_index.npy", program_groups_np)
    if morphology_prototypes_np is not None:
        np.save(
            args.output_dir / "train_only_morphology_program_prototypes.npy",
            morphology_prototypes_np,
        )
    if morphology_similarity_np is not None:
        np.save(
            args.output_dir / "train_only_morphology_program_similarity.npy",
            morphology_similarity_np,
        )

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=float(args.learning_rate), weight_decay=1.0e-4
    )
    train_slides = [slide for slide in slides if slide.split == "train"]
    rng = np.random.default_rng(int(args.seed))
    best_state = copy.deepcopy(model.state_dict())
    best_pcc = -math.inf
    best_loss = math.inf
    stale = 0
    history: list[dict[str, float]] = []
    for epoch in range(1, int(args.epochs) + 1):
        model.train()
        weights = _scheduled_weights(args, epoch)
        loss_values: list[float] = []
        part_values: dict[str, list[float]] = {}
        for index in rng.permutation(len(train_slides)):
            slide = train_slides[int(index)]
            features = torch.from_numpy(
                np.asarray(slide.features, dtype=np.float32)
            ).to(device)
            graphs = _to_device_graph_pair(slide.graphs, device)
            optimizer.zero_grad(set_to_none=True)
            rna_logits, composition_logits, _, auxiliary = model(
                features,
                graphs,
                support_gamma=float(weights["support_gamma"]),
                domain_strength=float(weights["domain_ramp"]),
            )
            local_expression = (
                torch.from_numpy(slide.local_theta_rna).to(device)
                @ torch.from_numpy(slide.local_basis_panel).to(device)
            )
            slide_consensus_basis_panel = torch.from_numpy(
                slide.consensus_basis_panel
            ).to(device)
            loss, parts = _loss(
                rna_logits=rna_logits,
                composition_logits=composition_logits,
                auxiliary=auxiliary,
                rna_target=torch.from_numpy(slide.theta_rna).to(device),
                composition_target=torch.from_numpy(slide.theta_composition).to(device),
                local_expression_target=local_expression,
                consensus_basis_panel=slide_consensus_basis_panel,
                gene_weight=gene_weight,
                program_weight=program_weight,
                rna_per_cell=rna_per_cell,
                program_group_index=model.program_group_index,
                graphs=graphs,
                type_lift_alpha=float(args.type_lift_alpha),
                source_target=torch.tensor(slide.source_index, device=device),
                technology_target=torch.tensor(slide.technology_index, device=device),
                support_mass=float(args.support_mass),
                weights=weights,
            )
            loss = loss * float(source_weight[slide.source])
            if not bool(torch.isfinite(loss)):
                raise FloatingPointError(f"non-finite loss for {slide.key}")
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            loss_values.append(float(loss.detach().cpu()))
            for name, value in parts.items():
                part_values.setdefault(name, []).append(float(value))
        if bool(args.final_train_all):
            record = {
                "epoch": float(epoch),
                "pcc_ramp": float(weights["pcc_ramp"]),
                "train_loss": float(np.mean(loss_values)),
                "support_gamma": float(weights["support_gamma"]),
                "domain_ramp": float(weights["domain_ramp"]),
            }
            record.update(
                {
                    f"train_{key}": float(np.mean(value))
                    for key, value in part_values.items()
                }
            )
            history.append(record)
            best_state = copy.deepcopy(model.state_dict())
            print(json.dumps({"final_train_all": record}), flush=True)
            continue
        validation_loss, validation, _ = _evaluate(
            model,
            None,
            slides,
            split="val",
            consensus_basis_panel=consensus_basis_panel,
            gene_weight=gene_weight,
            program_weight=program_weight,
            rna_per_cell=rna_per_cell,
            device=device,
            weights=weights,
        )
        validation_pcc = float(
            validation["patient_source_balanced_loss_panel_gene_pcc_median"]
        )
        record: dict[str, float] = {
            "epoch": float(epoch),
            "pcc_ramp": float(weights["pcc_ramp"]),
            "train_loss": float(np.mean(loss_values)),
            "validation_loss": validation_loss,
            "validation_gene_pcc": validation_pcc,
            "support_gamma": float(weights["support_gamma"]),
            "domain_ramp": float(weights["domain_ramp"]),
        }
        record.update(
            {f"train_{key}": float(np.mean(value)) for key, value in part_values.items()}
        )
        history.append(record)
        print(json.dumps(record), flush=True)
        improved = validation_pcc > best_pcc + 1.0e-4 or (
            abs(validation_pcc - best_pcc) <= 1.0e-4
            and validation_loss < best_loss
        )
        if improved:
            best_pcc = validation_pcc
            best_loss = validation_loss
            best_state = copy.deepcopy(model.state_dict())
            stale = 0
        else:
            stale += 1
            if stale >= int(args.patience):
                break
    model.load_state_dict(best_state)
    final_weights = _scheduled_weights(args, int(args.epochs))
    if bool(args.final_train_all):
        pd.DataFrame(history).to_csv(
            args.output_dir / "training_history.csv", index=False
        )
        pd.DataFrame().to_csv(args.output_dir / "inr_history.csv", index=False)
        program_diagnostics = pd.DataFrame(
            {
                "program_index": np.arange(programs, dtype=np.int64),
                "program_group": program_groups_np,
                "loss_weight": program_weight_np,
                "rna_local_graph_gate": torch.sigmoid(model.rna_local_refine)
                .detach()
                .cpu()
                .numpy(),
                "rna_context_graph_gate": torch.sigmoid(model.rna_context_refine)
                .detach()
                .cpu()
                .numpy(),
                "composition_local_graph_gate": torch.sigmoid(
                    model.composition_local_refine
                )
                .detach()
                .cpu()
                .numpy(),
                "composition_context_graph_gate": torch.sigmoid(
                    model.composition_context_refine
                )
                .detach()
                .cpu()
                .numpy(),
            }
        )
        program_diagnostics.to_csv(
            args.output_dir / "learned_program_graph_gates.csv", index=False
        )
        torch.save(
            {
                "model": model.state_dict(),
                "inr": None,
                "inr_accepted": False,
                "input_dim": input_dim,
                "programs": programs,
                "hidden": int(args.hidden),
                "graph_layers": int(args.graph_layers),
                "context_graph_layers": int(args.context_graph_layers),
                "graph_heads": int(args.graph_heads),
                "context_neighbors": int(args.context_neighbors),
                "program_groups": int(program_groups_np.max()) + 1,
                "type_lift_alpha": float(args.type_lift_alpha),
                "sources": len(source_levels),
                "technologies": len(technology_levels),
                "source_levels": source_levels,
                "technology_levels": technology_levels,
                "support_mass": float(args.support_mass),
                "hierarchical_logits": bool(args.hierarchical_logits),
                "program_group_method": str(args.program_group_method),
                "support_gamma": float(final_weights["support_gamma"]),
                "neighbors": int(args.neighbors),
                "consensus_reference_sha256": _sha256(
                    args.consensus_dir / "consensus_reference_probability_10k.npy"
                ),
            },
            args.output_dir / "model_checkpoint.pt",
        )
        payload = {
            "method": f"virchow2_pan_atlas_g{programs}_full_development_final_mapper",
            "purpose": "deployment_checkpoint_not_used_for_oof_reporting",
            "selection_rule": "fixed epoch count from median inner-validation selected epoch",
            "fixed_epochs": int(args.epochs),
            "train_slides": int(len(train_slides)),
            "train_spots": int(sum(len(slide.features) for slide in train_slides)),
            "target_programs": programs,
            "expression_loss_genes": int(len(panel_np)),
            "gene_pcc_loss_genes": int(len(panel_np)),
            "program_groups": int(program_groups_np.max()) + 1,
            "sources": source_levels,
            "technologies": technology_levels,
            "source_training_weights": source_weight,
            "training_configuration": {
                key: str(value) if isinstance(value, Path) else value
                for key, value in vars(args).items()
            },
            "history": history,
            "runtime_seconds": time.time() - started,
            "peak_cuda_memory_bytes": int(torch.cuda.max_memory_allocated()),
        }
        (args.output_dir / "metrics.json").write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(json.dumps(payload, indent=2, ensure_ascii=False), flush=True)
        return
    base_val_loss, base_val, _ = _evaluate(
        model,
        None,
        slides,
        split="val",
        consensus_basis_panel=consensus_basis_panel,
        gene_weight=gene_weight,
        program_weight=program_weight,
        rna_per_cell=rna_per_cell,
        device=device,
        weights=final_weights,
    )
    base_val_pcc = float(
        base_val["patient_source_balanced_loss_panel_gene_pcc_median"]
    )

    selected_inr: BoundedINR | None = None
    inr_history: list[dict[str, float]] = []
    inr_accepted = False
    if int(args.inr_epochs) > 0:
        inr = BoundedINR(int(args.hidden), programs * 2).to(device)
        inr_optimizer = torch.optim.Adam(
            inr.parameters(), lr=float(args.inr_learning_rate)
        )
        best_inr_state = copy.deepcopy(inr.state_dict())
        best_inr_pcc = base_val_pcc
        best_inr_loss = base_val_loss
        for epoch in range(1, int(args.inr_epochs) + 1):
            model.eval()
            inr.train()
            values = []
            for slide in train_slides:
                features = torch.from_numpy(
                    np.asarray(slide.features, dtype=np.float32)
                ).to(device)
                graphs = _to_device_graph_pair(slide.graphs, device)
                coordinates = torch.from_numpy(
                    normalize_slide_coordinates(slide.coordinates)
                ).to(device)
                with torch.no_grad():
                    rna_logits, composition_logits, context, auxiliary = model(
                        features,
                        graphs,
                        support_gamma=float(final_weights["support_gamma"]),
                        domain_strength=0.0,
                    )
                inr_optimizer.zero_grad(set_to_none=True)
                correction = inr(coordinates, context)
                rna_correction, composition_correction = correction.split(programs, dim=1)
                local_expression = (
                    torch.from_numpy(slide.local_theta_rna).to(device)
                    @ torch.from_numpy(slide.local_basis_panel).to(device)
                )
                slide_consensus_basis_panel = torch.from_numpy(
                    slide.consensus_basis_panel
                ).to(device)
                loss, _ = _loss(
                    rna_logits=rna_logits + rna_correction,
                    composition_logits=composition_logits + composition_correction,
                    auxiliary=auxiliary,
                    rna_target=torch.from_numpy(slide.theta_rna).to(device),
                    composition_target=torch.from_numpy(slide.theta_composition).to(device),
                    local_expression_target=local_expression,
                    consensus_basis_panel=slide_consensus_basis_panel,
                    gene_weight=gene_weight,
                    program_weight=program_weight,
                    rna_per_cell=rna_per_cell,
                    program_group_index=model.program_group_index,
                    graphs=graphs,
                    type_lift_alpha=float(args.type_lift_alpha),
                    source_target=torch.tensor(slide.source_index, device=device),
                    technology_target=torch.tensor(slide.technology_index, device=device),
                    support_mass=float(args.support_mass),
                    weights=final_weights,
                )
                loss = loss * float(source_weight[slide.source])
                loss.backward()
                torch.nn.utils.clip_grad_norm_(inr.parameters(), 5.0)
                inr_optimizer.step()
                values.append(float(loss.detach().cpu()))
            validation_loss, validation, _ = _evaluate(
                model,
                inr,
                slides,
                split="val",
                consensus_basis_panel=consensus_basis_panel,
                gene_weight=gene_weight,
                program_weight=program_weight,
                rna_per_cell=rna_per_cell,
                device=device,
                weights=final_weights,
            )
            validation_pcc = float(
                validation["patient_source_balanced_loss_panel_gene_pcc_median"]
            )
            record = {
                "epoch": float(epoch),
                "train_loss": float(np.mean(values)),
                "validation_loss": validation_loss,
                "validation_gene_pcc": validation_pcc,
                "gate": float(torch.sigmoid(inr.gate_logit).detach().cpu()),
            }
            inr_history.append(record)
            print(json.dumps({"inr": record}), flush=True)
            if validation_pcc > best_inr_pcc + 1.0e-4 and validation_loss <= base_val_loss * 1.01:
                best_inr_pcc = validation_pcc
                best_inr_loss = validation_loss
                best_inr_state = copy.deepcopy(inr.state_dict())
        inr.load_state_dict(best_inr_state)
        inr_accepted = bool(best_inr_pcc > base_val_pcc + 1.0e-4)
        selected_inr = inr if inr_accepted else None

    _, validation_final, validation_predictions = _evaluate(
        model,
        selected_inr,
        slides,
        split="val",
        consensus_basis_panel=consensus_basis_panel,
        gene_weight=gene_weight,
        program_weight=program_weight,
        rna_per_cell=rna_per_cell,
        device=device,
        weights=final_weights,
    )
    _, test_final, test_predictions = _evaluate(
        model,
        selected_inr,
        slides,
        split="test",
        consensus_basis_panel=consensus_basis_panel,
        gene_weight=gene_weight,
        program_weight=program_weight,
        rna_per_cell=rna_per_cell,
        device=device,
        weights=final_weights,
    )
    validation_full, validation_gene_table = _full_10k_evaluation(
        slides=slides,
        split="val",
        predictions=validation_predictions,
        consensus_reference=consensus_reference_np,
        genes=genes,
        block_size=int(args.full_evaluation_block),
    )
    test_full, test_gene_table = _full_10k_evaluation(
        slides=slides,
        split="test",
        predictions=test_predictions,
        consensus_reference=consensus_reference_np,
        genes=genes,
        block_size=int(args.full_evaluation_block),
    )
    validation_gene_table.to_csv(
        args.output_dir / "validation_gene_pcc_10k.csv", index=False
    )
    test_gene_table.to_csv(args.output_dir / "test_gene_pcc_10k.csv", index=False)
    pd.DataFrame(history).to_csv(args.output_dir / "training_history.csv", index=False)
    pd.DataFrame(inr_history).to_csv(args.output_dir / "inr_history.csv", index=False)
    program_diagnostics = pd.DataFrame(
        {
            "program_index": np.arange(programs, dtype=np.int64),
            "program_group": program_groups_np,
            "loss_weight": program_weight_np,
            "rna_local_graph_gate": torch.sigmoid(model.rna_local_refine)
            .detach()
            .cpu()
            .numpy(),
            "rna_context_graph_gate": torch.sigmoid(model.rna_context_refine)
            .detach()
            .cpu()
            .numpy(),
            "composition_local_graph_gate": torch.sigmoid(
                model.composition_local_refine
            )
            .detach()
            .cpu()
            .numpy(),
            "composition_context_graph_gate": torch.sigmoid(
                model.composition_context_refine
            )
            .detach()
            .cpu()
            .numpy(),
        }
    )
    program_diagnostics.to_csv(
        args.output_dir / "learned_program_graph_gates.csv", index=False
    )

    test_slides = [slide for slide in slides if slide.split == "test"]
    test_spots = pd.concat(
        [
            slide.spot_index.assign(
                image_model_slide_key=slide.key,
                image_model_split="test",
                reference_row=slide.reference_rows,
            )
            for slide in test_slides
        ],
        ignore_index=True,
    )
    np.save(
        args.output_dir / "theta_rna_predicted_test.npy",
        np.concatenate([test_predictions[slide.key][0] for slide in test_slides], axis=0),
    )
    np.save(
        args.output_dir / "theta_composition_predicted_test.npy",
        np.concatenate([test_predictions[slide.key][1] for slide in test_slides], axis=0),
    )
    test_spots.to_csv(args.output_dir / "test_spot_index.csv", index=False)
    torch.save(
        {
            "model": model.state_dict(),
            "inr": selected_inr.state_dict() if selected_inr is not None else None,
            "inr_accepted": inr_accepted,
            "input_dim": input_dim,
            "programs": programs,
            "hidden": int(args.hidden),
            "graph_layers": int(args.graph_layers),
            "context_graph_layers": int(args.context_graph_layers),
            "graph_heads": int(args.graph_heads),
            "context_neighbors": int(args.context_neighbors),
            "program_groups": int(program_groups_np.max()) + 1,
            "type_lift_alpha": float(args.type_lift_alpha),
            "sources": len(source_levels),
            "technologies": len(technology_levels),
            "source_levels": source_levels,
            "technology_levels": technology_levels,
            "support_mass": float(args.support_mass),
            "hierarchical_logits": bool(args.hierarchical_logits),
            "program_group_method": str(args.program_group_method),
            "support_gamma": float(final_weights["support_gamma"]),
            "neighbors": int(args.neighbors),
            "consensus_reference_sha256": _sha256(
                args.consensus_dir / "consensus_reference_probability_10k.npy"
            ),
        },
        args.output_dir / "model_checkpoint.pt",
    )
    payload = {
        "method": f"virchow2_pan_atlas_g{programs}_full10k_mapper",
        "target_programs": programs,
        "target_rna": "theta_rna_consensus.npy",
        "target_composition": "theta_composition_consensus.npy",
        "expression_decoder": "predicted_theta_rna @ slide_matched_batch_adapted_consensus_W",
        "pcc_reference": "original_per_slide_top10k_theta_rna @ basis_probability",
        "loss_updates": {
            "dual_composition_heads": True,
            "bounded_jsd_for_both_theta_targets": True,
            "support_weighted_theta_pcc": True,
            "logcp10k_huber_against_original_local_inverse": True,
            "gene_pcc_against_original_local_inverse": True,
            "pcc_warmup_and_ramp": True,
            "bulk_and_spatial_gradient_matching": True,
            "rna_composition_consistency": True,
            "entropy_matching": True,
            "full_10k_expression_log_loss": int(len(panel_np)) == 10000,
            "full_10k_gene_pcc_loss": int(len(panel_np)) == 10000,
            "dual_scale_morphology_graph": True,
            "shared_sparse_program_support": True,
            "conditional_rna_and_composition_abundance": True,
            "raw_plus_slide_centered_virchow2": True,
            "train_only_morphology_program_groups": str(
                args.program_group_method
            ) == "morphology",
            "hierarchical_group_and_type_logits": bool(args.hierarchical_logits),
            "learned_edge_boundary_and_homogeneous_channels": True,
            "source_technology_gradient_reversal": True,
            "train_only_spatial_reliability_gene_weight": not bool(
                args.disable_spatial_gene_weighting
            ),
        },
        "sources": list(args.sources),
        "source_balance_power": float(args.source_balance_power),
        "source_training_weights": source_weight,
        "train_slides": sum(slide.split == "train" for slide in slides),
        "validation_slides": sum(slide.split == "val" for slide in slides),
        "test_slides": sum(slide.split == "test" for slide in slides),
        "expression_loss_genes": int(len(panel_np)),
        "gene_pcc_loss_genes": int(len(panel_np)),
        "program_groups": int(program_groups_np.max()) + 1,
        "training_configuration": {
            key: str(value) if isinstance(value, Path) else value
            for key, value in vars(args).items()
        },
        "learned_architecture_gates": {
            "centered_feature": float(
                torch.sigmoid(model.centered_feature_gate).detach().cpu()
            ),
            "support_prototype": float(
                torch.sigmoid(model.support_type_gate).detach().cpu()
            ),
            "support_temperature": float(
                model.log_support_temperature.exp().clamp(0.05, 1.0).detach().cpu()
            ),
            "local_graph_trunk": float(torch.sigmoid(model.local_graph_gate).detach().cpu()),
            "context_graph_trunk": float(
                torch.sigmoid(model.context_graph_gate).detach().cpu()
            ),
            "rna_type": float(torch.sigmoid(model.rna_type_gate).detach().cpu()),
            "composition_type": float(
                torch.sigmoid(model.composition_type_gate).detach().cpu()
            ),
            "rna_group": float(torch.sigmoid(model.rna_group_gate).detach().cpu()),
            "composition_group": float(
                torch.sigmoid(model.composition_group_gate).detach().cpu()
            ),
            "type_temperature": float(
                model.log_type_temperature.exp().clamp(0.05, 1.0).detach().cpu()
            ),
        },
        "history": history,
        "inr_history": inr_history,
        "inr_accepted": inr_accepted,
        "validation_loss_panel": validation_final,
        "test_loss_panel": test_final,
        "validation_full_10k": validation_full,
        "test_full_10k": test_full,
        "runtime_seconds": time.time() - started,
        "peak_cuda_memory_bytes": int(torch.cuda.max_memory_allocated()),
    }
    (args.output_dir / "metrics.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
