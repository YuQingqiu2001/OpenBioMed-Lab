"""Frozen mathematical definitions used by CoVarST; historical CLI omitted."""
from __future__ import annotations

import argparse

import copy

import json

import math

import time

from dataclasses import dataclass

from pathlib import Path

import numpy as np

import pandas as pd

import torch

from scipy.spatial import cKDTree

from torch import nn

from torch.nn import functional as F

from openst_final.he_uni_mamba_bnn_inr import SineLayer, normalize_slide_coordinates

EPS = 1.0e-8

@dataclass
class SlideData:
    key: str
    sample_id: str
    source: str
    patient: str
    split: str
    reference_rows: np.ndarray
    features: np.ndarray
    fseg_histogram: np.ndarray
    coordinates: np.ndarray
    theta: np.ndarray
    target: np.ndarray
    presence: np.ndarray
    spot_index: pd.DataFrame
    graph: tuple[np.ndarray, np.ndarray, np.ndarray]

class FSegGATv2(nn.Module):
    def __init__(self, hidden: int, heads: int, dropout: float) -> None:
        super().__init__()
        if hidden % heads:
            raise ValueError("hidden must be divisible by heads")
        self.heads = int(heads)
        self.head_dim = int(hidden // heads)
        self.norm = nn.LayerNorm(hidden)
        self.left = nn.Linear(hidden, hidden, bias=False)
        self.right = nn.Linear(hidden, hidden, bias=False)
        self.message = nn.Linear(hidden, hidden, bias=False)
        self.attention = nn.Parameter(torch.empty(heads, self.head_dim))
        self.output = nn.Linear(hidden, hidden, bias=False)
        self.dropout = nn.Dropout(dropout)
        nn.init.xavier_uniform_(self.attention)

    def forward(
        self,
        value: torch.Tensor,
        graph: tuple[torch.Tensor, torch.Tensor, torch.Tensor],
    ) -> torch.Tensor:
        source, target, prior = graph
        n = value.shape[0]
        x = self.norm(value)
        left = self.left(x).reshape(n, self.heads, self.head_dim)
        right = self.right(x).reshape(n, self.heads, self.head_dim)
        message = self.message(x).reshape(n, self.heads, self.head_dim)
        dynamic = F.leaky_relu(left[source] + right[target], 0.2)
        score = (dynamic * self.attention[None, :, :]).sum(dim=2)
        score = score.float() + torch.log(prior.clamp_min(1.0e-6))[:, None]
        maximum = torch.full(
            (n, self.heads),
            -torch.inf,
            dtype=score.dtype,
            device=score.device,
        )
        maximum.scatter_reduce_(
            0,
            source[:, None].expand(-1, self.heads),
            score,
            reduce="amax",
            include_self=True,
        )
        attention = torch.exp(score - maximum[source])
        denominator = torch.zeros(
            (n, self.heads), dtype=attention.dtype, device=score.device
        )
        denominator.index_add_(0, source, attention)
        # Keep this normalization out-of-place: ``attention`` is also an input
        # to ``denominator.index_add_`` and autograd needs its original version
        # when differentiating the normalization.
        attention = attention / denominator[source].clamp_min(EPS)
        aggregate = torch.zeros(
            (n, self.heads, self.head_dim),
            dtype=message.dtype,
            device=value.device,
        )
        aggregate.index_add_(
            0,
            source,
            attention.to(message.dtype)[:, :, None] * message[target],
        )
        return value + self.dropout(self.output(aggregate.reshape(n, -1)))

class SlideContextK64(nn.Module):
    def __init__(
        self,
        input_dim: int,
        *,
        hidden: int,
        k: int,
        graph_layers: int,
        graph_heads: int,
        dropout: float,
    ) -> None:
        super().__init__()
        self.encoder = nn.Sequential(
            nn.LayerNorm(input_dim),
            nn.Linear(input_dim, hidden * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden * 2, hidden),
            nn.LayerNorm(hidden),
        )
        self.graph = nn.ModuleList(
            [
                FSegGATv2(hidden, graph_heads, dropout)
                for _ in range(graph_layers)
            ]
        )
        self.graph_gate = nn.Parameter(torch.tensor(-2.0))
        self.slide_encoder = nn.Sequential(
            nn.Linear(hidden * 2, hidden),
            nn.GELU(),
            nn.LayerNorm(hidden),
        )
        self.local_head = nn.Sequential(
            nn.Linear(hidden * 3, hidden * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden * 2, k),
        )
        self.slide_head = nn.Linear(hidden, k)

    def forward(
        self,
        features: torch.Tensor,
        graph: tuple[torch.Tensor, torch.Tensor, torch.Tensor],
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        local = self.encoder(features)
        spatial = local
        for block in self.graph:
            spatial = block(spatial, graph)
        gate = torch.sigmoid(self.graph_gate)
        spatial = local + gate * (spatial - local)
        pooled = torch.cat(
            [
                spatial.mean(dim=0),
                spatial.std(dim=0, unbiased=False),
            ],
            dim=0,
        )
        slide = self.slide_encoder(pooled)
        repeated = slide[None, :].expand(len(spatial), -1)
        logits = self.local_head(
            torch.cat([spatial, repeated, spatial * repeated], dim=1)
        )
        logits = logits + self.slide_head(slide)[None, :]
        logits = logits - logits.mean(dim=1, keepdim=True)
        return logits, spatial, slide

class BoundedINR(nn.Module):
    def __init__(self, hidden: int, k: int) -> None:
        super().__init__()
        self.coordinate = nn.Sequential(
            SineLayer(2, 64, 8.0),
            SineLayer(64, 64, 1.0),
            SineLayer(64, 64, 1.0),
        )
        self.context = nn.Sequential(
            nn.LayerNorm(hidden),
            nn.Linear(hidden, 64),
            nn.GELU(),
        )
        self.output = nn.Linear(128, k)
        self.gate_logit = nn.Parameter(torch.tensor(-2.0))
        nn.init.zeros_(self.output.weight)
        nn.init.zeros_(self.output.bias)

    def forward(
        self, coordinates: torch.Tensor, context: torch.Tensor
    ) -> torch.Tensor:
        value = self.output(
            torch.cat(
                [self.coordinate(coordinates), self.context(context)], dim=1
            )
        )
        value = 0.10 * torch.tanh(value)
        value = value - value.mean(dim=0, keepdim=True)
        return torch.sigmoid(self.gate_logit) * value

def _to_device_graph(
    graph: tuple[np.ndarray, np.ndarray, np.ndarray],
    device: torch.device,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    return (
        torch.from_numpy(graph[0]).to(device=device, dtype=torch.long),
        torch.from_numpy(graph[1]).to(device=device, dtype=torch.long),
        torch.from_numpy(graph[2]).to(device=device, dtype=torch.float32),
    )

def _loss(
    *,
    logits: torch.Tensor,
    theta_target: torch.Tensor,
    expression_target: torch.Tensor,
    basis: torch.Tensor,
    presence: torch.Tensor,
    panel: torch.Tensor,
    graph: tuple[torch.Tensor, torch.Tensor, torch.Tensor],
    expression_weight: float,
    gene_pcc_weight: float,
    theta_weight: float,
    bulk_weight: float,
    gradient_weight: float,
) -> tuple[torch.Tensor, dict[str, float]]:
    theta = torch.softmax(logits.float(), dim=1)
    expression = theta @ basis
    expression = expression * presence[None, :]
    expression = expression / expression.sum(dim=1, keepdim=True).clamp_min(EPS)
    target = expression_target.clamp_min(0.0)
    target = target / target.sum(dim=1, keepdim=True).clamp_min(EPS)

    expression_ce = -(
        target * torch.log(expression.clamp_min(1.0e-8))
    ).sum(dim=1).mean()
    theta_hellinger = 0.5 * (
        torch.sqrt(theta.clamp_min(EPS))
        - torch.sqrt(theta_target.clamp_min(EPS))
    ).square().sum(dim=1).mean()

    target_log = torch.log1p(1.0e4 * target[:, panel])
    predicted_log = torch.log1p(1.0e4 * expression[:, panel])
    target_centered = target_log - target_log.mean(dim=0, keepdim=True)
    predicted_centered = predicted_log - predicted_log.mean(
        dim=0, keepdim=True
    )
    # A raw ``sqrt(0)`` has an infinite derivative.  Even when zero-variance
    # genes are masked below, the backward pass can otherwise produce 0 * inf
    # and poison the first optimizer step with NaNs.
    target_sum_square = target_centered.square().sum(dim=0)
    predicted_sum_square = predicted_centered.square().sum(dim=0)
    target_norm = torch.sqrt(target_sum_square.clamp_min(1.0e-12))
    predicted_norm = torch.sqrt(predicted_sum_square.clamp_min(1.0e-12))
    valid = (target_norm > 1.0e-3) & (predicted_norm > 1.0e-5)
    if bool(valid.any()):
        correlation = (
            target_centered[:, valid] * predicted_centered[:, valid]
        ).sum(dim=0) / (
            target_norm[valid] * predicted_norm[valid]
        ).clamp_min(EPS)
        gene_pcc = 1.0 - correlation.clamp(-1.0, 1.0).mean()
    else:
        gene_pcc = expression_ce.new_tensor(1.0)

    bulk = 0.5 * (
        torch.sqrt(theta.mean(dim=0).clamp_min(EPS))
        - torch.sqrt(theta_target.mean(dim=0).clamp_min(EPS))
    ).square().sum()

    source, target_edge, weight = graph
    nonself = source != target_edge
    if bool(nonself.any()):
        predicted_gradient = theta[source[nonself]] - theta[target_edge[nonself]]
        target_gradient = (
            theta_target[source[nonself]] - theta_target[target_edge[nonself]]
        )
        gradient_per_edge = F.smooth_l1_loss(
            predicted_gradient,
            target_gradient,
            reduction="none",
        ).mean(dim=1)
        gradient = (
            gradient_per_edge * weight[nonself]
        ).sum() / weight[nonself].sum().clamp_min(EPS)
    else:
        gradient = expression_ce.new_zeros(())

    total = (
        expression_weight * expression_ce
        + gene_pcc_weight * gene_pcc
        + theta_weight * theta_hellinger
        + bulk_weight * bulk
        + gradient_weight * gradient
    )
    return total, {
        "expression_ce": float(expression_ce.detach().cpu()),
        "gene_pcc_loss": float(gene_pcc.detach().cpu()),
        "theta_hellinger": float(theta_hellinger.detach().cpu()),
        "bulk_hellinger": float(bulk.detach().cpu()),
        "gradient": float(gradient.detach().cpu()),
    }

@torch.no_grad()
def _evaluate_split(
    model: SlideContextK64,
    inr: BoundedINR | None,
    slides: list[SlideData],
    *,
    split: str,
    basis: torch.Tensor,
    panel: torch.Tensor,
    device: torch.device,
    weights: dict[str, float],
) -> tuple[float, dict[str, object], list[np.ndarray]]:
    model.eval()
    if inr is not None:
        inr.eval()
    losses: list[float] = []
    records: list[dict[str, object]] = []
    pcc_records: list[tuple[str, str, np.ndarray]] = []
    theta_outputs: list[np.ndarray] = []
    for slide in slides:
        if slide.split != split:
            continue
        feature = torch.from_numpy(
            np.asarray(slide.features, dtype=np.float32)
        ).to(device)
        graph = _to_device_graph(slide.graph, device)
        logits, context, _ = model(feature, graph)
        correction = None
        if inr is not None:
            coordinates = torch.from_numpy(
                normalize_slide_coordinates(slide.coordinates)
            ).to(device)
            correction = inr(coordinates, context)
            logits = logits + correction
        theta_target = torch.from_numpy(slide.theta).to(device)
        expression_target = torch.from_numpy(slide.target).to(device)
        presence = torch.from_numpy(slide.presence).to(device)
        loss, parts = _loss(
            logits=logits,
            theta_target=theta_target,
            expression_target=expression_target,
            basis=basis,
            presence=presence,
            panel=panel,
            graph=graph,
            **weights,
        )
        predicted_theta = torch.softmax(logits, dim=1)
        predicted_expression = predicted_theta @ basis
        predicted_expression *= presence[None, :]
        predicted_expression /= predicted_expression.sum(
            dim=1, keepdim=True
        ).clamp_min(EPS)
        target_log = torch.log1p(1.0e4 * expression_target[:, panel])
        predicted_log = torch.log1p(1.0e4 * predicted_expression[:, panel])
        x = target_log - target_log.mean(dim=0, keepdim=True)
        y = predicted_log - predicted_log.mean(dim=0, keepdim=True)
        denominator = torch.sqrt(x.square().sum(0) * y.square().sum(0))
        valid = denominator > 1.0e-12
        pcc = (
            (x[:, valid] * y[:, valid]).sum(0) / denominator[valid]
        ).cpu().numpy()
        pcc_full = np.full(int(panel.numel()), np.nan, dtype=np.float32)
        pcc_full[valid.cpu().numpy()] = pcc
        losses.append(float(loss.cpu()))
        record = {
            "slide_key": slide.key,
            "source": slide.source,
            "patient": slide.patient,
            "spots": len(slide.features),
            "loss": float(loss.cpu()),
            "panel_gene_pcc_median": float(np.median(pcc)),
            "panel_gene_pcc_mean": float(np.mean(pcc)),
            **parts,
            "inr_mean_abs_correction": (
                float(correction.abs().mean().cpu())
                if correction is not None
                else 0.0
            ),
        }
        records.append(record)
        pcc_records.append((slide.source, slide.patient, pcc_full))
        theta_outputs.append(predicted_theta.cpu().numpy().astype(np.float32))
    by_source = {}
    for source in sorted({record["source"] for record in records}):
        selected = [
            record["panel_gene_pcc_median"]
            for record in records
            if record["source"] == source
        ]
        by_source[source] = float(np.mean(selected))
    patient_vectors: list[np.ndarray] = []
    by_source_patient_balanced: dict[str, float] = {}
    for source in sorted({item[0] for item in pcc_records}):
        source_patients: list[np.ndarray] = []
        for patient in sorted(
            {item[1] for item in pcc_records if item[0] == source}
        ):
            source_patients.append(
                np.nanmean(
                    np.stack(
                        [
                            value
                            for current_source, current_patient, value
                            in pcc_records
                            if current_source == source
                            and current_patient == patient
                        ],
                        axis=0,
                    ),
                    axis=0,
                )
            )
        source_balanced = np.nanmean(
            np.stack(source_patients, axis=0),
            axis=0,
        )
        patient_vectors.append(source_balanced)
        by_source_patient_balanced[source] = float(
            np.nanmedian(source_balanced)
        )
    patient_balanced = np.nanmean(
        np.stack(patient_vectors, axis=0),
        axis=0,
    )
    summary = {
        "split": split,
        "slides": records,
        "mean_loss": float(np.mean(losses)),
        "source_balanced_slide_median_gene_pcc": float(
            np.mean(list(by_source.values()))
        ),
        "by_source_slide_median_gene_pcc": by_source,
        "patient_balanced_gene_pcc_median": float(
            np.nanmedian(patient_balanced)
        ),
        "by_source_patient_balanced_gene_pcc_median": (
            by_source_patient_balanced
        ),
    }
    return float(np.mean(losses)), summary, theta_outputs
