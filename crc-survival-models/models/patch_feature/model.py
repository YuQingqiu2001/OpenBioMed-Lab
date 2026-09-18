"""Portable five-fold inference for the patch-feature CRC survival model."""

from __future__ import annotations

import json
import math
import pickle
from pathlib import Path

import h5py
import numpy as np
import torch
import torch.nn as nn


PHYSICAL_SCALES_UM = (250.0, 500.0, 1000.0)
MAXIMUM_NEIGHBOURS = 12
MAX_PATIENT_TOKENS = 768
MAX_SLIDE_TOKENS = 192


def _torch_load(path: Path):
    try:
        return torch.load(path, map_location="cpu", weights_only=False)
    except TypeError:
        return torch.load(path, map_location="cpu")


class PatchEncoder(nn.Module):
    def __init__(self, image_dim: int, topology_dim: int, guided_dim: int, free_dim: int):
        super().__init__()
        hidden = 192
        self.shared = nn.Sequential(nn.Linear(image_dim, hidden), nn.LayerNorm(hidden), nn.GELU(), nn.Dropout(0.10))
        self.guided = nn.Sequential(nn.Linear(hidden, guided_dim), nn.LayerNorm(guided_dim), nn.GELU())
        self.free = nn.Sequential(nn.Linear(hidden, free_dim), nn.LayerNorm(free_dim), nn.GELU())
        self.topology_head = nn.Sequential(nn.Linear(guided_dim, 96), nn.GELU(), nn.Dropout(0.05), nn.Linear(96, topology_dim))
        self.free_reconstruction_head = nn.Sequential(nn.Linear(free_dim, hidden), nn.GELU(), nn.Linear(hidden, image_dim))

    def forward(self, image: torch.Tensor):
        shared = self.shared(image)
        guided = self.guided(shared)
        free = self.free(shared)
        return guided, free, self.topology_head(guided), self.free_reconstruction_head(free)


class SpatialField(nn.Module):
    def __init__(self, input_dim: int, image_dim: int, topology_dim: int, hidden: int, free_dim: int):
        super().__init__()
        self.field = nn.Sequential(
            nn.Linear(input_dim, hidden), nn.LayerNorm(hidden), nn.GELU(), nn.Dropout(0.12),
            nn.Linear(hidden, hidden), nn.LayerNorm(hidden), nn.GELU(), nn.Dropout(0.08),
        )
        self.topology_mu = nn.Linear(hidden, topology_dim)
        self.topology_logvar = nn.Linear(hidden, topology_dim)
        self.free = nn.Sequential(nn.Linear(image_dim, 96), nn.LayerNorm(96), nn.GELU(), nn.Dropout(0.10), nn.Linear(96, free_dim))
        self.reconstruct = nn.Sequential(nn.GELU(), nn.Linear(free_dim, image_dim))

    def forward(self, field_x: torch.Tensor, image_x: torch.Tensor):
        hidden = self.field(field_x)
        mu = self.topology_mu(hidden)
        logvar = self.topology_logvar(hidden).clamp(-4.0, 3.0)
        free = self.free(image_x)
        return mu, logvar, free, self.reconstruct(free)


class SurvivalAggregator(nn.Module):
    def __init__(self, input_dim: int, hidden: int = 64):
        super().__init__()
        self.project = nn.Sequential(nn.Linear(input_dim, hidden), nn.LayerNorm(hidden), nn.GELU(), nn.Dropout(0.18))
        self.attention_v = nn.Linear(hidden, hidden // 2)
        self.attention_u = nn.Linear(hidden, hidden // 2)
        self.attention_w = nn.Linear(hidden // 2, 1, bias=False)
        self.instance_risk = nn.Linear(hidden, 1)
        slide_dim = hidden * 3 + 6
        self.slide_project = nn.Sequential(nn.LayerNorm(slide_dim), nn.Linear(slide_dim, hidden), nn.GELU(), nn.Dropout(0.16))
        self.slide_attention = nn.Linear(hidden, 1, bias=False)
        self.slide_risk = nn.Linear(hidden, 1)
        self.patient_risk = nn.Sequential(nn.LayerNorm(hidden + 3), nn.Dropout(0.15), nn.Linear(hidden + 3, 1))
        self.horizon_head = nn.Sequential(nn.LayerNorm(hidden), nn.Linear(hidden, 3))

    def forward(self, token, compartment_weight, valid_mask, slide_index, max_slides):
        hidden = self.project(token)
        base_attention = self.attention_w(torch.tanh(self.attention_v(hidden)) * torch.sigmoid(self.attention_u(hidden))).squeeze(-1)
        instance = self.instance_risk(hidden).squeeze(-1)
        pooled_by_compartment, hotspot_by_compartment, coverage_by_compartment = [], [], []
        slide_masks = [valid_mask & (slide_index == value) for value in range(max_slides)]
        slide_valid = torch.stack([mask.any(dim=1) for mask in slide_masks], dim=1)
        for compartment in range(3):
            pooled_slides, hotspot_slides, coverage_slides = [], [], []
            for mask in slide_masks:
                weight = compartment_weight[:, :, compartment].clamp_min(1e-6)
                logits = (base_attention + torch.log(weight)).masked_fill(~mask, -1e4)
                attention = torch.softmax(logits, dim=1) * mask
                attention = attention / attention.sum(dim=1, keepdim=True).clamp_min(1e-8)
                pooled_slides.append((attention[:, :, None] * hidden).sum(dim=1))
                count = mask.sum(dim=1).clamp_min(1).float()
                weighted_instance = (instance + 0.20 * torch.log(weight)).masked_fill(~mask, -1e4)
                hotspot = 0.20 * (torch.logsumexp(weighted_instance / 0.20, dim=1) - torch.log(count))
                hotspot_slides.append(torch.where(mask.any(dim=1), hotspot, torch.zeros_like(hotspot)))
                coverage_slides.append((weight * mask).sum(dim=1) / count)
            pooled_by_compartment.append(torch.stack(pooled_slides, dim=1))
            hotspot_by_compartment.append(torch.stack(hotspot_slides, dim=1))
            coverage_by_compartment.append(torch.stack(coverage_slides, dim=1))
        slide_summary = torch.cat(
            [*pooled_by_compartment, torch.stack(hotspot_by_compartment, dim=2), torch.stack(coverage_by_compartment, dim=2)], dim=2
        )
        slide_hidden = self.slide_project(slide_summary)
        slide_risk = self.slide_risk(slide_hidden).squeeze(-1).masked_fill(~slide_valid, -1e4)
        slide_logits = self.slide_attention(slide_hidden).squeeze(-1).masked_fill(~slide_valid, -1e4)
        slide_attention = torch.softmax(slide_logits, dim=1) * slide_valid
        slide_attention = slide_attention / slide_attention.sum(dim=1, keepdim=True).clamp_min(1e-8)
        patient_hidden = (slide_attention[:, :, None] * slide_hidden).sum(dim=1)
        mean_risk = (slide_attention * slide_risk.masked_fill(~slide_valid, 0.0)).sum(dim=1)
        variance = (slide_attention * (slide_risk.masked_fill(~slide_valid, 0.0) - mean_risk[:, None]).square()).sum(dim=1)
        count = slide_valid.sum(dim=1).clamp_min(1).float()
        worst_risk = 0.25 * (torch.logsumexp(slide_risk / 0.25, dim=1) - torch.log(count))
        patient_summary = torch.cat([patient_hidden, mean_risk[:, None], worst_risk[:, None], variance.clamp_min(1e-8).sqrt()[:, None]], dim=1)
        return self.patient_risk(patient_summary).squeeze(1), self.horizon_head(patient_hidden)


def _spatial_design(image_x: np.ndarray, probability: np.ndarray, xy_um: np.ndarray):
    n = len(xy_um)
    if n == 1:
        gradients = [np.zeros_like(image_x) for _ in PHYSICAL_SCALES_UM]
        densities = np.zeros((1, len(gradients)), np.float32)
        local_probability = probability.copy()
    else:
        difference = xy_um[:, None, :] - xy_um[None, :, :]
        distance2 = np.square(difference).sum(axis=2)
        np.fill_diagonal(distance2, np.inf)
        k = min(MAXIMUM_NEIGHBOURS, n - 1)
        neighbour = np.argpartition(distance2, kth=k - 1, axis=1)[:, :k]
        row = np.arange(n)[:, None]
        neighbour_distance2 = distance2[row, neighbour]
        gradients, density_columns, contexts = [], [], []
        for scale in PHYSICAL_SCALES_UM:
            weight = np.exp(-0.5 * neighbour_distance2 / scale**2).astype(np.float32)
            weight *= neighbour_distance2 <= (3.0 * scale) ** 2
            denominator = weight.sum(axis=1, keepdims=True)
            empty = denominator[:, 0] <= 1e-8
            denominator[empty] = 1.0
            context = (weight[..., None] * image_x[neighbour]).sum(axis=1) / denominator
            p_context = (weight[..., None] * probability[neighbour]).sum(axis=1) / denominator
            context[empty], p_context[empty] = image_x[empty], probability[empty]
            gradients.append(context - image_x)
            contexts.append(p_context)
            density_columns.append(np.log1p(weight.sum(axis=1)) / math.log1p(max(k, 1)))
        densities = np.stack(density_columns, axis=1).astype(np.float32)
        local_probability = contexts[min(1, len(contexts) - 1)]
    entropy = -(probability * np.log(np.clip(probability, 1e-6, 1.0))).sum(axis=1, keepdims=True) / math.log(probability.shape[1])
    design = np.concatenate([image_x, *gradients, probability, local_probability - probability, densities, entropy.astype(np.float32)], axis=1)
    return design.astype(np.float32), local_probability.astype(np.float32)


def _compartment_weights(slide: dict) -> dict[str, np.ndarray]:
    probability, local = slide["probability"], slide["local_probability"]
    p_tum, local_tum = np.clip(probability[:, 8], 0.0, 1.0), np.clip(local[:, 8], 0.0, 1.0)
    confidence = np.exp(-0.35 * np.clip(slide["uncertainty"], 0.0, 8.0))
    result = {
        "core": p_tum**2 * confidence,
        "front": np.clip(0.65 * (p_tum * (1.0 - local_tum) + (1.0 - p_tum) * local_tum) + 0.35 * 4.0 * p_tum * (1.0 - p_tum), 0.0, 1.0) * confidence,
        "peri": np.clip((1.0 - p_tum) * local_tum, 0.0, 1.0) * confidence,
    }
    for key, weight in result.items():
        if float(weight.sum()) <= 1e-6:
            result[key] = np.ones_like(weight) * confidence
    return result


def _encode_cache(path: Path, transforms, field, point, device):
    with h5py.File(path, "r") as handle:
        if not bool(handle.attrs.get("complete", 0)):
            raise ValueError(f"Incomplete patch cache: {path}")
        embedding = handle["uni/embedding"][:]
        delta = handle["uni/neighbor_delta_center"][:]
        probability = handle["uni/coarse_probability"][:, 0].astype(np.float32)
        xy = handle["regions/xy_um"][:].astype(np.float32)
        tissue = handle["regions/rgb_tissue_fraction"][:].astype(np.float32)
    valid = np.isfinite(embedding).all(axis=(1, 2)) & np.isfinite(delta).all(axis=1) & np.isfinite(probability).all(axis=1) & (tissue >= 0.03)
    if not valid.any():
        raise ValueError(f"No valid tissue regions in {path}")
    image = np.concatenate([embedding[valid, 0], delta[valid]], axis=1).astype(np.float32)
    image_x = transforms["image_pca"].transform(image).astype(np.float32)
    probability = probability[valid]
    design, local_probability = _spatial_design(image_x, probability, xy[valid])
    mu_parts, logvar_parts, free_parts, point_topology_parts, point_free_parts = [], [], [], [], []
    for start in range(0, len(image_x), 2048):
        field_x = torch.from_numpy(transforms["field_scaler"].transform(design[start : start + 2048]).astype(np.float32)).to(device)
        image_t = torch.from_numpy(image_x[start : start + 2048]).to(device)
        with torch.inference_mode():
            mu, logvar, free, _ = field(field_x, image_t)
            _, point_free, point_topology, _ = point(image_t)
        mu_parts.append(mu.cpu().numpy()); logvar_parts.append(logvar.cpu().numpy()); free_parts.append(free.cpu().numpy())
        point_topology_parts.append(point_topology.cpu().numpy()); point_free_parts.append(point_free.cpu().numpy())
    topology = np.concatenate(mu_parts).astype(np.float32)
    point_topology = np.concatenate(point_topology_parts).astype(np.float32)
    return {
        "topology": topology,
        "free": np.concatenate(free_parts).astype(np.float32),
        "uncertainty": np.exp(np.concatenate(logvar_parts)).mean(axis=1).astype(np.float32),
        "probability": probability,
        "local_probability": local_probability,
        "point_topology": point_topology,
        "point_free": np.concatenate(point_free_parts).astype(np.float32),
        "field_residual": topology - point_topology,
    }


def _make_token(slide: dict) -> np.ndarray:
    uncertainty = np.log1p(np.clip(slide["uncertainty"], 0.0, 20.0))[:, None]
    return np.concatenate([slide["point_topology"], slide["topology"], slide["field_residual"], slide["point_free"], slide["probability"], uncertainty], axis=1).astype(np.float32)


def _select(slide: dict, cap: int) -> np.ndarray:
    n = len(slide["topology"])
    if n <= cap:
        return np.arange(n, dtype=np.int64)
    current = _compartment_weights(slide)
    quotas = {"front": int(round(cap * 0.50)), "core": int(round(cap * 0.25))}
    quotas["peri"] = cap - quotas["front"] - quotas["core"]
    selected, used = [], set()
    for name in ("front", "core", "peri"):
        taken = 0
        for index in np.argsort(-np.asarray(current[name], float), kind="stable"):
            if int(index) not in used:
                selected.append(int(index)); used.add(int(index)); taken += 1
            if taken >= quotas[name]:
                break
    if len(selected) < cap:
        support = np.maximum.reduce([np.asarray(current[name], float) for name in ("front", "core", "peri")])
        for index in np.argsort(-support, kind="stable"):
            if int(index) not in used:
                selected.append(int(index)); used.add(int(index))
            if len(selected) >= cap:
                break
    return np.sort(np.asarray(selected[:cap], dtype=np.int64))


def _pack(slides: list[dict], scaler, device):
    cap = min(MAX_SLIDE_TOKENS, max(48, MAX_PATIENT_TOKENS // max(len(slides), 1)))
    token_parts, weight_parts, slide_parts = [], [], []
    for slide_number, slide in enumerate(slides):
        token = _make_token(slide)
        current = _compartment_weights(slide)
        weight = np.stack([current[name] for name in ("core", "front", "peri")], axis=1).astype(np.float32)
        index = _select(slide, cap)
        token_parts.append(token[index]); weight_parts.append(weight[index]); slide_parts.append(np.full(len(index), slide_number, np.int64))
    token, weight, slide_index = np.concatenate(token_parts), np.concatenate(weight_parts), np.concatenate(slide_parts)
    if len(token) > MAX_PATIENT_TOKENS:
        keep = []
        for value in np.unique(slide_index):
            candidates = np.flatnonzero(slide_index == value)
            keep.extend(candidates[: min(24, len(candidates))].tolist())
        used = set(keep)
        for index in np.argsort(-weight.max(axis=1), kind="stable"):
            if int(index) not in used:
                keep.append(int(index)); used.add(int(index))
            if len(keep) >= MAX_PATIENT_TOKENS:
                break
        keep = np.sort(np.asarray(keep[:MAX_PATIENT_TOKENS], dtype=np.int64))
        token, weight, slide_index = token[keep], weight[keep], slide_index[keep]
    token = scaler.transform(token).astype(np.float32)
    valid = np.ones((1, len(token)), dtype=bool)
    return (
        torch.from_numpy(token[None]).to(device), torch.from_numpy(weight[None]).to(device),
        torch.from_numpy(valid).to(device), torch.from_numpy(slide_index[None]).to(device), len(slides),
    )


def predict(cache_files: list[str | Path], weights_dir: str | Path | None = None, device: str | None = None, case_id: str = "external_case") -> dict:
    weights_dir = Path(weights_dir) if weights_dir else Path(__file__).resolve().parent / "weights"
    target = torch.device(device or ("cuda:0" if torch.cuda.is_available() else "cpu"))
    fold_rows = []
    for fold in range(5):
        root = weights_dir / f"fold_{fold}"
        point_checkpoint = _torch_load(root / "patch_encoder.pt")
        point = PatchEncoder(int(point_checkpoint["image_dim"]), int(point_checkpoint["topology_dim"]), int(point_checkpoint["guided_dim"]), int(point_checkpoint["free_dim"]))
        point.load_state_dict(point_checkpoint["state_dict"], strict=True); point.to(target).eval()
        with (root / "spatial_transforms.pkl").open("rb") as stream:
            transforms = pickle.load(stream)
        field_checkpoint = _torch_load(root / "spatial_field.pt")
        field = SpatialField(field_checkpoint["input_dim"], field_checkpoint["image_dim"], field_checkpoint["topology_dim"], field_checkpoint["hidden_dim"], field_checkpoint["free_dim"])
        field.load_state_dict(field_checkpoint["state_dict"], strict=True); field.to(target).eval()
        risk_checkpoint = _torch_load(root / "survival_aggregator.pt")
        survival = SurvivalAggregator(int(risk_checkpoint["input_dim"]), int(risk_checkpoint["hidden_dim"]))
        survival.load_state_dict(risk_checkpoint["state_dict"], strict=True); survival.to(target).eval()
        with (root / "inference_contract.pkl").open("rb") as stream:
            contract = pickle.load(stream)
        slides = [_encode_cache(Path(path), transforms, field, point, target) for path in cache_files]
        packed = _pack(slides, contract["token_scaler"], target)
        with torch.inference_mode():
            raw, _ = survival(*packed)
        raw_value = float(raw.cpu().item())
        standardized = (raw_value - float(contract["train_risk_median"])) / max(float(contract["train_risk_iqr"]), 1e-8)
        fold_rows.append({"fold": fold, "raw_risk": raw_value, "risk": standardized})
        del point, field, survival
        if target.type == "cuda":
            torch.cuda.empty_cache()
    values = np.asarray([row["risk"] for row in fold_rows], dtype=float)
    risk = float(np.median(values))
    return {
        "model": "patch_feature",
        "case_id": case_id,
        "n_slides": len(cache_files),
        "risk": risk,
        "risk_group": "high" if risk >= 0.0 else "low",
        "threshold": 0.0,
        "fold_mean": float(values.mean()),
        "fold_sd": float(values.std(ddof=1)),
        "fold_predictions": fold_rows,
    }


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Predict CRC OS risk from prepared patch-feature caches")
    parser.add_argument("caches", nargs="+")
    parser.add_argument("--weights-dir", type=Path)
    parser.add_argument("--device")
    parser.add_argument("--case-id", default="external_case")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = predict(args.caches, args.weights_dir, args.device, args.case_id)
    text = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
