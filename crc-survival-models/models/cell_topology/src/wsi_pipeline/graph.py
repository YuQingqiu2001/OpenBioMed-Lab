from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree


PROBABILITY_NAMES = [
    "p_tumour",
    "p_non_tumour_epithelial",
    "p_fibroblast_like",
    "p_lymphocyte",
    "p_plasma",
    "p_neutrophil",
    "p_eosinophil",
    "p_dead",
    "p_spindle_candidate",
    "p_macrophage_candidate",
    "p_mitotic_tumour_candidate",
]

MORPHOLOGY_NAMES = [
    "log1p_area_um2",
    "log1p_perimeter_um",
    "circularity",
    "log1p_convex_area_um2",
    "solidity",
    "log1p_bbox_width_um",
    "log1p_bbox_height_um",
    "log_aspect_ratio",
    "log1p_equivalent_diameter_um",
    "log1p_major_axis_um",
    "log1p_minor_axis_um",
    "eccentricity",
]

FEATURE_NAMES = (
    PROBABILITY_NAMES
    + MORPHOLOGY_NAMES
    + ["level1_confidence", "level2_confidence"]
    + [f"consensus_status_{i}" for i in range(4)]
    + ["candidate_macrophage", "candidate_mitotic", "candidate_spindle"]
    + [f"missing_{name}" for name in MORPHOLOGY_NAMES]
    + ["full_embedding_available"]
)

MAIN_TYPE_INDICES = np.arange(0, 8, dtype=np.int64)
CANDIDATE_TYPE_INDICES = np.arange(8, 11, dtype=np.int64)
MORPH_INDICES = np.arange(11, 23, dtype=np.int64)
NORMALIZE_INDICES = np.asarray(list(range(11, 25)), dtype=np.int64)


def _decode_json_attr(value) -> dict | list:
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    return json.loads(str(value))


def _physical_morphology(morph: np.ndarray, target_mpp: float) -> tuple[np.ndarray, np.ndarray]:
    """Convert pixel measurements to physical, rotation-safe morphology.

    The orientation angle is intentionally excluded from the primary feature vector.
    It would otherwise break the rotation-invariance requirement unless every
    augmentation transformed the angle consistently.
    """
    morph = np.asarray(morph, dtype=np.float32)
    out = np.empty((len(morph), 12), dtype=np.float32)
    area_scale = float(target_mpp) ** 2
    length_scale = float(target_mpp)
    out[:, 0] = np.log1p(np.clip(morph[:, 0] * area_scale, 0, None))
    out[:, 1] = np.log1p(np.clip(morph[:, 1] * length_scale, 0, None))
    out[:, 2] = morph[:, 2]
    out[:, 3] = np.log1p(np.clip(morph[:, 3] * area_scale, 0, None))
    out[:, 4] = morph[:, 4]
    out[:, 5] = np.log1p(np.clip(morph[:, 5] * length_scale, 0, None))
    out[:, 6] = np.log1p(np.clip(morph[:, 6] * length_scale, 0, None))
    out[:, 7] = np.log(np.clip(morph[:, 7], 1e-4, None))
    out[:, 8] = np.log1p(np.clip(morph[:, 8] * length_scale, 0, None))
    out[:, 9] = np.log1p(np.clip(morph[:, 9] * length_scale, 0, None))
    out[:, 10] = np.log1p(np.clip(morph[:, 10] * length_scale, 0, None))
    out[:, 11] = morph[:, 11]
    missing = ~np.isfinite(out)
    out[missing] = 0.0
    return out, missing.astype(np.float32)


def load_slide_arrays(h5_path: str | Path) -> dict:
    h5_path = Path(h5_path)
    with h5py.File(h5_path, "r") as h5:
        coords = h5["nuclei/centroid_global_um"][:].astype(np.float32, copy=False)
        probabilities = np.stack(
            [h5[f"nuclei/{name}"][:] for name in PROBABILITY_NAMES], axis=1
        ).astype(np.float32, copy=False)
        roi_record = _decode_json_attr(h5.attrs["roi_record_json"])
        target_mpp = float(roi_record["target_mpp"])
        morph, morph_missing = _physical_morphology(h5["nuclei/morphology"][:], target_mpp)
        l1_conf = h5["nuclei/level1_confidence"][:].astype(np.float32)[:, None]
        l2_conf = h5["nuclei/level2_confidence"][:].astype(np.float32)[:, None]
        status = h5["nuclei/consensus_status"][:].astype(np.int64)
        status_onehot = np.eye(4, dtype=np.float32)[np.clip(status, 0, 3)]
        flags = np.stack(
            [
                h5["nuclei/candidate_macrophage_flag"][:],
                h5["nuclei/candidate_mitotic_flag"][:],
                h5["nuclei/candidate_spindle_flag"][:],
            ],
            axis=1,
        ).astype(np.float32)
        embedding_available = np.zeros((len(coords), 1), dtype=np.float32)
        features = np.concatenate(
            [
                probabilities,
                morph,
                l1_conf,
                l2_conf,
                status_onehot,
                flags,
                morph_missing,
                embedding_available,
            ],
            axis=1,
        ).astype(np.float32, copy=False)
        if features.shape[1] != len(FEATURE_NAMES):
            raise RuntimeError(f"Feature contract mismatch: {features.shape[1]} != {len(FEATURE_NAMES)}")
        return {
            "coords_um": coords,
            "features_raw": features,
            "probabilities": probabilities,
            "target_mpp": target_mpp,
            "tissue_mask_manifest": str(h5.attrs["tissue_mask_manifest"]),
            "schema": str(h5.attrs["schema"]),
        }


def map_tissue_components(
    coords_um: np.ndarray, manifest_path: str | Path
) -> tuple[np.ndarray, dict, np.ndarray]:
    manifest_path = Path(manifest_path)
    with manifest_path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    mask_npz = Path(manifest["mask_npz"])
    with np.load(mask_npz) as payload:
        mask = payload["mask"].astype(bool)
        mask_mpp = float(np.asarray(payload["mpp"]).ravel()[0])
    labels, n_components = ndimage.label(mask, structure=np.ones((3, 3), dtype=np.uint8))
    origin_x, origin_y = map(float, manifest.get("mask_origin_um", [0.0, 0.0]))
    col = np.floor((coords_um[:, 0] - origin_x) / mask_mpp).astype(np.int64)
    row = np.floor((coords_um[:, 1] - origin_y) / mask_mpp).astype(np.int64)
    row = np.clip(row, 0, mask.shape[0] - 1)
    col = np.clip(col, 0, mask.shape[1] - 1)
    component = labels[row, col].astype(np.int32)
    remapped = component == 0
    if remapped.any() and mask.any():
        nearest = ndimage.distance_transform_edt(~mask, return_distances=False, return_indices=True)
        component[remapped] = labels[
            nearest[0, row[remapped], col[remapped]],
            nearest[1, row[remapped], col[remapped]],
        ]
    metadata = {
        "mask_npz": str(mask_npz),
        "mask_mpp": mask_mpp,
        "mask_origin_um": [origin_x, origin_y],
        "mask_shape": list(mask.shape),
        "n_components": int(n_components),
        "n_boundary_cells_remapped": int(remapped.sum()),
    }
    return component, metadata, mask


def build_symmetric_knn_union(
    coords_um: np.ndarray,
    component: np.ndarray,
    k: int = 8,
    max_edge_um: float = 50.0,
    query_chunk: int = 200_000,
    query_workers: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    n = len(coords_um)
    if n < 2:
        return np.empty((2, 0), np.int64), np.empty(0, np.float32)
    tree = cKDTree(coords_um)
    all_u: list[np.ndarray] = []
    all_v: list[np.ndarray] = []
    all_d: list[np.ndarray] = []
    n_query = min(k + 1, n)
    for start in range(0, n, query_chunk):
        stop = min(n, start + query_chunk)
        distance, neighbor = tree.query(coords_um[start:stop], k=n_query, workers=query_workers)
        if n_query == 1:
            distance = distance[:, None]
            neighbor = neighbor[:, None]
        source = np.repeat(np.arange(start, stop, dtype=np.int64), n_query)
        target = neighbor.reshape(-1).astype(np.int64)
        dist = distance.reshape(-1).astype(np.float32)
        valid = (
            (source != target)
            & np.isfinite(dist)
            & (dist <= max_edge_um)
            & (component[source] > 0)
            & (component[source] == component[target])
        )
        source, target, dist = source[valid], target[valid], dist[valid]
        all_u.append(np.minimum(source, target))
        all_v.append(np.maximum(source, target))
        all_d.append(dist)
    if not all_u:
        return np.empty((2, 0), np.int64), np.empty(0, np.float32)
    u = np.concatenate(all_u)
    v = np.concatenate(all_v)
    d = np.concatenate(all_d)
    # A single int64 key is materially faster and less memory-hungry than a
    # two-column lexsort on million-cell slides.
    key = u * np.int64(n) + v
    _, first = np.unique(key, return_index=True)
    edge_index = np.stack([u[first], v[first]], axis=0)
    return edge_index, d[first]


def assign_regions(
    coords_um: np.ndarray,
    mask: np.ndarray,
    mask_meta: dict,
    core_um: float,
) -> dict:
    origin_x, origin_y = mask_meta["mask_origin_um"]
    ix = np.floor((coords_um[:, 0] - origin_x) / core_um).astype(np.int32)
    iy = np.floor((coords_um[:, 1] - origin_y) / core_um).astype(np.int32)
    cell_pairs = np.stack([iy, ix], axis=1)

    mask_mpp = float(mask_meta["mask_mpp"])
    rows, cols = np.nonzero(mask)
    mask_labels, _ = ndimage.label(mask, structure=np.ones((3, 3), dtype=np.uint8))
    mask_component = mask_labels[rows, cols].astype(np.int32)
    mask_x = origin_x + (cols + 0.5) * mask_mpp
    mask_y = origin_y + (rows + 0.5) * mask_mpp
    mask_ix = np.floor((mask_x - origin_x) / core_um).astype(np.int32)
    mask_iy = np.floor((mask_y - origin_y) / core_um).astype(np.int32)
    unique_mask_pairs, mask_counts = np.unique(
        np.stack([mask_iy, mask_ix], axis=1), axis=0, return_counts=True
    )
    # Preserve tissue tiles with zero detected nuclei. They carry an explicit
    # empty-tissue flag downstream, so a detection failure cannot silently
    # drop tissue from the region graph.
    region_pairs, combined_inverse = np.unique(
        np.concatenate([cell_pairs, unique_mask_pairs], axis=0),
        axis=0,
        return_inverse=True,
    )
    inverse = combined_inverse[: len(cell_pairs)]
    counts = np.bincount(inverse, minlength=len(region_pairs)).astype(np.int32)
    tissue_counter: dict[tuple[int, int], int] = {}
    for pair, count in zip(unique_mask_pairs, mask_counts):
        tissue_counter[(int(pair[0]), int(pair[1]))] = int(count)
    tissue_area = np.asarray(
        [tissue_counter.get((int(y), int(x)), 0) * mask_mpp**2 for y, x in region_pairs],
        dtype=np.float32,
    )
    component_counter: dict[tuple[int, int], tuple[int, int]] = {}
    component_triplets, component_counts = np.unique(
        np.stack([mask_iy, mask_ix, mask_component], axis=1), axis=0, return_counts=True
    )
    for triplet, count in zip(component_triplets, component_counts):
        key = (int(triplet[0]), int(triplet[1]))
        previous = component_counter.get(key)
        if previous is None or int(count) > previous[1]:
            component_counter[key] = (int(triplet[2]), int(count))
    region_component = np.asarray(
        [component_counter.get((int(y), int(x)), (0, 0))[0] for y, x in region_pairs],
        dtype=np.int32,
    )
    tissue_fraction = np.clip(tissue_area / (core_um**2), 0.0, 1.0)
    return {
        "region_index": inverse.astype(np.int32),
        "region_iy_ix": region_pairs.astype(np.int32),
        "region_n_core": counts,
        "region_tissue_area_um2": tissue_area,
        "region_tissue_fraction": tissue_fraction.astype(np.float32),
        "region_component_id": region_component,
    }


def build_slide_cache(
    source_h5: str | Path,
    destination_h5: str | Path,
    *,
    k: int = 8,
    max_edge_um: float = 50.0,
    core_um: float = 250.0,
) -> dict:
    source_h5, destination_h5 = Path(source_h5), Path(destination_h5)
    payload = load_slide_arrays(source_h5)
    component, mask_meta, mask = map_tissue_components(
        payload["coords_um"], payload["tissue_mask_manifest"]
    )
    edge_index, edge_distance = build_symmetric_knn_union(
        payload["coords_um"], component, k=k, max_edge_um=max_edge_um
    )
    regions = assign_regions(payload["coords_um"], mask, mask_meta, core_um)
    destination_h5.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination_h5.with_suffix(destination_h5.suffix + ".tmp")
    with h5py.File(tmp, "w") as out:
        out.attrs["schema"] = "whole_slide_cell_graph_cache_v1_1"
        out.attrs["complete"] = 0
        out.attrs["source_h5"] = str(source_h5)
        out.attrs["source_schema"] = payload["schema"]
        out.attrs["feature_names_json"] = json.dumps(FEATURE_NAMES)
        out.attrs["normalize_indices_json"] = json.dumps(NORMALIZE_INDICES.tolist())
        out.attrs["main_type_indices_json"] = json.dumps(MAIN_TYPE_INDICES.tolist())
        out.attrs["candidate_type_indices_json"] = json.dumps(CANDIDATE_TYPE_INDICES.tolist())
        out.attrs["morph_indices_json"] = json.dumps(MORPH_INDICES.tolist())
        out.attrs["graph_k"] = int(k)
        out.attrs["max_edge_um"] = float(max_edge_um)
        out.attrs["region_core_um"] = float(core_um)
        out.attrs["tissue_metadata_json"] = json.dumps(mask_meta)
        out.create_dataset("cells/coords_um", data=payload["coords_um"], compression="gzip", chunks=True)
        out.create_dataset("cells/features_raw", data=payload["features_raw"], compression="gzip", chunks=True)
        out.create_dataset("cells/component_id", data=component, compression="gzip", chunks=True)
        out.create_dataset("cells/region_index", data=regions["region_index"], compression="gzip", chunks=True)
        out.create_dataset("graph/edge_index_undirected", data=edge_index, compression="gzip", chunks=True)
        out.create_dataset("graph/edge_distance_um", data=edge_distance, compression="gzip", chunks=True)
        for name, values in regions.items():
            if name == "region_index":
                continue
            out.create_dataset(f"regions/{name}", data=values, compression="gzip", chunks=True)
        out.attrs["complete"] = 1
    tmp.replace(destination_h5)
    return {
        "source_h5": str(source_h5),
        "cache_h5": str(destination_h5),
        "n_cells": int(len(payload["coords_um"])),
        "n_edges_undirected": int(edge_index.shape[1]),
        "n_regions": int(len(regions["region_iy_ix"])),
        "n_tissue_components": int(mask_meta["n_components"]),
        "n_boundary_cells_remapped": int(mask_meta["n_boundary_cells_remapped"]),
        "core_assignment_exact": bool(len(regions["region_index"]) == len(payload["coords_um"])),
    }


@dataclass
class CachedSlide:
    path: Path
    coords_um: np.ndarray
    features_raw: np.ndarray
    component_id: np.ndarray
    region_index: np.ndarray
    edge_index: np.ndarray
    edge_distance_um: np.ndarray
    region_iy_ix: np.ndarray
    region_n_core: np.ndarray
    region_tissue_area_um2: np.ndarray
    region_tissue_fraction: np.ndarray
    region_component_id: np.ndarray
    core_um: float
    origin_um: tuple[float, float]

    @classmethod
    def load(cls, path: str | Path) -> "CachedSlide":
        path = Path(path)
        with h5py.File(path, "r") as h5:
            if not bool(h5.attrs.get("complete", 0)):
                raise RuntimeError(f"Incomplete graph cache: {path}")
            tissue_metadata = _decode_json_attr(h5.attrs["tissue_metadata_json"])
            return cls(
                path=path,
                coords_um=h5["cells/coords_um"][:],
                features_raw=h5["cells/features_raw"][:],
                component_id=h5["cells/component_id"][:],
                region_index=h5["cells/region_index"][:],
                edge_index=h5["graph/edge_index_undirected"][:].astype(np.int64),
                edge_distance_um=h5["graph/edge_distance_um"][:],
                region_iy_ix=h5["regions/region_iy_ix"][:],
                region_n_core=h5["regions/region_n_core"][:],
                region_tissue_area_um2=h5["regions/region_tissue_area_um2"][:],
                region_tissue_fraction=h5["regions/region_tissue_fraction"][:],
                region_component_id=(
                    h5["regions/region_component_id"][:]
                    if "region_component_id" in h5["regions"]
                    else np.zeros(len(h5["regions/region_iy_ix"]), dtype=np.int32)
                ),
                core_um=float(h5.attrs["region_core_um"]),
                origin_um=tuple(map(float, tissue_metadata["mask_origin_um"])),
            )


@dataclass
class RobustFeatureNormalizer:
    center: np.ndarray
    scale: np.ndarray
    normalize_indices: np.ndarray

    @classmethod
    def fit(
        cls,
        slides: Iterable[CachedSlide],
        normalize_indices: np.ndarray = NORMALIZE_INDICES,
        maximum_rows: int = 500_000,
        seed: int = 20260907,
    ) -> "RobustFeatureNormalizer":
        slides = list(slides)
        rng = np.random.default_rng(seed)
        per_slide = max(1, maximum_rows // max(1, len(slides)))
        sampled = []
        for slide in slides:
            n = len(slide.features_raw)
            idx = rng.choice(n, size=min(n, per_slide), replace=False)
            sampled.append(slide.features_raw[idx][:, normalize_indices])
        values = np.concatenate(sampled, axis=0)
        center = np.nanmedian(values, axis=0).astype(np.float32)
        q25, q75 = np.nanquantile(values, [0.25, 0.75], axis=0)
        scale = np.asarray(q75 - q25, dtype=np.float32)
        scale[~np.isfinite(scale) | (scale < 1e-5)] = 1.0
        return cls(center=center, scale=scale, normalize_indices=np.asarray(normalize_indices))

    @classmethod
    def fit_paths(
        cls,
        cache_paths: Iterable[str | Path],
        normalize_indices: np.ndarray = NORMALIZE_INDICES,
        maximum_rows: int = 500_000,
        seed: int = 20260907,
    ) -> "RobustFeatureNormalizer":
        paths = [Path(path) for path in cache_paths]
        rng = np.random.default_rng(seed)
        per_slide = max(1, maximum_rows // max(1, len(paths)))
        sampled = []
        for path in paths:
            with h5py.File(path, "r") as h5:
                dataset = h5["cells/features_raw"]
                n = len(dataset)
                idx = np.sort(rng.choice(n, size=min(n, per_slide), replace=False))
                sampled.append(dataset[idx][:, normalize_indices])
        values = np.concatenate(sampled, axis=0)
        center = np.nanmedian(values, axis=0).astype(np.float32)
        q25, q75 = np.nanquantile(values, [0.25, 0.75], axis=0)
        scale = np.asarray(q75 - q25, dtype=np.float32)
        scale[~np.isfinite(scale) | (scale < 1e-5)] = 1.0
        return cls(center=center, scale=scale, normalize_indices=np.asarray(normalize_indices))

    def transform(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=np.float32).copy()
        x[:, self.normalize_indices] = (
            x[:, self.normalize_indices] - self.center[None, :]
        ) / self.scale[None, :]
        x[:, self.normalize_indices] = np.clip(x[:, self.normalize_indices], -8.0, 8.0)
        return x

    def to_dict(self) -> dict:
        return {
            "center": self.center.tolist(),
            "scale": self.scale.tolist(),
            "normalize_indices": self.normalize_indices.tolist(),
            "feature_names": FEATURE_NAMES,
        }


def extract_region_subgraph(
    slide: CachedSlide,
    region_id: int,
    halo_um: float,
    normalizer: RobustFeatureNormalizer,
) -> dict:
    iy, ix = map(int, slide.region_iy_ix[region_id])
    origin_x, origin_y = slide.origin_um
    x0, x1 = origin_x + ix * slide.core_um, origin_x + (ix + 1) * slide.core_um
    y0, y1 = origin_y + iy * slide.core_um, origin_y + (iy + 1) * slide.core_um
    c = slide.coords_um
    selected = (
        (c[:, 0] >= x0 - halo_um)
        & (c[:, 0] < x1 + halo_um)
        & (c[:, 1] >= y0 - halo_um)
        & (c[:, 1] < y1 + halo_um)
    )
    node_global = np.flatnonzero(selected)
    core_global = np.flatnonzero(slide.region_index == region_id)
    lookup = np.full(len(c), -1, dtype=np.int64)
    lookup[node_global] = np.arange(len(node_global), dtype=np.int64)
    u, v = slide.edge_index
    edge_keep = selected[u] & selected[v]
    local_u = lookup[u[edge_keep]]
    local_v = lookup[v[edge_keep]]
    undirected = np.stack([local_u, local_v], axis=0)
    directed = np.concatenate([undirected, undirected[::-1]], axis=1)
    distance = np.concatenate(
        [slide.edge_distance_um[edge_keep], slide.edge_distance_um[edge_keep]], axis=0
    ).astype(np.float32)
    core_local = lookup[core_global]
    if (core_local < 0).any():
        raise RuntimeError("Core-to-halo assignment failure")
    return {
        "x": normalizer.transform(slide.features_raw[node_global]),
        "edge_index": directed,
        "edge_distance_um": distance,
        "core_local_index": core_local,
        "global_node_index": node_global,
        "region_id": int(region_id),
        "tissue_area_um2": float(slide.region_tissue_area_um2[region_id]),
        "tissue_fraction": float(slide.region_tissue_fraction[region_id]),
    }


def build_region_graph(slide: CachedSlide) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create an undirected adjacent-region graph with observed cross-cell support."""
    pair_to_id = {
        (int(iy), int(ix)): region_id
        for region_id, (iy, ix) in enumerate(slide.region_iy_ix)
    }
    candidate_u: list[int] = []
    candidate_v: list[int] = []
    for region_id, (iy_raw, ix_raw) in enumerate(slide.region_iy_ix):
        iy, ix = int(iy_raw), int(ix_raw)
        component = int(slide.region_component_id[region_id])
        if component <= 0:
            continue
        for dy, dx in [(0, 1), (1, -1), (1, 0), (1, 1)]:
            neighbor = pair_to_id.get((iy + dy, ix + dx))
            if neighbor is None or int(slide.region_component_id[neighbor]) != component:
                continue
            candidate_u.append(region_id)
            candidate_v.append(neighbor)
    if not candidate_u:
        return np.empty((2, 0), np.int64), np.empty(0, np.float32), np.empty(0, np.float32)
    edge_index = np.asarray([candidate_u, candidate_v], dtype=np.int64)
    delta = slide.region_iy_ix[edge_index[0]] - slide.region_iy_ix[edge_index[1]]
    distance = np.sqrt((delta.astype(np.float32) ** 2).sum(axis=1)) * slide.core_um

    cell_u, cell_v = slide.edge_index
    region_u = slide.region_index[cell_u]
    region_v = slide.region_index[cell_v]
    cross = region_u != region_v
    ru = np.minimum(region_u[cross], region_v[cross]).astype(np.int64)
    rv = np.maximum(region_u[cross], region_v[cross]).astype(np.int64)
    if len(ru):
        observed_pairs, observed_counts = np.unique(
            np.stack([ru, rv], axis=1), axis=0, return_counts=True
        )
        count_map = {
            (int(pair[0]), int(pair[1])): int(count)
            for pair, count in zip(observed_pairs, observed_counts)
        }
    else:
        count_map = {}
    cross_count = np.asarray(
        [count_map.get((int(u), int(v)), 0) for u, v in edge_index.T], dtype=np.float32
    )
    return edge_index, distance.astype(np.float32), cross_count
