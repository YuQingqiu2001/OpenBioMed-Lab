"""Extract outcome-blind pure cell-topology region tokens from sparse WSI graphs."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


HERE = Path(__file__).resolve().parent
CONFIG_PATH = HERE / "token_config.json"


def load_config() -> dict:
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def config_sha256(cfg: dict) -> str:
    payload = json.dumps(cfg, sort_keys=True, ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    os.replace(temporary, path)


def region_adjacency(iy_ix: np.ndarray, component: np.ndarray) -> np.ndarray:
    lookup = {
        (int(comp), int(iy), int(ix)): index
        for index, ((iy, ix), comp) in enumerate(zip(iy_ix, component))
    }
    source: list[int] = []
    target: list[int] = []
    for index, ((iy, ix), comp) in enumerate(zip(iy_ix, component)):
        for dy, dx in ((0, 1), (1, 0)):
            other = lookup.get((int(comp), int(iy + dy), int(ix + dx)))
            if other is not None:
                source.extend((index, other))
                target.extend((other, index))
    if not source:
        return np.empty((2, 0), dtype=np.int64)
    return np.asarray([source, target], dtype=np.int64)


def aggregate_mean_std(values: np.ndarray, region: np.ndarray, n_regions: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    good = (region >= 0) & (region < n_regions) & np.isfinite(values).all(axis=1)
    r = region[good]
    x = values[good].astype(np.float64, copy=False)
    counts = np.bincount(r, minlength=n_regions).astype(np.float64)
    sums = np.zeros((n_regions, x.shape[1]), dtype=np.float64)
    sums2 = np.zeros_like(sums)
    np.add.at(sums, r, x)
    np.add.at(sums2, r, x * x)
    denominator = np.maximum(counts[:, None], 1.0)
    mean = sums / denominator
    variance = np.maximum(sums2 / denominator - mean * mean, 0.0)
    return mean.astype(np.float32), np.sqrt(variance).astype(np.float32), counts.astype(np.int64)


def pair_index_names(cell_types: list[str]) -> tuple[list[tuple[int, int]], list[str]]:
    pairs: list[tuple[int, int]] = []
    names: list[str] = []
    for first in range(len(cell_types)):
        for second in range(first, len(cell_types)):
            pairs.append((first, second))
            names.append(f"edge_log2_enrichment.{cell_types[first]}__{cell_types[second]}")
    return pairs, names


def pair_ids(first: np.ndarray, second: np.ndarray, n_types: int) -> np.ndarray:
    low = np.minimum(first, second)
    high = np.maximum(first, second)
    return low * n_types - (low * (low - 1)) // 2 + (high - low)


def expected_pair_proportions(q: np.ndarray, pairs: list[tuple[int, int]]) -> np.ndarray:
    expected = np.empty((q.shape[0], len(pairs)), dtype=np.float64)
    for index, (first, second) in enumerate(pairs):
        value = q[:, first] * q[:, second]
        if first != second:
            value *= 2.0
        expected[:, index] = value
    return expected


def typed_edge_enrichment(
    edge_region: np.ndarray,
    first_label: np.ndarray,
    second_label: np.ndarray,
    first_confidence: np.ndarray,
    second_confidence: np.ndarray,
    edge_weight: np.ndarray,
    n_regions: int,
    n_types: int,
    pairs: list[tuple[int, int]],
) -> np.ndarray:
    n_pairs = len(pairs)
    pids = pair_ids(first_label, second_label, n_types)
    confidence_weight = edge_weight * first_confidence * second_confidence
    observed = np.bincount(
        edge_region * n_pairs + pids,
        weights=confidence_weight,
        minlength=n_regions * n_pairs,
    ).reshape(n_regions, n_pairs)
    stub = np.zeros((n_regions, n_types), dtype=np.float64)
    np.add.at(stub, (edge_region, first_label), edge_weight * first_confidence)
    np.add.at(stub, (edge_region, second_label), edge_weight * second_confidence)
    q = stub / np.maximum(stub.sum(axis=1, keepdims=True), 1e-8)
    expected = expected_pair_proportions(q, pairs) * observed.sum(axis=1, keepdims=True)
    enrichment = np.log2((observed + 1.0) / (expected + 1.0))
    return np.clip(enrichment, -6.0, 6.0).astype(np.float32)


def complete_output(path: Path, expected_hash: str) -> bool:
    if not path.exists():
        return False
    try:
        with h5py.File(path, "r") as handle:
            return bool(handle.attrs.get("complete", 0)) and str(handle.attrs.get("config_sha256", "")) == expected_hash
    except OSError:
        return False


def extract_one(record: dict, cfg: dict, expected_hash: str) -> dict:
    started = time.time()
    sample_id = str(record["sample_id"])
    case_id = str(record["case_id"])
    graph_path = Path(cfg["input"]["graph_cache_dir"]) / f"{sample_id}_cell_graph.h5"
    output_dir = Path(cfg["output_root"]) / "01_region_tokens"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{sample_id}.topology_tokens.h5"
    if complete_output(output_path, expected_hash):
        return {"sample_id": sample_id, "case_id": case_id, "status": "skipped_complete", "runtime_seconds": 0.0}

    with h5py.File(graph_path, "r") as handle:
        n_cells = int(handle["cells/coords_um"].shape[0])
        region_index = handle["cells/region_index"][:].astype(np.int64, copy=False)
        edge_index = handle["graph/edge_index_undirected"][:].astype(np.int64, copy=False)
        edge_distance = handle["graph/edge_distance_um"][:].astype(np.float32, copy=False)
        iy_ix = handle["regions/region_iy_ix"][:].astype(np.int32, copy=False)
        component = handle["regions/region_component_id"][:].astype(np.int32, copy=False)
        tissue_fraction = handle["regions/region_tissue_fraction"][:].astype(np.float32, copy=False)
        source_h5 = Path(str(handle.attrs["source_h5"]))
        region_core_um = float(handle.attrs.get("region_core_um", 250.0))

    with h5py.File(source_h5, "r") as source:
        label_raw = source["nuclei/level2_class_id"][:].astype(np.int16, copy=False)
        confidence_raw = source["nuclei/level2_confidence"][:].astype(np.float32, copy=False)
        consensus_status = source["nuclei/consensus_status"][:].astype(np.uint8, copy=False)

    if len(label_raw) != n_cells:
        raise ValueError(f"cell count mismatch for {sample_id}: graph={n_cells}, nuclei={len(label_raw)}")

    cell_types = list(cfg["cell_types"])
    n_types = len(cell_types)
    uncertain = n_types - 1
    minimum_confidence = float(cfg["minimum_main_class_confidence"])
    resolved_main = (label_raw >= 0) & (label_raw < uncertain) & (confidence_raw >= minimum_confidence) & (consensus_status == 0)
    label = np.where(resolved_main, label_raw, uncertain).astype(np.int16)
    confidence = np.where(resolved_main, confidence_raw, 1.0).astype(np.float32)
    x0 = np.zeros((n_cells, n_types), dtype=np.float32)
    x0[np.arange(n_cells), label] = confidence
    x0[:, uncertain] += np.where(resolved_main, 1.0 - confidence, 0.0).astype(np.float32)

    u = edge_index[0]
    v = edge_index[1]
    kernel_scale = float(cfg["edge_kernel_scale_um"])
    edge_kernel = np.exp(-np.square(edge_distance / kernel_scale)).astype(np.float32)
    adjacency = sparse.csr_matrix((edge_kernel, (u, v)), shape=(n_cells, n_cells), dtype=np.float32)
    adjacency = adjacency + adjacency.T
    adjacency.setdiag(1.0)
    adjacency.sum_duplicates()
    degree = np.asarray(adjacency.sum(axis=1)).ravel().astype(np.float32)
    inverse_degree = 1.0 / np.maximum(degree, 1e-8)
    adjacency.data *= np.repeat(inverse_degree, np.diff(adjacency.indptr))

    h1 = adjacency @ x0
    h2 = adjacency @ h1
    global_composition = x0.mean(axis=0, dtype=np.float64).astype(np.float32)
    subtree = np.concatenate((h1 - global_composition, h2 - global_composition, h2 - h1), axis=1)
    subtree_mean, subtree_std, region_counts = aggregate_mean_std(subtree, region_index, len(iy_ix))

    entropy1 = -(np.clip(h1, 1e-8, 1.0) * np.log(np.clip(h1, 1e-8, 1.0))).sum(axis=1) / np.log(n_types)
    entropy2 = -(np.clip(h2, 1e-8, 1.0) * np.log(np.clip(h2, 1e-8, 1.0))).sum(axis=1) / np.log(n_types)
    agreement = (x0 * h1).sum(axis=1)
    scalar_node = np.stack((entropy1, entropy2, agreement), axis=1).astype(np.float32)
    scalar_mean, scalar_std, _ = aggregate_mean_std(scalar_node, region_index, len(iy_ix))

    valid_edge = (
        (region_index[u] >= 0)
        & (region_index[u] == region_index[v])
        & (region_index[u] < len(iy_ix))
    )
    edge_region = region_index[u[valid_edge]].astype(np.int64)
    edge_u = u[valid_edge]
    edge_v = v[valid_edge]
    within_distance = edge_distance[valid_edge]
    within_kernel = edge_kernel[valid_edge]
    pairs, pair_names = pair_index_names(cell_types)
    pair_enrichment = typed_edge_enrichment(
        edge_region,
        label[edge_u],
        label[edge_v],
        confidence[edge_u],
        confidence[edge_v],
        within_kernel,
        len(iy_ix),
        n_types,
        pairs,
    )

    selected_pairs = [(0, 2), (0, 3), (2, 3), (0, 4), (0, 5), (3, 4)]
    selected_names = [f"{cell_types[a]}__{cell_types[b]}" for a, b in selected_pairs]
    multiscale: list[np.ndarray] = []
    multiscale_names: list[str] = []
    for threshold in cfg["filtration_thresholds_um"]:
        keep = within_distance <= float(threshold)
        if keep.any():
            enrichment = typed_edge_enrichment(
                edge_region[keep],
                label[edge_u[keep]],
                label[edge_v[keep]],
                confidence[edge_u[keep]],
                confidence[edge_v[keep]],
                np.ones(int(keep.sum()), dtype=np.float32),
                len(iy_ix),
                n_types,
                pairs,
            )
        else:
            enrichment = np.zeros((len(iy_ix), len(pairs)), dtype=np.float32)
        pair_lookup = {pair: index for index, pair in enumerate(pairs)}
        multiscale.append(enrichment[:, [pair_lookup[pair] for pair in selected_pairs]])
        multiscale_names.extend([f"filtration_{int(threshold)}um.{name}" for name in selected_names])

    total_edges = np.bincount(edge_region, minlength=len(iy_ix)).astype(np.float64)
    edge_stats: list[np.ndarray] = []
    edge_stat_names: list[str] = []
    for threshold in cfg["filtration_thresholds_um"]:
        numerator = np.bincount(edge_region[within_distance <= float(threshold)], minlength=len(iy_ix))
        edge_stats.append((numerator / np.maximum(total_edges, 1.0)).astype(np.float32))
        edge_stat_names.append(f"filtration_edge_fraction_le_{int(threshold)}um")
    heterotypic = label[edge_u] != label[edge_v]
    typed = (label[edge_u] != uncertain) & (label[edge_v] != uncertain)
    edge_stats.append(
        (np.bincount(edge_region[heterotypic], minlength=len(iy_ix)) / np.maximum(total_edges, 1.0)).astype(np.float32)
    )
    edge_stat_names.append("heterotypic_edge_fraction")
    edge_stats.append(
        (np.bincount(edge_region[typed], minlength=len(iy_ix)) / np.maximum(total_edges, 1.0)).astype(np.float32)
    )
    edge_stat_names.append("resolved_typed_edge_fraction")

    base = np.concatenate(
        (
            subtree_mean,
            subtree_std,
            scalar_mean,
            scalar_std,
            pair_enrichment,
            *multiscale,
            np.stack(edge_stats, axis=1),
        ),
        axis=1,
    ).astype(np.float32)
    subtree_names = []
    for statistic in ("mean", "std"):
        for block in ("softwl1_minus_slide", "softwl2_minus_slide", "softwl2_minus_softwl1"):
            subtree_names.extend([f"{block}.{cell_type}.{statistic}" for cell_type in cell_types])
    scalar_names = [f"{name}.{statistic}" for statistic in ("mean", "std") for name in ("neighbor_entropy_h1", "neighbor_entropy_h2", "self_neighbor_agreement")]
    base_names = subtree_names + scalar_names + pair_names + multiscale_names + edge_stat_names
    if base.shape[1] != len(base_names):
        raise AssertionError((base.shape, len(base_names)))

    valid_region = (
        (region_counts >= int(cfg["minimum_region_cells"]))
        & (tissue_fraction >= float(cfg["minimum_region_tissue_fraction"]))
        & np.isfinite(base).all(axis=1)
    )
    region_edges = region_adjacency(iy_ix, component)
    if region_edges.shape[1]:
        edge_keep = valid_region[region_edges[0]] & valid_region[region_edges[1]]
        region_edges = region_edges[:, edge_keep]
    region_graph = sparse.csr_matrix(
        (np.ones(region_edges.shape[1], dtype=np.float32), (region_edges[0], region_edges[1])),
        shape=(len(iy_ix), len(iy_ix)),
    )
    region_graph.setdiag(1.0)
    region_graph.sum_duplicates()
    region_degree = np.asarray(region_graph.sum(axis=1)).ravel()
    region_graph.data *= np.repeat(1.0 / np.maximum(region_degree, 1.0), np.diff(region_graph.indptr))
    neighbor = region_graph @ np.where(valid_region[:, None], base, 0.0)
    region_tokens = np.concatenate((base, neighbor - base), axis=1)[valid_region]
    token_names = base_names + [f"region_neighbor_delta.{name}" for name in base_names]
    selected_old = np.flatnonzero(valid_region)
    old_to_new = np.full(len(iy_ix), -1, dtype=np.int64)
    old_to_new[selected_old] = np.arange(len(selected_old))
    if region_edges.shape[1]:
        selected_edges = old_to_new[region_edges]
    else:
        selected_edges = np.empty((2, 0), dtype=np.int64)
    region_xy = iy_ix[selected_old][:, ::-1].astype(np.float32) * region_core_um + region_core_um / 2.0

    composition_mean, _, _ = aggregate_mean_std(x0, region_index, len(iy_ix))
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    if temporary.exists():
        temporary.unlink()
    with h5py.File(temporary, "w") as output:
        output.attrs["schema"] = cfg["schema"] + "_region_tokens"
        output.attrs["complete"] = 1
        output.attrs["config_sha256"] = expected_hash
        output.attrs["sample_id"] = sample_id
        output.attrs["case_id"] = case_id
        output.attrs["source_graph_cache"] = str(graph_path)
        output.attrs["source_nuclei_h5"] = str(source_h5)
        output.attrs["n_source_cells"] = n_cells
        output.attrs["n_source_edges"] = edge_index.shape[1]
        output.attrs["region_core_um"] = region_core_um
        output.attrs["feature_semantics"] = "pure typed-cell topology; no morphology, colour, UNI, scBPS, stage or clinical covariates"
        output.create_dataset("tokens/features", data=region_tokens, compression="gzip", compression_opts=4)
        output.create_dataset("tokens/feature_names", data=np.asarray(token_names, dtype=h5py.string_dtype("utf-8")))
        output.create_dataset("tokens/region_xy_um", data=region_xy, compression="gzip", compression_opts=4)
        output.create_dataset("tokens/region_old_index", data=selected_old.astype(np.int32), compression="gzip", compression_opts=4)
        output.create_dataset("tokens/n_cells", data=region_counts[valid_region].astype(np.int32), compression="gzip", compression_opts=4)
        output.create_dataset("tokens/tissue_fraction", data=tissue_fraction[valid_region], compression="gzip", compression_opts=4)
        output.create_dataset("tokens/nuisance_cell_composition", data=composition_mean[valid_region], compression="gzip", compression_opts=4)
        output.create_dataset("tokens/nuisance_names", data=np.asarray([f"composition.{name}" for name in cell_types], dtype=h5py.string_dtype("utf-8")))
        output.create_dataset("graph/edge_index", data=selected_edges, compression="gzip", compression_opts=4)
    os.replace(temporary, output_path)
    return {
        "sample_id": sample_id,
        "case_id": case_id,
        "status": "complete",
        "n_cells": n_cells,
        "n_regions": len(iy_ix),
        "n_valid_tokens": int(valid_region.sum()),
        "runtime_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()
    cfg = load_config()
    expected_hash = config_sha256(cfg)
    output_root = Path(cfg["output_root"])
    manifest_dir = output_root / "00_manifests"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    audit = pd.read_csv(cfg["input"]["slide_audit_tsv"], sep="\t", dtype=str)
    audit = audit[(audit["valid_for_graph"] == "True") & (audit["clinical_linked"] == "True")].copy()
    records = audit[["sample_id", "case_id"]].drop_duplicates().to_dict("records")
    if args.limit is not None:
        records = records[: args.limit]
    workers = args.workers or int(cfg["extract_workers"])
    started = time.time()
    results: list[dict] = []
    failures: list[dict] = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(extract_one, record, cfg, expected_hash): record for record in records}
        for index, future in enumerate(as_completed(futures), start=1):
            record = futures[future]
            try:
                result = future.result()
                results.append(result)
            except Exception as exc:
                result = {"sample_id": record["sample_id"], "case_id": record["case_id"], "status": "failed", "error": repr(exc)}
                failures.append(result)
                results.append(result)
            if index == 1 or index % 10 == 0 or index == len(records):
                status = {
                    "schema": cfg["schema"] + "_extraction_status",
                    "processed": index,
                    "total": len(records),
                    "complete_or_skipped": sum(row["status"] in {"complete", "skipped_complete"} for row in results),
                    "failed": len(failures),
                    "last": result,
                    "elapsed_seconds": time.time() - started,
                }
                atomic_json(manifest_dir / "extraction_status_live.json", status)
                print(json.dumps(status, ensure_ascii=False), flush=True)
    table = pd.DataFrame(results).sort_values("sample_id")
    table.to_csv(manifest_dir / "region_token_status.tsv", sep="\t", index=False)
    summary = {
        "schema": cfg["schema"] + "_extraction_summary",
        "config_sha256": expected_hash,
        "n_requested": len(records),
        "n_complete_or_skipped": int(table["status"].isin(["complete", "skipped_complete"]).sum()),
        "n_failed": int((table["status"] == "failed").sum()),
        "elapsed_seconds": time.time() - started,
        "locked_constraints": cfg["locked_constraints"],
    }
    atomic_json(manifest_dir / "extraction_summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
