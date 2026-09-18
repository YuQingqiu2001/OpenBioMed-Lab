#!/usr/bin/env python3
"""Create compact multi-head CellViT++ nucleus HDF5 with CRC consensus labels."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Dict, Mapping, Tuple

import cv2
import h5py
import numpy as np
import snappy
import torch
import ujson

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from common import atomic_json, sha256_file
from config import (
    CLASSIFIER_FILES,
    CLASSIFIER_ROOT,
    CELLVIT_ROOT,
    FROZEN,
    LEVEL1_LABELS,
    LEVEL2_LABELS,
)

sys.path.insert(0, str(CELLVIT_ROOT))
from cellvit.models.classifier.linear_classifier import LinearClassifier
from tissue_mask_utils import load_tissue_mask, sample_points


PANNUKE_LABELS = {
    0: "Background",
    1: "Neoplastic",
    2: "Inflammatory",
    3: "Connective",
    4: "Dead",
    5: "Epithelial",
}

FALLBACK_HEAD_LABELS = {
    "nucls_main": {
        0: "tumor_nonMitotic",
        1: "tumor_mitotic",
        2: "nonTILnonMQ_stromal",
        3: "macrophage",
        4: "lymphocyte",
        5: "plasma_cell",
        6: "other_nucleus",
    },
    "nucls_super": {
        0: "tumor_any",
        1: "nonTIL_stromal",
        2: "sTIL",
        3: "other_nucleus",
    },
}

MORPHOLOGY_FEATURES = (
    "area_px2",
    "perimeter_px",
    "circularity",
    "convex_area_px2",
    "solidity",
    "bbox_width_px",
    "bbox_height_px",
    "aspect_ratio",
    "equivalent_diameter_px",
    "major_axis_px",
    "minor_axis_px",
    "eccentricity",
    "orientation_deg",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", type=Path, required=True)
    parser.add_argument("--cells-json-snappy", type=Path, required=True)
    parser.add_argument("--roi-record", type=Path, required=True)
    parser.add_argument("--roi-uid", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tissue-mask-manifest", type=Path)
    parser.add_argument("--batch-size", type=int, default=32768)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def normalize_label(label: str) -> str:
    return "".join(character.lower() for character in str(label) if character.isalnum())


def load_classifier(path: Path, device: torch.device, head_name: str | None = None):
    checkpoint = torch.load(path, map_location="cpu")
    config = checkpoint["config"]
    state = checkpoint["model_state_dict"]
    labels = {
        int(key.rsplit(".", 1)[1]): str(value)
        for key, value in config.items()
        if key.startswith("data.label_map.")
    }
    if not labels and head_name in FALLBACK_HEAD_LABELS:
        labels = FALLBACK_HEAD_LABELS[str(head_name)]
    model = LinearClassifier(
        embed_dim=int(state["fc1.weight"].shape[1]),
        hidden_dim=int(config["model.hidden_dim"]),
        num_classes=int(config["data.num_classes"]),
        drop_rate=0,
    )
    model.load_state_dict(state)
    model.eval().to(device)
    return model, labels


def probability_by_label(
    probabilities: np.ndarray, labels: Mapping[int, str], candidates: Tuple[str, ...]
) -> np.ndarray:
    normalized = {normalize_label(label): index for index, label in labels.items()}
    columns = [normalized[name] for name in candidates if name in normalized]
    if not columns:
        return np.zeros(probabilities.shape[0], dtype=np.float32)
    return probabilities[:, columns].astype(np.float32).sum(axis=1)


def read_indexed_rows_dense(
    dataset: h5py.Dataset,
    indices: np.ndarray,
    *,
    source_span_rows: int = 65536,
    output_dtype: np.dtype | type | None = None,
    progress_label: str | None = None,
) -> np.ndarray:
    """Read sorted HDF5 rows through bounded dense spans.

    h5py point selection is extremely slow for large, sparse index arrays. The
    CellViT tissue filter and review set both produce sorted unique indices, so
    reading bounded dense source spans and selecting in RAM is equivalent while
    avoiding millions of HDF5 point-selection operations.
    """
    rows = np.asarray(indices, dtype=np.int64)
    if rows.ndim != 1:
        raise ValueError("indices must be one-dimensional")
    if rows.size and (rows[0] < 0 or rows[-1] >= dataset.shape[0]):
        raise IndexError("row index outside HDF5 dataset")
    if rows.size > 1 and np.any(np.diff(rows) <= 0):
        raise ValueError("indices must be sorted and unique")
    dtype = np.dtype(output_dtype) if output_dtype is not None else dataset.dtype
    output = np.empty((rows.size, *dataset.shape[1:]), dtype=dtype)
    if not rows.size:
        return output

    source_span_rows = max(1, int(source_span_rows))
    cursor = 0
    last_reported = -1
    while cursor < rows.size:
        source_start = int(rows[cursor])
        stop = int(np.searchsorted(rows, source_start + source_span_rows, side="left"))
        stop = max(stop, cursor + 1)
        source_stop = int(rows[stop - 1]) + 1
        dense = dataset[source_start:source_stop]
        selected = dense[rows[cursor:stop] - source_start]
        output[cursor:stop] = selected.astype(dtype, copy=False)
        cursor = stop
        if progress_label:
            percent = int(round(100.0 * cursor / rows.size))
            bucket = min(10, percent // 10)
            if bucket != last_reported or cursor == rows.size:
                print(
                    f"{progress_label} {cursor}/{rows.size} ({percent}%)",
                    flush=True,
                )
                last_reported = bucket
    return output


def pannuke_signal(types: np.ndarray, confidence: np.ndarray, class_id: int) -> np.ndarray:
    return np.where(types == class_id, confidence, 0.0).astype(np.float32)


def contour_features(contour_xy: np.ndarray, bbox: np.ndarray) -> np.ndarray:
    result = np.full(len(MORPHOLOGY_FEATURES), np.nan, dtype=np.float32)
    if contour_xy.shape[0] < 3:
        return result
    contour = contour_xy.astype(np.float32).reshape((-1, 1, 2))
    area = float(abs(cv2.contourArea(contour)))
    perimeter = float(cv2.arcLength(contour, True))
    hull = cv2.convexHull(contour)
    convex_area = float(abs(cv2.contourArea(hull)))
    width = float(abs(bbox[1, 0] - bbox[0, 0]))
    height = float(abs(bbox[1, 1] - bbox[0, 1]))
    circularity = 4.0 * math.pi * area / perimeter**2 if perimeter else np.nan
    solidity = area / convex_area if convex_area else np.nan
    aspect_ratio = width / height if height else np.nan
    equivalent_diameter = math.sqrt(4.0 * area / math.pi) if area else 0.0
    major_axis = minor_axis = eccentricity = orientation = np.nan
    if contour_xy.shape[0] >= 5:
        (_, _), (axis_a, axis_b), angle = cv2.fitEllipse(contour)
        major_axis, minor_axis = float(max(axis_a, axis_b)), float(min(axis_a, axis_b))
        if major_axis > 0:
            eccentricity = math.sqrt(max(0.0, 1.0 - (minor_axis / major_axis) ** 2))
        orientation = float(angle if axis_a >= axis_b else (angle + 90.0) % 180.0)
    result[:] = (
        area,
        perimeter,
        circularity,
        convex_area,
        solidity,
        width,
        height,
        aspect_ratio,
        equivalent_diameter,
        major_axis,
        minor_axis,
        eccentricity,
        orientation,
    )
    return result


def consensus(
    pannuke_type: np.ndarray,
    pannuke_conf: np.ndarray,
    head_prob: Dict[str, np.ndarray],
    head_labels: Dict[str, Dict[int, str]],
) -> Dict[str, np.ndarray]:
    n = pannuke_type.shape[0]
    zero = np.zeros(n, dtype=np.float32)

    def hp(head: str, *labels: str) -> np.ndarray:
        if head not in head_prob:
            return zero.copy()
        return probability_by_label(head_prob[head], head_labels[head], tuple(labels))

    pan_neoplastic = pannuke_signal(pannuke_type, pannuke_conf, 1)
    pan_inflammatory = pannuke_signal(pannuke_type, pannuke_conf, 2)
    pan_connective = pannuke_signal(pannuke_type, pannuke_conf, 3)
    pan_dead = pannuke_signal(pannuke_type, pannuke_conf, 4)
    pan_epithelial = pannuke_signal(pannuke_type, pannuke_conf, 5)

    ocelot_tumour = hp("ocelot", "tumorcell", "tumourcell")
    nucls_tumour = hp("nucls_super", "tumorany", "tumourany")
    lizard_epithelial = hp("lizard", "epithelial")
    panop_epithelial = hp("panoptils", "epithelialcells", "epithelial")
    lizard_connective = hp("lizard", "connectivetissue", "connective")
    panop_stromal = hp("panoptils", "stromalcells", "stromal")
    consep_spindle = hp("consep", "spindleshaped", "spindle")
    lizard_lymph = hp("lizard", "lymphocyte")
    panop_til = hp("panoptils", "tils", "til")
    lizard_plasma = hp("lizard", "plasma", "plasmacell")
    lizard_neut = hp("lizard", "neutrophil")
    lizard_eos = hp("lizard", "eosinophil")
    nucls_mac = hp("nucls_main", "macrophage")
    midog_mitotic = hp("midog", "mitotic")

    p_tumour = np.clip(0.50 * pan_neoplastic + 0.35 * ocelot_tumour + 0.15 * nucls_tumour, 0, 1)
    p_non_tumour_epithelial = np.clip(
        (0.50 * pan_epithelial + 0.30 * lizard_epithelial + 0.20 * panop_epithelial)
        * (1.0 - 0.50 * p_tumour),
        0,
        1,
    )
    p_fib = np.clip(
        0.45 * pan_connective + 0.30 * lizard_connective + 0.15 * panop_stromal + 0.10 * consep_spindle,
        0,
        1,
    )
    p_lymph = np.clip(0.55 * lizard_lymph + 0.25 * panop_til + 0.20 * pan_inflammatory, 0, 1)
    p_plasma = lizard_plasma
    p_neut = lizard_neut
    p_eos = lizard_eos
    p_mac = nucls_mac
    p_dead = pan_dead
    p_spindle = consep_spindle
    p_mitotic_tumour = np.minimum(p_tumour, midog_mitotic)

    p_immune = np.maximum.reduce([p_lymph, p_plasma, p_neut, p_eos, 0.75 * p_mac])
    p_epithelial = np.maximum(p_tumour, p_non_tumour_epithelial)
    p_other = p_dead
    broad = np.column_stack([p_epithelial, p_fib, p_immune, p_other])
    broad_id_raw = broad.argmax(axis=1)
    broad_conf = broad.max(axis=1).astype(np.float32)
    level1_id = np.where(broad_conf >= FROZEN.level1_probability_min, broad_id_raw, 4).astype(np.uint8)

    category_evidence = np.column_stack(
        [
            np.maximum.reduce([pan_neoplastic, pan_epithelial, ocelot_tumour, lizard_epithelial, panop_epithelial]),
            np.maximum.reduce([pan_connective, lizard_connective, panop_stromal, consep_spindle]),
            np.maximum.reduce([pan_inflammatory, lizard_lymph, lizard_plasma, lizard_neut, lizard_eos, panop_til]),
            pan_dead,
        ]
    )
    # Discordance means two incompatible biological categories have strong
    # evidence, not merely that two heads agree on the same category.
    strong_vote_count = (category_evidence >= FROZEN.level1_probability_min).sum(axis=1)
    discordant = strong_vote_count >= 2

    fine = np.column_stack(
        [
            p_tumour,
            p_non_tumour_epithelial,
            p_fib,
            p_lymph,
            p_plasma,
            p_neut,
            p_eos,
            p_dead,
            p_spindle,
            p_mac,
            p_mitotic_tumour,
        ]
    )
    fine_id = fine.argmax(axis=1).astype(np.uint8)
    fine_conf = fine.max(axis=1).astype(np.float32)
    level2_id = np.where(fine_conf >= FROZEN.level2_probability_min, fine_id, 11).astype(np.uint8)

    mitotic_gate = (p_mitotic_tumour >= FROZEN.level2_probability_min) & (
        p_tumour >= FROZEN.level2_probability_min
    )
    level2_id[mitotic_gate] = 10

    # Out-of-domain candidate heads cannot silently replace primary CRC labels.
    candidate_ids = np.isin(level2_id, [8, 9, 10])
    consensus_status = np.full(n, 0, dtype=np.uint8)  # 0 resolved, 1 unresolved, 2 discordant, 3 candidate
    consensus_status[level1_id == 4] = 1
    consensus_status[discordant] = 2
    consensus_status[candidate_ids & ~discordant] = 3

    return {
        "level1_class_id": level1_id,
        "level1_confidence": broad_conf,
        "level2_class_id": level2_id,
        "level2_confidence": fine_conf,
        "consensus_status": consensus_status,
        "p_tumour": p_tumour,
        "p_non_tumour_epithelial": p_non_tumour_epithelial,
        "p_fibroblast_like": p_fib,
        "p_lymphocyte": p_lymph,
        "p_plasma": p_plasma,
        "p_neutrophil": p_neut,
        "p_eosinophil": p_eos,
        "p_macrophage_candidate": p_mac,
        "p_dead": p_dead,
        "p_spindle_candidate": p_spindle,
        "p_mitotic_tumour_candidate": p_mitotic_tumour,
        "candidate_spindle_flag": (p_spindle >= FROZEN.level2_probability_min).astype(np.uint8),
        "candidate_macrophage_flag": (p_mac >= FROZEN.level2_probability_min).astype(np.uint8),
        "candidate_mitotic_flag": mitotic_gate.astype(np.uint8),
    }


def main() -> None:
    args = parse_args()
    if args.output.exists() and not args.overwrite:
        print(f"Already complete: {args.output}")
        return
    roi_manifest = json.loads(args.roi_record.read_text(encoding="utf-8"))
    records = {record["roi_uid"]: record for record in roi_manifest["rois"]}
    if args.roi_uid not in records:
        raise KeyError(f"Unknown ROI uid {args.roi_uid}")
    roi = records[args.roi_uid]

    tissue_manifest = None
    tissue_mask = None
    if args.tissue_mask_manifest is not None:
        tissue_mask, tissue_manifest = load_tissue_mask(args.tissue_mask_manifest)

    x0_um, y0_um, _, _ = map(float, roi["read_bbox_um"])
    mpp = float(roi["target_mpp"])
    with h5py.File(args.graph, "r") as graph:
        if int(graph.attrs.get("complete", 0)) != 1:
            raise RuntimeError(f"Incomplete graph: {args.graph}")
        n_source_cells, embedding_dim = graph["x"].shape
        positions_source = graph["positions"][:].astype(np.float32)
        if len(positions_source) != n_source_cells:
            raise RuntimeError("Graph positions and embeddings differ in length")
        source_cell_index = np.arange(n_source_cells, dtype=np.int64)
        if tissue_mask is not None:
            centroid_source_um = positions_source * mpp + np.asarray(
                [x0_um, y0_um], dtype=np.float32
            )
            inside = sample_points(
                tissue_mask,
                float(tissue_manifest["mask_mpp"]),
                centroid_source_um,
            )
            source_cell_index = np.flatnonzero(inside).astype(np.int64)
        positions = positions_source[source_cell_index]
        n_cells = int(len(source_cell_index))
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        models, head_labels = {}, {}
        for head, filename in CLASSIFIER_FILES.items():
            model, labels = load_classifier(CLASSIFIER_ROOT / filename, device, head_name=head)
            if int(model.fc1.in_features) != embedding_dim:
                raise RuntimeError(f"Embedding mismatch for {head}")
            models[head], head_labels[head] = model, labels
        head_prob = {
            head: np.empty((n_cells, len(labels)), dtype=np.float16)
            for head, labels in head_labels.items()
        }
        n_batches = int(math.ceil(n_cells / args.batch_size)) if n_cells else 0
        with torch.inference_mode():
            for batch_index, start in enumerate(range(0, n_cells, args.batch_size), start=1):
                end = min(start + args.batch_size, n_cells)
                token_rows = read_indexed_rows_dense(
                    graph["x"],
                    source_cell_index[start:end],
                    source_span_rows=max(65536, args.batch_size * 2),
                )
                tokens = torch.from_numpy(token_rows).to(
                    device=device, dtype=torch.float32
                )
                for head, model in models.items():
                    head_prob[head][start:end] = torch.softmax(model(tokens), dim=1).cpu().numpy().astype(np.float16)
                print(
                    f"MULTIHEAD_GPU_PROGRESS batch={batch_index}/{n_batches} "
                    f"cells={end}/{n_cells}",
                    flush=True,
                )

    if device.type == "cuda":
        if n_cells:
            del tokens
        del models
        torch.cuda.empty_cache()
    print(f"MULTIHEAD_GPU_COMPLETE cells={n_cells}", flush=True)

    cell_payload = ujson.loads(snappy.decompress(args.cells_json_snappy.read_bytes()))
    source_cells = cell_payload["cells"]
    if len(source_cells) != n_source_cells:
        raise RuntimeError(f"Cell/token mismatch: {len(source_cells)} vs {n_source_cells}")
    cells = [source_cells[int(index)] for index in source_cell_index]
    bbox_local = np.empty((n_cells, 2, 2), dtype=np.int32)
    pannuke_type = np.empty(n_cells, dtype=np.uint8)
    pannuke_conf = np.empty(n_cells, dtype=np.float32)
    morphology = np.empty((n_cells, len(MORPHOLOGY_FEATURES)), dtype=np.float32)
    contour_offsets = np.zeros(n_cells + 1, dtype=np.int64)
    contour_parts = []
    morphology_report_every = max(1, min(100000, int(math.ceil(n_cells / 10))))
    for index, cell in enumerate(cells):
        bbox = np.asarray(cell["bbox"], dtype=np.int32)
        contour = np.asarray(cell["contour"], dtype=np.int32)
        bbox_local[index] = bbox
        pannuke_type[index] = int(cell["type"])
        pannuke_conf[index] = float(cell.get("type_prob", 0.0))
        morphology[index] = contour_features(contour, bbox)
        contour_parts.append(contour)
        contour_offsets[index + 1] = contour_offsets[index] + len(contour)
        if (index + 1) % morphology_report_every == 0 or index + 1 == n_cells:
            print(
                f"MORPHOLOGY_PROGRESS {index + 1}/{n_cells} "
                f"({100.0 * (index + 1) / max(n_cells, 1):.1f}%)",
                flush=True,
            )
    contour_local = np.concatenate(contour_parts, axis=0) if contour_parts else np.empty((0, 2), dtype=np.int32)

    centroid_global_um = positions * mpp + np.asarray([x0_um, y0_um], dtype=np.float32)
    contour_global_um = contour_local.astype(np.float32) * mpp + np.asarray([x0_um, y0_um], dtype=np.float32)
    geometry = roi_manifest["slide_geometry"]
    centroid_native_px = np.column_stack(
        [
            centroid_global_um[:, 0] / float(geometry["native_mpp_x"]),
            centroid_global_um[:, 1] / float(geometry["native_mpp_y"]),
        ]
    ).astype(np.float32)
    derived = consensus(pannuke_type, pannuke_conf, head_prob, head_labels)
    low_confidence = (
        (derived["level1_class_id"] == 4)
        | (derived["level2_class_id"] == 11)
        | (derived["consensus_status"] != 0)
    )
    rng = np.random.default_rng(FROZEN.seed)
    qc_size = min(n_cells, max(100, min(2000, int(math.ceil(0.01 * n_cells)))))
    qc_index = (
        rng.choice(n_cells, size=qc_size, replace=False).astype(np.int64)
        if qc_size
        else np.empty(0, dtype=np.int64)
    )
    review_index = np.union1d(np.flatnonzero(low_confidence), qc_index).astype(np.int64)
    review_flags = np.zeros(review_index.size, dtype=np.uint8)
    review_flags[np.isin(review_index, qc_index)] |= 1
    review_flags[low_confidence[review_index]] |= 2
    review_flags[derived["consensus_status"][review_index] == 2] |= 4
    review_flags[derived["consensus_status"][review_index] == 3] |= 8
    with h5py.File(args.graph, "r") as graph:
        review_embeddings = read_indexed_rows_dense(
            graph["x"],
            source_cell_index[review_index],
            source_span_rows=65536,
            output_dtype=np.float16,
            progress_label="REVIEW_EMBEDDING_PROGRESS",
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".incomplete")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["schema"] = "cellular_interdigitating_topology_nuclei_v3"
        handle.attrs["complete"] = 0
        handle.attrs["roi_uid"] = args.roi_uid
        handle.attrs["roi_record_json"] = json.dumps(roi)
        handle.attrs["slide_geometry_json"] = json.dumps(geometry)
        handle.attrs["level1_labels_json"] = json.dumps(LEVEL1_LABELS)
        handle.attrs["level2_labels_json"] = json.dumps(LEVEL2_LABELS)
        handle.attrs["consensus_status_json"] = json.dumps({0: "resolved", 1: "unresolved", 2: "discordant", 3: "candidate"})
        handle.attrs["morphology_features_json"] = json.dumps(MORPHOLOGY_FEATURES)
        handle.attrs["graph_source"] = str(args.graph)
        handle.attrs["cell_source"] = str(args.cells_json_snappy)
        handle.attrs["n_source_cells_before_tissue_filter"] = int(n_source_cells)
        handle.attrs["n_cells_after_tissue_filter"] = int(n_cells)
        handle.attrs["whole_tissue_cellvit"] = int(roi.get("whole_tissue_cellvit", False))
        if tissue_manifest is not None:
            handle.attrs["tissue_first"] = 1
            handle.attrs["tissue_mask_manifest"] = str(args.tissue_mask_manifest)
            handle.attrs["tissue_mask_sha256"] = str(tissue_manifest["mask_sha256"])
        nuclei = handle.create_group("nuclei")
        nuclei.create_dataset("cell_id", data=np.arange(n_cells, dtype=np.int64), compression="lzf")
        nuclei.create_dataset("source_cell_id", data=source_cell_index, compression="lzf")
        nuclei.create_dataset("centroid_local_px", data=positions, compression="lzf")
        nuclei.create_dataset("centroid_global_um", data=centroid_global_um, compression="lzf")
        nuclei.create_dataset("centroid_native_px", data=centroid_native_px, compression="lzf")
        nuclei.create_dataset("bbox_local_px", data=bbox_local, compression="lzf")
        nuclei.create_dataset("contour_offsets", data=contour_offsets, compression="lzf")
        nuclei.create_dataset("contour_local_px", data=contour_local, compression="lzf")
        nuclei.create_dataset("contour_global_um", data=contour_global_um, compression="lzf")
        nuclei.create_dataset("pannuke_class_id", data=pannuke_type, compression="lzf")
        nuclei.create_dataset("pannuke_confidence", data=pannuke_conf, compression="lzf")
        nuclei.create_dataset("morphology", data=morphology, compression="lzf")
        for name, array in derived.items():
            nuclei.create_dataset(name, data=array, compression="lzf")
        probability_group = handle.create_group("probabilities")
        for head, probabilities in head_prob.items():
            dataset = probability_group.create_dataset(head, data=probabilities, compression="lzf")
            dataset.attrs["labels_json"] = json.dumps(head_labels[head])
            dataset.attrs["checkpoint"] = str(CLASSIFIER_ROOT / CLASSIFIER_FILES[head])
            dataset.attrs["checkpoint_sha256"] = sha256_file(CLASSIFIER_ROOT / CLASSIFIER_FILES[head])
        review = handle.create_group("embedding_review")
        review.create_dataset("cell_id", data=review_index, compression="lzf")
        review.create_dataset("reason_flags", data=review_flags, compression="lzf")
        review.create_dataset("embedding_float16", data=review_embeddings, compression="lzf")
        review.attrs["reason_flags_json"] = json.dumps(
            {1: "deterministic_qc_sample", 2: "low_confidence_or_unresolved", 4: "discordant", 8: "candidate"}
        )
        review.attrs["full_embedding_persisted"] = 0
        handle.attrs["complete"] = 1
        handle.flush()
    os.replace(temporary, args.output)

    summary = {
        "status": "complete",
        "output": str(args.output),
        "roi_uid": args.roi_uid,
        "n_cells": n_cells,
        "n_source_cells_before_tissue_filter": int(n_source_cells),
        "n_cells_removed_outside_explicit_tissue": int(n_source_cells - n_cells),
        "embedding_dim": embedding_dim,
        "n_review_embeddings": int(review_index.size),
        "n_low_confidence_review_embeddings": int(low_confidence.sum()),
        "full_embedding_persisted": False,
        "level1_counts": {
            LEVEL1_LABELS[int(key)]: int(value)
            for key, value in zip(*np.unique(derived["level1_class_id"], return_counts=True))
        },
        "level2_counts": {
            LEVEL2_LABELS[int(key)]: int(value)
            for key, value in zip(*np.unique(derived["level2_class_id"], return_counts=True))
        },
    }
    atomic_json(args.output.with_suffix(".summary.json"), summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
