#!/usr/bin/env python3
"""Run one persistent CellViT++ model across arbitrary CRC whole-slide images.

The explicit 8-um tissue mask is converted to a QuPath-compatible GeoJSON and
is supplied to PathoPatch before window creation. CellViT++ therefore sees
the complete pathology tissue domain once per WSI, without reloading the model
for each UNI front piece.  Nuclei outside the exact frozen raster mask are
removed again by ``build_multihead_nuclei.py`` before the final HDF5 is sealed.
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable

import cv2
import h5py
import numpy as np
import ray
import torch

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from config import CELLVIT_ROOT

sys.path.insert(0, str(CELLVIT_ROOT))

from cellvit.data.dataclass.cell_graph import CellGraphDataWSI
from cellvit.inference.inference_disk import CellViTInference
from cellvit.inference.inference_memory import CellViTInferenceMemory
import pathopatch.patch_extraction.dataset as patch_dataset

from common import atomic_json, sha256_file
from config import CELLVIT_CHECKPOINT, FROZEN, OUTPUT_ROOT


SCHEMA = "cellvitpp_whole_tissue_coad_batch_v3"
CELLVIT_MAX_SOURCE_MPP = 0.75
TERMINAL = {"complete", "no_nuclei_detected", "excluded_low_resolution"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument(
        "--status",
        type=Path,
        default=OUTPUT_ROOT / "batch/coad_cellvitpp_whole_tissue_status.json",
    )
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--head-batch-size", type=int, default=32768)
    parser.add_argument("--tissue-workers", type=int, default=8)
    parser.add_argument("--sample-id", action="append")
    parser.add_argument("--max-slides", type=int)
    parser.add_argument("--retry-failed", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--keep-full-embedding-graph", action="store_true")
    parser.add_argument("--minimum-free-gb", type=float, default=20.0)
    return parser.parse_args()


def tissue_manifest_path(record: dict) -> Path:
    return OUTPUT_ROOT / "tissue_masks" / record["cohort"] / record["sample_id"] / "tissue_mask_manifest.json"


def final_nuclei_path(record: dict) -> Path:
    return OUTPUT_ROOT / "nuclei_whole_tissue" / f"{record['sample_id']}_nuclei.h5"


def sample_complete_manifest(record: dict) -> Path:
    return OUTPUT_ROOT / "cellvit_whole_tissue" / record["cohort"] / record["sample_id"] / "complete.json"


def source_mpp_values(record: dict) -> list[float]:
    """Return finite positive source MPP values recorded for a slide."""
    values = []
    for key in ("mpp_x", "mpp_y"):
        value = record.get(key)
        if value is None:
            continue
        numeric = float(value)
        if np.isfinite(numeric) and numeric > 0:
            values.append(numeric)
    return values


def requires_low_resolution_exclusion(
    record: dict, max_source_mpp: float = CELLVIT_MAX_SOURCE_MPP
) -> bool:
    """Identify slides whose native sampling cannot support CellViT nuclei inference.

    Upsampling a low-resolution source does not recreate nuclear detail, so these
    slides must be excluded rather than relabelled with a fabricated MPP.
    """
    values = source_mpp_values(record)
    return bool(values) and max(values) > float(max_source_mpp)


def low_resolution_exclusion_result(record: dict, tissue_path: Path, tissue: dict) -> dict:
    values = source_mpp_values(record)
    observed = max(values)
    return {
        "schema": "cellvitpp_whole_tissue_slide_v3",
        "status": "excluded_low_resolution",
        "stage": "source_resolution_preflight",
        "sample_id": record["sample_id"],
        "case_id": record["case_id"],
        "source_slide": str(record["slide"]),
        "tissue_first": True,
        "whole_tissue_cellvit": True,
        "tissue_mask_manifest": str(tissue_path),
        "tissue_mask_sha256": tissue["mask_sha256"],
        "source_mpp_x": record.get("mpp_x"),
        "source_mpp_y": record.get("mpp_y"),
        "observed_max_source_mpp": observed,
        "cellvit_max_source_mpp": CELLVIT_MAX_SOURCE_MPP,
        "exclusion_rule": "max(mpp_x, mpp_y) > cellvit_max_source_mpp",
        "requested_cellvit_target_mpp": float(FROZEN.cellvit_target_mpp),
        "reason_code": "native_source_mpp_above_cellvit_limit",
        "reason": (
            f"Native source MPP {observed:.4f} exceeds the CellViT++ input limit "
            f"of {CELLVIT_MAX_SOURCE_MPP:.2f} um/px."
        ),
        "scientific_rationale": (
            "Upsampling a low-resolution source cannot recover missing nuclear morphology; "
            "overriding MPP would create invalid cell-level measurements."
        ),
        "inference_performed": False,
        "silent_slide_substitution_performed": False,
        "formal_eligible": bool(record.get("formal_eligible", False)),
        "sensitivity_only": bool(record.get("sensitivity_only", False)),
    }


def is_complete_h5(path: Path) -> bool:
    if not path.exists():
        return False
    try:
        with h5py.File(path, "r") as handle:
            return int(handle.attrs.get("complete", 0)) == 1
    except Exception:
        return False


def load_complete_tissue_manifest(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("status") != "complete" or payload.get("tissue_first") is not True:
        raise RuntimeError(f"Incomplete tissue-first manifest: {path}")
    mask_path = Path(payload["mask_npz"])
    if not mask_path.exists() or sha256_file(mask_path) != payload["mask_sha256"]:
        raise RuntimeError(f"Tissue mask is missing or has a SHA256 mismatch: {mask_path}")
    return payload


def build_one_tissue_mask(record: dict) -> tuple[str, str]:
    destination = tissue_manifest_path(record)
    try:
        if destination.exists():
            load_complete_tissue_manifest(destination)
            return record["sample_id"], "existing_complete"
        command = [
            sys.executable,
            str(SCRIPT_DIR / "build_tissue_mask.py"),
            "--cohort",
            record["cohort"],
            "--slide",
            record["slide"],
            "--sample-id",
            record["sample_id"],
            "--case-id",
            record["case_id"],
        ]
        if record.get("mpp_x") is not None:
            command.extend(["--mpp-x", str(record["mpp_x"]), "--mpp-y", str(record["mpp_y"])])
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        load_complete_tissue_manifest(destination)
        return record["sample_id"], "built_complete"
    except Exception as error:
        return record["sample_id"], f"failed: {type(error).__name__}: {error}"


def prepare_tissue_masks(records: Iterable[dict], workers: int) -> Dict[str, str]:
    records = list(records)
    results: Dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = {pool.submit(build_one_tissue_mask, record): record for record in records}
        for completed, future in enumerate(as_completed(futures), start=1):
            sample_id, status = future.result()
            results[sample_id] = status
            if completed == 1 or completed % 25 == 0 or completed == len(records):
                failures = sum(value.startswith("failed") for value in results.values())
                print(
                    f"TISSUE_MASK_PROGRESS {completed}/{len(records)} failures={failures}",
                    flush=True,
                )
    return results


def close_ring(points: np.ndarray) -> list[list[float]]:
    coordinates = points.reshape(-1, 2).astype(float).tolist()
    if coordinates and coordinates[0] != coordinates[-1]:
        coordinates.append(coordinates[0])
    return coordinates


def export_tissue_geojson(tissue_manifest: dict, destination: Path) -> dict:
    mask_path = Path(tissue_manifest["mask_npz"])
    with np.load(mask_path) as bundle:
        mask = bundle["mask"].astype(np.uint8)
        mask_mpp = float(bundle["mpp"][0])
    geometry = tissue_manifest["slide_geometry"]
    scale = np.asarray(
        [mask_mpp / float(geometry["native_mpp_x"]), mask_mpp / float(geometry["native_mpp_y"])],
        dtype=np.float64,
    )
    contours, hierarchy = cv2.findContours(mask, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
    if hierarchy is None:
        raise RuntimeError("Explicit tissue mask contains no contour")
    hierarchy = hierarchy[0]
    features = []
    for index, contour in enumerate(contours):
        if hierarchy[index][3] != -1 or len(contour) < 3:
            continue
        exterior = cv2.approxPolyDP(contour, epsilon=0.5, closed=True).reshape(-1, 2) * scale
        if len(exterior) < 3:
            continue
        rings = [close_ring(exterior)]
        child = int(hierarchy[index][2])
        while child != -1:
            hole = cv2.approxPolyDP(contours[child], epsilon=0.5, closed=True).reshape(-1, 2) * scale
            if len(hole) >= 3:
                rings.append(close_ring(hole))
            child = int(hierarchy[child][0])
        features.append(
            {
                "type": "Feature",
                "geometry": {"type": "Polygon", "coordinates": rings},
                "properties": {"classification": {"name": "tissue"}},
            }
        )
    if not features:
        raise RuntimeError("No valid tissue polygon could be exported")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".incomplete")
    temporary.write_text(json.dumps(features, separators=(",", ":")), encoding="utf-8")
    os.replace(temporary, destination)
    return {"n_polygons": len(features), "geojson_sha256": sha256_file(destination)}


def safe_livepatch_next(loader):
    """Return the next non-empty live WSI batch.

    PathoPatch 1.0.2 indexes ``patches[0]`` when all remaining candidates are
    rejected by its final RGB-background check.  Exhausting only rejected
    candidates is a normal iterator terminal state, not an inference error.
    The loop also remains correct if a future PathoPatch version limits a
    single scan to fewer candidates: an empty scan advances to the next one.
    """

    while loader.i < len(loader.element_list):
        patches = []
        metadata = []
        masks = []
        batch_item_count = 0
        while batch_item_count < loader.batch_size and loader.i < len(loader.element_list):
            patch, meta, mask = loader.dataset[loader.element_list[loader.i]]
            loader.i += 1
            if patch is None and meta["discard_patch"]:
                loader.discard_count += 1
                continue
            if loader.dataset.config.filter_patches:
                output = loader.dataset.detector_model(
                    loader.dataset.detector_transforms(patch)[None, ...]
                )
                output_prob = torch.softmax(output, dim=-1)
                prediction = torch.argmax(output_prob, dim=-1)
                if int(prediction) != 0:
                    loader.discard_count += 1
                    continue
            patches.append(patch)
            metadata.append(meta)
            masks.append(mask)
            batch_item_count += 1
        if not patches:
            continue
        if len(patches) > 1:
            patch_tensor = torch.stack([torch.tensor(value) for value in patches])
        else:
            patch_tensor = torch.tensor(patches[0][None, ...])
        return patch_tensor, metadata, masks
    raise StopIteration


def install_upstream_safety_patches() -> None:
    original_module_exists = patch_dataset.module_exists

    def module_exists(module_name: str, error: str = "ignore") -> bool:
        if module_name == "cucim":
            return False
        return original_module_exists(module_name, error=error)

    patch_dataset.module_exists = module_exists

    # PathoPatch 1.0.2 calls ``self.annotation_path.suffix`` although the
    # Pydantic field is declared as ``str``.  It also discards a directly
    # supplied label map whenever no label-map file is given.  Correct those
    # two upstream issues locally so the frozen tissue annotation can be used.
    def fixed_config_post_init(self) -> None:
        if self.label_map_file is None and self.label_map is None:
            self.label_map = {"background": 0}
        if self.otsu_annotation is not None:
            self.otsu_annotation = self.otsu_annotation.lower()
        if self.tissue_annotation is not None:
            self.tissue_annotation = self.tissue_annotation.lower()
        if len(self.exclude_classes) > 0:
            self.exclude_classes = [value.lower() for value in self.exclude_classes]
        if self.tissue_annotation_intersection_ratio is None:
            self.tissue_annotation_intersection_ratio = self.min_intersection_ratio
        elif not 0 <= self.tissue_annotation_intersection_ratio <= 1:
            raise RuntimeError("Tissue annotation intersection ratio must be between 0 and 1")
        if self.annotation_path is not None:
            annotation_path = Path(self.annotation_path)
            if not annotation_path.exists():
                raise FileNotFoundError(f"Annotation path {annotation_path} does not exist")
            if annotation_path.suffix.lower() != ".json":
                raise ValueError("Only JSON annotations are supported")

    patch_dataset.LivePatchWSIConfig.__post_init_post_parse__ = fixed_config_post_init
    patch_dataset.LivePatchWSIDataloader.__next__ = safe_livepatch_next

    original_setup_worker = CellViTInference._setup_worker

    def setup_worker(self) -> None:
        original_setup_worker(self)
        schedulable = max(1, int(ray.cluster_resources().get("CPU", 1) // 8))
        if self.ray_actors > schedulable:
            self.logger.info(
                "Adjusting Ray post-processing actors from %d to %d",
                self.ray_actors,
                schedulable,
            )
            self.ray_actors = schedulable

    CellViTInference._setup_worker = setup_worker

    original_torch_save = torch.save

    def safe_graph_save(obj, target, *args, **kwargs):
        if isinstance(obj, CellGraphDataWSI) and str(target).endswith("_cells.pt"):
            target_path = Path(target).with_suffix(".h5")
            temporary_path = target_path.with_suffix(".h5.incomplete")
            x = obj.x.detach().cpu()
            positions = obj.positions.detach().cpu()
            n_cells = int(x.shape[0])
            chunk_rows = min(4096, max(1, n_cells))
            print(f"WRITE_GRAPH_H5 {target_path} cells={n_cells}", flush=True)
            with h5py.File(temporary_path, "w") as handle:
                handle.attrs["schema"] = "cellvit_graph_h5_v1"
                handle.attrs["metadata_json"] = json.dumps(obj.metadata, default=str)
                x_dataset = handle.create_dataset(
                    "x", shape=tuple(x.shape), dtype=x.numpy().dtype, chunks=(chunk_rows, int(x.shape[1]))
                )
                position_dataset = handle.create_dataset(
                    "positions",
                    shape=tuple(positions.shape),
                    dtype=positions.numpy().dtype,
                    chunks=(chunk_rows, int(positions.shape[1])),
                )
                for start in range(0, n_cells, chunk_rows):
                    end = min(start + chunk_rows, n_cells)
                    x_dataset[start:end] = x[start:end].numpy()
                    position_dataset[start:end] = positions[start:end].numpy()
                handle.attrs["complete"] = 1
                handle.flush()
            os.replace(temporary_path, target_path)
            return None
        return original_torch_save(obj, target, *args, **kwargs)

    torch.save = safe_graph_save


def update_status(path: Path, payload: dict, total: int) -> None:
    counts: Dict[str, int] = {}
    for record in payload["records"].values():
        status = record.get("status", "unknown")
        counts[status] = counts.get(status, 0) + 1
    payload["summary"] = {
        "n_manifest": total,
        "n_status_records": len(payload["records"]),
        "status_counts": counts,
        "n_complete": counts.get("complete", 0),
        "n_no_nuclei_detected": counts.get("no_nuclei_detected", 0),
        "n_excluded_low_resolution": counts.get("excluded_low_resolution", 0),
        "n_failed": counts.get("failed", 0),
        "n_running": counts.get("running", 0),
        "n_terminal": sum(counts.get(status, 0) for status in TERMINAL),
        "updated_at_epoch": time.time(),
    }
    atomic_json(path, payload)


def graph_metadata(path: Path) -> dict:
    with h5py.File(path, "r") as handle:
        if int(handle.attrs.get("complete", 0)) != 1:
            raise RuntimeError(f"Incomplete CellViT++ graph: {path}")
        return json.loads(handle.attrs["metadata_json"])


def main() -> None:
    args = parse_args()
    if not CELLVIT_CHECKPOINT.exists():
        raise FileNotFoundError(CELLVIT_CHECKPOINT)
    source = json.loads(args.manifest.read_text(encoding="utf-8"))
    records = list(source["records"])
    if args.sample_id:
        wanted = set(args.sample_id)
        records = [record for record in records if record["sample_id"] in wanted]
        missing = wanted - {record["sample_id"] for record in records}
        if missing:
            raise ValueError(f"Samples absent from manifest: {sorted(missing)}")
    if args.max_slides is not None:
        records = records[: args.max_slides]
    if not records:
        raise RuntimeError("No slides selected")

    status_payload = (
        json.loads(args.status.read_text(encoding="utf-8"))
        if args.status.exists()
        else {
            "schema": SCHEMA,
            "manifest": str(args.manifest),
            "cellvit_checkpoint": str(CELLVIT_CHECKPOINT),
            "cellvit_checkpoint_sha256": sha256_file(CELLVIT_CHECKPOINT),
            "tissue_first": True,
            "whole_tissue_cellvit": True,
            "records": {},
        }
    )
    if status_payload.get("schema") != SCHEMA:
        raise RuntimeError("Refusing to mix a non-whole-tissue batch status file")
    status_payload["runner_pid"] = os.getpid()
    status_payload["runner_started_at_epoch"] = time.time()
    update_status(args.status, status_payload, len(source["records"]))

    print(
        f"PREPARE_TISSUE_MASKS slides={len(records)} workers={args.tissue_workers}",
        flush=True,
    )
    tissue_results = prepare_tissue_masks(records, args.tissue_workers)
    failures = {key: value for key, value in tissue_results.items() if value.startswith("failed")}
    if failures:
        for sample_id, error in failures.items():
            status_payload["records"][sample_id] = {"status": "failed", "stage": "tissue_mask", "error": error}
        update_status(args.status, status_payload, len(source["records"]))
        records = [record for record in records if record["sample_id"] not in failures]
    if not records:
        raise RuntimeError("No slide has a valid explicit tissue mask")

    # CellViT++ refuses native MPP above 0.75 um/px. Treat genuine low-resolution
    # slides as a documented terminal data exclusion, not a retryable model error.
    # This preflight happens before model construction so a low-resolution-only
    # targeted repair does not allocate GPU memory unnecessarily.
    processable_records = []
    for selected_ordinal, record in enumerate(records, start=1):
        if not requires_low_resolution_exclusion(record):
            processable_records.append(record)
            continue
        sample_id = record["sample_id"]
        complete_path = sample_complete_manifest(record)
        previous = status_payload["records"].get(sample_id, {})
        if previous.get("status") == "excluded_low_resolution" and complete_path.exists() and not args.overwrite:
            print(f"SKIP_EXCLUDED_LOW_RESOLUTION {selected_ordinal}/{len(records)} {sample_id}", flush=True)
            continue
        final_h5 = final_nuclei_path(record)
        sample_root = complete_path.parent
        raw_exists = any((sample_root / "pannuke").glob("*_cells.h5")) or any(
            (sample_root / "pannuke").glob("*_cells.json.snappy")
        )
        if final_h5.exists() or raw_exists:
            raise RuntimeError(
                f"Low-resolution slide {sample_id} has existing CellViT outputs; "
                "refusing to classify them as valid or delete them automatically"
            )
        tissue_path = tissue_manifest_path(record)
        tissue = load_complete_tissue_manifest(tissue_path)
        complete_path.parent.mkdir(parents=True, exist_ok=True)
        result = low_resolution_exclusion_result(record, tissue_path, tissue)
        result.update(
            {
                "ordinal": previous.get("ordinal", selected_ordinal),
                "selected_total": len(records),
                "excluded_at_epoch": time.time(),
            }
        )
        atomic_json(complete_path, result)
        state = dict(record)
        state.update(result)
        state["error"] = ""
        state.pop("traceback", None)
        state.pop("failed_at_epoch", None)
        status_payload["records"][sample_id] = state
        update_status(args.status, status_payload, len(source["records"]))
        print(
            f"EXCLUDED_LOW_RESOLUTION {selected_ordinal}/{len(records)} {sample_id} "
            f"mpp={result['observed_max_source_mpp']:.4f}",
            flush=True,
        )

    records = processable_records
    if not records:
        update_status(args.status, status_payload, len(source["records"]))
        print(json.dumps(status_payload["summary"], indent=2), flush=True)
        return

    install_upstream_safety_patches()
    initial_outdir = OUTPUT_ROOT / "cellvit_whole_tissue" / "_model_initialization"
    initial_outdir.mkdir(parents=True, exist_ok=True)
    detector = CellViTInferenceMemory(
        model_path=CELLVIT_CHECKPOINT,
        classifier_path=None,
        binary=False,
        gpu=0,
        outdir=initial_outdir,
        geojson=False,
        graph=True,
        compression=True,
        batch_size=args.batch_size,
        enforce_mixed_precision=True,
    )

    for ordinal, record in enumerate(records, start=1):
        sample_id = record["sample_id"]
        final_h5 = final_nuclei_path(record)
        complete_path = sample_complete_manifest(record)
        previous = status_payload["records"].get(sample_id, {})
        previous_sample_root = sample_complete_manifest(record).parent
        previous_raw_exists = any((previous_sample_root / "pannuke").glob("*_cells.h5")) or any(
            (previous_sample_root / "pannuke").glob("*_cells.json.snappy")
        )
        if (
            not args.overwrite
            and previous.get("status") in TERMINAL
            and complete_path.exists()
            and (
                is_complete_h5(final_h5)
                or (
                    previous.get("status") in {"no_nuclei_detected", "excluded_low_resolution"}
                    and not previous_raw_exists
                )
            )
        ):
            print(f"SKIP_COMPLETE {ordinal}/{len(records)} {sample_id}", flush=True)
            continue
        if previous.get("status") == "failed" and not args.retry_failed and not args.overwrite:
            print(f"SKIP_FAILED {ordinal}/{len(records)} {sample_id}", flush=True)
            continue
        free_gb = shutil.disk_usage(OUTPUT_ROOT).free / (1024**3)
        if free_gb < args.minimum_free_gb:
            raise RuntimeError(f"Free-space guard: {free_gb:.1f} GiB < {args.minimum_free_gb:.1f} GiB")

        started = time.time()
        state = dict(record)
        state.update(
            {
                "status": "running",
                "stage": "whole_tissue_cellvitpp",
                "ordinal": ordinal,
                "selected_total": len(records),
                "started_at_epoch": started,
                "error": "",
            }
        )
        status_payload["records"][sample_id] = state
        update_status(args.status, status_payload, len(source["records"]))
        try:
            tissue_path = tissue_manifest_path(record)
            tissue = load_complete_tissue_manifest(tissue_path)
            sample_root = complete_path.parent
            pannuke_dir = sample_root / "pannuke"
            pannuke_dir.mkdir(parents=True, exist_ok=True)
            # PathoPatch accepts QuPath GeoJSON content but requires a .json suffix.
            tissue_geojson = sample_root / "explicit_tissue_domain.json"
            geojson_metrics = export_tissue_geojson(tissue, tissue_geojson)

            slide = Path(record["slide"])
            # Upstream names outputs from everything before the first period,
            # not pathlib's suffix-stripped stem (TCGA SVS names contain a UUID
            # after the first period).
            stem = slide.name.split(".", 1)[0]
            graph = pannuke_dir / f"{stem}_cells.h5"
            cells = pannuke_dir / f"{stem}_cells.json.snappy"
            if args.overwrite:
                for stale in (graph, graph.with_suffix(".h5.incomplete"), cells, final_h5):
                    if stale.exists():
                        stale.unlink()
            if not (graph.exists() and cells.exists()):
                detector.outdir = pannuke_dir
                wsi_properties = {}
                if record.get("mpp_x") is not None:
                    wsi_properties = {"slide_mpp": float(record["mpp_x"]), "magnification": 40}
                print(
                    f"CELLVIT_WSI_START {ordinal}/{len(records)} {sample_id} "
                    f"tissue_mm2={tissue['metrics']['tissue_area_mm2']:.3f}",
                    flush=True,
                )
                detector.process_wsi(
                    wsi_path=slide,
                    wsi_properties=wsi_properties,
                    resolution=FROZEN.cellvit_target_mpp,
                    annotation_path=str(tissue_geojson),
                    label_map={"background": 0, "tissue": 1},
                    tissue_annotation="tissue",
                    tissue_annotation_intersection_ratio=FROZEN.tissue_candidate_fov_fraction_min,
                    min_intersection_ratio=FROZEN.tissue_candidate_fov_fraction_min,
                    apply_prefilter=False,
                    filter_patches=False,
                )
            else:
                print(f"RESUME_VERIFIED_RAW_CELLVIT {sample_id}", flush=True)
            if not graph.exists() and not cells.exists():
                result = {
                    "schema": "cellvitpp_whole_tissue_slide_v3",
                    "status": "no_nuclei_detected",
                    "sample_id": sample_id,
                    "source_slide": str(slide),
                    "tissue_first": True,
                    "tissue_mask_manifest": str(tissue_path),
                    "tissue_mask_sha256": tissue["mask_sha256"],
                    "elapsed_seconds": time.time() - started,
                }
                atomic_json(complete_path, result)
                state.update(result)
                update_status(args.status, status_payload, len(source["records"]))
                continue
            if not graph.exists() or not cells.exists():
                raise RuntimeError("CellViT++ wrote only one of graph/cell-table outputs")

            metadata = graph_metadata(graph)
            wsi_metadata = metadata["wsi_metadata"]
            effective_mpp = float(wsi_metadata["target_patch_mpp"])
            synthetic_record = sample_root / "whole_tissue_record.json"
            atomic_json(
                synthetic_record,
                {
                    "schema": "whole_tissue_cellvit_roi_record_v3",
                    "tissue_first": True,
                    "cohort": record["cohort"],
                    "sample_id": sample_id,
                    "slide_geometry": tissue["slide_geometry"],
                    "tissue_mask_manifest": str(tissue_path),
                    "tissue_mask_sha256": tissue["mask_sha256"],
                    "rois": [
                        {
                            "roi_uid": "whole_tissue",
                            "whole_tissue_cellvit": True,
                            "read_bbox_um": [
                                0.0,
                                0.0,
                                float(tissue["slide_geometry"]["width_um"]),
                                float(tissue["slide_geometry"]["height_um"]),
                            ],
                            "target_mpp": effective_mpp,
                            "image_path": str(slide),
                        }
                    ],
                },
            )
            state["stage"] = "whole_tissue_multihead_consensus"
            update_status(args.status, status_payload, len(source["records"]))
            command = [
                sys.executable,
                str(SCRIPT_DIR / "build_multihead_nuclei.py"),
                "--graph",
                str(graph),
                "--cells-json-snappy",
                str(cells),
                "--roi-record",
                str(synthetic_record),
                "--roi-uid",
                "whole_tissue",
                "--tissue-mask-manifest",
                str(tissue_path),
                "--output",
                str(final_h5),
                "--batch-size",
                str(args.head_batch_size),
            ]
            if args.overwrite:
                command.append("--overwrite")
            subprocess.run(command, check=True)
            if not is_complete_h5(final_h5):
                raise RuntimeError("Final whole-tissue nucleus HDF5 is not complete")
            summary_path = final_h5.with_suffix(".summary.json")
            nucleus_summary = json.loads(summary_path.read_text(encoding="utf-8"))
            result = {
                "schema": "cellvitpp_whole_tissue_slide_v3",
                "status": "complete",
                "sample_id": sample_id,
                "case_id": record["case_id"],
                "source_slide": str(slide),
                "tissue_first": True,
                "whole_tissue_cellvit": True,
                "cellvit_checkpoint": str(CELLVIT_CHECKPOINT),
                "cellvit_checkpoint_sha256": status_payload["cellvit_checkpoint_sha256"],
                "tissue_mask_manifest": str(tissue_path),
                "tissue_mask_sha256": tissue["mask_sha256"],
                "tissue_geojson": str(tissue_geojson),
                "tissue_geojson_metrics": geojson_metrics,
                "effective_cellvit_mpp": effective_mpp,
                "pannuke_cells": str(cells),
                "intermediate_embedding_graph": str(graph) if args.keep_full_embedding_graph else None,
                "final_nuclei_h5": str(final_h5),
                "nuclei_summary": nucleus_summary,
                "elapsed_seconds": time.time() - started,
            }
            atomic_json(complete_path, result)
            if not args.keep_full_embedding_graph and graph.exists():
                graph.unlink()
                result["intermediate_embedding_graph_deleted_after_verified_conversion"] = True
                atomic_json(complete_path, result)
            state.update(result)
            print(
                f"CELLVIT_WSI_COMPLETE {ordinal}/{len(records)} {sample_id} "
                f"cells={nucleus_summary['n_cells']} elapsed_min={(time.time()-started)/60:.1f}",
                flush=True,
            )
        except Exception as error:
            state.update(
                {
                    "status": "failed",
                    "stage": state.get("stage", "unknown"),
                    "error": f"{type(error).__name__}: {error}",
                    "traceback": traceback.format_exc(),
                    "failed_at_epoch": time.time(),
                }
            )
            print(f"CELLVIT_WSI_FAILED {ordinal}/{len(records)} {sample_id}: {state['error']}", flush=True)
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            gc.collect()
        update_status(args.status, status_payload, len(source["records"]))

    update_status(args.status, status_payload, len(source["records"]))
    print(json.dumps(status_payload["summary"], indent=2), flush=True)


if __name__ == "__main__":
    main()
