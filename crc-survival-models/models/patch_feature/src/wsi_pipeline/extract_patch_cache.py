from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
import cv2


SCRIPT_ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(SCRIPT_ROOT))
from slide_io import SlideGeometry, open_physical_slide, tissue_fraction  # noqa: E402
from uni_model import IMAGENET_MEAN, IMAGENET_STD, load_uni_classifier  # noqa: E402


UNI_CHECKPOINT = Path(
    os.environ.get(
        "CRC_SURVIVAL_UNI_CHECKPOINT",
        str(REPOSITORY_ROOT / "external/uni/patch_encoder_backbone.pth"),
    )
).expanduser().resolve()
SEED = 20260908
GRID_UM = 250.0
VIEW_MPP = 1.0
IMAGE_SIZE = 224
MIN_COARSE_TISSUE = 0.05
MIN_RGB_TISSUE = 0.03
MAX_REGIONS = 192
MAX_CANDIDATES_TO_DECODE = 4096


@dataclass
class PreparedSlide:
    row: object
    output: Path
    temporary: Path
    tensor: torch.Tensor
    xy: np.ndarray
    grid_ij: np.ndarray
    coarse_tissue: np.ndarray
    rgb_tissue: np.ndarray
    geometry: dict
    n_coarse_candidates: int
    n_candidates_decoded: int
    started: float
    io_seconds: float


class CuCIMSVSReader:
    """Fast SVS reader with the same physical-coordinate contract as slide_io."""

    def __init__(self, path: Path):
        from cucim import CuImage

        self.path = Path(path)
        self.slide = CuImage(str(path))
        metadata = self.slide.metadata
        aperio = metadata.get("aperio", {})
        mpp = float(aperio.get("MPP", "nan"))
        if not np.isfinite(mpp):
            raise ValueError(f"cuCIM could not resolve MPP for {path}")
        width, height = self.slide.resolutions["level_dimensions"][0]
        self.level_downsamples = tuple(float(value) for value in self.slide.resolutions["level_downsamples"])
        self.geometry = SlideGeometry(
            path=str(path),
            reader="cucim",
            native_mpp_x=mpp,
            native_mpp_y=mpp,
            native_origin_x_px=0,
            native_origin_y_px=0,
            width_px=int(width),
            height_px=int(height),
        )

    def read_bbox_um(self, bbox_um, target_mpp: float) -> np.ndarray:
        x0_um, y0_um, x1_um, y1_um = bbox_um
        x0_um = max(0.0, float(x0_um))
        y0_um = max(0.0, float(y0_um))
        x1_um = min(self.geometry.width_um, float(x1_um))
        y1_um = min(self.geometry.height_um, float(y1_um))
        out_w = max(1, int(round((x1_um - x0_um) / target_mpp)))
        out_h = max(1, int(round((y1_um - y0_um) / target_mpp)))
        desired_downsample = max(1.0, target_mpp / self.geometry.native_mpp_x)
        valid_levels = [index for index, value in enumerate(self.level_downsamples) if value <= desired_downsample]
        level = max(valid_levels) if valid_levels else 0
        downsample = self.level_downsamples[level]
        level_w = max(1, int(math.ceil((x1_um - x0_um) / (self.geometry.native_mpp_x * downsample))))
        level_h = max(1, int(math.ceil((y1_um - y0_um) / (self.geometry.native_mpp_y * downsample))))
        location = (int(math.floor(x0_um / self.geometry.native_mpp_x)), int(math.floor(y0_um / self.geometry.native_mpp_y)))
        region = np.asarray(self.slide.read_region(location=location, size=(level_w, level_h), level=level))
        rgb = np.asarray(region[..., :3], dtype=np.uint8)
        if rgb.shape[:2] != (out_h, out_w):
            rgb = cv2.resize(rgb, (out_w, out_h), interpolation=cv2.INTER_AREA)
        return rgb

    def coarse_tissue_presence(self, target_mpp: float = 8.0) -> np.ndarray:
        rgb = self.read_bbox_um((0.0, 0.0, self.geometry.width_um, self.geometry.height_um), target_mpp)
        array = rgb.astype(np.float32)
        brightness = array.mean(axis=2)
        chroma = array.max(axis=2) - array.min(axis=2)
        return (brightness < 238.0) & (chroma > 6.0)

    def read_center_patches(
        self,
        centers_um: np.ndarray,
        target_mpp: float,
        size: int,
        num_workers: int,
        batch_size: int = 32,
    ) -> np.ndarray:
        """Ordered multi-thread reads; executor.map preserves centre-to-patch identity."""
        def read_one(center) -> np.ndarray:
            cx, cy = map(float, center)
            fov = float(target_mpp) * int(size)
            return self.read_bbox_um(
                (cx - fov / 2.0, cy - fov / 2.0, cx + fov / 2.0, cy + fov / 2.0),
                target_mpp,
            )

        with ThreadPoolExecutor(max_workers=max(1, int(num_workers))) as pool:
            return np.stack(list(pool.map(read_one, centers_um)))

    def close(self) -> None:
        self.slide = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()


def open_extraction_slide(path: Path):
    if path.suffix.lower() == ".svs":
        try:
            return CuCIMSVSReader(path)
        except Exception:
            return open_physical_slide(path)
    return open_physical_slide(path)


def stable_seed(text: str) -> int:
    digest = hashlib.sha256(f"{SEED}|{text}".encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "little")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(8 * 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def atomic_tsv(path: Path, table: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    table.to_csv(temporary, sep="\t", index=False)
    os.replace(temporary, path)


def integral_image(mask: np.ndarray) -> np.ndarray:
    return np.pad(mask.astype(np.int64), ((1, 0), (1, 0))).cumsum(0).cumsum(1)


def rect_sum(integral: np.ndarray, x0: int, y0: int, x1: int, y1: int) -> int:
    return int(integral[y1, x1] - integral[y0, x1] - integral[y1, x0] + integral[y0, x0])


def candidate_grid(reader, coarse_mpp: float = 8.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mask = reader.coarse_tissue_presence(target_mpp=coarse_mpp)
    integral = integral_image(mask)
    height, width = mask.shape
    nx = int(math.floor(reader.geometry.width_um / GRID_UM))
    ny = int(math.floor(reader.geometry.height_um / GRID_UM))
    xy, ij, frac = [], [], []
    half = GRID_UM / 2.0
    for iy in range(ny):
        cy = half + iy * GRID_UM
        y0 = max(0, int(math.floor((cy - half) / coarse_mpp)))
        y1 = min(height, int(math.ceil((cy + half) / coarse_mpp)))
        if y1 <= y0:
            continue
        for ix in range(nx):
            cx = half + ix * GRID_UM
            x0 = max(0, int(math.floor((cx - half) / coarse_mpp)))
            x1 = min(width, int(math.ceil((cx + half) / coarse_mpp)))
            if x1 <= x0:
                continue
            value = rect_sum(integral, x0, y0, x1, y1) / float((x1 - x0) * (y1 - y0))
            if value >= MIN_COARSE_TISSUE:
                xy.append((cx, cy))
                ij.append((iy, ix))
                frac.append(value)
    return np.asarray(xy, np.float32), np.asarray(ij, np.int32), np.asarray(frac, np.float32)


def read_center_patch(reader, cx: float, cy: float) -> np.ndarray:
    fov = VIEW_MPP * IMAGE_SIZE
    requested = (cx - fov / 2, cy - fov / 2, cx + fov / 2, cy + fov / 2)
    clipped = (
        max(0.0, requested[0]),
        max(0.0, requested[1]),
        min(reader.geometry.width_um, requested[2]),
        min(reader.geometry.height_um, requested[3]),
    )
    canvas = np.full((IMAGE_SIZE, IMAGE_SIZE, 3), 255, dtype=np.uint8)
    if clipped[2] <= clipped[0] or clipped[3] <= clipped[1]:
        return canvas
    rgb = reader.read_bbox_um(clipped, target_mpp=VIEW_MPP)
    x0 = int(round((clipped[0] - requested[0]) / VIEW_MPP))
    y0 = int(round((clipped[1] - requested[1]) / VIEW_MPP))
    x1 = min(IMAGE_SIZE, x0 + rgb.shape[1])
    y1 = min(IMAGE_SIZE, y0 + rgb.shape[0])
    if x1 > x0 and y1 > y0:
        canvas[y0:y1, x0:x1] = rgb[: y1 - y0, : x1 - x0]
    return canvas


def graph_edges(grid_ij: np.ndarray) -> np.ndarray:
    lookup = {tuple(value): index for index, value in enumerate(grid_ij.tolist())}
    edges = []
    for index, (iy, ix) in enumerate(grid_ij.tolist()):
        for key in ((iy + 1, ix), (iy, ix + 1)):
            other = lookup.get(key)
            if other is not None:
                edges.append((index, other))
    return np.asarray(edges, dtype=np.int64).T if edges else np.empty((2, 0), dtype=np.int64)


def cache_complete(path: Path) -> bool:
    if not path.exists():
        return False
    try:
        with h5py.File(path, "r") as handle:
            return bool(handle.attrs.get("complete", 0)) and handle["uni/embedding"].shape[1:] == (1, 1024)
    except Exception:
        return False


def prepare_slide(row, output_dir: Path, overwrite: bool, cucim_workers: int) -> PreparedSlide | dict:
    """CPU/I/O stage; safe to run concurrently because each reader owns one WSI."""
    output = output_dir / f"{row.sample_id}.patch_cache.h5"
    temporary = output.with_suffix(output.suffix + ".incomplete")
    if cache_complete(output) and not overwrite:
        with h5py.File(output, "r") as handle:
            n = int(handle["regions/xy_um"].shape[0])
        return {"cohort": row.cohort, "case_id": row.case_id, "sample_id": row.sample_id, "status": "skipped_complete", "n_regions": n, "seconds": 0.0, "output_h5": str(output)}
    started = time.time()
    if temporary.exists():
        temporary.unlink()
    with open_extraction_slide(Path(row.slide_path)) as reader:
        xy_all, ij_all, coarse_all = candidate_grid(reader)
        if len(xy_all) == 0:
            raise RuntimeError("no_coarse_tissue_candidates")
        rng = np.random.default_rng(stable_seed(str(row.sample_id)))
        order = rng.permutation(len(xy_all))[: min(len(xy_all), MAX_CANDIDATES_TO_DECODE)]
        selected_xy, selected_ij, selected_coarse, selected_rgb, images = [], [], [], [], []
        decoded_count = 0
        if hasattr(reader, "read_center_patches"):
            # Decode a modest candidate block using the reader's optimized path
            # (cuCIM worker threads or ordered IBL SQLite BLOB prefetch),
            # preserve the seeded candidate order, and stop once 192 tissue
            # patches have been accepted. The spatial/statistical contract is
            # identical to the scalar reader.
            for block_start in range(0, len(order), 256):
                block_index = order[block_start : block_start + 256]
                block = reader.read_center_patches(
                    xy_all[block_index], VIEW_MPP, IMAGE_SIZE, cucim_workers, batch_size=32
                )
                array = block.astype(np.float32)
                brightness = array.mean(axis=3)
                chroma = array.max(axis=3) - array.min(axis=3)
                fractions = ((brightness < 232.0) & (chroma > 8.0)).mean(axis=(1, 2))
                for source_index, image, rgb_fraction in zip(block_index, block, fractions):
                    decoded_count += 1
                    if float(rgb_fraction) < MIN_RGB_TISSUE:
                        continue
                    selected_xy.append(xy_all[source_index])
                    selected_ij.append(ij_all[source_index])
                    selected_coarse.append(coarse_all[source_index])
                    selected_rgb.append(float(rgb_fraction))
                    images.append(image)
                    if len(images) >= MAX_REGIONS:
                        break
                if len(images) >= MAX_REGIONS:
                    break
        else:
            for source_index in order:
                decoded_count += 1
                image = read_center_patch(reader, *xy_all[source_index].tolist())
                rgb_fraction = tissue_fraction(image)
                if rgb_fraction < MIN_RGB_TISSUE:
                    continue
                selected_xy.append(xy_all[source_index])
                selected_ij.append(ij_all[source_index])
                selected_coarse.append(coarse_all[source_index])
                selected_rgb.append(rgb_fraction)
                images.append(image)
                if len(images) >= MAX_REGIONS:
                    break
        geometry = reader.geometry.to_dict()
    if not images:
        raise RuntimeError("no_rgb_tissue_regions")
    array = np.stack(images).astype(np.float32) / 255.0
    array = (array - IMAGENET_MEAN) / IMAGENET_STD
    tensor = torch.from_numpy(np.ascontiguousarray(array.transpose(0, 3, 1, 2)))
    return PreparedSlide(
        row=row,
        output=output,
        temporary=temporary,
        tensor=tensor,
        xy=np.asarray(selected_xy, np.float32),
        grid_ij=np.asarray(selected_ij, np.int32),
        coarse_tissue=np.asarray(selected_coarse, np.float16),
        rgb_tissue=np.asarray(selected_rgb, np.float16),
        geometry=geometry,
        n_coarse_candidates=int(len(xy_all)),
        n_candidates_decoded=int(decoded_count),
        started=started,
        io_seconds=time.time() - started,
    )


@torch.inference_mode()
def infer_and_write(prepared: PreparedSlide, model, device: torch.device, batch_size: int, checkpoint_hash: str) -> dict:
    """Single-GPU stage; preparation of later slides continues in worker threads."""
    gpu_started = time.time()
    embedding_parts, probability_parts = [], []
    source = prepared.tensor
    for start in range(0, len(source), batch_size):
        tensor = source[start : start + batch_size]
        if device.type == "cuda":
            tensor = tensor.pin_memory().to(device, non_blocking=True)
        else:
            tensor = tensor.to(device)
        with torch.autocast(device_type="cuda", dtype=torch.float16, enabled=device.type == "cuda"):
            features = model.backbone(tensor)
            logits = model.classifier(features)
        embedding_parts.append(features.float().cpu().numpy().astype(np.float16))
        probability_parts.append(torch.softmax(logits.float(), dim=1).cpu().numpy().astype(np.float16))
    embedding = np.concatenate(embedding_parts)
    probability = np.concatenate(probability_parts)
    edge_index = graph_edges(prepared.grid_ij)
    center = embedding.astype(np.float32)
    sums = center.copy()
    degree = np.ones(len(center), dtype=np.float32)
    if edge_index.shape[1]:
        u, v = edge_index
        np.add.at(sums, u, center[v])
        np.add.at(sums, v, center[u])
        np.add.at(degree, u, 1.0)
        np.add.at(degree, v, 1.0)
    delta = (sums / degree[:, None] - center).astype(np.float16)
    prepared.output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(prepared.temporary, "w") as handle:
        handle.attrs["schema"] = "patch_aligned_topology_survival_v2_external_uni_cache_v1"
        handle.attrs["complete"] = 0
        handle.attrs["cohort"] = prepared.row.cohort
        handle.attrs["case_id"] = prepared.row.case_id
        handle.attrs["sample_id"] = prepared.row.sample_id
        handle.attrs["source_slide"] = prepared.row.slide_path
        handle.attrs["selection"] = "full_tissue_250um_grid_outcome_blind_deterministic_sample"
        handle.attrs["external_domain_difference"] = "RGB_tissue_QC_only; no CellViT n_cells>=20 filter"
        handle.attrs["uni_checkpoint_sha256"] = checkpoint_hash
        handle.attrs["view_mpp"] = VIEW_MPP
        handle.attrs["image_size_px"] = IMAGE_SIZE
        handle.attrs["grid_um"] = GRID_UM
        handle.attrs["slide_geometry_json"] = json.dumps(prepared.geometry)
        handle.attrs["loader_mode"] = "threaded_prefetch_single_gpu"
        handle.create_dataset("regions/xy_um", data=prepared.xy, compression="gzip", compression_opts=4)
        handle.create_dataset("regions/grid_ij", data=prepared.grid_ij, compression="gzip", compression_opts=4)
        handle.create_dataset("regions/coarse_tissue_fraction", data=prepared.coarse_tissue, compression="gzip", compression_opts=4)
        handle.create_dataset("regions/rgb_tissue_fraction", data=prepared.rgb_tissue, compression="gzip", compression_opts=4)
        handle.create_dataset("graph/edge_index", data=edge_index, compression="gzip", compression_opts=4)
        handle.create_dataset("uni/embedding", data=embedding[:, None, :], compression="gzip", compression_opts=4, chunks=(min(64, len(embedding)), 1, 1024))
        handle.create_dataset("uni/neighbor_delta_center", data=delta, compression="gzip", compression_opts=4, chunks=(min(128, len(delta)), 1024))
        handle.create_dataset("uni/coarse_probability", data=probability[:, None, :], compression="gzip", compression_opts=4)
        handle.attrs["complete"] = 1
        handle.flush()
    os.replace(prepared.temporary, prepared.output)
    elapsed = time.time() - prepared.started
    gpu_seconds = time.time() - gpu_started
    return {
        "cohort": prepared.row.cohort,
        "case_id": prepared.row.case_id,
        "sample_id": prepared.row.sample_id,
        "status": "complete",
        "n_coarse_candidates": prepared.n_coarse_candidates,
        "n_candidates_decoded": prepared.n_candidates_decoded,
        "n_regions": int(len(prepared.xy)),
        "mean_rgb_tissue_fraction": float(prepared.rgb_tissue.astype(np.float32).mean()),
        "n_graph_edges": int(edge_index.shape[1]),
        "io_seconds": prepared.io_seconds,
        "gpu_write_seconds": gpu_seconds,
        "seconds": elapsed,
        "regions_per_second": len(prepared.xy) / max(elapsed, 1e-8),
        "output_h5": str(prepared.output),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--batch-size", type=int, default=160)
    parser.add_argument("--loader-workers", type=int, default=6)
    parser.add_argument("--prefetch-slides", type=int, default=12)
    parser.add_argument("--cucim-workers", type=int, default=8)
    parser.add_argument("--limit-slides", type=int, default=0)
    parser.add_argument(
        "--sample-ids",
        nargs="*",
        default=(),
        help="Optional exact sample IDs for targeted repair/resume runs.",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    manifest = pd.read_csv(args.manifest, sep="\t")
    required = {"cohort", "case_id", "sample_id", "slide_path", "eligible"}
    missing_columns = sorted(required.difference(manifest.columns))
    if missing_columns:
        raise ValueError(f"Manifest is missing required columns: {missing_columns}")
    manifest = manifest[manifest.eligible.astype(str).str.lower().isin({"true", "1"})].copy()
    if args.sample_ids:
        requested = {str(value) for value in args.sample_ids}
        manifest = manifest[manifest.sample_id.astype(str).isin(requested)].copy()
        missing = requested.difference(set(manifest.sample_id.astype(str)))
        if missing:
            raise ValueError(f"Requested sample IDs are absent or ineligible: {sorted(missing)}")
    manifest = manifest.sort_values(["case_id", "sample_id"]).reset_index(drop=True)
    if args.limit_slides > 0:
        manifest = manifest.head(args.limit_slides).copy()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    audit_path = output_dir / "patch_cache_audit.tsv"
    existing = pd.read_csv(audit_path, sep="\t") if audit_path.exists() else pd.DataFrame()
    records = existing.to_dict("records") if not existing.empty else []
    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    model = load_uni_classifier(UNI_CHECKPOINT, device)
    checkpoint_hash = sha256_file(UNI_CHECKPOINT)
    rows = list(manifest.itertuples(index=False))
    next_row = 0
    completed_now = 0
    pending = {}
    with ThreadPoolExecutor(max_workers=max(1, args.loader_workers), thread_name_prefix="wsi-loader") as pool:
        while next_row < len(rows) and len(pending) < max(1, args.prefetch_slides):
            row = rows[next_row]
            pending[pool.submit(prepare_slide, row, output_dir, args.overwrite, args.cucim_workers)] = row
            next_row += 1
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                row = pending.pop(future)
                try:
                    prepared = future.result()
                    record = prepared if isinstance(prepared, dict) else infer_and_write(prepared, model, device, args.batch_size, checkpoint_hash)
                except Exception as exc:
                    record = {"cohort": row.cohort, "case_id": row.case_id, "sample_id": row.sample_id, "status": "error", "error": repr(exc), "n_regions": 0, "seconds": 0.0, "output_h5": ""}
                records = [value for value in records if str(value.get("sample_id")) != str(row.sample_id)]
                records.append(record)
                completed_now += 1
                atomic_tsv(audit_path, pd.DataFrame(records).sort_values(["case_id", "sample_id"]))
                print(json.dumps({"completed_this_run": completed_now, "total": len(manifest), "loader_workers": args.loader_workers, "prefetch_slides": args.prefetch_slides, "cucim_workers": args.cucim_workers, "batch_size": args.batch_size, **record}, ensure_ascii=False), flush=True)
                if next_row < len(rows):
                    new_row = rows[next_row]
                    pending[pool.submit(prepare_slide, new_row, output_dir, args.overwrite, args.cucim_workers)] = new_row
                    next_row += 1


if __name__ == "__main__":
    main()
