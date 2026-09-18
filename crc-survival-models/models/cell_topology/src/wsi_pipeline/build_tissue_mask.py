#!/usr/bin/env python3
"""Build the mandatory whole-slide tissue domain before any UNI analysis."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from common import atomic_json, sha256_file
from config import FROZEN, OUTPUT_ROOT
from slide_io import open_physical_slide
from tissue_mask_utils import segment_tissue


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", required=True, choices=["tcga", "sr386", "custom"])
    parser.add_argument("--slide", type=Path, required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--case-id", default="")
    parser.add_argument("--mpp-x", type=float)
    parser.add_argument("--mpp-y", type=float)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT / "tissue_masks")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def write_preview(rgb: np.ndarray, mask: np.ndarray, path: Path) -> None:
    boundary = mask ^ ndi.binary_erosion(mask)
    overlay = rgb.astype(np.float32).copy()
    overlay[mask] = 0.72 * overlay[mask] + 0.28 * np.asarray([44, 146, 160], dtype=np.float32)
    overlay[boundary] = np.asarray([205, 58, 64], dtype=np.float32)
    image = Image.fromarray(np.clip(overlay, 0, 255).astype(np.uint8), mode="RGB")
    maximum = 1600
    if max(image.size) > maximum:
        scale = maximum / max(image.size)
        image = image.resize(
            (max(1, int(round(image.width * scale))), max(1, int(round(image.height * scale)))),
            Image.Resampling.LANCZOS,
        )
    temporary = path.with_suffix(path.suffix + ".incomplete")
    image.save(temporary, format="PNG", optimize=True)
    os.replace(temporary, path)


def main() -> None:
    args = parse_args()
    if not args.slide.exists():
        raise FileNotFoundError(args.slide)
    output_dir = args.output_root / args.cohort / args.sample_id
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "tissue_mask_manifest.json"
    if manifest_path.exists() and not args.overwrite:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        if payload.get("tissue_first") is True and Path(payload["mask_npz"]).exists():
            print(f"Already complete: {manifest_path}")
            return

    started = time.time()
    with open_physical_slide(args.slide, args.mpp_x, args.mpp_y) as reader:
        rgb = reader.read_bbox_um(
            (0.0, 0.0, reader.geometry.width_um, reader.geometry.height_um),
            target_mpp=FROZEN.tissue_mask_mpp,
        )
        mask, metrics = segment_tissue(rgb, FROZEN.tissue_mask_mpp)
        if not mask.any():
            raise RuntimeError("Whole-slide tissue segmentation produced an empty tissue mask")

        mask_path = output_dir / "tissue_mask_8um.npz"
        temporary_mask = mask_path.with_suffix(mask_path.suffix + ".incomplete")
        with temporary_mask.open("wb") as handle:
            np.savez_compressed(
                handle,
                mask=mask.astype(np.uint8),
                mpp=np.asarray([FROZEN.tissue_mask_mpp], dtype=np.float64),
            )
        os.replace(temporary_mask, mask_path)
        preview_path = output_dir / "tissue_mask_qc_preview.png"
        write_preview(rgb, mask, preview_path)
        slide_stat = args.slide.stat()
        payload = {
            "schema": "whole_slide_tissue_mask_v3_tissue_first",
            "status": "complete",
            "tissue_first": True,
            "sample_id": args.sample_id,
            "case_id": args.case_id or args.sample_id,
            "cohort": args.cohort,
            "source_slide": str(args.slide),
            "source_slide_size_bytes": int(slide_stat.st_size),
            "source_slide_mtime_ns": int(slide_stat.st_mtime_ns),
            "slide_geometry": reader.geometry.to_dict(),
            "mask_npz": str(mask_path),
            "mask_sha256": sha256_file(mask_path),
            "mask_mpp": FROZEN.tissue_mask_mpp,
            "mask_shape": list(mask.shape),
            "mask_origin_um": [0.0, 0.0],
            "qc_preview": str(preview_path),
            "segmentation_definition": {
                "colour_space": "adaptive grayscale plus HSV saturation",
                "hard_domain_before_uni": True,
                "config": FROZEN.to_dict(),
            },
            "metrics": metrics,
            "elapsed_seconds": time.time() - started,
        }
        atomic_json(manifest_path, payload)
        print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
