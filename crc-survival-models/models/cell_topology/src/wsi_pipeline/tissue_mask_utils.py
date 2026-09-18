#!/usr/bin/env python3
"""Whole-slide tissue mask construction and coordinate projection utilities."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np
from scipy import ndimage as ndi
from skimage.color import rgb2gray, rgb2hsv
from skimage.filters import threshold_otsu
from skimage.morphology import binary_closing, binary_opening, disk, remove_small_holes, remove_small_objects

from common import sha256_file
from config import FROZEN


def segment_tissue(rgb: np.ndarray, mpp: float) -> tuple[np.ndarray, dict]:
    """Return a cleaned H&E tissue mask at low resolution.

    The adaptive grayscale term captures pale eosin tissue, while the
    saturation term preserves lightly stained mucosa that an intensity-only
    rule can miss. Tiny coloured annotations and scanner specks are removed by
    a physical-area component threshold.
    """

    if rgb.ndim != 3 or rgb.shape[2] != 3:
        raise ValueError(f"Expected RGB image, received shape={rgb.shape}")
    unit = np.clip(rgb.astype(np.float32) / 255.0, 0.0, 1.0)
    gray = rgb2gray(unit)
    hsv = rgb2hsv(unit)
    saturation = hsv[..., 1]
    value = hsv[..., 2]
    usable = gray[np.isfinite(gray) & (gray < 0.995)]
    if usable.size >= 256 and float(np.ptp(usable)) > 1e-4:
        adaptive_gray = float(np.clip(threshold_otsu(usable), 0.68, 0.94))
    else:
        adaptive_gray = 0.90
    non_black = value > 0.06
    raw = non_black & (
        (gray < adaptive_gray)
        | ((saturation > 0.045) & (value < 0.985))
    )

    close_radius_px = max(1, int(round(FROZEN.tissue_closing_radius_um / mpp)))
    minimum_pixels = max(
        1,
        int(math.ceil(FROZEN.tissue_component_min_mm2 * 1e6 / (mpp * mpp))),
    )
    maximum_hole_pixels = max(
        1,
        int(math.ceil(FROZEN.tissue_small_hole_max_mm2 * 1e6 / (mpp * mpp))),
    )
    cleaned = binary_closing(raw, footprint=disk(close_radius_px))
    cleaned = binary_opening(cleaned, footprint=disk(1))
    cleaned = remove_small_holes(cleaned, area_threshold=maximum_hole_pixels)
    cleaned = remove_small_objects(cleaned, min_size=minimum_pixels)
    cleaned = np.asarray(cleaned, dtype=bool)

    labels, n_components = ndi.label(cleaned)
    component_pixels = np.bincount(labels.ravel())[1:]
    component_areas_mm2 = np.sort(component_pixels.astype(float) * mpp * mpp / 1e6)[::-1]
    metrics = {
        "adaptive_gray_threshold": adaptive_gray,
        "raw_tissue_fraction": float(raw.mean()),
        "clean_tissue_fraction": float(cleaned.mean()),
        "tissue_area_mm2": float(cleaned.sum() * mpp * mpp / 1e6),
        "n_tissue_components": int(n_components),
        "largest_component_area_mm2": float(component_areas_mm2[0]) if component_areas_mm2.size else 0.0,
        "component_areas_mm2": component_areas_mm2.tolist(),
        "minimum_component_area_mm2": FROZEN.tissue_component_min_mm2,
        "maximum_filled_hole_area_mm2": FROZEN.tissue_small_hole_max_mm2,
        "closing_radius_um": FROZEN.tissue_closing_radius_um,
    }
    return cleaned, metrics


def load_tissue_mask(manifest_path: Path) -> tuple[np.ndarray, dict]:
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    if manifest.get("status") != "complete" or manifest.get("tissue_first") is not True:
        raise RuntimeError(f"Not a complete tissue-first manifest: {manifest_path}")
    mask_path = Path(manifest["mask_npz"])
    if sha256_file(mask_path) != manifest["mask_sha256"]:
        raise RuntimeError(f"Tissue mask SHA256 mismatch: {mask_path}")
    with np.load(mask_path) as bundle:
        mask = bundle["mask"].astype(bool)
        stored_mpp = float(bundle["mpp"][0])
    if not np.isclose(stored_mpp, float(manifest["mask_mpp"])):
        raise RuntimeError("Tissue mask MPP differs between NPZ and manifest")
    expected_shape = tuple(int(value) for value in manifest["mask_shape"])
    if mask.shape != expected_shape:
        raise RuntimeError(f"Tissue mask shape mismatch: {mask.shape} != {expected_shape}")
    return mask, manifest


def bbox_fraction(mask: np.ndarray, mpp: float, bbox_um: Iterable[float]) -> float:
    x0, y0, x1, y1 = map(float, bbox_um)
    col0 = max(0, int(math.floor(x0 / mpp)))
    row0 = max(0, int(math.floor(y0 / mpp)))
    col1 = min(mask.shape[1], int(math.ceil(x1 / mpp)))
    row1 = min(mask.shape[0], int(math.ceil(y1 / mpp)))
    if col1 <= col0 or row1 <= row0:
        return 0.0
    return float(mask[row0:row1, col0:col1].mean())


def sample_points(mask: np.ndarray, mpp: float, points_xy_um: np.ndarray) -> np.ndarray:
    points = np.asarray(points_xy_um, dtype=float)
    col = np.floor(points[:, 0] / mpp).astype(np.int64)
    row = np.floor(points[:, 1] / mpp).astype(np.int64)
    valid = (row >= 0) & (row < mask.shape[0]) & (col >= 0) & (col < mask.shape[1])
    result = np.zeros(len(points), dtype=bool)
    result[valid] = mask[row[valid], col[valid]]
    return result


def rasterize_bbox(
    mask: np.ndarray,
    mask_mpp: float,
    bbox_um: Tuple[float, float, float, float],
    shape: Tuple[int, int],
    raster_mpp: float,
) -> np.ndarray:
    x0, y0, _, _ = bbox_um
    xs = x0 + (np.arange(shape[1], dtype=float) + 0.5) * raster_mpp
    ys = y0 + (np.arange(shape[0], dtype=float) + 0.5) * raster_mpp
    cols = np.floor(xs / mask_mpp).astype(np.int64)
    rows = np.floor(ys / mask_mpp).astype(np.int64)
    valid_cols = (cols >= 0) & (cols < mask.shape[1])
    valid_rows = (rows >= 0) & (rows < mask.shape[0])
    result = np.zeros(shape, dtype=bool)
    if valid_rows.any() and valid_cols.any():
        result[np.ix_(valid_rows, valid_cols)] = mask[np.ix_(rows[valid_rows], cols[valid_cols])]
    return result
