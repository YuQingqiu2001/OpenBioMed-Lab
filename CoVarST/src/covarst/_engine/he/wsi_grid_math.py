"""Frozen mathematical definitions used by CoVarST; historical CLI omitted."""
from __future__ import annotations

import argparse

import json

import math

import re

import time

from concurrent.futures import ProcessPoolExecutor, as_completed

from datetime import datetime, timezone

from pathlib import Path

from typing import Any

import cv2

import numpy as np

import openslide

import pandas as pd

from PIL import Image, ImageDraw

from scipy.ndimage import binary_closing, binary_fill_holes, binary_opening, label

def property_float(properties: Any, names: tuple[str, ...]) -> float | None:
    for name in names:
        value = properties.get(name)
        if value is None:
            continue
        try:
            parsed = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(parsed) and parsed > 0:
            return parsed
    return None

def slide_mpp(slide: openslide.OpenSlide) -> tuple[float, float, str]:
    properties = slide.properties
    mpp_x = property_float(properties, (openslide.PROPERTY_NAME_MPP_X, "aperio.MPP"))
    mpp_y = property_float(properties, (openslide.PROPERTY_NAME_MPP_Y, "aperio.MPP"))
    source = "openslide_mpp"
    if mpp_x is None or mpp_y is None:
        objective = property_float(properties, (openslide.PROPERTY_NAME_OBJECTIVE_POWER, "aperio.AppMag"))
        if objective is None:
            raise ValueError("missing MPP and objective-power metadata")
        inferred = 10.0 / objective
        mpp_x = inferred if mpp_x is None else mpp_x
        mpp_y = inferred if mpp_y is None else mpp_y
        source = "objective_power_fallback"
    if not (0.10 <= mpp_x <= 4.00 and 0.10 <= mpp_y <= 4.00):
        raise ValueError(f"implausible MPP ({mpp_x}, {mpp_y})")
    anisotropy = abs(mpp_x - mpp_y) / max((mpp_x + mpp_y) / 2.0, 1.0e-8)
    if anisotropy > 0.05:
        raise ValueError(f"anisotropic MPP ({mpp_x}, {mpp_y})")
    return float(mpp_x), float(mpp_y), source

def segment_tissue(thumbnail: np.ndarray) -> tuple[np.ndarray, dict[str, float | int]]:
    rgb = np.asarray(thumbnail, dtype=np.uint8)
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    saturation = hsv[..., 1].astype(np.float32) / 255.0
    value = hsv[..., 2].astype(np.float32) / 255.0
    blurred = cv2.GaussianBlur(gray, (5, 5), 0)
    otsu_threshold, _ = cv2.threshold(
        blurred, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU
    )
    chromatic_tissue = (saturation >= 0.045) & (value <= 0.985)
    dark_tissue = (gray <= min(float(otsu_threshold) + 18.0, 225.0)) & (value >= 0.10)
    initial = (chromatic_tissue | dark_tissue) & (value >= 0.08)
    mask = binary_opening(initial, structure=np.ones((3, 3), dtype=bool))
    mask = binary_closing(mask, structure=np.ones((7, 7), dtype=bool))
    mask = binary_fill_holes(mask)
    components, count = label(mask)
    sizes = np.bincount(components.ravel(), minlength=count + 1)
    minimum_component = max(64, int(mask.size * 2.0e-5))
    retained = np.flatnonzero(sizes >= minimum_component)
    retained = retained[retained != 0]
    mask = np.isin(components, retained)
    return mask.astype(bool), {
        "otsu_gray_threshold": float(otsu_threshold),
        "minimum_component_pixels": int(minimum_component),
        "component_count_initial": int(count),
        "component_count_retained": int(len(retained)),
        "tissue_fraction_thumbnail": float(mask.mean()),
    }

def mask_fraction_for_boxes(
    mask: np.ndarray,
    left_x: np.ndarray,
    left_y: np.ndarray,
    patch_width: int,
    patch_height: int,
    scale_x: float,
    scale_y: float,
) -> tuple[np.ndarray, np.ndarray]:
    height, width = mask.shape
    x0 = np.clip(np.floor(left_x * scale_x).astype(int), 0, width - 1)
    y0 = np.clip(np.floor(left_y * scale_y).astype(int), 0, height - 1)
    x1 = np.clip(np.ceil((left_x + patch_width) * scale_x).astype(int), x0 + 1, width)
    y1 = np.clip(np.ceil((left_y + patch_height) * scale_y).astype(int), y0 + 1, height)
    integral = cv2.integral(mask.astype(np.uint8), sdepth=cv2.CV_64F)
    tissue = integral[y1, x1] - integral[y0, x1] - integral[y1, x0] + integral[y0, x0]
    area = np.maximum((x1 - x0) * (y1 - y0), 1)
    fraction = tissue / area
    centre_x = np.clip(np.rint((left_x + patch_width / 2.0) * scale_x).astype(int), 0, width - 1)
    centre_y = np.clip(np.rint((left_y + patch_height / 2.0) * scale_y).astype(int), 0, height - 1)
    return fraction.astype(np.float32), mask[centre_y, centre_x]

def build_patch_grid(
    width: int,
    height: int,
    mpp_x: float,
    mpp_y: float,
    mask: np.ndarray,
    scale_x: float,
    scale_y: float,
    patch_field_um: float,
    stride_um: float,
    minimum_tissue_fraction: float,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    patch_width = max(1, int(round(patch_field_um / mpp_x)))
    patch_height = max(1, int(round(patch_field_um / mpp_y)))
    horizontal_pitch_px = float(stride_um / mpp_x)
    vertical_pitch_px = float((math.sqrt(3.0) / 2.0) * stride_um / mpp_y)
    centre_y_values = np.arange(
        patch_height / 2.0,
        max(height - patch_height / 2.0 + 1.0e-6, patch_height / 2.0 + 1.0e-6),
        vertical_pitch_px,
        dtype=np.float64,
    )
    left_x_blocks: list[np.ndarray] = []
    left_y_blocks: list[np.ndarray] = []
    centre_x_blocks: list[np.ndarray] = []
    centre_y_blocks: list[np.ndarray] = []
    row_blocks: list[np.ndarray] = []
    column_blocks: list[np.ndarray] = []
    for grid_row, centre_y in enumerate(centre_y_values):
        offset = (horizontal_pitch_px / 2.0) if grid_row % 2 else 0.0
        centre_x_values = np.arange(
            patch_width / 2.0 + offset,
            max(width - patch_width / 2.0 + 1.0e-6, patch_width / 2.0 + offset + 1.0e-6),
            horizontal_pitch_px,
            dtype=np.float64,
        )
        centre_x_values = centre_x_values[
            centre_x_values <= width - patch_width / 2.0 + 1.0e-6
        ]
        if not len(centre_x_values):
            continue
        centre_y_block = np.full(len(centre_x_values), centre_y, dtype=np.float64)
        left_x_blocks.append(
            np.rint(centre_x_values - patch_width / 2.0).astype(np.int64)
        )
        left_y_blocks.append(
            np.rint(centre_y_block - patch_height / 2.0).astype(np.int64)
        )
        centre_x_blocks.append(centre_x_values)
        centre_y_blocks.append(centre_y_block)
        row_blocks.append(np.full(len(centre_x_values), grid_row, dtype=np.int64))
        column_blocks.append(np.arange(len(centre_x_values), dtype=np.int64))
    if not left_x_blocks:
        raise RuntimeError("empty physical hexagonal patch grid")
    left_x = np.concatenate(left_x_blocks)
    left_y = np.concatenate(left_y_blocks)
    centre_x = np.concatenate(centre_x_blocks)
    centre_y = np.concatenate(centre_y_blocks)
    grid_rows = np.concatenate(row_blocks)
    grid_columns = np.concatenate(column_blocks)
    if not len(left_x):
        raise RuntimeError("empty physical patch grid")
    fraction, centre_tissue = mask_fraction_for_boxes(
        mask,
        left_x,
        left_y,
        patch_width,
        patch_height,
        scale_x,
        scale_y,
    )
    keep = centre_tissue & (fraction >= minimum_tissue_fraction)
    left_x = left_x[keep]
    left_y = left_y[keep]
    centre_x = centre_x[keep]
    centre_y = centre_y[keep]
    grid_rows = grid_rows[keep]
    grid_columns = grid_columns[keep]
    fraction = fraction[keep]
    if not len(left_x):
        raise RuntimeError("tissue-filtered HEST patch grid is empty")
    order = np.lexsort((left_x, left_y))
    left_x = left_x[order]
    left_y = left_y[order]
    centre_x = centre_x[order]
    centre_y = centre_y[order]
    grid_rows = grid_rows[order]
    grid_columns = grid_columns[order]
    fraction = fraction[order]
    frame = pd.DataFrame(
        {
            "patch_index": np.arange(len(left_x), dtype=np.int64),
            "patch_left_x": left_x,
            "patch_left_y": left_y,
            "patch_center_x": centre_x,
            "patch_center_y": centre_y,
            "patch_width_src": patch_width,
            "patch_height_src": patch_height,
            "grid_row": grid_rows,
            "grid_column": grid_columns,
            "spot_hex_radius_um": float(stride_um / math.sqrt(3.0)),
            "tissue_fraction": fraction,
        }
    )
    return frame, {
        "patch_width_src_px": int(patch_width),
        "patch_height_src_px": int(patch_height),
        "lattice": "triangular_centres_regular_hexagonal_voronoi",
        "nearest_neighbour_spacing_um": float(stride_um),
        "horizontal_pitch_src_px": float(horizontal_pitch_px),
        "vertical_pitch_src_px": float(vertical_pitch_px),
        "alternate_row_offset_src_px": float(horizontal_pitch_px / 2.0),
        "hex_cell_circumradius_um": float(stride_um / math.sqrt(3.0)),
        "candidate_grid_rows": int(len(keep)),
        "retained_patch_rows": int(len(frame)),
        "retained_fraction": float(keep.mean()),
    }
