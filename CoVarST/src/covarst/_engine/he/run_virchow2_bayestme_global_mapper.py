"""Frozen mathematical definitions used by CoVarST; historical CLI omitted."""
from __future__ import annotations

import argparse

import copy

import json

import math

import time

from pathlib import Path

import numpy as np

import pandas as pd

import torch

from scipy.spatial import cKDTree

from openst_final.he_uni_mamba_bnn_inr import (
    normalize_slide_coordinates,
)

from run_virchow2_fseg_k64_mapper import (
    BoundedINR,
    SlideContextK64,
    SlideData,
    _evaluate_split,
    _loss,
    _to_device_graph,
)

EPS = 1.0e-8

def _morphology_graph(
    coordinates: np.ndarray,
    features: np.ndarray,
    *,
    neighbors: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    coordinates = np.asarray(coordinates, dtype=np.float32)
    features = np.asarray(features, dtype=np.float32)
    feature_norm = features / np.maximum(
        np.linalg.norm(features, axis=1, keepdims=True),
        EPS,
    )
    k = min(int(neighbors) + 1, len(coordinates))
    distance, index = cKDTree(coordinates).query(coordinates, k=k)
    if k == 1:
        distance = distance[:, None]
        index = index[:, None]
    source = np.repeat(np.arange(len(coordinates), dtype=np.int64), k)
    target = index.reshape(-1).astype(np.int64)
    flat_distance = distance.reshape(-1).astype(np.float32)
    positive = flat_distance[flat_distance > 0]
    scale = float(np.median(positive)) if len(positive) else 1.0
    cosine = np.sum(
        feature_norm[source] * feature_norm[target],
        axis=1,
    )
    similarity = np.clip(0.5 + 0.5 * cosine, 0.0, 1.0)
    spatial = np.exp(
        -flat_distance / max(scale, EPS)
    ).astype(np.float32)
    weight = np.maximum(
        spatial * np.square(similarity).astype(np.float32),
        1.0e-4,
    )
    return source, target, weight
