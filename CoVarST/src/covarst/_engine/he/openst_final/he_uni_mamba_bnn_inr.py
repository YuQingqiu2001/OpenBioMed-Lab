"""Frozen mathematical definitions used by CoVarST; historical CLI omitted."""
from __future__ import annotations

import hashlib

import json

import math

import os

import random

import shutil

import tempfile

import time

from dataclasses import asdict, dataclass, field, replace

from pathlib import Path

from typing import Any, Iterable, Mapping, Sequence

import numpy as np

import pandas as pd

import torch

from scipy import sparse as scipy_sparse

from scipy.spatial import cKDTree

from scipy.stats import spearmanr

from torch import Tensor, nn

from torch.nn import functional as F

from torch.utils.checkpoint import checkpoint

def normalize_slide_coordinates(coordinates: np.ndarray) -> np.ndarray:
    coordinates = np.asarray(coordinates, dtype=np.float32)
    lower = np.min(coordinates, axis=0)
    upper = np.max(coordinates, axis=0)
    center = 0.5 * (lower + upper)
    scale = np.maximum(0.5 * (upper - lower), 1.0)
    return np.clip((coordinates - center) / scale, -1.0, 1.0).astype(np.float32)

class SineLayer(nn.Module):
    def __init__(self, input_dim: int, output_dim: int, omega0: float):
        super().__init__()
        self.linear = nn.Linear(input_dim, output_dim)
        self.omega0 = float(omega0)
        bound = 1.0 / max(1, input_dim)
        nn.init.uniform_(self.linear.weight, -bound, bound)
        nn.init.zeros_(self.linear.bias)

    def forward(self, inputs: Tensor) -> Tensor:
        return torch.sin(self.omega0 * self.linear(inputs))
