"""Frozen mathematical definitions used by CoVarST; historical CLI omitted."""
from __future__ import annotations

import argparse

import csv

import gc

import json

import sys

import time

from pathlib import Path

import h5py

import numpy as np

import pandas as pd

import torch

from PIL import Image, ImageDraw

from tqdm import tqdm

def load_paramnet(paramnet_root: Path, device: torch.device) -> torch.nn.Module:
    sys.path.insert(0, str(paramnet_root))
    from source.model import ParamNet

    ckpt_path = paramnet_root / "checkpoints" / "ParamNet-Uni.pt"
    if not ckpt_path.exists():
        raise FileNotFoundError(f"ParamNet-Uni checkpoint not found: {ckpt_path}")
    model = ParamNet().to(device)
    ckpt = torch.load(ckpt_path, map_location=device, weights_only=False)
    state = ckpt["net_G_A"] if isinstance(ckpt, dict) and "net_G_A" in ckpt else ckpt
    model.load_state_dict(state)
    model.eval()
    for param in model.parameters():
        param.requires_grad_(False)
    return model

def normalize_batch(paramnet: torch.nn.Module, batch: np.ndarray, device: torch.device) -> np.ndarray:
    arr = np.asarray(batch, dtype=np.float32)
    x = ((arr / 255.0) - 0.5) / 0.5
    x_t = torch.from_numpy(x.transpose(0, 3, 1, 2)).to(device=device, dtype=torch.float32)
    with torch.no_grad():
        y = paramnet(x_t).detach().cpu().numpy()
    y = ((y * 0.5 + 0.5) * 255.0).clip(0, 255).astype(np.uint8)
    return y.transpose(0, 2, 3, 1)
