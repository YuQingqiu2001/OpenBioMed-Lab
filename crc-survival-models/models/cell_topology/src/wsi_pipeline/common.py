"""Shared deterministic and atomic I/O utilities."""

from __future__ import annotations

import hashlib
import json
import os
import random
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass


def sha256_file(path: Path, chunk_bytes: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_bytes)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".incomplete")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    os.replace(temporary, path)


def decode_if_bytes(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not valid.any():
        return float("nan")
    return float(np.average(values[valid], weights=weights[valid]))


def weighted_upper_cvar(
    values: np.ndarray, weights: np.ndarray, upper_fraction: float = 0.20
) -> float:
    """Weighted mean of the upper tail without duplicating observations."""

    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not valid.any():
        return float("nan")
    values = values[valid]
    weights = weights[valid]
    order = np.argsort(values)[::-1]
    values = values[order]
    weights = weights[order]
    target = float(weights.sum() * upper_fraction)
    remaining = target
    total = 0.0
    used = 0.0
    for value, weight in zip(values, weights):
        take = min(float(weight), remaining)
        total += float(value) * take
        used += take
        remaining -= take
        if remaining <= 1e-12:
            break
    return float(total / used) if used > 0 else float("nan")


def geometric_mean_01(values: Sequence[float], eps: float = 1e-6) -> float:
    array = np.clip(np.asarray(values, dtype=float), 0.0, 1.0)
    if not np.isfinite(array).all():
        return float("nan")
    return float(np.exp(np.mean(np.log(np.maximum(array, eps)))))


def canonical_json_hash(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, ensure_ascii=False, default=str).encode(
        "utf-8"
    )
    return hashlib.sha256(encoded).hexdigest()


def require_paths(paths: Iterable[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required path(s):\n" + "\n".join(missing))


def dataframe_to_hdf5_group(group: Any, frame: Any) -> None:
    """Write a pandas DataFrame as compressed column datasets."""
    group.attrs["columns_json"] = json.dumps(list(frame.columns))
    group.attrs["n_rows"] = int(len(frame))
    for column in frame.columns:
        values = frame[column].to_numpy()
        if values.dtype.kind in {"O", "U"}:
            encoded_values = [
                ("" if value is None else str(value)).encode("utf-8") for value in values
            ]
            maximum = max((len(value) for value in encoded_values), default=1)
            dataset = group.create_dataset(
                column,
                data=np.asarray(encoded_values, dtype=f"S{maximum}"),
                compression="lzf",
            )
            dataset.attrs["text_encoding"] = "utf-8"
        else:
            group.create_dataset(column, data=values, compression="lzf")
