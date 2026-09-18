"""Portable inference for the pure cell-topology survival model.

The public checkpoint contains numeric arrays only.  Input files are the
topology-token HDF5 files created by the bundled feature pipeline.
"""

from __future__ import annotations

import hashlib
import json
from collections import deque
from pathlib import Path

import h5py
import numpy as np


FRONT = {
    "minimum_resolved_fraction": 0.50,
    "minimum_tumour_membership": 0.15,
    "minimum_tumour_purity_among_resolved": 0.35,
    "minimum_core_component_regions": 2,
    "front_expansion_hops": 1,
    "minimum_front_tokens_per_slide": 4,
    "minimum_front_tokens_per_patient": 16,
}
MAXIMUM_TOKENS_PER_SLIDE = 512
MAXIMUM_TOKENS_PER_PATIENT = 1536
SAMPLING_SEED = 20260907


def _decode(values: np.ndarray) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def _stable_seed(text: str, seed: int) -> int:
    digest = hashlib.sha256(f"{text}|{seed}".encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "little")


def _components(mask: np.ndarray, edges: np.ndarray) -> list[np.ndarray]:
    adjacency = [[] for _ in range(len(mask))]
    for left, right in edges.T:
        adjacency[int(left)].append(int(right))
        adjacency[int(right)].append(int(left))
    seen = np.zeros(len(mask), dtype=bool)
    groups: list[np.ndarray] = []
    for start in np.flatnonzero(mask):
        if seen[start]:
            continue
        queue = deque([int(start)])
        seen[start] = True
        group: list[int] = []
        while queue:
            node = queue.popleft()
            group.append(node)
            for neighbour in adjacency[node]:
                if mask[neighbour] and not seen[neighbour]:
                    seen[neighbour] = True
                    queue.append(neighbour)
        groups.append(np.asarray(group, dtype=np.int64))
    return groups


def _front_mask(composition: np.ndarray, edges: np.ndarray) -> np.ndarray:
    resolved = np.clip(1.0 - composition[:, -1], 0.0, 1.0)
    tumour = composition[:, 0]
    purity = tumour / np.maximum(resolved, 1e-6)
    core = (
        (resolved >= FRONT["minimum_resolved_fraction"])
        & (tumour >= FRONT["minimum_tumour_membership"])
        & (purity >= FRONT["minimum_tumour_purity_among_resolved"])
    )
    retained = np.zeros(len(core), dtype=bool)
    for group in _components(core, edges):
        if len(group) >= FRONT["minimum_core_component_regions"]:
            retained[group] = True
    core = retained
    adjacency = [[] for _ in range(len(core))]
    boundary = np.zeros(len(core), dtype=bool)
    for left, right in edges.T:
        left, right = int(left), int(right)
        adjacency[left].append(right)
        adjacency[right].append(left)
        if core[left] != core[right]:
            boundary[left] = True
            boundary[right] = True
    front = boundary.copy()
    frontier = np.flatnonzero(boundary).tolist()
    for _ in range(FRONT["front_expansion_hops"]):
        expanded = set(frontier)
        for node in frontier:
            expanded.update(adjacency[node])
        frontier = sorted(expanded)
        front[frontier] = True
    return front


def _balanced_indices(lengths: list[int], seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    quota = max(1, min(MAXIMUM_TOKENS_PER_SLIDE, MAXIMUM_TOKENS_PER_PATIENT // len(lengths)))
    chosen, remainder = [], []
    for length in lengths:
        order = rng.permutation(length)
        take = min(length, quota)
        chosen.append(order[:take])
        remainder.append(order[take : min(length, MAXIMUM_TOKENS_PER_SLIDE)])
    capacity = MAXIMUM_TOKENS_PER_PATIENT - sum(map(len, chosen))
    pool = [(slide, int(value)) for slide, values in enumerate(remainder) for value in values]
    rng.shuffle(pool)
    for slide, value in pool[: max(0, capacity)]:
        chosen[slide] = np.append(chosen[slide], value)
    return [np.sort(values.astype(np.int64, copy=False)) for values in chosen]


def _load_slide(path: Path, selected_names: list[str]) -> dict:
    with h5py.File(path, "r") as handle:
        if not bool(handle.attrs.get("complete", 0)):
            raise ValueError(f"Incomplete topology-token file: {path}")
        names = _decode(handle["tokens/feature_names"][:])
        lookup = {name: index for index, name in enumerate(names)}
        missing = [name for name in selected_names if name not in lookup]
        if missing:
            raise ValueError(f"Feature contract mismatch in {path}; missing {missing[:5]}")
        columns = np.asarray([lookup[name] for name in selected_names], dtype=np.int64)
        return {
            "sample_id": str(handle.attrs.get("sample_id", path.stem)),
            "case_id": str(handle.attrs.get("case_id", handle.attrs.get("sample_id", path.stem))),
            "features": handle["tokens/features"][:, columns].astype(np.float32),
            "composition": handle["tokens/nuisance_cell_composition"][:].astype(np.float32),
            "edges": handle["graph/edge_index"][:].astype(np.int64),
        }


def _patient_bag(slides: list[dict], scope: str, case_id: str) -> tuple[np.ndarray, np.ndarray]:
    pools, retained = [], []
    for slide in slides:
        index = (
            np.flatnonzero(_front_mask(slide["composition"], slide["edges"]))
            if scope == "front"
            else np.arange(len(slide["features"]), dtype=np.int64)
        )
        if scope == "front" and len(index) < FRONT["minimum_front_tokens_per_slide"]:
            continue
        if len(index):
            pools.append(index)
            retained.append(slide)
    if not pools:
        raise ValueError(f"No eligible {scope} regions for {case_id}")
    selected = _balanced_indices(
        [len(pool) for pool in pools], _stable_seed(f"{case_id}|{scope}", SAMPLING_SEED)
    )
    blocks, slide_ids = [], []
    for number, (slide, pool, local) in enumerate(zip(retained, pools, selected)):
        if len(local):
            blocks.append(slide["features"][pool[local]])
            slide_ids.append(np.full(len(local), number, dtype=np.int16))
    values = np.concatenate(blocks)
    if scope == "front" and len(values) < FRONT["minimum_front_tokens_per_patient"]:
        raise ValueError(f"Only {len(values)} tumour-front regions for {case_id}; at least 16 are required")
    return values, np.concatenate(slide_ids)


def _fingerprint(values: np.ndarray, slide_index: np.ndarray, data, branch: str) -> np.ndarray:
    token_mean = data[f"{branch}__token_mean"]
    token_scale = data[f"{branch}__token_scale"]
    centres = data[f"{branch}__prototype_centres"]
    tau = float(data[f"{branch}__tau"][0])
    z = np.clip((values - token_mean) / token_scale, -8.0, 8.0)
    distances = np.sqrt(np.maximum(((z[:, None, :] - centres[None, :, :]) ** 2).sum(axis=2), 0.0))
    logits = -(distances**2) / (2.0 * tau * tau)
    logits -= logits.max(axis=1, keepdims=True)
    membership = np.exp(np.clip(logits, -40.0, 0.0))
    membership /= np.maximum(membership.sum(axis=1, keepdims=True), 1e-8)
    blocks = [membership.mean(axis=0), membership.std(axis=0), np.quantile(membership, 0.90, axis=0)]
    slide_means, slide_q90 = [], []
    for slide in np.unique(slide_index):
        current = membership[slide_index == slide]
        slide_means.append(current.mean(axis=0))
        slide_q90.append(np.quantile(current, 0.90, axis=0))
    slide_means = np.asarray(slide_means)
    blocks.append(np.asarray(slide_q90).max(axis=0))
    blocks.append(slide_means.std(axis=0) if len(slide_means) > 1 else np.zeros(membership.shape[1]))
    return np.concatenate(blocks)


def predict(token_files: list[str | Path], model_dir: str | Path | None = None, case_id: str | None = None) -> dict:
    model_dir = Path(model_dir) if model_dir else Path(__file__).resolve().parent / "weights"
    with np.load(model_dir / "cell_topology_model.npz", allow_pickle=False) as data:
        selected_names = _decode(data["selected_feature_names"])
        slides = [_load_slide(Path(path), selected_names) for path in token_files]
        inferred_ids = {slide["case_id"] for slide in slides}
        resolved_case = case_id or (next(iter(inferred_ids)) if len(inferred_ids) == 1 else "external_case")
        result = {"model": "cell_topology", "case_id": resolved_case, "n_slides": len(slides), "branches": {}}
        branch_risks = []
        for branch, scope in (("whole_tissue", "all"), ("tumour_front", "front")):
            values, slide_index = _patient_bag(slides, scope, resolved_case)
            fingerprint = _fingerprint(values, slide_index, data, branch)
            scaled = (fingerprint - data[f"{branch}__fingerprint_mean"]) / data[f"{branch}__fingerprint_scale"]
            raw = float(np.dot(scaled, data[f"{branch}__cox_coefficients"]))
            risk = (raw - float(data[f"{branch}__risk_median"][0])) / float(data[f"{branch}__risk_iqr"][0])
            result["branches"][branch] = {"risk": risk, "n_regions": int(len(values))}
            branch_risks.append(risk)
        result["risk"] = float(np.mean(branch_risks))
        result["risk_group"] = "high" if result["risk"] >= 0.0 else "low"
        result["threshold"] = 0.0
        return result


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Predict CRC OS risk from cell-topology token files")
    parser.add_argument("tokens", nargs="+")
    parser.add_argument("--model-dir", type=Path)
    parser.add_argument("--case-id")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = predict(args.tokens, args.model_dir, args.case_id)
    text = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
