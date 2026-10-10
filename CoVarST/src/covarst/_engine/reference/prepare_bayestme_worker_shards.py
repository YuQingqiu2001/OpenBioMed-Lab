"""Create memory-bounded, exact-row cache shards for parallel BayesTME workers."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def build_tissue_shards(root: Path, workers: int) -> None:
    prepared = root / "00_prepared_top10k"
    contract_path = prepared / "prepared_top10k_contract.json"
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    expected_method = "hest_tissue_marker_fixed_union_slide_balanced_hvg10k_raw_umi_v2"
    if contract.get("method") != expected_method:
        raise RuntimeError(
            f"{root.name}: expected marker-fixed v2 preparation, found {contract.get('method')}"
        )
    if int(contract["genes"]) != 10_000 or not bool(contract["raw_integer_umi"]):
        raise RuntimeError(f"{root.name}: invalid 10K raw-UMI preparation contract")
    panel = pd.read_csv(prepared / "gene_panel_10k.csv")
    if len(panel) != 10_000 or panel["gene"].astype(str).nunique() != 10_000:
        raise RuntimeError(f"{root.name}: panel is not 10,000 unique genes")
    mandatory_rows = panel["selection_reason"].astype(str).eq("mandatory_marker")
    if int(mandatory_rows.sum()) != int(contract["mandatory_markers_retained"]):
        raise RuntimeError(f"{root.name}: marker-first panel contract mismatch")
    if not {"measured_slides", "assay_missing_slides"}.issubset(panel.columns):
        raise RuntimeError(f"{root.name}: marker assay-coverage audit is missing")
    output = root / "00_worker_shards"
    completion = output / "shards_complete.json"
    if completion.is_file():
        existing = json.loads(completion.read_text(encoding="utf-8"))
        if (
            int(existing["workers"]) == workers
            and existing["source_fit_sha256"] == contract["fit_core_sha256"]
        ):
            print(json.dumps({"event": "shards_skip", "tissue": root.name}), flush=True)
            return
        raise RuntimeError(f"stale shard completion contract: {completion}")
    output.mkdir(parents=True, exist_ok=True)

    matrices = {
        name: sparse.load_npz(prepared / f"{name}_core.npz").tocsr()
        for name in ("fit", "tune", "audit")
    }
    shape = matrices["fit"].shape
    if any(matrix.shape != shape for matrix in matrices.values()):
        raise RuntimeError(f"{root.name}: cache matrix shapes differ")
    raw = (matrices["fit"] + matrices["tune"] + matrices["audit"]).tocsr()
    if raw.data.size and not np.allclose(raw.data, np.rint(raw.data)):
        raise RuntimeError(f"{root.name}: non-integer raw molecule cache")
    del raw

    # The first prepared-cache build used a pandas object array for sample IDs.
    # Load that local trusted metadata once, then write every production shard as
    # fixed-width Unicode so the BayesTME runner can keep allow_pickle=False.
    metadata = np.load(prepared / "metadata.npz", allow_pickle=True)
    slide_index = metadata["slide_index"].astype(np.int32)
    sample_ids = np.asarray(metadata["sample_ids"].astype(str).tolist(), dtype=str)
    sources = np.asarray(metadata["sources"].astype(str).tolist(), dtype=str)
    coordinates = np.load(prepared / "coordinates.npy", allow_pickle=False)
    genes = np.load(prepared / "discovery_gene_order.npy", allow_pickle=False)
    if not np.array_equal(genes.astype(str), panel["gene"].astype(str).to_numpy()):
        raise RuntimeError(f"{root.name}: panel CSV and model gene order differ")
    spots = pd.read_csv(prepared / "spot_index.csv")
    if len(spots) != shape[0] or len(slide_index) != shape[0]:
        raise RuntimeError(f"{root.name}: spot metadata row mismatch")

    shard_rows: list[dict[str, object]] = []
    for worker in range(workers):
        final = output / f"worker_{worker:02d}"
        final_contract = final / "shard_contract.json"
        assigned = np.flatnonzero(np.arange(len(sample_ids)) % workers == worker)
        row_mask = np.isin(slide_index, assigned)
        rows = np.flatnonzero(row_mask)
        expected_slides = int(len(assigned))
        if final_contract.is_file():
            prior = json.loads(final_contract.read_text(encoding="utf-8"))
            if prior["source_fit_sha256"] != contract["fit_core_sha256"]:
                raise RuntimeError(f"stale shard: {final}")
            shard_rows.append(prior)
            continue
        partial = output / f"worker_{worker:02d}.partial"
        if partial.exists():
            shutil.rmtree(partial)
        partial.mkdir()
        for name, matrix in matrices.items():
            sparse.save_npz(partial / f"{name}_core.npz", matrix[rows], compressed=True)
        np.save(partial / "coordinates.npy", coordinates[rows])
        np.save(partial / "discovery_gene_order.npy", genes)
        np.savez_compressed(
            partial / "metadata.npz",
            slide_index=slide_index[rows],
            sample_ids=sample_ids,
            sources=sources,
        )
        shard_spots = spots.iloc[rows].copy()
        shard_spots.to_csv(partial / "spot_index.csv", index=False)
        payload = {
            "method": "exact_row_partition_for_parallel_bayestme_v1",
            "tissue_key": root.name,
            "worker_id": worker,
            "workers": workers,
            "assigned_global_slide_indices": assigned.astype(int).tolist(),
            "slides": expected_slides,
            "spots": int(len(rows)),
            "genes": int(shape[1]),
            "preserves_global_slide_index": True,
            "preserves_global_spot_index": True,
            "source_fit_sha256": contract["fit_core_sha256"],
            "fit_core_sha256": sha256(partial / "fit_core.npz"),
        }
        write_json(partial / "shard_contract.json", payload)
        os.replace(partial, final)
        shard_rows.append(payload)
        print(
            json.dumps(
                {
                    "event": "shard_complete",
                    "tissue": root.name,
                    "worker": worker,
                    "slides": expected_slides,
                    "spots": int(len(rows)),
                }
            ),
            flush=True,
        )
    if sum(int(row["spots"]) for row in shard_rows) != shape[0]:
        raise RuntimeError(f"{root.name}: shard rows do not exactly partition the cache")
    write_json(
        completion,
        {
            "method": "exact_row_partition_for_parallel_bayestme_v1",
            "tissue_key": root.name,
            "workers": workers,
            "source_fit_sha256": contract["fit_core_sha256"],
            "slides": int(len(sample_ids)),
            "spots": int(shape[0]),
            "genes": int(shape[1]),
            "shards": shard_rows,
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--tissue", action="append", default=[])
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    workers = int(config["bayestme"]["workers"])
    requested = set(args.tissue)
    for tissue in config["tissues"]:
        if requested and str(tissue["key"]) not in requested:
            continue
        root = Path(str(config["output_root"])) / str(tissue["key"])
        build_tissue_shards(root, workers)


if __name__ == "__main__":
    main()
