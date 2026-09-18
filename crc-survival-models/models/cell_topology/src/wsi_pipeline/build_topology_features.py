"""Convert one typed-nucleus HDF5 file into a topology-token HDF5 file."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from extract_tokens import config_sha256, extract_one
from graph import build_slide_cache


def build(nuclei_h5: Path, sample_id: str, case_id: str, output_root: Path) -> Path:
    graph_dir = output_root / "graph_cache"
    graph_dir.mkdir(parents=True, exist_ok=True)
    graph_path = graph_dir / f"{sample_id}_cell_graph.h5"
    if not graph_path.exists():
        build_slide_cache(nuclei_h5, graph_path, k=8, max_edge_um=50.0, core_um=250.0)

    config_path = Path(__file__).resolve().parent / "token_config.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    config["input"]["graph_cache_dir"] = str(graph_dir)
    config["output_root"] = str(output_root / "topology_features")
    result = extract_one(
        {"sample_id": sample_id, "case_id": case_id},
        config,
        config_sha256(config),
    )
    if result["status"] not in {"complete", "skipped_complete"}:
        raise RuntimeError(f"Topology feature extraction failed: {result}")
    token_path = Path(config["output_root"]) / "01_region_tokens" / f"{sample_id}.topology_tokens.h5"
    if not token_path.exists():
        raise FileNotFoundError(token_path)
    return token_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Build cell-topology features from a typed-nucleus HDF5")
    parser.add_argument("--nuclei-h5", type=Path, required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    path = build(args.nuclei_h5, args.sample_id, args.case_id, args.output_root)
    print(json.dumps({"status": "complete", "topology_tokens": str(path)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
