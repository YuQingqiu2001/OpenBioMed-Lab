"""Standalone HE-only inference for the pan-atlas G384 v4 mapper.

This program deliberately does not load consensus theta, local BayesTME theta,
local expression bases, counts, or any expression reference.  Those inputs are
reserved for the separate evaluation program.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

WORKSPACE = Path(__file__).resolve().parents[3]
PIPELINE_SRC = Path(__file__).resolve().parents[1] / "src"
for module_root in (WORKSPACE, PIPELINE_SRC):
    if str(module_root) not in sys.path:
        sys.path.insert(0, str(module_root))

from openst_final.he_uni_mamba_bnn_inr import normalize_slide_coordinates
from run_virchow2_bayestme_global_mapper import _morphology_graph
from run_virchow2_fseg_k64_mapper import BoundedINR
from train_virchow2_pan_atlas_consensus_mapper import (
    PanAtlasDualMapper,
    _resolve_feature_path,
    _sha256,
    _to_device_graph_pair,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--program-group-index", type=Path, required=True)
    parser.add_argument("--master-spot-index", type=Path, required=True)
    parser.add_argument("--split-manifest", type=Path, required=True)
    parser.add_argument("--feature-root", type=Path)
    parser.add_argument("--minimum-feature-coverage", type=float, default=0.90)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--split", default="test")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.time()
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is required")
    device = torch.device("cuda")
    torch.set_float32_matmul_precision("high")
    checkpoint = torch.load(args.checkpoint, map_location="cpu")
    program_groups = np.load(args.program_group_index, allow_pickle=False)
    model = PanAtlasDualMapper(
        int(checkpoint["input_dim"]),
        hidden=int(checkpoint["hidden"]),
        programs=int(checkpoint["programs"]),
        program_groups=program_groups,
        graph_layers=int(checkpoint["graph_layers"]),
        context_graph_layers=int(checkpoint["context_graph_layers"]),
        graph_heads=int(checkpoint["graph_heads"]),
        dropout=0.10,
        type_lift_alpha=float(checkpoint["type_lift_alpha"]),
        sources=int(checkpoint["sources"]),
        technologies=int(checkpoint["technologies"]),
        support_mass=float(checkpoint["support_mass"]),
        hierarchical_logits=bool(checkpoint.get("hierarchical_logits", False)),
    ).to(device)
    model.load_state_dict(checkpoint["model"], strict=True)
    model.eval()

    inr = None
    if bool(checkpoint.get("inr_accepted", False)) and checkpoint.get("inr") is not None:
        inr = BoundedINR(
            int(checkpoint["hidden"]), int(checkpoint["programs"]) * 2
        ).to(device)
        inr.load_state_dict(checkpoint["inr"], strict=True)
        inr.eval()

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(f"refusing non-empty output {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    master = pd.read_csv(args.master_spot_index)
    split = pd.read_csv(args.split_manifest)
    rows = split.loc[split["split"].astype(str).eq(str(args.split))]
    if not len(rows):
        raise RuntimeError(f"no rows assigned to split {args.split}")

    predicted_rna: list[np.ndarray] = []
    predicted_composition: list[np.ndarray] = []
    predicted_support: list[np.ndarray] = []
    spot_tables: list[pd.DataFrame] = []
    slide_records: list[dict[str, object]] = []
    prediction_offset = 0

    with torch.no_grad():
        for row in rows.itertuples(index=False):
            source = str(row.source)
            sample_id = str(row.reference_sample_id)
            block = master.loc[
                master["source"].astype(str).eq(source)
                & master["sample_id"].astype(str).eq(sample_id)
            ].sort_values("local_index", kind="stable")
            if not len(block):
                raise RuntimeError(f"master spots absent for {sample_id}")
            feature_path = _resolve_feature_path(str(row.feature_path), args.feature_root)
            with np.load(feature_path, allow_pickle=False) as archive:
                features_all = np.asarray(archive["features"], dtype=np.float16)
                barcodes = np.asarray(archive["barcode"]).astype(str)
                coordinates_all = np.asarray(archive["coords"], dtype=np.float32)
            lookup = {barcode: index for index, barcode in enumerate(barcodes)}
            matched = block["barcode"].astype(str).isin(lookup)
            coverage = float(matched.mean())
            if coverage < float(args.minimum_feature_coverage):
                raise RuntimeError(
                    f"Virchow2/spot coverage below {args.minimum_feature_coverage:.1%} "
                    f"for {sample_id}: {coverage:.3%}"
                )
            block = block.loc[matched].copy()
            feature_rows = np.asarray(
                [lookup[value] for value in block["barcode"].astype(str)],
                dtype=np.int64,
            )
            features_np = features_all[feature_rows]
            coordinates_np = coordinates_all[feature_rows]
            local_graph = _morphology_graph(
                coordinates_np,
                np.asarray(features_np, dtype=np.float32),
                neighbors=int(checkpoint["neighbors"]),
            )
            context_graph = _morphology_graph(
                coordinates_np,
                np.asarray(features_np, dtype=np.float32),
                neighbors=int(checkpoint["context_neighbors"]),
            )
            features = torch.from_numpy(
                np.asarray(features_np, dtype=np.float32)
            ).to(device)
            graphs = _to_device_graph_pair((local_graph, context_graph), device)
            rna_logits, composition_logits, context, auxiliary = model(
                features,
                graphs,
                support_gamma=float(checkpoint.get("support_gamma", 1.0)),
                domain_strength=0.0,
            )
            if inr is not None:
                coordinates = torch.from_numpy(
                    normalize_slide_coordinates(coordinates_np)
                ).to(device)
                correction = inr(coordinates, context)
                rna_correction, composition_correction = correction.split(
                    int(checkpoint["programs"]), dim=1
                )
                rna_logits = rna_logits + rna_correction
                composition_logits = composition_logits + composition_correction

            rna_np = torch.softmax(rna_logits, dim=1).cpu().numpy().astype(np.float32)
            composition_np = (
                torch.softmax(composition_logits, dim=1)
                .cpu()
                .numpy()
                .astype(np.float32)
            )
            support_np = (
                torch.sigmoid(auxiliary["support_logits"])
                .cpu()
                .numpy()
                .astype(np.float32)
            )
            predicted_rna.append(rna_np)
            predicted_composition.append(composition_np)
            predicted_support.append(support_np)
            slide_key = f"{source}:{sample_id.split(':', 1)[-1]}"
            block = block.assign(
                inference_slide_key=slide_key,
                inference_split=str(args.split),
                prediction_row=np.arange(
                    prediction_offset, prediction_offset + len(block), dtype=np.int64
                ),
            )
            prediction_offset += len(block)
            spot_tables.append(block)
            slide_records.append(
                {
                    "slide_key": slide_key,
                    "sample_id": sample_id,
                    "source": source,
                    "patient": str(row.patient),
                    "spots": int(len(block)),
                    "coverage": coverage,
                    "feature_path": str(feature_path),
                }
            )

    np.save(
        args.output_dir / "theta_rna_predicted.npy",
        np.concatenate(predicted_rna, axis=0),
    )
    np.save(
        args.output_dir / "theta_composition_predicted.npy",
        np.concatenate(predicted_composition, axis=0),
    )
    np.save(
        args.output_dir / "program_support_probability.npy",
        np.concatenate(predicted_support, axis=0),
    )
    pd.concat(spot_tables, ignore_index=True).to_csv(
        args.output_dir / "inference_spot_index.csv", index=False
    )
    pd.DataFrame(slide_records).to_csv(
        args.output_dir / "inference_slides.csv", index=False
    )
    payload = {
        "method": f"standalone_he_only_pan_atlas_g{int(checkpoint['programs'])}_inference",
        "split": str(args.split),
        "checkpoint": str(args.checkpoint),
        "checkpoint_sha256": _sha256(args.checkpoint),
        "spots": int(prediction_offset),
        "slides": int(len(slide_records)),
        "programs": int(checkpoint["programs"]),
        "hierarchical_logits": bool(
            checkpoint.get("hierarchical_logits", False)
        ),
        "program_group_method": str(
            checkpoint.get("program_group_method", "expression")
        ),
        "inr_accepted": bool(checkpoint.get("inr_accepted", False)),
        "transcriptomic_inputs_loaded": False,
        "runtime_seconds": time.time() - started,
        "peak_cuda_memory_bytes": int(torch.cuda.max_memory_allocated()),
    }
    (args.output_dir / "inference_manifest.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
