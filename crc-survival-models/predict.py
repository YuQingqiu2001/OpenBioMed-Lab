"""Unified command-line interface for the two CRC survival models."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CLASSIFIER_FILES = (
    "consep.pth",
    "lizard.pth",
    "midog.pth",
    "nucls_main.pth",
    "nucls_super.pth",
    "ocelot.pth",
    "panoptils.pth",
)


def configured_path(variable: str, relative_default: str) -> Path:
    value = os.environ.get(variable)
    return Path(value).expanduser().resolve() if value else (ROOT / relative_default).resolve()


def display_path(path: Path) -> str:
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def write(result: dict, output: Path | None) -> None:
    text = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text, encoding="utf-8")
    print(text, end="")


def main() -> None:
    parser = argparse.ArgumentParser(description="Predict overall-survival risk for CRC pathology slides")
    sub = parser.add_subparsers(dest="command", required=True)

    topology = sub.add_parser("topology", help="predict from prepared cell-topology token HDF5 files")
    topology.add_argument("inputs", nargs="+", type=Path)
    topology.add_argument("--case-id")
    topology.add_argument("--output", type=Path)

    patch = sub.add_parser("patch", help="predict from prepared patch-cache HDF5 files")
    patch.add_argument("inputs", nargs="+", type=Path)
    patch.add_argument("--case-id", default="external_case")
    patch.add_argument("--device")
    patch.add_argument("--output", type=Path)

    validate = sub.add_parser("validate-install", help="check required files without running a WSI")
    validate.add_argument("--json", action="store_true")
    validate.add_argument(
        "--raw-wsi",
        action="store_true",
        help="also check user-supplied CellViT++ and UNI paths",
    )

    args = parser.parse_args()
    if args.command == "topology":
        from models.cell_topology.model import predict

        write(predict(args.inputs, case_id=args.case_id), args.output)
    elif args.command == "patch":
        from models.patch_feature.model import predict

        write(predict(args.inputs, device=args.device, case_id=args.case_id), args.output)
    else:
        required = [
            ROOT / "models/cell_topology/weights/cell_topology_model.npz",
            *[
                ROOT / f"models/patch_feature/weights/fold_{fold}/{name}"
                for fold in range(5)
                for name in (
                    "patch_encoder.pt",
                    "spatial_field.pt",
                    "spatial_transforms.pkl",
                    "survival_aggregator.pt",
                    "inference_contract.pkl",
                    "training_manifest.json",
                )
            ],
        ]
        checks = [(path, "project", "file") for path in required]
        if args.raw_wsi:
            cellvit_root = configured_path("CRC_SURVIVAL_CELLVIT_ROOT", "external/cellvit-plus-plus")
            cellvit_checkpoint = configured_path(
                "CRC_SURVIVAL_CELLVIT_CHECKPOINT", "external/cellvit/cell_segmentation_backbone.pth"
            )
            classifiers = configured_path(
                "CRC_SURVIVAL_CELLVIT_CLASSIFIERS", "external/cellvit/classifiers"
            )
            uni_checkpoint = configured_path(
                "CRC_SURVIVAL_UNI_CHECKPOINT", "external/uni/patch_encoder_backbone.pth"
            )
            checks.extend(
                [
                    (cellvit_root, "third_party", "dir"),
                    (cellvit_checkpoint, "third_party", "file"),
                    (uni_checkpoint, "third_party", "file"),
                    *[(classifiers / name, "third_party", "file") for name in CLASSIFIER_FILES],
                ]
            )
        rows = []
        for path, scope, kind in checks:
            exists = path.is_dir() if kind == "dir" else path.is_file()
            rows.append(
                {
                    "path": display_path(path),
                    "scope": scope,
                    "kind": kind,
                    "exists": exists,
                    "bytes": path.stat().st_size if exists and kind == "file" else None,
                }
            )
        result = {"ok": all(row["exists"] for row in rows), "raw_wsi_checked": args.raw_wsi, "files": rows}
        if args.json:
            print(json.dumps(result, ensure_ascii=False, indent=2))
        else:
            for row in rows:
                print(("OK   " if row["exists"] else "MISS ") + row["path"])
            raise SystemExit(0 if result["ok"] else 2)


if __name__ == "__main__":
    main()
