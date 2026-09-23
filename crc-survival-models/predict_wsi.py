"""Run CRC survival inference from one or more whole-slide images.

The two branches share only the input slide list. Each branch keeps its own
input contract and writes intermediate files under the output directory.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
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


def _external_path(cli_value: Path | None, variable: str, relative_default: str) -> Path:
    if cli_value is not None:
        return cli_value.expanduser().resolve()
    value = os.environ.get(variable)
    return Path(value).expanduser().resolve() if value else (ROOT / relative_default).resolve()


def _external_paths(args: argparse.Namespace) -> dict[str, Path]:
    return {
        "cellvit_root": _external_path(
            args.cellvit_root, "CRC_SURVIVAL_CELLVIT_ROOT", "external/cellvit-plus-plus"
        ),
        "cellvit_checkpoint": _external_path(
            args.cellvit_checkpoint,
            "CRC_SURVIVAL_CELLVIT_CHECKPOINT",
            "external/cellvit/cell_segmentation_backbone.pth",
        ),
        "cellvit_classifiers": _external_path(
            args.cellvit_classifiers,
            "CRC_SURVIVAL_CELLVIT_CLASSIFIERS",
            "external/cellvit/classifiers",
        ),
        "uni_checkpoint": _external_path(
            args.uni_checkpoint,
            "CRC_SURVIVAL_UNI_CHECKPOINT",
            "external/uni/patch_encoder_backbone.pth",
        ),
    }


def _sample_id(path: Path, index: int) -> str:
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", path.stem).strip("._-") or "slide"
    return f"slide_{index:02d}_{stem}"


def _run(command: list[str], env: dict[str, str] | None = None) -> None:
    print("RUN", " ".join(command), flush=True)
    subprocess.run(command, check=True, env=env)


def _preflight(slides: list[Path], models: set[str], external: dict[str, Path]) -> dict:
    missing = [str(path) for path in slides if not path.is_file()]
    patch_formats = {".svs", ".mrxs", ".ndpi", ".tif", ".tiff", ".scn", ".bif", ".ibl"}
    topology_formats = {".svs", ".mrxs", ".ndpi", ".tif", ".tiff", ".scn", ".bif"}
    supported = topology_formats if "topology" in models else patch_formats
    unsupported = [str(path) for path in slides if path.suffix.lower() not in supported]
    requirements = []
    if "topology" in models:
        requirements.extend(
            [
                external["cellvit_root"],
                external["cellvit_checkpoint"],
                *[external["cellvit_classifiers"] / name for name in CLASSIFIER_FILES],
                ROOT / "models/cell_topology/weights/cell_topology_model.npz",
            ]
        )
    if "patch" in models:
        requirements.extend(
            [
                external["uni_checkpoint"],
                ROOT / "models/patch_feature/weights/fold_0/survival_aggregator.pt",
            ]
        )
    absent = [str(path) for path in requirements if not path.exists()]
    return {
        "ok": not missing and not unsupported and not absent,
        "slides": [str(path) for path in slides],
        "missing_slides": missing,
        "unsupported_slide_formats": unsupported,
        "missing_model_files": absent,
        "models": sorted(models),
    }


def _topology(
    slides: list[Path],
    case_id: str,
    output: Path,
    mpp: float | None,
    extra: argparse.Namespace,
    external: dict[str, Path],
) -> dict:
    work = output / "cell_topology"
    work.mkdir(parents=True, exist_ok=True)
    records = []
    for index, slide in enumerate(slides, start=1):
        record = {
            "cohort": "custom",
            "case_id": case_id,
            "sample_id": _sample_id(slide, index),
            "slide": str(slide.resolve()),
            "formal_eligible": True,
            "sensitivity_only": False,
        }
        if mpp is not None:
            record.update({"mpp_x": mpp, "mpp_y": mpp})
        records.append(record)
    manifest = work / "slide_manifest.json"
    manifest.write_text(json.dumps({"records": records}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    env = os.environ.copy()
    env["CRC_SURVIVAL_OUTPUT"] = str(work.resolve())
    env["CRC_SURVIVAL_CELLVIT_ROOT"] = str(external["cellvit_root"])
    env["CRC_SURVIVAL_CELLVIT_CHECKPOINT"] = str(external["cellvit_checkpoint"])
    env["CRC_SURVIVAL_CELLVIT_CLASSIFIERS"] = str(external["cellvit_classifiers"])
    command = [
        sys.executable,
        str(ROOT / "models/cell_topology/src/wsi_pipeline/run_cellvit.py"),
        "--manifest",
        str(manifest),
        "--status",
        str(work / "cell_detection_status.json"),
        "--batch-size",
        str(extra.cell_batch_size),
        "--head-batch-size",
        str(extra.cell_head_batch_size),
        "--tissue-workers",
        str(extra.workers),
        "--minimum-free-gb",
        str(extra.minimum_free_gb),
        "--retry-failed",
    ]
    _run(command, env=env)
    tokens = []
    builder = ROOT / "models/cell_topology/src/wsi_pipeline/build_topology_features.py"
    for record in records:
        nuclei = work / "nuclei_whole_tissue" / f"{record['sample_id']}_nuclei.h5"
        if not nuclei.is_file():
            raise RuntimeError(f"Cell-topology feature generation did not produce {nuclei}")
        _run(
            [
                sys.executable,
                str(builder),
                "--nuclei-h5",
                str(nuclei),
                "--sample-id",
                record["sample_id"],
                "--case-id",
                case_id,
                "--output-root",
                str(work),
            ]
        )
        tokens.append(work / "topology_features/01_region_tokens" / f"{record['sample_id']}.topology_tokens.h5")
    from models.cell_topology.model import predict

    return predict(tokens, case_id=case_id)


def _patch(
    slides: list[Path],
    case_id: str,
    output: Path,
    extra: argparse.Namespace,
    external: dict[str, Path],
) -> dict:
    work = output / "patch_feature"
    cache = work / "patch_cache"
    work.mkdir(parents=True, exist_ok=True)
    rows = ["cohort\tcase_id\tsample_id\tslide_path\teligible"]
    for index, slide in enumerate(slides, start=1):
        rows.append(f"custom\t{case_id}\t{_sample_id(slide, index)}\t{slide.resolve()}\tTrue")
    manifest = work / "slide_manifest.tsv"
    manifest.write_text("\n".join(rows) + "\n", encoding="utf-8")
    env = os.environ.copy()
    env["CRC_SURVIVAL_UNI_CHECKPOINT"] = str(external["uni_checkpoint"])
    _run(
        [
            sys.executable,
            str(ROOT / "models/patch_feature/src/wsi_pipeline/extract_patch_cache.py"),
            "--manifest",
            str(manifest),
            "--output-dir",
            str(cache),
            "--device",
            extra.device,
            "--batch-size",
            str(extra.patch_batch_size),
            "--loader-workers",
            str(extra.workers),
            "--cucim-workers",
            str(extra.workers),
        ],
        env=env,
    )
    caches = [cache / f"{_sample_id(slide, index)}.patch_cache.h5" for index, slide in enumerate(slides, start=1)]
    from models.patch_feature.model import predict

    return predict(caches, device=extra.device, case_id=case_id)


def main() -> None:
    parser = argparse.ArgumentParser(description="CRC WSI survival-risk inference")
    parser.add_argument("slides", nargs="+", type=Path)
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--model", choices=["topology", "patch", "both"], default="both")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mpp", type=float, help="native micrometres per pixel when slide metadata is absent")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--cell-batch-size", type=int, default=8)
    parser.add_argument("--cell-head-batch-size", type=int, default=32768)
    parser.add_argument("--patch-batch-size", type=int, default=160)
    parser.add_argument("--minimum-free-gb", type=float, default=20.0)
    parser.add_argument("--cellvit-root", type=Path, help="local clone of the official CellViT++ repository")
    parser.add_argument("--cellvit-checkpoint", type=Path, help="local compatible CellViT++ backbone checkpoint")
    parser.add_argument("--cellvit-classifiers", type=Path, help="directory containing the seven compatible classifier heads")
    parser.add_argument("--uni-checkpoint", type=Path, help="local compatible UNI-derived coarse-classifier checkpoint")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()

    slides = [path.resolve() for path in args.slides]
    models = {"topology", "patch"} if args.model == "both" else {args.model}
    external = _external_paths(args)
    preflight = _preflight(slides, models, external)
    if args.check_only or not preflight["ok"]:
        print(json.dumps(preflight, ensure_ascii=False, indent=2))
        raise SystemExit(0 if preflight["ok"] else 2)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = {"case_id": args.case_id, "slides": [str(path) for path in slides], "models": {}}
    if "topology" in models:
        results["models"]["cell_topology"] = _topology(
            slides, args.case_id, args.output_dir, args.mpp, args, external
        )
    if "patch" in models:
        results["models"]["patch_feature"] = _patch(
            slides, args.case_id, args.output_dir, args, external
        )
    destination = args.output_dir / "survival_predictions.json"
    destination.write_text(json.dumps(results, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
