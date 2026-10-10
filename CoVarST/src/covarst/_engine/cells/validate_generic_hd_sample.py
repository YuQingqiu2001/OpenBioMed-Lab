"""Independent read-only whole-count integrity of a new cohort, CPU 28-31."""
import argparse
import datetime as dt
import json
import os
from pathlib import Path

os.sched_setaffinity(0, set(range(28, 32)))
for option in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[option] = "2"
from validate_acquired_hd import validate_bin, validate_tiff

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cohort", required=True)
    p.add_argument("--image", type=Path, required=True)
    args = p.parse_args()
    directory = ROOT / "data" / args.cohort
    if not directory.resolve().is_relative_to((ROOT / "data").resolve()):
        raise ValueError("Invalid cohort path")
    report = {"sample": args.cohort, "checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "completed": False, "independent_cell_expression_truth": False, "errors": [], "matrices": []}
    try:
        report["image"] = validate_tiff(args.image)
        if report["image"]["image_role"] != "original_full_resolution":
            raise ValueError("Input is not a full-resolution tissue image")
        for size in (2, 16):
            matrix = validate_bin(directory / "binned_outputs" / f"square_{size:03d}um")
            if matrix["positions_in_tissue"] != matrix["barcodes"]:
                raise ValueError(f"Tissue position count differs from filtered matrix at {size}um")
            report["matrices"].append(matrix)
        report["completed"] = True
    except Exception as error:
        report["errors"].append(str(error))
    output = directory / "semantic_validation.json"
    temporary = output.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    os.replace(temporary, output)
    print(json.dumps({"sample": args.cohort, "completed": report["completed"], "errors": report["errors"]}))
    if not report["completed"]:
        raise RuntimeError(report["errors"])


if __name__ == "__main__":
    main()
