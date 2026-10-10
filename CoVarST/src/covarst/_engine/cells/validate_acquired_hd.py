"""Read-only scientific file validation in existing WSL sc_spatial_env.

No images are emitted. TIFF dimensions/tags and decoded crops are reported as
numbers. Original files remain unchanged; output is a separate JSON report.
"""
from __future__ import annotations
import argparse
import datetime as dt
import json
import os
from pathlib import Path
if hasattr(os, "sched_getaffinity"):
    os.sched_setaffinity(0, set(range(28, 32)))
for option in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[option] = "2"
import h5py
import numpy as np
import pyarrow.parquet as pq
import tifffile
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
NATIVE_BREAST = Path('DATA/binned_outputs')


def validate_tiff(path):
    size = path.stat().st_size
    with tifffile.TiffFile(path) as tf:
        result = {"path": str(path), "bytes": size, "pages": [], "decoded_samples": []}
        for page in tf.pages:
            info = {"shape": list(page.shape), "dtype": str(page.dtype), "compression": str(page.compression),
                    "tiled": bool(page.is_tiled), "segment_count": len(page.dataoffsets)}
            if any(offset < 0 or offset + count > size for offset, count in zip(page.dataoffsets, page.databytecounts)):
                raise ValueError("TIFF segment points outside completed file")
            result["pages"].append(info)
        page = tf.pages[0]
        # tifffile exposes JPEG/LZW segment decoding without loading a whole WSI.
        wanted = sorted(set([0, len(page.dataoffsets) // 2, len(page.dataoffsets) - 1]))
        with path.open("rb") as handle:
            for index in wanted:
                handle.seek(page.dataoffsets[index])
                payload = handle.read(page.databytecounts[index])
                options = {"jpegtables": page.jpegtables} if page.jpegtables is not None else {}
                decoded, indices, shape = page.decode(payload, index, **options)
                if decoded is None or not decoded.size:
                    raise ValueError(f"TIFF segment {index} did not decode")
                result["decoded_samples"].append({"segment": index, "shape": list(decoded.shape),
                                                   "min": int(decoded.min()), "max": int(decoded.max())})
        result["valid"] = True
        result["pixel_read_scope"] = "first_middle_last_compressed_segments; all_segment_bounds_checked"
        result["image_role"] = "original_full_resolution" if min(page.shape[:2]) >= 10000 else "downsampled_auxiliary"
    return result


def validate_bin(directory):
    path = directory / "filtered_feature_bc_matrix.h5"
    spatial = directory / "spatial"
    report = {"path": str(path), "spatial": str(spatial)}
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        shape = m["shape"][:].tolist()
        n_features, n_barcodes = shape
        if len(m["barcodes"]) != n_barcodes or len(m["indptr"]) != n_barcodes + 1:
            raise ValueError("Sparse matrix dimensions do not match barcode/indptr arrays")
        if int(m["indptr"][0]) != 0 or int(m["indptr"][-1]) != len(m["data"]):
            raise ValueError("Sparse matrix indptr endpoints invalid")
        total = 0
        maximum = 0
        for start in range(0, len(m["data"]), 4_000_000):
            values = m["data"][start:start + 4_000_000]
            if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
                raise ValueError("Nonfinite, negative or noninteger gene counts")
            total += int(values.sum(dtype=np.uint64))
            maximum = max(maximum, int(values.max(initial=0)))
        probe = set(x.decode() for x in list(m["barcodes"][:25]) + list(m["barcodes"][-25:]))
        report.update(genes=n_features, barcodes=n_barcodes, nonzero_entries=len(m["data"]),
                      total_umis=total, maximum_entry=maximum, all_count_values_checked=True)
    parquet = pq.ParquetFile(spatial / "tissue_positions.parquet")
    columns = parquet.schema.names
    selected = [x for x in ["barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"] if x in columns]
    remaining = set(probe)
    tissue = 0
    limits = {"pxl_row_in_fullres": [float("inf"), float("-inf")], "pxl_col_in_fullres": [float("inf"), float("-inf")]}
    for batch in parquet.iter_batches(batch_size=262144, columns=selected):
        frame = batch.to_pandas()
        if "barcode" in frame:
            remaining -= set(frame["barcode"].isin(remaining).pipe(lambda found: frame.loc[found, "barcode"]))
        if "in_tissue" in frame:
            tissue += int(frame["in_tissue"].sum())
        for name in limits:
            if name in frame:
                limits[name][0] = min(limits[name][0], float(frame[name].min()))
                limits[name][1] = max(limits[name][1], float(frame[name].max()))
    if remaining:
        raise ValueError(f"Sampled matrix barcodes missing in positions: {len(remaining)}")
    scales = json.loads((spatial / "scalefactors_json.json").read_text())
    report.update(position_rows=parquet.metadata.num_rows, positions_in_tissue=tissue, coordinate_ranges=limits,
                  scalefactors=scales, sampled_barcode_alignment_count=len(probe), valid=True)
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--images-only", action="store_true")
    args = parser.parse_args()
    report = {"checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "images": [], "matrices": [],
              "errors": [], "independent_cell_ground_truth": False}
    original_batch_images = [path for cohort in ("BRCA_FF", "LUAD_FX")
                             for path in (ROOT / "data" / cohort).glob("*image.*")]
    for path in sorted(original_batch_images):
        if path.suffix not in (".tif", ".btf"):
            continue
        try:
            report["images"].append(validate_tiff(path))
        except Exception as e:
            report["errors"].append({"path": str(path), "error": str(e)})
    if not args.images_only:
        for root in [NATIVE_BREAST, ROOT / "data/LUAD_FX/binned_outputs"]:
            for bin_name in ["square_002um", "square_016um"]:
                directory = root / bin_name
                if not (directory / "filtered_feature_bc_matrix.h5").exists():
                    report["errors"].append({"path": str(directory), "error": "not_downloaded_or_extracted_yet"})
                    continue
                try:
                    report["matrices"].append(validate_bin(directory))
                except Exception as e:
                    report["errors"].append({"path": str(directory), "error": str(e)})
    outpath = ROOT / "data" / ("image_validation.json" if args.images_only else "semantic_validation.json")
    outpath.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"report": str(outpath), "images": len(report["images"]), "matrices": len(report["matrices"]),
                      "errors": report["errors"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
