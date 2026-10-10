"""Independently verify the completed clustered top10K BayesTME deliverables."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np


def _atomic_json(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    os.replace(temporary, path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _probability_rows(path: Path) -> tuple[int, float]:
    matrix = np.load(path, mmap_mode="r")
    if matrix.ndim != 2 or not np.isfinite(matrix).all() or np.any(matrix < 0):
        raise RuntimeError(f"invalid probability matrix: {path.name}")
    maximum_error = 0.0
    for start in range(0, matrix.shape[0], 1024):
        stop = min(start + 1024, matrix.shape[0])
        maximum_error = max(
            maximum_error,
            float(np.max(np.abs(np.asarray(matrix[start:stop]).sum(axis=1) - 1.0))),
        )
    return int(matrix.shape[0]), maximum_error


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--integration-dir", type=Path, required=True)
    parser.add_argument("--figures-dir", type=Path, required=True)
    args = parser.parse_args()

    integration = args.integration_dir
    figures = args.figures_dir
    integration_contract = json.loads(
        (integration / "top10k_integration_contract.json").read_text(encoding="utf-8")
    )
    reconstruction_contract = json.loads(
        (integration / "top10k_reconstruction_contract.json").read_text(
            encoding="utf-8"
        )
    )
    figures_contract = json.loads(
        (figures / "program_figures_contract.json").read_text(encoding="utf-8")
    )

    if not bool(integration_contract["production_ready"]):
        raise RuntimeError("integration is not production-ready")
    if int(integration_contract["completed_slides"]) != 177:
        raise RuntimeError("integration does not include all 177 slides")
    if not bool(integration_contract["all_local_programs_assigned_once"]):
        raise RuntimeError("not every local program is assigned exactly once")
    if max(
        float(integration_contract["basis_row_sum_maximum_error"]),
        float(integration_contract["theta_row_sum_maximum_error"]),
    ) > 1.0e-6:
        raise RuntimeError("integration contract has invalid probability normalization")
    if reconstruction_contract["standardization"] != (
        "log1p(1e4 * row-normalized expression probability)"
    ):
        raise RuntimeError("standardization contract differs from the required logCP10K")
    if int(reconstruction_contract["local_slides"]) != 177:
        raise RuntimeError("reconstruction omits slides")
    if max(
        float(reconstruction_contract["global_reload_maximum_absolute_error"]),
        float(reconstruction_contract["local_reload_maximum_absolute_error"]),
        float(reconstruction_contract["local_factorization_maximum_absolute_error"]),
    ) > 1.0e-6:
        raise RuntimeError("reconstruction contract has failed a reload/factorization check")

    paths = {
        "basis_sha256": integration / "reference_probability_global.npy",
        "theta_sha256": integration / "theta_rna_global.npy",
        "gene_order_sha256": integration / "gene_order.npy",
        "local_logcp10k_sha256": integration / "reconstructed_logcp10k_local_10k.npy",
        "global_logcp10k_sha256": integration / "reconstructed_logcp10k_global_10k.npy",
    }
    hashes = {name: _sha256(path) for name, path in paths.items()}
    mismatches = [
        name
        for name, value in hashes.items()
        if value != reconstruction_contract[name]
    ]
    if mismatches:
        raise RuntimeError(f"reconstruction hash mismatch: {', '.join(mismatches)}")

    basis_rows, basis_error = _probability_rows(paths["basis_sha256"])
    theta_rows, theta_error = _probability_rows(paths["theta_sha256"])
    if max(basis_error, theta_error) > 1.0e-6:
        raise RuntimeError("reloaded basis or theta is not row-normalized")

    required_figures = [str(value) for value in figures_contract["figures"]]
    for name in required_figures:
        pdf = figures / name
        png = figures / name.replace(".pdf", ".png")
        if not pdf.is_file() or pdf.stat().st_size < 1024:
            raise RuntimeError(f"missing or empty PDF figure: {name}")
        if pdf.read_bytes()[:4] != b"%PDF":
            raise RuntimeError(f"invalid PDF signature: {name}")
        if not png.is_file() or png.stat().st_size < 1024:
            raise RuntimeError(f"missing or empty PNG preview: {png.name}")

    audit = {
        "status": "complete",
        "method": "independent_pure_bayestme_top10k_completion_audit_v1",
        "completed_slides": 177,
        "global_programs": int(integration_contract["global_programs"]),
        "local_programs": int(integration_contract["local_programs"]),
        "spots": int(reconstruction_contract["spots"]),
        "standardization": reconstruction_contract["standardization"],
        "basis_rows": basis_rows,
        "basis_rowsum_maximum_error": basis_error,
        "theta_rows": theta_rows,
        "theta_rowsum_maximum_error": theta_error,
        "hashes_verified": sorted(hashes),
        "figures_verified": required_figures,
    }
    _atomic_json(integration / "top10k_completion_audit.json", audit)
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == "__main__":
    main()
