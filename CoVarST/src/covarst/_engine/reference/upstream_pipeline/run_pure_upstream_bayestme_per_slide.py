"""Run the installed, unmodified BayesTME 1.0.0 independently per slide.

This is deliberately a thin orchestration layer.  It does not patch BayesTME,
replace its variational guide, initialise its latent variables externally, or
merge/prune local components.  Every candidate K and every final fit calls the
public upstream functions directly.

The input count matrix is the exact sum of the previously serialized
fit/tune/audit molecule splits.  Those files are used only as a lossless raw
UMI cache; no old basis, theta, program labels, or checkpoints are read.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import json
import os
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyro
import scipy
import torch
from bayestme import cv_likelihoods
from bayestme.common import Layout
from bayestme.data import SpatialExpressionDataset
from bayestme.deconvolution import sample_from_posterior
from bayestme.phenotype_selection import (
    create_folds,
    run_phenotype_selection_single_job,
)
from bayestme.svi.deconvolution import BayesTME_VI
from bayestme.utils import get_edges
from scipy import sparse


EPS = 1.0e-12
EXPECTED_BAYESTME_VERSION = "1.0.0"


def _safe_name(value: str) -> str:
    return (
        str(value)
        .replace(":", "__")
        .replace("/", "_")
        .replace("\\", "_")
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, allow_nan=False),
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _normalise_rows(values: np.ndarray) -> np.ndarray:
    values = np.maximum(np.asarray(values, dtype=np.float64), 0.0)
    return values / np.maximum(values.sum(axis=1, keepdims=True), EPS)


def _gene_pcc(target: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    target = np.log1p(1.0e4 * _normalise_rows(target))
    prediction = np.log1p(1.0e4 * _normalise_rows(prediction))
    target -= target.mean(axis=0, keepdims=True)
    prediction -= prediction.mean(axis=0, keepdims=True)
    denominator = np.sqrt(
        np.square(target).sum(axis=0)
        * np.square(prediction).sum(axis=0)
    )
    correlation = np.full(target.shape[1], np.nan, dtype=np.float64)
    valid = denominator > EPS
    correlation[valid] = (
        target[:, valid] * prediction[:, valid]
    ).sum(axis=0) / denominator[valid]
    finite = correlation[np.isfinite(correlation)]
    if len(finite) == 0:
        return {
            "median": float("nan"),
            "mean": float("nan"),
            "p10": float("nan"),
            "evaluable_genes": 0,
        }
    return {
        "median": float(np.median(finite)),
        "mean": float(np.mean(finite)),
        "p10": float(np.quantile(finite, 0.10)),
        "evaluable_genes": int(len(finite)),
    }


def _select_shared_genes(
    raw: sparse.csr_matrix,
    slide_index: np.ndarray,
    sources: np.ndarray,
    maximum_genes: int,
) -> np.ndarray:
    """Select a source/slide-balanced count-derived common gene panel."""
    n_genes = raw.shape[1]
    if maximum_genes >= n_genes:
        return np.arange(n_genes, dtype=np.int64)
    source_names = np.unique(sources.astype(str))
    source_scores: list[np.ndarray] = []
    detection = np.asarray((raw > 0).sum(axis=0)).ravel()
    for source in source_names:
        slide_scores: list[np.ndarray] = []
        source_slides = np.flatnonzero(sources.astype(str) == source)
        for slide in source_slides:
            rows = np.flatnonzero(slide_index == slide)
            block = raw[rows].tocsr().astype(np.float64)
            library = np.asarray(block.sum(axis=1)).ravel()
            transformed = block.multiply(
                (1.0e4 / np.maximum(library, 1.0))[:, None]
            ).tocsr()
            transformed.data = np.log1p(transformed.data)
            mean = np.asarray(transformed.mean(axis=0)).ravel()
            second = np.asarray(
                transformed.multiply(transformed).mean(axis=0)
            ).ravel()
            variance = np.maximum(second - np.square(mean), 0.0)
            order = np.argsort(variance, kind="stable")
            rank = np.empty(n_genes, dtype=np.float64)
            rank[order] = np.linspace(0.0, 1.0, n_genes)
            slide_scores.append(rank)
        source_scores.append(np.mean(slide_scores, axis=0))
    score = np.mean(source_scores, axis=0)
    score[detection < max(10, int(0.001 * raw.shape[0]))] = -np.inf
    selected = np.argsort(-score, kind="stable")[:maximum_genes]
    return np.sort(selected.astype(np.int64))


def _audit_upstream_install() -> dict[str, object]:
    version = importlib.metadata.version("bayestme")
    if version != EXPECTED_BAYESTME_VERSION:
        raise RuntimeError(
            f"expected bayestme {EXPECTED_BAYESTME_VERSION}, found {version}"
        )
    module_path = Path(inspect.getsourcefile(BayesTME_VI) or "").resolve()
    guide_source = inspect.getsource(BayesTME_VI.spatial_guide)
    symmetric_guide = (
        "torch.ones(self.n_celltypes, self.n_genes)" in guide_source
        and "torch.ones(self.N, self.n_celltypes)" in guide_source
    )
    if not symmetric_guide:
        raise RuntimeError(
            "BayesTME spatial guide no longer matches the unmodified "
            "upstream all-ones guide"
        )
    forbidden_loaded = sorted(
        name
        for name in tuple(__import__("sys").modules)
        if "bayestme_asymmetric_init" in name
    )
    if forbidden_loaded:
        raise RuntimeError(
            f"forbidden BayesTME patch module loaded: {forbidden_loaded}"
        )
    return {
        "bayestme_version": version,
        "bayestme_svi_source": str(module_path),
        "bayestme_svi_source_sha256": _sha256(module_path),
        "symmetric_all_ones_upstream_guide": True,
        "runtime_patch_modules_loaded": forbidden_loaded,
        "python": __import__("sys").version,
        "torch": torch.__version__,
        "pyro": pyro.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "cuda_visible_to_torch": bool(torch.cuda.is_available()),
        "upstream_bayestme_device": (
            "cpu; BayesTME 1.0.0 constructs CPU tensors internally"
        ),
    }


def _build_dataset(
    counts: np.ndarray,
    coordinates: np.ndarray,
    genes: np.ndarray,
    barcodes: np.ndarray,
) -> SpatialExpressionDataset:
    positions = np.rint(np.asarray(coordinates)).astype(np.int64)
    edges = get_edges(positions, layout=Layout.IRREGULAR)
    return SpatialExpressionDataset.from_arrays(
        raw_counts=np.asarray(counts, dtype=np.float32),
        positions=positions,
        tissue_mask=np.ones(len(counts), dtype=bool),
        gene_names=np.asarray(genes, dtype=str),
        layout=Layout.IRREGULAR,
        edges=np.asarray(edges, dtype=np.int64),
        barcodes=np.asarray(barcodes, dtype=str),
    )


def _run_official_k_selection(
    dataset: SpatialExpressionDataset,
    output_dir: Path,
    *,
    k_min: int,
    k_max: int,
    spatial_smoothing: float,
    n_splits: int,
    n_svi_steps: int,
    n_samples: int,
    seed: int,
) -> tuple[int, list[dict[str, object]]]:
    state_path = output_dir / "phenotype_selection.json"
    if state_path.exists():
        state = json.loads(state_path.read_text(encoding="utf-8"))
        records = list(state.get("candidates", []))
    else:
        records = []
    completed_k = {int(record["k"]) for record in records}
    fold_mask = next(
        create_folds(
            dataset,
            n_fold=1,
            n_splits=int(n_splits),
        )
    )
    np.save(output_dir / "phenotype_selection_holdout_mask.npy", fold_mask)

    for k in range(int(k_min), int(k_max) + 1):
        if k in completed_k:
            continue
        started = time.time()
        result = run_phenotype_selection_single_job(
            spatial_smoothing_parameter=float(spatial_smoothing),
            n_components=int(k),
            mask=fold_mask,
            fold_number=0,
            stdata=dataset,
            n_samples=int(n_samples),
            n_svi_steps=int(n_svi_steps),
            use_spatial_guide=True,
            rng=np.random.default_rng(int(seed + 1009 * k)),
        )
        record = {
            "k": int(k),
            "lambda": float(spatial_smoothing),
            "fold": 0,
            "heldout_log_likelihood_mean": float(
                np.nanmean(result.log_lh_test_trace)
            ),
            "train_log_likelihood_mean": float(
                np.nanmean(result.log_lh_train_trace)
            ),
            "heldout_log_likelihood_sd": float(
                np.nanstd(result.log_lh_test_trace)
            ),
            "posterior_draws": int(result.cell_prob_trace.shape[0]),
            "runtime_seconds": float(time.time() - started),
            "seed": int(seed + 1009 * k),
            "n_svi_steps": int(n_svi_steps),
            "n_samples": int(n_samples),
        }
        records.append(record)
        records.sort(key=lambda item: int(item["k"]))
        _atomic_json(
            state_path,
            {
                "method": (
                    "unmodified BayesTME "
                    "run_phenotype_selection_single_job"
                ),
                "selection_rule": (
                    "upstream mean heldout likelihood maximum"
                ),
                "k_min": int(k_min),
                "k_max": int(k_max),
                "lambda_values": [float(spatial_smoothing)],
                "n_fold": 1,
                "n_splits": int(n_splits),
                "candidates": records,
                "complete": False,
            },
        )
        print(json.dumps(record), flush=True)

    k_values = list(range(int(k_min), int(k_max) + 1))
    likelihoods = np.full((2, len(k_values), 1, 1), np.nan)
    record_by_k = {int(record["k"]): record for record in records}
    for index, k in enumerate(k_values):
        likelihoods[0, index, 0, 0] = float(
            record_by_k[k]["train_log_likelihood_mean"]
        )
        likelihoods[1, index, 0, 0] = float(
            record_by_k[k]["heldout_log_likelihood_mean"]
        )
    selected_k = int(
        cv_likelihoods.get_max_likelihood_n_components(
            likelihoods,
            k_values,
        )
    )
    _atomic_json(
        state_path,
        {
            "method": (
                "unmodified BayesTME "
                "run_phenotype_selection_single_job"
            ),
            "selection_rule": "upstream mean heldout likelihood maximum",
            "selected_k": selected_k,
            "selected_lambda": float(spatial_smoothing),
            "k_min": int(k_min),
            "k_max": int(k_max),
            "lambda_values": [float(spatial_smoothing)],
            "n_fold": 1,
            "n_splits": int(n_splits),
            "boundary_maximum": bool(selected_k in {k_min, k_max}),
            "candidates": records,
            "complete": True,
        },
    )
    np.save(output_dir / "phenotype_selection_likelihoods.npy", likelihoods)
    return selected_k, records


def _run_slide(
    *,
    slide: int,
    rows: np.ndarray,
    raw: sparse.csr_matrix,
    coordinates: np.ndarray,
    spot_index: pd.DataFrame,
    genes: np.ndarray,
    sample_id: str,
    source: str,
    output_dir: Path,
    k_min: int,
    k_max: int,
    selection_steps: int,
    selection_samples: int,
    selection_splits: int,
    final_steps: int,
    final_samples: int,
    spatial_smoothing: float,
    seed: int,
) -> dict[str, object]:
    started = time.time()
    dense_counts = raw[rows].toarray().astype(np.float32)
    block_index = spot_index.iloc[rows].copy()
    expected_local = np.arange(len(rows), dtype=np.int64)
    observed_local = block_index["local_index"].to_numpy(dtype=np.int64)
    if not np.array_equal(observed_local, expected_local):
        raise RuntimeError(f"{sample_id}: non-canonical local spot order")
    dataset = _build_dataset(
        dense_counts,
        coordinates[rows],
        genes,
        block_index["barcode"].astype(str).to_numpy(),
    )
    selected_k, candidates = _run_official_k_selection(
        dataset,
        output_dir,
        k_min=int(k_min),
        k_max=int(k_max),
        spatial_smoothing=float(spatial_smoothing),
        n_splits=int(selection_splits),
        n_svi_steps=int(selection_steps),
        n_samples=int(selection_samples),
        seed=int(seed),
    )

    final_started = time.time()
    result = sample_from_posterior(
        data=dataset,
        n_components=int(selected_k),
        spatial_smoothing_parameter=float(spatial_smoothing),
        n_samples=int(final_samples),
        n_svi_steps=int(final_steps),
        use_spatial_guide=True,
        rng=np.random.default_rng(int(seed + 9_999_991)),
    )
    psi_draws = np.asarray(result.cell_prob_trace, dtype=np.float64)
    beta_draws = np.asarray(result.beta_trace, dtype=np.float64)
    basis_draws = np.asarray(result.expression_trace, dtype=np.float64)
    theta_composition_draws = psi_draws / np.maximum(
        psi_draws.sum(axis=2, keepdims=True),
        EPS,
    )
    theta_rna_draws = psi_draws * beta_draws[:, None, :]
    theta_rna_draws /= np.maximum(
        theta_rna_draws.sum(axis=2, keepdims=True),
        EPS,
    )
    basis_probability_draws = basis_draws / np.maximum(
        basis_draws.sum(axis=2, keepdims=True),
        EPS,
    )
    theta_composition = _normalise_rows(
        theta_composition_draws.mean(axis=0)
    ).astype(np.float32)
    theta_rna = _normalise_rows(
        theta_rna_draws.mean(axis=0)
    ).astype(np.float32)
    basis_probability = _normalise_rows(
        basis_probability_draws.mean(axis=0)
    ).astype(np.float32)
    upstream_reconstruction = np.asarray(
        result.nb_probs,
        dtype=np.float64,
    ).mean(axis=0)
    upstream_reconstruction = _normalise_rows(
        upstream_reconstruction
    ).astype(np.float32)
    factorized_reconstruction = _normalise_rows(
        theta_rna @ basis_probability
    ).astype(np.float32)
    count_probability = _normalise_rows(dense_counts).astype(np.float32)
    library_size = dense_counts.sum(axis=1).astype(np.int64)
    component_mass = (
        theta_rna.astype(np.float64)
        * np.maximum(library_size, 1)[:, None]
    ).sum(axis=0)
    component_mass /= np.maximum(component_mass.sum(), EPS)

    np.save(
        output_dir / "cell_prob_trace.npy",
        np.asarray(result.cell_prob_trace, dtype=np.float16),
    )
    np.save(
        output_dir / "expression_trace.npy",
        np.asarray(result.expression_trace, dtype=np.float16),
    )
    np.save(
        output_dir / "beta_trace.npy",
        np.asarray(result.beta_trace, dtype=np.float32),
    )
    np.save(
        output_dir / "cell_num_total_trace.npy",
        np.asarray(result.cell_num_total_trace, dtype=np.float32),
    )
    np.save(
        output_dir / "losses.npy",
        np.asarray(result.losses, dtype=np.float64),
    )
    np.save(output_dir / "theta_composition_mean.npy", theta_composition)
    np.save(output_dir / "theta_rna_mean.npy", theta_rna)
    np.save(output_dir / "basis_probability_mean.npy", basis_probability)
    np.save(
        output_dir / "reconstructed_probability_upstream.npy",
        upstream_reconstruction,
    )
    np.save(
        output_dir / "reconstructed_probability_factorized.npy",
        factorized_reconstruction,
    )
    np.save(output_dir / "library_size.npy", library_size)
    np.save(output_dir / "coordinates.npy", coordinates[rows].astype(np.float32))
    block_index.to_csv(output_dir / "spot_index.csv", index=False)
    pd.DataFrame(
        {
            "slide_index": int(slide),
            "sample_id": sample_id,
            "source": source,
            "local_program_index": np.arange(selected_k, dtype=np.int32),
            "rna_mass_fraction": component_mass,
        }
    ).to_csv(output_dir / "local_program_index.csv", index=False)

    metrics = {
        "method": "unmodified_upstream_bayestme_1.0.0_per_slide",
        "slide_index": int(slide),
        "sample_id": sample_id,
        "source": source,
        "spots": int(len(rows)),
        "genes": int(len(genes)),
        "selected_k": int(selected_k),
        "spatial_smoothing_parameter": float(spatial_smoothing),
        "phenotype_selection_candidates": int(len(candidates)),
        "phenotype_selection_steps_per_candidate": int(selection_steps),
        "final_svi_steps": int(final_steps),
        "final_posterior_samples": int(final_samples),
        "final_fit_runtime_seconds": float(time.time() - final_started),
        "runtime_seconds": float(time.time() - started),
        "raw_count_reconstruction_upstream": _gene_pcc(
            count_probability,
            upstream_reconstruction,
        ),
        "raw_count_reconstruction_factorized_mean": _gene_pcc(
            count_probability,
            factorized_reconstruction,
        ),
        "upstream_vs_factorized_gene_pcc": _gene_pcc(
            upstream_reconstruction,
            factorized_reconstruction,
        ),
        "component_mass_fraction": component_mass.tolist(),
        "local_components_pruned": 0,
        "local_components_merged": 0,
        "old_basis_or_theta_read": 0,
        "images_read": 0,
        "partitions_read": 0,
        "upstream_source_or_guide_modified": False,
    }
    _atomic_json(output_dir / "complete.json", metrics)
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared-cache", type=Path, required=True)
    parser.add_argument("--spot-index", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--maximum-genes", type=int, default=1000)
    parser.add_argument("--k-min", type=int, default=2)
    parser.add_argument("--k-max", type=int, default=12)
    parser.add_argument("--selection-steps", type=int, default=1000)
    parser.add_argument("--selection-samples", type=int, default=4)
    parser.add_argument("--selection-splits", type=int, default=5)
    parser.add_argument("--final-steps", type=int, default=10000)
    parser.add_argument("--final-samples", type=int, default=32)
    parser.add_argument("--spatial-smoothing", type=float, default=100.0)
    parser.add_argument("--slide-indices", default="")
    parser.add_argument("--worker-id", default="main")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument(
        "--skip-queue-metadata-write",
        action="store_true",
        help=(
            "Reuse queue-level panel/runtime files created by a serial preflight. "
            "This avoids concurrent metadata writes when multiple cache shards share one queue."
        ),
    )
    parser.add_argument("--seed", type=int, default=260729)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    slides_dir = args.output_dir / "slides"
    slides_dir.mkdir(exist_ok=True)
    upstream_audit = _audit_upstream_install()
    if not args.skip_queue_metadata_write:
        _atomic_json(args.output_dir / "upstream_runtime_audit.json", upstream_audit)

    fit = sparse.load_npz(args.prepared_cache / "fit_core.npz").tocsr()
    tune = sparse.load_npz(args.prepared_cache / "tune_core.npz").tocsr()
    audit = sparse.load_npz(args.prepared_cache / "audit_core.npz").tocsr()
    raw = (fit + tune + audit).tocsr()
    if not np.allclose(raw.data, np.rint(raw.data)):
        raise RuntimeError("raw molecule add-back is not integer")
    metadata = np.load(args.prepared_cache / "metadata.npz", allow_pickle=False)
    slide_index = metadata["slide_index"].astype(np.int32)
    sample_ids = metadata["sample_ids"].astype(str)
    sources = metadata["sources"].astype(str)
    coordinates = np.load(
        args.prepared_cache / "coordinates.npy",
        allow_pickle=False,
    )
    all_genes = np.load(
        args.prepared_cache / "discovery_gene_order.npy",
        allow_pickle=False,
    ).astype(str)
    spot_index = pd.read_csv(args.spot_index)
    if len(spot_index) != raw.shape[0]:
        raise RuntimeError("spot index and raw cache row counts differ")
    if not np.array_equal(
        spot_index["slide_index"].to_numpy(dtype=np.int32),
        slide_index,
    ):
        raise RuntimeError("spot index and raw cache slide order differ")

    selected_path = args.output_dir / "selected_gene_indices.npy"
    gene_order_path = args.output_dir / "gene_order.npy"
    if selected_path.exists() and gene_order_path.exists():
        selected = np.load(selected_path, allow_pickle=False).astype(np.int64)
        genes = np.load(gene_order_path, allow_pickle=False).astype(str)
        if not np.array_equal(genes, all_genes[selected]):
            raise RuntimeError("saved common gene panel is inconsistent")
    else:
        selected = _select_shared_genes(
            raw,
            slide_index,
            sources,
            maximum_genes=int(args.maximum_genes),
        )
        genes = all_genes[selected]
        np.save(selected_path, selected)
        np.save(gene_order_path, genes)
    raw = raw[:, selected].tocsr()

    manifest = {
        "method": "unmodified_upstream_bayestme_1.0.0_per_slide_v1",
        "slides": int(len(sample_ids)),
        "spots": int(raw.shape[0]),
        "shared_count_derived_genes": int(raw.shape[1]),
        "input": "fit_core+tune_core+audit_core exact raw UMI add-back",
        "per_slide_k_selection": {
            "implementation": (
                "upstream run_phenotype_selection_single_job"
            ),
            "rule": "upstream mean heldout likelihood maximum",
            "k_min": int(args.k_min),
            "k_max": int(args.k_max),
            "n_fold": 1,
            "n_splits": int(args.selection_splits),
            "n_svi_steps": int(args.selection_steps),
            "posterior_samples": int(args.selection_samples),
            "spatial_smoothing_parameter": float(args.spatial_smoothing),
        },
        "final_fit": {
            "n_svi_steps": int(args.final_steps),
            "posterior_samples": int(args.final_samples),
        },
        "local_component_postprocessing": "none",
        "old_basis_or_theta_read": 0,
        "images_read": 0,
        "partitions_read": 0,
        "upstream_source_or_guide_modified": False,
        "runtime_audit": upstream_audit,
    }
    if not args.skip_queue_metadata_write:
        _atomic_json(args.output_dir / "queue_manifest.json", manifest)
    else:
        for required_queue_file in (
            args.output_dir / "upstream_runtime_audit.json",
            args.output_dir / "queue_manifest.json",
            args.output_dir / "selected_gene_indices.npy",
            args.output_dir / "gene_order.npy",
        ):
            if not required_queue_file.is_file():
                raise FileNotFoundError(
                    f"serial queue preflight file is missing: {required_queue_file}"
                )
    if args.prepare_only:
        print(json.dumps(manifest, indent=2), flush=True)
        return

    if args.slide_indices.strip():
        requested = [
            int(value)
            for value in args.slide_indices.split(",")
            if value.strip()
        ]
    else:
        requested = list(range(len(sample_ids)))

    records: list[dict[str, object]] = []
    started_queue = time.time()
    for position, slide in enumerate(requested):
        if slide < 0 or slide >= len(sample_ids):
            raise ValueError(f"invalid slide index {slide}")
        sample_id = str(sample_ids[slide])
        safe = f"{slide:03d}__{_safe_name(sample_id)}"
        final_dir = slides_dir / safe
        complete_path = final_dir / "complete.json"
        if complete_path.exists():
            records.append(json.loads(complete_path.read_text()))
            continue
        partial_dir = slides_dir / f"{safe}.partial"
        if partial_dir.exists() and not (
            partial_dir / "phenotype_selection.json"
        ).exists():
            shutil.rmtree(partial_dir)
        partial_dir.mkdir(exist_ok=True)
        rows = np.flatnonzero(slide_index == slide)
        metrics = _run_slide(
            slide=slide,
            rows=rows,
            raw=raw,
            coordinates=coordinates,
            spot_index=spot_index,
            genes=genes,
            sample_id=sample_id,
            source=str(sources[slide]),
            output_dir=partial_dir,
            k_min=int(args.k_min),
            k_max=int(args.k_max),
            selection_steps=int(args.selection_steps),
            selection_samples=int(args.selection_samples),
            selection_splits=int(args.selection_splits),
            final_steps=int(args.final_steps),
            final_samples=int(args.final_samples),
            spatial_smoothing=float(args.spatial_smoothing),
            seed=int(args.seed + slide * 100_003),
        )
        os.replace(partial_dir, final_dir)
        records.append(metrics)
        print(
            json.dumps(
                {
                    "worker_id": args.worker_id,
                    "completed_in_worker": position + 1,
                    "requested_in_worker": len(requested),
                    "sample_id": sample_id,
                    "selected_k": metrics["selected_k"],
                    "slide_runtime_seconds": metrics["runtime_seconds"],
                    "queue_elapsed_seconds": time.time() - started_queue,
                }
            ),
            flush=True,
        )

    worker_summary = pd.DataFrame(records)
    if not worker_summary.empty:
        worker_summary.sort_values("slide_index").to_csv(
            args.output_dir / f"worker_{_safe_name(args.worker_id)}_summary.csv",
            index=False,
        )


if __name__ == "__main__":
    main()
