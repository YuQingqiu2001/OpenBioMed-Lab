"""Fit a shared ~40-program dictionary directly against per-slide BayesTME inverses.

The shared reference is learned on a locked patient-level training split.  Each
held-out slide may fit only its local-program-to-shared-program map; validation
selects the dictionary size and test is evaluated once after selection.

The reconstruction target is always the original per-slide inverse::

    theta_local @ basis_local ~= theta_local @ local_map @ W_batch

No H128/G128 reconstruction term is used.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from fit_bayestme_top10k_rccd_consensus import (
    EPS,
    _fit_mapping_only,
    _fit_model,
    _initial_reference,
    _mapping_logits,
    _model_arrays,
    _save_final,
    _sha256,
    _summary,
    evaluate_rows,
    load_bundle,
)
from train_virchow2_pan_atlas_consensus_mapper import (
    _balanced_pcc_summary,
    _correlation_vector,
)


METHOD = "bayestme_top10k_direct_local_inverse_locked_split_consensus_v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue-dir", type=Path, required=True)
    parser.add_argument("--master-spot-index", type=Path, required=True)
    parser.add_argument("--split-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--candidates", type=int, nargs="+", default=[36, 40, 44, 48])
    parser.add_argument("--target-programs", type=int, default=40)
    parser.add_argument("--selection-tolerance", type=float, default=0.002)
    parser.add_argument("--candidate-epochs", type=int, default=180)
    parser.add_argument("--mapping-epochs", type=int, default=140)
    parser.add_argument("--final-epochs", type=int, default=260)
    parser.add_argument("--warmup-epochs", type=int, default=40)
    parser.add_argument("--learning-rate", type=float, default=0.025)
    parser.add_argument("--technology-bound", type=float, default=0.35)
    parser.add_argument("--batch-bound", type=float, default=0.20)
    parser.add_argument("--lambda-basis", type=float, default=0.20)
    parser.add_argument("--lambda-log-basis", type=float, default=0.25)
    parser.add_argument("--lambda-entropy", type=float, default=0.010)
    parser.add_argument("--lambda-collision", type=float, default=0.05)
    parser.add_argument("--lambda-diversity", type=float, default=0.02)
    parser.add_argument("--lambda-technology", type=float, default=0.05)
    parser.add_argument("--lambda-batch", type=float, default=0.10)
    parser.add_argument("--gate-balanced-gene-pcc", type=float, default=0.75)
    parser.add_argument("--gate-relative-median", type=float, default=0.20)
    parser.add_argument("--gate-relative-q95", type=float, default=0.35)
    parser.add_argument("--gate-spot-cosine", type=float, default=0.97)
    parser.add_argument("--gate-gene-pcc", type=float, default=0.75)
    parser.add_argument("--evaluation-block", type=int, default=512)
    parser.add_argument("--seed", type=int, default=260803)
    parser.add_argument(
        "--posterior-stability", action=argparse.BooleanOptionalAction, default=True
    )
    parser.add_argument("--max-slides", type=int)
    # Compatibility attributes consumed by the shared artifact writer.
    parser.set_defaults(holdout_fraction=0.0)
    return parser.parse_args()


def _locked_rows(bundle, manifest: pd.DataFrame) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    required = {"sample_id", "split"}
    if not required.issubset(manifest.columns):
        raise ValueError(f"split manifest missing {sorted(required - set(manifest.columns))}")
    split_by_sample = (
        manifest[["sample_id", "split"]]
        .drop_duplicates("sample_id")
        .assign(sample_id=lambda x: x["sample_id"].astype(str))
        .set_index("sample_id")["split"]
        .astype(str)
        .to_dict()
    )
    labels = bundle.slides["sample_id"].astype(str).map(split_by_sample).fillna("unassigned")
    allowed = {"train", "val", "test", "unassigned"}
    unexpected = sorted(set(labels) - allowed)
    if unexpected:
        raise ValueError(f"unexpected split labels: {unexpected}")
    rows = {
        name: np.flatnonzero(labels.to_numpy() == name).astype(np.int64)
        for name in ["train", "val", "test", "unassigned"]
    }
    for name in ["train", "val", "test"]:
        if not len(rows[name]):
            raise RuntimeError(f"locked split {name!r} contains no slides")
    split = bundle.slides.copy()
    split["split"] = labels.to_numpy()
    return rows, split


def _direct_gene_pcc(
    bundle,
    slide_rows: np.ndarray,
    mappings: np.ndarray,
    reference_by_batch: np.ndarray,
    *,
    block: int,
) -> tuple[dict[str, object], np.ndarray]:
    records: list[tuple[str, str, np.ndarray]] = []
    slide_medians: list[float] = []
    spot_count = 0
    for local_position, slide_row in enumerate(slide_rows):
        k = int(bundle.valid_mask[slide_row].sum())
        theta = np.asarray(bundle.theta[slide_row], dtype=np.float32)
        basis = np.asarray(bundle.basis[slide_row, :k], dtype=np.float32)
        mapping = np.asarray(mappings[local_position, :k], dtype=np.float32)
        reference = np.asarray(
            reference_by_batch[int(bundle.batch_index[slide_row])], dtype=np.float32
        )
        predicted_basis = mapping @ reference
        gene_pcc = np.full(bundle.n_genes, np.nan, dtype=np.float32)
        for start in range(0, bundle.n_genes, int(block)):
            stop = min(start + int(block), bundle.n_genes)
            target = np.log1p(1.0e4 * (theta @ basis[:, start:stop]))
            predicted = np.log1p(1.0e4 * (theta @ predicted_basis[:, start:stop]))
            gene_pcc[start:stop] = _correlation_vector(target, predicted)
        row = bundle.slides.iloc[int(slide_row)]
        records.append((str(row["source"]), str(row["patient"]), gene_pcc))
        slide_medians.append(float(np.nanmedian(gene_pcc)))
        spot_count += len(theta)
    median, by_source, vector = _balanced_pcc_summary(records)
    return (
        {
            "slides": int(len(slide_rows)),
            "patients": int(bundle.slides.iloc[slide_rows]["patient"].nunique()),
            "spots": int(spot_count),
            "genes": int(bundle.n_genes),
            "patient_source_balanced_gene_pcc_median": float(median),
            "by_source_patient_balanced_gene_pcc_median": by_source,
            "per_slide_gene_pcc_median": float(np.median(slide_medians)),
        },
        vector,
    )


def _fit_heldout_mapping(bundle, rows, neutral, by_batch, args, seed, device):
    if not len(rows):
        return np.zeros((0, bundle.k_max, len(neutral)), dtype=np.float32)
    return _fit_mapping_only(
        bundle,
        rows,
        neutral,
        by_batch,
        device=device,
        epochs=int(args.mapping_epochs),
        learning_rate=float(args.learning_rate),
        seed=int(seed),
        args=args,
    )


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"refusing existing output directory: {args.output_dir}")
    random.seed(int(args.seed))
    np.random.seed(int(args.seed))
    torch.manual_seed(int(args.seed))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    bundle = load_bundle(
        args.queue_dir,
        args.master_spot_index,
        use_posterior_stability=bool(args.posterior_stability),
        max_slides=args.max_slides,
    )
    manifest = pd.read_csv(args.split_manifest)
    rows, split = _locked_rows(bundle, manifest)
    print(
        json.dumps(
            {
                "device": str(device),
                "slides": bundle.n_slides,
                "local_programs": int(bundle.valid_mask.sum()),
                "genes": bundle.n_genes,
                "locked_split_slides": {name: int(len(value)) for name, value in rows.items()},
                "objective": "theta_local @ basis_local ~= theta_local @ map_local @ W_batch",
                "h128_reconstruction_target_used": False,
            },
            ensure_ascii=False,
            indent=2,
        ),
        flush=True,
    )

    candidate_records: list[dict[str, object]] = []
    fitted: dict[int, dict[str, object]] = {}
    all_rows = np.arange(bundle.n_slides, dtype=np.int64)
    train_local_programs = int(bundle.valid_mask[rows["train"]].sum())
    for programs in sorted(set(map(int, args.candidates))):
        print(json.dumps({"candidate": programs, "stage": "train_reference"}), flush=True)
        reference_init, mapping_init = _initial_reference(
            bundle,
            rows["train"],
            programs,
            int(args.seed) + programs,
            fallback_slide_rows=all_rows,
        )
        model, history = _fit_model(
            bundle,
            rows["train"],
            reference_init,
            mapping_init,
            device=device,
            epochs=int(args.candidate_epochs),
            warmup_epochs=min(int(args.warmup_epochs), int(args.candidate_epochs) // 2),
            learning_rate=float(args.learning_rate),
            seed=int(args.seed) + programs,
            args=args,
        )
        neutral, by_batch, technology_delta, batch_delta, train_mapping = _model_arrays(model)
        validation_mapping = _fit_heldout_mapping(
            bundle,
            rows["val"],
            neutral,
            by_batch,
            args,
            int(args.seed) + 10_000 + programs,
            device,
        )
        validation_slide_metrics = evaluate_rows(
            bundle,
            rows["val"],
            validation_mapping,
            by_batch,
            block=int(args.evaluation_block),
        )
        validation_summary = _summary(validation_slide_metrics)
        balanced, vector = _direct_gene_pcc(
            bundle,
            rows["val"],
            validation_mapping,
            by_batch,
            block=int(args.evaluation_block),
        )
        record = {
            "consensus_programs": int(programs),
            "reference_initialization_scope": (
                "train"
                if int(programs) <= train_local_programs
                else "all_slides_initialization_only"
            ),
            "validation_balanced_gene_pcc_median": balanced[
                "patient_source_balanced_gene_pcc_median"
            ],
            "validation_per_slide_gene_pcc_median": balanced["per_slide_gene_pcc_median"],
            "validation_spot_cosine_median": validation_summary[
                "median_spot_logcp10k_cosine"
            ]["median"],
            "validation_relative_error_median": validation_summary[
                "relative_probability_error"
            ]["median"],
            "validation_relative_error_q95": validation_summary[
                "relative_probability_error"
            ]["q95"],
        }
        for source, value in balanced["by_source_patient_balanced_gene_pcc_median"].items():
            record[f"validation_pcc_{source}"] = value
        candidate_records.append(record)
        fitted[programs] = {
            "neutral": neutral,
            "by_batch": by_batch,
            "technology_delta": technology_delta,
            "batch_delta": batch_delta,
            "train_mapping": train_mapping,
            "validation_mapping": validation_mapping,
            "history": history,
            "validation_balanced": balanced,
            "validation_gene_vector": vector,
        }
        print(json.dumps({"candidate_result": record}, ensure_ascii=False), flush=True)
        del model
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    candidate_table = pd.DataFrame(candidate_records).sort_values("consensus_programs")
    best = float(candidate_table["validation_balanced_gene_pcc_median"].max())
    eligible = candidate_table.loc[
        candidate_table["validation_balanced_gene_pcc_median"]
        >= best - float(args.selection_tolerance)
    ].copy()
    eligible["distance_to_target"] = (
        eligible["consensus_programs"] - int(args.target_programs)
    ).abs()
    selected_row = eligible.sort_values(
        ["distance_to_target", "consensus_programs"], ascending=[True, True]
    ).iloc[0]
    selected_programs = int(selected_row["consensus_programs"])
    selected = fitted[selected_programs]
    print(
        json.dumps(
            {
                "selected_programs": selected_programs,
                "selection_split": "val",
                "best_validation_balanced_gene_pcc": best,
                "selection_tolerance": float(args.selection_tolerance),
                "test_used_for_selection": False,
            },
            ensure_ascii=False,
        ),
        flush=True,
    )

    # Test mapping and test recovery are computed only after K and W are frozen.
    test_mapping = _fit_heldout_mapping(
        bundle,
        rows["test"],
        selected["neutral"],
        selected["by_batch"],
        args,
        int(args.seed) + 20_000 + selected_programs,
        device,
    )
    unassigned_mapping = _fit_heldout_mapping(
        bundle,
        rows["unassigned"],
        selected["neutral"],
        selected["by_batch"],
        args,
        int(args.seed) + 30_000 + selected_programs,
        device,
    )
    all_mapping = np.zeros(
        (bundle.n_slides, bundle.k_max, selected_programs), dtype=np.float32
    )
    all_mapping[rows["train"]] = selected["train_mapping"]
    all_mapping[rows["val"]] = selected["validation_mapping"]
    all_mapping[rows["test"]] = test_mapping
    if len(rows["unassigned"]):
        all_mapping[rows["unassigned"]] = unassigned_mapping

    recovery: dict[str, dict[str, object]] = {}
    gene_vectors: dict[str, np.ndarray] = {}
    for split_name, mapping in [
        ("train", selected["train_mapping"]),
        ("val", selected["validation_mapping"]),
        ("test", test_mapping),
    ]:
        recovery[split_name], gene_vectors[split_name] = _direct_gene_pcc(
            bundle,
            rows[split_name],
            mapping,
            selected["by_batch"],
            block=int(args.evaluation_block),
        )

    # After K is locked and the untouched test audit has been recorded, refit
    # the production dictionary on all 177 slides as explicitly requested.
    full_mapping_init = _mapping_logits(
        bundle.basis, bundle.valid_mask, selected["neutral"]
    )
    print(
        json.dumps(
            {
                "stage": "production_refit_all_slides",
                "slides": int(bundle.n_slides),
                "programs": int(selected_programs),
            }
        ),
        flush=True,
    )
    final_model, final_history = _fit_model(
        bundle,
        all_rows,
        selected["neutral"],
        full_mapping_init,
        device=device,
        epochs=int(args.final_epochs),
        warmup_epochs=min(int(args.warmup_epochs), int(args.final_epochs) // 2),
        learning_rate=float(args.learning_rate),
        seed=int(args.seed) + 40_000 + selected_programs,
        args=args,
    )
    (
        final_neutral,
        final_by_batch,
        final_technology_delta,
        final_batch_delta,
        final_mapping,
    ) = _model_arrays(final_model)
    full_recovery, full_gene_vector = _direct_gene_pcc(
        bundle,
        all_rows,
        final_mapping,
        final_by_batch,
        block=int(args.evaluation_block),
    )
    all_metrics = evaluate_rows(
        bundle,
        all_rows,
        final_mapping,
        final_by_batch,
        block=int(args.evaluation_block),
    )
    _save_final(
        bundle,
        args.output_dir,
        final_neutral,
        final_by_batch,
        final_technology_delta,
        final_batch_delta,
        final_mapping,
        all_metrics,
        final_history,
        candidate_table,
        split,
        selected_programs,
        "validation_balanced_gene_pcc_then_nearest_to_target_within_tolerance",
        args,
    )
    for split_name, vector in gene_vectors.items():
        pd.DataFrame(
            {
                "gene_index": np.arange(bundle.n_genes),
                "gene": bundle.genes,
                "oracle_pcc": vector,
            }
        ).to_csv(args.output_dir / f"gene_pcc_10k_{split_name}.csv", index=False)
    pd.DataFrame(
        {
            "gene_index": np.arange(bundle.n_genes),
            "gene": bundle.genes,
            "oracle_pcc": full_gene_vector,
        }
    ).to_csv(args.output_dir / "gene_pcc_10k_all177_full_fit.csv", index=False)

    reference_path = args.output_dir / "consensus_reference_probability_10k.npy"
    theta_path = args.output_dir / "theta_rna_consensus.npy"
    integrity = {
        "reference_row_sum_maximum_error": float(
            np.max(np.abs(np.asarray(final_neutral).sum(axis=1) - 1.0))
        ),
        "reference_nonnegative": bool(np.all(np.asarray(final_neutral) >= 0)),
        "mapping_row_sum_maximum_error": float(
            np.max(np.abs(final_mapping[bundle.valid_mask].sum(axis=1) - 1.0))
        ),
    }
    contract = {
        "method": METHOD,
        "objective": "theta_local @ W_local ~= theta_local @ M_local @ W_shared_batch",
        "gold_standard": "pre-integration per-slide BayesTME inverse expression",
        "h128_reconstruction_target_used": False,
        "candidate_reference_fit_split": "train_only",
        "candidate_selection_split": "val_only",
        "test_used_for_candidate_selection": False,
        "heldout_mapping_rule": "fit M_local only; W_shared and batch adapters frozen",
        "production_reference_fit": "all_177_slides_after_K_was_locked",
        "production_reference_fit_slides": int(bundle.n_slides),
        "batch_calibration": "technology plus centered source-by-technology logit adapters; train-only for heldout audit and all-177 for production refit",
        "consensus_programs": selected_programs,
        "target_programs": int(args.target_programs),
        "genes": int(bundle.n_genes),
        "candidate_selection_tolerance": float(args.selection_tolerance),
        "initialization_fallback_policy": (
            "When candidate K exceeds the number of local programs in the training "
            "partition, all-tissue local programs are used only for KMeans initialization; "
            "candidate parameter optimization remains train-only."
        ),
        "candidates": candidate_records,
        "selected_candidate": selected_row.to_dict(),
        "heldout_recovery_before_all177_refit": recovery,
        "production_all177_fit_recovery_against_original_local_inverse": full_recovery,
        "gate_balanced_gene_pcc": float(args.gate_balanced_gene_pcc),
        "production_ready": bool(
            full_recovery["patient_source_balanced_gene_pcc_median"]
            >= float(args.gate_balanced_gene_pcc)
            and integrity["reference_nonnegative"]
            and integrity["reference_row_sum_maximum_error"] <= 1.0e-5
            and integrity["mapping_row_sum_maximum_error"] <= 1.0e-5
        ),
        "integrity": integrity,
        "reference_sha256": _sha256(reference_path),
        "theta_sha256": _sha256(theta_path),
    }
    (args.output_dir / "consensus_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(contract, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
