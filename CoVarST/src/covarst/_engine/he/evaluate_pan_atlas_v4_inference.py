"""Evaluate standalone v4 HE inference against held-out local BayesTME inverses."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from train_virchow2_pan_atlas_consensus_mapper import (
    EPS,
    _balanced_pcc_summary,
    _correlation_vector,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inference-dir", type=Path, required=True)
    parser.add_argument("--consensus-dir", type=Path, required=True)
    parser.add_argument("--per-slide-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--block-size", type=int, default=512)
    parser.add_argument("--support-mass", type=float, default=0.90)
    return parser.parse_args()


def _support_target(
    rna: np.ndarray, composition: np.ndarray, *, mass: float
) -> np.ndarray:
    value = 0.5 * (rna + composition)
    value /= np.maximum(value.sum(axis=1, keepdims=True), EPS)
    order = np.argsort(-value, axis=1, kind="stable")
    sorted_value = np.take_along_axis(value, order, axis=1)
    previous = np.cumsum(sorted_value, axis=1) - sorted_value
    keep = previous < float(mass)
    result = np.zeros_like(value, dtype=bool)
    np.put_along_axis(result, order, keep, axis=1)
    return result


def _balanced_scalar(records: list[tuple[str, str, float]]) -> float:
    source_values = []
    for source in sorted({value[0] for value in records}):
        patient_values = []
        for patient in sorted({value[1] for value in records if value[0] == source}):
            patient_values.append(
                float(
                    np.mean(
                        [
                            value[2]
                            for value in records
                            if value[0] == source and value[1] == patient
                        ]
                    )
                )
            )
        source_values.append(float(np.mean(patient_values)))
    return float(np.mean(source_values))


def _patient_balanced_vectors_by_source(
    records: list[tuple[str, str, np.ndarray]],
) -> dict[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    for source in sorted({value[0] for value in records}):
        patient_vectors = []
        for patient in sorted({value[1] for value in records if value[0] == source}):
            patient_vectors.append(
                np.nanmean(
                    np.stack(
                        [
                            vector
                            for current_source, current_patient, vector in records
                            if current_source == source and current_patient == patient
                        ],
                        axis=0,
                    ),
                    axis=0,
                )
            )
        result[source] = np.nanmean(np.stack(patient_vectors, axis=0), axis=0)
    return result


def main() -> None:
    args = parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(f"refusing non-empty output {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    inference_manifest = json.loads(
        (args.inference_dir / "inference_manifest.json").read_text(encoding="utf-8")
    )
    if bool(inference_manifest.get("transcriptomic_inputs_loaded", True)):
        raise RuntimeError("inference manifest does not certify HE-only inference")
    spots = pd.read_csv(args.inference_dir / "inference_spot_index.csv")
    predicted_rna = np.load(
        args.inference_dir / "theta_rna_predicted.npy", mmap_mode="r"
    )
    predicted_composition = np.load(
        args.inference_dir / "theta_composition_predicted.npy", mmap_mode="r"
    )
    predicted_support = np.load(
        args.inference_dir / "program_support_probability.npy", mmap_mode="r"
    )
    if not (
        len(spots)
        == len(predicted_rna)
        == len(predicted_composition)
        == len(predicted_support)
    ):
        raise RuntimeError("inference outputs have inconsistent spot counts")

    reference = np.load(
        args.consensus_dir / "consensus_reference_probability_10k.npy",
        mmap_mode="r",
    )
    batch_reference_path = (
        args.consensus_dir / "consensus_reference_probability_by_batch_10k.npy"
    )
    batch_metadata_path = args.consensus_dir / "batch_adapter_metadata.csv"
    batch_references = None
    batch_lookup: dict[str, int] = {}
    if batch_reference_path.is_file() and batch_metadata_path.is_file():
        batch_references = np.load(batch_reference_path, mmap_mode="r")
        batch_metadata = pd.read_csv(batch_metadata_path)
        batch_lookup = {
            str(row.batch): int(row.batch_index)
            for row in batch_metadata.itertuples(index=False)
        }
    genes = np.load(
        args.consensus_dir / "gene_order_10k.npy", allow_pickle=False
    ).astype(str)
    theta_rna_all = np.load(
        args.consensus_dir / "theta_rna_consensus.npy", mmap_mode="r"
    )
    theta_composition_all = np.load(
        args.consensus_dir / "theta_composition_consensus.npy", mmap_mode="r"
    )
    assignments = pd.read_csv(args.consensus_dir / "local_program_assignment.csv")
    slide_directory = (
        assignments[["sample_id", "slide_dir"]]
        .drop_duplicates()
        .set_index("sample_id")["slide_dir"]
        .astype(str)
        .to_dict()
    )

    predicted_gene_records: list[tuple[str, str, np.ndarray]] = []
    oracle_gene_records: list[tuple[str, str, np.ndarray]] = []
    rna_program_records: list[tuple[str, str, np.ndarray]] = []
    composition_program_records: list[tuple[str, str, np.ndarray]] = []
    support_recall_records: list[tuple[str, str, float]] = []
    support_targets: list[np.ndarray] = []
    support_predictions: list[np.ndarray] = []
    slide_records: list[dict[str, object]] = []

    for slide_key, block in spots.groupby("inference_slide_key", sort=False):
        block = block.sort_values("prediction_row", kind="stable")
        prediction_rows = block["prediction_row"].to_numpy(dtype=np.int64)
        global_rows = block["global_index"].to_numpy(dtype=np.int64)
        local_rows = block["local_index"].to_numpy(dtype=np.int64)
        source = str(block["source"].iloc[0])
        technology = str(block["technology"].iloc[0])
        patient = str(block["patient"].iloc[0])
        sample_id = str(block["sample_id"].iloc[0])
        slide_reference = reference
        batch_key = f"{source}|{technology}"
        if batch_references is not None:
            if batch_key not in batch_lookup:
                raise KeyError(f"no batch-adapted consensus reference for {batch_key}")
            slide_reference = batch_references[batch_lookup[batch_key]]
        raw_slide_dir = Path(slide_directory[sample_id])
        slide_dir = raw_slide_dir
        if not slide_dir.is_dir():
            slide_dir = args.per_slide_dir / "slides" / raw_slide_dir.name
        if not slide_dir.is_dir():
            raise FileNotFoundError(raw_slide_dir)

        local_theta = np.load(
            slide_dir / "theta_rna_mean.npy", allow_pickle=False
        ).astype(np.float32)[local_rows]
        local_basis = np.load(
            slide_dir / "basis_probability_mean.npy", mmap_mode="r"
        )
        predicted_theta = np.asarray(predicted_rna[prediction_rows], dtype=np.float32)
        predicted_comp = np.asarray(
            predicted_composition[prediction_rows], dtype=np.float32
        )
        target_theta = np.asarray(theta_rna_all[global_rows], dtype=np.float32)
        target_theta /= np.maximum(target_theta.sum(axis=1, keepdims=True), EPS)
        target_comp = np.asarray(
            theta_composition_all[global_rows], dtype=np.float32
        )
        target_comp /= np.maximum(target_comp.sum(axis=1, keepdims=True), EPS)

        predicted_gene_pcc = np.full(len(genes), np.nan, dtype=np.float32)
        oracle_gene_pcc = np.full(len(genes), np.nan, dtype=np.float32)
        for start in range(0, len(genes), int(args.block_size)):
            stop = min(start + int(args.block_size), len(genes))
            target_log = np.log1p(
                1.0e4 * (local_theta @ np.asarray(local_basis[:, start:stop]))
            )
            predicted_log = np.log1p(
                1.0e4
                * (predicted_theta @ np.asarray(slide_reference[:, start:stop]))
            )
            oracle_log = np.log1p(
                1.0e4 * (target_theta @ np.asarray(slide_reference[:, start:stop]))
            )
            predicted_gene_pcc[start:stop] = _correlation_vector(
                target_log, predicted_log
            )
            oracle_gene_pcc[start:stop] = _correlation_vector(target_log, oracle_log)

        rna_program_pcc = _correlation_vector(
            np.sqrt(np.maximum(target_theta, EPS)),
            np.sqrt(np.maximum(predicted_theta, EPS)),
        )
        composition_program_pcc = _correlation_vector(
            np.sqrt(np.maximum(target_comp, EPS)),
            np.sqrt(np.maximum(predicted_comp, EPS)),
        )
        target_support = _support_target(
            target_theta, target_comp, mass=float(args.support_mass)
        )
        support_probability = np.asarray(
            predicted_support[prediction_rows], dtype=np.float32
        )
        support_count = target_support.sum(axis=1)
        predicted_order = np.argsort(-support_probability, axis=1, kind="stable")
        rank = np.argsort(predicted_order, axis=1, kind="stable")
        predicted_top = rank < support_count[:, None]
        spot_recall = (
            (predicted_top & target_support).sum(axis=1)
            / np.maximum(support_count, 1)
        )
        support_recall = float(np.mean(spot_recall))

        predicted_gene_records.append((source, patient, predicted_gene_pcc))
        oracle_gene_records.append((source, patient, oracle_gene_pcc))
        rna_program_records.append((source, patient, rna_program_pcc))
        composition_program_records.append(
            (source, patient, composition_program_pcc)
        )
        support_recall_records.append((source, patient, support_recall))
        support_targets.append(target_support.astype(np.uint8))
        support_predictions.append(support_probability)
        slide_records.append(
            {
                "slide_key": str(slide_key),
                "source": source,
                "patient": patient,
                "spots": int(len(block)),
                "gene_pcc_median": float(np.nanmedian(predicted_gene_pcc)),
                "oracle_gene_pcc_median": float(np.nanmedian(oracle_gene_pcc)),
                "rna_program_pcc_median": float(np.nanmedian(rna_program_pcc)),
                "composition_program_pcc_median": float(
                    np.nanmedian(composition_program_pcc)
                ),
                "top_program_recall": support_recall,
            }
        )

    predicted_median, predicted_by_source, predicted_vector = _balanced_pcc_summary(
        predicted_gene_records
    )
    oracle_median, oracle_by_source, oracle_vector = _balanced_pcc_summary(
        oracle_gene_records
    )
    rna_program_median, rna_program_by_source, rna_program_vector = (
        _balanced_pcc_summary(rna_program_records)
    )
    composition_program_median, composition_program_by_source, composition_program_vector = (
        _balanced_pcc_summary(composition_program_records)
    )
    predicted_gene_vectors_by_source = _patient_balanced_vectors_by_source(
        predicted_gene_records
    )
    oracle_gene_vectors_by_source = _patient_balanced_vectors_by_source(
        oracle_gene_records
    )

    macro_auprc = float("nan")
    try:
        from sklearn.metrics import average_precision_score

        target_all = np.concatenate(support_targets, axis=0)
        prediction_all = np.concatenate(support_predictions, axis=0)
        valid = (target_all.sum(axis=0) > 0) & (
            target_all.sum(axis=0) < len(target_all)
        )
        if np.any(valid):
            macro_auprc = float(
                np.nanmean(
                    average_precision_score(
                        target_all[:, valid], prediction_all[:, valid], average=None
                    )
                )
            )
    except ImportError:
        pass

    pd.DataFrame(slide_records).to_csv(
        args.output_dir / "metrics_by_slide.csv", index=False
    )
    gene_table = {
        "gene_index": np.arange(len(genes), dtype=np.int64),
        "gene": genes,
        "predicted_patient_source_balanced_pcc": predicted_vector,
        "oracle_patient_source_balanced_pcc": oracle_vector,
    }
    for source, vector in predicted_gene_vectors_by_source.items():
        gene_table[f"patient_balanced_pcc__{source}"] = vector
    for source, vector in oracle_gene_vectors_by_source.items():
        gene_table[f"oracle_patient_balanced_pcc__{source}"] = vector
    pd.DataFrame(gene_table).to_csv(
        args.output_dir / "gene_pcc_10k.csv", index=False
    )
    pd.DataFrame(
        {
            "program_index": np.arange(reference.shape[0], dtype=np.int64),
            "rna_patient_source_balanced_pcc": rna_program_vector,
            "composition_patient_source_balanced_pcc": composition_program_vector,
        }
    ).to_csv(
        args.output_dir / f"program_pcc_g{int(reference.shape[0])}.csv",
        index=False,
    )

    payload = {
        "method": f"standalone_evaluation_of_he_only_pan_atlas_g{int(reference.shape[0])}_inference",
        "inference_manifest": inference_manifest,
        "reference_for_pcc": "original_per_slide_top10k_theta_rna_times_basis_probability",
        "prediction_decoder": (
            "standalone_predicted_theta_rna_times_slide_matched_batch_adapted_pan_atlas_W"
            if batch_references is not None
            else "standalone_predicted_theta_rna_times_batch_neutral_pan_atlas_W"
        ),
        "patient_source_balanced_gene_pcc_median": predicted_median,
        "by_source_patient_balanced_gene_pcc_median": predicted_by_source,
        "oracle_consensus_gene_pcc_median": oracle_median,
        "oracle_by_source_gene_pcc_median": oracle_by_source,
        "rna_program_pcc_median": rna_program_median,
        "rna_program_pcc_by_source": rna_program_by_source,
        "composition_program_pcc_median": composition_program_median,
        "composition_program_pcc_by_source": composition_program_by_source,
        "patient_source_balanced_top_program_recall": _balanced_scalar(
            support_recall_records
        ),
        "support_macro_auprc": macro_auprc,
        "slides": slide_records,
    }
    (args.output_dir / "metrics.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
