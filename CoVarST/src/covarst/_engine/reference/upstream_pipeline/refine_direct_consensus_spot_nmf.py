"""Spot-level refinement of a direct local-inverse consensus dictionary.

This stage keeps the selected program count fixed, but removes the restrictive
theta_shared = theta_local @ M_local constraint.  All slides jointly refine a
free nonnegative spot-by-program matrix and the shared/batch-adapted reference
against each slide's original BayesTME inverse expression.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from fit_bayestme_top10k_rccd_consensus import EPS, RCCD, _sha256, load_bundle
from train_virchow2_pan_atlas_consensus_mapper import (
    _balanced_pcc_summary,
    _correlation_vector,
)


METHOD = "bayestme_top10k_direct_local_inverse_spot_nmf_all177_v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue-dir", type=Path, required=True)
    parser.add_argument("--master-spot-index", type=Path, required=True)
    parser.add_argument("--split-manifest", type=Path, required=True)
    parser.add_argument("--input-consensus-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--steps", type=int, default=3200)
    parser.add_argument("--spots-per-step", type=int, default=192)
    parser.add_argument("--genes-per-step", type=int, default=2048)
    parser.add_argument("--reference-learning-rate", type=float, default=0.006)
    parser.add_argument("--theta-learning-rate", type=float, default=0.025)
    parser.add_argument("--gene-pcc-weight", type=float, default=0.75)
    parser.add_argument("--spot-cosine-weight", type=float, default=0.20)
    parser.add_argument("--theta-anchor-weight", type=float, default=0.01)
    parser.add_argument("--adapter-weight", type=float, default=0.03)
    parser.add_argument("--technology-bound", type=float, default=0.45)
    parser.add_argument("--batch-bound", type=float, default=0.30)
    parser.add_argument("--evaluation-block", type=int, default=512)
    parser.add_argument("--gate-balanced-gene-pcc", type=float, default=0.75)
    parser.add_argument("--seed", type=int, default=260804)
    return parser.parse_args()


def _initialise_rccd(bundle, input_dir: Path, args, device):
    reference = np.load(
        input_dir / "consensus_reference_probability_10k.npy", allow_pickle=False
    ).astype(np.float32)
    dummy_mapping = np.zeros((1, 1, len(reference)), dtype=np.float32)
    model = RCCD(
        reference,
        dummy_mapping,
        len(bundle.technology_names),
        len(bundle.batch_names),
        bundle.batch_to_technology,
        bundle.technology_prior,
        bundle.batch_prior_within_technology,
        technology_bound=float(args.technology_bound),
        batch_bound=float(args.batch_bound),
    ).to(device)
    model.a_logits.requires_grad_(False)
    technology = np.load(
        input_dir / "technology_gene_adapters.npy", allow_pickle=False
    ).astype(np.float32)
    batch = np.load(input_dir / "batch_gene_adapters.npy", allow_pickle=False).astype(
        np.float32
    )
    with torch.no_grad():
        model.technology_raw.copy_(
            torch.from_numpy(
                np.arctanh(
                    np.clip(
                        technology / max(float(args.technology_bound), EPS), -0.95, 0.95
                    )
                )
            ).to(device)
        )
        model.batch_raw.copy_(
            torch.from_numpy(
                np.arctanh(
                    np.clip(batch / max(float(args.batch_bound), EPS), -0.95, 0.95)
                )
            ).to(device)
        )
    return model, reference


def _sample_genes(rng, gene_weight: np.ndarray, count: int) -> np.ndarray:
    count = min(int(count), len(gene_weight))
    uniform_count = count // 2
    weighted_count = count - uniform_count
    uniform = rng.choice(len(gene_weight), size=uniform_count, replace=False)
    remaining = np.ones(len(gene_weight), dtype=bool)
    remaining[uniform] = False
    candidates = np.flatnonzero(remaining)
    probability = np.maximum(gene_weight[candidates].astype(np.float64), EPS)
    probability /= probability.sum()
    weighted = rng.choice(
        candidates, size=weighted_count, replace=False, p=probability
    )
    return np.sort(np.concatenate([uniform, weighted])).astype(np.int64)


def _fit(bundle, input_dir: Path, args, device):
    rng = np.random.default_rng(int(args.seed))
    model, initial_reference = _initialise_rccd(bundle, input_dir, args, device)
    theta_initial = np.load(
        input_dir / "theta_rna_consensus.npy", allow_pickle=False
    ).astype(np.float32)
    programs = int(theta_initial.shape[1])
    theta_embedding = torch.nn.Embedding(
        len(theta_initial), programs, sparse=True, device=device
    )
    with torch.no_grad():
        theta_embedding.weight.copy_(
            torch.log(torch.from_numpy(theta_initial).to(device).clamp_min(1.0e-7))
        )
    dense_parameters = [model.w_logits, model.technology_raw, model.batch_raw]
    dense_optimizer = torch.optim.Adam(
        dense_parameters, lr=float(args.reference_learning_rate)
    )
    theta_optimizer = torch.optim.SparseAdam(
        theta_embedding.parameters(), lr=float(args.theta_learning_rate)
    )
    slide_probability = np.asarray(bundle.slide_weight, dtype=np.float64)
    slide_probability /= slide_probability.sum()
    gene_weight = np.asarray(bundle.gene_weight, dtype=np.float64)
    history: list[dict[str, float]] = []
    for step in range(1, int(args.steps) + 1):
        slide = int(rng.choice(bundle.n_slides, p=slide_probability))
        local_theta_np = bundle.theta[slide]
        take = min(int(args.spots_per_step), len(local_theta_np))
        local_rows = rng.choice(len(local_theta_np), size=take, replace=False)
        genes_np = _sample_genes(rng, gene_weight, int(args.genes_per_step))
        global_rows_np = bundle.global_rows[slide][local_rows]
        global_rows = torch.from_numpy(global_rows_np).long().to(device)
        genes = torch.from_numpy(genes_np).long().to(device)
        local_theta = torch.from_numpy(local_theta_np[local_rows]).to(device)
        local_basis = torch.from_numpy(bundle.basis[slide, : local_theta_np.shape[1]]).to(
            device
        )
        target = local_theta @ local_basis[:, genes]

        dense_optimizer.zero_grad(set_to_none=True)
        theta_optimizer.zero_grad(set_to_none=True)
        theta_logits = theta_embedding(global_rows)
        theta = torch.softmax(theta_logits, dim=1)
        technology_delta, batch_delta = model.deltas()
        batch_index = int(bundle.batch_index[slide])
        technology_index = int(bundle.batch_to_technology[batch_index])
        reference = torch.softmax(
            model.w_logits
            + technology_delta[technology_index][None, :]
            + batch_delta[batch_index][None, :],
            dim=1,
        )
        predicted = theta @ reference[:, genes]
        target_log = torch.log1p(1.0e4 * target)
        predicted_log = torch.log1p(1.0e4 * predicted)
        log_huber = torch.nn.functional.smooth_l1_loss(
            predicted_log, target_log, beta=0.25
        )
        target_center = target_log - target_log.mean(dim=0, keepdim=True)
        predicted_center = predicted_log - predicted_log.mean(dim=0, keepdim=True)
        denominator = torch.sqrt(
            torch.sum(target_center.square(), dim=0).clamp_min(1.0e-8)
            * torch.sum(predicted_center.square(), dim=0).clamp_min(1.0e-8)
        )
        gene_pcc = torch.sum(target_center * predicted_center, dim=0) / denominator
        gene_pcc_loss = 1.0 - torch.mean(gene_pcc)
        spot_cosine_loss = 1.0 - torch.nn.functional.cosine_similarity(
            target_log, predicted_log, dim=1
        ).mean()
        theta_anchor = torch.from_numpy(theta_initial[global_rows_np]).to(device)
        theta_kl = torch.mean(
            torch.sum(
                theta_anchor
                * (
                    torch.log(theta_anchor.clamp_min(1.0e-8))
                    - torch.log(theta.clamp_min(1.0e-8))
                ),
                dim=1,
            )
        )
        adapter_penalty = torch.mean(technology_delta.square()) + torch.mean(
            batch_delta.square()
        )
        total = (
            log_huber
            + float(args.gene_pcc_weight) * gene_pcc_loss
            + float(args.spot_cosine_weight) * spot_cosine_loss
            + float(args.theta_anchor_weight) * theta_kl
            + float(args.adapter_weight) * adapter_penalty
        )
        if not torch.isfinite(total):
            raise FloatingPointError(f"non-finite loss at step {step}")
        total.backward()
        torch.nn.utils.clip_grad_norm_(dense_parameters, 5.0)
        dense_optimizer.step()
        theta_optimizer.step()
        if step == 1 or step % 100 == 0 or step == int(args.steps):
            row = {
                "step": float(step),
                "total": float(total.detach().cpu()),
                "log_huber": float(log_huber.detach().cpu()),
                "gene_pcc": float((1.0 - gene_pcc_loss).detach().cpu()),
                "spot_cosine": float((1.0 - spot_cosine_loss).detach().cpu()),
                "theta_kl": float(theta_kl.detach().cpu()),
                "adapter_penalty": float(adapter_penalty.detach().cpu()),
            }
            history.append(row)
            print(json.dumps({"spot_refine": row}), flush=True)

    with torch.no_grad():
        neutral = model.neutral_reference().cpu().numpy().astype(np.float32)
        by_batch = model.reference_by_batch().cpu().numpy().astype(np.float32)
        technology_delta, batch_delta = model.deltas()
        theta = np.zeros_like(theta_initial)
        for start in range(0, len(theta), 4096):
            indices = torch.arange(start, min(start + 4096, len(theta)), device=device)
            theta[start : start + len(indices)] = (
                torch.softmax(theta_embedding(indices), dim=1).cpu().numpy()
            )
    return (
        neutral,
        by_batch,
        technology_delta.cpu().numpy().astype(np.float32),
        batch_delta.cpu().numpy().astype(np.float32),
        theta,
        history,
        initial_reference,
        theta_initial,
    )


def _evaluate(bundle, theta, by_batch, split_manifest: Path, block: int):
    manifest = pd.read_csv(split_manifest)
    split_by_sample = (
        manifest[["sample_id", "split"]]
        .drop_duplicates("sample_id")
        .assign(sample_id=lambda x: x["sample_id"].astype(str))
        .set_index("sample_id")["split"]
        .astype(str)
        .to_dict()
    )
    records_by_split: dict[str, list[tuple[str, str, np.ndarray]]] = {
        "train": [],
        "val": [],
        "test": [],
        "all177": [],
    }
    slide_rows: list[dict[str, object]] = []
    spots_by_split = {name: 0 for name in records_by_split}
    for slide in range(bundle.n_slides):
        k = int(bundle.valid_mask[slide].sum())
        local_theta = bundle.theta[slide]
        basis = bundle.basis[slide, :k]
        global_rows = bundle.global_rows[slide]
        shared_theta = theta[global_rows]
        reference = by_batch[int(bundle.batch_index[slide])]
        gene_pcc = np.full(bundle.n_genes, np.nan, dtype=np.float32)
        spot_dot = np.zeros(len(local_theta), dtype=np.float64)
        target_sq = np.zeros(len(local_theta), dtype=np.float64)
        predicted_sq = np.zeros(len(local_theta), dtype=np.float64)
        for start in range(0, bundle.n_genes, int(block)):
            stop = min(start + int(block), bundle.n_genes)
            target = np.log1p(1.0e4 * (local_theta @ basis[:, start:stop]))
            predicted = np.log1p(
                1.0e4 * (shared_theta @ reference[:, start:stop])
            )
            gene_pcc[start:stop] = _correlation_vector(target, predicted)
            spot_dot += np.sum(target * predicted, axis=1)
            target_sq += np.sum(np.square(target), axis=1)
            predicted_sq += np.sum(np.square(predicted), axis=1)
        row = bundle.slides.iloc[slide]
        sample_id = str(row["sample_id"])
        split = split_by_sample.get(sample_id, "unassigned")
        record = (str(row["source"]), str(row["patient"]), gene_pcc)
        records_by_split["all177"].append(record)
        spots_by_split["all177"] += len(local_theta)
        if split in records_by_split:
            records_by_split[split].append(record)
            spots_by_split[split] += len(local_theta)
        spot_cosine = spot_dot / np.maximum(
            np.sqrt(target_sq * predicted_sq), EPS
        )
        slide_rows.append(
            {
                **row.to_dict(),
                "split": split,
                "gene_pcc_median": float(np.nanmedian(gene_pcc)),
                "spot_logcp10k_cosine_median": float(np.median(spot_cosine)),
                "spot_logcp10k_cosine_q05": float(np.quantile(spot_cosine, 0.05)),
            }
        )
    metrics = pd.DataFrame(slide_rows)
    summaries: dict[str, dict[str, object]] = {}
    vectors: dict[str, np.ndarray] = {}
    for name, records in records_by_split.items():
        median, by_source, vector = _balanced_pcc_summary(records)
        vectors[name] = vector
        block_rows = metrics if name == "all177" else metrics.loc[metrics["split"] == name]
        summaries[name] = {
            "slides": int(len(records)),
            "spots": int(spots_by_split[name]),
            "patient_source_balanced_gene_pcc_median": float(median),
            "by_source_patient_balanced_gene_pcc_median": by_source,
            "per_slide_gene_pcc_median": float(block_rows["gene_pcc_median"].median()),
            "per_slide_spot_cosine_median": float(
                block_rows["spot_logcp10k_cosine_median"].median()
            ),
        }
    return summaries, vectors, metrics


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"refusing existing output {args.output_dir}")
    np.random.seed(int(args.seed))
    torch.manual_seed(int(args.seed))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    bundle = load_bundle(
        args.queue_dir,
        args.master_spot_index,
        use_posterior_stability=True,
        max_slides=None,
    )
    print(
        json.dumps(
            {
                "device": str(device),
                "slides_entering_fit": int(bundle.n_slides),
                "spots": int(sum(len(value) for value in bundle.theta)),
                "genes": int(bundle.n_genes),
                "input_consensus": str(args.input_consensus_dir),
            },
            ensure_ascii=False,
        ),
        flush=True,
    )
    (
        neutral,
        by_batch,
        technology_delta,
        batch_delta,
        theta,
        history,
        initial_reference,
        theta_initial,
    ) = _fit(bundle, args.input_consensus_dir, args, device)
    summaries, vectors, metrics = _evaluate(
        bundle, theta, by_batch, args.split_manifest, int(args.evaluation_block)
    )
    shutil.copytree(args.input_consensus_dir, args.output_dir)
    np.save(args.output_dir / "consensus_reference_probability_10k.npy", neutral)
    np.save(
        args.output_dir / "consensus_reference_logcp10k_10k.npy",
        np.log1p(1.0e4 * neutral).astype(np.float32),
    )
    np.save(args.output_dir / "consensus_reference_probability_by_batch_10k.npy", by_batch)
    np.save(args.output_dir / "technology_gene_adapters.npy", technology_delta)
    np.save(args.output_dir / "batch_gene_adapters.npy", batch_delta)
    np.save(args.output_dir / "theta_rna_consensus.npy", theta)
    pd.DataFrame(history).to_csv(
        args.output_dir / "spot_nmf_refinement_history.csv", index=False
    )
    metrics.to_csv(
        args.output_dir / "spot_nmf_reconstruction_metrics_by_slide.csv", index=False
    )
    genes = bundle.genes
    for name, vector in vectors.items():
        pd.DataFrame(
            {"gene_index": np.arange(len(genes)), "gene": genes, "oracle_pcc": vector}
        ).to_csv(args.output_dir / f"gene_pcc_10k_spot_nmf_{name}.csv", index=False)

    metadata_path = args.output_dir / "consensus_program_metadata.csv"
    if metadata_path.exists():
        metadata = pd.read_csv(metadata_path)
        mass = theta.astype(np.float64).sum(axis=0)
        mass /= mass.sum()
        metadata["rna_mass_fraction"] = mass
        metadata["top_genes"] = [
            ";".join(genes[np.argsort(-neutral[row], kind="stable")[:50]])
            for row in range(len(neutral))
        ]
        metadata.to_csv(metadata_path, index=False)

    reference_drift = np.mean(
        np.sum(
            (
                np.sqrt(np.maximum(initial_reference, 0.0))
                - np.sqrt(np.maximum(neutral, 0.0))
            )
            ** 2,
            axis=1,
        )
    )
    theta_drift = np.mean(
        np.sum(
            theta_initial
            * (
                np.log(np.maximum(theta_initial, EPS))
                - np.log(np.maximum(theta, EPS))
            ),
            axis=1,
        )
    )
    reference_path = args.output_dir / "consensus_reference_probability_10k.npy"
    theta_path = args.output_dir / "theta_rna_consensus.npy"
    contract = {
        "method": METHOD,
        "objective": "free theta44 and shared W44 directly reconstruct every pre-integration theta_local @ W_local",
        "input_consensus_dir": str(args.input_consensus_dir),
        "program_count_locked_before_refinement": int(neutral.shape[0]),
        "slides_entering_production_fit": int(bundle.n_slides),
        "spots_entering_production_fit": int(sum(len(value) for value in bundle.theta)),
        "genes": int(bundle.n_genes),
        "batch_calibration": "bounded technology and centered source-by-technology gene adapters",
        "spot_level_theta_is_free": True,
        "theta_composition_semantics": "retained local-composition mapping from the K44 parent; RNA theta is spot-refined",
        "recovery_against_original_local_inverse": summaries,
        "gate_balanced_gene_pcc": float(args.gate_balanced_gene_pcc),
        "production_ready": bool(
            summaries["all177"]["patient_source_balanced_gene_pcc_median"]
            >= float(args.gate_balanced_gene_pcc)
        ),
        "reference_hellinger_squared_drift_mean": float(reference_drift),
        "theta_kl_drift_mean": float(theta_drift),
        "reference_row_sum_maximum_error": float(
            np.max(np.abs(neutral.sum(axis=1) - 1.0))
        ),
        "theta_row_sum_maximum_error": float(
            np.max(np.abs(theta.sum(axis=1) - 1.0))
        ),
        "reference_nonnegative": bool(np.all(neutral >= 0)),
        "theta_nonnegative": bool(np.all(theta >= 0)),
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
