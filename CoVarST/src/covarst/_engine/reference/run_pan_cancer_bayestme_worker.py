"""Run one resumable worker across all prepared tissue-specific BayesTME queues."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--worker-id", type=int, required=True)
    parser.add_argument("--workers", type=int, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    if args.worker_id < 0 or args.worker_id >= args.workers:
        raise ValueError("invalid worker id")
    runner = Path(str(config["code_root"])) / "upstream_pipeline" / "run_pure_upstream_bayestme_per_slide.py"
    bayes = config["bayestme"]
    environment = os.environ.copy()
    threads = str(int(bayes["threads_per_worker"]))
    environment.update(
        {
            "OMP_NUM_THREADS": threads,
            "MKL_NUM_THREADS": threads,
            "OPENBLAS_NUM_THREADS": threads,
            "TQDM_DISABLE": "1",
        }
    )
    started = time.time()
    completed_tissues = 0
    for tissue_position, tissue in enumerate(config["tissues"]):
        key = str(tissue["key"])
        root = Path(str(config["output_root"])) / key
        prepared = root / "00_prepared_top10k"
        contract = prepared / "prepared_top10k_contract.json"
        if not contract.is_file():
            raise FileNotFoundError(contract)
        shard = root / "00_worker_shards" / f"worker_{args.worker_id:02d}"
        shard_contract = shard / "shard_contract.json"
        if not shard_contract.is_file():
            raise FileNotFoundError(shard_contract)
        samples = pd.read_csv(prepared / "patient_split.csv").sort_values("slide_index")
        assigned = [
            int(value)
            for value in samples["slide_index"].astype(int).tolist()
            if int(value) % int(args.workers) == int(args.worker_id)
        ]
        if not assigned:
            continue
        output = root / "01_per_slide_bayestme"
        command = [
            str(config["python"]["bayestme"]),
            str(runner),
            "--prepared-cache",
            str(shard),
            "--spot-index",
            str(shard / "spot_index.csv"),
            "--output-dir",
            str(output),
            "--maximum-genes",
            str(int(config["panel"]["genes"])),
            "--k-min",
            str(int(bayes["k_min"])),
            "--k-max",
            str(int(bayes["k_max"])),
            "--selection-steps",
            str(int(bayes["selection_steps"])),
            "--selection-samples",
            str(int(bayes["selection_samples"])),
            "--selection-splits",
            str(int(bayes["selection_splits"])),
            "--final-steps",
            str(int(bayes["final_steps"])),
            "--final-samples",
            str(int(bayes["final_samples"])),
            "--spatial-smoothing",
            str(float(bayes["spatial_smoothing"])),
            "--slide-indices",
            ",".join(str(value) for value in assigned),
            "--worker-id",
            f"pan_cancer_w{args.worker_id}_{key}",
            "--skip-queue-metadata-write",
            "--seed",
            str(int(bayes["seed"]) + tissue_position * 1_000_003),
        ]
        print(
            json.dumps(
                {
                    "event": "worker_tissue_start",
                    "worker": args.worker_id,
                    "tissue": key,
                    "slides": assigned,
                }
            ),
            flush=True,
        )
        subprocess.run(command, check=True, env=environment)
        completed_tissues += 1
        print(
            json.dumps(
                {
                    "event": "worker_tissue_complete",
                    "worker": args.worker_id,
                    "tissue": key,
                    "elapsed_seconds": time.time() - started,
                }
            ),
            flush=True,
        )
    print(
        json.dumps(
            {
                "event": "worker_complete",
                "worker": args.worker_id,
                "tissues_with_assigned_slides": completed_tissues,
                "runtime_seconds": time.time() - started,
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
