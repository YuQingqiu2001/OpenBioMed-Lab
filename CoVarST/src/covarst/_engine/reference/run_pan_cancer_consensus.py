"""Fit one independent approximately-40-program consensus per eligible tissue."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd


def run(command: list[str]) -> None:
    print(json.dumps({"command": command}), flush=True)
    subprocess.run(command, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    project_code_root = Path(str(config["code_root"]))
    code_root = project_code_root / "upstream_pipeline"
    source_root = project_code_root / "src"
    existing_pythonpath = os.environ.get("PYTHONPATH", "")
    os.environ["PYTHONPATH"] = os.pathsep.join(
        value for value in (str(source_root), existing_pythonpath) if value
    )
    consensus = config["consensus"]
    status_rows = []
    for tissue_position, tissue in enumerate(config["tissues"]):
        key = str(tissue["key"])
        root = Path(str(config["output_root"])) / key
        prepared = root / "00_prepared_top10k"
        queue = root / "01_per_slide_bayestme"
        complete_paths = sorted((queue / "slides").glob("*/complete.json"))
        records = [json.loads(path.read_text(encoding="utf-8")) for path in complete_paths]
        expected = len(pd.read_csv(prepared / "patient_split.csv"))
        local_programs = int(sum(int(record["selected_k"]) for record in records))
        if len(records) != expected:
            raise RuntimeError(f"{key}: per-slide BayesTME is incomplete")
        if len(records) < 3 or local_programs < int(consensus["minimum_programs"]):
            status = {
                "tissue_key": key,
                "status": "insufficient_independent_slides_or_local_programs",
                "slides": len(records),
                "total_local_programs": local_programs,
                "minimum_programs": int(consensus["minimum_programs"]),
                "reason": tissue.get("known_consensus_limitation", "insufficient local rank for a stable approximately-40-program reference"),
            }
            (root / "consensus_status.json").write_text(
                json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
            )
            status_rows.append(status)
            print(json.dumps(status, ensure_ascii=False), flush=True)
            continue

        candidates = [
            value
            for value in (32, 36, 40, 44, 48)
            if value <= local_programs and value <= int(consensus["maximum_programs"])
        ]
        if not candidates:
            raise RuntimeError(f"{key}: no identifiable consensus candidate")
        target = min(candidates, key=lambda value: (abs(value - int(consensus["preferred_programs"])), value))
        locked = root / "02a_locked_consensus"
        refined = root / "02b_final_consensus"
        locked_contract = locked / "consensus_contract.json"
        if not locked_contract.is_file():
            if locked.exists() and any(locked.iterdir()):
                raise FileExistsError(f"{key}: incomplete locked consensus output {locked}")
            if locked.exists():
                locked.rmdir()
            run(
                [
                    sys.executable,
                    str(code_root / "fit_bayestme_top10k_direct_consensus_locked_split.py"),
                    "--queue-dir",
                    str(queue),
                    "--master-spot-index",
                    str(prepared / "spot_index.csv"),
                    "--split-manifest",
                    str(prepared / "patient_split.csv"),
                    "--output-dir",
                    str(locked),
                    "--candidates",
                    *[str(value) for value in candidates],
                    "--target-programs",
                    str(target),
                    "--selection-tolerance",
                    str(float(consensus["selection_tolerance"])),
                    "--candidate-epochs",
                    str(int(consensus["candidate_epochs"])),
                    "--mapping-epochs",
                    str(int(consensus["mapping_epochs"])),
                    "--final-epochs",
                    str(int(consensus["final_epochs"])),
                    "--gate-balanced-gene-pcc",
                    str(float(consensus["gate_balanced_gene_pcc"])),
                    "--seed",
                    str(int(consensus["seed_locked"]) + tissue_position * 1009),
                ]
            )
        refined_contract = refined / "consensus_contract.json"
        if not refined_contract.is_file():
            if refined.exists() and any(refined.iterdir()):
                raise FileExistsError(f"{key}: incomplete refined consensus output {refined}")
            if refined.exists():
                refined.rmdir()
            run(
                [
                    sys.executable,
                    str(code_root / "refine_direct_consensus_spot_nmf.py"),
                    "--queue-dir",
                    str(queue),
                    "--master-spot-index",
                    str(prepared / "spot_index.csv"),
                    "--split-manifest",
                    str(prepared / "patient_split.csv"),
                    "--input-consensus-dir",
                    str(locked),
                    "--output-dir",
                    str(refined),
                    "--steps",
                    str(int(consensus["spot_nmf_steps"])),
                    "--spots-per-step",
                    str(int(consensus["spots_per_step"])),
                    "--genes-per-step",
                    str(int(consensus["genes_per_step"])),
                    "--gate-balanced-gene-pcc",
                    str(float(consensus["gate_balanced_gene_pcc"])),
                    "--seed",
                    str(int(consensus["seed_refine"]) + tissue_position * 1009),
                ]
            )
        contract = json.loads(refined_contract.read_text(encoding="utf-8"))
        status = {
            "tissue_key": key,
            "status": "complete",
            "slides": len(records),
            "total_local_programs": local_programs,
            "programs": int(contract["program_count_locked_before_refinement"]),
            "balanced_gene_pcc_all": float(
                contract["recovery_against_original_local_inverse"]["all177"][
                    "patient_source_balanced_gene_pcc_median"
                ]
            ),
            "reference": str(refined / "consensus_reference_probability_10k.npy"),
            "theta_rna": str(refined / "theta_rna_consensus.npy"),
            "theta_composition": str(refined / "theta_composition_consensus.npy"),
        }
        (root / "consensus_status.json").write_text(
            json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        status_rows.append(status)
        print(json.dumps(status, ensure_ascii=False), flush=True)
    pd.DataFrame(status_rows).to_csv(
        Path(str(config["output_root"])) / "consensus_completion.csv", index=False
    )


if __name__ == "__main__":
    main()
