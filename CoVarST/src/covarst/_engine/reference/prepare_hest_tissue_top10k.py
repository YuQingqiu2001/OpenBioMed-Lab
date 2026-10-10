"""Prepare independent marker-first/HVG-fill 10K HEST caches by cancer type.

Only original nonnegative integer counts are serialized.  No BayesTME basis,
theta, image feature, or prior reference is read.  Each cancer type receives a
separate gene panel and a separate spot index.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


EPS = 1.0e-12


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--tissue", action="append", default=[])
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def atomic_npy(path: Path, values: np.ndarray) -> None:
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("wb") as handle:
        np.save(handle, values)
    os.replace(temporary, path)


def atomic_sparse(path: Path, values: sparse.csr_matrix) -> None:
    temporary = path.with_name(path.stem + ".partial" + path.suffix)
    sparse.save_npz(temporary, values, compressed=True)
    os.replace(temporary, path)


def as_csr(adata: ad.AnnData) -> sparse.csr_matrix:
    values = adata.X[:]
    if sparse.issparse(values):
        matrix = values.tocsr()
    else:
        matrix = sparse.csr_matrix(np.asarray(values))
    matrix.sum_duplicates()
    matrix.sort_indices()
    if np.any(matrix.data < 0):
        raise ValueError("negative expression values are not raw counts")
    if not np.allclose(matrix.data, np.rint(matrix.data)):
        raise ValueError("non-integer expression values are not raw UMI counts")
    return matrix


def clean_gene_names(values: pd.Index) -> np.ndarray:
    return np.asarray([str(value).strip() for value in values], dtype=str)


def subset_and_sum_duplicates(
    matrix: sparse.csr_matrix,
    source_genes: np.ndarray,
    target_genes: list[str],
    *,
    allow_missing: bool = False,
) -> sparse.csr_matrix:
    target = {gene: index for index, gene in enumerate(target_genes)}
    source_columns: list[int] = []
    output_columns: list[int] = []
    for column, gene in enumerate(source_genes):
        if gene in target:
            source_columns.append(column)
            output_columns.append(target[gene])
    observed = set(source_genes[source_columns].tolist())
    missing = sorted(set(target_genes) - observed)
    if missing and not allow_missing:
        raise RuntimeError(f"selected genes absent from slide: {missing[:20]}")
    selected = matrix[:, np.asarray(source_columns, dtype=np.int64)]
    projection = sparse.csr_matrix(
        (
            np.ones(len(source_columns), dtype=np.float32),
            (
                np.arange(len(source_columns), dtype=np.int64),
                np.asarray(output_columns, dtype=np.int64),
            ),
        ),
        shape=(len(source_columns), len(target_genes)),
    )
    result = (selected @ projection).tocsr()
    result.sum_duplicates()
    result.sort_indices()
    return result


def choose_coordinates(adata: ad.AnnData) -> tuple[np.ndarray, str]:
    if {"array_col", "array_row"}.issubset(adata.obs.columns):
        coordinates = adata.obs[["array_col", "array_row"]].to_numpy(dtype=np.float32)
        if np.isfinite(coordinates).all() and len(
            np.unique(np.rint(coordinates).astype(np.int64), axis=0)
        ) == len(coordinates):
            return coordinates, "obs.array_col,array_row"
    coordinates = np.asarray(adata.obsm["spatial"], dtype=np.float32)
    if coordinates.shape != (adata.n_obs, 2) or not np.isfinite(coordinates).all():
        raise ValueError("invalid spatial coordinates")
    return coordinates, "obsm.spatial"


def split_manifest(samples: pd.DataFrame, seed: int) -> tuple[pd.DataFrame, str]:
    frame = samples.copy()
    patient = frame["patient"].fillna("").astype(str).str.strip()
    patient = patient.mask(patient.isin({"", "nan", "NA", "Unknown", "UNKNOWN"}))
    if patient.nunique(dropna=True) < 3:
        patient = pd.Series(
            [f"sample:{sample_id}" for sample_id in frame["sample_id"].astype(str)],
            index=frame.index,
        )
        grouping = "sample_fallback_because_fewer_than_three_patient_ids"
    else:
        patient = patient.fillna(
            pd.Series(
                [f"sample:{sample_id}" for sample_id in frame["sample_id"].astype(str)],
                index=frame.index,
            )
        )
        grouping = "patient_when_available_with_sample_fallback"
    frame["patient"] = patient.astype(str)
    groups = np.asarray(sorted(frame["patient"].unique()), dtype=object)
    rng = np.random.default_rng(int(seed))
    rng.shuffle(groups)
    labels: dict[str, str] = {}
    if len(groups) >= 3:
        test_count = max(1, int(round(0.20 * len(groups))))
        val_count = max(1, int(round(0.20 * len(groups))))
        if test_count + val_count >= len(groups):
            test_count = 1
            val_count = 1
        for value in groups[:test_count]:
            labels[str(value)] = "test"
        for value in groups[test_count : test_count + val_count]:
            labels[str(value)] = "val"
        for value in groups[test_count + val_count :]:
            labels[str(value)] = "train"
    elif len(groups) == 2:
        labels[str(groups[0])] = "val"
        labels[str(groups[1])] = "train"
    elif len(groups) == 1:
        labels[str(groups[0])] = "train"
    frame["split"] = frame["patient"].map(labels)
    return frame, grouping


def prepare_tissue(
    config: dict[str, object],
    tissue: dict[str, object],
    metadata: pd.DataFrame,
    marker_catalog: dict[str, object],
    tissue_position: int,
) -> dict[str, object]:
    key = str(tissue["key"])
    code = str(tissue["oncotree_code"])
    output_root = Path(str(config["output_root"]))
    output_dir = output_root / key / "00_prepared_top10k"
    output_dir.mkdir(parents=True, exist_ok=True)
    contract_path = output_dir / "prepared_top10k_contract.json"
    if contract_path.exists():
        return json.loads(contract_path.read_text(encoding="utf-8"))
    unexpected = list(output_dir.iterdir())
    if unexpected:
        raise FileExistsError(
            f"refusing partially populated prepared directory without contract: {output_dir}"
        )

    panel_genes = int(config["panel"]["genes"])
    species = str(config["panel"]["species"])
    technologies = {str(value) for value in tissue["include_technologies"]}
    rows = metadata.loc[
        metadata["species"].astype(str).eq(species)
        & metadata["oncotree_code"].astype(str).eq(code)
    ].copy()
    rows["sample_id"] = rows["id"].astype(str)
    rows["h5ad_path"] = rows["sample_id"].map(
        lambda value: str(Path(str(config["h5ad_root"])) / f"{value}.h5ad")
    )
    rows["included"] = rows["st_technology"].astype(str).isin(technologies)
    rows["exclusion_reason"] = np.where(
        rows["included"], "", "technology_not_selected_for_true_10k_panel"
    )
    minimum_assay_genes = int(tissue.get("minimum_assay_genes", panel_genes))
    too_small = rows["nb_genes"].fillna(0).astype(float) < minimum_assay_genes
    rows.loc[rows["included"] & too_small, "included"] = False
    rows.loc[
        rows["exclusion_reason"].eq("") & too_small,
        "exclusion_reason",
    ] = "assay_has_fewer_than_10000_genes"
    missing_file = ~rows["h5ad_path"].map(lambda value: Path(value).is_file())
    rows.loc[rows["included"] & missing_file, "included"] = False
    rows.loc[
        rows["exclusion_reason"].eq("") & missing_file,
        "exclusion_reason",
    ] = "h5ad_missing"
    included = rows.loc[rows["included"]].sort_values("sample_id", kind="stable").copy()
    if included.empty:
        raise RuntimeError(f"{key}: no eligible samples")

    gene_sets: list[set[str]] = []
    observed_shapes: list[dict[str, object]] = []
    for row in included.itertuples(index=False):
        adata = ad.read_h5ad(row.h5ad_path, backed="r")
        genes = clean_gene_names(adata.var_names)
        gene_sets.append(set(genes.tolist()))
        observed_shapes.append(
            {
                "sample_id": str(row.sample_id),
                "spots": int(adata.n_obs),
                "genes": int(adata.n_vars),
                "unique_genes": int(len(set(genes.tolist()))),
            }
        )
        adata.file.close()
        print(
            json.dumps(
                {
                    "stage": "gene_universe_scan",
                    "tissue": key,
                    "sample_id": str(row.sample_id),
                    "completed": len(observed_shapes),
                    "total": int(len(included)),
                }
            ),
            flush=True,
        )
    common = set.intersection(*gene_sets)
    union = set.union(*gene_sets)
    exclude_patterns = [
        re.compile(str(pattern)) for pattern in config["panel"]["exclude_gene_patterns"]
    ]
    marker_entry = marker_catalog["tissues"][code]
    mandatory_order = [str(value) for value in marker_entry["mandatory_markers"]]
    mandatory_set = set(mandatory_order)
    forced = list(
        dict.fromkeys(
            gene
            for gene in mandatory_order
            if gene in union and gene and not gene.startswith("__")
        )
    )
    forced_set = set(forced)
    hvg_genes = sorted(
        gene
        for gene in common
        if gene
        and not gene.startswith("__")
        and not any(pattern.match(gene) for pattern in exclude_patterns)
    )
    eligible_panel_universe = set(hvg_genes) | forced_set
    if len(eligible_panel_universe) < panel_genes:
        raise RuntimeError(
            f"{key}: only {len(eligible_panel_universe)} marker/HVG eligible genes "
            f"for a {panel_genes} panel"
        )

    gene_to_index = {gene: index for index, gene in enumerate(hvg_genes)}
    forced_to_index = {gene: index for index, gene in enumerate(forced)}
    score = np.zeros(len(hvg_genes), dtype=np.float64)
    detected_spots = np.zeros(len(hvg_genes), dtype=np.int64)
    detected_slides = np.zeros(len(hvg_genes), dtype=np.int32)
    marker_detected_spots = np.zeros(len(forced), dtype=np.int64)
    marker_detected_slides = np.zeros(len(forced), dtype=np.int32)
    total_spots = 0
    for score_position, row in enumerate(included.itertuples(index=False), start=1):
        adata = ad.read_h5ad(row.h5ad_path)
        matrix = as_csr(adata)
        source_genes = clean_gene_names(adata.var_names)
        block = subset_and_sum_duplicates(matrix, source_genes, hvg_genes).astype(
            np.float32
        )
        marker_block = subset_and_sum_duplicates(
            matrix, source_genes, forced, allow_missing=True
        )
        library = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64)
        transformed = block.multiply(
            (1.0e4 / np.maximum(library, 1.0))[:, None]
        ).tocsr()
        transformed.data = np.log1p(transformed.data)
        mean = np.asarray(transformed.mean(axis=0)).ravel()
        second = np.asarray(transformed.multiply(transformed).mean(axis=0)).ravel()
        variance = np.maximum(second - np.square(mean), 0.0)
        order = np.argsort(variance, kind="stable")
        rank = np.empty(len(hvg_genes), dtype=np.float64)
        rank[order] = np.linspace(0.0, 1.0, len(hvg_genes))
        score += rank
        detection = np.asarray((block > 0).sum(axis=0)).ravel().astype(np.int64)
        detected_spots += detection
        detected_slides += (detection > 0).astype(np.int32)
        marker_detection = np.asarray((marker_block > 0).sum(axis=0)).ravel().astype(
            np.int64
        )
        marker_detected_spots += marker_detection
        marker_detected_slides += (marker_detection > 0).astype(np.int32)
        total_spots += int(block.shape[0])
        del adata, matrix, block, marker_block, transformed
        print(
            json.dumps(
                {
                    "stage": "hvg_scoring",
                    "tissue": key,
                    "sample_id": str(row.sample_id),
                    "completed": score_position,
                    "total": int(len(included)),
                }
            ),
            flush=True,
        )
    score /= len(included)

    absent_markers = [gene for gene in mandatory_order if gene not in union]
    remaining = [gene for gene in hvg_genes if gene not in forced_set]
    remaining.sort(
        key=lambda gene: (
            -score[gene_to_index[gene]],
            -detected_slides[gene_to_index[gene]],
            -detected_spots[gene_to_index[gene]],
            gene,
        )
    )
    selected = forced + remaining[: panel_genes - len(forced)]
    if len(selected) != panel_genes or len(set(selected)) != panel_genes:
        raise RuntimeError(f"{key}: failed to construct a unique 10K panel")

    minimum_core = set(str(value) for value in marker_entry["minimum_core_markers"])
    tissue_markers = set(str(value) for value in marker_entry["tissue_specific_markers"])
    universal_tme = set(str(value) for value in marker_catalog["universal_tme_markers"])
    measured_slides = {
        gene: int(sum(gene in genes for genes in gene_sets))
        for gene in eligible_panel_universe
    }
    selection_rows = []
    selected_set = set(selected)
    selected_index = {gene: index for index, gene in enumerate(selected)}
    for gene in sorted(eligible_panel_universe):
        scopes = []
        if gene in minimum_core:
            scopes.append("minimum_core")
        if gene in tissue_markers:
            scopes.append("tissue_specific")
        if gene in universal_tme:
            scopes.append("universal_tme")
        selection_rows.append(
            {
                "gene": gene,
                "selected": gene in selected_set,
                "panel_index": selected_index.get(gene, -1),
                "selection_reason": "mandatory_marker" if gene in forced_set else (
                    "slide_balanced_hvg" if gene in selected_set else "not_selected"
                ),
                "marker_scopes": ";".join(scopes),
                "slide_balanced_hvg_score": (
                    float(score[gene_to_index[gene]]) if gene in gene_to_index else np.nan
                ),
                "detected_slides": int(
                    marker_detected_slides[forced_to_index[gene]]
                    if gene in forced_to_index
                    else detected_slides[gene_to_index[gene]]
                ),
                "detected_spots": int(
                    marker_detected_spots[forced_to_index[gene]]
                    if gene in forced_to_index
                    else detected_spots[gene_to_index[gene]]
                ),
                "measured_slides": measured_slides[gene],
                "assay_missing_slides": int(len(included) - measured_slides[gene]),
            }
        )
    selection_frame = pd.DataFrame(selection_rows).sort_values(
        ["selected", "panel_index", "slide_balanced_hvg_score", "gene"],
        ascending=[False, True, False, True],
        kind="stable",
    )
    selection_frame.to_csv(output_dir / "gene_selection_top10k.csv", index=False)
    pd.DataFrame(
        {
            "panel_index": np.arange(panel_genes, dtype=np.int64),
            "gene": selected,
            "selection_reason": [
                "mandatory_marker" if gene in forced_set else "slide_balanced_hvg"
                for gene in selected
            ],
            "measured_slides": [measured_slides[gene] for gene in selected],
            "assay_missing_slides": [
                int(len(included) - measured_slides[gene]) for gene in selected
            ],
        }
    ).to_csv(output_dir / "gene_panel_10k.csv", index=False)
    pd.DataFrame({"gene": absent_markers}).to_csv(
        output_dir / "mandatory_markers_absent_from_all_eligible_slides.csv", index=False
    )

    count_blocks: list[sparse.csr_matrix] = []
    coordinate_blocks: list[np.ndarray] = []
    spot_tables: list[pd.DataFrame] = []
    sample_records: list[dict[str, object]] = []
    global_offset = 0
    for slide_index, row in enumerate(included.itertuples(index=False)):
        adata = ad.read_h5ad(row.h5ad_path)
        matrix = as_csr(adata)
        source_genes = clean_gene_names(adata.var_names)
        counts = subset_and_sum_duplicates(
            matrix, source_genes, selected, allow_missing=True
        )
        coordinates, coordinate_source = choose_coordinates(adata)
        barcodes = adata.obs_names.astype(str).to_numpy()
        keep = np.asarray(counts.sum(axis=1)).ravel() > 0
        removed = int(np.sum(~keep))
        counts = counts[keep].tocsr()
        coordinates = coordinates[keep]
        barcodes = barcodes[keep]
        patient = str(row.patient).strip()
        if patient in {"", "nan", "NA", "Unknown", "UNKNOWN"}:
            patient = f"sample:{row.sample_id}"
        count_blocks.append(counts)
        coordinate_blocks.append(coordinates.astype(np.float32))
        spot_tables.append(
            pd.DataFrame(
                {
                    "global_index": np.arange(
                        global_offset, global_offset + counts.shape[0], dtype=np.int64
                    ),
                    "slide_index": int(slide_index),
                    "local_index": np.arange(counts.shape[0], dtype=np.int64),
                    "source": "HEST",
                    "tissue_key": key,
                    "oncotree_code": code,
                    "sample_id": str(row.sample_id),
                    "patient": patient,
                    "technology": str(row.st_technology),
                    "barcode": barcodes,
                }
            )
        )
        sample_records.append(
            {
                "slide_index": int(slide_index),
                "sample_id": str(row.sample_id),
                "source": "HEST",
                "patient": patient,
                "technology": str(row.st_technology),
                "spots": int(counts.shape[0]),
                "removed_zero_panel_spots": removed,
                "coordinate_source": coordinate_source,
                "h5ad_path": str(row.h5ad_path),
                "original_genes": int(adata.n_vars),
            }
        )
        global_offset += counts.shape[0]
        del adata, matrix, counts
        print(
            json.dumps(
                {
                    "stage": "cache_build",
                    "tissue": key,
                    "sample_id": str(row.sample_id),
                    "completed": slide_index + 1,
                    "total": int(len(included)),
                }
            ),
            flush=True,
        )

    counts10k = sparse.vstack(count_blocks, format="csr")
    coordinates = np.concatenate(coordinate_blocks, axis=0).astype(np.float32)
    spot_index = pd.concat(spot_tables, ignore_index=True)
    samples = pd.DataFrame(sample_records)
    samples, split_grouping = split_manifest(
        samples, seed=int(config["bayestme"]["seed"]) + tissue_position * 1009
    )
    split_lookup = samples.set_index("sample_id")[["patient", "split"]]
    spot_index["patient"] = spot_index["sample_id"].map(split_lookup["patient"])
    if len(spot_index) != counts10k.shape[0] or len(coordinates) != len(spot_index):
        raise RuntimeError(f"{key}: combined row contracts disagree")
    if not np.allclose(counts10k.data, np.rint(counts10k.data)):
        raise RuntimeError(f"{key}: prepared counts are not integer UMI")

    zero = sparse.csr_matrix(counts10k.shape, dtype=counts10k.dtype)
    atomic_sparse(output_dir / "fit_core.npz", counts10k)
    atomic_sparse(output_dir / "tune_core.npz", zero)
    atomic_sparse(output_dir / "audit_core.npz", zero)
    atomic_npy(output_dir / "discovery_gene_order.npy", np.asarray(selected, dtype=str))
    atomic_npy(output_dir / "coordinates.npy", coordinates)
    with (output_dir / "metadata.npz.partial").open("wb") as handle:
        np.savez(
            handle,
            slide_index=spot_index["slide_index"].to_numpy(dtype=np.int32),
            sample_ids=np.asarray(
                samples.sort_values("slide_index")["sample_id"].astype(str).tolist(),
                dtype=str,
            ),
            sources=np.asarray(["HEST"] * len(samples), dtype=str),
        )
    os.replace(output_dir / "metadata.npz.partial", output_dir / "metadata.npz")
    spot_index.to_csv(output_dir / "spot_index.csv", index=False)
    samples.to_csv(output_dir / "patient_split.csv", index=False)
    rows.to_csv(output_dir / "sample_eligibility.csv", index=False)
    pd.DataFrame(observed_shapes).to_csv(output_dir / "input_h5ad_shapes.csv", index=False)

    contract = {
        "method": "hest_tissue_marker_fixed_union_slide_balanced_hvg10k_raw_umi_v2",
        "tissue_key": key,
        "oncotree_code": code,
        "species": species,
        "genes": int(counts10k.shape[1]),
        "slides": int(len(samples)),
        "spots": int(counts10k.shape[0]),
        "patients_or_sample_fallback_groups": int(samples["patient"].nunique()),
        "split_counts": {
            str(name): int(value) for name, value in samples["split"].value_counts().items()
        },
        "split_grouping": split_grouping,
        "technologies": {
            str(name): int(value)
            for name, value in samples["technology"].value_counts().items()
        },
        "common_gene_universe": int(len(common)),
        "union_gene_universe": int(len(union)),
        "eligible_common_hvg_universe": int(len(hvg_genes)),
        "eligible_marker_plus_hvg_universe": int(len(eligible_panel_universe)),
        "mandatory_marker_catalog_count": int(len(mandatory_order)),
        "mandatory_markers_retained": int(len(forced)),
        "mandatory_markers_absent_from_all_eligible_slides": absent_markers,
        "mandatory_markers_not_measured_in_every_slide": int(
            sum(measured_slides[gene] < len(included) for gene in forced)
        ),
        "assay_missing_marker_values_are_explicit_zero": True,
        "hvg_fill_count": int(panel_genes - len(forced)),
        "selection": str(config["panel"]["selection"]),
        "raw_integer_umi": True,
        "old_basis_theta_checkpoint_read": 0,
        "images_read": 0,
        "input_h5ad_count": int(len(samples)),
        "excluded_metadata_rows": int(np.sum(~rows["included"])),
        "total_umi": int(counts10k.sum()),
        "nnz": int(counts10k.nnz),
        "marker_catalog_sha256": sha256(Path(str(config["marker_catalog"]))),
        "fit_core_sha256": sha256(output_dir / "fit_core.npz"),
        "gene_order_sha256": sha256(output_dir / "discovery_gene_order.npy"),
        "spot_index_sha256": sha256(output_dir / "spot_index.csv"),
        "known_consensus_limitation": tissue.get("known_consensus_limitation"),
    }
    atomic_json(contract_path, contract)
    print(json.dumps(contract, ensure_ascii=False), flush=True)
    return contract


def main() -> None:
    args = parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    metadata = pd.read_csv(config["metadata"])
    marker_catalog = json.loads(
        Path(str(config["marker_catalog"])).read_text(encoding="utf-8")
    )
    requested = set(args.tissue)
    contracts = []
    for position, tissue in enumerate(config["tissues"]):
        if requested and str(tissue["key"]) not in requested:
            continue
        contracts.append(
            prepare_tissue(config, tissue, metadata, marker_catalog, position)
        )
    print(
        json.dumps(
            {
                "prepared": [contract["tissue_key"] for contract in contracts],
                "count": len(contracts),
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
