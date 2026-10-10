#!/usr/bin/env python3
"""Run the P1 MV4.2 BayesGraph 278-gene gated prototype.

V55 and the 2 um direct-cell reference are fitted independently.  V55 uses
buffered contiguous spatial folds on fractional spot mass; 2 um uses exact
integer split-half validation followed by empirical-Bayes shrinkage.  No 2 um
expression is read by the V55 fit.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace
import time

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import torch

import cell_resolved_visium55_mv3 as mv3
import p1_visium55_mv4_panel_prototype as proto
import spot2cell_bayesgraph_core as bg
import spot2cell_mv4_core as mv4
import spot2cell_stinr_core as mass_core


SCHEMA = "P1_MV42_BayesGraph_panel_v1"
CLASS_NAMES = mv3.CLASS_NAMES
MAX_PROGRAMS = 8
EPS = 1.0e-10


def log(message: str) -> None:
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), message, flush=True)


def parse_args() -> argparse.Namespace:
    root = Path('DATA/Visium_HD')
    mv2 = root / "cell_resolved_spot2cell_mv2/P1/visium55_v2"
    mv3_root = root / "cell_resolved_spot2cell_mv3/P1/visium55_v3"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-2um", type=Path, default=root / "cell_resolved_002um_voronoi/P1/P1.Visium2_direct.full_support_input.h5ad")
    parser.add_argument("--mv2-h5ad", type=Path, default=mv2 / "P1.visium55.spot2cell_mv2.h5ad")
    parser.add_argument("--mv2-allocation", type=Path, default=mv2 / "P1.visium55.spot2cell_mv2.allocation.h5")
    parser.add_argument("--spot-positions", type=Path, default=mv2 / "P1.visium55.spots.tsv.gz")
    parser.add_argument("--mv3-state", type=Path, default=mv3_root / "P1.visium55.spot2cell_mv3.state.h5")
    parser.add_argument("--panel-template", type=Path, default=mv3_root / "P1.visium55.spot2cell_mv3.P_state_278.h5ad")
    parser.add_argument("--output-dir", type=Path, default=root / "cell_resolved_spot2cell_mv42/P1/bayesgraph_panel_v1")
    parser.add_argument("--ranks", default="0,2,4")
    parser.add_argument("--trend-lambdas", default="0.05,0.2")
    parser.add_argument("--tau-grid", default="0,2,5,10,20,50")
    parser.add_argument("--spatial-folds", type=int, default=5)
    parser.add_argument("--spatial-buffer", type=float, default=1.5)
    parser.add_argument("--cv-steps-v55", type=int, default=100)
    parser.add_argument("--final-steps-v55", type=int, default=320)
    parser.add_argument("--cv-steps-2um", type=int, default=160)
    parser.add_argument("--final-steps-2um", type=int, default=360)
    parser.add_argument("--permutations", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260828)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--force-rank-v55", type=int, default=None,
                        help="development-only override used to exercise the nonzero-rank path")
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def parse_numbers(value: str, kind=float):
    return [kind(item.strip()) for item in value.split(",") if item.strip()]


def sha256(path: Path, block: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(block)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(value: dict, path: Path) -> None:
    partial = path.with_suffix(path.suffix + ".partial")
    with partial.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, default=json_default)
        handle.write("\n")
    os.replace(partial, path)


def json_default(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def write_tsv(rows, path: Path) -> None:
    partial = path.with_suffix(path.suffix + ".partial")
    pd.DataFrame(rows).to_csv(partial, sep="\t", index=False)
    os.replace(partial, path)


def write_h5ad_atomic(data: ad.AnnData, path: Path) -> None:
    partial = path.with_suffix(path.suffix + ".partial")
    if path.exists() or partial.exists():
        raise FileExistsError(path)
    data.write_h5ad(partial, compression="gzip", compression_opts=4)
    os.replace(partial, path)


def load_v55(args):
    data, old_summary = proto.load_data(args)
    with h5py.File(args.panel_template, "r") as handle:
        panel_genes = mv3.read_var_names(handle).astype(str)
    lookup = {str(gene): index for index, gene in enumerate(data["gene_names"])}
    panel_index = np.asarray([lookup[gene] for gene in panel_genes], np.int32)
    counts = np.asarray(data["matrix"][panel_index].T.toarray(), np.float64)
    profiles = np.maximum(data["profiles"][:, panel_index], EPS)
    profiles /= profiles.sum(axis=1, keepdims=True)
    geometry = data["geometry"].tocsr()
    capacity = data["capacity"].astype(np.float64)
    type_mass = np.zeros((geometry.shape[0], len(CLASS_NAMES)), np.float64)
    for value in range(len(CLASS_NAMES)):
        type_mass[:, value] = np.asarray(
            geometry @ (capacity * (data["class_id"] == value))
        ).ravel()
    total = type_mass.sum(axis=1)
    composition = np.divide(type_mass, total[:, None], out=np.zeros_like(type_mass), where=total[:, None] > 0)
    return data, old_summary, panel_genes, panel_index, counts, profiles.astype(np.float32), composition.astype(np.float32), total


def candidate_grid(ranks, lambdas):
    result = [(0, 0.0)]
    result.extend((rank, lam) for rank in ranks if rank > 0 for lam in lambdas)
    return result


def select_one_se(rows: list[dict]) -> dict:
    best = min(rows, key=lambda row: row["mean_nll"])
    threshold = best["mean_nll"] + best.get("se_nll", 0.0)
    eligible = [row for row in rows if row["mean_nll"] <= threshold]
    return sorted(eligible, key=lambda row: (row["rank"], -row["lambda_trend"]))[0]


def fit_v55(args, counts, composition, profiles, coords, component, graph):
    library = counts.sum(axis=1)
    occupied = composition.sum(axis=1) > 0
    folds = mv4.contiguous_spatial_folds(coords, args.spatial_folds)
    rows = []
    for rank, lam in candidate_grid(parse_numbers(args.ranks, int), parse_numbers(args.trend_lambdas, float)):
        fold_losses = []
        audits = []
        if rank == 0:
            probability = composition @ profiles
            probability /= np.maximum(probability.sum(axis=1, keepdims=True), EPS)
            for fold in range(args.spatial_folds):
                test = np.flatnonzero((folds == fold) & occupied & (library > 0))
                fold_losses.append(bg.dense_multinomial_nll(counts, probability, test))
        else:
            for fold in range(args.spatial_folds):
                train, test_mask = mv4._training_mask_with_buffer(
                    coords, folds, fold, args.spatial_buffer,
                )
                eligible = np.flatnonzero(train & occupied & (library > 0))
                test = np.flatnonzero(test_mask & occupied & (library > 0))
                log(f"V55 CV rank={rank} lambda={lam:g} fold={fold + 1}/{args.spatial_folds}")
                fit = bg.fit_spot_program(
                    counts, composition, profiles, graph,
                    rank=rank, lambda_trend=lam, steps=args.cv_steps_v55,
                    eligible=eligible, seed=args.seed + rank * 1000 + int(lam * 100) * 10 + fold,
                    device=args.device,
                )
                fold_losses.append(bg.dense_multinomial_nll(counts, fit.probabilities, test))
                audits.append(fit.audit)
                if fit.model is not None:
                    fit.model.to("cpu")
                del fit
                torch.cuda.empty_cache(); gc.collect()
        rows.append({
            "rank": int(rank), "lambda_trend": float(lam),
            "fold_nll": fold_losses, "mean_nll": float(np.mean(fold_losses)),
            "se_nll": float(np.std(fold_losses, ddof=1) / np.sqrt(len(fold_losses))),
            "fit_audits": audits,
        })
    selected = select_one_se(rows)
    if args.force_rank_v55 is not None:
        forced = [row for row in rows if row["rank"] == args.force_rank_v55]
        if not forced:
            raise ValueError("forced V55 rank was not evaluated")
        selected = min(forced, key=lambda row: row["mean_nll"])
    log(f"V55 selected rank={selected['rank']} lambda={selected['lambda_trend']:g}")
    final = bg.fit_spot_program(
        counts, composition, profiles, graph,
        rank=selected["rank"], lambda_trend=selected["lambda_trend"],
        steps=args.final_steps_v55, seed=args.seed + 9000,
        device=args.device,
    )
    probability = final.probabilities.copy()
    zero = probability.sum(axis=1) <= 0
    observed_probability = np.divide(counts, library[:, None], out=np.zeros_like(counts), where=library[:, None] > 0)
    probability[zero] = observed_probability[zero]
    remaining = probability.sum(axis=1) <= 0
    probability[remaining] = np.mean(profiles, axis=0)
    probability /= np.maximum(probability.sum(axis=1, keepdims=True), EPS)
    imputed = library[:, None] * probability
    return final, imputed.astype(np.float32), rows, selected, folds


def direct_tau_nll(train_dense, train_library, test_coo, prior, tau):
    denominator = train_library[test_coo.row] + tau
    numerator = train_dense[test_coo.row, test_coo.col] + tau * prior[test_coo.row, test_coo.col]
    fallback = denominator <= 0
    probability = np.divide(numerator, np.maximum(denominator, EPS))
    if np.any(fallback):
        probability[fallback] = prior[test_coo.row[fallback], test_coo.col[fallback]]
    return float(-(test_coo.data * np.log(np.maximum(probability, EPS))).sum() / max(float(test_coo.data.sum()), 1.0))


def fit_2um(args, counts, coords, class_id, graph):
    train, test = bg.split_sparse_counts(counts, args.seed + 101)
    profiles_train = bg.class_profiles(train, class_id, len(CLASS_NAMES))
    train_dense = train.toarray().astype(np.float32)
    train_library = np.asarray(train.sum(axis=1)).ravel().astype(np.float64)
    test_coo = test.tocoo()
    tau_grid = parse_numbers(args.tau_grid, float)
    rows = []
    for rank, lam in candidate_grid(parse_numbers(args.ranks, int), parse_numbers(args.trend_lambdas, float)):
        log(f"2 um split-half rank={rank} lambda={lam:g}")
        fit = bg.fit_direct_program(
            train, class_id, profiles_train, graph,
            rank=rank, lambda_trend=lam, steps=args.cv_steps_2um,
            seed=args.seed + 2000 + rank * 100 + int(lam * 100), device=args.device,
        )
        tau_rows = [{"tau": tau, "heldout_nll": direct_tau_nll(
            train_dense, train_library, test_coo, fit.probabilities, tau,
        )} for tau in tau_grid]
        chosen_tau = min(tau_rows, key=lambda row: row["heldout_nll"])
        rows.append({"rank": int(rank), "lambda_trend": float(lam),
                     "mean_nll": chosen_tau["heldout_nll"], "se_nll": 0.0,
                     "selected_tau": chosen_tau["tau"], "tau_grid": tau_rows,
                     "fit_audit": fit.audit})
        if fit.model is not None:
            fit.model.to("cpu")
        del fit
        torch.cuda.empty_cache(); gc.collect()
    best_nll = min(row["mean_nll"] for row in rows)
    eligible = [row for row in rows if row["mean_nll"] <= best_nll + 1.0e-3]
    selected = sorted(eligible, key=lambda row: (row["rank"], -row["lambda_trend"], row["selected_tau"]))[0]
    log(f"2 um selected rank={selected['rank']} lambda={selected['lambda_trend']:g} tau={selected['selected_tau']:g}")
    profiles_full = bg.class_profiles(counts, class_id, len(CLASS_NAMES))
    final = bg.fit_direct_program(
        counts, class_id, profiles_full, graph,
        rank=selected["rank"], lambda_trend=selected["lambda_trend"],
        steps=args.final_steps_2um, seed=args.seed + 9900, device=args.device,
    )
    posterior = bg.posterior_shrinkage(counts, final.probabilities, selected["selected_tau"])
    return final, posterior, rows, selected


def install_new_programs(data, panel_index, composition, fit):
    rank = fit.v.shape[1]
    support = np.count_nonzero(composition > 1.0e-5, axis=0)
    program_count = np.where(support >= 40, rank, 0).astype(np.int32)
    loadings = np.zeros((len(CLASS_NAMES), MAX_PROGRAMS, len(data["gene_names"])), np.float64)
    if rank:
        loadings[:, :rank, panel_index] = fit.v
    group_bins = []
    group_classes = []
    group_scores = []
    for value in range(len(CLASS_NAMES)):
        if program_count[value] <= 0:
            continue
        bins = np.flatnonzero(composition[:, value] > 1.0e-5)
        group_bins.append(bins)
        group_classes.append(np.full(len(bins), value, np.int32))
        block = np.zeros((len(bins), MAX_PROGRAMS), np.float64)
        block[:, :rank] = fit.w[bins, value]
        group_scores.append(block)
    data["program_count"] = program_count
    data["loadings"] = loadings
    data["group_bin"] = np.concatenate(group_bins) if group_bins else np.empty(0, np.int32)
    data["group_class"] = np.concatenate(group_classes) if group_classes else np.empty(0, np.int32)
    data["group_scores"] = np.concatenate(group_scores) if group_scores else np.empty((0, MAX_PROGRAMS), np.float64)
    return program_count, support


def overlap_weighted_base_state(data, fit):
    rank = fit.v.shape[1]
    state = np.zeros((len(data["class_id"]), MAX_PROGRAMS), np.float32)
    if rank == 0:
        return state
    geometry = data["geometry"].tocsr()
    denominator = np.asarray(geometry.sum(axis=0)).ravel()
    for value in range(len(CLASS_NAMES)):
        cells = np.flatnonzero(data["class_id"] == value)
        if not len(cells):
            continue
        numerator = np.asarray(geometry[:, cells].T @ fit.w[:, value, :rank])
        state[cells, :rank] = np.divide(
            numerator, denominator[cells, None], out=np.zeros_like(numerator),
            where=denominator[cells, None] > 0,
        )
    return state


def refine_cells(args, data, panel_index, composition, fit, library):
    program_count, support = install_new_programs(data, panel_index, composition, fit)
    context = proto.robust_standardize(np.column_stack((
        np.sqrt(np.maximum(composition, 0.0)),
        np.log1p(library)[:, None],
        np.log1p(np.asarray(data["geometry"].getnnz(axis=1)))[:, None],
    )))
    state_args = SimpleNamespace(
        spatial_folds=args.spatial_folds, permutations=args.permutations,
        spatial_buffer=args.spatial_buffer, graph_neighbors=10, graph_radius=3.5,
        boundary_strength=3.0, lambda_prior=0.25, lambda_graph=0.50,
        seed=args.seed,
    )
    graph_contrast, morph_contrast, edges, audits = proto.fit_states(
        state_args, data, context,
    )
    base_state = overlap_weighted_base_state(data, fit)
    state = base_state + graph_contrast + morph_contrast
    return state.astype(np.float32), graph_contrast, morph_contrast, edges, audits, program_count, support


def cell_probability(state, class_id, profiles, loadings, program_count):
    result = np.empty((len(class_id), profiles.shape[1]), np.float32)
    for value in range(len(CLASS_NAMES)):
        rows = np.flatnonzero(class_id == value)
        if not len(rows):
            continue
        k = int(program_count[value])
        logits = np.tile(np.log(np.maximum(profiles[value], EPS)), (len(rows), 1))
        if k:
            logits += np.clip(state[rows, :k] @ loadings[value, :k], -6.0, 6.0)
        logits -= logits.max(axis=1, keepdims=True)
        prob = np.exp(logits)
        prob /= np.maximum(prob.sum(axis=1, keepdims=True), EPS)
        result[rows] = prob.astype(np.float32)
    return result


def v55_obs_var(template: Path):
    backed = ad.read_h5ad(template, backed="r")
    obs = backed.obs.copy(); var = backed.var.copy()
    obsm = {key: np.asarray(backed.obsm[key]).copy() for key in backed.obsm.keys()}
    backed.file.close()
    return obs, var, obsm


def write_cell_output(path, template, matrix, semantics, audit, probability=None):
    obs, var, obsm = v55_obs_var(template)
    data = ad.AnnData(X=matrix, obs=obs, var=var)
    for key, value in obsm.items():
        data.obsm[key] = value
    if probability is not None:
        data.layers["log1p_cp10k"] = np.log1p(1.0e4 * probability).astype(np.float32)
    data.uns["schema"] = SCHEMA
    data.uns["X_semantics"] = semantics
    data.uns["audit_json"] = json.dumps(audit, ensure_ascii=False, default=json_default)
    write_h5ad_atomic(data, path)


def save_state(path, fit, state, graph_contrast, morph_contrast, edges, graph, audits):
    partial = path.with_suffix(path.suffix + ".partial")
    with h5py.File(partial, "w") as handle:
        handle.attrs["schema"] = SCHEMA
        handle.attrs["complete"] = 0
        handle.create_dataset("spot_program_activity", data=fit.w, compression="lzf")
        handle.create_dataset("program_loadings", data=fit.v, compression="lzf")
        handle.create_dataset("cell_program_activity", data=state, compression="lzf")
        handle.create_dataset("graph_contrast", data=graph_contrast, compression="lzf")
        handle.create_dataset("morph_contrast", data=morph_contrast, compression="lzf")
        handle.create_dataset("spot_graph_edges", data=graph.edges, compression="lzf")
        handle.create_dataset("spot_graph_edge_weight", data=graph.edge_weight, compression="lzf")
        handle.create_dataset("cell_graph_edges", data=edges, compression="lzf")
        handle.attrs["class_audit_json"] = json.dumps(audits, ensure_ascii=False, default=json_default)
        handle.attrs["complete"] = 1
        handle.flush()
    os.replace(partial, path)


def main() -> None:
    args = parse_args()
    for path in (args.reference_2um, args.mv2_h5ad, args.mv2_allocation,
                 args.spot_positions, args.mv3_state, args.panel_template):
        if not path.exists():
            raise FileNotFoundError(path)
    if args.output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite: {args.output_dir}")
    args.output_dir.mkdir(parents=True)
    if args.smoke:
        args.ranks = "0,2"; args.trend_lambdas = "0.2"; args.spatial_folds = 2
        args.cv_steps_v55 = 4; args.final_steps_v55 = 6
        args.cv_steps_2um = 4; args.final_steps_2um = 6; args.permutations = 2
    started = time.time()

    log("Loading immutable 2 um direct-cell observations")
    reference = ad.read_h5ad(args.reference_2um)
    counts2 = sparse.csr_matrix(reference.layers["measured_2um_sum_counts"], dtype=np.int32)
    coords2 = np.asarray(reference.obsm["spatial"], np.float64)
    class2 = reference.obs["cellvit_class_id"].to_numpy(np.int64)
    genes2 = reference.var["gene"].astype(str).to_numpy() if "gene" in reference.var else reference.var_names.astype(str).to_numpy()
    log("Building independent same-class 2 um graph")
    graph2 = bg.build_spatial_graph(coords2, group=class2, neighbors=6, radius_factor=4.0)
    final2, posterior2, cv2, selected2 = fit_2um(args, counts2, coords2, class2, graph2)
    output2 = ad.AnnData(X=posterior2, obs=reference.obs.copy(), var=reference.var.copy())
    for key in reference.obsm.keys():
        output2.obsm[key] = np.asarray(reference.obsm[key]).copy()
    output2.layers["measured_2um_sum_counts"] = counts2
    output2.layers["log1p_cp10k"] = np.log1p(1.0e4 * posterior2).astype(np.float32)
    output2.uns["schema"] = SCHEMA
    output2.uns["X_semantics"] = "independent BayesGraph empirical-Bayes latent gene probability; not UMI"
    output2.uns["raw_source_sha256"] = sha256(args.reference_2um)
    output2.uns["fit_audit_json"] = json.dumps({"selected": selected2, "cv": cv2, "graph": graph2.audit}, ensure_ascii=False, default=json_default)
    write_h5ad_atomic(output2, args.output_dir / "P1.2um.BayesGraph_latent_278.h5ad")
    del output2, posterior2, reference
    if final2.model is not None:
        final2.model.to("cpu")
    torch.cuda.empty_cache(); gc.collect()

    log("Loading V55 fractional spots, CellViT census and frozen geometry")
    data55, old_summary, panel_genes, panel_index, counts55, profiles55, composition55, cell_mass55 = load_v55(args)
    if not np.array_equal(genes2.astype(str), panel_genes.astype(str)):
        raise RuntimeError("2 um and V55 panel gene order differs")
    coords55 = data55["spots"][["hex_col", "hex_row"]].to_numpy(np.float64)
    context55 = np.column_stack((np.sqrt(np.maximum(composition55, 0.0)), np.log1p(cell_mass55)))
    log("Building V55 composition-aware spot graph")
    graph55 = bg.build_spatial_graph(
        coords55, component=data55["spot_component"], context=context55,
        neighbors=6, radius_factor=2.5,
    )
    final55, imputed55, cv55, selected55, folds55 = fit_v55(
        args, counts55, composition55, profiles55, coords55,
        data55["spot_component"], graph55,
    )
    log("Mapping newly learned spot programs to CellViT nuclei")
    state, graph_contrast, morph_contrast, cell_edges, class_audits, program_count, support = refine_cells(
        args, data55, panel_index, composition55, final55, counts55.sum(axis=1),
    )
    probability55 = cell_probability(state, data55["class_id"], profiles55, final55.v, program_count)
    log("Projecting raw and imputed V55 ledgers to cells")
    raw_mass, raw_audit = mass_core.exact_panel_mass(
        counts55, data55["geometry"], data55["capacity"], probability55,
    )
    imputed_mass, imputed_audit = mass_core.exact_panel_mass(
        imputed55, data55["geometry"], data55["capacity"], probability55,
    )
    audit55 = {
        "schema": SCHEMA, "selected": selected55, "cv": cv55,
        "spot_graph": graph55.audit, "program_count_by_class": program_count,
        "support_spots_by_class": support, "class_refinement": class_audits,
        "raw_mass_conservation": raw_audit, "imputed_mass_conservation": imputed_audit,
        "raw_spot_mass": float(counts55.sum()), "imputed_spot_mass": float(imputed55.sum()),
        "fixed_44_programs_used": False, "2um_expression_used_for_v55_fit": False,
        "external_reference_used": False,
    }

    spots = data55["spots"].copy()
    if "barcode" in spots:
        spots.index = spots["barcode"].astype(str)
    _, var55, _ = v55_obs_var(args.panel_template)
    spot_data = ad.AnnData(X=imputed55, obs=spots, var=var55.copy())
    spot_data.layers["raw_fractional_counts"] = counts55.astype(np.float32)
    spot_data.layers["imputed_probability"] = final55.probabilities.astype(np.float32)
    spot_data.uns["schema"] = SCHEMA
    spot_data.uns["X_semantics"] = "BayesGraph imputed V55 spot mass; per-spot library conserved"
    spot_data.uns["audit_json"] = json.dumps(audit55, ensure_ascii=False, default=json_default)
    write_h5ad_atomic(spot_data, args.output_dir / "P1.V55.BayesGraph_spot_imputed_278.h5ad")
    write_cell_output(
        args.output_dir / "P1.V55.BayesGraph_Lambda_278.h5ad", args.panel_template,
        probability55, "BayesGraph+nucleus cell latent probability; not UMI", audit55,
        probability=probability55,
    )
    write_cell_output(
        args.output_dir / "P1.V55.BayesGraph_X_raw_mass_278.h5ad", args.panel_template,
        sparse.csr_matrix(raw_mass), "exact raw V55 spot-gene-conserved fractional mass; audit", audit55,
    )
    write_cell_output(
        args.output_dir / "P1.V55.BayesGraph_X_imputed_mass_278.h5ad", args.panel_template,
        imputed_mass, "BayesGraph imputed spot mass allocated to cells; per-spot imputed ledger conserved", audit55,
        probability=probability55,
    )
    save_state(
        args.output_dir / "P1.V55.BayesGraph_state.h5", final55, state,
        graph_contrast, morph_contrast, cell_edges, graph55, class_audits,
    )
    write_tsv([
        {"modality": "2um", **{k: v for k, v in row.items() if k not in ("fit_audit", "tau_grid")}}
        for row in cv2
    ], args.output_dir / "P1.2um.BayesGraph_cv.tsv")
    write_tsv([
        {"modality": "V55", "rank": row["rank"], "lambda_trend": row["lambda_trend"],
         "mean_nll": row["mean_nll"], "se_nll": row["se_nll"],
         **{f"fold_{i}": value for i, value in enumerate(row["fold_nll"])} }
        for row in cv55
    ], args.output_dir / "P1.V55.BayesGraph_spatial_cv.tsv")
    summary = {
        "schema": SCHEMA, "status": "complete_panel_prototype_not_frozen_for_full_gene",
        "elapsed_seconds": time.time() - started, "n_genes": len(panel_genes),
        "n_2um_cells": len(class2), "n_v55_cells": len(data55["class_id"]),
        "n_v55_spots": len(counts55), "2um": {"selected": selected2, "cv": cv2, "graph": graph2.audit},
        "V55": audit55,
        "parameters": vars(args),
    }
    atomic_json(summary, args.output_dir / "P1.dual_BayesGraph.training_qc.json")
    log(f"Completed P1 MV4.2 BayesGraph prototype: {args.output_dir}")


if __name__ == "__main__":
    main()
