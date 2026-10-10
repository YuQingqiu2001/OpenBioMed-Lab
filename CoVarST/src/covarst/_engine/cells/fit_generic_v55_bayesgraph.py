#!/usr/bin/env python3
"""Fit the original BayesGraph functions to complete generic virtual55 inputs.

The sole expression input is a prepared coarse55 HDF5: exact circle/16um-square
fractional mass, independent coarse55 class profiles, and real image features.
No fine-resolution input argument or CRC panel fallback exists. Category names
and dimensions come from input metadata; original numerical functions remain
unchanged. Use --preflight-only before fitting into a fresh output directory.
"""
from __future__ import annotations

import argparse
import copy
import gc
from functools import partial
import importlib
import json
import os
from pathlib import Path
import re
import sys
import time

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import torch

import fit_native16_bayesgraph as native

SCHEMA = "generic_virtual55_BayesGraph_frozen_function_adapter_v1"
SOURCE_NAMES = ["p1_dual_bayesgraph_panel.py", "p1_bayesgraph_fullgene_rank_sweep.py",
                "p1_visium55_mv4_panel_prototype.py", "cell_resolved_visium55_mv3.py",
                "spot2cell_bayesgraph_core.py", "spot2cell_mv4_core.py"]
CELL_EXTRAS = ["source_cell_index", "domain_component", "direct_capture_overlap_fraction",
               "cross_component_fallback", "spatial_evidence_level"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True)
    p.add_argument("--prepared", type=Path, required=True)
    p.add_argument("--pipeline-dir", type=Path, default=native.SOURCE)
    p.add_argument("--output-dir", type=Path)
    p.add_argument("--checkpoint-dir", type=Path)
    p.add_argument("--preflight-only", action="store_true")
    p.add_argument("--preflight-json", type=Path)
    p.add_argument("--ranks", default="0,8")
    p.add_argument("--trend-lambdas", default="0.05,0.2")
    p.add_argument("--spatial-folds", type=int, default=5)
    p.add_argument("--spatial-buffer-um", type=float, default=150.0)
    p.add_argument("--cv-steps-v55", type=int, default=100)
    p.add_argument("--final-steps-v55", type=int, default=320)
    p.add_argument("--fit-batch-size", type=int, default=384)
    p.add_argument("--trend-batch-size", type=int, default=4096)
    p.add_argument("--cpu-threads", type=int, default=8)
    p.add_argument("--max-cuda-gib", type=float, default=16.0)
    p.add_argument("--permutations", type=int, default=100)
    p.add_argument("--gene-block", type=int, default=512)
    p.add_argument("--seed", type=int, default=20260828)
    p.add_argument("--device", default="cuda")
    a = p.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_.+-]+", a.sample):
        p.error("Sample ID must be a safe filename identifier")
    a.prepared = native.safe_input(a.prepared)
    ranks = [int(x) for x in a.ranks.split(",")]
    if not ranks or any(x < 0 or x > 8 for x in ranks) or len(ranks) != len(set(ranks)):
        p.error("Ranks must be unique integers in [0,8]")
    if a.spatial_buffer_um < 0 or a.spatial_folds < 2:
        p.error("Invalid spatial folds or physical buffer")
    if min(a.fit_batch_size, a.trend_batch_size, a.cpu_threads, a.cv_steps_v55,
           a.final_steps_v55, a.gene_block, a.permutations) < 1 or a.cpu_threads > 8:
        p.error("Positive sizes/steps required; CPU threads must be in [1,8]")
    if not 0 < a.max_cuda_gib <= 16:
        p.error("Additional Torch CUDA memory is limited to at most 16 GiB")
    if not a.preflight_only and a.output_dir is None:
        p.error("A fresh --output-dir is required for fitting")
    return a


def inspect_inputs(a):
    """Validate metadata, complete marker, and coordinate vectors before counts."""
    with h5py.File(a.prepared, "r") as h:
        if int(h.attrs.get("complete", 0)) != 1:
            raise ValueError("Prepared input is not complete; no data have been read")
        if h.attrs.get("schema") != "generic_hd_virtual55_independent_inputs_v1":
            raise ValueError("Unexpected prepared55 schema")
        provenance = json.loads(h.attrs["input_provenance_json"])
        if provenance.get("sample") != a.sample:
            raise ValueError("Prepared input sample does not match --sample")
        if provenance.get("2um_expression_used_for_fit") is not False:
            raise ValueError("Prepared provenance must explicitly exclude fine expression")
        if provenance.get("source_resolution_um") != 16 or provenance.get("virtual_spot_diameter_um") != 55:
            raise ValueError("Prepared input must be virtual55 from native16 square mass")
        a.coordinate_unit_um = float(h.attrs["coordinate_unit_um"])
        if a.coordinate_unit_um != 100.0:
            raise ValueError("Current adapter contract requires physical 100um coordinate units")
        a.spatial_buffer = a.spatial_buffer_um / a.coordinate_unit_um
        mpp = float(h.attrs["microns_per_pixel"])
        if not np.isfinite(mpp) or mpp <= 0 or not provenance.get("physical_affine"):
            raise ValueError("Missing physical pixel scale or native array affine provenance")
        genes = native.decode(native.read(h, "gene_name"))
        names = native.decode(native.read(h, "class_names"))
        panel = native.decode(native.read(h, "panel_gene_name"))
        if len(names) != 5 or len(set(names)) != 5:
            raise ValueError("Generic PanNuke contract requires five real input classes")
        if not 20 <= len(panel) <= 1000 or len(set(panel)) != len(panel):
            raise ValueError("Expected an explicit unique coarse-only panel of 20-1000 genes")
        if panel.tolist() != provenance.get("panel_genes"):
            raise ValueError("Panel identifiers disagree with prepared selection provenance")
        lookup = {gene: i for i, gene in enumerate(genes)}
        missing = [gene for gene in panel if gene not in lookup]
        if missing:
            raise ValueError(f"Panel genes missing from prepared matrix: {missing[:12]}")
        cell_id = native.read(h, "cells/cell_id").astype(np.int64)
        class_id = native.read(h, "cells/class_id").astype(np.int32)
        owner = native.read(h, "cells/owner_bin").astype(np.int32)
        cell_coord_epsilon = np.finfo(h["cells/coords"].dtype).eps
        spot_coord_epsilon = np.finfo(h["spots/coords"].dtype).eps
        coords = native.read(h, "cells/coords").astype(np.float64)
        cell_uv = native.read(h, "cells/native16_array_uv").astype(np.float64)
        spot_coords = native.read(h, "spots/coords").astype(np.float64)
        spot_uv = native.read(h, "spots/native16_array_colrow").astype(np.float64)
        component = native.read(h, "spots/component").astype(np.int32)
        profiles = native.read(h, "model/class_gene_probability").astype(np.float64)
        shape = tuple(int(v) for v in native.read(h, "matrix/shape"))
        geom_shape = tuple(int(v) for v in native.read(h, "geometry/shape"))
        if h["matrix"].attrs.get("format") != "csc" or h["geometry"].attrs.get("format") != "csr":
            raise ValueError("Prepared sparse matrix format mismatch")
        n, s, g = len(cell_id), len(spot_coords), len(genes)
        if n == 0 or s < a.spatial_folds or shape != (g, s) or geom_shape != (s, n):
            raise ValueError("Invalid matrix/geometry/census shapes")
        if len(set(cell_id)) != n or class_id.shape != (n,) or class_id.min() < 0 or class_id.max() >= len(names):
            raise ValueError("Invalid cell IDs or class schema")
        if owner.shape != (n,) or owner.min() < 0 or owner.max() >= s:
            raise ValueError("Invalid virtual55 cell owner indices")
        coordinate_errors = {}
        for actual, uv, expected, label, epsilon in [(coords, cell_uv, (n, 2), "cells", cell_coord_epsilon),
                                           (spot_coords, spot_uv, (s, 2), "spots", spot_coord_epsilon)]:
            if actual.shape != expected or uv.shape != expected or not np.isfinite(actual).all():
                raise ValueError(f"Invalid {label} coordinate dimensions/values")
            physical = uv * 16.0 / a.coordinate_unit_um
            # Prepared cell coordinates may be float32: tolerate only rounding
            # implied by their stored dtype, not a geometric re-registration.
            if not np.allclose(actual, physical, rtol=max(1e-10, 4 * epsilon), atol=1e-8):
                raise ValueError(f"{label} coordinates are not physical native16 array units")
            coordinate_errors[label] = float(np.max(np.abs(actual - physical)))
        if component.shape != (s,) or np.any(component < 0):
            raise ValueError("Invalid measured-domain spot components")
        if profiles.shape != (len(names), g) or not np.isfinite(profiles).all() or np.any(profiles < 0):
            raise ValueError("Invalid independent coarse55 class profiles")
        if not np.allclose(profiles.sum(1), 1., atol=2e-6):
            raise ValueError("Coarse55 class profiles must sum to one")
        extras = {key: native.read(h, "cells/" + key) for key in CELL_EXTRAS if "cells/" + key in h}
        for required in ["domain_component", "direct_capture_overlap_fraction", "cross_component_fallback"]:
            if required not in extras or extras[required].shape != (n,):
                raise ValueError(f"Missing cell evaluation-stratum vector: {required}")
        if not np.isfinite(extras["direct_capture_overlap_fraction"]).all() or np.any(extras["direct_capture_overlap_fraction"] < 0):
            raise ValueError("Invalid direct capture overlap fractions")
        feature_shape = h["cells/features"].shape
        if len(feature_shape) != 2 or feature_shape[0] != n or feature_shape[1] == 0:
            raise ValueError("Invalid real nucleus feature dimensions")
        if h["cells/spatial_px"].shape != (n, 2) or h["cells/rna_capacity"].shape != (n,):
            raise ValueError("Invalid pixel coordinates or RNA capacity dimensions")
        matrix_nnz = h["matrix/data"].size
    spots = pd.DataFrame({"hex_col": spot_coords[:, 0], "hex_row": spot_coords[:, 1]})
    data = {"gene_names": genes, "panel_genes": panel,
            "panel_index": np.asarray([lookup[g] for g in panel], np.int32),
            "cell_id": cell_id, "class_id": class_id, "owner_bin": owner,
            "coords": coords, "spots": spots, "spot_component": component,
            "cell_component": extras["domain_component"].astype(np.int32),
            "profiles": profiles, "cell_extras": extras}
    report = {"schema": SCHEMA, "sample": a.sample, "preflight_passed": True,
              "n_cells": n, "n_virtual55_spots": s, "n_genes": g,
              "n_panel_genes": len(panel), "class_names": names.tolist(),
              "class_schema": "input_metadata_five_PanNuke_classes", "ledger_nnz": matrix_nnz,
              "nucleus_features_shape": list(feature_shape),
              "coordinate_unit_um": a.coordinate_unit_um, "cv_buffer_um": a.spatial_buffer_um,
              "cell_graph_radius_um": 3.5 * a.coordinate_unit_um,
              "physical_coordinate_check_passed": True, "physical_coordinate_max_roundoff_model_units": coordinate_errors,
              "microns_per_pixel": mpp,
              "2um_expression_used_for_fit": False, "old_cell_expression_read": False,
              "old_programs_or_offsets_read": False, "CRC_panel_or_Lizard_labels_used": False,
              "input_provenance": provenance, "inputs": {"prepared": native.fingerprint(a.prepared)},
              "CV_interpretation": "conditional/transductive coarse55 model selection with fixed global profiles",
              "approx_bin_draws_per_virtual55_spot": {
                  "per_cv_fit": a.fit_batch_size * a.cv_steps_v55 / s,
                  "per_final_fit": a.fit_batch_size * a.final_steps_v55 / s},
              "environment": {"python": sys.version, "torch": torch.__version__, "numpy": np.__version__}}
    return report, data


def load_large_inputs(a, data):
    with h5py.File(a.prepared, "r") as h:
        if int(h.attrs.get("complete", 0)) != 1:
            raise ValueError("Prepared input complete marker changed")
        data["capacity"] = native.read(h, "cells/rna_capacity").astype(np.float64)
        data["spatial_px"] = native.read(h, "cells/spatial_px").astype(np.float64)
        features = native.read(h, "cells/features")
        if not np.isfinite(features).all():
            raise ValueError("Nucleus features contain nonfinite entries")
        data["features"] = np.clip(features, -6., 6.).astype(np.float64)
        for key, constructor in [("matrix", sparse.csc_matrix), ("geometry", sparse.csr_matrix)]:
            values = native.read(h, key + "/data").astype(np.float64)
            if not np.isfinite(values).all() or np.any(values < 0):
                raise ValueError(f"Invalid sparse {key} values")
            # Deliberately preserve fractional capture mass; no rounding or integer cast.
            data[key] = constructor((values, native.read(h, key + "/indices").astype(np.int32),
                                     native.read(h, key + "/indptr").astype(np.int64)),
                                    shape=tuple(native.read(h, key + "/shape")))
            data[key].check_format(full_check=True)
    if not np.isfinite(data["capacity"]).all() or np.any(data["capacity"] <= 0) or not np.isfinite(data["spatial_px"]).all():
        raise ValueError("Invalid capacity or pixel coordinates")
    mass = np.column_stack([np.asarray(data["geometry"] @ (data["capacity"] * (data["class_id"] == c))).ravel()
                            for c in range(len(data["profiles"]))])
    total = mass.sum(1)
    composition = np.divide(mass, total[:, None], out=np.zeros_like(mass), where=total[:, None] > 0).astype(np.float32)
    panel_counts = data["matrix"][data["panel_index"]].T.toarray().astype(np.float64)
    if np.count_nonzero((panel_counts.sum(1) > 0) & (total > 0)) < a.spatial_folds:
        raise ValueError("Insufficient positive occupied panel measurements")
    profiles = np.maximum(data["profiles"][:, data["panel_index"]], native.EPS)
    profiles /= profiles.sum(1, keepdims=True)
    return composition, total, panel_counts, profiles.astype(np.float32)


def main():
    a = parse_args()
    torch.set_num_threads(a.cpu_threads)
    torch.set_num_interop_threads(1)
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:a.cpu_threads])
    report, data = inspect_inputs(a)
    sys.path.insert(0, str(a.pipeline_dir))
    dual = importlib.import_module("p1_dual_bayesgraph_panel")
    sweep = importlib.import_module("p1_bayesgraph_fullgene_rank_sweep")
    bg = importlib.import_module("spot2cell_bayesgraph_core")
    report["runtime_compatibility"] = native.install_runtime_compatibility(dual.mv4)
    names = tuple(report["class_names"])
    dual.CLASS_NAMES = dual.mv3.CLASS_NAMES = dual.proto.CLASS_NAMES = sweep.mv3.CLASS_NAMES = names
    bg.fit_spot_program = partial(bg.fit_spot_program, batch_size=a.fit_batch_size,
                                 trend_batch_size=a.trend_batch_size)
    report["original_sources_sha256"] = {name: native.sha256(a.pipeline_dir / name) for name in SOURCE_NAMES}
    report["adapter_sha256"] = native.sha256(Path(__file__))
    report["shared_helper_sha256"] = native.sha256(Path(native.__file__))
    report["read_audit"] = native.READ_LOG.copy()
    if a.preflight_json:
        native.write_json(a.preflight_json, report)
    if a.preflight_only:
        print(json.dumps({k: v for k, v in report.items() if k != "read_audit"},
                         ensure_ascii=False, default=native.json_default, indent=2))
        return
    if a.output_dir.exists():
        raise FileExistsError(f"Output must be a fresh directory: {a.output_dir}")
    if a.device.startswith("cuda"):
        torch.cuda.set_per_process_memory_fraction(
            a.max_cuda_gib * (1 << 30) / torch.cuda.get_device_properties(0).total_memory)
    composition, total, panel_counts, profiles = load_large_inputs(a, data)
    a.output_dir.mkdir(parents=True)
    import bayesgraph_runtime
    bayesgraph_runtime.install_checkpoints(dual.mv4,
        a.checkpoint_dir or a.output_dir / "numerical_checkpoints", report["original_sources_sha256"])
    report["checkpoint_helper_sha256"] = native.sha256(Path(bayesgraph_runtime.__file__))
    native.write_json(a.output_dir / f"{a.sample}.virtual55.preflight.json", report)
    coords = data["spots"][["hex_col", "hex_row"]].to_numpy(np.float64)
    context = np.column_stack((np.sqrt(np.maximum(composition, 0)), np.log1p(total)))
    graph = bg.build_spatial_graph(coords, component=data["spot_component"], context=context,
                                   neighbors=6, radius_factor=2.5)
    started, results = time.time(), {}
    for rank in dual.parse_numbers(a.ranks, int):
        local = copy.copy(a)
        local.ranks, local.force_rank_v55 = str(rank), rank
        stage_started = time.time()
        final, _, cv, selected, folds = dual.fit_v55(local, panel_counts, composition,
            profiles, coords, data["spot_component"], graph)
        fit_finished = time.time()
        print(f"[{a.sample}] virtual55 rank{rank}: starting original cell refinement", flush=True)
        state, contrast, morph, edges, audits, program_count, support = dual.refine_cells(
            local, data, data["panel_index"], composition, final, panel_counts.sum(1))
        refine_finished = time.time()
        loadings = sweep.extend_v55_loadings(data, composition, final, data["panel_index"], a.gene_block)
        loadings[program_count == 0] = 0
        summary = {k: v for k, v in report.items() if k != "read_audit"}
        summary.update({"rank": rank, "selected": selected, "cv": cv, "spot_graph": graph.audit,
            "program_count_by_class": program_count, "support_spots_by_class": support,
            "class_refinement": audits, "parameters": vars(a),
            "factor_formula": "softmax(log(class_gene_probability)+cell_state@program_loadings)",
            "factor_is_count_conserving_allocation": False, "fixed_44_programs_used": False,
            "external_reference_used": False, "stage_seconds": {
                "CV_and_final_fit": fit_finished - stage_started,
                "cell_refinement": refine_finished - fit_finished,
                "all_gene_extension": time.time() - refine_finished}})
        output = a.output_dir / f"{a.sample}.V55.BayesGraph_rank{rank}_all_gene.factorized.h5"
        staging = output.with_name(output.name + ".adding_strata")
        sweep.write_factor(staging, modality=f"generic_V55_BG_rank{rank}", genes=data["gene_names"],
            cell_id=data["cell_id"], class_id=data["class_id"], spatial=data["spatial_px"],
            profiles=data["profiles"], loadings=loadings, state=state, program_count=program_count,
            summary=summary, sample=a.sample, extra={"spot_program_activity": final.w.astype(np.float32),
                "graph_contrast": contrast, "morph_contrast": morph, "spatial_fold": folds.astype(np.int8),
                "cell_graph_edges": edges.astype(np.float32), "panel_gene_index": data["panel_index"]})
        with h5py.File(staging, "r+") as h:
            h.attrs["complete"] = 0
            for key, values in data["cell_extras"].items():
                h["cells"].create_dataset(key, data=values, compression="lzf")
            h.create_dataset("class_names", data=np.asarray(names, object), dtype=h5py.string_dtype("utf-8"))
            h.attrs["complete"] = 1
        os.replace(staging, output)
        results[str(rank)] = {"path": str(output), "selected": selected}
        print(f"[{a.sample}] Completed virtual55 rank{rank}: {output}", flush=True)
        if final.model is not None:
            final.model.to("cpu")
        del final, state, contrast, morph, loadings
        torch.cuda.empty_cache(); gc.collect()
    native.write_json(a.output_dir / f"{a.sample}.virtual55.completed.json", {
        "schema": SCHEMA, "status": "complete", "results": results,
        "elapsed_seconds": time.time() - started, "parameters": vars(a),
        "read_audit": native.READ_LOG, "2um_expression_used_for_fit": False})


if __name__ == "__main__":
    main()
