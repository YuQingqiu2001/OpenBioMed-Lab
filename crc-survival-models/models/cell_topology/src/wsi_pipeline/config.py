"""Frozen configuration for the outcome-blind CellViT++ topology pipeline."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import os
from pathlib import Path
from typing import Dict, Tuple


REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PACKAGE_ROOT = Path(__file__).resolve().parent
OUTPUT_ROOT = Path(
    os.environ.get("CRC_SURVIVAL_OUTPUT", str(Path.cwd() / "crc_survival_output"))
).resolve()


def external_path(variable: str, relative_default: str) -> Path:
    """Resolve a user-supplied third-party path without bundling its contents."""
    value = os.environ.get(variable)
    return Path(value).expanduser().resolve() if value else (REPOSITORY_ROOT / relative_default).resolve()


CELLVIT_ROOT = external_path("CRC_SURVIVAL_CELLVIT_ROOT", "external/cellvit-plus-plus")
CELLVIT_CHECKPOINT = external_path(
    "CRC_SURVIVAL_CELLVIT_CHECKPOINT", "external/cellvit/cell_segmentation_backbone.pth"
)
CELLVIT_WRAPPER = PACKAGE_ROOT / "run_cellvit.py"
CLASSIFIER_ROOT = external_path("CRC_SURVIVAL_CELLVIT_CLASSIFIERS", "external/cellvit/classifiers")

TCGA_INVENTORY = Path()
SR386_WSI_ROOT = Path()

UNI_CLASS_NAMES: Tuple[str, ...] = (
    "ADI",
    "BACK",
    "DEB",
    "LYM",
    "MUC",
    "MUS",
    "NORM",
    "STR",
    "TUM",
)

CLASSIFIER_FILES: Dict[str, str] = {
    "consep": "consep.pth",
    "lizard": "lizard.pth",
    "midog": "midog.pth",
    "nucls_main": "nucls_main.pth",
    "nucls_super": "nucls_super.pth",
    "ocelot": "ocelot.pth",
    "panoptils": "panoptils.pth",
}

LEVEL1_LABELS: Tuple[str, ...] = (
    "Epithelial",
    "Stromal",
    "Immune",
    "Dead_or_other",
    "Unresolved",
)

LEVEL2_LABELS: Tuple[str, ...] = (
    "Tumor_epithelial",
    "Non_tumor_epithelial",
    "Fibroblast_like",
    "Lymphocyte_TIL",
    "Plasma_cell",
    "Neutrophil",
    "Eosinophil",
    "Dead_apoptotic",
    "Spindle_stromal_candidate",
    "Macrophage_candidate",
    "Mitotic_tumor_candidate",
    "Other_unresolved",
)

TOPOLOGY_FEATURES: Tuple[str, ...] = ("E", "I", "P_L", "X", "R")
TOPOLOGY_LABELS: Tuple[str, ...] = ("A", "B", "C", "D")


@dataclass(frozen=True)
class FrozenConfig:
    schema_version: str = "cellular_interdigitating_topology_v3_tissue_first"
    seed: int = 20260830

    # Whole-slide tissue is identified and frozen before UNI or CellViT++.
    tissue_mask_mpp: float = 8.0
    tissue_candidate_fov_fraction_min: float = 0.20
    tissue_component_min_mm2: float = 0.05
    tissue_small_hole_max_mm2: float = 0.25
    tissue_closing_radius_um: float = 32.0

    # UNI is trained on 224 pixels representing approximately 112 micrometres.
    uni_target_mpp: float = 0.5
    uni_image_size_px: int = 224
    uni_stride_px: int = 112
    uni_tumor_threshold: float = 0.50
    uni_sensitivity_threshold: float = 0.35
    uni_tissue_fraction_min: float = 0.35
    uni_min_component_mm2: float = 0.10

    front_band_each_side_um: float = 500.0
    roi_read_halo_um: float = 100.0

    cellvit_target_mpp: float = 0.25
    cellvit_patch_size_px: int = 1024
    cellvit_overlap_px: int = 64
    duplicate_centroid_um: float = 3.0
    duplicate_contour_iou: float = 0.50

    level1_probability_min: float = 0.60
    level2_probability_min: float = 0.70

    local_knn_k: int = 6
    local_scale_min_um: float = 5.0
    local_scale_max_um: float = 30.0
    delaunay_scale_multiplier: float = 2.5
    delaunay_edge_max_um: float = 50.0
    alpha_primary_lambda: float = 1.6
    alpha_lambda_min: float = 0.8
    alpha_lambda_max: float = 2.4
    alpha_lambda_step: float = 0.2
    topology_raster_um: float = 5.0

    front_segment_arc_um: float = 100.0
    front_neighborhood_um: float = 250.0
    front_tissue_coverage_min: float = 0.80
    front_nuclei_min: int = 200

    topology_k: int = 4
    topology_cvar_fraction: float = 0.20
    topology_membership_threshold: float = 0.50
    high_confidence_pi: float = 0.40
    high_confidence_margin: float = 0.10
    prototype_bootstrap_n: int = 100

    def to_dict(self) -> dict:
        return asdict(self)


FROZEN = FrozenConfig()
