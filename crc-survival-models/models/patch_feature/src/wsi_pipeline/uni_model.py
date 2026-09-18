"""Exact local UNI coarse classifier architecture and preprocessing."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import timm
import torch
from torch import nn

from config import UNI_CLASS_NAMES


IMAGENET_MEAN = np.asarray([0.485, 0.456, 0.406], dtype=np.float32)
IMAGENET_STD = np.asarray([0.229, 0.224, 0.225], dtype=np.float32)


class UNIHistologyClassifier(nn.Module):
    def __init__(self):
        super().__init__()
        self.backbone = timm.create_model(
            "vit_large_patch16_224",
            pretrained=False,
            init_values=1e-5,
            dynamic_img_size=True,
            num_classes=0,
        )
        self.classifier = nn.Sequential(
            nn.LayerNorm(1024),
            nn.Dropout(0.2),
            nn.Linear(1024, 512),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(512, len(UNI_CLASS_NAMES)),
        )

    def forward(self, images: torch.Tensor) -> torch.Tensor:
        return self.classifier(self.backbone(images))


def load_uni_classifier(checkpoint: Path, device: torch.device) -> UNIHistologyClassifier:
    model = UNIHistologyClassifier()
    payload = torch.load(checkpoint, map_location="cpu")
    state = payload.get("model_state_dict", payload.get("state_dict", payload))
    model.load_state_dict(state, strict=True)
    model.to(device).eval()
    return model


def preprocess_rgb(rgb: np.ndarray) -> torch.Tensor:
    if rgb.shape != (224, 224, 3):
        raise ValueError(f"UNI input must be 224x224 RGB, received {rgb.shape}")
    image = rgb.astype(np.float32) / 255.0
    image = (image - IMAGENET_MEAN) / IMAGENET_STD
    return torch.from_numpy(np.ascontiguousarray(image)).permute(2, 0, 1)
