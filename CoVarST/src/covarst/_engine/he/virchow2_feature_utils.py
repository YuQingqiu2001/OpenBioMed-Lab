from __future__ import annotations

from pathlib import Path

import numpy as np
import torch
import timm
from safetensors.torch import load_file as load_safetensors
from timm.layers import SwiGLUPacked


def feature_subdir_name(embedding: str = "cls_mean") -> str:
    return f"virchow2_{embedding}_paramnet_fp16"


def load_virchow2_model(
    root: str | Path,
    device: torch.device,
    checkpoint: str = "model.safetensors",
    embedding: str = "cls_mean",
) -> tuple[torch.nn.Module, dict]:
    root = Path(root)
    if not root.exists():
        raise FileNotFoundError(root)
    ckpt_path = root / checkpoint
    if not ckpt_path.exists() and checkpoint == "model.safetensors":
        ckpt_path = root / "pytorch_model.bin"
    if not ckpt_path.exists():
        raise FileNotFoundError(ckpt_path)
    if embedding != "cls_mean":
        raise ValueError("Only embedding='cls_mean' is currently supported.")

    model = timm.create_model(
        "vit_huge_patch14_224",
        pretrained=False,
        mlp_layer=SwiGLUPacked,
        act_layer=torch.nn.SiLU,
        img_size=224,
        init_values=1e-5,
        num_classes=0,
        reg_tokens=4,
        mlp_ratio=5.3375,
        global_pool="",
        dynamic_img_size=True,
    )
    if ckpt_path.suffix == ".safetensors":
        state = load_safetensors(str(ckpt_path), device="cpu")
    else:
        state = torch.load(str(ckpt_path), map_location="cpu")
        if isinstance(state, dict) and "state_dict" in state:
            state = state["state_dict"]
    missing, unexpected = model.load_state_dict(state, strict=False)
    if missing or unexpected:
        raise RuntimeError(
            f"Virchow2 state_dict mismatch: missing={missing[:8]}, unexpected={unexpected[:8]}"
        )
    model.eval().to(device)
    meta = {
        "feature_model": f"Virchow2(local={root}, embedding={embedding})",
        "variant": "virchow2",
        "embedding": embedding,
        "checkpoint": str(ckpt_path),
        "token_dim": 1280,
        "feature_dim": 2560,
        "license": "CC-BY-NC-ND-4.0",
    }
    return model, meta


def virchow2_features(
    model: torch.nn.Module,
    imgs_uint8: np.ndarray,
    device: torch.device,
    mean: torch.Tensor,
    std: torch.Tensor,
) -> np.ndarray:
    arr = np.asarray(imgs_uint8)
    if arr.ndim != 4 or arr.shape[-1] != 3:
        raise ValueError(f"Expected NHWC uint8 RGB images, got shape={arr.shape}")
    x = torch.from_numpy(arr).to(device=device, dtype=torch.float32)
    x = x.permute(0, 3, 1, 2).contiguous().div_(255.0)
    x = (x - mean) / std
    amp_enabled = device.type == "cuda"
    with torch.inference_mode(), torch.autocast(device_type="cuda", dtype=torch.float16, enabled=amp_enabled):
        tokens = model(x)
        cls = tokens[:, 0]
        patch_mean = tokens[:, 5:].mean(dim=1)
        emb = torch.cat([cls, patch_mean], dim=-1)
    return emb.detach().to(torch.float16).cpu().numpy()
