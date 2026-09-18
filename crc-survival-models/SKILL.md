---
name: crc-survival-models
description: Run the two frozen CRC pathology survival-risk models from prepared features or locally configured WSI dependencies.
---

# CRC Survival Models agent protocol

1. Read `README.md`, `docs/QUICKSTART.md`, and both model cards before inference.
2. Treat all slides supplied under one `case-id` as one patient; never score them as separate patients.
3. Prefer prepared-feature mode when topology-token or patch-cache HDF5 files already exist.
4. Before raw-WSI inference, run `python predict_wsi.py ... --check-only` and resolve every missing dependency.
5. Never guess micrometres-per-pixel. Use `--mpp` only when the scanner or acquisition record supplies it.
6. Do not download or redistribute UNI or CellViT++ assets automatically. Direct the user to the official links in `docs/THIRD_PARTY_NOTICES.md` and require local paths.
7. Report the continuous relative-risk score first. The frozen zero-threshold group is descriptive and is not a clinical decision boundary.
8. Never reuse Fig3 outcome-adaptive cutpoints for a new patient.

Prepared topology inference:

```bash
python predict.py topology slide_a.topology_tokens.h5 --case-id CRC_001 --output results/topology.json
```

Prepared patch inference:

```bash
python predict.py patch slide_a.patch_cache.h5 --case-id CRC_001 --output results/patch.json
```

本协议要求：同一患者多张切片必须合并输入；原始 WSI 运行前必须预检；禁止猜测 MPP；
禁止自动下载或再分发 UNI、CellViT++；结果优先解释连续相对风险，不把零阈值分组当作
临床决策界值，也不得把 Fig3 的结局自适应切点用于新患者。
