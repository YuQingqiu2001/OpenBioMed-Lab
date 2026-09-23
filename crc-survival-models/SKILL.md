---
name: crc-survival-models
description: Score CRC overall-survival risk from prepared topology tokens, patch caches, or locally configured WSI dependencies.
---

# CRC Survival Models

Read `README.md`, `docs/QUICKSTART.md`, and the relevant model card before running a case.

## Rules

- Group every slide from one patient under the same `case-id` and score them together.
- Use prepared topology tokens or patch caches when they already exist.
- Run the WSI command with `--check-only` before starting feature extraction.
- Supply `--mpp` only when the scanner or acquisition record gives the value.
- Do not download or redistribute UNI or CellViT++ files. Point the user to
  `docs/THIRD_PARTY_NOTICES.md` and use local paths.
- Report continuous risk before the optional zero-threshold group.
- Do not apply the outcome-adaptive Fig3 cutpoints to a new patient.

Topology tokens:

```bash
python predict.py topology slide_a.topology_tokens.h5 \
  --case-id CRC_001 --output results/topology.json
```

Patch cache:

```bash
python predict.py patch slide_a.patch_cache.h5 \
  --case-id CRC_001 --output results/patch.json
```

执行要点：同一患者的切片一起输入；MPP 只能来自扫描记录；Fig3 的结局自适应
切点不能用于新患者。
