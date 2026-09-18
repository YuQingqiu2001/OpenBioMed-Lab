# CRC Survival Models / 结直肠癌生存风险模型

This repository packages the two pathology models used for the survival curves in `Fig3.png`:

1. **Cell Topology Survival Model** — a morphology-free model built from typed-cell neighbourhood topology.
2. **Patch Feature Survival Model** — a topology-guided spatial neural-field model built from histology patch features.

本仓库打包了 `Fig3.png` 生存曲线实际使用的两套病理模型：

1. **细胞拓扑生存模型**：只使用细胞类型及邻接拓扑，不使用图像颜色、形态、临床变量或微生物数据。
2. **Patch 特征生存模型**：将病理 Patch 表征、局部空间梯度及肿瘤核心/边界/瘤周结构联合建模。

Public model names and folders intentionally contain no historical development labels. Immutable internal HDF5 schema fields may retain technical schema identifiers for backward compatibility. / 对外模型名和文件夹不再包含历史开发标签；个别 HDF5 内部 schema 标识为兼容旧特征文件而保留。

## What a prediction means / 预测结果含义

Both models output a **continuous relative risk score** for overall survival and a frozen descriptive group:

- `risk >= 0`: high risk
- `risk < 0`: low risk

两套模型均输出总体生存的连续相对风险分数，并按冻结阈值 `0` 给出高/低风险组。分数不是生存概率、剩余寿命或临床诊断，不能直接替代病理分期或治疗决策。

The highly significant curves shown in Fig3 are preserved under `reference_results`. Some Patch-model figure panels used outcome-adaptive cohort-specific cutpoints to visualize the strongest separation. Those cutpoints are **figure-reproduction statistics only** and are not used for a new patient. New slides always receive the continuous risk and the frozen zero threshold.

Fig3 中高度显著的曲线统计量保存在 `reference_results`。其中部分 Patch 模型曲线使用了各队列基于结局选择的最优切点，这些切点仅用于复现图形，不能作为新患者的通用阈值。新切片只输出连续风险和冻结的零阈值分组。

## Quick use / 快速使用

```bash
git clone https://github.com/YuQingqiu2001/OpenBioMed-Lab.git
cd OpenBioMed-Lab/crc-survival-models
conda env create -f environment.yml
conda activate crc-survival
python predict.py validate-install
```

The repository includes both frozen CRC survival models. It intentionally does **not**
redistribute UNI or CellViT++ source code and checkpoints. Prepared-feature inference
works with the files in this subproject; raw-WSI inference additionally requires locally
obtained compatible third-party checkpoints. / 本子项目包含两套冻结的 CRC 生存模型，
但有意不再分发 UNI 或 CellViT++ 的源码与权重。已有特征推理可直接使用；原始 WSI
推理还需用户从官方来源自行获取兼容的第三方文件。

Official upstream repositories / 官方原库：

- CellViT++: https://github.com/TIO-IKIM/CellViT-plus-plus
- UNI: https://github.com/mahmoodlab/UNI

Prepared-feature inference / 已有特征文件直接推理：

```bash
python predict.py topology slide_a.topology_tokens.h5 slide_b.topology_tokens.h5 \
  --case-id CRC_001 --output results/topology.json

python predict.py patch slide_a.patch_cache.h5 slide_b.patch_cache.h5 \
  --case-id CRC_001 --device cuda:0 --output results/patch.json
```

Raw WSI inference / 原始切片端到端推理：

```bash
python predict_wsi.py patient_slide.svs \
  --case-id CRC_001 --model both --output-dir results/CRC_001 \
  --cellvit-root /path/to/CellViT-plus-plus \
  --cellvit-checkpoint /path/to/cellvit_checkpoint.pth \
  --cellvit-classifiers /path/to/cellvit_classifiers \
  --uni-checkpoint /path/to/compatible_uni_classifier.pth
```

For multiple slides from the same patient, pass all slide paths in one command. Never run them as independent patients. / 同一患者多张切片必须在同一条命令中输入，不能拆成多个患者。

Both branches accept standard OpenSlide-compatible formats (`SVS`, `MRXS`, `NDPI`, `TIF/TIFF`, `SCN`, `BIF`). The custom `IBL` reader is available for the Patch branch only. Unsupported or QC-failing slides stop with an explicit error rather than receiving a fabricated score. / 两条模型链均支持 OpenSlide 兼容格式；自定义 `IBL` 读取器仅用于 Patch 分支。不支持或质控失败的切片会明确停止，不会强行生成分数。

Run `--check-only` before a long WSI job. The cell-topology branch requires a CUDA GPU with approximately 24 GB VRAM, at least 32 GB RAM, native slide resolution no worse than 0.75 µm/pixel, and adequate free disk space. The Patch branch also benefits strongly from CUDA. Windows users should run the WSI pipeline in WSL2/Linux; prepared-feature inference is platform neutral.

长任务前可先运行 `--check-only`。细胞拓扑分支建议使用约 24 GB 显存的 CUDA GPU、至少 32 GB 内存，并要求原始切片分辨率不低于 0.75 µm/pixel；Patch 分支也建议使用 CUDA。Windows 用户建议在 WSL2/Linux 中运行原始 WSI 流程。

## Repository map / 目录说明

- `models/cell_topology/`: topology inference code, numeric deployment weights, feature pipeline and Fig3 reference statistics.
- `models/patch_feature/`: five frozen folds, patch-cache extractor, inference code and multi-cohort reference statistics.
- `docs/`: bilingual usage, provenance, validation and release documentation.
- `release/`: package inventory, checksums and automated audit.
- `predict.py`: prepared-feature and installation-check CLI.
- `predict_wsi.py`: end-to-end WSI CLI.
- `docs/FIG3_PROVENANCE.md`: exact curve-to-model mapping and statistical caveats.
- `docs/THIRD_PARTY_NOTICES.md`: official upstream links and local setup contract.

`external/`, `vendor/`, `weights/` and all `*.pth` files are ignored deliberately so
third-party assets cannot be committed accidentally. / 上述外部目录和权重扩展名已加入
忽略规则，防止第三方资产被误提交。

## Validation status / 验证状态

- Cell-topology prepared-feature inference: smoke-tested on a two-slide CRC patient.
- Patch prepared-feature inference: smoke-tested on a real external cache; risk-group assignment matched the historical result and continuous risk differed by less than `0.001` because of CPU/GPU numerical variation.
- Raw-WSI CLI: end-to-end validated from a fresh 304 MB CRC SVS through both feature pipelines and both survival models. The run produced 236,169 consensus-typed cells and 192 patch regions.
- Public-repository audit confirms that no UNI or CellViT++ source/checkpoint is included.
- The package is for research use only and is not a medical device.

- 细胞拓扑特征级推理：已在一个双切片 CRC 病例上通过冒烟测试。
- Patch 特征级推理：已在真实外部缓存上通过测试；风险分组与历史结果一致，连续分数的 CPU/GPU 数值差异小于 `0.001`。
- 原始 WSI 入口：已使用一张 304 MB CRC SVS 从头跑完两条特征链和两套生存模型；共生成 236,169 个共识标注细胞和 192 个 Patch 区域。
- 公开仓库审计确认未包含 UNI 或 CellViT++ 源码及权重。
- 本仓库仅限科研用途，不构成医疗器械或临床建议。

See [Quickstart](docs/QUICKSTART.md), [Cell Topology model card](models/cell_topology/MODEL_CARD.md) and [Patch Feature model card](models/patch_feature/MODEL_CARD.md) before use. / 使用前请阅读快速开始和两份模型卡。
