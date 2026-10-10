# CoVarST

CoVarST 是 OpenBioMed Lab（安徽医科大学）的空间转录组方法。公开包展示最终方案的三部分：**独立表达程序（Programs）合并为癌种统一参考矩阵、H&E 推理为 spot 级空间转录组，以及同一切片低分辨率 ST 的细胞级重建**。方法原名为 FSegGATv2，公开名称统一为 CoVarST。

本地已完成**九个癌种各一个整合权重**，共约 189.1 MB，保留原网络容量，全部通过对原始教师集成的一致性门槛。详细结果见[发布状态](docs/RELEASE_STATUS.md)和[本轮核查报告](docs/LOCAL_AUDIT.md)。

## 快速开始

下列命令在 WSL 或 Linux 终端、已有兼容环境中运行。运行不会自动安装依赖或下载第三方权重。

```bash
cd CoVarST
export PYTHONPATH="$PWD/src"
python -m covarst --help
python -m covarst list-models
python -m covarst verify-assets
python -m covarst infer-spots --cancer colorectal_cancer \
  --features DATA/slide.features.npz --output WORK/slide.spot_expression.h5 --device cuda
```

`verify-assets` 校验发布清单中的权重和参考文件。`infer-spots` 接受已有的 2560 维图像特征；从原始 H&E 开始时，先执行 `extract-wsi`。完整命令见[使用指南](docs/QUICKSTART.md)。输出必须使用新路径。

本地检查环境为 Python 3.9、PyTorch 2.7.1、torch-geometric 2.6.1、numpy 1.26.4、scipy 1.10.1、timm 0.9.16；BayesTME 1.0.0 使用单独的兼容环境。可按需执行 `pip install -e .`，也可直接设置 `PYTHONPATH`。`models/` 和 `references/` 需放在同一资产根目录；从其他位置运行时，可通过 `--assets /path/to/CoVarST` 指定。

在 Python 中也可调用相同入口：

```python
from covarst.he import infer_spots
result = infer_spots(
    "colorectal_cancer", "DATA/slide.features.npz",
    "WORK/slide.spot_expression.h5", assets="/path/to/CoVarST", device="cuda",
)
```

## 三个公开模块

| 模块 | 实际提供的内容 | 方法说明 |
|---|---|---|
| 癌种参考构建 | 原始计数预处理、逐切片 BayesTME、局部到共识映射、规模选择、联合细化、九套冻结参考 | [参考矩阵构建](docs/REFERENCE_PROGRAMS.md) |
| H&E → spot 级空间转录组 | 每癌种一个权重、ParamNet 标准化、Virchow2 特征、双头图网络、全基因解码 | [H&E 推理](docs/HE_TO_SPOTS.md) |
| 粗 ST → 细胞级表达 | 同切片 RNA 谱、类内连续变化轴、真实细胞核几何与特征、空间交叉拟合、测量零空间约束、表达概率和守恒核账 | [细胞级重建](docs/SPOT_TO_CELL.md) |

细胞模块的唯一 RNA 来源是**当前切片实测的粗分辨率 ST**，不用外部单细胞参考；真实 H&E 细胞核提供位置、形态和图像特征。H&E 推理概率不能直接作为 `prepare-coarse` 的实测计数输入。

## 癌种与权重

| 癌种 | 命令中的标识 | 参考 Programs 数 | 权重数 |
|---|---|---:|---:|
| 膀胱癌 | `bladder_cancer` | 48 | 1 |
| 乳腺癌 | `breast_cancer` | 44 | 1 |
| 结直肠癌 | `colorectal_cancer` | 44 | 1 |
| 皮肤鳞状细胞癌 | `cutaneous_squamous_cell_carcinoma` | 36 | 1 |
| 室管膜瘤 | `ependymoma` | 32 | 1 |
| 肾透明细胞癌 | `kidney_clear_cell_carcinoma` | 48 | 1 |
| 肺腺癌 | `lung_adenocarcinoma` | 48 | 1 |
| 胰腺癌 | `pancreatic_adenocarcinoma` | 48 | 1 |
| 前列腺癌 | `prostate_adenocarcinoma` | 48 | 1 |

九套参考共 396 个 Programs，每套含固定的 10,000 个基因。[模型整合](docs/MODEL_CONSOLIDATION.md)使用同癌种教师的概率集成目标得到一个部署模型。166 个原始患者留一权重保留在研究目录，未附带于本包。ParamNet、Virchow2、CellViT++ 等权重只提供[官方链接](docs/THIRD_PARTY.md)。

## 输出与检验口径

spot 输出为固定参考上的表达概率和 `log1p CP10K` 估计值。细胞 `P_post` 是按全基因归一化的推断表达分布；`X_mass` 是带未分配质量记录的分数计数核账。计数守恒只能证明核账一致。

原始患者留一成绩在 `provenance/` 单独保存。参考矩阵由原始全队列 RNA 构建，原始成绩属于固定队列参考条件下的图像泛化评估。新学生的一致性衡量对教师集成的拟合，不能替代独立患者或临床验证。通用粗 ST 适配器已做软件检查，细胞表达真实性仍需独立验证。

检查命令：

```bash
python -m covarst verify-assets
python tests/check_software.py
python tests/check_generic_coarse.py
```

小型合成测试仅检查软件和数值约束。真实切片检查、来源哈希和逐切片一致性分别见 `validation/`、`provenance/` 及[核查报告](docs/LOCAL_AUDIT.md)。本轮未重跑完整 WSI 特征提取或从原始计数重新构建参考。

## 文件结构

```text
configs/       运行配置与历史方案记录
models/        每癌种一个权重、Program 分组索引和模型说明
references/    固定参考 W、基因顺序、批次适配器和共识映射
src/covarst/   命令入口与冻结数学实现
provenance/    来源哈希、原始教师注册表和患者留一成绩
validation/    整合一致性、真实切片与软件检查记录
docs/          中文方法说明、输入规范、使用方法与核查范围
```

本包未附带患者 RNA 矩阵、H&E 原图或原始教师集合。冻结实现内部保留历史 `FSegGATv2` 类名以兼容加载；癌种标识、命令、字段和文件名保持稳定。

代码遵循 [CC BY-NC-SA 4.0](LICENSE)。第三方组件遵循各自许可和访问条件。论文发表后补充引文。当前仅在本地构建和核查，尚未推送 GitHub。

[中文阅读导引](GUIDE_ZH.md) · [运行规范](SKILL.md)
