# 使用指南

在 `CoVarST` 目录的 WSL 或 Linux 终端、已有兼容环境中执行。先设置 `export PYTHONPATH="$PWD/src"`。输入输出显式指定，输出需新路径；程序不自动安装依赖或下载权重。

## H&E → spot 级空间转录组

```bash
python -m covarst verify-assets
python -m covarst extract-wsi --wsi DATA/slide.svs \
  --paramnet-root THIRD_PARTY/ParamNet --virchow2-root THIRD_PARTY/Virchow2 \
  --output WORK/slide.features.npz --device cuda --batch-size 16
python -m covarst infer-spots --cancer colorectal_cancer \
  --features WORK/slide.features.npz --output WORK/slide.spot_expression.h5 --device cuda
```

已有符合训练表示的 2560 维 NPZ 特征时可跳过提取。规范见 [H&E 推理](HE_TO_SPOTS.md)。

| 参数 | 用途 |
|---|---|
| `--cancer` | 首页表格中的准确癌种标识 |
| `--features` | 含 `features`、`coords`、`barcode` 的 NPZ |
| `--assets` | 模型与参考所在的 CoVarST 根目录，通常可省略 |
| `--device` | `cpu` 或 `cuda` |
| `--batch` | 精确登记的 `SOURCE\|TECHNOLOGY` 标识；省略时使用中性参考 |
| `--gene-block` | spot 全基因解码块大小，默认 512 |
| `--mpp` | 有实测依据的微米/像素标度 |

MPP 优先从切片元数据读取，历史物镜倍率回退会记录来源。有确认的实测标度时才用 `--mpp` 覆盖。批次参数不能随意填入患者名或未知标识。

## 实测粗 ST → 细胞级表达

按[输入规范](SPOT_TO_CELL.md)准备实测整数计数、真实细胞核、物理捕获重叠矩阵和五类 PanNuke 类别。H&E 预测概率不能直接替代实测计数。

```bash
python -m covarst prepare-coarse --input DATA/section.coarse_and_nuclei.h5 \
  --sample SECTION --output WORK/section.prepared.h5
python -m covarst fit-coarse -- --sample SECTION --prepared WORK/section.prepared.h5 \
  --output-dir WORK/section.factors --ranks 0,8 --device cuda --max-cuda-gib 6
python -m covarst materialize-cells \
  --factors WORK/section.factors/SECTION.V55.BayesGraph_rank8_all_gene.factorized.h5 \
  --output WORK/section.P_post.h5
python -m covarst allocate-cell-mass --prepared WORK/section.prepared.h5 \
  --factors WORK/section.factors/SECTION.V55.BayesGraph_rank8_all_gene.factorized.h5 \
  --output WORK/section.X_mass.h5
```

`--ranks 0,8` 产生两个候选，应按粗分辨率空间交叉验证记录选择。上例的 `rank8` 只演示文件读取，不能预先认定它最优。精细表达成绩不能用于选模，秩不代表细胞类型数。

历史文件名保留 `V55` 以兼容拟合器；输入和摘要会区分真实粗 ST 与 virtual55 降采样。完整拟合参数用 `python -m covarst fit-coarse -- --help` 查看。

## 受控 HD 降采样

```bash
python -m covarst prepare-virtual55 --help
python -m covarst fit-virtual55 --help
python -m covarst prepare-native16 --help
python -m covarst fit-native16 --help
```

`prepare-virtual55` 接收原始 native16 计数、坐标、标度和真实核特征，独立拟合 virtual55 分数捕获质量，不接收精细表达。`prepare-native16 --mode coarse` 准备 native16 核普查和粗 RNA 谱；严格拟合器使用有内存上限的秩揭示投影器。精细表达只用于另外执行的评估。

## 构建新癌种参考

```bash
python -m covarst prepare-reference-counts --help
python -m covarst deconvolve-local --help
python -m covarst fit-fixed-reference --help
python -m covarst refine-reference --help
python -m covarst train-he --help
```

BayesTME 使用单独的 1.0.0 兼容环境。最终候选规模与选择参数必须显式指定，见[参考构建](REFERENCE_PROGRAMS.md)；早期脚本默认值不等于最终方案。新癌种需重新构建参考并训练自己的图像模型。

## 本地检查

```bash
python -m covarst verify-assets
python tests/check_software.py
python tests/check_generic_coarse.py
```

前者校验资产哈希，后两者检查数值约束、输入边界及粗 ST 到细胞输出的软件连接，记录写入 `validation/`。小型合成测试不能作为生物学验证。冻结脚本的原始帮助可能保留英文，命令和数据字段保持原名。
