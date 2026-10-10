# H&E 推理为 spot 级空间转录组

每个癌种部署一个整合模型，架构保持原样：2560 维输入、256 隐藏宽度、6 邻居局部形态图、18 邻居上下文图、GATv2、RNA/组成双头，以及学习的支持门控。教师已经接受的有限 INR 修正被吸收进学生目标，部署时不再单独运行 INR。

## 图像特征

`extract-wsi` 沿用冻结的组织掩膜和物理六角网格：112 μm 图块视野、100 μm 最近邻间距、组织比例至少 0.20，且中心落在组织内。读取 RGB 后缩放到 224×224，经 ParamNet 标准化，再拼接 Virchow2 的 CLS 与图块 token 均值（排除四个 register token），得到 2560 维特征。

第三方权重目录见[官方来源](THIRD_PARTY.md)。本轮没有重跑完整的 WSI 主干提取。

| NPZ 字段 | 形状/类型 | 含义 |
|---|---|---|
| `features` | N × 2560，float16/float32 | 与训练一致的 ParamNet/Virchow2 表示 |
| `coords` | N × 2，float32 | 统一切片笛卡尔坐标 |
| `barcode` | N 个唯一字符串 | 精确特征行顺序 |
| `coords_um` | 可选 N × 2 | WSI 导出的微米位置 |

特征和坐标必须是有限值，条码唯一并逐行对齐。

## 校验与推理

先执行 `verify-assets`，校验清单内完整权重和参考文件的哈希。`infer-spots` 自己核对权重用途、癌种、一致性通过标记，以及权重登记的中性参考/基因顺序哈希，然后严格加载状态并重建冻结图网络；完整权重、分组索引和批次适配器的文件哈希由前置的 `verify-assets` 负责。

推理不读取样本 RNA、局部 theta 或实测表达。预测 Program 概率后，分块计算全部 10,000 个基因的 `theta_rna @ W`。默认使用中性 W；明确指定已登记的来源×技术标识时改用批次 W，遇到未知标识报错。

## 输出 HDF5

| 数据集 | 含义 |
|---|---|
| `barcode`、`spatial` | 对齐的标识和坐标 |
| `theta_rna` | 表达解码 Program 概率 |
| `theta_composition` | 组成分支 Program 概率 |
| `program_support_probability` | 支持概率 |
| `gene_name` | 固定基因标识顺序 |
| `expression_probability` | 全基因表达概率 |
| `expression_log1p_cp10k` | `log1p(10000 × expression_probability)` |

文件里记录了完成标记、参考/基因顺序哈希和 `transcriptomic_inputs_loaded=False`。这些概率估计的是表达组成，不能恢复实测文库总量或实测 UMI。

CellViT++ 用于独立的真实核普查和细胞模块，spot 解码不需要它。本包不附带 ParamNet、Virchow2 或 CellViT++ 权重，也不会自动下载。
