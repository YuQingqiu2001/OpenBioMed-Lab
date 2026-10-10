# H&E 推理为 spot 级空间转录组

每癌种部署一个整合模型，保留原架构：2560 维输入、256 隐藏宽度、6 邻居局部形态图、18 邻居上下文图、GATv2、RNA/组成双头和学习的支持门控。学生目标吸收教师已接受的有界 INR 效果，不单独部署 INR。

## 图像特征

`extract-wsi` 沿用冻结组织掩膜和物理六角网格：112 μm 图块视野、100 μm 最近邻间距、至少 0.20 组织比例且中心在组织内。读取 RGB、缩放到 224×224，经 ParamNet 标准化后，拼接 Virchow2 CLS 和图块 token 均值，排除四个 register token，得到 2560 维特征。

第三方权重目录见[官方来源](THIRD_PARTY.md)。本轮未重跑完整 WSI 主干提取。

| NPZ 字段 | 形状/类型 | 含义 |
|---|---|---|
| `features` | N × 2560，float16/float32 | 与训练一致的 ParamNet/Virchow2 表示 |
| `coords` | N × 2，float32 | 统一切片笛卡尔坐标 |
| `barcode` | N 个唯一字符串 | 精确特征行顺序 |
| `coords_um` | 可选 N × 2 | WSI 导出的微米位置 |

特征和坐标需为有限值，条码唯一并逐行对齐。

## 校验与推理

先执行 `verify-assets` 校验清单内完整权重和参考文件哈希。`infer-spots` 自身核对权重用途、癌种、一致性通过标记及权重登记的中性参考/基因顺序哈希，严格加载状态并重建冻结图网络。完整权重、分组索引和批次适配器的文件哈希由前置 `verify-assets` 校验。

推理不读取样本 RNA、局部 theta 或实测表达。预测 Program 概率后，分块计算全部 10,000 基因的 `theta_rna @ W`。默认用中性 W；明确指定登记的来源×技术标识时选择批次 W，未知标识报错。

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

文件记录完成标记、参考/基因顺序哈希和 `transcriptomic_inputs_loaded=False`。概率估计表达组成，不能恢复实测文库总量或实测 UMI。

CellViT++ 用于独立的真实核普查和细胞模块，spot 解码不需要它。本包不附带 ParamNet、Virchow2 或 CellViT++ 权重，也不自动下载。
