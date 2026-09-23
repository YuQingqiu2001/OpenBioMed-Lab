# Patch 特征生存模型

## 用途

模型从结直肠癌 H&E 全切片计算总体生存相对风险，仅供科研使用。同一患者的全部切片
需要作为一个层级病例包联合推理。

## 输入

WSI 提取器在 250 µm 物理网格上采样，每张切片最多保留 192 个组织区域。每个
224 x 224 RGB Patch 在 1.0 µm/pixel 下对应 224 µm 视野。Patch cache 保存：

- 1024 维 UNI 表征
- 局部邻域表征差值
- 九类粗组织概率
- 物理坐标和组织比例

已有 Patch cache 时可以跳过 WSI 特征提取，直接运行生存模型。

## 模型结构

1. Patch 编码器把 UNI 表征和邻域差值拆成拓扑引导流与自由残差流。前者学习细胞
   拓扑，后者保留不能由拓扑解释的图像信息。
2. 空间场在 250、500 和 1000 µm 三个尺度上汇总最多 12 个邻居。输入包括局部表征
   梯度、组织概率上下文、密度和熵，不使用 WSI 绝对坐标。
3. 空间场输出拓扑均值和不确定性，点编码器给出局部拓扑与自由纹理。两者的残差描述
   与空间上下文有关的组织结构。
4. 肿瘤概率和邻域肿瘤上下文定义肿瘤核心、侵袭前沿和瘤周的软权重，不确定性会降低
   相应区域的权重。
5. 每张切片保留区域中，约 50% 优先来自侵袭前沿，25% 来自肿瘤核心，25% 来自瘤周。
   单张切片最多 192 个区域，每名患者最多 768 个。
6. 门控注意力先按分区汇总区域，再将均值、热点和覆盖率组成切片表征，最后经过切片
   注意力和患者级生存头得到风险。
7. 五个冻结模型分别按本折训练风险的中位数和 IQR 标准化。患者最终分数取五折风险
   的中位数。

## 验证结果

同一套五折模型用于 TCGA-COAD、TCGA-READ、SR386 和自收集 CRC 队列。连续风险每
IQR 的随机效应 Meta 分析为 HR 1.743（95% CI 1.176 至 2.583），P=`0.0206`。

## 质控与适用范围

- 输入需要可靠的物理坐标和足够的组织区域。
- Patch 采样由样本 ID 确定，结果可重复，但不会覆盖每一个组织像素。
- 染色、扫描仪和制片差异可能改变 UNI 表征及风险分布。
- 风险分数不是校准后的生存概率，也不是治疗效应估计。
- `.pkl` 中的 sklearn 变换对象只能从本仓库读取，并应使用 `environment.yml` 锁定的版本。

<details>
<summary>English model card</summary>

## Intended use

This research model estimates relative overall-survival risk from colorectal cancer H&E
whole-slide images. All slides from one patient must be processed as one hierarchical bag.

## Inputs and method

The WSI extractor samples up to 192 regions per slide on a 250 µm grid. Each 224 x 224
crop covers 224 µm at 1.0 µm/pixel. The cache stores a 1024-dimensional UNI embedding,
local neighbour deltas, nine coarse tissue probabilities, physical coordinates, and
tissue fractions.

The patch encoder separates a topology-guided stream from a free residual stream. A
translation-invariant field summarizes up to 12 neighbours at 250, 500, and 1000 µm.
Soft tumour-core, invasive-front, and peritumour weights guide region sampling and gated
attention pooling. Slide representations are combined at patient level. Five frozen folds
produce median/IQR-standardized risks, and the final score is their median.

## Validation and limitations

The same frozen model was evaluated in TCGA-COAD, TCGA-READ, SR386, and a self-collected
CRC cohort. Continuous risk gave a random-effects HR of 1.743 per IQR
(95% CI 1.176 to 2.583), P=`0.0206`.

Valid physical coordinates and tissue-rich regions are required. Sampling is
deterministic for a sample identifier but does not cover every tissue pixel. Stain,
scanner, and preparation shifts may affect the features. The score is not a calibrated
survival probability or a treatment-effect estimate.

</details>
