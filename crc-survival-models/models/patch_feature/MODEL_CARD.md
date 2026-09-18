# Patch Feature Survival Model / Patch 特征生存模型

## Intended use / 预期用途

Research-only overall-survival risk scoring from colorectal-cancer H&E whole-slide images. All slides belonging to one patient must be processed as one hierarchical bag.

用于 CRC H&E 全切片总体生存风险评分的科研模型。同一患者的全部切片必须作为一个层级病例包联合推理。

## Input contract / 输入契约

The WSI extractor samples up to 192 tissue regions per slide on a 250 µm physical grid. Each 224×224 RGB crop represents a 224 µm field of view at 1.0 µm/pixel. The cache stores a 1024-dimensional UNI embedding, local neighbour delta, nine-class coarse tissue probabilities, physical coordinates and tissue fractions.

WSI 提取器在 250 µm 物理网格上每张切片最多采样 192 个组织区域；每个 224×224 RGB Patch 在 1.0 µm/pixel 下对应 224 µm 视野。缓存保存 1024 维 UNI 表征、局部邻域差值、九类粗组织概率、物理坐标和组织比例。

## Model logic / 模型逻辑

1. **Patch encoder.** UNI embeddings and neighbour deltas are projected into a topology-guided stream and a morphology-free residual stream. The guided stream was trained to distil cell-topology features; only the free stream reconstructs UNI appearance.
2. **Translation-invariant spatial field.** For every patch, isotropic kernels at 250, 500 and 1000 µm summarize at most 12 neighbours. The network receives local embedding gradients, probability context, density and entropy. Absolute WSI coordinates are prohibited.
3. **Field and point representations.** A spatial field predicts topology mean and uncertainty; a point encoder supplies local topology and free texture. Their residual captures context-dependent spatial organization.
4. **Compartment weights.** Tumour probability and local tumour context define soft tumour-core, invasive-front and peritumour weights, attenuated by uncertainty.
5. **Outcome-blind sampling.** Per slide, 50% of retained regions prioritize invasive front, 25% tumour core and 25% peritumour. The patient cap is 768 regions and the slide cap is 192.
6. **Hierarchical MIL.** Gated attention pools regions separately by compartment, combines mean/hotspot/coverage information into slide representations, then applies slide attention and a patient-level survival head.
7. **Five-fold ensemble.** Each frozen fold produces a training-median/IQR standardized risk. The final patient score is the median of the five fold risks.

1. **Patch 编码器。** UNI 表征与邻域差值被拆分为拓扑引导流和形态自由残差流；引导流学习细胞拓扑，只有自由流重建图像表征。
2. **平移不变空间场。** 在 250、500、1000 µm 三个尺度上汇总最多 12 个邻居，输入局部表征梯度、组织概率上下文、密度和熵，禁止使用 WSI 绝对坐标。
3. **空间场与点表征。** 空间场预测拓扑均值和不确定性；点编码器提供局部拓扑与自由纹理；二者残差描述依赖上下文的空间组织。
4. **软分区权重。** 根据肿瘤概率及邻域肿瘤上下文构建肿瘤核心、侵袭前沿和瘤周权重，并用不确定性衰减。
5. **无结局采样。** 每张切片保留区域中约 50% 优先来自侵袭前沿、25% 来自核心、25% 来自瘤周；患者最多 768 个区域、单张切片最多 192 个。
6. **层级 MIL。** 分区门控注意力汇总区域，将均值/热点/覆盖率组成切片表征，再由切片注意力和患者级生存头输出风险。
7. **五折集成。** 每折风险按训练集风险中位数和 IQR 标准化，最终患者分数取五折中位数。

## Validation and the Fig3 curves / 验证与 Fig3 曲线

The frozen model generalized across TCGA-COAD, TCGA-READ, SR386 and a self-collected CRC cohort. Fig3 displayed cohort-specific outcome-adaptive cutpoints, which produced highly significant curves but are selection-inflated. They are kept only for exact figure reproduction. New slides use the continuous five-fold median risk and the frozen descriptive threshold `0`.

冻结模型在 TCGA-COAD、TCGA-READ、SR386 和自收集 CRC 队列中进行了验证。Fig3 为展示最强分层使用了各队列结局自适应切点，因此显著性存在选择性放大；这些切点仅用于图形复现。新切片使用连续五折中位风险及冻结阈值 `0`。

The continuous per-IQR random-effects meta-analysis was HR 1.743 (95% CI 1.176–2.583), P=`0.0206`. Detailed fixed-threshold and figure statistics are in `reference_results`.

连续风险每 IQR 的随机效应 Meta 分析为 HR 1.743（95% CI 1.176–2.583），P=`0.0206`；冻结阈值和图形统计量见 `reference_results`。

## QC and limitations / 质控与局限

- Requires valid physical coordinates and tissue-rich regions.
- Patch sampling is deterministic for a sample identifier but does not cover every tissue pixel.
- Colour, scanner and tissue-preparation domain shifts can affect UNI features.
- The score is not a calibrated survival probability or treatment-effect estimate.
- The bundled sklearn transformation objects must be opened only from this trusted repository and with the pinned sklearn version.

- 需要有效物理坐标及足够的组织区域。
- Patch 采样对样本 ID 可重复，但并不覆盖全部组织像素。
- 染色、扫描仪和制片域偏移会影响 UNI 表征。
- 风险分数不是校准后的生存概率，也不是治疗效应估计。
- sklearn 变换对象只应从本可信仓库加载，并使用锁定版本。
