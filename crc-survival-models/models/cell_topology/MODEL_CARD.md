# Cell Topology Survival Model / 细胞拓扑生存模型

## Intended use / 预期用途

Research-only relative overall-survival risk stratification for colorectal-cancer resection histology. The model accepts all eligible slides from one patient as one bag. It must not be used for diagnosis, treatment selection or an absolute survival-time estimate.

用于结直肠癌切除标本 H&E 全切片的科研性总体生存相对风险分层。同一患者的全部合格切片组成一个病例包。不得用于诊断、治疗决策或绝对生存期预测。

## Input contract / 输入契约

Raw WSI mode first performs tissue masking, CellViT++ nucleus detection, seven-head phenotype consensus and physical-coordinate cell-graph construction. Prepared mode accepts topology-token HDF5 files containing:

- complete-file flag and sample/case identifiers;
- 256 outcome-blind topology features per 250 µm region;
- typed-cell composition with an unresolved class;
- region adjacency edges.

原始 WSI 模式依次完成组织掩膜、CellViT++ 细胞核检测、七分类头共识标注和物理坐标细胞图构建。已有特征模式读取每个 250 µm 区域的 256 维无结局拓扑特征、细胞组成及区域邻接关系。

## Model logic / 模型逻辑

1. **Typed cell graph.** Nuclei are nodes; edges are a symmetric physical k-nearest-neighbour union limited to 50 µm. Absolute canvas position, colour and raw morphology do not enter the survival model.
2. **Outcome-blind regional topology.** Soft graph propagation, typed edge enrichment, distance filtrations, cell-composition and neighbourhood-mixing summaries form each region token.
3. **Two spatial scopes.** The whole-tissue branch uses all valid regions. The tumour-front branch derives tumour core from resolved-cell composition, keeps connected core components, detects core/non-core adjacency boundaries, and expands the front by one graph hop.
4. **Prototype dictionary.** Each branch contains 24 topology prototypes. Region features are standardized and converted into soft radial-basis memberships.
5. **Patient fingerprint.** For every prototype the model summarizes mean membership, standard deviation, 90th percentile, worst-slide 90th percentile and between-slide heterogeneity: 24 × 5 = 120 features per branch.
6. **Survival heads.** Two ridge-penalized Cox models produce branch risks. Each is normalized by its training median and IQR. Final risk is the equal mean of whole-tissue and tumour-front standardized risks.

1. **细胞图。** 细胞核为节点；边为最大 50 µm 的对称物理 k 近邻并集。生存模型不使用切片画布绝对位置、颜色或原始形态。
2. **无结局区域拓扑。** 通过软图传播、类型化边富集、距离滤过、细胞组成和邻域混合构建区域 token。
3. **双空间范围。** 全组织分支使用全部有效区域；肿瘤前沿分支由细胞组成确定肿瘤核心、保留连通核心、识别核心/非核心边界并扩展一跳。
4. **原型字典。** 每个分支包含 24 个拓扑原型，标准化区域特征被转换为软 RBF 隶属度。
5. **患者指纹。** 每个原型计算均值、标准差、90 分位数、最差切片 90 分位数和切片间异质性，每个分支共 120 维。
6. **生存头。** 两个岭惩罚 Cox 模型分别输出风险，经训练中位数和 IQR 标准化后等权平均。

## Validation / 验证

The Fig3 curve was based on leakage-controlled nested five-fold OOF predictions: n=431, events=92, C-index=0.6454, log-rank P=`8.75e-05`, HR per IQR=1.810 (95% CI 1.329–2.464). The deployment file is a full-development-cohort refit, so these OOF statistics—not apparent training separation—are the valid performance evidence.

Fig3 曲线来自嵌套五折无泄漏 OOF 预测：n=431、事件=92、C-index=0.6454、log-rank P=`8.75e-05`、每 IQR HR=1.810（95% CI 1.329–2.464）。部署权重为全开发队列重拟合，因此有效性能证据是上述 OOF 结果，而不是训练集表观分离。

## QC and limitations / 质控与局限

- Native MPP above 0.75 µm/pixel is rejected for nuclear analysis.
- At least four front regions per retained slide and sixteen per patient are required.
- Scanner, staining, specimen preparation and nucleus-phenotype domain shift can change risk distributions.
- The model was developed for CRC resection pathology; biopsy, metastasis and non-CRC tissue are out of scope unless independently validated.
- The score is relative and not calibrated to a patient-specific survival probability.

- 原始 MPP 大于 0.75 µm/pixel 时拒绝进行细胞核级分析。
- 每张保留切片至少需 4 个前沿区域，每名患者至少需 16 个前沿区域。
- 扫描仪、染色、制片及细胞类型识别的域偏移会改变风险分布。
- 开发对象为 CRC 切除病理；活检、转移灶和非 CRC 组织未经独立验证。
- 分数是相对风险，不是患者个体生存概率。
