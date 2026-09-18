# 细胞拓扑生存模型

## 用途

模型用于结直肠癌切除标本 H&E 全切片的总体生存相对风险研究。同一患者的所有合格
切片组成一个病例包。模型不适用于诊断、治疗选择或绝对生存期预测。

## 输入

原始 WSI 模式依次完成组织掩膜、CellViT++ 细胞核检测、七分类头共识标注和物理坐标
细胞图构建。也可以直接输入拓扑 token HDF5。每个文件包含：

- 文件完整性标记及样本、病例标识
- 每个 250 µm 区域的 256 维拓扑特征
- 含未解析类别的细胞组成
- 区域邻接边

这些区域特征在生成时不使用生存结局。

## 模型结构

1. 细胞核作为图节点。边取最大 50 µm 范围内的对称物理 k 近邻并集。切片绝对位置、
   颜色和细胞原始形态不进入生存模型。
2. 区域 token 由软图传播、类型化边富集、距离滤过、细胞组成和邻域混合构成。
3. 全组织分支读取所有有效区域。肿瘤前沿分支先根据细胞组成确定肿瘤核心，再识别
   核心与非核心的邻接边界并向外扩展一跳。
4. 每个分支有 24 个拓扑原型。标准化后的区域特征被转换为软 RBF 隶属度。
5. 每个原型计算均值、标准差、90 分位数、最差切片 90 分位数和切片间异质性，得到
   每个分支 120 维的患者表征。
6. 两个岭惩罚 Cox 模型分别输出分支风险。风险按训练集的中位数和 IQR 标准化，最终
   分数取全组织分支与肿瘤前沿分支的等权平均。

## Fig3 结果

Fig3 的生存曲线来自嵌套五折 OOF 预测，并在每一层训练中控制信息泄漏：

| 指标 | 数值 |
|---|---:|
| 患者数 | 431 |
| 事件数 | 92 |
| C-index | 0.6454 |
| Log-rank P | `8.75e-05` |
| 每 IQR 风险增量 HR | 1.810 |
| 95% CI | 1.329 至 2.464 |

仓库中的部署权重是在全部合格开发病例上重新拟合的。上表报告的是 OOF 性能，不是训练集
上的表观分离。

## 质控与适用范围

- 原始 MPP 大于 0.75 µm/pixel 时，不进行细胞核级分析。
- 每张保留切片至少需要 4 个前沿区域，每名患者至少需要 16 个。
- 扫描仪、染色、制片和细胞类型识别的域偏移可能改变风险分布。
- 模型开发对象是 CRC 切除标本。活检、转移灶和非 CRC 组织需要独立验证。
- 分数是相对风险，不是个体生存概率。

<details>
<summary>English model card</summary>

## Intended use

This model estimates relative overall-survival risk from colorectal cancer resection
histology for research use. All eligible slides from one patient form one bag. The score
must not be used for diagnosis, treatment selection, or absolute survival-time prediction.

## Inputs and method

Raw WSI mode performs tissue masking, CellViT++ nucleus detection, consensus labelling
with seven classifier heads, and physical-coordinate graph construction. Prepared mode
accepts topology-token HDF5 files with 256 outcome-blind features per 250 µm region,
typed-cell composition, and region adjacency.

Nuclei are graph nodes. Edges are the symmetric union of physical k-nearest neighbours
within 50 µm. Whole-tissue and tumour-front branches each map region features to 24 soft
RBF prototype memberships. Five summaries per prototype produce a 120-dimensional
patient fingerprint for each branch. Ridge-penalized Cox heads return standardized branch
risks, and the final score is their equal mean.

## Validation and limitations

Nested five-fold OOF evaluation included 431 patients and 92 events. C-index was 0.6454,
log-rank P was `8.75e-05`, and HR per IQR of risk was 1.810
(95% CI 1.329 to 2.464). The distributed checkpoint is a refit on all eligible development
patients; these OOF values remain the performance estimate.

The pipeline rejects native MPP above 0.75 µm/pixel for nuclear analysis. It requires at
least four front regions per retained slide and sixteen per patient. Scanner, stain,
preparation, and phenotype-domain shifts may alter the risk distribution. Biopsies,
metastases, and non-CRC tissue are outside the validated scope.

</details>
