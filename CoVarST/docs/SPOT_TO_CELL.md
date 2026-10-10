# 不依赖外部单细胞参考的细胞级重建

唯一 RNA 来源为当前切片的**实测粗分辨率 ST**，同时使用真实细胞核、物理捕获几何、图像特征、PanNuke 大类和正的形态容量。外部 scRNA-seq、外部细胞谱、既有细胞表达或精细表达不能进入拟合。标记基因名称提供弱锚定，表达值仍从同一粗 ST 估计。

## 拟合与约束

拟合粗 RNA 类别谱、收缩项和每个受支持图像大类的 0–8 个连续变化轴。秩表示类内表达维度，与细胞类型数和上游固定 32–48 个 Programs 不同。

空间分块交叉验证、缓冲区、图像交叉拟合和置换控制约束拟合。原始细化使用图对比和已验证的连通分量测量零空间投影，保留绝对 **1e-9** 阈值及全部原始测量行检查。弱支持类别按原规则保持秩零。全局谱来自当前切片，后续空间验证为条件式/传导式，不能证明独立细胞表达真实性。

## 通用输入 HDF5

属性要求：

```text
schema = covarst_coarse_counts_and_image_nuclei_v1
complete = 1
rna_source = measured_coarse_st
external_single_cell_reference_used = False
real_image_nuclei = True
```

核对实测非负整数计数、有效正总量、唯一标识、物理坐标和同分量几何。H&E 预测概率不满足实测计数约定。

| 数据集 | 形状与含义 |
|---|---|
| `gene_name` | G 个唯一基因符号，顺序与矩阵一致 |
| `matrix/{shape,data,indices,indptr}` | CSC G × S 实测粗计数 |
| `geometry/{shape,data,indices,indptr}` | CSR S × N 真实核/Voronoi 捕获重叠权重 |
| `spots/barcode` | S 个唯一 spot 标识 |
| `spots/coords_um` | S × 2 微米笛卡尔位置 |
| `spots/component` | S 个非负测量域连通分量编号 |
| `cells/cell_id` | N 个唯一整数核标识 |
| `cells/class_id` | N 个 0–4 编号，按下表顺序 |
| `cells/coords_um` | N × 2 微米位置 |
| `cells/spatial_px` | N × 2 配准图像坐标 |
| `cells/features` | N × D 标准化真实核图像特征 |
| `cells/rna_capacity` | N 个正形态容量，不来自外部 RNA |
| `cells/owner_bin` | N 个分量内归属/插值 spot 索引 |
| `cells/domain_component` | N 个与归属 spot 一致的分量编号 |

| 编号 | PanNuke 名称 | 中文含义 |
|---:|---|---|
| 0 | Neoplastic | 肿瘤性细胞 |
| 1 | Inflammatory | 炎症细胞 |
| 2 | Connective | 结缔组织细胞 |
| 3 | Dead | 死亡细胞 |
| 4 | Epithelial | 上皮细胞 |

矩阵列、几何行和 spot 元数据同序；几何列与细胞元数据同序。保留实际重叠和无核正计数 spot，最近 spot 分配不能替代捕获几何。

通用适配器负责打包并从粗 ST 估计谱，拟合交给冻结 BayesGraph 与同一严格投影器。native16/virtual55 保留原物理仿射、有限 Voronoi 域和分数捕获积分。

## 输出含义

因子文件保存类别概率、全基因载荷和细胞状态。主要表达输出为：

```text
P_post = softmax_all_genes(log(class_gene_probability)
                         + cell_state @ class_program_loadings)
```

`materialize-cells` 分块计算，归一化覆盖全部基因，不能将子面板重新归一化后当作全转录组。`log1p_cp10k` 为推断概率的分析表示。

`allocate-cell-mass` 对每个原始 spot–基因观测执行冻结责任权重/IPF 分配，保存 `X_mass` 分数计数、spot–细胞–基因连接和无核 spot 未分配质量。守恒只证明核账一致，输出保留粗捕获印记，不能直接证明生物学单细胞分辨率；表达结构以 `P_post` 为主。

## 验证范围

历史 HD 降采样和图像诊断与生物学接受分开。本包不声称跨癌种的真实细胞身份准确度已验证。2 μm 精细表达如可用，只能另行评估，与准备、拟合和选模隔离。

通用任意粗分辨率适配器是封装扩展，已做软件和输入规范检查，不能继承历史 HD 队列生物学验证。运行见[使用指南](QUICKSTART.md)，证据见[核查报告](LOCAL_AUDIT.md)。
