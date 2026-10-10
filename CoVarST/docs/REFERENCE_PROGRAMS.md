# 独立 Programs 如何形成癌种统一参考

最终产物是 `references/CANCER/consensus_reference_probability_10k.npy`，记作 W，形状为 K × 10,000，每行是一个 Program 的基因概率分布。`gene_order_10k.npy`、Program 顺序和批次适配器都是冻结资产。教师和学生都用这套参考解码 `theta_rna @ W`。

## 最终构建步骤

1. 按癌种整理实测原始计数。先覆盖标记基因，再按切片等权 logCP10K 方差挑选高变基因，拼成 10,000 基因面板。按最终配置排除线粒体、RPS 和 RPL 基因，并保存原始基因标识、顺序和切片来源。BayesTME 必须接收原始计数。
2. 每张切片独立运行原版 BayesTME 1.0.0，局部 K 取 2–12，保留后验 theta 和局部基因基底 W。不同切片的局部编号之间没有天然对应关系。
3. 拟合局部到共识的软映射和共享字典。最终候选数为 32、36、40、44、48，偏好 40。在验证基因 PCC 距最优不超过 0.002 的候选里，优先选最接近 40 的规模，保存成绩与映射。RCCD 提供映射机制，最终锁定划分的直接共识入口执行选择规则。
4. 联合细化共享 W、自由 spot theta，以及有界的技术/来源×技术适配器。目标函数结合 logCP10K Huber 损失、基因 PCC、余弦和弱 theta 锚定，不能简写成局部基底平均或新的 KL-NMF。适配器按平台/批次登记，不按患者身份设置。
5. 冻结 K、W、Program 与基因顺序、适配器数组及哈希，之后才训练和评估图像映射器，不逐个留一折重建参考。

## 已发布参考

| 癌种 | 癌种标识 | 固定 K | 权重数 |
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

共 396 个 Programs。乳腺癌保留历史的专用输入流程和 44 Program 资产；通用计数预处理示例覆盖后续八个癌种，只是输入模板，不能据此说九个癌种的原始路径完全相同。

## 入口与核对文件

| 命令 | 工作 |
|---|---|
| `prepare-reference-counts` | 原始计数与癌种基因面板 |
| `deconvolve-local` | 逐切片原版 BayesTME 分解 |
| `fit-fixed-reference` | 锁定划分下拟合并选择共识字典 |
| `refine-reference` | 最终联合细化 |
| `train-he` | 在固定参考下训练图像模型 |

各入口的 `--help` 会给出计数队列、spot 索引、划分清单和输出目录参数。最终入口需要显式设置：

```bash
python -m covarst fit-fixed-reference --queue-dir WORK/local_bayes \
  --master-spot-index WORK/master_spot_index.csv --split-manifest WORK/reference_split.csv \
  --output-dir WORK/shared_reference --candidates 32 36 40 44 48 \
  --target-programs 40 --selection-tolerance 0.002
```

`compress-programs` 是早期仅供 RCCD 使用的辅助入口，它的开发默认值不是最终方案；最终选择用上面的命令，再执行对应的联合细化。每份癌种参考都附带历史候选选择、细化摘要和配置。

| 文件 | 内容 |
|---|---|
| `local_to_consensus_mapping.npy` | 局部到共识映射 |
| `local_program_assignment.csv` | 局部分配记录 |
| `candidate_dictionary_selection.csv` | 候选规模选择依据 |
| `consensus_contract.json` | 参考构建约定 |
| `consensus_reference_probability_10k.npy` | 中性参考 W |
| `gene_order_10k.npy` | 精确基因顺序 |

患者 theta 目标没有公开，重建需要自行准备原始数据。

## 评估范围

历史固定参考在图像模型的患者留一评估之前就用上了全部可用的癌种 RNA，因此是队列条件下的传导式参考。原始成绩描述的是固定参考下的图像泛化，不代表参考对全新 RNA 队列的完全独立泛化。教师整合不改变这一范围。
