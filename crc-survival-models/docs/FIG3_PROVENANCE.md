# Fig3 model provenance / Fig3 模型溯源

`Fig3.png` is the authoritative figure used to identify the release candidates. The mapping below was verified against the saved prediction summaries and curve statistics.

`Fig3.png` 是本次确认发布模型的权威依据；以下映射已用保存的预测汇总和生存曲线统计量逐项核对。

## Panel C: Cell Topology Survival Model / C 面板：细胞拓扑生存模型

- Endpoint: overall survival.
- Patients/events: 431/92.
- Low/high: 215/216.
- Leakage-controlled evaluation: nested five-fold out-of-fold prediction.
- Harrell C-index: 0.6454.
- Log-rank P: `8.75e-05` (displayed as `P<0.001`).
- Hazard ratio per interquartile risk increment: 1.810, 95% CI 1.329–2.464.

This curve is produced by the equal-weight ensemble of the whole-tissue and tumour-front topology branches. It is not the separate strict-validation experiment that had previously been placed in an earlier staging folder.

该曲线来自“全组织拓扑分支”和“肿瘤前沿拓扑分支”的等权集成，并不是先前误放入旧暂存目录的另一套严格验证实验。

The public numeric checkpoint is an all-eligible-development-patient refit for deployment. The values above remain the leakage-controlled out-of-fold estimates; they are not training-fit performance.

公开数值权重是在全部符合条件的开发病例上重新拟合的部署模型；上述性能仍来自无泄漏 OOF 评估，不是训练集拟合结果。

## Panels K/L: Patch Feature Survival Model / K/L 面板：Patch 特征生存模型

The same five-fold frozen model was applied to all four cohorts. The figure used cohort-specific outcome-adaptive cutpoints to display the strongest survival separation:

同一套冻结的五折模型应用于四个队列；图中为了展示最强分层，使用了各队列基于结局选择的切点：

| Cohort / 队列 | n | Events / 事件 | HR (95% CI) | Log-rank P |
|---|---:|---:|---:|---:|
| TCGA-COAD | 439 | 97 | 2.262 (1.505–3.400) | `8.66e-05` |
| TCGA-READ | 161 | 28 | 3.363 (1.275–8.869) | `0.0142` |
| SR386 | 417 | 148 | 1.730 (1.250–2.395) | `0.000957` |
| Self-collected CRC / 自收集 CRC | 294 | 100 | 1.744 (1.159–2.625) | `0.00762` |

Random-effects meta-analysis: HR 1.922 (95% CI 1.365–2.706), P=`0.00894`.

These cutpoints are statistically selection-inflated because survival outcomes were used to choose them. They are retained only to reproduce Fig3. Portable inference uses continuous risk; the frozen descriptive threshold is zero. Across cohorts, the continuous per-IQR random-effects estimate was HR 1.743 (95% CI 1.176–2.583), P=`0.0206`.

这些切点因使用生存结局进行选择而存在选择性放大，只用于复现 Fig3。对新样本的可迁移推理使用连续风险，冻结描述性阈值为 0。跨队列连续风险每 IQR 的随机效应 HR 为 1.743（95% CI 1.176–2.583），P=`0.0206`。

## Reproducibility files / 复现文件

- `models/cell_topology/reference_results/figure_survival_summary.json`
- `models/cell_topology/reference_results/model_comparison.tsv`
- `models/patch_feature/reference_results/figure_survival_statistics.tsv`
- `models/patch_feature/reference_results/figure_meta_analysis.tsv`
- `models/patch_feature/reference_results/fixed_threshold_statistics.tsv`
- `models/patch_feature/reference_results/continuous_risk_meta_analysis.tsv`

Only aggregate results are included; patient-level predictions and identifiers are deliberately excluded. / 仓库只保留聚合统计量，不包含患者级预测或身份信息。
