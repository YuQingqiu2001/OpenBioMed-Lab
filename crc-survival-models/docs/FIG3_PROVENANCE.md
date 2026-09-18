# Fig3 模型与统计量

这份记录说明 `Fig3.png` 中的生存曲线对应哪一套模型，以及哪些统计量可以用于评价
新样本。模型身份已与保存的预测汇总和曲线统计表核对。

## C 面板：细胞拓扑模型

| 项目 | 数值 |
|---|---:|
| 终点 | 总体生存 |
| 患者 / 事件 | 431 / 92 |
| 低风险 / 高风险 | 215 / 216 |
| 评估方式 | 嵌套五折 OOF |
| Harrell C-index | 0.6454 |
| Log-rank P | `8.75e-05` |
| 每 IQR 风险增量 HR | 1.810 |
| 95% CI | 1.329 至 2.464 |

这条曲线使用全组织拓扑分支和肿瘤前沿拓扑分支的等权平均。它不是早期暂存目录中那套
独立的严格验证实验。

公开 checkpoint 是在全部合格开发病例上重新拟合的部署模型。上表仍以 OOF 结果作为
性能估计，不使用重拟合模型在训练病例上的分离程度。

## K/L 面板：Patch 特征模型

四个队列使用同一套冻结的五折模型。Fig3 中的分组切点由各队列的生存结局选择：

| 队列 | n | 事件 | HR (95% CI) | Log-rank P |
|---|---:|---:|---:|---:|
| TCGA-COAD | 439 | 97 | 2.262 (1.505 至 3.400) | `8.66e-05` |
| TCGA-READ | 161 | 28 | 3.363 (1.275 至 8.869) | `0.0142` |
| SR386 | 417 | 148 | 1.730 (1.250 至 2.395) | `0.000957` |
| 自收集 CRC | 294 | 100 | 1.744 (1.159 至 2.625) | `0.00762` |

图中分组结果的随机效应 Meta 分析为 HR 1.922（95% CI 1.365 至 2.706），
P=`0.00894`。

由于切点使用了生存结局，表中的显著性带有选择偏倚。这些结果只用于复现 Fig3。对新样本，
模型输出连续风险，并可用固定阈值 0 做描述性分组。连续风险每 IQR 的随机效应结果为
HR 1.743（95% CI 1.176 至 2.583），P=`0.0206`。

## 数据文件

- `models/cell_topology/reference_results/figure_survival_summary.json`
- `models/cell_topology/reference_results/model_comparison.tsv`
- `models/patch_feature/reference_results/figure_survival_statistics.tsv`
- `models/patch_feature/reference_results/figure_meta_analysis.tsv`
- `models/patch_feature/reference_results/fixed_threshold_statistics.tsv`
- `models/patch_feature/reference_results/continuous_risk_meta_analysis.tsv`

仓库只保留聚合统计量，不包含患者级预测或身份信息。

<details>
<summary>English summary</summary>

Panel C uses the equal-weight ensemble of the whole-tissue and tumour-front Cell Topology
branches. Nested five-fold OOF evaluation included 431 patients and 92 events. Harrell
C-index was 0.6454, log-rank P was `8.75e-05`, and HR per IQR of risk was 1.810
(95% CI 1.329 to 2.464). The public checkpoint is a refit on all eligible development
patients, while these OOF values remain the performance estimate.

Panels K and L use the same five frozen Patch Feature folds in all four cohorts. The
figure used cohort-specific outcome-adaptive cutpoints. Their random-effects estimate was
HR 1.922 (95% CI 1.365 to 2.706), P=`0.00894`, but this separation is inflated by cutpoint
selection. The portable continuous-risk estimate was HR 1.743 per IQR
(95% CI 1.176 to 2.583), P=`0.0206`. New patients should be evaluated with continuous
risk, with zero used only as a fixed descriptive split.

Only aggregate results are included in this repository.

</details>
