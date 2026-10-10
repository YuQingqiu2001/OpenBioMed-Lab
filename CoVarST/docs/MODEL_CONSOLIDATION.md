# 每癌种一个原容量整合模型

原研究有 166 个患者留一图像模型，本包整合为九癌种各一个部署模型。隐藏宽度、网络容量与数值精度保留，部署只加载一个模型。

## 整合过程

1. 全部同癌种教师在每张合格切片、同一冻结 Program 坐标系上执行 H&E 推理，有接受的有界 INR 修正时先应用修正。
2. 对 RNA 与组成 softmax 概率分别跨折等权算术平均。教师共用 W、Program 和 10,000 基因顺序，概率可以对应；不直接平均网络参数。
3. 从登记教师初始化原容量学生，匹配概率分布、每 Program 中心化变化和局部空间差分。所有合格切片和 spots 参与，特征/坐标按原始条码精确对齐。
4. 候选必须同时满足全部门槛，再按综合损失选择权重。失败候选不进入公开目录。

共 302 张切片、166 位患者、498,060 spots；历史 303 张中一张因原特征覆盖边界排除。封装没有新增训练 spot 抽样。蒸馏不重新读取实测 RNA，教师和固定参考的原始训练使用过 RNA。教师、特征和目标缓存保留本地，不公开分发。

## 一致性门槛

| 指标 | 门槛 |
|---|---:|
| RNA JSD | ≤ 0.05 |
| 组成 JSD | ≤ 0.05 |
| 解码基因 PCC | ≥ 0.90 |
| 解码 spot 余弦 | ≥ 0.98 |
| 空间 theta 梯度余弦 | ≥ 0.85 |

解码评估覆盖全 10,000 基因。逐切片记录后按患者/来源平衡汇总，详见 `validation/` 和各癌种 `model_card.json`。选模同时检查五项，包括空间梯度。

这是整合输入上的教师集成拟合一致性。部分教师训练见过对应患者 RNA，不能作为新留出患者评估、证明新患者非劣效或继承教师成绩。原始成绩见 `provenance/original_teacher_LOO_metrics.csv`，学生一致性见 `validation/student_fidelity_summary.csv`。权重记录参考/基因哈希、架构、教师数、训练/选模约定和结果，独立患者与临床验证仍需另行完成。

## 本地复现

需要完整原研究目录：教师注册表、全部教师权重、已完成患者留一推理/spot 索引和登记图像特征。这些原始输入不在公开包中，输出必须位于原研究目录之外。

```bash
python scripts/consolidate_checkpoints.py --study-root ORIGINAL_STUDY --work-root WORK/consolidation --phase cache
python scripts/consolidate_checkpoints.py --study-root ORIGINAL_STUDY --work-root WORK/consolidation --phase targets
python scripts/consolidate_checkpoints.py --study-root ORIGINAL_STUDY --work-root WORK/consolidation --phase student
```

| 阶段 | 工作 |
|---|---|
| `cache` | 整理精确教师、特征和合格切片 |
| `targets` | 同癌种全教师概率目标 |
| `student` | 原容量训练、门槛检查和选择 |

先完成 `targets` 再训练。可移植实现保留损失与架构，调整作者路径，保留满足全部门槛的候选。随机性和设备运算可能使重训权重字节不同，发布模型以实际哈希和检查结果为准。
