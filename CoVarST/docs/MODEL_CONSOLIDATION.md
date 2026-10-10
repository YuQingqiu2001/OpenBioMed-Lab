# 每癌种一个原容量整合模型

原研究有 166 个患者留一图像模型，本包把它们整合成九个癌种各一个部署模型。隐藏宽度、网络容量和数值精度都保留，部署时只加载一个模型。

## 整合过程

1. 全部同癌种教师都在每张合格切片、同一个冻结 Program 坐标系上跑 H&E 推理；有被接受的有界 INR 修正时先应用修正。
2. 对 RNA 和组成的 softmax 概率分别跨折做等权算术平均。教师共用同一套 W、Program 和 10,000 基因顺序，概率可以对应，但不直接平均网络参数。
3. 从登记教师初始化原容量学生，匹配概率分布、每 Program 的中心化变化和局部空间差分。所有合格切片和 spots 都参与，特征与坐标按原始条码精确对齐。
4. 候选必须同时满足全部门槛，再按综合损失选权重。失败候选不进入公开目录。

共 302 张切片、166 位患者、498,060 spots；历史记录里的 303 张中有一张因原特征覆盖边界被排除。封装过程没有新增训练 spot 抽样。蒸馏不重新读取实测 RNA；教师和固定参考的原始训练用过 RNA。教师、特征和目标缓存保留在本地，不公开分发。

## 一致性门槛

| 指标 | 门槛 |
|---|---:|
| RNA JSD | ≤ 0.05 |
| 组成 JSD | ≤ 0.05 |
| 解码基因 PCC | ≥ 0.90 |
| 解码 spot 余弦 | ≥ 0.98 |
| 空间 theta 梯度余弦 | ≥ 0.85 |

解码评估覆盖全部 10,000 基因。逐切片记录后按患者/来源平衡汇总，详见 `validation/` 和各癌种的 `model_card.json`。选模时五项一起检查，包括空间梯度。

这里衡量的是整合输入上的教师集成拟合一致性。部分教师训练时见过对应患者的 RNA，因此不能当作新留出患者的评估，也不能证明新患者上非劣效或继承教师成绩。原始成绩见 `provenance/original_teacher_LOO_metrics.csv`，学生一致性见 `validation/student_fidelity_summary.csv`。权重记录里保存了参考/基因哈希、架构、教师数、训练和选模约定及结果；独立患者与临床验证仍需另行完成。

## 本地复现

需要完整的原研究目录：教师注册表、全部教师权重、已完成的患者留一推理/spot 索引，以及登记图像特征。这些原始输入不在公开包中，输出必须放在原研究目录之外。

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

先跑完 `targets` 再训练。可移植实现保留了损失和架构，只调整了作者路径，并保留满足全部门槛的候选。随机性和设备运算会让重训权重字节不同，发布模型以实际哈希和检查结果为准。
