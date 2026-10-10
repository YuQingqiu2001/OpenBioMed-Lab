# CoVarST 本地发布包

这个包展示最终方案的三部分：独立 Programs 到癌种统一参考矩阵；H&E 到 spot 空转；仅以当前粗分辨率 ST 作为 RNA 来源的细胞精度化。

每个癌种只发布一个整合后的 H&E mapper checkpoint，共九个。它们保持原网络容量，在同一癌种所有折的预测概率均值上做集成蒸馏。166 个原始交叉验证权重不放进公开包。已通过整合一致性门槛的权重才进入 `models/`，状态见 [RELEASE_STATUS](docs/RELEASE_STATUS.md)。

在已有环境中运行：

```bash
cd CoVarST
export PYTHONPATH="$PWD/src"
python -m covarst verify-assets
python -m covarst infer-spots --cancer colorectal_cancer --features DATA/slide.npz --output WORK/slide.h5 --device cuda
```

`features` 输入必须是原流程的 ParamNet 标准化、Virchow2 CLS + patch mean 拼接所得 2560 维特征，NPZ 包含 `features`、`coords`、`barcode`。原始 H&E 输入可用 `extract-wsi`，第三方模型目录由用户明确提供，不自动下载。

细胞精度化允许使用真实 H&E 细胞核的位置、形态和特征，RNA profiles 和连续 Programs 均从当前低分辨率 ST 重新拟合，不能传入外部 scRNA-seq、已有细胞表达或用于评估的 2 μm 表达。细胞概率用于表达结构分析；守恒分配 `X_mass` 用于计数核账，两者的含义不同。

参考矩阵固定后才执行历史患者级 LOO；参考本身利用过全队列 RNA，因此这是固定参考下的预测评估，不能宣称严格 RNA 完全独立。新整合模型的统计是对教师集成的拟合一致性，不复用历史模型的 PCC 作为新模型成绩。细胞表达生物学真实性也需要独立验证。

[完整命令](docs/QUICKSTART.md) · [参考矩阵构建](docs/REFERENCE_PROGRAMS.md) · [细胞模块输入规范](docs/SPOT_TO_CELL.md)
