# 第三方软件与权重

源码和权重都从官方项目获取。本包不附带第三方权重，也不会自动下载。

| 组件 | 用途 | 官方来源 |
|---|---|---|
| BayesTME 1.0.0 | 逐切片原始计数 Program 分解 | [BayesTME](https://github.com/tansey-lab/bayestme) |
| ParamNet | H&E 染色标准化 | [ParamNet](https://github.com/khtao/ParamNet) |
| Virchow2 | CLS 与图块 token 均值特征 | [Virchow2](https://huggingface.co/paige-ai/Virchow2) |
| CellViT++ | 真实细胞核分割与特征 | [CellViT++](https://github.com/TIO-IKIM/CellViT-plus-plus) |
| torch-geometric | GATv2 图运算 | [PyTorch Geometric](https://github.com/pyg-team/pytorch_geometric) |

Virchow2 的获取和使用遵循官方模型说明与许可，本包许可不给第三方权重附加任何使用权。ParamNet 根目录下需要有 `source.model.ParamNet` 和 `checkpoints/ParamNet-Uni.pt`；Virchow2 根目录下需要放官方兼容权重。实际配准核普查按 CellViT++ 的文档生成，并把类别显式映射到本包的五类规范。

`_engine` 是本研究的数学实现和适配代码。`provenance/source_manifest.json` 记录了完整复制、精确语法树子集和路径调整的来源哈希；内部历史名称因兼容性保留，公开名称统一为 CoVarST。
