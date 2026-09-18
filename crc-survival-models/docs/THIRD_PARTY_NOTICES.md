# 外部依赖

本仓库不提供 UNI 或 CellViT++ 的源码、checkpoint、分类头和衍生主干权重。需要运行
原始 WSI 流程时，请从官方项目获取文件并遵守相应许可证：

- [CellViT++](https://github.com/TIO-IKIM/CellViT-plus-plus)
- [UNI](https://github.com/mahmoodlab/UNI)

官方仓库只说明文件来源。它们不一定包含本流程所需的项目特定分类头，用户仍需确认
checkpoint 与代码接口是否兼容。

## 本地路径

可以通过命令行参数提供文件，也可以设置环境变量：

| 环境变量 | 内容 |
|---|---|
| `CRC_SURVIVAL_CELLVIT_ROOT` | CellViT++ 本地仓库 |
| `CRC_SURVIVAL_CELLVIT_CHECKPOINT` | 兼容的 CellViT++ 分割 checkpoint |
| `CRC_SURVIVAL_CELLVIT_CLASSIFIERS` | 七个兼容表型分类头所在目录 |
| `CRC_SURVIVAL_UNI_CHECKPOINT` | 兼容的 UNI 衍生粗分类 checkpoint |

默认路径位于 `external/`。`.gitignore` 同时排除了根目录下的 `external/`、`vendor/`、
`weights/` 和所有 `.pth` 文件，避免把本地第三方文件提交到仓库。

本子项目的许可证只覆盖仓库原创代码、文档和随附的 CRC 生存模型，不改变用户另行取得
的 UNI 或 CellViT++ 材料的许可证。

<details>
<summary>English</summary>

This repository does not redistribute UNI or CellViT++ source, checkpoints, classifier
heads, or derived backbone files. Obtain the required files from the official
[CellViT++](https://github.com/TIO-IKIM/CellViT-plus-plus) and
[UNI](https://github.com/mahmoodlab/UNI) projects and follow their current licence and
access terms.

The official repositories identify the upstream sources but may not contain the exact
project-specific classifier heads expected by this pipeline. Provide compatible local
paths with the command-line options or the four environment variables listed above.

The licence in this subproject covers only the original code, documentation, and included
CRC survival models. It does not relicense separately obtained UNI or CellViT++ material.

</details>
